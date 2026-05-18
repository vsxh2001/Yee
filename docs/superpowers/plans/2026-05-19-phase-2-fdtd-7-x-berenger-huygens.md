# Phase 2.fdtd.7.x — Berenger Huygens-surface subgridding closure — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan track-by-track.

**Companion spec:** `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`
**Base SHA:** `618b49b` (`main` HEAD: Phase 2.fdtd.7 Q1–Q5 + Q4.1 merged; Q5 strict 0.5% gate still `#[ignore]`'d).
**Target phase:** 2.fdtd.7.x only — the closure-only rewrite. 7.1–7.5 stay deferred per spec `2026-05-18` §9.
**Tech-stack additions:** none. Reuses the existing `SubgridRegion` / `SubgriddedSolver` API surface, the Q3 snapshot + interpolation helpers, the Q4.1 `snapshot_fine_h_mid_step` helper, and the Q1 `WalkingSkeletonSolver` per-stage helpers (`update_h_only`, `update_e_only`, `apply_cpml_h`, `apply_cpml_e`, `advance_clock`). No new dependencies; the `#![forbid(unsafe_code)]` floor is preserved.

---

## Goal

Replace the spec `2026-05-18` §3 stage-6 (`average_fine_h_to_coarse`) and stage-7 (`overwrite_coarse_e_from_fine`) bidirectional direct-copy closure — diagnosed unstable by VVVVVV (`a2abb4c`) — with a Berenger 2006 Huygens-surface scheme that injects equivalent `J = +n̂ × H_tot` and `M = −n̂ × E_tot` surface currents from the fine subdomain into the coarse parent. Clear the Q5 strict 0.5%-of-peak plane-wave traversal gate (un-`#[ignore]` it), introduce a new Q6 10 000-step round-trip energy-drift gate (`|W(N) − W(0)| / W(0) ≤ 0.5%`), and unblock Q7 fdtd-007 Maloney-Smith production gate (forward-port from the original 7.0 plan unchanged).

## Pre-flight

No new dependencies. No refactor required — the Q1 step-refactor (commit `99db783`) already provides the composable `update_h_only` / `update_e_only` / `apply_cpml_*` / `advance_clock` surface that the Berenger pipeline needs. Step B1 below is the first behaviour-changing step.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-fdtd/src/subgrid.rs` | Modify | Add `inject_equivalent_currents_to_coarse` + Huygens-face index enumeration helpers. Mark `average_fine_h_to_coarse` / `overwrite_coarse_e_from_fine` `#[doc(hidden)]` (retained for posterity per ADR-0035). Rewrite `SubgriddedSolver::step` and `step_with_gaussian_source_ez` to the spec §3 stage list. |
| `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` | Modify | Remove the `#[ignore]` attribute on the strict 0.5%-of-peak test; tighten the assertion message to reference Berenger 2006. |
| `crates/yee-fdtd/tests/subgrid_energy_balance.rs` | Create | Q6 10 000-step round-trip energy-drift gate. |
| `crates/yee-fdtd/validation/README.md` | Modify | Add the `fdtd-007 (stability, Berenger closure)` row and forward-port the `fdtd-007 (Maloney-Smith)` row. |
| `crates/yee-validation/src/lib.rs` | Modify (Step B5 only) | `run_fdtd_007_maloney_smith_slot` driver — forward-ported unchanged from spec `2026-05-18` Q7. |
| `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` | Create (Step B5 only) | fdtd-007 production-gate test. |

No changes to `yee-core`, `yee-mesh`, `yee-cli`, `yee-py`, `yee-gui`, `yee-cuda`. Lane is `crates/yee-fdtd/**` for B1–B4 and adds `crates/yee-validation/**` for B5 only.

## Step ladder

### Step B1 — `inject_equivalent_currents_to_coarse` skeleton + Huygens-face index enumeration

- **Brief:** Add `pub fn inject_equivalent_currents_to_coarse(&self, parent: &mut YeeGrid)` to `SubgridRegion`. Body enumerates the six Huygens faces (`±x`, `±y`, `±z`) and computes the per-face equivalent surface currents `J = +n̂ × H_tot` (sampled from the Q4.1 `snapshot_fine_h_mid_step` mid-step `H` storage) and `M = −n̂ × E_tot` (sampled from the post-`update_fine_e` outer-layer fine `E_t`). Currents are added as RHS source terms to the coarse `E_t` (from `J`) and coarse `H_t` (from `M`) on the interface plane: `E_coarse += (dt / ε_0) · J_S`, `H_coarse += (dt / μ_0) · M_S`. Half-open tangential index range per face assigns each cuboid edge cell to exactly one of two adjacent faces (axis-ordering rule: lower-numbered axis owns the shared edge). Not yet wired into `SubgriddedSolver::step` — that is Step B2.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs`.
- **Base SHA dep:** `618b49b` (main HEAD).
- **DoD:** function present and `pub`; rustdoc cites Berenger 2006 §III by full reference; unit tests `closed_surface_constant_j_recovers_discrete_divergence` (constant `J` on the closed Huygens surface integrates to discrete-divergence-free) and `cuboid_edge_owned_by_one_face_only` (no edge cell receives a contribution from both adjacent faces) both pass; existing `average_fine_h_to_coarse` and `overwrite_coarse_e_from_fine` are marked `#[doc(hidden)]` with a doc-comment pointer to ADR-0035 but otherwise unchanged.
- **Verification:** `cargo clippy -p yee-fdtd --all-targets -- -D warnings && cargo test -p yee-fdtd --release subgrid::inject` exits 0.
- **Escape hatch:** blocked > 15 min on the per-face index map → write out the `+x` face only with three unit tests (`J_y`, `J_z`, edge-ownership), copy-paste to the remaining five with a sign-table comment. Do not attempt a general 6-face dispatch on first pass.
- **LOC:** ~280.

### Step B2 — `SubgriddedSolver::step` Berenger pipeline rewrite

- **Brief:** Rewrite `SubgriddedSolver::step` and `step_with_gaussian_source_ez` to the spec §3 stage list: keep stages 1–5 unchanged (Q3 forward injection + Q4.1 mid-step snapshot for `J` sourcing), replace stages 6+7 (`average_fine_h_to_coarse` + `overwrite_coarse_e_from_fine`) with the single call to `inject_equivalent_currents_to_coarse`. CPML and clock-advance stages are unchanged. The Q4.1 helper now feeds the `J = +n̂ × H` term; the new outer-layer post-`update_fine_e` sample feeds the `M = −n̂ × E` term.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs`.
- **Base SHA dep:** Step B1 merged.
- **DoD:** `SubgriddedSolver::step` body matches spec §3 stage list 1-by-1 (verify by inspection of the doc-comment / source mapping); every existing `subgrid::*` unit test under `crates/yee-fdtd/tests/` continues to pass (no regression in non-`#[ignore]`'d coverage); plane-wave traversal **loose** gate (the existing one without `#[ignore]`) still passes within its existing tolerance.
- **Verification:** `cargo clippy -p yee-fdtd --all-targets -- -D warnings && cargo test -p yee-fdtd --release` exits 0. The Q5 strict gate remains `#[ignore]`'d at this step (un-ignored in B3).
- **Escape hatch:** blocked > 15 min on a regression in the existing loose Q5 gate → first verify Step B1's two unit tests still pass, then bisect by reverting `step_with_gaussian_source_ez` to the old closure while keeping `step` on Berenger; existing loose Q5 uses `step_with_gaussian_source_ez` (source-driven). If the two solver methods diverge in behaviour, the Q3-Dirichlet → `update_fine_e` → `inject_equivalent_currents_to_coarse` chain is wrong; re-read Berenger 2006 §III stage diagram.
- **LOC:** ~120.

### Step B3 — Q5 strict 0.5%-of-peak plane-wave traversal — un-`#[ignore]`

- **Brief:** Remove the `#[ignore]` attribute on `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs:111`. Update the test docstring to cite the Berenger 2006 closure and ADR-0035. The test geometry is unchanged: `(64, 32, 32)` coarse + `(16, 16, 16)`-coarse-cell fine, plane-wave traversal across the nest, 5 probes downstream, 500 steps, `max|E_z_sub − E_z_ref| / peak ≤ 0.5%`.
- **Lane:** `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs`.
- **Base SHA dep:** Step B2 merged.
- **DoD:** test passes within the 0.5% bound without `#[ignore]`; the strict-gate docstring cites Berenger 2006 + ADR-0035; CI `cargo test -p yee-fdtd --release` includes this test in the default (non-ignored) set.
- **Verification:** `cargo test -p yee-fdtd --release subgrid_plane_wave_traversal -- --include-ignored=false` exits 0 with the strict gate counted as a non-ignored pass.
- **Escape hatch:** blocked > 15 min with `rel_err` between 0.5% and 5% after B2 → re-check sign conventions in Step B1's per-face table; the most likely failure modes (per spec §6 risks) are (a) inverted `n̂` on `−x`/`−y`/`−z` faces, (b) `J_S` time-centered at `n + 1/4` instead of `n + 1/2` (Q4.1 snapshot timing), (c) double-counting at cuboid edges. If `rel_err > 5%`, the Berenger pipeline is wired wrong end-to-end and B2 needs re-review.
- **LOC:** ~20 (mostly comment + `#[ignore]` removal).

### Step B4 — Q6 round-trip energy-drift gate (10 000 steps)

- **Brief:** New test `crates/yee-fdtd/tests/subgrid_energy_balance.rs`. Initialise a Gaussian-modulated sinusoidal pulse inside the fine region of a closed-PEC `(64, 64, 64)` coarse + `(16, 16, 16)`-coarse-cell fine domain (no CPML — closed cavity so round-trip is well-defined). Propagate forward through the Huygens interface into the coarse region, reflect off PEC walls, propagate back. Integrate `W(t) = ∫ [ε₀|E|² + μ₀|H|²] dV` over coarse + fine at `t = 0` and `t = 10 000 · dt_coarse`. Assert `|W(N) − W(0)| / W(0) ≤ 0.5%`. Wall-time budget `< 5 min` `--release`; hardware-gate (`#[ignore]`) if it overruns the runner.
- **Lane:** `crates/yee-fdtd/tests/subgrid_energy_balance.rs` (create), `crates/yee-fdtd/validation/README.md`.
- **Base SHA dep:** Step B3 merged.
- **DoD:** test passes within the 0.5% bound at 10 000 steps; `crates/yee-fdtd/validation/README.md` has the new `fdtd-007 (stability, Berenger closure)` row.
- **Verification:** `cargo test -p yee-fdtd --release subgrid_energy_balance` exits 0.
- **Escape hatch:** blocked > 15 min with drift between 0.5% and 5% → record drift as `// regression-tracked: <value>`, mark `#[ignore]`, surface as a Phase 2.fdtd.7.x finding. The fallback path here would be a closer reading of Berenger 2006 §IV on the Huygens-surface energy balance — the spec did **not** derive this in full detail (see §3 escape note). Do not weaken the 0.5% gate without a spec amendment.
- **LOC:** ~220.

### Step B5 — Q7 fdtd-007 Maloney-Smith production gate (forward-port)

- **Brief:** Forward-port the spec `2026-05-18` Q7 brief unchanged. Geometry: slot `w = 0.5 mm` × `L = 30 mm` in infinite PEC ground plane, dielectric backing (`ε_r = 2.2`, `h = 1.524 mm`), delta-gap voltage drive. Coarse `dx = 1 mm`, fine `dx = 0.5 mm` over a `(40 × 6 × 4) mm` box. Extract `S_11(f)` via the existing `LumpedRlcPort` infrastructure; compare to Maloney & Smith 1993 Fig. 9.
- **Lane:** `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` (create), `crates/yee-fdtd/validation/README.md`.
- **Base SHA dep:** Step B4 merged. (Independent of B4's energy-balance pass/fail in principle, but serialised here to keep the validation rollup coherent.)
- **DoD:** `f_res` within ±2% of Maloney-Smith Fig. 9; `|S_11(f_res)|` within ±1 dB; internal sanity check (same problem on a globally uniform `dx = 0.5 mm` grid) agrees with the subgridded result within 0.3% / 0.3 dB at five spot frequencies. Wall-time `< 30 min` `--release`; hardware-gate behind `#[ignore]` if it overruns.
- **Verification:** `cargo test -p yee-validation --release fdtd_007_maloney_smith_slot` exits 0.
- **Escape hatch:** unchanged from the original 7.0 Q7 escape hatch — if Maloney-Smith Fig. 9 numeric values are not at hand, digitise to ±5% by eye, run with the loosened tolerance, mark `// TBD: tighten when Fig. 9 digitisation verified against the journal`. Do not invent reference numbers.
- **LOC:** ~350.

## Track-letter sequencing

Strict serial: B1 → B2 → B3 → B4 → B5. Each step is a closure-layer change with whole-pipeline knock-on effects (B1 adds the helper, B2 wires it into the step, B3 verifies the short-time gate, B4 verifies the long-time gate, B5 verifies the published-benchmark gate); they do not parallelise. Critical-path depth 5, well within CLAUDE.md §5's parallel-agent envelope.

The track-letter assignment is handed off to the implementation dispatcher (queued as track ZZZZZZ for the spec/plan landing; sub-tracks within ZZZZZZ pick letters).

## Validation rollup

| Gate | Step | Tolerance | Run-time | Status before 7.x |
|------|------|-----------|----------|-------------------|
| Q5 strict plane-wave traversal | B3 | `max\|E_z_sub − E_z_ref\| / peak ≤ 0.5%` over 500 steps, 5 probes | `< 30 s` `--release` | `#[ignore]`'d |
| Q6 round-trip energy-drift | B4 | `\|W(10⁴) − W(0)\| / W(0) ≤ 0.5%` | `< 5 min` `--release` | not implemented (blocked) |
| **fdtd-007 Maloney-Smith** | B5 | `f_res ±2%`, `\|S_11(f_res)\| ±1 dB`, 0.3% / 0.3 dB sanity vs uniform-fine | `< 30 min` `--release`, hardware-gate if overrun | not implemented (blocked) |

Per CLAUDE.md §4 ("No solver feature ships without a published-benchmark validation case") — fdtd-007 (Maloney-Smith) is the published benchmark; Q5 + Q6 are the internal regressions that catch the closure-instability failure mode VVVVVV diagnosed. All three rows land in `crates/yee-fdtd/validation/README.md`.

## Lane / file inventory

| Step | Files |
|------|-------|
| B1 | `crates/yee-fdtd/src/subgrid.rs` |
| B2 | `crates/yee-fdtd/src/subgrid.rs` |
| B3 | `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` |
| B4 | `crates/yee-fdtd/tests/subgrid_energy_balance.rs` (create), `crates/yee-fdtd/validation/README.md` |
| B5 | `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` (create), `crates/yee-fdtd/validation/README.md` |

`yee-core`, `yee-mesh`, `yee-cli`, `yee-py`, `yee-gui`, `yee-cuda` are untouched. CLI / Python / GUI exposure remains deferred to Phase 2.fdtd.7.0.1 follow-up.

## Risk register

1. **TF/SF sign-convention bookkeeping at the Huygens surface** (spec §6 risk 1). The most likely first-cut failure mode is an inverted normal on a `−x`/`−y`/`−z` face or a swapped sign on `J` vs `M`. **Surfaces in Step B3** as the Q5 gate failing at the `~1-10%` level (sign error gives a clearly bounded miscoupling, not a divergence). Mitigation: B1's per-face index table is a single 6-row source-of-truth.
2. **Cuboid-edge double-counting** (spec §6 risk 2). If two adjacent faces both include the shared edge in their tangential index range, the corresponding edge cell receives `2×` the equivalent current. **Surfaces in Step B1** as the `cuboid_edge_owned_by_one_face_only` unit test failing — caught before B2 wiring.
3. **Time-centering mismatch** (spec §6 risk 3). `J` sampled at `n + 1/2` (correct, from Q4.1 helper), `M` sampled at `n + 1` (correct, from post-`update_fine_e` outer fine `E_t`). Silent on Q5 short-time; **surfaces in Step B4** as a slow energy drift between 0.5% and 5%. Mitigation: B4's escape hatch is to mark `#[ignore]` with the recorded drift, not to invent a phase shift; that would be a spec amendment.

## Out of scope

Explicit non-goals for this plan, unchanged from spec `2026-05-18` §2:

- **No dispersive ADE inside the fine region** (Drude / Lorentz / Debye). Phase 2.fdtd.7.2.
- **No multi-region nesting or nest-inside-nest.** Phase 2.fdtd.7.1.
- **No co-location with CPML or TF/SF faces.** Spec `2026-05-18` §2 documents this as a runtime error in 7.0; the Berenger closure does not change that. Phase 2.fdtd.7.3.
- **No GPU.** CPU-only, single-threaded, scalar FP64 throughout.
- **No refinement ratios other than 2×.** Phase 2.fdtd.7.4.
- **No higher-order spatial interpolation.** Linear in space and time only; cubic deferred to Phase 2.fdtd.7.4.
- **No CLI / Python / GUI exposure.** Direct Rust API only in 7.x. Phase 2.fdtd.7.0.1 follow-up.
- **No removal of the old `average_fine_h_to_coarse` / `overwrite_coarse_e_from_fine` helpers.** Per ADR-0035, retained `#[doc(hidden)]` for posterity. Removal is a separate spec amendment.

## Final verification

```bash
cargo build  -p yee-fdtd -p yee-validation
cargo clippy -p yee-fdtd -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-fdtd --release
cargo test   -p yee-validation --release
cargo fmt    --check --all
cargo doc    --no-deps -p yee-fdtd
mdbook build docs/
```

All seven must exit 0. Every existing `crates/yee-fdtd/tests/` non-`#[ignore]`'d regression stays green; the Q5 strict gate is now non-`#[ignore]`'d and passing; the new Q6 energy-balance gate passes; the new Q7 Maloney-Smith production gate passes (or is `#[ignore]`'d for wall-time only, per the escape hatch).

## Estimated total

- LOC: ~990 (B1 ~280, B2 ~120, B3 ~20, B4 ~220, B5 ~350).
- Wall-time per agent: 3–5 days end-to-end. Critical path B1→B3 is ~2 days; B4 + B5 add 1–2 days serial.
- Risk concentration: Step B4 (no in-spec mitigation beyond "Berenger 2006 §IV says it works, verify by test"). Escape hatch preserves merge throughput at the cost of deferring the energy-balance gate to a hypothetical Phase 2.fdtd.7.y.
