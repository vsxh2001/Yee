# Phase 2.fdtd.7 — FDTD subgridding (walking skeleton) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan track-by-track.

**Companion spec:** `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`
**Base SHA:** `6f75b5a` (Track OOOOO merge).
**Target phase:** 2.fdtd.7.0 only. 7.1–7.5 are explicitly deferred — see §9.
**Tech-stack additions:** none. Reuses `ndarray`, `update::update_h` / `update::update_e`, and the existing `grid_and_cpml_mut` split-borrow primitive at `crates/yee-fdtd/src/lib.rs:158`.

---

## Goal

Phase 2.fdtd.7.0 ships one axis-aligned cuboidal 2× fine sub-region nested inside a uniform coarse `YeeGrid`, time-subcycled 2:1, with linear spatial / temporal coarse → fine `E_t` interpolation and area-averaged fine → coarse `H_t` closure (Chevalier 1997 §III–IV). Two validation gates: a 10 000-step round-trip energy-drift regression (≤ 0.5%) and the external Maloney-Smith 1993 dielectric-loaded thin-slot reference (resonance ±2%, `|S_11|` ±1 dB). CPU-only, single-threaded, scalar FP64, non-dispersive isotropic materials in the fine region, no co-location with CPML / TF-SF — same execution model as today's `WalkingSkeletonSolver`.

## Pre-flight refactor — surfaced from Track MMMMM finding #1

Spec §7 calls the per-stage `WalkingSkeletonSolver::step` split a "quality-of-life improvement, not blocking." Track MMMMM's design review disagreed: every fine sub-step needs `update_h_only` / `update_e_only` *without* the CPML / PEC / clock-advance side effects that the existing `step` / `step_with_source` / `step_with_plane_wave` bundle in. Doing the refactor in-flight conflates "split the step" with "introduce subgridding" and makes the diff hostile to review. Step 1 lands the refactor on its own, proving behaviour preservation against every existing `crates/yee-fdtd/tests/` regression byte-for-byte. Every subsequent step builds on the post-refactor surface.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-fdtd/src/lib.rs` | Modify | Split `step` / `step_with_source` / `step_with_plane_wave` into helpers (`update_h_only`, `apply_cpml_h`, `update_e_only`, `apply_cpml_e`, `apply_gaussian_source_ez`); declare `pub mod subgrid`. |
| `crates/yee-fdtd/src/subgrid.rs` | Create | `SubgridRegion`, `SubgriddedSolver`, interface helpers, seven-stage `step`. |
| `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` | Create | Step 5 traversal integration test. |
| `crates/yee-fdtd/tests/subgrid_round_trip_energy.rs` | Create | Step 6 stability/reciprocity gate. |
| `crates/yee-fdtd/validation/README.md` | Modify | New `fdtd-007 (stability)` and `fdtd-007 (Maloney-Smith)` rows. |
| `crates/yee-validation/src/lib.rs` | Modify | Step 7 — `run_fdtd_007_maloney_smith_slot` driver. |
| `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` | Create | Step 7 production-gate test. |

No changes to `yee-core`, `yee-mesh`, `yee-cli`, `yee-py`, `yee-gui`, `yee-cuda`. The `yee-fdtd` `#![forbid(unsafe_code)]` floor is preserved throughout.

## Step ladder

### Step 1 (Track Q1) — `WalkingSkeletonSolver::step` refactor into composable helpers

- **Brief:** Mechanical extraction of `step`, `step_with_source`, `step_with_plane_wave`, `step_with_source_and_ntff` into five `pub fn` helpers: `update_h_only`, `apply_cpml_h`, `update_e_only`, `apply_cpml_e`, `apply_gaussian_source_ez`. Re-express the existing step methods as straight helper sequences plus `advance_clock`. No reordering, no "tidy-up."
- **Lane:** `crates/yee-fdtd/src/lib.rs`.
- **Base SHA dep:** none — branches off `6f75b5a` directly.
- **DoD:** five new `pub fn` helpers, documented under `#![warn(missing_docs)]`; all four existing `step*` methods re-expressed on top; every existing test under `crates/yee-fdtd/tests/` (9 files: cpml_reflection, dipole_pattern, dispersive, fdtd_propagation, lumped_resistor, ntff_dipole, plane_wave_finite_box, plane_wave_oblique, plane_wave_propagation) passes with byte-identical numerics.
- **Verification:** `cargo clippy -p yee-fdtd --all-targets -- -D warnings && cargo test -p yee-fdtd --release` exits 0.
- **Escape hatch:** blocked > 15 min on the byte-identical regression → surface and stop. Likely cause is reordering inside `step_with_source` (source injection runs *between* H-update + CPML-H and E-update; see lib.rs lines 199–214). Do not refactor beyond mechanical extraction.
- **LOC:** ~80.

### Step 2 (Track Q2) — `SubgridRegion` struct + axis-aligned 2× sub-Yee-grid scaffold

- **Brief:** Create `crates/yee-fdtd/src/subgrid.rs`. `SubgridRegion::new(parent: &YeeGrid, lo, hi)` constructs a fine `YeeGrid` with `dx_fine = parent.dx / 2`, `dy_fine = parent.dy / 2`, `dz_fine = parent.dz / 2`, `dt_fine = parent.dt / 2`, sized `(2·(hi.0-lo.0), 2·(hi.1-lo.1), 2·(hi.2-lo.2))` cells, inheriting scalar `eps_r` / `mu_r`. No coupling yet. Stub `SubgriddedSolver::step` as `inner.step()` (fine grid dormant) to unlock downstream parallel work against the type signatures. Constructor errors out if `lo`/`hi` overlap CPML thickness or any TF/SF box face (spec §6 — co-location with CPML/TF-SF is a documented 7.0 runtime error).
- **Lane:** `crates/yee-fdtd/src/lib.rs` (one `pub mod subgrid;` + re-export), `crates/yee-fdtd/src/subgrid.rs`.
- **Base SHA dep:** Step 1 merged (uses the new helper surface).
- **DoD:** `SubgridRegion::new` returns a fine grid with halved `dx`, `dy`, `dz`, `dt` (assert within `f64::EPSILON`); `fine_grid()`/`_mut()` getters present; `SubgriddedSolver::step` placeholder produces fields identical to `WalkingSkeletonSolver::step` baseline at step 10; `cargo doc -p yee-fdtd --no-deps` is `missing_docs`-clean.
- **Verification:** `cargo clippy -p yee-fdtd --all-targets -- -D warnings && cargo test -p yee-fdtd --release subgrid && cargo doc -p yee-fdtd --no-deps` exits 0.
- **Escape hatch:** blocked > 15 min on `YeeGrid` construction (material-array shape mismatch) → read `crates/yee-fdtd/src/grid.rs` first; do not refactor `YeeGrid::new`.
- **LOC:** ~180.

### Step 3 (Track Q3) — coarse → fine `E_t` spatial + temporal interpolation

- **Brief:** Implement spec §3 coarse → fine `E`-field driver on `SubgridRegion`. `snapshot_coarse_e_t(parent: &YeeGrid)` caches start- and end-of-coarse-step parent `E_t` on the six interface faces. `interpolate_coarse_e_to_fine(frac: f64)` writes Dirichlet fine boundary `E_t`: linear spatial interpolation between bracketing coarse edges (Chevalier 1997 §III) blended in time at `frac ∈ {0.25, 0.75}` per fine sub-step. Six face orientations × (`E_y`, `E_z` on `±x`; `E_x`, `E_z` on `±y`; `E_x`, `E_y` on `±z`). Not yet wired into `step` — that is Step 5.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs`.
- **Base SHA dep:** Step 2 merged.
- **DoD:** unit tests — `uniform_field_maps_to_uniform_fine_e_t` (constant parent → constant fine within `f64::EPSILON · 10`); `linear_gradient_preserved` (parent `E_y(x) = a + b·x` → exact at fine-edge midpoints within `f64::EPSILON · 100`); `temporal_blend_at_frac_one_quarter` (start=0, end=1 → fine=0.25). All six face orientations exercised.
- **Verification:** `cargo test -p yee-fdtd --release subgrid::interp` exits 0; `cargo clippy --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the coarse↔fine index map for a face → write out the index table for `+x` only, get its three unit tests green, copy-paste to the other five. Do not attempt the general 6-face dispatch on first pass.
- **LOC:** ~300.

### Step 4 (Track Q4) — fine → coarse `H_t` area-averaging closure

- **Brief:** Implement spec §3 fine → coarse closure. `SubgridRegion::average_fine_h_to_coarse(parent: &mut YeeGrid)` overwrites coarse `H_t` on the interface with the area-weighted average of the four fine `H_t` cells covering each coarse face (Chevalier 1997 §IV — the step that closes the discrete energy balance). Symmetric `overwrite_coarse_e_from_fine(parent: &mut YeeGrid)` for stage 7 of the time-step pattern (edge-average of two coincident fine `E_t` edges per coarse edge).
- **Lane:** `crates/yee-fdtd/src/subgrid.rs`.
- **Base SHA dep:** Step 3 merged (serialised after Q3 to keep `subgrid.rs` diffs reviewable; dependency-wise these are siblings of Q3).
- **DoD:** unit tests — `forward_reverse_round_trip_preserves_static_field` (snapshot → interpolate(0.5) → no fine update → overwrite ⇒ coarse matches original within `f64::EPSILON · 100`); `area_average_of_uniform_h_is_uniform` (constant fine `H_z` on a face → exact constant coarse `H_z`).
- **Verification:** `cargo test -p yee-fdtd --release subgrid::average` exits 0; `cargo clippy --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the 4-fine-cell coverage map → hard-code the four-index sum for the `+x` face first, generalise after.
- **LOC:** ~220.

### Step 5 (Track Q5) — time-subcycling loop in `SubgriddedSolver::step`

- **Brief:** Replace the Step-2 placeholder with the full seven-stage spec §3 sequence: `inner.update_h_only → inner.apply_cpml_h → snapshot_coarse_e_t → fine k=1 (interpolate(0.25), update_h_fine, update_e_fine) → inner.update_e_only → inner.apply_cpml_e → snapshot_coarse_e_t → fine k=2 (interpolate(0.75), update_h_fine, update_e_fine) → average_fine_h_to_coarse → overwrite_coarse_e_from_fine → inner.advance_clock`. Companion `step_with_gaussian_source_ez` injects on the coarse grid only (fine carries no source in 7.0). The Step-1 helpers earn their keep here — body is a direct transcription with no CPML reimplementation.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs`, `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` (create).
- **Base SHA dep:** Steps 2, 3, 4 all merged.
- **DoD:** integration test `subgrid_plane_wave_traversal` — drive a Gaussian-modulated plane wave across a vacuum fine region, compare to a uniform-fine-grid reference run (same `dx = dx_fine` everywhere, no nest). At 5 probe points downstream of the nest, time-domain `E_z` traces agree within **0.5%** of peak amplitude over the first 500 steps.
- **Verification:** `cargo test -p yee-fdtd --release subgrid_plane_wave_traversal` exits 0; `cargo clippy --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on > 0.5% error → first re-run the Step-2 passthrough test, then disable area-averaging (stage 7) to isolate which coupling direction is mis-wired. Pure forward-only coupling should still give < 5%; > 5% means spatial interpolation or temporal blend fraction is wrong.
- **LOC:** ~250.

### Step 6 (Track Q6) — stability / reciprocity gate (10 000-step round-trip energy drift)

- **Brief:** Implement spec §5 round-trip energy test. Initialise a Gaussian-modulated sinusoid in the fine region, propagate forward through the interface into the coarse region, reflect off PEC walls on all six outer coarse faces (no CPML — closed cavity so round-trip is well-defined), propagate back through the interface. Integrate `W(t) = ∫ [ε₀|E|² + µ₀|H|²] dV` (sum coarse + fine, using `dx³_coarse` / `dx³_fine` per region) at `t = 0` and `t = 10 000 · dt_coarse`. Assert `|W(N) − W(0)| / W(0) ≤ 0.5%`. Geometry: `(64,64,64)` coarse + `(16,16,16)`-coarse-cell fine region centred, yielding `(32,32,32)` fine cells. Wall-time `< 5 min` `--release`; if overrun, hardware-gate behind `#[ignore]` per the mom-001 / Phase 1.5 precedent (CLAUDE.md §4).
- **Lane:** `crates/yee-fdtd/tests/subgrid_round_trip_energy.rs` (create), `crates/yee-fdtd/validation/README.md`.
- **Base SHA dep:** Step 5 merged.
- **DoD:** test passes within the 0.5% bound at 10 000 steps; validation README has the `fdtd-007 (stability)` row.
- **Verification:** `cargo test -p yee-fdtd --release subgrid_round_trip_energy` exits 0.
- **Escape hatch:** blocked > 15 min with drift between 0.5% and 5% → record drift as `// regression-tracked: <value>`, mark `#[ignore]`, surface as Phase 2.fdtd.7.x finding. Do **not** weaken the 0.5% gate without a spec amendment — Berenger 2003 §IV says higher drift indicates asymmetric-coupling failure and that is a spec-level decision.
- **LOC:** ~200.

### Step 7 (Track Q7) — fdtd-007 production gate (Maloney-Smith dielectric-loaded thin slot)

- **Brief:** Implement spec §5 external validation. Geometry: slot `w = 0.5 mm` × `L = 30 mm` in infinite PEC ground plane, dielectric backing (`ε_r = 2.2`, `h = 1.524 mm`), delta-gap voltage drive. Coarse `dx = 1 mm`; fine `dx = 0.5 mm` over a `(40 × 6 × 4) mm` box centred on slot + substrate. Extract `S_11(f)` via the existing `LumpedRlcPort` infrastructure (`crates/yee-fdtd/src/lumped.rs`). Compare to Maloney & Smith 1993 Fig. 9.
- **Lane:** `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` (create), `crates/yee-fdtd/validation/README.md`.
- **Base SHA dep:** Step 5 merged. Independent of Step 6 — both gates can run in parallel post-Q5.
- **DoD:** `f_res` within ±2% of Maloney-Smith Fig. 9; `|S_11(f_res)|` within ±1 dB; internal sanity check (same problem on a globally uniform `dx = 0.5 mm` grid) agrees with the subgridded result within 0.3% / 0.3 dB at five spot frequencies. Wall-time `< 30 min` `--release`; hardware-gate behind `#[ignore]` if overrun.
- **Verification:** `cargo test -p yee-validation --release fdtd_007_maloney_smith_slot` exits 0.
- **Escape hatch:** blocked > 15 min on the Maloney-Smith Fig. 9 numeric values → digitise the curve to ±5% by eye, run with the loosened tolerance, mark `// TBD: tighten when Fig. 9 digitisation verified against the journal`, surface as a finding. Do not invent reference numbers.
- **LOC:** ~350.

## Track-letter sequencing

Q1 must land first (refactor base). Q2 must land before Q3, Q4, Q5 (defines type signatures). Q3 and Q4 are dependency-siblings on Q2; per CLAUDE.md §5 we serialise (Q3 → Q4) to keep individual `subgrid.rs` diffs reviewable. Q5 depends on Q1, Q2, Q3, Q4. Q6 and Q7 are independent gates that run in parallel once Q5 lands.

Critical path: `Q1 → Q2 → Q3 → Q4 → Q5` (5 sequential merges) + `Q6 ‖ Q7` (one parallel pair). Within CLAUDE.md §5's "up to 5 parallel agents" envelope.

## Validation rollup

| Gate | Step | Tolerance | Run-time |
|------|------|-----------|----------|
| Stability / reciprocity — 10 000-step round-trip energy drift | Step 6 (Q6) | `|W(10⁴) − W(0)| / W(0) ≤ 0.5%` | `< 5 min` `--release` |
| **fdtd-007 Maloney-Smith** — dielectric-loaded thin slot `S_11` | Step 7 (Q7) | `f_res` ±2%, `|S_11(f_res)|` ±1 dB, plus 0.3% / 0.3 dB sanity check vs uniform-fine | `< 30 min` `--release`, hardware-gate if overrun |

Both rows land in `crates/yee-fdtd/validation/README.md`. Per CLAUDE.md §4 "No solver feature ships without a published-benchmark validation case" — fdtd-007 is the published benchmark; the stability test is the internal regression catching late-time interface instability (Berenger 2003 §IV).

## Lane / file inventory

| Step | Files |
|------|-------|
| Q1 | `crates/yee-fdtd/src/lib.rs` |
| Q2 | `crates/yee-fdtd/src/lib.rs` (one `pub mod subgrid;`), `crates/yee-fdtd/src/subgrid.rs` (create) |
| Q3 | `crates/yee-fdtd/src/subgrid.rs` |
| Q4 | `crates/yee-fdtd/src/subgrid.rs` |
| Q5 | `crates/yee-fdtd/src/subgrid.rs`, `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` (create) |
| Q6 | `crates/yee-fdtd/tests/subgrid_round_trip_energy.rs` (create), `crates/yee-fdtd/validation/README.md` |
| Q7 | `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` (create), `crates/yee-fdtd/validation/README.md` |

Cross-lane consumers (`yee-cli`, `yee-py`, `yee-gui`) are not touched in 7.0. CLI / Python / GUI exposure lands as follow-up Phase 2.fdtd.7.0.1 once both gates are green.

## Risk register

1. **Late-time interface instability** (spec §6, Berenger 2003 §IV). Classical subgridding failure mode; Chevalier area-average is supposed to close the energy balance to second order in `dx_coarse` but asymmetric coupling can still grow over `O(10⁴)` steps. **Surfaces in Step 6 (Q6)** — the 10 000-step gate is the canary. Fallback is Berenger 2006's Huygens-surface scheme, out of scope here.
2. **CPML / TF-SF co-location with the fine region** (spec §6). Out of scope; documented as runtime error. **Materialises in Step 2 (Q2)** as a constructor-time check: `SubgridRegion::new` returns an `Error::Invalid` if `lo`/`hi` overlap the CPML thickness or any TF/SF box face.
3. **Fine CFL forcing the sub-cycle** (spec §6). The 2× time-subcycle *is* the mitigation; at 2× the fine update cost is roughly the coarse update cost. **Materialises in Step 5 (Q5)** — if the traversal test runs > 10× slower than uniform-coarse, the per-stage cost balance is off (likely a buffer reallocation in `snapshot_coarse_e_t`).
4. **NTFF surface intersecting the fine region** (spec §6). NTFF DFT bins go inconsistent if the Huygens surface crosses the nest. **Surfaces in Step 5** as a runtime error in `SubgriddedSolver::new` and a doc-comment limitation. Phase 2.fdtd.7.3 lifts this.
5. **Maloney-Smith reference values not digitised**. The reference is Fig. 9 of a 1993 IEEE T-AP paper. **Materialises in Step 7 (Q7)** as the explicit escape hatch — digitise to ±5%, mark TBD, do not invent numbers.

## Out of scope

Explicit non-goals for this plan, per spec §2:

- **No dispersive ADE inside the fine region** (Drude / Lorentz / Debye). Phase 2.fdtd.7.2.
- **No multi-region nesting or nest-inside-nest.** Phase 2.fdtd.7.1.
- **No co-location with CPML or TF/SF faces.** Spec §2 documents this as a runtime error in 7.0. Phase 2.fdtd.7.3.
- **No GPU.** CPU-only, single-threaded, scalar FP64 throughout; same execution model as `WalkingSkeletonSolver` today.
- **No refinement ratios other than 2×.** Phase 2.fdtd.7.4.
- **No higher-order spatial interpolation.** Linear in space and time only. Phase 2.fdtd.7.4.
- **No CLI / Python / GUI exposure.** Direct Rust API only in 7.0. Phase 2.fdtd.7.0.1 follow-up.
- **No conformal (Dey-Mittra) interaction.** Staircase geometry workspace-wide.

## Final verification

```bash
cargo build  -p yee-fdtd -p yee-validation
cargo clippy -p yee-fdtd -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-fdtd --release
cargo test   -p yee-validation --release
cargo fmt    --check --all
cargo doc    --no-deps -p yee-fdtd
```

All six must exit 0. Every existing `crates/yee-fdtd/tests/` regression stays green (Step 1's byte-identical floor); the new subgrid tests, the stability gate, and the Maloney-Smith production gate all pass.

## Estimated total

- LOC: ~1 580 (Q1 ~80, Q2 ~180, Q3 ~300, Q4 ~220, Q5 ~250, Q6 ~200, Q7 ~350).
- Wall-time per agent: 4–6 days end-to-end at one-engineer pace. Critical path Q1→Q5 is ~3 days; Q6 + Q7 add 1–2 days in parallel.
- Risk concentration: Step 6's stability gate (no in-spec mitigation beyond "Chevalier area-average, verify by test"). Escape hatch preserves merge throughput at the cost of deferring the gate to Phase 2.fdtd.7.x.
