# Phase 2.fdtd.7.y — Berenger M-side compensating-source amendment — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans`
> to drive this plan step-by-step.

**Companion spec:** `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md`
**Base SHA:** `2e54e6f` (`main` HEAD after Track OOOOOOO `464c7ba` merged: Phase
2.fdtd.7.x B1 + B2 + B2.1 + B2.2 wired; J-side coarse-ghost subtraction live; M
side B2.1 form `M = -n̂ × E_fine_end_of_step` (no ghost subtraction); 500-step
Berenger canary + Q5 strict + Q6 round-trip all `#[ignore]`'d.).
**Target phase:** 2.fdtd.7.y only — the M-side closure amendment. 7.1–7.5
stay deferred per Phase 2.fdtd.7.x §9.
**Tech-stack additions:** none.

---

## Goal

Amend the Phase 2.fdtd.7.x M-side closure to use **Option β** (compensating
source `M = -n̂ × (E_post − E_pre)`, where `E_pre` is the fine outer-layer
`E_t` immediately after the Q3 Dirichlet interpolation writes it and `E_post`
is the same `E_t` after sub-step 2's `update_fine_e` completes). Keep the
B2.2 J-side coarse-ghost subtraction unchanged. Un-`#[ignore]` the 500-step
Berenger canary, the strict Q5 0.5% gate, and the Q6 10 000-step round-trip
energy-drift gate.

## Pre-flight

No new dependencies. The B2.1 `snapshot_fine_e_end_of_step` helper provides
the array-clone pattern; the new `snapshot_fine_e_pre_update` and
`snapshot_fine_e_post_update` follow it.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-fdtd/src/subgrid.rs` | Modify | Add `snapshot_fine_e_pre_update` / `snapshot_fine_e_post_update`. Switch `inject_m_to_coarse_h` to the compensating-source form. Plumb the two new snapshot calls into `SubgriddedSolver::step` / `step_with_gaussian_source_ez`. Keep the legacy `snapshot_fine_e_end_of_step` for one release as `#[doc(hidden)]` (consumed by no caller after C1). |
| `crates/yee-fdtd/tests/berenger_traversal.rs` | Modify | Remove the `#[ignore]` on `berenger_step_propagates_without_divergence_500_steps`; update its docstring to cite Phase 2.fdtd.7.y Option β. |
| `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` | Modify | Remove the `#[ignore]` on `strict_05pct_peak_over_500_steps` (Phase 2.fdtd.7.x B3 deferred work); update docstring to cite Phase 2.fdtd.7.y Option β. |
| `crates/yee-fdtd/tests/subgrid_energy_balance.rs` | Modify or create | Phase 2.fdtd.7.x B4 deferred work — un-`#[ignore]` if it already exists, else create. |
| `crates/yee-fdtd/validation/README.md` | Modify | Bump the Berenger stability rollup from "B2.2 partial, M-side deferred" to "B2.2 + Phase 2.fdtd.7.y Option β, gates passing". |

No changes to `yee-core`, `yee-mesh`, `yee-cli`, `yee-py`, `yee-gui`, `yee-cuda`,
`yee-mom`, `yee-io`, `yee-plotters`, `yee-validation`. Lane is
`crates/yee-fdtd/**`.

## Step ladder

### Step C1 — `snapshot_fine_e_pre_update` + `snapshot_fine_e_post_update` + step plumbing

- **Brief:** Add two public methods on `SubgridRegion` that clone the fine
  `E_t` arrays (`ex`, `ey`, `ez`) into named cache slots
  (`fine_e_pre_snapshot` and `fine_e_post_snapshot`), mirroring the existing
  `snapshot_fine_e_end_of_step` helper's pattern. Rustdoc each method with a
  one-line spec reference and a one-paragraph time-level diagram (E_pre =
  after Q3 Dirichlet, before update_fine_e; E_post = after update_fine_e).
  Plumb the two calls into `SubgriddedSolver::step` and
  `step_with_gaussian_source_ez` per spec §3 sub-step 2 ordering:

  ```text
  fine sub-step 2:
    region.interpolate_coarse_e_to_fine(0.75)
    region.snapshot_fine_e_pre_update            ← NEW
    region.update_fine_e
    region.snapshot_fine_e_post_update           ← NEW
    region.update_fine_h                         (unchanged)
  ```

  Do **not** yet wire the new caches into `inject_m_to_coarse_h` — that is
  Step C2. The M source remains the B2.1 form for C1; the new caches are
  populated-but-unread.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs` only.
- **Base SHA dep:** `2e54e6f`.
- **DoD:** both methods present and `pub`; rustdoc cites
  `2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3 by relative path; a
  new unit test
  `fine_e_pre_and_post_snapshots_differ_after_fine_update`
  exercises the snapshot pair on a non-trivial fine-`E` state and asserts the
  two cached arrays are not bit-identical (sanity check that
  `update_fine_e` actually moves the fine E between the two snapshot
  points). `cargo clippy -p yee-fdtd --all-targets -- -D warnings` exits 0.
- **Verification:** `cargo test -p yee-fdtd --release subgrid::fine_e_pre`
  exits 0; existing `berenger_step_propagates_without_divergence` 100-step
  canary still passes (no behaviour change yet).
- **Escape hatch:** blocked > 15 min on cache-slot lifetime / borrow-checker
  issues → mirror the existing `FineESnapshot` struct verbatim and add a
  second field; do not invent a generic snapshot abstraction.
- **LOC:** ~80.

### Step C2 — `inject_m_to_coarse_h` compensating-source switch

- **Brief:** Replace `inject_m_to_coarse_h`'s M source from
  `M = -n̂ × E_fine_end_of_step` to `M = -n̂ × (E_post − E_pre)` per spec
  Option β. The implementation reads the two new snapshot caches and feeds
  the difference into the existing `inject_*_face` per-face helpers. Crucially:
  the new `use_compensating_source = true` flag must also force
  `do_ghost = false` on the M side (no double-differencing). Add a
  `debug_assert!(!(do_ghost && use_compensating_source))` guard in the
  per-face helper. Update the docstring on `inject_m_to_coarse_h` to cite
  Phase 2.fdtd.7.y Option β and the OOOOOOO empirical regression that
  motivated it.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs` only.
- **Base SHA dep:** Step C1 merged.
- **DoD:** `inject_m_to_coarse_h` body computes `E_post − E_pre` from the
  C1 caches before calling `inject_*_face`; the J-side
  `inject_j_to_coarse_e` is **byte-for-byte unchanged**; the existing
  100-step `berenger_step_propagates_without_divergence` canary still
  passes with peak `|E_z|_fine ≤ 2.75` V/m (the B2.2 baseline) — i.e. C2
  does not regress short-time behaviour. Verified by running both
  100-step and 200-step variants of the canary.
- **Verification:** `cargo test -p yee-fdtd --release berenger_traversal::berenger_step_propagates_without_divergence`
  exits 0 with peak `|E_z|_fine` reported `< 3.0` V/m.
- **Escape hatch:** if the 100-step canary regresses past 3.0 V/m
  (compensating source too noisy / wrong sign), revert the
  `inject_m_to_coarse_h` body to the B2.1 form, keep the C1 plumbing
  for diagnostic logging, and surface as a Step C5 trigger.
- **LOC:** ~100.

### Step C3 — Un-`#[ignore]` the 500-step Berenger canary and the strict Q5 gate

- **Brief:** Remove the `#[ignore]` attributes on
  `berenger_traversal::berenger_step_propagates_without_divergence_500_steps`
  and `subgrid_plane_wave_traversal::strict_05pct_peak_over_500_steps`.
  Update both docstrings to cite Phase 2.fdtd.7.y Option β and the spec /
  ADR locations. Verify each test passes under `cargo test --release` with
  the new closure live; if **either** still fails, that is a Step C5
  trigger.
- **Lane:** `crates/yee-fdtd/tests/berenger_traversal.rs`,
  `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs`.
- **Base SHA dep:** Step C2 merged.
- **DoD:** both tests pass; the 500-step canary reports peak `|E_z|_fine`
  bounded below `1e3` V/m at step 499; the strict Q5 gate reports
  `rel_err < 0.5%` of peak across all 5 downstream probes; CI default-set
  `cargo test -p yee-fdtd --release` counts both as passing.
- **Verification:** `cargo test -p yee-fdtd --release berenger_step_propagates_without_divergence_500_steps subgrid_plane_wave_traversal -- --include-ignored=false`
  exits 0.
- **Escape hatch:** see spec §6 risk 1 / risk 2 — if the 500-step canary
  passes but the strict Q5 gate sits at `0.5% < rel_err < 5%`, the
  compensating source is the right sign but undersized — pivot to
  Step C5 (Option α). If both fail simultaneously, Option β has
  degenerated to B2.2; Step C5.
- **LOC:** ~20.

### Step C4 — Q6 round-trip energy-drift gate sweep

- **Brief:** Phase 2.fdtd.7.x B4 forward-port. If
  `crates/yee-fdtd/tests/subgrid_energy_balance.rs` already exists (Phase
  2.fdtd.7.x B4 deferred), un-`#[ignore]` it and update its docstring; else
  create it per the Phase 2.fdtd.7.x B4 brief verbatim. Run the 10 000-step
  closed-PEC `(64, 64, 64)` coarse + `(16, 16, 16)`-coarse-cell fine cavity,
  assert `|W(10 000) − W(0)| / W(0) ≤ 0.5%`. Bump
  `crates/yee-fdtd/validation/README.md` Berenger stability rollup
  accordingly.
- **Lane:** `crates/yee-fdtd/tests/subgrid_energy_balance.rs`,
  `crates/yee-fdtd/validation/README.md`.
- **Base SHA dep:** Step C3 merged.
- **DoD:** Q6 gate passes within 0.5% bound at 10 000 steps; validation
  README rollup table reflects the Phase 2.fdtd.7.y status.
- **Verification:** `cargo test -p yee-fdtd --release subgrid_energy_balance`
  exits 0; wall-time < 5 min `--release`. If it overruns, hardware-gate
  per Phase 2.fdtd.7.x B4 escape hatch.
- **Escape hatch:** if drift between `0.5%` and `5%`, record the drift in a
  `// regression-tracked: <value>` comment and `#[ignore]` — Step C5
  trigger.
- **LOC:** ~30 (assuming the B4 test body already exists; else inherit the
  Phase 2.fdtd.7.x B4 ~220 LOC budget).

### Step C5 (optional escape hatch) — Option α absorbing-BC pivot

- **Brief:** Conditional on Step C3 or C4 failing. Replace the Q3
  coarse→fine `E_t` Dirichlet interpolation with a second-order Mur
  absorbing BC on the fine grid's outer `E_t` layer. Restore Berenger's
  canonical M source `M = -n̂ × (E_TF_fine − E_SF_coarse_ghost)` with
  coarse-ghost subtraction. Verify against the same Q5 / Q6 / canary
  gates.
- **Lane:** `crates/yee-fdtd/src/subgrid.rs`, `crates/yee-fdtd/tests/**`.
- **Base SHA dep:** Step C4 attempted-and-failed.
- **DoD:** spec amendment landed (new Phase 2.fdtd.7.y.α spec); all three
  gates pass with the absorbing-BC accuracy ceiling. Note that the
  strict Q5 gate's `0.5%` tolerance may need re-spec to `1%` per the
  Mur reflection floor analysis in spec §3 Option α trade-off; if so,
  amend the spec before merging C5.
- **Verification:** same as C3 + C4.
- **Escape hatch:** if Option α also fails to clear `0.5%`, stop and
  surface as a Phase 2.fdtd.8 (higher-order sub-gridding) trigger.
- **LOC:** ~250 (Mur BC + ghost-restoration + spec/ADR delta).

## Track-letter sequencing

Strict serial: C1 → C2 → C3 → C4, with C5 as a conditional fallback off
C3 or C4. Each step is a closure-layer change with whole-pipeline
knock-on effects (C1 adds plumbing, C2 wires the new source, C3 verifies
short-time, C4 verifies long-time, C5 is the escape pivot). Critical-path
depth 4 (5 with C5).

The track-letter assignment is handed off to the implementation
dispatcher.

## Validation rollup

| Gate | Step | Tolerance | Run-time | Status before 7.y |
|------|------|-----------|----------|-------------------|
| Berenger 500-step canary | C3 | `peak \|E_z\|_fine < 1e3` V/m | `< 30 s` `--release` | `#[ignore]`'d (peak ≈ 1.035e3 at step ~137) |
| Q5 strict plane-wave traversal | C3 | `max\|E_z_sub − E_z_ref\| / peak ≤ 0.5%` over 500 steps, 5 probes | `< 30 s` `--release` | `#[ignore]`'d |
| Q6 round-trip energy-drift | C4 | `\|W(10⁴) − W(0)\| / W(0) ≤ 0.5%` | `< 5 min` `--release` | `#[ignore]`'d (or not implemented if B4 deferred) |

Phase 2.fdtd.7.x's published-benchmark gate (fdtd-007 Maloney-Smith) is
unchanged; Phase 2.fdtd.7.y is an amendment to the closure mathematics
only.

## Lane / file inventory

| Step | Files |
|------|-------|
| C1 | `crates/yee-fdtd/src/subgrid.rs` |
| C2 | `crates/yee-fdtd/src/subgrid.rs` |
| C3 | `crates/yee-fdtd/tests/berenger_traversal.rs`, `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` |
| C4 | `crates/yee-fdtd/tests/subgrid_energy_balance.rs`, `crates/yee-fdtd/validation/README.md` |
| C5 | `crates/yee-fdtd/src/subgrid.rs`, `crates/yee-fdtd/tests/**`, plus Phase 2.fdtd.7.y.α spec amendment under `docs/superpowers/specs/**` and ADR-0039 under `docs/src/decisions/**` |

## Risk register

1. **Round-off floor swallows the compensating source.** Surfaces in
   Step C3 with strict Q5 rel_err in the `0.5%–5%` band. Mitigation:
   Step C5 (Option α).
2. **Compensating source degenerates to zero.** Surfaces in C3 with
   the 500-step canary still diverging at step ~137 (peak `|E_z|`
   matches B2.2 baseline). Mitigation: C5.
3. **Inadvertent ghost double-counting.** Surfaces as immediate
   100-step canary regression in C2; the `debug_assert!` guard
   catches it in debug builds. Mitigation: spec §6 risk 3 — docstring
   forbids ghost + compensating together; runtime guard backstops.

## Out of scope

Same non-goals as Phase 2.fdtd.7.x:

- No dispersive ADE in the fine region (Phase 2.fdtd.7.2).
- No multi-region nesting (Phase 2.fdtd.7.1).
- No co-location with CPML / TF/SF (Phase 2.fdtd.7.3).
- No GPU.
- No refinement ratios other than 2× (Phase 2.fdtd.7.4).
- No higher-order spatial interpolation (Phase 2.fdtd.7.4).
- No CLI / Python / GUI exposure (Phase 2.fdtd.7.0.1).
- No removal of B2.1 `snapshot_fine_e_end_of_step` — retained
  `#[doc(hidden)]` for one release cycle for rollback safety.

## Final verification

```bash
cargo build  -p yee-fdtd
cargo clippy -p yee-fdtd --all-targets -- -D warnings
cargo test   -p yee-fdtd --release
cargo fmt    --check --all
cargo doc    --no-deps -p yee-fdtd
mdbook build docs/
```

All six must exit 0. Every existing non-`#[ignore]`'d regression in
`crates/yee-fdtd/tests/` stays green; the three gates listed above
become non-`#[ignore]`'d and pass.

## Estimated total

- LOC: ~230 (C1 ~80, C2 ~100, C3 ~20, C4 ~30; C5 ~250 conditional).
- Wall-time per agent: 2–3 days end-to-end. Critical path C1→C3 is ~1
  day; C4 adds 1 day serial (energy-balance gate run-time
  dominates). C5 adds 2 days if triggered.
- Risk concentration: Step C3 (binary pass/fail; everything downstream
  depends on the compensating source being numerically non-trivial
  AND correctly signed AND of the right order of magnitude).
