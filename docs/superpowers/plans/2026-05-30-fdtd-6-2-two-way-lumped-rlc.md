# Phase 2.fdtd.6.2 — stable two-way lumped RLC port — Plan

**Spec:** `2026-05-30-fdtd-6-2-two-way-lumped-rlc-design.md` · **ADR:** ADR-0116

## Lane
`crates/yee-fdtd/**` ONLY (`src/lumped.rs` + a `tests/` gate; `ci.yml` for the
release job is allowed). Do NOT edit yee-voxel/yee-filter. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first). Branch `feature/fdtd-6-2-twoway-lumped`.

## Pattern files (READ)
- `crates/yee-fdtd/src/lumped.rs` — the CURRENT `LumpedRlcPort` (`series_rlc`,
  `pure_resistor`, `correct_e`, state fields `e_z_prev`/`inductor_current`/
  `capacitor_voltage`, the module docs that already note the one-way limitation +
  Phase 2.fdtd.6.2). This is what you rework.
- The existing lumped validation test(s) referenced in `lumped.rs` docs (the
  series-RLC / resistor gate) — must stay green.
- `crates/yee-fdtd/tests/cpml_reflection.rs` — the FDTD-test idiom (build a small
  grid, drive, measure a reflection/ratio) to mirror for `lumped_rlc_twoway_001`.
- `.github/workflows/ci.yml` `fdtd-coupling-gate` / `fdtd-lumped-gate` — the
  parallel `--release --ignored` release-gate idiom.
- Reference: Taflove & Hagness, *Computational Electrodynamics*, lumped-element
  (Piket-May 1994) chapter — the unconditionally-stable two-way semi-implicit
  R-L-C `E_z` update.

## Steps
1. Rework `correct_e` to solve the coupled `E_z^{n+1}` + branch-current/charge
   update implicitly (the stable two-way formulation) for general series R-L-C;
   verify the R-only, C-only (`l=0`), L-only (`c=∞`) limits + the Thévenin source.
   Keep constructors/signatures.
2. `tests/lumped_rlc_twoway_001.rs` (`#[ignore]`'d): stability (low-loss reactive
   element, no NaN over the record) + two-way correctness (single lumped load Γ
   vs analytic at a few f). Build + iterate IN THE CONTAINER
   (`YEE_BOX_DIR=worktrees/fdtd-6-2 scripts/yee-box.sh cargo test -p yee-fdtd
   --release -- --ignored lumped_rlc_twoway_001 --nocapture`).
3. `ci.yml`: a `fdtd-lumped-rlc-gate` release job (no `needs: lint-test`).

## Verify
- LOCAL light: fmt + `cargo clippy -p yee-fdtd --all-targets`. Existing lumped
  test still green (run `cargo test -p yee-fdtd <existing-lumped-test>` in the
  container if heavy).
- CI: `fdtd-lumped-rlc-gate` GREEN on the branch before merge.

## Escape hatch
Blocked > 60 min (the stable two-way derivation won't converge / still unstable /
the analytic-Γ cross-check won't match) → STOP + surface the update equations you
implemented, the stability behaviour, and the Γ mismatch. Do NOT weaken the gate;
do NOT regress `pure_resistor`; do NOT touch yee-voxel/yee-filter. This is
research-grade FDTD — surfacing a precise partial result is a good outcome.

## Done when
DoD 1–4; existing lumped tests green (no regression); the new two-way gate GREEN
in CI on the branch before merge; diff = `crates/yee-fdtd/**` + `ci.yml`. Then
F2.3's `fdtd_lumped_001` passes unchanged on top.
