# Phase 2.fdtd.6.6 — reactive lumped-port reformulation (sheet→mode coupling) — Plan

**Spec:** `2026-05-30-fdtd-6-6-reactive-port-reformulation-design.md` · **ADR:** ADR-0121

## Lane
`crates/yee-fdtd/**` ONLY (`src/lumped.rs` + `tests/reactive_deembed_001.rs`; may
touch `tests/lumped_rlc_twoway_001.rs` only if a shared change requires it). May
READ yee-voxel for context. Edit nothing else. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first). Branch `feature/fdtd-6-6-reactive-port`.
main already has the de-embed bench (`reactive_deembed_001`) + the canonical
per-element port (the foundation you reformulate).

## Pattern files (READ FIRST)
- `docs/src/decisions/0121-fdtd-6-6-reactive-port-reformulation.md` (this ADR — the
  hypotheses in order) and `0119-...` Outcome (the bench + the PORT-WRONG numbers:
  capacitor `Z_in≈94 Ω` vs `~3175 Ω` expected, the well-conditioned signal).
- `crates/yee-fdtd/tests/reactive_deembed_001.rs` — YOUR HARNESS. It already
  measures `Z_L(ω)` (V=∫E·dz, I=∮H·dl modal current, Z₀ from the incident wave,
  resistor anchor κ). Use it to validate every change; it counts the sheet cells
  and prints per-load `Z_L` tables.
- `crates/yee-fdtd/src/lumped.rs` — the canonical per-element `correct_e` (verified
  correct per-EDGE). The reformulation is the SHEET→MODE coupling, NOT the
  per-edge constitutive — do not "fix" the per-edge math.

## Steps (each validated against the bench, ~3 s/run)
1. Determine `N` (interior transverse `E_z` edges the full-width sheet spans) in
   the bench geometry. Implement **hypothesis 1 (value-normalization)**: shunt
   `C → C/N` per cell, shunt `L → N·L` per cell (decide where: a sheet-aware
   constructor in `lumped.rs`, or the bench's sheet-build helper — keep the public
   single-element API stable). Re-run the bench; check the well-conditioned
   capacitor `Z_L` moves from ~94 Ω toward ~3175 Ω and the resistor anchor still
   holds.
2. If insufficient, **hypothesis 2 (modal coupling factor)**: reconcile the
   per-cell back-action / source with `V=∫E·dz`, `I=∮H·dl` using the resistor's
   measured `κ`. Re-run.
3. If still wrong, **hypothesis 3 (modal lumped port)** — enforce `V=Z_L·I` on the
   measured modal V/I. Larger; only if (1)/(2) fail.
4. Once the reactive arms match: turn the bench's reactive prints into ASSERTS
   (shunt-C at minimum, within a loose Δ of `1/(jωC)`; shunt-L/series-RLC as
   conditioning allows), flip the module verdict to PORT-CORRECT, keep the
   resistor anchor asserted.
- Container loop:
  `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
  --release --test reactive_deembed_001 -- --ignored --nocapture`
  (cargo direct or `bash -c '…'`, NEVER `bash -lc`).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd` + `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container, `bash -c`) → exit 0.
- No regression: `cargo test -p yee-fdtd --release --test lumped_lc_resonance
  --test lumped_resistor --test lumped_rlc_twoway_001 -- --include-ignored` GREEN
  (resistor anchor exact).
- Gate: `reactive_deembed_001` GREEN with the shunt-C arm asserted to match.

## Escape hatch
Blocked > 60 min, OR hypotheses (1)→(3) cannot bring the well-conditioned
capacitor within tol while staying stable → STOP and surface: each coupling you
tried and the bench `Z_L(ω)` it produced (the bench quantifies the residual), and
whether the remaining scope is the multi-week modal port (3). Do NOT weaken the
resistor anchor, fake the reactive arms, or touch yee-voxel/yee-filter. A measured
partial (e.g. "shunt-C now matches, shunt-L still ill-conditioned") is acceptable
and is real progress — record it.

## Done when
DoD 1–4: shunt-C arm asserted GREEN on the bench (verdict PORT-CORRECT for at least
the well-conditioned case), resistor + fdtd-206 non-regressed; diff =
`crates/yee-fdtd/**`. Then F2.3 is re-run on top (separate follow-on).
