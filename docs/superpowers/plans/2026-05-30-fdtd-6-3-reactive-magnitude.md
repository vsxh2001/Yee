# Phase 2.fdtd.6.3 — reactive-magnitude correctness of the two-way lumped port — Plan

**Spec:** `2026-05-30-fdtd-6-3-reactive-magnitude-design.md` · **ADR:** ADR-0117

## Lane
`crates/yee-fdtd/**` ONLY (`src/lumped.rs` + `tests/lumped_rlc_twoway_001.rs`).
May READ `crates/yee-voxel/src/lumped_sim.rs` (the F2.3 consumer) for context but
do NOT edit it or yee-filter. Out of lane → finding, not fix.

## Base
New worktree off `main` (re-fetch first). Branch `feature/fdtd-6-3-reactive-mag`.
main has ADR-0116 (the two-way port) already merged.

## Pattern files (READ FIRST)
- `crates/yee-fdtd/src/lumped.rs` — the CURRENT two-way `correct_e` (the
  `K = R + L/dt + dt/(2C)`, `β`, `I^{n+1/2}`, `E_z^{n+1}`, `V_C` update). This is
  what you derive `Z_d(ω)` from and fix.
- `crates/yee-fdtd/tests/lumped_rlc_twoway_001.rs` — the gate. It already SWEEPS
  pure-L, pure-C, series-RLC and PRINTS |Γ|_fdtd vs |Γ|_anal with a scalar
  calibration `A`. You convert those prints into asserts; reuse its calibration.
- The diagnostic data (ADR-0117 table): R exact, L transparent (|Γ|≈0.013), C
  near-open (|Γ|≈1.0). This tells you the sign of the bug: L under-couples, C
  over-couples → the reactive coefficients carry a wrong `dz`/`dA`/`ε₀`/`dt`
  factor relative to the (correct) resistor + field-coupling terms.

## Steps
1. **Derive `Z_d(ω)`.** Z-transform the discrete branch recurrences (`I`, `V_C`,
   `E_z`) to get the discrete impedance the port presents; evaluate at
   `z = e^{jωdt}`. In the low-`ωdt` limit it must → `R + jωL + 1/(jωC)`. The R
   term matches; find the mis-scaled L and/or C factor. WRITE the derivation in a
   comment or the report — this is derivation-first, not parameter-fishing.
2. **Fix** the reactive coefficient(s) in `correct_e`. Keep `K + β > 0`
   (unconditional stability) and the exact resistor limit. Keep public API.
3. **Strengthen the gate**: turn the reactive |Γ| prints into asserts (pure-L,
   pure-C, series-RLC within Δ|Γ| ≤ 0.15 after the existing scalar calibration,
   at 4/6/9 GHz). Keep resistor-exact + stability asserts.
4. Iterate IN THE CONTAINER (fast, ~3 s/run):
   `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
   --release --test lumped_rlc_twoway_001 -- --ignored --nocapture`
   (run cargo via `bash -c '…'` or directly — NOT `bash -lc`, which drops cargo
   from PATH in-container).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd`, `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container, `bash -c`).
- No regression: `cargo test -p yee-fdtd --release --test lumped_lc_resonance
  --test lumped_resistor -- --include-ignored` GREEN (container).
- Gate: `lumped_rlc_twoway_001` GREEN with the reactive asserts (container).

## Escape hatch
Blocked > 60 min, OR the reactive |Γ| cannot be brought within the loose tol
without destabilising the update after deriving `Z_d(ω)` → STOP and surface: the
derived `Z_d(ω)`, the coefficient you changed and why, and the residual |Γ| table.
Do NOT weaken the gate back to a print; do NOT relax the resistor-exact tol; do
NOT fake a pass; do NOT touch yee-voxel/yee-filter. A precise partial (the
derivation + residual) is a good outcome — this is research-grade FDTD.

## Done when
DoD 1–4; resistor + fdtd-206 + resistor gates non-regressed; the reactive
asserts GREEN in the container; diff = `crates/yee-fdtd/**` only. Then F2.3's
`fdtd_lumped_001` is re-run on top (separate follow-up) and should acquire its
band-pass shape.
