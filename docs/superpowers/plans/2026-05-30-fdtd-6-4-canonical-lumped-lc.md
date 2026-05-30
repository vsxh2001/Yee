# Phase 2.fdtd.6.4 — canonical per-element Taflove lumped L/C updates — Plan

**Spec:** `2026-05-30-fdtd-6-4-canonical-lumped-lc-design.md` · **ADR:** ADR-0118

## Lane
`crates/yee-fdtd/**` ONLY (`src/lumped.rs` + `tests/lumped_rlc_twoway_001.rs`).
May READ `crates/yee-voxel/src/lumped_sim.rs` for context; do NOT edit it or
yee-filter. Out of lane → finding, not fix.

## Base
New worktree off `main` (re-fetch first). Branch `feature/fdtd-6-4-canonical-lc`.
main has ADR-0116/0117/0118 docs merged.

## Pattern files (READ FIRST)
- `docs/src/decisions/0117-fdtd-6-3-reactive-magnitude.md` — the **Outcome**
  section: WHY the current update fails (loads by instantaneous `K`, not the
  physical reactance; opposite-direction L/C errors). Do NOT re-attempt a
  coefficient rescale — that path is closed.
- `crates/yee-fdtd/src/lumped.rs` — the CURRENT `correct_e`. Keep the resistor
  path; replace the reactive arms with the canonical per-element updates.
- `crates/yee-fdtd/tests/lumped_rlc_twoway_001.rs` — the gate. It already sweeps
  shunt-L, shunt-C, series-RLC, prints |Γ|_fdtd vs the analytic `−Z₀/(2Z_L+Z₀)`
  shunt law with a scalar calibration `A`, and confirms the resistor exactly.
  Reuse that harness; turn the reactive prints into asserts.
- Taflove & Hagness Ch. on lumped elements (Piket-May 1994): canonical
  capacitor (effective permittivity `ε_eff = ε₀ + C·dz/dA`), inductor
  (accumulated current `I_L += (dt·dz/L)·E_z`, subtracted from the `E_z` update),
  series-RLC combined `E` update.

## Steps
1. **Capacitor arm** (`l=0`): implement `ε_eff = ε₀ + C·dz/dA` (equivalently the
   `I_C = C·dz·dE_z/dt` displacement-current injection). Verify shunt-C |Γ| →
   analytic in the container.
2. **Inductor arm** (`c=∞`): implement the auxiliary accumulated current
   `I_L^{n+1/2} = I_L^{n−1/2} + (dt·dz/L)·E_z^n`, subtract `(dt/(ε₀·dA))·I_L` from
   the `E_z` update. Verify shunt-L |Γ| → analytic.
3. **Series R-L-C**: the canonical combined update (R + inductor accumulator +
   capacitor-voltage state in series). Verify series-RLC |Γ|. If it resists the
   tol after the shunt cases pass, invoke the escape hatch (ship shunt, defer
   series-RLC with an in-test note).
4. **Resistor**: confirm byte-identical / numerically-identical to the validated
   path (no regression).
5. **Strengthen the gate**: reactive prints → asserts (Δ|Γ| ≤ 0.15 at 4/6/9 GHz,
   after the test's scalar calibration). Keep resistor-exact + stability asserts.
6. Iterate IN THE CONTAINER (fast, ~3 s/run):
   `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
   --release --test lumped_rlc_twoway_001 -- --ignored --nocapture`
   (run cargo directly or via `bash -c '…'`, NEVER `bash -lc` — it drops cargo
   from PATH in-container).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd`, `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container, `bash -c`).
- No regression: `cargo test -p yee-fdtd --release --test lumped_lc_resonance
  --test lumped_resistor -- --include-ignored` GREEN.
- Gate: `lumped_rlc_twoway_001` GREEN with reactive asserts.

## Escape hatch
Blocked > 60 min, OR a reactive arm cannot reach the loose tol while staying
stable → STOP and surface: the update you implemented, its measured |Γ| table,
and which arm fails. Do NOT weaken the gate, relax the resistor tol, fake, or
touch yee-voxel/yee-filter. Shipping shunt-L + shunt-C and deferring series-RLC
(documented) is an acceptable partial — that already unblocks F2.3's selectivity.

## Done when
DoD 1–4; resistor + fdtd-206 non-regressed; reactive asserts GREEN in the
container (series-RLC green or explicitly deferred); diff = `crates/yee-fdtd/**`.
Then F2.3's `fdtd_lumped_001` is re-run on top (separate follow-up).
