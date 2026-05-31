# Filter Phase F2.3-g — PEC-box 2-point standing-wave CW de-embed — Plan

**Spec:** `2026-05-31-f2-3-g-pec-2point-deembed-design.md` · **ADR:** ADR-0132

## Lane
`crates/yee-voxel/**` ONLY (`src/lumped_sim.rs` + module doc). Do NOT edit
yee-fdtd/yee-filter. Out of lane → finding. (Cargo.lock/ci.yml from the merge.)

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(tip `7efd173` — the F2.3-f matched-CPML attempt). Merge current `main` first
(Cargo.lock `--theirs`, keep all CI jobs). `cargo check -p yee-voxel` green in
container.

## Pattern files (READ FIRST)
- `crates/yee-voxel/src/lib.rs` `run_line_eeff` (ADR-0108) — the STABLE PEC +
  CW/forward-wave-on-a-line pattern for this microstrip geometry (NO CPML — CPML
  into the substrate is unstable here, the whole reason F2.3-f failed). Mirror its
  PEC box + standing-wave handling.
- `crates/yee-voxel/src/lumped_sim.rs` `run_board_solve` — the current CW drive +
  the F2.3-f CPML termination you REPLACE with a PEC box + 2-point standing-wave
  de-embed. Keep the aperture-port placement + the CW drive.
- `crates/yee-fdtd/tests/cap_cw_001.rs` — the CW steady-state phasor measurement
  idiom (single-bin DFT over a settled window) for sampling V(x) phasors.
- ADR-0132 (the 2-point method) + ADR-0131 Outcome (why CPML failed; over-unity =
  bad de-embed) + ADR-0108 (PEC-not-CPML for microstrip).

## Steps
1. Merge `main`; `cargo check -p yee-voxel` green (container).
2. PEC box: drop the CPML termination; run the microstrip in a PEC-bounded grid,
   line long enough past each port for a developed standing wave + element
   clearance. (Stable, unlike CPML-into-substrate.)
3. CW drive (Hann-ramped, settle to steady state). At each port reference region,
   single-bin-DFT the line-voltage phasor at ≥2 points of spacing `d`. Solve
   `V(x)=a e^{−jβx}+b e^{+jβx}` for forward `a` / backward `b` (β from a thru ε_eff
   calibration or a 3-point fit).
4. `S21(f) = (b₂/a₁)_dut / (b₂/a₁)_thru`. Frequency set = gate points + a few.
5. Re-run `fdtd_lumped_001`; capture the |S21| (physical? band-pass? notch depth).

## Verify (bounded container — heavy FDTD, generous timeout; keep dx=0.4mm, bounded freq set)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- `... scripts/yee-box.sh cargo test -p yee-voxel --release --test fdtd_lumped_001
  -- --ignored --nocapture` → REPORT the |S21| + verdict. (cargo direct or
  `bash -c`, NEVER `bash -lc`. Keep the run bounded — dx=0.4mm, few freqs — to avoid
  the multi-hour runs that bit F2.3-e/f.)

## Escape hatch
Do NOT weaken `fdtd_lumped_001`. Three honest outcomes (all valuable): band-pass
≥20 dB (ship), physical-but-shallow band-pass (→ port accuracy next), or
still-monotone/no-band-pass (→ the board integration genuinely doesn't resonate, a
deeper finding to surface). If the 2-point β-extraction is unstable or the PEC box
has a cavity mode at a gate freq, surface the precise issue. Blocked > 90 min OR a
single run > ~45 min → surface (don't burn multi-hour runs). Do NOT touch
yee-fdtd/yee-filter; do NOT merge/push.

## Done when
fmt/clippy clean; `fdtd_lumped_001` re-run with the PEC-box 2-point de-embed + the
|S21| reported. Either GREEN at 20 dB (→ review + F2.3 merge, EM-sim ships 6/6) OR a
precise verdict (shallow band-pass → port accuracy; or no band-pass → board
integration). diff = `crates/yee-voxel/**` (+ merge artifacts).
