# Filter Phase F2.3-e — finer-grid lumped EM sim — Plan

**Spec:** `2026-05-31-f2-3-e-finer-grid-emsim-design.md` · **ADR:** ADR-0129

## Lane
`crates/yee-voxel/**` ONLY (`src/lumped_sim.rs` + a scratch sweep test + the gate's
dx if changed). Do NOT edit yee-fdtd/yee-filter. Out of lane → finding.

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(tip `6373bc7` — aperture port + CW drive). Merge current `main` first (Cargo.lock
`--theirs`, keep all CI jobs). `cargo check -p yee-voxel` green in container.

## Pattern files (READ FIRST)
- `crates/yee-voxel/src/lumped_sim.rs` — the CW drive (`run_board_solve`,
  `cw_settle_cycles`, `cw_freqs_hz`) + `LumpedSimConfig.dx_m`. The settle is in
  *carrier cycles*; confirm it scales physically as dt shrinks with dx (so finer dx
  doesn't under-settle).
- `crates/yee-voxel/tests/fdtd_lumped_001.rs` — the gate (its config + the 2.0/2.4
  GHz checks).
- ADR-0129 (the dx-refinement decision) + ADR-0128 Outcome (the coarse-grid 5 dB
  saturation + the over-unity passband artifact).

## Steps
1. Merge `main` into the branch. `cargo check -p yee-voxel` green (container).
2. Scratch dx-sweep test: `simulate_lumped_board` at dx = 0.4 / 0.2 mm (0.1 mm if
   runtime allows) at {2.0, 2.4 GHz}, CW drive (settle scaled to hold physical
   settle time as dt shrinks). Print the 2.4 GHz rejection + 2.0 GHz |S21| + the
   wall-clock per dx. Run in the BOUNDED CONTAINER (generous timeout; finer dx is
   N⁴-expensive — if 0.1 mm is too slow, do 0.2 mm + extrapolate and `log` it).
3. **Decide:** rejection climbs toward 20 dB ⇒ find the dx meeting the gate, set the
   F2.3 default dx, re-run `fdtd_lumped_001` (strict 20 dB). Caps shallow ⇒ record
   the trend + conclude a higher-accuracy port is needed (next ADR).
4. Remove the scratch test before finishing (untracked); keep only a principled
   default-dx change if the gate passes.

## Verify (bounded container — heavy)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- The dx-sweep + (if a dx passes) `fdtd_lumped_001` at the strict gate, in the
  container. (cargo direct or `bash -c`, NEVER `bash -lc`.)

## Escape hatch
Do NOT weaken `fdtd_lumped_001`. If finer dx caps the rejection shallow (doesn't
climb toward 20 dB) OR the runtime is prohibitive even in the container, that is an
honest result: report the rejection-vs-dx trend + runtime + the extrapolated dx (or
"caps at X dB → higher-accuracy port needed"). Blocked > 90 min → surface. Do NOT
touch yee-fdtd/yee-filter.

## Done when
Either: `fdtd_lumped_001` GREEN at the strict 20 dB at a feasible finer dx (→ branch
ready for review + the F2.3 merge, EM-sim ships) OR a precise "rejection vs dx +
runtime, caps/extrapolation" finding driving the next sub-increment. diff =
`crates/yee-voxel/**` (+ merge artifacts).
