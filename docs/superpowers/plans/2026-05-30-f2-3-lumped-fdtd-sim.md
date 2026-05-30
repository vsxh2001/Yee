# Filter Phase F2.3 — lumped-LC FDTD EM simulation — Plan

**Spec:** `2026-05-30-f2-3-lumped-fdtd-sim-design.md` · **ADR:** ADR-0115

## Lane
`crates/yee-voxel/**` (driver + `#[ignore]`'d gate) + `.github/workflows/ci.yml`
(release job). Add `yee-filter` to yee-voxel's deps. Consume yee-fdtd
`LumpedRlcPort` + yee-layout + yee-filter `lumped_board`/`LumpedLadder` public
APIs — do NOT edit them. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first). Branch `feature/filter-f2-3-lumped-fdtd`.

## Pattern files (READ)
- `crates/yee-voxel/src/lib.rs` — `voxelize_microstrip`, `run_line_eeff` (the
  F1.1b.1 FDTD driver: solver build, step loop, lumped ports, DFT, time-gating).
  MIRROR the driver shape.
- `crates/yee-fdtd/src/lumped.rs` — `LumpedRlcPort::series_rlc(cell,r,l,c,src)`
  (+ `pure_resistor`, `correct_e`, `SourceWaveform`); the series/shunt decomposition
  in the spec.
- `crates/yee-filter/src/board.rs` — `lumped_board(&LumpedLadder,&Substrate,Footprint)
  -> LumpedBoard{layout,placements}`; `Placement{ref_des,footprint,kind,center_m}` →
  map each placement's `center_m` to a grid `(i,j,k_top)` cell (mirror how
  voxelize maps ports).
- `crates/yee-filter/src/lumped.rs` — `LumpedLadder`/`LcResonator`/`LcBranch`.
- `.github/workflows/ci.yml` `fdtd-coupling-gate` (the parallel `--release
  --ignored` release-gate idiom from F1.1b.1) — mirror for `fdtd-lumped-gate`.
- CLAUDE.md §10 — the voxelizer z-stack gotcha + the bounded-container usage.

## Steps
1. yee-voxel Cargo.toml: add yee-filter.
2. `simulate_lumped_board` + `LumpedSimConfig`: build board → voxelize → for each
   placement, map to a cell + add the per-branch lumped element(s) (series: 1
   series_rlc(L,C); shunt: pure-L ‖ pure-C two elements) → drive input port,
   matched output → step loop (correct_e for every port after update_e) → DFT
   input+output → S21(f). Document public items.
3. `tests/fdtd_lumped_001.rs`: `#[ignore]`'d (reason string); build the cheb N=5
   ladder, `simulate_lumped_board`, assert |S21| in-band ≈ 0 dB (within a few dB)
   and ≥ ~20 dB rejection at 2.4 GHz, cross-checked vs `ladder_s21`.
4. ci.yml: `fdtd-lumped-gate` job (release, no `needs: lint-test`) running
   `cargo test -p yee-voxel --release -- --ignored fdtd_lumped_001 --nocapture`.

## Verify
- LOCAL (light): fmt + `cargo clippy -p yee-voxel --all-targets`. Build + run the
  gate IN THE CONTAINER (`scripts/yee-box.sh`, host-safe) until green — paste the
  FDTD vs analytic S21.
- CI: push branch; the `fdtd-lumped-gate` job GREEN before merge.

## Escape hatch
Blocked > 45 min (the shunt parallel-LC via two elements doesn't behave;
S21 extraction wrong; coarse FDTD can't hit even a loose tol; cell-mapping a
placement to a gap is ambiguous) → STOP + surface the measured vs analytic S21,
the grid, and the lumped-element placement. Do NOT weaken the gate to a no-op; do
NOT run FDTD on the host (container only); do NOT edit yee-filter/yee-fdtd.

## Done when
DoD 1–4; the CI `fdtd-lumped-gate` GREEN on the branch before merge; diff =
`crates/yee-voxel/**` + `ci.yml` (+ Cargo.lock); WASM-safe crates untouched.
