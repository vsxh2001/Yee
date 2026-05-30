# Filter Phase F2.3-b — full-width-sheet lumped-element placement — Plan

**Spec:** `2026-05-30-f2-3-full-width-sheet-placement-design.md` · **ADR:** ADR-0124

## Lane
`crates/yee-voxel/**` ONLY (`src/lumped_sim.rs` + its module doc). Do NOT edit
`yee-fdtd` / `yee-filter`. Out of lane → finding. (`Cargo.lock` may change from the
main-merge — that's expected.)

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(F2.3 driver + gate, with `.with_two_way()` already applied; tip `3c417e6`+).
**FIRST: merge current `main` into the branch** (it predates the canonical port /
per-axis CPML / the 6.x work). `Cargo.lock` conflict → `git checkout --theirs
Cargo.lock`, `cargo check` (in container), commit (CLAUDE.md §5). ci.yml may
conflict (the F2.3 gate job vs the new fdtd gate jobs) — keep BOTH jobs.

## Pattern files (READ FIRST)
- `crates/yee-voxel/src/lumped_sim.rs` — the CURRENT driver. The element-placement
  block (`for (ri, res) in ladder.resonators...`, the `Series`/`Shunt` arms calling
  `cell_for(cx,cy,k_elem)` + `LumpedRlcPort::series_rlc(...).with_two_way()`). This
  is what you change from single-cell to a value-distributed sheet.
- `crates/yee-fdtd/tests/reactive_deembed_001.rs` (READ-ONLY) — how the bench
  builds a value-distributed full-width SHEET of lumped ports across the guide
  cross-section (the pattern to mirror for the trace width).
- The voxel model / layout: find the trace WIDTH in cells (the `N` transverse `E_z`
  edges across the trace at `k_elem`) — from `model.dims` + the layout trace width
  + the `cell_for`/`x0,y0` mapping already in the driver.
- `docs/src/decisions/0124-...md` (the value-distribution rule) + 0123 Outcome.

## Steps
1. Merge `main` into the F2.3 branch (resolve Cargo.lock `--theirs` + ci.yml keep
   both gate jobs); `cargo check -p yee-voxel` in the container green.
2. Compute `N` = trace-width-in-cells at the element x (the y-span of trace `E_z`
   edges at `k_elem`). 
3. Replace each single-cell element with a loop over the `N` transverse cells,
   emitting a value-distributed `LumpedRlcPort` per cell (shunt C → C/N, shunt
   L → N·L, both `.with_two_way()`; series arm split consistently — document).
4. Re-run `fdtd_lumped_001` in the container; capture the full |S21| sweep.
5. Update the `lumped_sim.rs` module doc (sheet placement).

## Verify (bounded container — heavy FDTD)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- `... scripts/yee-box.sh cargo test -p yee-voxel --release --test fdtd_lumped_001
  -- --ignored --nocapture` → REPORT the |S21| sweep + whether it meets the loose
  tol (GREEN) or how close it gets.
- (cargo direct or `bash -c`, NEVER `bash -lc`.)

## Escape hatch
Do NOT weaken `fdtd_lumped_001` to force GREEN. If, after value-distributed sheet
placement, the |S21| still does not meet the loose tol, that is an HONEST, valuable
result: record the achieved sweep (how close — e.g. "band-pass visible, stopband
only 8 dB") and conclude the ≈0.37 single-cell port accuracy gates it → the
multi-cell aperture port is required. Blocked > 60 min on the merge or the build →
surface. Do NOT touch yee-fdtd/yee-filter.

## Done when
fmt/clippy clean; `fdtd_lumped_001` re-run with sheet placement and the |S21| sweep
reported. Either: GREEN (→ this branch is then ready for review + the F2.3 merge,
EM-sim ships) OR a precise "how close + multi-cell port needed" finding. diff =
`crates/yee-voxel/**` (+ the main-merge artifacts).
