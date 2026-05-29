# Filter Phase F1.1a — `yee-voxel` — Implementation Plan

**Spec:** `2026-05-29-filter-f1-1a-yee-voxel-design.md` · **ADR:** ADR-0091

## Lane
`crates/yee-voxel/**` (new) + root `Cargo.toml` (add member). Out of lane
(yee-layout, yee-fdtd, any other crate — consume read-only) → finding, not fix.
**Do NOT modify `yee-layout` or add a yee-fdtd dep to it** (WASM-safety, ADR-0089).

## Base
Worktree `worktrees/voxel`, branch `feature/filter-f1-1a-yee-voxel`, base `3b63e87`.

## Pattern files
- `crates/yee-synth/Cargo.toml` — new-crate manifest shape (`[lints.rust]` form).
- `crates/yee-fdtd/src/grid.rs` — READ `YeeGrid::vacuum` (line 131),
  `with_eps_r_cells` (199), `with_pec_mask_ez` (335) for the EXACT `Array3`
  shapes + any internal `assert!` on shapes; build arrays to match.
- `crates/yee-layout/src/lib.rs` — `Layout`/`Polygon`/`PortRef`/`Substrate`/
  `BBox`/`Point2` field names; `Polygon.verts` for point-in-polygon.
- `crates/yee-fdtd/tests/` (cavity_resonance / lumped_lc_resonance) — how a grid
  is constructed + materials set, for idiom.

## Steps
1. `crates/yee-voxel/Cargo.toml` per spec; root `Cargo.toml` add `"crates/yee-voxel"`.
2. `src/lib.rs`: `VoxelOptions`, `MicrostripModel`, `voxelize_microstrip` per spec
   §API/§z-stack/§material. Point-in-polygon helper (ray-cast). Doc every public item.
3. `tests/voxel_001_microstrip_line.rs` per spec §DoD 4 — build a 1-trace
   microstrip `Layout` (use `yee-layout`'s public types directly, or its
   `edge_coupled_bpf` with a single section if simpler), voxelize, assert.

## Verify (exit 0; nice -n 19, --jobs 2; NO FDTD run, NO --workspace)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-voxel --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-voxel --jobs 2
nice -n 19 cargo build -p yee-voxel --jobs 2
```
`yee-voxel` pulls `yee-fdtd` (already compiled in the workspace target) — mostly
link. The voxelizer does NOT time-step; tests are ms-scale.

## Escape hatch
Blocked >15 min — `YeeGrid` array-shape mismatch (the `with_*` methods assert a
specific shape you can't satisfy), or `Layout` lacks a field you need → STOP,
commit what compiles, surface the exact shape/field gap. Do NOT modify yee-layout
or yee-fdtd to work around it (surface as a finding).

## Done when
DoD 1–5 pass; `git diff --stat 3b63e87..HEAD` shows only `crates/yee-voxel/**`
+ root `Cargo.toml`/`Cargo.lock` + the 3 committed docs; yee-layout untouched.
