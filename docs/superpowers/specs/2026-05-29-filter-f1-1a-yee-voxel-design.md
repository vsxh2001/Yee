# Filter Phase F1.1a — `yee-voxel` Layout→FDTD voxelizer — Design Spec

**Phase:** F1.1a · **ADR:** ADR-0091 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
Convert a planar microstrip `yee-layout::Layout` into a material-assigned
`yee-fdtd::YeeGrid` (ε_r slab + PEC ground + PEC traces). The one missing block
before FDTD k/Qe extraction (F1.1b). Pure rasterization + array assignment — no
time-stepping, so gateable in milliseconds. Native crate (deps yee-fdtd).

## Crate `yee-voxel`
`Cargo.toml`: deps `yee-layout`, `yee-fdtd`, `ndarray` (all `{ workspace = true }`).
`[lints.rust] unsafe_code="forbid"`, `missing_docs="warn"`.

### API
```rust
pub struct VoxelOptions {
    pub dx_m: f64,               // isotropic cell size
    pub xy_margin_cells: usize,  // air margin around the bbox in x and y
    pub air_above_cells: usize,  // air layers above the top metal
}
pub struct MicrostripModel {
    pub grid: yee_fdtd::YeeGrid,
    pub dims: (usize, usize, usize),
    pub dx_m: f64,
    pub port_cells: Vec<(usize, usize, usize)>,
}
pub fn voxelize_microstrip(layout: &yee_layout::Layout, opts: &VoxelOptions) -> MicrostripModel;
```

### Z-stack (z increasing upward)
- `k = 0`: ground plane — PEC over the whole x-y extent.
- `k = 1 .. 1+n_sub`: substrate, `n_sub = (substrate.height_m / dx).round()` (≥1),
  `ε_r = substrate.eps_r`.
- `k = k_top = 1 + n_sub`: top-metal layer — PEC where a trace `Polygon` covers
  the cell centre (point-in-polygon); ε_r = 1 elsewhere.
- `k_top+1 .. k_top+1+air_above_cells`: air (ε_r = 1).
- `nz = k_top + 1 + air_above_cells`.

### X-Y extent
`x0 = bbox.min.x − margin·dx`, `x1 = bbox.max.x + margin·dx` (same in y);
`nx = ceil((x1−x0)/dx)`, `ny = ceil((y1−y0)/dx)`. Cell (i,j) centre =
`(x0 + (i+0.5)dx, y0 + (j+0.5)dx)`.

### Material assignment (match YeeGrid's exact array shapes)
- Read `YeeGrid::with_eps_r_cells` / `with_pec_mask_ez` in `yee-fdtd/src/grid.rs`
  for the EXACT `Array3` shapes they require (eps is `[nx+1,ny+1,nz+1]`-style;
  pec_mask_ez matches the staggered Ez component) — build arrays of those shapes,
  fill per the z-stack, then `YeeGrid::vacuum(nx,ny,nz,dx).with_eps_r_cells(eps)
  .with_pec_mask_ez(pec)`. Point-in-polygon: standard ray-cast / winding test.
- `port_cells`: map each `layout.ports[p].at` (x,y) to its `(i,j,k_top)` cell.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-voxel --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-voxel` exit 0 (fast; NO FDTD time-stepping).
4. `voxel_001_microstrip_line` (`tests/`): build a single straight microstrip
   line `Layout` (one rectangle trace `w×l`, FR-4 `ε_r=4.4`, `h=1.6 mm`),
   `voxelize_microstrip` with a chosen `dx`, and assert:
   - `dims` match the hand-computed `nx,ny,nz` (bbox+margin/dx; `n_sub` substrate
     + 1 ground + 1 trace-layer + `air_above` layers).
   - a substrate cell has `ε_r` within 1e-9 of 4.4; an air cell `1.0`.
   - the ground-plane layer (k=0) is fully PEC in `pec_mask_ez`.
   - the trace-layer PEC-cell count ≈ `(w·l)/dx²` within ±1 row/col of rounding.
   - a cell whose centre is under the trace is PEC; one clearly off the trace
     (in the margin) is not.
   - `port_cells.len() == layout.ports.len()` and each is in `[0,nx)×[0,ny)×[0,nz)`.
5. `cargo build -p yee-voxel` exit 0.

## Out of scope
Any FDTD run / k/Qe extraction (F1.1b); waveguide or lumped voxelization;
sub-cell/conformal metal. Cubic-cell occupancy only.
