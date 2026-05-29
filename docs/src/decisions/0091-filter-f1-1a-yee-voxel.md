# ADR-0091: Filter Phase F1.1a — `yee-voxel` Layout→FDTD voxelizer

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0086 (`yee-layout`), ADR-0089 (app architecture / WASM-safety),
`FILTER-DESIGN-ROADMAP.md` §5/§5a (F1.1)

---

## Context

`FILTER-DESIGN-ROADMAP` F1.1 (FDTD coupling/Qe extraction) is the EM-in-the-loop
building block for planar dimensional synthesis. A read-only recon of `yee-fdtd`
(2026-05-29) found **every FDTD primitive already exists** — per-cell ε_r
(`YeeGrid::with_eps_r_cells`), PEC masks (`with_pec_mask_ez`), `set_sigma_box`,
`LumpedRlcPort` (drive + `inductor_current`/`capacitor_voltage`), the single-bin
DFT pattern (`cavity_resonance.rs`), and decay-fit Q (`cavity_q.rs`). The **one
missing block** is a converter from a `yee-layout::Layout` (polygons on a
substrate) to a material-assigned `YeeGrid`. Nothing rasterizes geometry to a
voxel grid in-tree.

So F1.1 splits: **F1.1a = the voxelizer** (this ADR; bounded, gateable with NO
FDTD run) → **F1.1b = k/Qe extraction** (composes the existing primitives; heavy
runs; separate increment).

## Decision

New native crate **`yee-voxel`** (deps `yee-layout`, `yee-fdtd`, `ndarray`):
turn a planar microstrip `Layout` into a material-assigned `YeeGrid` for FDTD.

```rust
pub struct MicrostripModel {
    pub grid: yee_fdtd::YeeGrid,            // ε_r + PEC assigned
    pub dims: (usize, usize, usize),        // (nx, ny, nz)
    pub dx_m: f64,
    pub port_cells: Vec<(usize, usize, usize)>,  // layout ports → grid cells
}
pub struct VoxelOptions { pub dx_m: f64, pub xy_margin_cells: usize, pub air_above_cells: usize }
pub fn voxelize_microstrip(layout: &yee_layout::Layout, opts: &VoxelOptions) -> MicrostripModel;
```

The microstrip z-stack (z = up): ground plane (PEC) at the bottom cell layer;
substrate slab of `ε_r = layout.substrate.eps_r` for `round(height_m/dx)` cell
layers; the top-metal traces (PEC, one cell thick) at the substrate-top layer,
rasterized by point-in-polygon over each trace `Polygon`; air (ε_r = 1) above.
X-Y extent from `layout.bbox` + `xy_margin_cells`. The `eps_r_cells`
(`Array3<f64>`) and the PEC masks are built at the exact shapes the `YeeGrid`
builders require, then fed to `YeeGrid::vacuum(..).with_eps_r_cells(..)…`.

> **Review correction (shipped impl).** This ADR's prose said `pec_mask_ez`. That
> is the WRONG component: a horizontal PEC sheet (ground plane, traces) zeroes the
> **tangential** field — `Ex` and `Ey` — on its plane, not the normal `Ez`. The
> shipped voxelizer masks `pec_mask_ex` (shape `(nx, ny+1, nz+1)`) **and**
> `pec_mask_ey` (`(nx+1, ny, nz+1)`) at the ground (k=0) and trace (k_top) planes,
> at their staggered node positions (`Ex` at `((i+0.5)dx, j·dx)`, `Ey` at
> `(i·dx, (j+0.5)dx)`), and chains `.with_pec_mask_ex(..).with_pec_mask_ey(..)`.
> Read `pec_mask_ez` below as `pec_mask_ex`+`pec_mask_ey`.

**No EM:** building the grid assigns materials only — there is no time-stepping,
so the gate runs in milliseconds.

## Consequences

**Ships:** `yee-voxel` crate + `voxelize_microstrip`. Gate (crate test, §4, NO
FDTD run): voxelize a single straight microstrip line (one rectangle trace, FR-4
`ε_r=4.4`, `h=1.6 mm`) and assert — grid dims match the bbox+margin / dx +
substrate+air layer counts; substrate cells carry `ε_r≈4.4`, air cells `1.0`;
the ground-plane layer is all PEC; the trace-layer PEC-cell count ≈
`trace_area/dx²` (within rounding); a cell under the trace is PEC, one clearly
off it is not; `port_cells` map the layout's ports to in-range cell indices.

**Constraint honored (ADR-0089):** `yee-layout` is NOT modified and gains no
`yee-fdtd` dep — it stays WASM-safe. `yee-voxel` is the separate native crate
that bridges the two.

**Not in scope:** the k/Qe extraction + any FDTD run (F1.1b); waveguide/lumped
voxelization. New dep: `ndarray` (already a `yee-fdtd` dep — in the lockfile).

---

## References
- `yee-fdtd/src/grid.rs` (`vacuum`/`with_eps_r_cells`/`with_pec_mask_ez`);
  the F1.1 recon; `yee-layout` `Layout`/`Polygon`/`Substrate`.
- `docs/superpowers/specs/2026-05-29-filter-f1-1a-yee-voxel-design.md`;
  `docs/superpowers/plans/2026-05-29-filter-f1-1a-yee-voxel.md`.
