//! # yee-voxel
//!
//! Native bridge from a planar microstrip [`yee_layout::Layout`] to a
//! material-assigned [`yee_fdtd::YeeGrid`] (Filter Phase F1.1a, ADR-0091).
//!
//! [`voxelize_microstrip`] rasterizes the layout's top-metal polygons onto a
//! cubic Yee grid: a PEC ground plane at the bottom cell layer, a dielectric
//! substrate slab of `ε_r = layout.substrate.eps_r`, a one-cell-thick PEC
//! top-metal layer where a trace polygon covers the cell centre
//! (point-in-polygon ray-cast), and air above. The result is a `YeeGrid` with
//! per-cell `ε_r` and tangential `Ex`+`Ey` PEC masks already attached (a
//! horizontal PEC sheet zeroes the in-plane field, not the normal `Ez`), ready
//! for the F1.1b k/Q_e extraction step.
//!
//! This crate does **no** EM time-stepping — building the grid assigns
//! materials only, so it runs in milliseconds.
//!
//! ## WASM-safety boundary (ADR-0089)
//!
//! `yee-layout` deliberately has no `yee-fdtd` dependency so it stays
//! WASM-safe. `yee-voxel` is the separate **native** crate that bridges the
//! two; it depends on both.
//!
//! ## Z-stack (`z` increasing upward)
//!
//! The ground plane sits on the plane `z = 0` and the trace on `z = k_top·dx`,
//! so the dielectric fills the `n_sub` `E_z` edges spanning the ground-to-trace
//! gap (`k = 0 .. n_sub`) with **no air series gap** at the ground.
//!
//! - `k = 0` (plane `z = 0`): ground plane — tangential-E PEC over the whole
//!   x-y extent.
//! - `E_z` edges `k = 0 .. n_sub` (`z ∈ [0, n_sub·dx]`): substrate,
//!   `n_sub = round(height_m / dx)` (≥ 1), `ε_r = layout.substrate.eps_r`.
//! - `k = k_top = n_sub` (plane `z = n_sub·dx ≈ h`): top-metal layer —
//!   tangential-E PEC where a trace polygon covers the cell centre.
//! - `E_z` edges `k = k_top .. k_top + air_above_cells`: air (`ε_r = 1`).
//! - `nz = k_top + 1 + air_above_cells`.
//!
//! ## Example
//!
//! ```
//! use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
//! use yee_voxel::{voxelize_microstrip, VoxelOptions};
//!
//! let substrate = Substrate {
//!     eps_r: 4.4,
//!     height_m: 1.6e-3,
//!     loss_tangent: 0.0,
//!     metal_thickness_m: 35e-6,
//! };
//! let trace = Polygon::rect(0.0, 0.0, 3.0e-3, 20.0e-3);
//! let traces = vec![trace];
//! let bbox = BBox::from_polygons(&traces);
//! let layout = Layout {
//!     substrate,
//!     traces,
//!     ports: vec![PortRef {
//!         at: Point2::new(1.5e-3, 0.0),
//!         width_m: 3.0e-3,
//!         ref_impedance_ohm: 50.0,
//!     }],
//!     bbox,
//! };
//! let opts = VoxelOptions { dx_m: 0.5e-3, xy_margin_cells: 4, air_above_cells: 8 };
//! let model = voxelize_microstrip(&layout, &opts);
//! assert_eq!(model.port_cells.len(), 1);
//! ```

use ndarray::Array3;
use yee_fdtd::{LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver, YeeGrid};
use yee_layout::{Layout, Point2, Polygon, Stackup};

/// Lumped-LC FDTD EM simulation of a synthesized filter board (Filter Phase
/// F2.3, ADR-0115): place each ladder L/C as a [`yee_fdtd::LumpedRlcPort`] on
/// the voxelized board, drive/sense two ports, and extract `|S21|(f)`.
pub mod lumped_sim;
pub use lumped_sim::{LumpedSimConfig, SERIES_ESR_OHM, simulate_lumped_board};

/// Voxelization parameters for [`voxelize_microstrip`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoxelOptions {
    /// Isotropic cell size `dx = dy = dz`, metres.
    pub dx_m: f64,
    /// Air margin (in cells) added around the layout bounding box in `x` and `y`.
    pub xy_margin_cells: usize,
    /// Number of air layers placed above the top-metal layer.
    pub air_above_cells: usize,
}

/// A material-assigned FDTD model produced from a planar microstrip layout.
#[derive(Debug)]
pub struct MicrostripModel {
    /// Ground-sheet layer index (`E_z` edge below the substrate): `0` for
    /// the classic floor-ground stack ([`voxelize_microstrip`]); the
    /// `air_below_cells` of a lifted stack ([`voxelize_microstrip_open`]).
    pub k_gnd: usize,
    /// The Yee grid with per-cell `ε_r` and tangential `Ex`+`Ey` PEC masks attached.
    pub grid: YeeGrid,
    /// Grid cell dimensions `(nx, ny, nz)`.
    pub dims: (usize, usize, usize),
    /// Cell size used to build the grid, metres (echoes [`VoxelOptions::dx_m`]).
    pub dx_m: f64,
    /// Each layout port mapped to its `(i, j, k_top)` grid cell.
    pub port_cells: Vec<(usize, usize, usize)>,
}

/// Voxelize a planar microstrip [`Layout`] into a material-assigned
/// [`YeeGrid`].
///
/// Builds a cubic-cell grid sized from `layout.bbox` plus
/// [`VoxelOptions::xy_margin_cells`] of air in `x`/`y`, with the z-stack
/// described in the [crate] docs. Returns the assigned grid together with its
/// dimensions, the cell size, and the per-port cell indices.
///
/// # Panics
///
/// Panics if [`VoxelOptions::dx_m`] is not positive and finite, or if the
/// resulting grid would have a zero dimension (an empty / degenerate bounding
/// box with no margin).
pub fn voxelize_microstrip(layout: &Layout, opts: &VoxelOptions) -> MicrostripModel {
    voxelize_with_ground_level(layout, opts, 0)
}

/// Voxelize with a **lifted** z-stack (FS.1a.1b, ADR-0205): `air_below_cells`
/// air layers under the ground sheet, so the domain floor is free air (an
/// absorber can sit there) and the ground is a finite PEC sheet mid-domain —
/// what a truncated-ground antenna (quasi-Yagi) needs, where the classic
/// floor-ground stack lets the PEC bottom boundary act as an infinite image
/// plane no mask truncation can remove. `air_below_cells = 0` is exactly
/// [`voxelize_microstrip`]. Pair with `AperturePortSpec::k_lo = model.k_gnd`.
pub fn voxelize_microstrip_open(
    layout: &Layout,
    opts: &VoxelOptions,
    air_below_cells: usize,
) -> MicrostripModel {
    voxelize_with_ground_level(layout, opts, air_below_cells)
}

/// Voxelize a **finite board** (FS.2b.1, ADR-0207): the dielectric slab
/// AND the ground sheet extend only `board_margin_m` beyond the layout
/// bbox instead of filling the whole domain — a real PCB floating in
/// air, lifted `air_below_cells` above the domain floor. This is the
/// fixture absolute far-field products need: with the substrate bounded,
/// the NTFF equivalence box encloses the entire antenna and passes
/// through homogeneous air on every face (the measured motivation: the
/// whole-domain slab forced the box through dielectric and the patch
/// read a non-physical 22 dBi). Pair with all-six-face CPML and
/// `AperturePortSpec::k_lo = model.k_gnd`.
pub fn voxelize_finite_board(
    layout: &Layout,
    opts: &VoxelOptions,
    board_margin_m: f64,
    air_below_cells: usize,
) -> MicrostripModel {
    voxelize_inner(layout, opts, air_below_cells, Some(board_margin_m))
}

fn voxelize_with_ground_level(
    layout: &Layout,
    opts: &VoxelOptions,
    k_gnd: usize,
) -> MicrostripModel {
    voxelize_inner(layout, opts, k_gnd, None)
}

fn voxelize_inner(
    layout: &Layout,
    opts: &VoxelOptions,
    k_gnd: usize,
    board_margin_m: Option<f64>,
) -> MicrostripModel {
    let dx = opts.dx_m;
    assert!(
        dx.is_finite() && dx > 0.0,
        "VoxelOptions::dx_m must be positive and finite"
    );

    // --- X-Y extent: bbox padded by `margin` cells of air on every side. ---
    let margin = opts.xy_margin_cells as f64 * dx;
    let x0 = layout.bbox.min.x - margin;
    let x1 = layout.bbox.max.x + margin;
    let y0 = layout.bbox.min.y - margin;
    let y1 = layout.bbox.max.y + margin;

    let nx = ((x1 - x0) / dx).ceil() as usize;
    let ny = ((y1 - y0) / dx).ceil() as usize;
    assert!(
        nx > 0 && ny > 0,
        "voxelize_microstrip: degenerate x-y extent (nx={nx}, ny={ny}); \
         increase xy_margin_cells or check the layout bbox"
    );

    // --- Z-stack. ---
    //
    // The ground plane (tangential-E PEC) sits on the plane `z = 0` (the `k = 0`
    // staggered node). The trace (tangential-E PEC) sits on the plane
    // `z = k_top·dx`. The dielectric must fill the *entire* ground-to-trace gap,
    // which is `n_sub` cubic cells thick, so the trace is at `k_top = n_sub`
    // (giving a ground-to-trace spacing of exactly `n_sub·dx ≈ h`) and the
    // dielectric fills the `E_z` edges `k = 0 .. n_sub`.
    //
    // (Earlier the trace sat at `k_top = 1 + n_sub` with the dielectric filling
    // only `k = 1 ..= n_sub`, which left a one-cell *air* gap between the ground
    // and the substrate — a series air capacitance that drove the FDTD-measured
    // ε_eff ~20 % too low: the propagation gate fdtd-line-eeff-001 measured
    // ε_eff ≈ 2.5 vs analytic 3.33 until this gap was closed, after which it
    // measured 3.31, ADR-0108.)
    let n_sub = ((layout.substrate.height_m / dx).round() as usize).max(1);
    // Ground sheet at k_gnd (0 in the classic stack; air_below_cells in the
    // lifted stack) + n_sub dielectric cells -> top-metal layer.
    let k_top = k_gnd + n_sub;
    let nz = k_top + 1 + opts.air_above_cells;

    let eps_r_sub = layout.substrate.eps_r;

    // --- Material arrays at YeeGrid's exact required shapes. ---
    // `with_eps_r_cells` requires `(nx+1, ny+1, nz+1)`.
    // A horizontal PEC sheet (the ground plane and the metal traces) zeroes the
    // TANGENTIAL field — `Ex` and `Ey` — on its plane, NOT the normal `Ez`. So
    // mask the two in-plane components at their staggered node positions:
    //   `with_pec_mask_ex` requires `(nx, ny+1, nz+1)`; `Ex` node at ((i+0.5)dx, j·dx).
    //   `with_pec_mask_ey` requires `(nx+1, ny, nz+1)`; `Ey` node at (i·dx, (j+0.5)dx).
    let mut eps = Array3::<f64>::from_elem((nx + 1, ny + 1, nz + 1), 1.0);
    let mut pec_ex = Array3::<bool>::from_elem((nx, ny + 1, nz + 1), false);
    let mut pec_ey = Array3::<bool>::from_elem((nx + 1, ny, nz + 1), false);

    let in_trace = |x: f64, y: f64| {
        layout
            .traces
            .iter()
            .any(|p| point_in_polygon(Point2 { x, y }, p))
    };

    // Substrate dielectric on the `E_z` edges spanning the ground-to-trace gap,
    // `k = 0 .. n_sub` (i.e. `0 ..= k_top − 1`); air `ε_r = 1.0` elsewhere is the
    // array default. Filling from `k = 0` (the cell directly above the ground
    // plane) leaves NO air series gap at the ground — see the z-stack note above.
    // Finite board: dielectric + ground live only inside the board rect
    // (cell-centre / node-position tests against the physical bounds).
    let board = board_margin_m.map(|m| {
        (
            layout.bbox.min.x - m,
            layout.bbox.min.y - m,
            layout.bbox.max.x + m,
            layout.bbox.max.y + m,
        )
    });
    let on_board = |x: f64, y: f64| match board {
        None => true,
        Some((bx0, by0, bx1, by1)) => x >= bx0 && x <= bx1 && y >= by0 && y <= by1,
    };
    for i in 0..nx {
        for j in 0..ny {
            let (xc, yc) = (x0 + (i as f64 + 0.5) * dx, y0 + (j as f64 + 0.5) * dx);
            if !on_board(xc, yc) {
                continue;
            }
            for k in k_gnd..k_top {
                eps[(i, j, k)] = eps_r_sub;
            }
        }
    }

    // Tangential `Ex` PEC: ground plane (k = 0, whole layer) + traces
    // (k = k_top, where the `Ex` node ((i+0.5)dx, j·dx) lies under a trace).
    for i in 0..nx {
        for j in 0..=ny {
            if !on_board(x0 + (i as f64 + 0.5) * dx, y0 + j as f64 * dx) {
                continue;
            }
            pec_ex[(i, j, k_gnd)] = true;
            let x = x0 + (i as f64 + 0.5) * dx;
            let y = y0 + j as f64 * dx;
            if in_trace(x, y) {
                pec_ex[(i, j, k_top)] = true;
            }
        }
    }

    // Tangential `Ey` PEC: ground plane (k = 0) + traces (k = k_top), with the
    // `Ey` node at (i·dx, (j+0.5)dx).
    for i in 0..=nx {
        for j in 0..ny {
            if !on_board(x0 + i as f64 * dx, y0 + (j as f64 + 0.5) * dx) {
                continue;
            }
            pec_ey[(i, j, k_gnd)] = true;
            let x = x0 + i as f64 * dx;
            let y = y0 + (j as f64 + 0.5) * dx;
            if in_trace(x, y) {
                pec_ey[(i, j, k_top)] = true;
            }
        }
    }

    // --- Map layout ports to grid cells at the top-metal layer. ---
    let port_cells = layout
        .ports
        .iter()
        .map(|port| {
            let i = (((port.at.x - x0) / dx).floor() as isize).clamp(0, nx as isize - 1) as usize;
            let j = (((port.at.y - y0) / dx).floor() as isize).clamp(0, ny as isize - 1) as usize;
            (i, j, k_top)
        })
        .collect();

    let grid = YeeGrid::vacuum(nx, ny, nz, dx)
        .with_eps_r_cells(eps)
        .with_pec_mask_ex(pec_ex)
        .with_pec_mask_ey(pec_ey);

    MicrostripModel {
        k_gnd,
        grid,
        dims: (nx, ny, nz),
        dx_m: dx,
        port_cells,
    }
}

// ===========================================================================
// Graded (nonuniform) voxelization — FS.0b.1, ADR-0210
// ===========================================================================

/// Per-axis graded grid description for [`voxelize_microstrip_graded`]
/// (FS.0b.1, ADR-0210): primal cell widths (the `yee_compute` /
/// `yee_engine` `GradedSpacings` convention, ADR-0208) plus the domain
/// origin and the z-stack layer indices a rule generator
/// (`yee_engine::automesh::auto_spacings`) computed. The voxelizer
/// rasterizes against the **true** cell centres and node planes these
/// arrays imply — cumulative sums from `x0_m`/`y0_m`/`z = 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct GradedVoxelGrid {
    /// Primal cell widths along x (length `nx`, metres).
    pub dx_m: Vec<f64>,
    /// Primal cell widths along y (length `ny`, metres).
    pub dy_m: Vec<f64>,
    /// Primal cell widths along z (length `nz`, metres).
    pub dz_m: Vec<f64>,
    /// Layout-frame x of grid node 0 (typically `bbox.min.x − margin`).
    pub x0_m: f64,
    /// Layout-frame y of grid node 0.
    pub y0_m: f64,
    /// Ground-sheet layer index (`0` for the classic floor-ground stack).
    pub k_gnd: usize,
    /// Trace-plane layer index; the substrate fills cells
    /// `k = k_gnd .. k_top` (the ADR-0108 no-air-gap z-stack).
    pub k_top: usize,
}

/// A material-assigned graded model produced from a planar microstrip
/// layout by [`voxelize_microstrip_graded`].
///
/// Deliberately **not** a `yee_fdtd::YeeGrid`: that type carries a scalar
/// `dx` and a dt derived from it, both meaningless on a graded grid. The
/// raw arrays feed `yee_engine::MaterialsSpec` directly and the node
/// coordinates support probe/port placement by physical position.
#[derive(Debug, Clone)]
pub struct GradedMicrostripModel {
    /// Grid cell dimensions `(nx, ny, nz)`.
    pub dims: (usize, usize, usize),
    /// Per-cell relative permittivity, shape `[nx+1, ny+1, nz+1]`.
    pub eps_r_cells: Array3<f64>,
    /// Tangential `E_x` PEC mask (ground + traces), shape `[nx, ny+1, nz+1]`.
    pub pec_mask_ex: Array3<bool>,
    /// Tangential `E_y` PEC mask (ground + traces), shape `[nx+1, ny, nz+1]`.
    pub pec_mask_ey: Array3<bool>,
    /// Each layout port mapped to its `(i, j, k_top)` grid cell by
    /// coordinate lookup.
    pub port_cells: Vec<(usize, usize, usize)>,
    /// Ground-sheet layer index (echoes [`GradedVoxelGrid::k_gnd`]).
    pub k_gnd: usize,
    /// Trace-plane layer index (echoes [`GradedVoxelGrid::k_top`]).
    pub k_top: usize,
    /// Layout-frame x of every grid node plane (length `nx + 1`).
    pub x_nodes_m: Vec<f64>,
    /// Layout-frame y of every grid node plane (length `ny + 1`).
    pub y_nodes_m: Vec<f64>,
    /// z of every grid node plane, ground plane at the `k_gnd` node
    /// (length `nz + 1`, starting at `0.0`).
    pub z_nodes_m: Vec<f64>,
}

impl GradedMicrostripModel {
    /// The x-index of the primal cell containing layout-frame `x_m`
    /// (largest `i` with `x_nodes_m[i] ≤ x_m`, clamped to the grid) —
    /// the graded generalization of `floor((x − x0)/dx)`.
    pub fn cell_at_x(&self, x_m: f64) -> usize {
        cell_containing(&self.x_nodes_m, x_m)
    }

    /// The y-index of the primal cell containing layout-frame `y_m`.
    pub fn cell_at_y(&self, y_m: f64) -> usize {
        cell_containing(&self.y_nodes_m, y_m)
    }
}

/// Largest cell index `i` with `nodes[i] ≤ v`, clamped to `[0, n − 1]`
/// (`nodes` has length `n + 1`).
fn cell_containing(nodes: &[f64], v: f64) -> usize {
    let n_cells = nodes.len() - 1;
    nodes[..n_cells].partition_point(|&x| x <= v).max(1) - 1
}

/// Node coordinates (length `n + 1`) and cell centres (length `n`) from an
/// origin and primal widths.
///
/// Computed per **maximal run of identical widths** as
/// `origin + (base + m·d)` / `origin + (base + (m + 0.5)·d)`, so a constant
/// array reproduces the uniform voxelizer's `x0 + i·dx` and
/// `x0 + (i + 0.5)·dx` coordinates **bit-exactly** (the FS.0b.0
/// degenerate-by-construction discipline; a naive running sum differs by an
/// ulp, which flips point-in-polygon at trace edges that sit exactly on
/// node planes — measured, gate `voxel-graded-001`).
fn axis_coords(origin: f64, widths: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = widths.len();
    let mut nodes = Vec::with_capacity(n + 1);
    let mut centres = Vec::with_capacity(n);
    let mut base = 0.0_f64;
    nodes.push(origin + base);
    let mut s = 0;
    while s < n {
        let d = widths[s];
        let mut e = s + 1;
        while e < n && widths[e] == d {
            e += 1;
        }
        for m in 0..(e - s) {
            centres.push(origin + (base + (m as f64 + 0.5) * d));
            nodes.push(origin + (base + (m as f64 + 1.0) * d));
        }
        base += (e - s) as f64 * d;
        s = e;
    }
    (nodes, centres)
}

/// Voxelize a planar microstrip [`Layout`] onto a **graded** grid
/// (FS.0b.1, ADR-0210): the same rasterization as [`voxelize_microstrip`]
/// — PEC ground sheet at `k_gnd`, dielectric cells `k_gnd .. k_top`,
/// trace PEC at `k_top` where a polygon covers the cell centre, air above
/// — but against per-axis coordinate arrays instead of a scalar `dx`:
/// cell centres are `node[i] + d[i]/2`, `E_x`/`E_y` node positions use the
/// true staggered coordinates, and port cells are found by coordinate
/// lookup. With constant spacing arrays the produced masks are
/// bit-identical to the uniform voxelizer's (gate `voxel-graded-001`).
///
/// # Panics
///
/// Panics if any spacing array is empty or contains a non-positive /
/// non-finite width, or if `k_top` does not satisfy
/// `k_gnd < k_top < nz` (the trace plane needs an air layer above it).
pub fn voxelize_microstrip_graded(
    layout: &Layout,
    grid: &GradedVoxelGrid,
) -> GradedMicrostripModel {
    for (axis, arr) in [
        ("dx_m", &grid.dx_m),
        ("dy_m", &grid.dy_m),
        ("dz_m", &grid.dz_m),
    ] {
        assert!(!arr.is_empty(), "GradedVoxelGrid::{axis} must be non-empty");
        assert!(
            arr.iter().all(|d| d.is_finite() && *d > 0.0),
            "GradedVoxelGrid::{axis} widths must be positive and finite"
        );
    }
    let (nx, ny, nz) = (grid.dx_m.len(), grid.dy_m.len(), grid.dz_m.len());
    assert!(
        grid.k_gnd < grid.k_top && grid.k_top < nz,
        "GradedVoxelGrid: need k_gnd ({}) < k_top ({}) < nz ({nz})",
        grid.k_gnd,
        grid.k_top
    );

    let (x_nodes, xc) = axis_coords(grid.x0_m, &grid.dx_m);
    let (y_nodes, yc) = axis_coords(grid.y0_m, &grid.dy_m);
    let (z_nodes, _) = axis_coords(0.0, &grid.dz_m);

    let mut eps = Array3::<f64>::from_elem((nx + 1, ny + 1, nz + 1), 1.0);
    let mut pec_ex = Array3::<bool>::from_elem((nx, ny + 1, nz + 1), false);
    let mut pec_ey = Array3::<bool>::from_elem((nx + 1, ny, nz + 1), false);

    let in_trace = |x: f64, y: f64| {
        layout
            .traces
            .iter()
            .any(|p| point_in_polygon(Point2 { x, y }, p))
    };

    // Substrate dielectric on the E_z edges spanning the ground-to-trace
    // gap, `k = k_gnd .. k_top` — the ADR-0108 no-air-gap z-stack.
    let eps_r_sub = layout.substrate.eps_r;
    for i in 0..nx {
        for j in 0..ny {
            for k in grid.k_gnd..grid.k_top {
                eps[(i, j, k)] = eps_r_sub;
            }
        }
    }

    // Tangential Ex PEC: ground (whole k_gnd layer) + traces at k_top,
    // Ex node at (cell-centre x, node y).
    for i in 0..nx {
        for j in 0..=ny {
            pec_ex[(i, j, grid.k_gnd)] = true;
            if in_trace(xc[i], y_nodes[j]) {
                pec_ex[(i, j, grid.k_top)] = true;
            }
        }
    }

    // Tangential Ey PEC: ground + traces, Ey node at (node x, cell-centre y).
    for i in 0..=nx {
        for j in 0..ny {
            pec_ey[(i, j, grid.k_gnd)] = true;
            if in_trace(x_nodes[i], yc[j]) {
                pec_ey[(i, j, grid.k_top)] = true;
            }
        }
    }

    let port_cells = layout
        .ports
        .iter()
        .map(|port| {
            (
                cell_containing(&x_nodes, port.at.x),
                cell_containing(&y_nodes, port.at.y),
                grid.k_top,
            )
        })
        .collect();

    GradedMicrostripModel {
        dims: (nx, ny, nz),
        eps_r_cells: eps,
        pec_mask_ex: pec_ex,
        pec_mask_ey: pec_ey,
        port_cells,
        k_gnd: grid.k_gnd,
        k_top: grid.k_top,
        x_nodes_m: x_nodes,
        y_nodes_m: y_nodes,
        z_nodes_m: z_nodes,
    }
}

/// Surface resistance of a good conductor at frequency `f_hz` (R.0b,
/// ADR-0202): `R_s = sqrt(pi f mu0 / sigma)` ohms per square - the value
/// the resistive-sheet trace boundary (`MaterialsSpec::sheet_r_ohm`) takes
/// at the design frequency. Copper at 5 GHz: ~18.4 milli-ohm/sq.
pub fn surface_resistance_ohm(f_hz: f64, sigma_s_m: f64) -> f64 {
    const MU0: f64 = 4.0e-7 * std::f64::consts::PI;
    (std::f64::consts::PI * f_hz * MU0 / sigma_s_m).sqrt()
}

/// Map a substrate `tan δ` to per-cell conductivity for the engine's
/// lossy CA/CB update (R.0, ADR-0194): `σ = 2π f_ref ε₀ ε_r tan δ` on
/// every cell whose ε_r exceeds 1 (the substrate), zero in air. FDTD σ
/// is frequency-flat, so the map is exact at `f_ref` (the design
/// frequency) and the standard single-frequency approximation elsewhere.
/// Returns a flat `[nx+1, ny+1, nz+1]` row-major vector matching the
/// `MaterialsSpec::sigma_cells` convention.
///
/// # Panics
///
/// Panics if the model carries no ε_r map, or on non-physical inputs.
pub fn substrate_sigma_cells(model: &MicrostripModel, tan_d: f64, f_ref_hz: f64) -> Vec<f64> {
    assert!(
        tan_d >= 0.0 && tan_d.is_finite(),
        "tan_d must be finite and >= 0 (got {tan_d})"
    );
    assert!(
        f_ref_hz > 0.0 && f_ref_hz.is_finite(),
        "f_ref_hz must be positive and finite (got {f_ref_hz})"
    );
    const EPS0: f64 = 8.854_187_817e-12;
    let eps = model
        .grid
        .eps_r_cells
        .as_ref()
        .expect("substrate_sigma_cells: model has no eps_r map");
    let omega = std::f64::consts::TAU * f_ref_hz;
    eps.as_slice()
        .unwrap()
        .iter()
        .map(|&er| {
            if er > 1.0 {
                omega * EPS0 * er * tan_d
            } else {
                0.0
            }
        })
        .collect()
}

/// Map an N-layer [`Stackup`]'s per-layer `loss_tangent` to per-cell
/// conductivity (FS.4.2b, ADR-0226): σ = 2π f_ref ε₀ ε_r(layer) tan
/// δ(layer) for every `E_z` edge inside that layer's k-band, `0`
/// elsewhere (air, the lid plane, any plane above the stack). Re-derives
/// the k-bands from `stackup.layers` heights and `model.dx_m` with the
/// exact quantization [`voxelize_stackup`] used (`round(height/dx).max(1)`,
/// contiguous bands from `k = 0`) rather than inferring layer identity
/// from the ε value — two layers can share `ε_r`. Returns a flat
/// `[nx+1, ny+1, nz+1]` row-major vector matching the
/// `MaterialsSpec::sigma_cells` convention — the multilayer
/// generalization of [`substrate_sigma_cells`].
///
/// All-zero `tan δ` across every layer is a provable no-op: each band's
/// σ is `... * 0.0` exactly, so the returned vector is all-zero.
///
/// # Panics
///
/// Panics if `model` carries no ε_r map, if any layer height is not
/// positive, or on a non-physical (negative/non-finite) `loss_tangent`
/// or `f_ref_hz`.
pub fn stackup_sigma_cells(model: &MicrostripModel, stackup: &Stackup, f_ref_hz: f64) -> Vec<f64> {
    assert!(
        f_ref_hz > 0.0 && f_ref_hz.is_finite(),
        "f_ref_hz must be positive and finite (got {f_ref_hz})"
    );
    const EPS0: f64 = 8.854_187_817e-12;
    let omega = std::f64::consts::TAU * f_ref_hz;
    let dx = model.dx_m;

    // Same k-band bookkeeping voxelize_stackup used to fill ε: contiguous
    // cell bands from k = 0, height quantized against the model's dx.
    // sigma_by_k[k] is the layer sigma covering that cell; entries beyond
    // the stack (air, lid plane) are never resized in, so a plain
    // out-of-range lookup below naturally returns 0.
    let mut sigma_by_k = Vec::new();
    let mut k = 0usize;
    for layer in &stackup.layers {
        assert!(
            layer.height_m > 0.0,
            "stackup layer height must be positive"
        );
        assert!(
            layer.loss_tangent >= 0.0 && layer.loss_tangent.is_finite(),
            "stackup layer loss_tangent must be finite and >= 0 (got {})",
            layer.loss_tangent
        );
        let n_cells = ((layer.height_m / dx).round() as usize).max(1);
        let sigma = omega * EPS0 * layer.eps_r * layer.loss_tangent;
        sigma_by_k.resize(k + n_cells, sigma);
        k += n_cells;
    }

    let eps = model
        .grid
        .eps_r_cells
        .as_ref()
        .expect("stackup_sigma_cells: model has no eps_r map");
    // The dielectric fill (both voxelize_stackup and voxelize_microstrip)
    // only ever writes i in 0..nx, j in 0..ny — the eps array's `nx+1`/
    // `ny+1` shape carries one padding plane per axis (shared with the
    // other field components' node counts) that is never touched and
    // stays at the `1.0` (air) default. Mirror that exact footprint here
    // so a lookup at the padding plane reads 0, matching the real ε map
    // instead of extending the layer's σ past where its ε was ever set.
    let (nx, ny, _) = model.dims;
    eps.indexed_iter()
        .map(|((ii, jj, kk), _)| {
            if ii < nx && jj < ny {
                sigma_by_k.get(kk).copied().unwrap_or(0.0)
            } else {
                0.0
            }
        })
        .collect()
}

/// Punch a **via** through the substrate at grid column `(i, j)`
/// (R.1, ADR-0194): the `E_z` edges `k = 0 .. k_top` become PEC — a ground-to-trace
/// shorting post when `k_top` is the trace plane.
///
/// Equivalent to [`with_via_between`]`(model, i, j, 0, k_top)` (FS.4.1,
/// ADR-0221 kept it as the single-layer back-compat spelling).
///
/// # Panics
///
/// Panics if `(i, j)` or `k_top` exceeds the `E_z` grid.
pub fn with_via_at_cell(model: &mut MicrostripModel, i: usize, j: usize, k_top: usize) {
    with_via_between(model, i, j, 0, k_top);
}

/// Punch a **blind via** at grid column `(i, j)` (FS.4.1, ADR-0221): the
/// `E_z` edges `k = k_lo .. k_hi` become PEC — a metal post spanning the
/// z node-plane `k_lo` to the z node-plane `k_hi`, touching neither
/// outer plate unless the range does.
///
/// `k_lo`/`k_hi` are **grid cell indices**, not stackup layer indices:
/// callers on a [`voxelize_stackup`] model quantize layer heights the
/// same way the voxelizer did (`round(h/dx).max(1)` cumulative bands)
/// and follow the crate's cell-index post-processing idiom
/// ([`with_via_at_cell`] / [`truncate_ground_at_cell`]).
///
/// # Panics
///
/// Panics if `(i, j)` exceeds the `E_z` grid or the range is malformed
/// (`k_lo > k_hi` or `k_hi > nz`) — a malformed via is a caller bug,
/// not data.
pub fn with_via_between(model: &mut MicrostripModel, i: usize, j: usize, k_lo: usize, k_hi: usize) {
    let (nx, ny, nz) = model.dims;
    assert!(
        i <= nx && j <= ny && k_lo <= k_hi && k_hi <= nz,
        "via at ({i}, {j}) spanning k {k_lo}..{k_hi} outside grid {nx}x{ny}x{nz}"
    );
    let mut mask = model
        .grid
        .pec_mask_ez
        .take()
        .unwrap_or_else(|| ndarray::Array3::from_elem((nx + 1, ny + 1, nz), false));
    for k in k_lo..k_hi {
        mask[(i, j, k)] = true;
    }
    model.grid.pec_mask_ez = Some(mask);
}

/// Punch a **through-via** at grid column `(i, j)` (FS.4.1, ADR-0221):
/// the full-stack case of [`with_via_between`] — every `E_z` edge
/// `k = 0 .. nz` becomes PEC, connecting the ground plane to the domain
/// top. On a lidded [`voxelize_stackup`] model the node plane `nz` *is*
/// the lid, so this is a ground-to-lid barrel through every layer (and
/// through the trace plane, shorting any trace it passes).
pub fn with_through_via_at_cell(model: &mut MicrostripModel, i: usize, j: usize) {
    let nz = model.dims.2;
    with_via_between(model, i, j, 0, nz);
}

/// Truncate the ground plane to the strip of the first `i_ground_end` grid
/// columns (FS.1a.0, ADR-0205): the `k = 0` tangential-E PEC then covers
/// exactly `x ∈ [grid x₀, grid x₀ + i_ground_end·dx]` and the rest of the
/// bottom layer is open substrate. This is the quasi-Yagi's defining
/// structural feature — the truncated ground *is* the reflector element —
/// and follows the [`with_via_at_cell`] post-processing idiom (cell-index
/// frame; callers map layout-frame x through the same
/// `x₀ = bbox.min.x − margin` origin the voxelizer used).
///
/// Edge rule, documented so the boundary is exact: the ground edge sits on
/// the plane `x = x₀ + i_ground_end·dx`. `Ex` nodes (centred at
/// `(i + 0.5)·dx`) stay PEC for `i < i_ground_end`; `Ey` nodes (at `i·dx`,
/// on-or-inside the edge) stay PEC for `i ≤ i_ground_end`. Trace masks at
/// `k_top` are untouched. `i_ground_end ≥ nx` is a no-op (full ground) so
/// callers can pass a clamped layout coordinate without special-casing.
///
/// # Panics
///
/// Panics if the model was built without the `Ex`/`Ey` PEC masks (every
/// [`voxelize_microstrip`] model has them).
pub fn truncate_ground_at_cell(model: &mut MicrostripModel, i_ground_end: usize) {
    let (nx, ny, _) = model.dims;
    if i_ground_end >= nx {
        return;
    }
    let kg = model.k_gnd;
    let ex = model
        .grid
        .pec_mask_ex
        .as_mut()
        .expect("truncate_ground_at_cell: model has no Ex PEC mask");
    for i in i_ground_end..nx {
        for j in 0..=ny {
            ex[(i, j, kg)] = false;
        }
    }
    let ey = model
        .grid
        .pec_mask_ey
        .as_mut()
        .expect("truncate_ground_at_cell: model has no Ey PEC mask");
    for i in (i_ground_end + 1)..=nx {
        for j in 0..ny {
            ey[(i, j, kg)] = false;
        }
    }
}

/// Test whether point `p` lies inside polygon `poly` via the standard
/// even-odd ray-cast (crossing-number) rule.
///
/// A horizontal ray is cast to `+x`; the point is inside when it crosses an odd
/// number of polygon edges. The polygon is treated as implicitly closed
/// (last vertex → first). Robust for the axis-aligned rectangular traces this
/// crate targets; points exactly on an edge are reported consistently (no
/// special handling) which is acceptable for cell-centre sampling.
fn point_in_polygon(p: Point2, poly: &Polygon) -> bool {
    let verts = &poly.verts;
    let n = verts.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = verts[i];
        let vj = verts[j];
        // Does edge (vj -> vi) straddle the horizontal line y = p.y, and is the
        // edge's x at y = p.y to the right of p.x?
        if (vi.y > p.y) != (vj.y > p.y) {
            let x_cross = (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y) + vi.x;
            if p.x < x_cross {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

// ===========================================================================
// FDTD propagation-based microstrip ε_eff driver (Filter Phase F1.1b.1,
// ADR-0108)
// ===========================================================================
//
// `run_line_eeff` is the SHIPPED full-wave gate for F1.1b.1. It supersedes the
// resonant-split method above, which PR #1 (7 CI/local iterations) proved
// unworkable for a microstrip on an *open* domain:
//
//   * a small hard-PEC box CONFINES the fringing / air-gap fields that set the
//     even/odd ε_eff split → the split collapses (k ≈ 0.02);
//   * a large hard-PEC box becomes a resonant CAVITY whose box modes swamp the
//     microstrip resonances → argmax picks box modes (k ≈ 0.01);
//   * open CPML walls remove both pathologies but KILL the resonator Q (the
//     λ/2 line radiates into the absorber) → zero detectable peaks.
//
// There is no box that is simultaneously high-Q and non-confining/
// non-resonant, so a *resonant* split is the wrong observable. The robust,
// textbook FDTD coupled-line characterization is a PROPAGATION measurement:
// drive a long, NON-resonant line, terminate both propagation ends in matched
// (CPML) loads, and read the phase velocity of the traveling wave directly off
// two probe planes a known distance apart. No Q, no cavity modes, no
// peak-picking.
//
// `run_line_eeff` measures the single-line ε_eff this way (the validated
// walking-skeleton gate `fdtd-line-eeff-001`). The same machinery generalizes
// to the even/odd coupled split via `run_coupled_line_eeff` (drive two coupled
// strips in-phase → even, anti-phase → odd; gate `fdtd-line-eeff-coupled-001`).

/// Configuration for the propagation-based ε_eff drivers [`run_line_eeff`] and
/// [`run_coupled_line_eeff`].
///
/// The defaults target a single straight FR-4 microstrip line driven at 5 GHz
/// (a short guided wavelength → fewer cells per λ_g → a cheap grid) with a
/// modulated-Gaussian lumped port and a long line terminated in stable hard-PEC
/// walls (the DFT is time-gated to the forward pulse, before the far-wall
/// reflection returns — see [`Self::gate_steps`]). Two probe planes a known
/// distance apart sample `E_z`; the single-bin DFT phase advance between them at
/// the drive centre frequency gives the phase velocity → `ε_eff`. The defaults
/// are tuned to converge the `fdtd-line-eeff-001` gate within ≤ 15 % on a
/// tractable grid (validated locally in the bounded dev container, ADR-0108).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineRunConfig {
    /// Isotropic Yee cell size `dx = dy = dz`, metres.
    pub dx_m: f64,
    /// Air margin (in cells) around the layout bounding box in `x`/`y`
    /// (forwarded to [`VoxelOptions::xy_margin_cells`]). Must keep the strip a
    /// few cells clear of the CPML region.
    pub xy_margin_cells: usize,
    /// Air layers above the top metal (forwarded to
    /// [`VoxelOptions::air_above_cells`]).
    pub air_above_cells: usize,
    /// Drive centre frequency `f0` (Hz). The single-bin DFT phase is read at
    /// this frequency; the guided wavelength `λ_g = c / (f0·√ε_eff)` sets the
    /// probe spacing in wavelengths.
    pub f0_hz: f64,
    /// Fractional FWHM bandwidth of the modulated-Gaussian drive, as a fraction
    /// of `f0` (`bw = freq_span · f0`). Broadband enough to launch a clean
    /// traveling pulse but narrow enough that `f0` dominates the spectrum at the
    /// probes.
    pub freq_span: f64,
    /// Total number of FDTD time steps. Must be long enough for the launched
    /// pulse to fully transit both probe planes and clear the far CPML, so the
    /// single-bin DFT at each probe integrates the whole passage.
    pub n_steps: usize,
    /// Series resistance (Ω) of the drive port. A moderate value (≈ the line
    /// impedance) launches the wave without a strong reflection at the feed.
    pub port_resistance_ohm: f64,
    /// Peak drive voltage (V) of the Gaussian-modulated pulse.
    pub drive_v0: f64,
    /// Number of time steps over which the single-bin DFT at each probe is
    /// **integrated**, counted from `t = 0`. This *time-gates* the measurement
    /// to the FORWARD traveling pulse, before the reflection off the far line
    /// end returns to the probes.
    ///
    /// The line is terminated in hard-PEC outer walls (stable — see below), so
    /// the far end reflects. A long enough line (several λ_g of clearance past
    /// the downstream probe) keeps that reflection out of the gate window: the
    /// forward pulse fully transits both probes, and the DFT integrates only it.
    /// `None` integrates the whole `n_steps` record (use only when the line is
    /// long enough that no reflection ever returns within `n_steps`).
    ///
    /// **Why PEC walls + a time gate, not CPML:** the determined PR #1 method
    /// was open CPML terminations, but `CpmlParams::for_grid` applies CPML on
    /// all six faces and is **late-time unstable** for a microstrip whose PEC
    /// ground plane and high-ε substrate run *into* the boundary region
    /// (container measurement, ADR-0108: fields diverge ~1e13 by ~2.5 k steps).
    /// A PEC box is unconditionally stable; gating the DFT to the forward
    /// passage recovers a clean, reflection-free traveling-wave phase — the
    /// standard "long line, short look" FDTD line characterization.
    pub gate_steps: Option<usize>,
    /// Absolute `x` position (metres, in layout coordinates) of probe plane A —
    /// the *upstream* probe. Placed a fraction of a wavelength past the feed so
    /// the launch transient has settled into a clean traveling wave.
    pub probe_a_x_m: f64,
    /// Absolute `x` position (metres, in layout coordinates) of probe plane B —
    /// the *downstream* probe. `probe_b_x_m > probe_a_x_m`; the phase advance is
    /// read between A and B, so the spacing `probe_b_x_m − probe_a_x_m` should
    /// be ≈ λ_g/4 … λ_g/2 — large enough to resolve cleanly, small enough that
    /// the true advance is `< 2π` (unambiguous; no phase wrap).
    pub probe_b_x_m: f64,
}

impl Default for LineRunConfig {
    /// Walking-skeleton defaults (single FR-4 line at 5 GHz, cheap grid). See
    /// the per-field docs; container-validated to pass `fdtd-line-eeff-001`
    /// within ≤ 15 % (ADR-0108).
    fn default() -> Self {
        Self {
            // 0.4 mm cells: ~4 substrate cells through a 1.6 mm FR-4 slab and
            // ~80 cells per guided wavelength at 5 GHz (λ_g ≈ 33 mm). Fine
            // enough to resolve the phase advance between the two probe planes.
            dx_m: 0.4e-3,
            xy_margin_cells: 14,
            air_above_cells: 16,
            // 5 GHz: short λ_g → fewer cells per wavelength → a cheaper grid for
            // a given probe spacing (and well within FR-4's quasi-static band).
            f0_hz: 5.0e9,
            // 80 % FWHM: a SHORT launched pulse (≈ 810 steps) so it fully
            // transits both probes before the far-wall reflection returns —
            // wide enough to fit inside the time gate, narrow enough that f0
            // still dominates the probe spectra.
            freq_span: 0.8,
            // The forward pulse reaches the downstream probe within a few
            // hundred steps and is ~800 steps long; 2.5 k steps capture its full
            // passage well before the gate cutoff.
            n_steps: 2_500,
            // ≈ line impedance: launch the wave with a weak feed reflection.
            port_resistance_ohm: 50.0,
            drive_v0: 1.0,
            // Time-gate the DFT to the forward passage (the test sets this to
            // just before the far-PEC reflection returns to the downstream
            // probe — see `gate_steps` docs). `None` here; the test always sets
            // it explicitly to match its line length.
            gate_steps: None,
            // Probe planes in absolute layout-x metres. The default geometry is
            // a single ~6 λ_g line (λ_g ≈ 33 mm at 5 GHz); place the probes near
            // the middle, a third of a wavelength apart → Δx ≈ λ_g/3, a ~120°
            // advance that is unambiguous (< 2π) and well-resolved. The test
            // overrides these to match its line length.
            probe_a_x_m: 82.0e-3,
            probe_b_x_m: 93.0e-3,
        }
    }
}

/// Result of a propagation-based ε_eff measurement ([`run_line_eeff`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineRunResult {
    /// FDTD-extracted effective relative permittivity
    /// `ε_eff = (c / v_p)²`, from the measured phase velocity.
    pub eps_eff: f64,
    /// Measured phase velocity `v_p = ω·Δx / Δφ` (m/s).
    pub v_p: f64,
    /// Unwrapped phase advance `Δφ = φ_A − φ_B` between the two probe planes at
    /// the drive centre frequency (radians; positive for a `+x`-traveling
    /// wave).
    pub delta_phi: f64,
    /// Probe-plane separation `Δx` along the line (metres).
    pub delta_x: f64,
}

/// Run a full-wave FDTD solve on a single straight microstrip line and extract
/// its effective relative permittivity `ε_eff` from the measured **phase
/// velocity** of a traveling wave.
///
/// This is the SHIPPED walking-skeleton gate for Filter Phase F1.1b.1
/// (`fdtd-line-eeff-001`, ADR-0108): the first full-wave EM solve in the
/// filter-design pipeline. It composes shipped primitives —
/// [`voxelize_microstrip`] (layout→grid) and `yee-fdtd`'s [`LumpedRlcPort`]
/// drive + [`WalkingSkeletonSolver`] time-stepping.
///
/// # Method
///
/// A *non-resonant* propagation measurement (the textbook FDTD line
/// characterization; replaces the unworkable resonant-split method — see the
/// module comment and PR #1):
///
/// 1. Voxelize the line with hard-PEC outer walls (unconditionally stable; CPML
///    is late-time unstable for a microstrip whose PEC ground / high-ε substrate
///    run into the boundary — ADR-0108). The line is made several λ_g long so
///    the far-wall reflection returns *after* the forward pulse has cleared the
///    probes.
/// 2. Drive `E_z` at one end with a modulated-Gaussian lumped port, launching a
///    `+x`-traveling pulse.
/// 3. Record `E_z` in the substrate under the strip centre at two probe planes A
///    and B a known distance `Δx` apart along the line.
/// 4. Time-gated single-bin DFT at `f0` at each probe (gated to the forward
///    passage via [`LineRunConfig::gate_steps`], before the reflection returns)
///    → complex phasors; the phase advance `Δφ = φ_A − φ_B` (positive,
///    downstream lags) gives the phase velocity `v_p = ω·Δx / Δφ` and hence
///    `ε_eff = (c / v_p)²`.
///
/// # Panics
///
/// Panics if the layout has no ports (a drive port is required), if the probe
/// positions do not satisfy `probe_b_x_m > probe_a_x_m`, or if the resulting
/// probe planes collapse to the same grid column.
pub fn run_line_eeff(layout: &Layout, cfg: &LineRunConfig) -> LineRunResult {
    assert!(
        !layout.ports.is_empty(),
        "run_line_eeff: need ≥ 1 port (a drive port)"
    );
    assert!(
        cfg.probe_b_x_m > cfg.probe_a_x_m,
        "run_line_eeff: require probe_b_x_m ({}) > probe_a_x_m ({})",
        cfg.probe_b_x_m,
        cfg.probe_a_x_m
    );

    // --- 1. Voxelize: layout -> material-assigned grid + port cells. --------
    let opts = VoxelOptions {
        dx_m: cfg.dx_m,
        xy_margin_cells: cfg.xy_margin_cells,
        air_above_cells: cfg.air_above_cells,
    };
    let model = voxelize_microstrip(layout, &opts);
    let (nx, _ny, _nz) = model.dims;
    let drive_cell = model.port_cells[0];
    let (_i_drive, j_strip, k_top) = drive_cell;
    let dt = model.grid.dt;
    let dx = model.dx_m;

    // Probe `E_z` in the SUBSTRATE, under the strip — that is where the
    // quasi-TEM mode's dominant vertical field lives (between the trace at
    // `k_top` and the ground at `k = 0`). The `E_z` node at `(i, j, k)` spans
    // `z ∈ [k·dx, (k+1)·dx]`, so the substrate cell directly beneath the trace
    // is `k_top − 1`. (The node at `k_top` itself sits in the air just above
    // the metal, where the field is far weaker.)
    let k_probe = k_top.saturating_sub(1).max(1);

    // Map the two probe planes from absolute layout-x to grid columns. The
    // voxelizer's grid origin is `x0 = bbox.min.x − xy_margin_cells·dx` (see
    // `voxelize_microstrip`), so the column for layout-x `xp` is
    // `round((xp − x0)/dx)`. Sample `E_z` in the substrate under the strip
    // centre column `j_strip` at `k_probe`. The phase advance is read between
    // these two interior columns; the probe spacing `Δx = probe_b_x_m −
    // probe_a_x_m` is chosen (by the caller) < λ_g so the advance is
    // unambiguous.
    let x0 = layout.bbox.min.x - cfg.xy_margin_cells as f64 * dx;
    let i_for = |xp: f64| -> usize {
        (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize
    };
    let i_a = i_for(cfg.probe_a_x_m);
    let i_b = i_for(cfg.probe_b_x_m);
    assert!(
        i_b > i_a,
        "run_line_eeff: probe planes collapsed to the same column (i_a = i_b = {i_a}); \
         widen the probe spacing or refine dx"
    );
    let probe_a = (i_a, j_strip, k_probe);
    let probe_b = (i_b, j_strip, k_probe);
    let delta_x = (i_b - i_a) as f64 * dx;

    // --- 2. Hard-PEC outer walls (unconditionally stable). The CPML path was
    //        found late-time unstable for a microstrip whose PEC ground / high-ε
    //        substrate run into the boundary (ADR-0108); a PEC box plus a
    //        time-gated DFT (`gate_steps`) on the FORWARD pulse recovers a
    //        clean, reflection-free traveling-wave phase instead. --------------
    let mut solver = WalkingSkeletonSolver::new(model.grid);

    // --- 3. Drive port: modulated-Gaussian launch at the −x end. ------------
    let bw = cfg.freq_span * cfg.f0_hz;
    // Centre the pulse ~3.5 time-constants in so its t = 0 tail is negligible.
    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw))
        / dt)
        .ceil() as usize;
    let wave = SourceWaveform::GaussianPulse {
        v0: cfg.drive_v0,
        f0: cfg.f0_hz,
        bw,
        t0_steps,
    };
    let mut port = LumpedRlcPort::pure_resistor(drive_cell, cfg.port_resistance_ohm, wave);

    // --- 4. Time-step; integrate the single-bin DFT at each probe over the
    //        time gate (the forward passage, before the far-wall reflection
    //        returns). The step body matches the canonical `step_with_source`
    //        order: H + boundary, E + boundary, then the lumped-port correction
    //        (after the E boundary, per `LumpedRlcPort::correct_e`'s call site),
    //        then advance the clock. -------------------------------------------
    let mut acc = [0.0_f64; 4]; // [reA, imA, reB, imB] single-bin DFT
    let omega = 2.0 * std::f64::consts::PI * cfg.f0_hz;
    let gate = cfg.gate_steps.unwrap_or(cfg.n_steps).min(cfg.n_steps);
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        solver.apply_cpml_e();
        port.correct_e(solver.grid_mut(), n, dt);

        solver.advance_clock();

        if n < gate {
            let grid = solver.grid();
            let ez_a = grid.ez[probe_a];
            let ez_b = grid.ez[probe_b];
            // Single-bin DFT accumulation at f0. The gate confines the record to
            // the forward traveling pulse (finite support, fully contained), so
            // the rectangular window introduces no high-Q sidelobe comb.
            let phase = omega * n as f64 * dt;
            let (c, s) = (phase.cos(), phase.sin());
            acc[0] += ez_a * c;
            acc[1] -= ez_a * s;
            acc[2] += ez_b * c;
            acc[3] -= ez_b * s;
        }
    }

    // --- 5. Phase advance A → B → phase velocity → ε_eff. -------------------
    let phi_a = acc[1].atan2(acc[0]);
    let phi_b = acc[3].atan2(acc[2]);
    // A +x-traveling wave lags downstream, so φ decreases from A to B; the
    // phase advance Δφ = φ_A − φ_B is positive. Reduce into (0, 2π) so a wrap
    // across the atan2 branch cut does not flip the sign (the probe spacing is
    // chosen < λ_g so the true advance is < 2π).
    let mut delta_phi = phi_a - phi_b;
    while delta_phi <= 0.0 {
        delta_phi += 2.0 * std::f64::consts::PI;
    }
    while delta_phi > 2.0 * std::f64::consts::PI {
        delta_phi -= 2.0 * std::f64::consts::PI;
    }

    let v_p = omega * delta_x / delta_phi;
    let eps_eff = (C0_M_S / v_p).powi(2);

    let mag_a = (acc[0] * acc[0] + acc[1] * acc[1]).sqrt();
    let mag_b = (acc[2] * acc[2] + acc[3] * acc[3]).sqrt();
    eprintln!(
        "[fdtd-line-eeff DIAG] probes i_a={i_a} i_b={i_b} (j={j_strip}, k={k_probe}), \
         Δx={:.3} mm | φ_A={:.4} (|A|={:.3e}) φ_B={:.4} (|B|={:.3e}) | \
         Δφ={:.4} rad, v_p={:.4e} m/s, ε_eff={:.4}",
        delta_x * 1e3,
        phi_a,
        mag_a,
        phi_b,
        mag_b,
        delta_phi,
        v_p,
        eps_eff,
    );

    LineRunResult {
        eps_eff,
        v_p,
        delta_phi,
        delta_x,
    }
}

/// Speed of light in vacuum (m/s), for the `ε_eff = (c / v_p)²` conversion.
const C0_M_S: f64 = 299_792_458.0;

/// Result of a coupled-line even/odd ε_eff measurement
/// ([`run_coupled_line_eeff`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupledLineResult {
    /// Even-mode (in-phase drive) effective permittivity `ε_eff,e`.
    pub eps_eff_e: f64,
    /// Odd-mode (anti-phase drive) effective permittivity `ε_eff,o`.
    pub eps_eff_o: f64,
    /// The even/odd ε_eff-split `k = (ε_eff,e − ε_eff,o)/(ε_eff,e + ε_eff,o)`.
    /// (The even mode concentrates more field in the substrate, so normally
    /// `ε_eff,e ≥ ε_eff,o` and `k ≥ 0`.)
    pub k_split: f64,
}

/// Run two full-wave FDTD propagation solves on a *coupled* microstrip pair and
/// extract the even- and odd-mode effective permittivities `ε_eff,e` /
/// `ε_eff,o` from their phase velocities (Filter Phase F1.1b.1 coupled
/// follow-on, gate `fdtd-line-eeff-coupled-001`, ADR-0108).
///
/// # Method
///
/// The same reflection-free propagation measurement as [`run_line_eeff`],
/// applied to two parallel edge-coupled strips, run twice:
///
/// 1. **Even mode** — drive *both* strips in phase (`+v0`, `+v0`). The symmetric
///    excitation launches only the even supermode; the parity-matched probe
///    `E_z(strip1) + E_z(strip2)` reconstructs it. Its phase velocity →
///    `ε_eff,e`.
/// 2. **Odd mode** — drive the strips anti-phase (`+v0`, `−v0`). The
///    antisymmetric excitation launches only the odd supermode; the
///    parity-matched probe `E_z(strip1) − E_z(strip2)` reconstructs it (a SUM
///    would cancel). Its phase velocity → `ε_eff,o`.
///
/// Each run is the stable PEC-box + time-gated-DFT propagation measurement of
/// [`run_line_eeff`]; the split `k = (ε_eff,e − ε_eff,o)/(ε_eff,e + ε_eff,o)`
/// is the physically-correct even/odd ε_eff-split (NOT the impedance coupling —
/// see PR #1 root cause, ADR-0108).
///
/// # Panics
///
/// Panics if the layout has fewer than two ports (one drive per strip), or for
/// the same probe-geometry preconditions as [`run_line_eeff`].
pub fn run_coupled_line_eeff(
    layout: &Layout,
    cfg: &LineRunConfig,
    probe_b_offset_m: f64,
) -> CoupledLineResult {
    assert!(
        layout.ports.len() >= 2,
        "run_coupled_line_eeff: need ≥ 2 ports (one drive per strip); got {}",
        layout.ports.len()
    );
    let eps_eff_e = coupled_mode_eeff(layout, cfg, probe_b_offset_m, false);
    let eps_eff_o = coupled_mode_eeff(layout, cfg, probe_b_offset_m, true);
    // Order so the higher-ε (even) mode is first; the split is then ≥ 0.
    let (eps_eff_e, eps_eff_o) = if eps_eff_e >= eps_eff_o {
        (eps_eff_e, eps_eff_o)
    } else {
        (eps_eff_o, eps_eff_e)
    };
    let k_split = (eps_eff_e - eps_eff_o) / (eps_eff_e + eps_eff_o);
    CoupledLineResult {
        eps_eff_e,
        eps_eff_o,
        k_split,
    }
}

/// Drive a coupled microstrip pair into one supermode (even if `anti_phase ==
/// false`, odd if `true`) and return that mode's propagation ε_eff via the
/// time-gated phase-velocity measurement.
///
/// `probe_b_offset_m` is the downstream-probe offset from probe A (≈ λ_g/3); the
/// upstream probe A sits at `cfg.probe_a_x_m`. Both strips are probed in the
/// substrate (`k_top − 1`); the parity-matched combination
/// `E_z(strip1) ± E_z(strip2)` isolates the supermode.
fn coupled_mode_eeff(
    layout: &Layout,
    cfg: &LineRunConfig,
    probe_b_offset_m: f64,
    anti_phase: bool,
) -> f64 {
    let opts = VoxelOptions {
        dx_m: cfg.dx_m,
        xy_margin_cells: cfg.xy_margin_cells,
        air_above_cells: cfg.air_above_cells,
    };
    let model = voxelize_microstrip(layout, &opts);
    let (nx, _ny, _nz) = model.dims;
    let cell0 = model.port_cells[0];
    let cell1 = model.port_cells[1];
    let (_i0, j0, k_top) = cell0;
    let (_i1, j1, _k1) = cell1;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - cfg.xy_margin_cells as f64 * dx;
    let i_for = |xp: f64| -> usize {
        (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize
    };
    let i_a = i_for(cfg.probe_a_x_m);
    let i_b = i_for(cfg.probe_a_x_m + probe_b_offset_m);
    assert!(
        i_b > i_a,
        "run_coupled_line_eeff: probe planes collapsed (i_a = i_b = {i_a})"
    );
    let delta_x = (i_b - i_a) as f64 * dx;

    // Stable PEC outer walls; the lateral PEC walls are far from the strips.
    let mut solver = WalkingSkeletonSolver::new(model.grid);

    let bw = cfg.freq_span * cfg.f0_hz;
    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw))
        / dt)
        .ceil() as usize;
    let make_wave = |v0: f64| SourceWaveform::GaussianPulse {
        v0,
        f0: cfg.f0_hz,
        bw,
        t0_steps,
    };
    let v1 = if anti_phase {
        -cfg.drive_v0
    } else {
        cfg.drive_v0
    };
    let mut port0 =
        LumpedRlcPort::pure_resistor(cell0, cfg.port_resistance_ohm, make_wave(cfg.drive_v0));
    let mut port1 = LumpedRlcPort::pure_resistor(cell1, cfg.port_resistance_ohm, make_wave(v1));

    // Probe both strips at each plane; the parity-matched combination
    // (strip0 + strip1 for even, strip0 − strip1 for odd) reconstructs the
    // supermode, mirroring the drive parity (a mismatched combination would
    // cancel the mode of interest).
    let probe0_a = (i_a, j0, k_probe);
    let probe1_a = (i_a, j1, k_probe);
    let probe0_b = (i_b, j0, k_probe);
    let probe1_b = (i_b, j1, k_probe);

    let mut acc = [0.0_f64; 4]; // [reA, imA, reB, imB]
    let omega = 2.0 * std::f64::consts::PI * cfg.f0_hz;
    let gate = cfg.gate_steps.unwrap_or(cfg.n_steps).min(cfg.n_steps);
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        solver.apply_cpml_e();
        port0.correct_e(solver.grid_mut(), n, dt);
        port1.correct_e(solver.grid_mut(), n, dt);

        solver.advance_clock();

        if n < gate {
            let grid = solver.grid();
            let comb = |p0: (usize, usize, usize), p1: (usize, usize, usize)| {
                if anti_phase {
                    grid.ez[p0] - grid.ez[p1]
                } else {
                    grid.ez[p0] + grid.ez[p1]
                }
            };
            let ez_a = comb(probe0_a, probe1_a);
            let ez_b = comb(probe0_b, probe1_b);
            let phase = omega * n as f64 * dt;
            let (c, s) = (phase.cos(), phase.sin());
            acc[0] += ez_a * c;
            acc[1] -= ez_a * s;
            acc[2] += ez_b * c;
            acc[3] -= ez_b * s;
        }
    }

    let phi_a = acc[1].atan2(acc[0]);
    let phi_b = acc[3].atan2(acc[2]);
    let mut delta_phi = phi_a - phi_b;
    while delta_phi <= 0.0 {
        delta_phi += 2.0 * std::f64::consts::PI;
    }
    while delta_phi > 2.0 * std::f64::consts::PI {
        delta_phi -= 2.0 * std::f64::consts::PI;
    }
    let v_p = omega * delta_x / delta_phi;
    let eps_eff = (C0_M_S / v_p).powi(2);

    let mode = if anti_phase { "odd " } else { "even" };
    eprintln!(
        "[fdtd-line-eeff-coupled DIAG] {mode}: Δx={:.3} mm, Δφ={:.4} rad, \
         v_p={:.4e} m/s, ε_eff={:.4}",
        delta_x * 1e3,
        delta_phi,
        v_p,
        eps_eff,
    );
    eps_eff
}

// ===========================================================================
// Multilayer stackup voxelization — FS.4.0, ADR-0215
// ===========================================================================

/// Voxelize a planar layout on an N-layer [`Stackup`] (FS.4.0,
/// ADR-0215): ground plane at `k = 0`, `stackup.layers` filled bottom-up
/// with **no air gap anywhere in the stack** (the ADR-0108 lesson
/// generalized — an accidental series gap between layers poisons every
/// downstream result), the trace PEC on the plane at the TOP of
/// `layers[trace_layer]`, and — when `stackup.lid` — a whole-plane PEC
/// lid directly above the last layer (`air_above_cells` then ignored:
/// the domain ends at the lid).
///
/// Each layer quantizes to `round(height/dx).max(1)` cells (height error
/// ≤ dx/2, the walking-skeleton trade recorded in the ADR). A
/// single-layer open stackup with `trace_layer = 0` reproduces
/// [`voxelize_microstrip`] bit-identically (gate `voxel-stackup-001`).
///
/// Returns the standard [`MicrostripModel`] (ports mapped to the trace
/// plane) so every downstream fixture works unchanged.
///
/// # Panics
///
/// Panics if `trace_layer` is out of range or a layer height is not
/// positive — malformed stackups are caller bugs, not data.
pub fn voxelize_stackup(
    layout: &Layout,
    stackup: &Stackup,
    trace_layer: usize,
    opts: &VoxelOptions,
) -> MicrostripModel {
    let dx = opts.dx_m;
    assert!(
        dx.is_finite() && dx > 0.0,
        "VoxelOptions::dx_m must be positive and finite"
    );
    assert!(
        trace_layer < stackup.layers.len(),
        "trace_layer {trace_layer} out of range ({} layers)",
        stackup.layers.len()
    );

    // --- Layer quantization: contiguous cell bands, no gaps. ---
    let mut k_starts = Vec::with_capacity(stackup.layers.len());
    let mut k = 0usize;
    for layer in &stackup.layers {
        assert!(
            layer.height_m > 0.0,
            "stackup layer height must be positive"
        );
        k_starts.push(k);
        k += ((layer.height_m / dx).round() as usize).max(1);
    }
    let n_stack = k;
    let k_trace = k_starts[trace_layer]
        + ((stackup.layers[trace_layer].height_m / dx).round() as usize).max(1);
    // Open top: one guaranteed cell layer above the stack (mirroring the
    // microstrip voxelizer's top-metal layer, which makes the
    // single-layer case bit-identical) plus the air margin. Lidded: the
    // domain ends exactly at the lid plane.
    let nz = if stackup.lid {
        n_stack
    } else {
        n_stack + 1 + opts.air_above_cells
    };

    // --- X-Y extent: identical to the microstrip voxelizer. ---
    let margin = opts.xy_margin_cells as f64 * dx;
    let x0 = layout.bbox.min.x - margin;
    let x1 = layout.bbox.max.x + margin;
    let y0 = layout.bbox.min.y - margin;
    let y1 = layout.bbox.max.y + margin;
    let nx = ((x1 - x0) / dx).ceil() as usize;
    let ny = ((y1 - y0) / dx).ceil() as usize;
    assert!(
        nx > 0 && ny > 0,
        "voxelize_stackup: degenerate x-y extent (nx={nx}, ny={ny})"
    );

    let mut eps = Array3::<f64>::from_elem((nx + 1, ny + 1, nz + 1), 1.0);
    let mut pec_ex = Array3::<bool>::from_elem((nx, ny + 1, nz + 1), false);
    let mut pec_ey = Array3::<bool>::from_elem((nx + 1, ny, nz + 1), false);

    let in_trace = |x: f64, y: f64| {
        layout
            .traces
            .iter()
            .any(|p| point_in_polygon(Point2 { x, y }, p))
    };

    // Dielectric: each layer's E_z-edge band, contiguous from k = 0.
    for i in 0..nx {
        for j in 0..ny {
            for (layer, &ks) in stackup.layers.iter().zip(&k_starts) {
                let ke = ks + ((layer.height_m / dx).round() as usize).max(1);
                for kk in ks..ke {
                    eps[(i, j, kk)] = layer.eps_r;
                }
            }
        }
    }

    // Tangential PEC masks: ground plane (k = 0), trace plane (k_trace,
    // under the traces), lid plane (k = n_stack, whole plane) when lidded.
    for i in 0..nx {
        for j in 0..=ny {
            pec_ex[(i, j, 0)] = true;
            let (x, y) = (x0 + (i as f64 + 0.5) * dx, y0 + j as f64 * dx);
            if in_trace(x, y) {
                pec_ex[(i, j, k_trace)] = true;
            }
            if stackup.lid {
                pec_ex[(i, j, n_stack)] = true;
            }
        }
    }
    for i in 0..=nx {
        for j in 0..ny {
            pec_ey[(i, j, 0)] = true;
            let (x, y) = (x0 + i as f64 * dx, y0 + (j as f64 + 0.5) * dx);
            if in_trace(x, y) {
                pec_ey[(i, j, k_trace)] = true;
            }
            if stackup.lid {
                pec_ey[(i, j, n_stack)] = true;
            }
        }
    }

    let port_cells = layout
        .ports
        .iter()
        .map(|port| {
            let i = (((port.at.x - x0) / dx).floor() as isize).clamp(0, nx as isize - 1) as usize;
            let j = (((port.at.y - y0) / dx).floor() as isize).clamp(0, ny as isize - 1) as usize;
            (i, j, k_trace)
        })
        .collect();

    let grid = YeeGrid::vacuum(nx, ny, nz, dx)
        .with_eps_r_cells(eps)
        .with_pec_mask_ex(pec_ex)
        .with_pec_mask_ey(pec_ey);

    MicrostripModel {
        k_gnd: 0,
        grid,
        dims: (nx, ny, nz),
        dx_m: dx,
        port_cells,
    }
}

#[cfg(test)]
mod rf_tool_tests {
    use super::*;
    use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, StackupLayer, Substrate};

    fn tiny_layout() -> Layout {
        let traces = vec![Polygon::rect(0.0, 0.0, 6.0e-3, 3.0e-3)];
        let bbox = BBox::from_polygons(&traces);
        Layout {
            substrate: Substrate {
                eps_r: 4.4,
                height_m: 1.6e-3,
                loss_tangent: 0.0,
                metal_thickness_m: 35e-6,
            },
            traces,
            ports: vec![PortRef {
                at: Point2::new(0.5e-3, 1.5e-3),
                width_m: 3.0e-3,
                ref_impedance_ohm: 50.0,
            }],
            bbox,
        }
    }

    fn tiny_opts() -> VoxelOptions {
        VoxelOptions {
            dx_m: 0.4e-3,
            xy_margin_cells: 4,
            air_above_cells: 6,
        }
    }

    fn tiny_model() -> MicrostripModel {
        voxelize_microstrip(&tiny_layout(), &tiny_opts())
    }

    #[test]
    fn sigma_cells_map_substrate_only_at_the_pozar_value() {
        let model = tiny_model();
        let tan_d = 0.02;
        let f = 5.0e9;
        let sigma = substrate_sigma_cells(&model, tan_d, f);
        let eps = model.grid.eps_r_cells.as_ref().unwrap();
        let eps_flat = eps.as_slice().unwrap();
        assert_eq!(sigma.len(), eps_flat.len());
        // σ = 2π f ε₀ ε_r tanδ = 2π·5e9·8.854e-12·4.4·0.02 = 0.02448 S/m.
        let expect = std::f64::consts::TAU * f * 8.854_187_817e-12 * 4.4 * tan_d;
        let mut in_sub = 0usize;
        for (s, e) in sigma.iter().zip(eps_flat) {
            if *e > 1.0 {
                assert!((s - expect).abs() < 1e-12, "σ = {s}, want {expect}");
                in_sub += 1;
            } else {
                assert_eq!(*s, 0.0, "air cell must be lossless");
            }
        }
        assert!(in_sub > 0, "no substrate cells found");
    }

    /// Two-layer stack: 2 cells of ε_r=2.2/tanδ=0.02 under 3 cells of
    /// ε_r=4.4/tanδ=0.005, open top (FS.4.2b fixture).
    fn two_layer_stack() -> Stackup {
        Stackup {
            layers: vec![
                StackupLayer {
                    eps_r: 2.2,
                    height_m: 0.8e-3, // 2 cells @ dx = 0.4e-3
                    loss_tangent: 0.02,
                },
                StackupLayer {
                    eps_r: 4.4,
                    height_m: 1.2e-3, // 3 cells
                    loss_tangent: 0.005,
                },
            ],
            lid: false,
        }
    }

    #[test]
    fn stackup_sigma_cells_matches_each_layer_band_exactly() {
        let stack = two_layer_stack();
        let model = voxelize_stackup(&tiny_layout(), &stack, 0, &tiny_opts());
        let f = 5.0e9;
        let sigma = stackup_sigma_cells(&model, &stack, f);
        let eps = model.grid.eps_r_cells.as_ref().unwrap();
        assert_eq!(sigma.len(), eps.len());

        const EPS0: f64 = 8.854_187_817e-12;
        let omega = std::f64::consts::TAU * f;
        let sigma_l0 = omega * EPS0 * 2.2 * 0.02;
        let sigma_l1 = omega * EPS0 * 4.4 * 0.005;
        assert!(sigma_l0 > sigma_l1, "fixture drift: bands must differ");

        let (nx, ny, nz) = model.dims;
        assert_eq!(nz, 2 + 3 + 1 + 6, "fixture drift: open-top stack height");
        let sigma3 = Array3::from_shape_vec((nx + 1, ny + 1, nz + 1), sigma).unwrap();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..5 {
                    let want = if k < 2 { sigma_l0 } else { sigma_l1 };
                    assert!(
                        (sigma3[(i, j, k)] - want).abs() < 1e-12,
                        "band mismatch at k={k}: got {}, want {want}",
                        sigma3[(i, j, k)]
                    );
                }
                // Boundary k=1 (last of layer 0) / k=2 (first of layer 1)
                // are exact, distinct band values — not blended.
                assert!((sigma3[(i, j, 1)] - sigma_l0).abs() < 1e-12);
                assert!((sigma3[(i, j, 2)] - sigma_l1).abs() < 1e-12);
                // Above the stack (air + lid-less top): lossless.
                for k in 5..=nz {
                    assert_eq!(sigma3[(i, j, k)], 0.0, "air cell at k={k} must be lossless");
                }
            }
        }
    }

    #[test]
    fn all_zero_loss_tangent_is_a_provable_no_op() {
        let stack = Stackup {
            layers: vec![
                StackupLayer {
                    eps_r: 2.2,
                    height_m: 0.8e-3,
                    loss_tangent: 0.0,
                },
                StackupLayer {
                    eps_r: 4.4,
                    height_m: 1.2e-3,
                    loss_tangent: 0.0,
                },
            ],
            lid: false,
        };
        let model = voxelize_stackup(&tiny_layout(), &stack, 0, &tiny_opts());
        let sigma = stackup_sigma_cells(&model, &stack, 5.0e9);
        assert!(!sigma.is_empty());
        assert!(sigma.iter().all(|&s| s == 0.0), "loss-off must be all-zero");
    }

    #[test]
    fn single_layer_stackup_sigma_matches_substrate_sigma_cells() {
        let tan_d = 0.02;
        let f = 5.0e9;
        let micro = tiny_model();
        let micro_sigma = substrate_sigma_cells(&micro, tan_d, f);

        let stack = Stackup {
            layers: vec![StackupLayer {
                eps_r: 4.4,
                height_m: 1.6e-3,
                loss_tangent: tan_d,
            }],
            lid: false,
        };
        let multi = voxelize_stackup(&tiny_layout(), &stack, 0, &tiny_opts());
        let multi_sigma = stackup_sigma_cells(&multi, &stack, f);

        // Same k-band formula, same input eps map (voxel-stackup-001 pins
        // the single-layer stackup bit-identical to voxelize_microstrip):
        // the two code paths must agree bit-for-bit, not just numerically.
        assert_eq!(
            micro_sigma, multi_sigma,
            "single-layer stackup_sigma_cells must match substrate_sigma_cells exactly"
        );
    }

    #[test]
    fn via_sets_a_full_ez_column_and_attaches_the_mask() {
        let mut model = tiny_model();
        assert!(model.grid.pec_mask_ez.is_none());
        let k_top = model.port_cells[0].2;
        with_via_at_cell(&mut model, 7, 5, k_top);
        let mask = model.grid.pec_mask_ez.as_ref().unwrap();
        for k in 0..k_top {
            assert!(mask[(7, 5, k)], "via edge k = {k} not set");
        }
        assert!(!mask[(7, 5, k_top.min(mask.dim().2 - 1))] || k_top == mask.dim().2);
        assert!(!mask[(6, 5, 0)], "neighbour column must stay free");
    }
}
