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
//! per-cell `ε_r` and a per-component `E_z` PEC mask already attached, ready
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
//! - `k = 0`: ground plane — PEC over the whole x-y extent.
//! - `k = 1 ..= n_sub`: substrate, `n_sub = round(height_m / dx)` (≥ 1),
//!   `ε_r = layout.substrate.eps_r`.
//! - `k = k_top = 1 + n_sub`: top-metal layer — PEC where a trace polygon
//!   covers the cell centre; `ε_r = 1` elsewhere.
//! - `k = k_top + 1 ..= k_top + air_above_cells`: air (`ε_r = 1`).
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
use yee_fdtd::YeeGrid;
use yee_layout::{Layout, Point2, Polygon};

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
    /// The Yee grid with per-cell `ε_r` and an `E_z` PEC mask attached.
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
    let n_sub = ((layout.substrate.height_m / dx).round() as usize).max(1);
    let k_top = 1 + n_sub; // ground (k=0) + n_sub substrate layers -> top-metal layer
    let nz = k_top + 1 + opts.air_above_cells;

    let eps_r_sub = layout.substrate.eps_r;

    // --- Material arrays at YeeGrid's exact required shapes. ---
    // `with_eps_r_cells` requires `(nx+1, ny+1, nz+1)`.
    // `with_pec_mask_ez`  requires `(nx+1, ny+1, nz)`  (the `ez` component shape).
    let mut eps = Array3::<f64>::from_elem((nx + 1, ny + 1, nz + 1), 1.0);
    let mut pec = Array3::<bool>::from_elem((nx + 1, ny + 1, nz), false);

    // Cell (i, j) centre in metres.
    let cell_center = |i: usize, j: usize| Point2 {
        x: x0 + (i as f64 + 0.5) * dx,
        y: y0 + (j as f64 + 0.5) * dx,
    };

    // Logical cell (i, j, k) is represented at node index (i, j, k) of each
    // array. Both arrays are oversized by +1 in x and y relative to the cell
    // count, so writing node (i, j, *) for i in 0..nx, j in 0..ny is always in
    // bounds. `eps` is additionally +1 in z (covers k in 0..=nz); `pec` has
    // exactly `nz` z-samples (covers k in 0..nz), matching the `Ez` staggering.
    for i in 0..nx {
        for j in 0..ny {
            let c = cell_center(i, j);

            // Ground plane: PEC across the whole layer (k = 0).
            pec[(i, j, 0)] = true;

            // Substrate dielectric for k = 1 ..= n_sub.
            for k in 1..=n_sub {
                eps[(i, j, k)] = eps_r_sub;
            }
            // (Air ε_r = 1.0 above k_top is already the array default.)

            // Top-metal layer (k = k_top): PEC where a trace covers the centre.
            if layout.traces.iter().any(|p| point_in_polygon(c, p)) {
                pec[(i, j, k_top)] = true;
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
        .with_pec_mask_ez(pec);

    MicrostripModel {
        grid,
        dims: (nx, ny, nz),
        dx_m: dx,
        port_cells,
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
