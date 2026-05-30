//! `voxel_001` — straight microstrip line voxelization gate (Filter F1.1a,
//! ADR-0091 §Consequences / spec §DoD 4).
//!
//! Builds a single straight microstrip line `Layout` (one rectangle trace,
//! FR-4 `ε_r = 4.4`, `h = 1.6 mm`), voxelizes it, and asserts the grid
//! dimensions, the substrate/air `ε_r`, the all-PEC ground layer, the
//! trace-layer PEC-cell count, on/off-trace cell PEC state, and the port-cell
//! mapping. No FDTD time-stepping — this runs in milliseconds.

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

/// Trace width along `x`, metres (3.0 mm).
const W: f64 = 3.0e-3;
/// Trace length along `y`, metres (20.0 mm).
const L: f64 = 20.0e-3;
/// Substrate height, metres (1.6 mm).
const H: f64 = 1.6e-3;
/// Cell size, metres (0.5 mm).
const DX: f64 = 0.5e-3;
/// Air margin in cells.
const MARGIN: usize = 4;
/// Air layers above the top metal.
const AIR_ABOVE: usize = 8;
/// Substrate relative permittivity (FR-4).
const EPS_R: f64 = 4.4;

fn straight_line_layout() -> Layout {
    let substrate = Substrate {
        eps_r: EPS_R,
        height_m: H,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    // Lower-left at the origin; W along x, L along y.
    let trace = Polygon::rect(0.0, 0.0, W, L);
    let traces = vec![trace];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate,
        traces,
        // One port at the bottom-centre end of the line, one at the top-centre.
        ports: vec![
            PortRef {
                at: Point2::new(W / 2.0, 0.0),
                width_m: W,
                ref_impedance_ohm: 50.0,
            },
            PortRef {
                at: Point2::new(W / 2.0, L),
                width_m: W,
                ref_impedance_ohm: 50.0,
            },
        ],
        bbox,
    }
}

#[test]
fn voxel_001_microstrip_line() {
    let layout = straight_line_layout();
    let opts = VoxelOptions {
        dx_m: DX,
        xy_margin_cells: MARGIN,
        air_above_cells: AIR_ABOVE,
    };
    let model = voxelize_microstrip(&layout, &opts);

    // --- Hand-computed dims. ---
    // x extent: bbox.x ∈ [0, 3mm], padded ±2mm -> [-2, 5]mm = 7mm -> 14 cells.
    // y extent: bbox.y ∈ [0, 20mm], padded ±2mm -> [-2, 22]mm = 24mm -> 48 cells.
    // n_sub = round(1.6/0.5) = 3; trace at k_top = n_sub = 3 (so the ground-to-
    // trace gap is n_sub·dx = h, all dielectric — no air series gap, ADR-0108);
    // nz = k_top + 1 + 8 = 12.
    let exp_nx = 14;
    let exp_ny = 48;
    let n_sub = 3;
    let k_top = n_sub; // 3
    let exp_nz = k_top + 1 + AIR_ABOVE; // 12
    assert_eq!(
        model.dims,
        (exp_nx, exp_ny, exp_nz),
        "grid dims mismatch (computed nx,ny,nz)"
    );
    assert_eq!(model.dx_m, DX);

    let (nx, ny, nz) = model.dims;
    let grid = &model.grid;

    // --- Substrate vs air ε_r. ---
    let eps = grid
        .eps_r_cells
        .as_ref()
        .expect("eps_r_cells must be attached");
    assert_eq!(eps.dim(), (nx + 1, ny + 1, nz + 1));
    // A substrate cell (mid-slab, k=2) carries ε_r ≈ 4.4.
    assert!(
        (eps[(nx / 2, ny / 2, 2)] - EPS_R).abs() < 1e-9,
        "substrate cell ε_r should be {EPS_R}, got {}",
        eps[(nx / 2, ny / 2, 2)]
    );
    // An air cell well above the top metal carries ε_r = 1.0.
    assert!(
        (eps[(nx / 2, ny / 2, nz - 1)] - 1.0).abs() < 1e-9,
        "air cell ε_r should be 1.0, got {}",
        eps[(nx / 2, ny / 2, nz - 1)]
    );

    // --- Ground plane (k=0): tangential Ex AND Ey fully PEC. ---
    // A horizontal PEC plane zeroes the in-plane (tangential) field, so both
    // Ex and Ey are masked at k=0; the normal Ez is NOT (ADR-0091 review fix).
    let pec_ex = grid
        .pec_mask_ex
        .as_ref()
        .expect("pec_mask_ex must be attached");
    let pec_ey = grid
        .pec_mask_ey
        .as_ref()
        .expect("pec_mask_ey must be attached");
    assert_eq!(pec_ex.dim(), (nx, ny + 1, nz + 1));
    assert_eq!(pec_ey.dim(), (nx + 1, ny, nz + 1));
    for i in 0..nx {
        for j in 0..=ny {
            assert!(pec_ex[(i, j, 0)], "ground Ex node ({i},{j},0) must be PEC");
        }
    }
    for i in 0..=nx {
        for j in 0..ny {
            assert!(pec_ey[(i, j, 0)], "ground Ey node ({i},{j},0) must be PEC");
        }
    }

    // --- Trace-layer (k=k_top) tangential-PEC count (Ex nodes) ≈ (w·l)/dx². ---
    let expected_trace_cells = (W * L) / (DX * DX); // 240.0
    let mut trace_count = 0usize;
    for i in 0..nx {
        for j in 0..=ny {
            if pec_ex[(i, j, k_top)] {
                trace_count += 1;
            }
        }
    }
    // Staggered Ex nodes near the trace boundary add a row/col of rounding slack.
    let nrow = (W / DX).round() as usize; // 6
    let ncol = (L / DX).round() as usize; // 40
    let slack = (nrow + ncol + 2) as f64;
    assert!(
        (trace_count as f64 - expected_trace_cells).abs() <= slack,
        "trace Ex PEC-node count {trace_count} should be ≈ {expected_trace_cells} (±{slack})"
    );

    // --- An Ex node under the trace is PEC; one in the margin is not. ---
    let ci = (((W / 2.0 - (layout.bbox.min.x - MARGIN as f64 * DX)) / DX).floor()) as usize;
    let cj = (((L / 2.0 - (layout.bbox.min.y - MARGIN as f64 * DX)) / DX).floor()) as usize;
    assert!(
        pec_ex[(ci, cj, k_top)],
        "Ex node ({ci},{cj},{k_top}) under the trace centre must be PEC"
    );
    // A corner node (i=0, j=0) is in the air margin, well off the trace.
    assert!(
        !pec_ex[(0, 0, k_top)],
        "margin Ex node (0,0,{k_top}) must not be PEC"
    );

    // --- Ports map to in-range cells at the top-metal layer. ---
    assert_eq!(model.port_cells.len(), layout.ports.len());
    for &(i, j, k) in &model.port_cells {
        assert!(i < nx, "port i={i} out of range (nx={nx})");
        assert!(j < ny, "port j={j} out of range (ny={ny})");
        assert!(k < nz, "port k={k} out of range (nz={nz})");
        assert_eq!(k, k_top, "ports should sit at the top-metal layer");
    }
}
