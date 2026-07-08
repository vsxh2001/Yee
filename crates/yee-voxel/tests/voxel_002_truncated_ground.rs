//! `voxel_002` — truncated-ground gate (FS.1a.0, ADR-0205).
//!
//! The quasi-Yagi's reflector is the ground plane's truncation edge, so the
//! truncation must be *exact*: this gate hand-computes the `k = 0` PEC mask
//! population before and after [`yee_voxel::truncate_ground_at_cell`] and
//! pins the edge rule (`Ex` kept for `i < i_ground_end`, `Ey` for
//! `i ≤ i_ground_end`), the untouched trace layer, and the full-ground
//! default. Milliseconds, no time-stepping.

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
use yee_voxel::{VoxelOptions, truncate_ground_at_cell, voxelize_microstrip};

const W: f64 = 3.0e-3;
const L: f64 = 20.0e-3;
const H: f64 = 1.6e-3;
const DX: f64 = 0.5e-3;

fn line_layout() -> Layout {
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: H,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let trace = Polygon::rect(0.0, 0.0, L, W); // along +x
    let traces = vec![trace];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate,
        traces,
        ports: vec![PortRef {
            at: Point2::new(0.0, W / 2.0),
            width_m: W,
            ref_impedance_ohm: 50.0,
        }],
        bbox,
    }
}

fn opts() -> VoxelOptions {
    VoxelOptions {
        dx_m: DX,
        xy_margin_cells: 4,
        air_above_cells: 8,
    }
}

/// Count `true` cells in the `k = 0` layer of a mask.
fn k0_count(mask: &ndarray::Array3<bool>) -> usize {
    let (d0, d1, _) = mask.dim();
    (0..d0)
        .flat_map(|i| (0..d1).map(move |j| (i, j)))
        .filter(|&(i, j)| mask[(i, j, 0)])
        .count()
}

#[test]
fn voxel_002_truncated_ground_edge_is_exact() {
    let layout = line_layout();
    let mut model = voxelize_microstrip(&layout, &opts());
    let (nx, ny, _) = model.dims;

    // Full ground before truncation: every k = 0 node PEC.
    let ex_full = k0_count(model.grid.pec_mask_ex.as_ref().unwrap());
    let ey_full = k0_count(model.grid.pec_mask_ey.as_ref().unwrap());
    assert_eq!(ex_full, nx * (ny + 1), "full-ground Ex layer");
    assert_eq!(ey_full, (nx + 1) * ny, "full-ground Ey layer");
    let k_top = model.port_cells[0].2;
    let ex_trace_before = {
        let ex = model.grid.pec_mask_ex.as_ref().unwrap();
        let (d0, d1, _) = ex.dim();
        (0..d0)
            .flat_map(|i| (0..d1).map(move |j| (i, j)))
            .filter(|&(i, j)| ex[(i, j, k_top)])
            .count()
    };

    // Truncate mid-grid.
    let g = nx / 2;
    truncate_ground_at_cell(&mut model, g);

    let ex = model.grid.pec_mask_ex.as_ref().unwrap();
    let ey = model.grid.pec_mask_ey.as_ref().unwrap();

    // Population: Ex keeps i < g, Ey keeps i <= g.
    assert_eq!(k0_count(ex), g * (ny + 1), "truncated Ex population");
    assert_eq!(k0_count(ey), (g + 1) * ny, "truncated Ey population");

    // The edge itself, both sides, every j.
    for j in 0..=ny {
        assert!(
            ex[(g - 1, j, 0)],
            "Ex just inside the edge (i = g-1, j = {j})"
        );
        assert!(!ex[(g, j, 0)], "Ex just outside the edge (i = g, j = {j})");
    }
    for j in 0..ny {
        assert!(ey[(g, j, 0)], "Ey on the edge plane (i = g, j = {j})");
        assert!(!ey[(g + 1, j, 0)], "Ey past the edge (i = g+1, j = {j})");
    }

    // Trace layer untouched.
    let ex_trace_after = {
        let (d0, d1, _) = ex.dim();
        (0..d0)
            .flat_map(|i| (0..d1).map(move |j| (i, j)))
            .filter(|&(i, j)| ex[(i, j, k_top)])
            .count()
    };
    assert_eq!(ex_trace_after, ex_trace_before, "trace masks untouched");
}

#[test]
fn voxel_002_full_ground_is_the_untruncated_default() {
    let layout = line_layout();
    let reference = voxelize_microstrip(&layout, &opts());
    let mut noop = voxelize_microstrip(&layout, &opts());
    let (nx, _, _) = noop.dims;
    // i_ground_end >= nx is the documented no-op.
    truncate_ground_at_cell(&mut noop, nx);
    assert_eq!(
        reference.grid.pec_mask_ex.as_ref().unwrap(),
        noop.grid.pec_mask_ex.as_ref().unwrap(),
        "no-op truncation must leave Ex bit-identical"
    );
    assert_eq!(
        reference.grid.pec_mask_ey.as_ref().unwrap(),
        noop.grid.pec_mask_ey.as_ref().unwrap(),
        "no-op truncation must leave Ey bit-identical"
    );
}
