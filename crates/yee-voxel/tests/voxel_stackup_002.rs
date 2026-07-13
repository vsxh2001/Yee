//! Gate `voxel-stackup-002` (FS.4.1, ADR-0221): vias through multilayer
//! stackups —
//!
//! 1. a **through-via** on the FS.4.0 3-layer lidded stack masks exactly
//!    the full `E_z` column `k = 0..nz` at its grid column (ground →
//!    trace plane → lid);
//! 2. a **blind via** (`with_via_between`) masks exactly its
//!    node-plane-to-node-plane band and nothing above or below;
//! 3. the **whole-mask set-cell count equals the sum of the two
//!    columns** — nothing else anywhere was touched (stronger than
//!    spot-checking neighbours, which we also do);
//! 4. `with_via_at_cell` (R.1) is **bit-identical** to
//!    `with_via_between(…, 0, k_top)` — the back-compat delegation.
//!
//! Instant, non-ignored.

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, StackupLayer, Substrate};
use yee_voxel::{
    MicrostripModel, VoxelOptions, voxelize_microstrip, voxelize_stackup, with_through_via_at_cell,
    with_via_at_cell, with_via_between,
};

const DX: f64 = 0.4e-3;

fn line_layout(h_m: f64) -> Layout {
    let traces = vec![Polygon::rect(0.0, 0.0, 20.0e-3, 3.0e-3)];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: Substrate {
            eps_r: 4.4,
            height_m: h_m,
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

fn opts() -> VoxelOptions {
    VoxelOptions {
        dx_m: DX,
        xy_margin_cells: 6,
        air_above_cells: 8,
    }
}

/// The FS.4.0 3-layer lidded stack (voxel-stackup-001's fixture):
/// 2 + 3 + 2 cells, trace at the top of layer 1 (k = 5), lid at k = 7.
fn three_layer_model() -> MicrostripModel {
    let layout = line_layout(1.6e-3); // substrate field unused by the stackup path
    let stack = Stackup {
        layers: vec![
            StackupLayer {
                eps_r: 2.2,
                height_m: 0.8e-3, // 2 cells
                loss_tangent: 0.0,
            },
            StackupLayer {
                eps_r: 4.4,
                height_m: 1.2e-3, // 3 cells
                loss_tangent: 0.0,
            },
            StackupLayer {
                eps_r: 3.0,
                height_m: 0.8e-3, // 2 cells
                loss_tangent: 0.0,
            },
        ],
        lid: true,
    };
    voxelize_stackup(&layout, &stack, 1, &opts())
}

#[test]
fn through_and_blind_vias_mask_exactly_their_columns() {
    let mut model = three_layer_model();
    let (nx, ny, nz) = model.dims;
    assert_eq!(nz, 7, "fixture drift: the FS.4.0 stack is 7 cells");
    assert_eq!(model.port_cells[0].2, 5, "fixture drift: trace at k = 5");
    assert!(model.grid.pec_mask_ez.is_none());

    let (i_t, j_t) = (nx / 2, ny / 2);
    // Through-via: ground (plane 0) -> lid (plane nz = 7), all 7 edges.
    with_through_via_at_cell(&mut model, i_t, j_t);
    // Blind via 3 cells away in x: layer-0 top (plane 2) -> trace plane
    // (plane 5), edges k = 2, 3, 4 only.
    let (i_b, j_b) = (i_t + 3, j_t);
    with_via_between(&mut model, i_b, j_b, 2, 5);

    let mask = model.grid.pec_mask_ez.as_ref().expect("mask attached");
    assert_eq!(mask.dim(), (nx + 1, ny + 1, nz));

    // 1. Through column: every edge 0..7.
    for k in 0..nz {
        assert!(mask[(i_t, j_t, k)], "through-via edge k = {k} not set");
    }
    // 2. Blind column: exactly 2..5.
    for k in 0..nz {
        assert_eq!(
            mask[(i_b, j_b, k)],
            (2..5).contains(&k),
            "blind-via edge k = {k} wrong"
        );
    }
    // 3. Nothing else anywhere: exact whole-mask set count.
    let set = mask.iter().filter(|&&b| b).count();
    assert_eq!(
        set,
        nz + 3,
        "stray E_z PEC cells: expected {} (through) + 3 (blind), found {set}",
        nz
    );
    // Neighbour columns spot-checked too (the voxel-001 idiom).
    for (i, j) in [
        (i_t - 1, j_t),
        (i_t + 1, j_t),
        (i_t, j_t - 1),
        (i_t, j_t + 1),
        (i_b + 1, j_b),
    ] {
        for k in 0..nz {
            assert!(!mask[(i, j, k)], "neighbour ({i},{j}) touched at k = {k}");
        }
    }
}

#[test]
fn via_at_cell_is_bit_identical_to_via_between_from_ground() {
    let layout = line_layout(1.6e-3);
    let mut a = voxelize_microstrip(&layout, &opts());
    let mut b = voxelize_microstrip(&layout, &opts());
    let (i, j, k_top) = a.port_cells[0];
    with_via_at_cell(&mut a, i, j, k_top);
    with_via_between(&mut b, i, j, 0, k_top);
    assert_eq!(
        a.grid.pec_mask_ez.as_ref().unwrap(),
        b.grid.pec_mask_ez.as_ref().unwrap(),
        "R.1 back-compat delegation not bit-identical"
    );
}
