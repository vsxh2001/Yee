//! Gate `voxel-stackup-001` (FS.4.0, ADR-0215): the stackup voxelizer is
//! a strict generalization of the microstrip one —
//!
//! 1. a single-layer open stackup reproduces [`voxelize_microstrip`]
//!    **bit-identically** (ε and both PEC masks — the FS.0b idiom);
//! 2. a 3-layer stack puts each ε_r in exactly its k-band with **no air
//!    gap anywhere** (the ADR-0108 lesson generalized), the trace on the
//!    chosen interface, and the lid masked across the whole plane;
//! 3. the symmetric-stripline convenience yields a mid-stack trace under
//!    a lidded, homogeneous fill.
//!
//! Instant, non-ignored.

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, StackupLayer, Substrate};
use yee_voxel::{VoxelOptions, voxelize_microstrip, voxelize_stackup};

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

#[test]
fn single_layer_stackup_is_bit_identical_to_microstrip() {
    let layout = line_layout(1.6e-3);
    let micro = voxelize_microstrip(&layout, &opts());
    let stack = Stackup {
        layers: vec![StackupLayer {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
        }],
        lid: false,
    };
    let multi = voxelize_stackup(&layout, &stack, 0, &opts());

    assert_eq!(micro.dims, multi.dims);
    assert_eq!(micro.port_cells, multi.port_cells);
    assert_eq!(micro.k_gnd, multi.k_gnd);
    assert_eq!(
        micro.grid.eps_r_cells.as_ref().unwrap(),
        multi.grid.eps_r_cells.as_ref().unwrap(),
        "eps not bit-identical"
    );
    assert_eq!(
        micro.grid.pec_mask_ex.as_ref().unwrap(),
        multi.grid.pec_mask_ex.as_ref().unwrap(),
        "Ex mask not bit-identical"
    );
    assert_eq!(
        micro.grid.pec_mask_ey.as_ref().unwrap(),
        multi.grid.pec_mask_ey.as_ref().unwrap(),
        "Ey mask not bit-identical"
    );
}

#[test]
fn three_layer_bands_are_contiguous_and_correct() {
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
    // Trace buried at the top of layer 1 (k = 2 + 3 = 5).
    let model = voxelize_stackup(&layout, &stack, 1, &opts());
    let (nx, ny, nz) = model.dims;
    // Lidded: the domain ends at the stack top (2 + 3 + 2 = 7 cells).
    assert_eq!(nz, 7);
    assert_eq!(model.port_cells[0].2, 5, "trace plane at the interface");

    let eps = model.grid.eps_r_cells.as_ref().unwrap();
    let (ic, jc) = (nx / 2, ny / 2);
    let expect = [2.2, 2.2, 4.4, 4.4, 4.4, 3.0, 3.0];
    for (k, &e) in expect.iter().enumerate() {
        assert_eq!(
            eps[(ic, jc, k)],
            e,
            "eps band wrong at k={k}: no-gap contiguity broken"
        );
    }

    // Lid: whole-plane tangential PEC at the stack top.
    let pec_ex = model.grid.pec_mask_ex.as_ref().unwrap();
    let pec_ey = model.grid.pec_mask_ey.as_ref().unwrap();
    for i in 0..nx {
        for j in 0..=ny {
            assert!(pec_ex[(i, j, 7)], "lid Ex hole at ({i},{j})");
        }
    }
    for i in 0..=nx {
        for j in 0..ny {
            assert!(pec_ey[(i, j, 7)], "lid Ey hole at ({i},{j})");
        }
    }
    // Trace: masked at k=5 under the strip, NOT the whole plane.
    assert!(pec_ex[(ic, jc, 5)], "trace missing at the interface");
    assert!(!pec_ex[(ic, 0, 5)], "trace mask leaked to the margin");
}

#[test]
fn symmetric_stripline_has_mid_stack_trace_under_a_lid() {
    let layout = line_layout(3.2e-3);
    let stack = Stackup::symmetric_stripline(4.4, 3.2e-3);
    let model = voxelize_stackup(&layout, &stack, 0, &opts());
    let (_, _, nz) = model.dims;
    // b = 3.2 mm at dx = 0.4 mm: 4 cells per half, 8 total, trace at 4.
    assert_eq!(nz, 8);
    assert_eq!(model.port_cells[0].2, 4);
    // Homogeneous fill, ground-to-lid.
    let eps = model.grid.eps_r_cells.as_ref().unwrap();
    for k in 0..8 {
        assert_eq!(eps[(model.dims.0 / 2, model.dims.1 / 2, k)], 4.4);
    }
}
