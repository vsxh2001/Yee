//! `voxel_003` — finite-board gate (FS.2b.1, ADR-0207): the dielectric
//! slab and the ground sheet end `board_margin_m` beyond the layout bbox
//! (a real PCB in air) instead of filling the domain — the fixture that
//! lets an NTFF equivalence box pass through homogeneous air. Bounds
//! pinned cell-exactly; the infinite path stays bit-identical.

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
use yee_voxel::{VoxelOptions, voxelize_finite_board, voxelize_microstrip_open};

fn line_layout() -> Layout {
    let traces = vec![Polygon::rect(0.0, 0.0, 20.0e-3, 3.0e-3)];
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
            at: Point2::new(0.0, 1.5e-3),
            width_m: 3.0e-3,
            ref_impedance_ohm: 50.0,
        }],
        bbox,
    }
}

#[test]
fn voxel_003_board_bounds_are_exact_and_outside_is_air() {
    let layout = line_layout();
    let opts = VoxelOptions {
        dx_m: 0.5e-3,
        xy_margin_cells: 20,
        air_above_cells: 8,
    };
    let air_below = 6;
    let board_margin = 3.0e-3; // 6 cells beyond the bbox
    let model = voxelize_finite_board(&layout, &opts, board_margin, air_below);
    let (nx, ny, _) = model.dims;
    let kg = model.k_gnd;
    assert_eq!(kg, air_below);

    let eps = model.grid.eps_r_cells.as_ref().unwrap();
    let ex = model.grid.pec_mask_ex.as_ref().unwrap();
    let dx = model.dx_m;
    let x0 = layout.bbox.min.x - 20.0 * dx;
    let y0 = layout.bbox.min.y - 20.0 * dx;
    // Cell centre inside the board rect (bbox + 3 mm) → dielectric + PEC;
    // outside → air, no ground.
    let inside = |x: f64, y: f64| {
        x >= -board_margin
            && x <= 20.0e-3 + board_margin
            && y >= -board_margin
            && y <= 3.0e-3 + board_margin
    };
    for &(i, j) in &[
        (5usize, 5usize), // deep in the air margin
        (20, 20),         // on the board (bbox origin area)
        (nx - 5, ny - 5), // far corner, air
        (nx / 2, ny / 2), // board centre
    ] {
        let (xc, yc) = (x0 + (i as f64 + 0.5) * dx, y0 + (j as f64 + 0.5) * dx);
        if inside(xc, yc) {
            assert_eq!(eps[(i, j, kg)], 4.4, "dielectric at ({i}, {j})");
            assert!(ex[(i, j, kg)], "ground PEC at ({i}, {j})");
        } else {
            assert_eq!(eps[(i, j, kg)], 1.0, "air at ({i}, {j})");
            assert!(!ex[(i, j, kg)], "no ground at ({i}, {j})");
        }
    }
    // Domain floor open (lifted stack).
    assert!(!ex[(nx / 2, ny / 2, 0)], "floor must be open");

    // Infinite path unchanged: finite with a huge margin == open stack.
    let a = voxelize_microstrip_open(&layout, &opts, air_below);
    let b = voxelize_finite_board(&layout, &opts, 10.0, air_below);
    assert_eq!(
        a.grid.eps_r_cells.as_ref().unwrap(),
        b.grid.eps_r_cells.as_ref().unwrap(),
        "huge board margin must reproduce the infinite slab"
    );
    assert_eq!(
        a.grid.pec_mask_ex.as_ref().unwrap(),
        b.grid.pec_mask_ex.as_ref().unwrap()
    );
}
