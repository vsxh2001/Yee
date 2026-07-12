//! Gate `voxel-poly-001` (FS.3.2a, ADR-0217): non-axis-aligned polygon
//! edges rasterize correctly through the general even-odd
//! `point_in_polygon` path — pinned on a hand-computable 45°-cut square
//! and on the [`double_jog`] mitered-vs-square pair.

use yee_layout::{BBox, Layout, MiterStyle, Point2, Polygon, PortRef, Substrate, double_jog};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

fn substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

fn mask_count(l: &Layout, dx: f64) -> (usize, yee_voxel::MicrostripModel) {
    let model = voxelize_microstrip(
        l,
        &VoxelOptions {
            dx_m: dx,
            xy_margin_cells: 2,
            air_above_cells: 3,
        },
    );
    let n = model
        .grid
        .pec_mask_ex
        .as_ref()
        .unwrap()
        .iter()
        .filter(|&&m| m)
        .count()
        + model
            .grid
            .pec_mask_ey
            .as_ref()
            .unwrap()
            .iter()
            .filter(|&&m| m)
            .count();
    (n, model)
}

/// A 10×10 mm square with its top-right corner cut at 45° (legs 4.5 mm,
/// chosen so the cut line x+y = 15.5 mm never passes through a cell
/// centre — centre coordinate sums are integers at dx = 1 mm, and a
/// boundary-exact centre is even-odd-implementation-defined), voxelized
/// at dx = 1 mm: centres above the line must be air, the rest metal.
#[test]
fn diagonal_cut_staircase_is_exact() {
    let poly = Polygon {
        verts: vec![
            Point2::new(0.0, 0.0),
            Point2::new(10.0e-3, 0.0),
            Point2::new(10.0e-3, 5.5e-3),
            Point2::new(5.5e-3, 10.0e-3),
            Point2::new(0.0, 10.0e-3),
        ],
    };
    let traces = vec![poly];
    let bbox = BBox::from_polygons(&traces);
    let l = Layout {
        substrate: substrate(),
        traces,
        ports: vec![PortRef {
            at: Point2::new(0.5e-3, 5.0e-3),
            width_m: 3.0e-3,
            ref_impedance_ohm: 50.0,
        }],
        bbox,
    };
    let (_, model) = mask_count(&l, 1.0e-3);
    let k_top = model.port_cells[0].2;
    let mask = model.grid.pec_mask_ex.as_ref().unwrap();
    // Cell (i, j) covers [i·dx, (i+1)·dx] from the bbox min minus margin
    // (2 cells): trace cell (u, v) in layout frame is grid (u+2, v+2).
    for u in 0..10 {
        for v in 0..10 {
            let (cx, cy) = (u as f64 + 0.5, v as f64 + 0.5); // mm
            let expected = cx + cy < 15.5; // under the 45° cut line
            assert_eq!(
                mask[(u + 2, v + 2, k_top)],
                expected,
                "cell centre ({cx}, {cy}) mm: expected metal = {expected}"
            );
        }
    }
}

/// The mitered double-jog masks strictly fewer PEC edges than the square
/// one (four corner triangles removed), on the identical grid.
#[test]
fn mitered_jog_masks_strictly_fewer_cells_than_square() {
    let s = substrate();
    let sq = double_jog(&s, 3.0e-3, 24.0e-3, 9.0e-3, 9.0e-3, MiterStyle::Square);
    let mi = double_jog(
        &s,
        3.0e-3,
        24.0e-3,
        9.0e-3,
        9.0e-3,
        MiterStyle::Mitered { f: 0.7 },
    );
    let (n_sq, m_sq) = mask_count(&sq, 0.3e-3);
    let (n_mi, m_mi) = mask_count(&mi, 0.3e-3);
    assert_eq!(m_sq.dims, m_mi.dims, "variants must share the grid");
    assert!(
        n_mi < n_sq,
        "mitered ({n_mi}) must mask fewer PEC edges than square ({n_sq})"
    );
    // Each cut removes a triangle of area (f·w)²/2 = 2.205 mm²; four cuts
    // ≈ 8.82 mm² ≈ 98 cells at dx = 0.3 mm. Allow generous staircase
    // slack, but the removal must be the right order (per mask array).
    let removed = n_sq - n_mi;
    assert!(
        (100..400).contains(&removed),
        "removed {removed} mask edges; expected ~196 (2 arrays × ~98 cells)"
    );
}
