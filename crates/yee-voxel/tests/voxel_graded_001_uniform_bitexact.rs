//! Gate `voxel-graded-001` (FS.0b.1, ADR-0210): the graded voxelizer with
//! **constant** spacing arrays equal to `dx` must produce eps / PEC masks,
//! dims, and port cells **bit-identical** to [`voxelize_microstrip`]'s —
//! the FS.0b.0 bit-exact-on-uniform discipline (compute-018) applied to
//! rasterization. Plus graded-specific sanity: a refined substrate z-stack
//! fills the right cells and coordinate lookups handle nonuniform axes.

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
use yee_voxel::{GradedVoxelGrid, VoxelOptions, voxelize_microstrip, voxelize_microstrip_graded};

const DX: f64 = 0.4e-3;
const MARGIN_CELLS: usize = 4;
const AIR_ABOVE: usize = 6;

/// A small two-trace layout: a line plus an offset stub, so the trace mask
/// has interesting x AND y structure (edges off the domain boundary).
fn layout() -> Layout {
    let traces = vec![
        Polygon::rect(0.0, 0.0, 8.0e-3, 3.0e-3),
        Polygon::rect(3.2e-3, 3.0e-3, 1.6e-3, 2.4e-3),
    ];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, 1.5e-3),
                width_m: 3.0e-3,
                ref_impedance_ohm: 50.0,
            },
            PortRef {
                at: Point2::new(7.5e-3, 1.5e-3),
                width_m: 3.0e-3,
                ref_impedance_ohm: 50.0,
            },
        ],
        bbox,
    }
}

#[test]
fn constant_spacings_match_uniform_voxelizer_bit_exactly() {
    let layout = layout();
    let uniform = voxelize_microstrip(
        &layout,
        &VoxelOptions {
            dx_m: DX,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE,
        },
    );
    let (nx, ny, nz) = uniform.dims;
    let k_top = uniform.port_cells[0].2;

    let graded = voxelize_microstrip_graded(
        &layout,
        &GradedVoxelGrid {
            dx_m: vec![DX; nx],
            dy_m: vec![DX; ny],
            dz_m: vec![DX; nz],
            x0_m: layout.bbox.min.x - MARGIN_CELLS as f64 * DX,
            y0_m: layout.bbox.min.y - MARGIN_CELLS as f64 * DX,
            k_gnd: 0,
            k_top,
        },
    );

    assert_eq!(graded.dims, uniform.dims);
    assert_eq!(graded.port_cells, uniform.port_cells);
    assert_eq!(graded.k_gnd, 0);
    assert_eq!(graded.k_top, k_top);

    // Bit-identical material + PEC arrays (exact equality, no tolerance).
    let u_eps = uniform.grid.eps_r_cells.as_ref().expect("uniform eps map");
    assert_eq!(u_eps.dim(), graded.eps_r_cells.dim());
    assert!(
        u_eps
            .as_slice()
            .unwrap()
            .iter()
            .zip(graded.eps_r_cells.as_slice().unwrap())
            .all(|(a, b)| a == b),
        "eps arrays diverged"
    );
    let u_ex = uniform.grid.pec_mask_ex.as_ref().expect("uniform ex mask");
    assert_eq!(u_ex, &graded.pec_mask_ex, "Ex PEC masks diverged");
    let u_ey = uniform.grid.pec_mask_ey.as_ref().expect("uniform ey mask");
    assert_eq!(u_ey, &graded.pec_mask_ey, "Ey PEC masks diverged");

    // The trace mask is non-trivial (some PEC beyond the ground layer).
    let trace_ex = graded
        .pec_mask_ex
        .indexed_iter()
        .filter(|((_, _, k), v)| *k == k_top && **v)
        .count();
    assert!(trace_ex > 0, "no trace PEC rasterized");
}

#[test]
fn graded_z_stack_fills_the_refined_substrate() {
    let layout = layout();
    // Substrate refined to 4 × 0.4 mm; air grows 0.4 → 0.52 → 0.6760 →
    // then three coarser layers. k_top = 4, nz = 8.
    let dz = vec![
        0.4e-3, 0.4e-3, 0.4e-3, 0.4e-3, 0.52e-3, 0.676e-3, 0.8e-3, 0.8e-3,
    ];
    let graded = voxelize_microstrip_graded(
        &layout,
        &GradedVoxelGrid {
            dx_m: vec![DX; 30],
            dy_m: vec![DX; 22],
            dz_m: dz.clone(),
            x0_m: layout.bbox.min.x - MARGIN_CELLS as f64 * DX,
            y0_m: layout.bbox.min.y - MARGIN_CELLS as f64 * DX,
            k_gnd: 0,
            k_top: 4,
        },
    );
    // Dielectric fills k = 0..4 (no air gap at the ground, ADR-0108), air above.
    for k in 0..8 {
        let want = if k < 4 { 4.4 } else { 1.0 };
        assert_eq!(graded.eps_r_cells[(10, 8, k)], want, "eps at k = {k}");
    }
    // Ground sheet PEC over the whole k = 0 layer; trace PEC only at k = 4.
    assert!(graded.pec_mask_ex[(10, 8, 0)]);
    assert!(graded.pec_mask_ex[(10, 8, 4)]); // (10, 8) sits under the line
    assert!(!graded.pec_mask_ex[(10, 8, 2)]);
    // Node coordinates are the cumulative sums.
    assert_eq!(graded.z_nodes_m.len(), 9);
    let z_top: f64 = dz[..4].iter().sum();
    assert!((graded.z_nodes_m[4] - z_top).abs() < 1e-15);
}

#[test]
fn coordinate_lookup_handles_nonuniform_axes() {
    let layout = layout();
    // x axis: 10 fine cells (0.2 mm) then 20 coarse (0.5 mm), origin −1.6 mm.
    let mut dx = vec![0.2e-3; 10];
    dx.extend(vec![0.5e-3; 20]);
    let graded = voxelize_microstrip_graded(
        &layout,
        &GradedVoxelGrid {
            dx_m: dx,
            dy_m: vec![DX; 22],
            dz_m: vec![DX; 10],
            x0_m: -1.6e-3,
            y0_m: layout.bbox.min.y - MARGIN_CELLS as f64 * DX,
            k_gnd: 0,
            k_top: 4,
        },
    );
    // Fine region: x = −1.6 + i·0.2 mm; x = −0.5 mm is inside cell 5
    // ([−0.6, −0.4) mm). Coarse region starts at node 10 (x = 0.4 mm):
    // x = 1.0 mm is inside cell 11 ([0.9, 1.4) mm).
    assert_eq!(graded.cell_at_x(-0.5e-3), 5);
    assert_eq!(graded.cell_at_x(1.0e-3), 11);
    // Clamping at both ends.
    assert_eq!(graded.cell_at_x(-99.0), 0);
    assert_eq!(graded.cell_at_x(99.0), 29);
    // Port 0 at x = 0.5 mm lands in coarse cell 10 ([0.4, 0.9) mm).
    assert_eq!(graded.port_cells[0].0, 10);
}
