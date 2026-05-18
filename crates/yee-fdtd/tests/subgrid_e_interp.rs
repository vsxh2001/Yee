//! Phase 2.fdtd.7.0 Q3 — coarse → fine `E_t` spatial + temporal interpolation.
//!
//! Covers the three behavioural invariants of [`SubgridRegion`]'s Q3 surface:
//!
//! 1. A constant parent `E_t` maps to a constant fine boundary `E_t` (within
//!    `f64::EPSILON · 10`) on **all six** interface faces.
//! 2. A linear-in-x parent gradient reproduces exactly (within
//!    `f64::EPSILON · 100`) at fine-edge midpoints, including the
//!    extrapolation cells at the very edge of the subgrid.
//! 3. The temporal blend is linear in `frac` so a (0, 1) snapshot pair
//!    sampled at `frac = 0.25` and `frac = 0.75` reads back as `0.25` and
//!    `0.75` respectively.
//!
//! Plus a dispatch test that proves no face is silently skipped: build a
//! parent grid with a non-symmetric `E_t` pattern, run the interp, and
//! confirm every face's tangential `E_t` carries a non-trivial value.

use yee_fdtd::{SubgridRegion, YeeGrid};

const N: usize = 16;
const DX: f64 = 1.0e-3;
const TOL_UNIFORM: f64 = f64::EPSILON * 10.0;
const TOL_LINEAR: f64 = f64::EPSILON * 100.0;

fn parent() -> YeeGrid {
    YeeGrid::vacuum(N, N, N, DX)
}

/// Drive every coarse `E` array to a uniform constant.
fn set_parent_uniform(grid: &mut YeeGrid, val: f64) {
    grid.ex.fill(val);
    grid.ey.fill(val);
    grid.ez.fill(val);
}

// -----------------------------------------------------------------------
// Test 1 — uniform parent field → uniform fine boundary, all six faces.
// -----------------------------------------------------------------------

#[test]
fn uniform_field_maps_to_uniform_fine_e_t() {
    let mut p = parent();
    let val = 0.42;
    set_parent_uniform(&mut p, val);

    let lo = (4, 4, 4);
    let hi = (10, 10, 10);
    let mut region = SubgridRegion::new(&p, lo, hi).expect("valid subgrid");

    // Start == end snapshot: temporal blend collapses out, leaving the
    // pure spatial interp. A constant field interpolates to a constant.
    region.snapshot_coarse_e_t(&p);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.5);

    let fine = region.fine_grid();
    let fine_nx = fine.nx;
    let fine_ny = fine.ny;
    let fine_nz = fine.nz;

    // ±x faces — check E_y and E_z.
    for face_i in [0usize, fine_nx] {
        for j in 0..fine_ny {
            for k in 0..=fine_nz {
                let got = fine.ey[(face_i, j, k)];
                assert!(
                    (got - val).abs() <= TOL_UNIFORM,
                    "+/-x face E_y at ({face_i}, {j}, {k}) = {got}, expected {val}"
                );
            }
        }
        for j in 0..=fine_ny {
            for k in 0..fine_nz {
                let got = fine.ez[(face_i, j, k)];
                assert!(
                    (got - val).abs() <= TOL_UNIFORM,
                    "+/-x face E_z at ({face_i}, {j}, {k}) = {got}, expected {val}"
                );
            }
        }
    }

    // ±y faces — check E_x and E_z.
    for face_j in [0usize, fine_ny] {
        for i in 0..fine_nx {
            for k in 0..=fine_nz {
                let got = fine.ex[(i, face_j, k)];
                assert!(
                    (got - val).abs() <= TOL_UNIFORM,
                    "+/-y face E_x at ({i}, {face_j}, {k}) = {got}, expected {val}"
                );
            }
        }
        for i in 0..=fine_nx {
            for k in 0..fine_nz {
                let got = fine.ez[(i, face_j, k)];
                assert!(
                    (got - val).abs() <= TOL_UNIFORM,
                    "+/-y face E_z at ({i}, {face_j}, {k}) = {got}, expected {val}"
                );
            }
        }
    }

    // ±z faces — check E_x and E_y.
    for face_k in [0usize, fine_nz] {
        for i in 0..fine_nx {
            for j in 0..=fine_ny {
                let got = fine.ex[(i, j, face_k)];
                assert!(
                    (got - val).abs() <= TOL_UNIFORM,
                    "+/-z face E_x at ({i}, {j}, {face_k}) = {got}, expected {val}"
                );
            }
        }
        for i in 0..=fine_nx {
            for j in 0..fine_ny {
                let got = fine.ey[(i, j, face_k)];
                assert!(
                    (got - val).abs() <= TOL_UNIFORM,
                    "+/-z face E_y at ({i}, {j}, {face_k}) = {got}, expected {val}"
                );
            }
        }
    }
}

// -----------------------------------------------------------------------
// Test 2 — linear gradient in x preserved on the +z face's E_x.
// -----------------------------------------------------------------------

#[test]
fn linear_gradient_preserved() {
    let mut p = parent();
    let a = 1.0_f64;
    let b = 3.0_f64;

    // Drive parent.ex with E_x(i, j, k) = a + b · x_coarse where the
    // E_x edge at (i, j, k) sits at coarse-x = (i + 0.5)·dx (Taflove
    // staggering). Other E components do not matter for this face.
    for ((i, j, k), v) in p.ex.indexed_iter_mut() {
        let x = (i as f64 + 0.5) * DX;
        let _ = (j, k);
        *v = a + b * x;
    }

    let lo = (4, 4, 4);
    let hi = (10, 10, 10);
    let mut region = SubgridRegion::new(&p, lo, hi).expect("valid subgrid");

    region.snapshot_coarse_e_t(&p);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.5);

    let fine = region.fine_grid();
    let dx_fine = DX / 2.0;

    // +z face: fine E_x at (i_f, j_f, fine.nz) sits at fine coordinates
    // x_fine = (i_f + 0.5) · dx_fine + lo.0 · dx_coarse.
    let face_k = fine.nz;
    for j_f in 0..=fine.ny {
        for i_f in 0..fine.nx {
            let x_fine = (i_f as f64 + 0.5) * dx_fine + (lo.0 as f64) * DX;
            let expected = a + b * x_fine;
            let got = fine.ex[(i_f, j_f, face_k)];
            assert!(
                (got - expected).abs() <= TOL_LINEAR * expected.abs().max(1.0),
                "+z face E_x at ({i_f}, {j_f}, {face_k}) = {got}, expected {expected}"
            );
        }
    }
}

// -----------------------------------------------------------------------
// Test 3 — temporal blend at frac = 0.25.
// -----------------------------------------------------------------------

#[test]
fn temporal_blend_at_frac_one_quarter() {
    let mut p = parent();
    let lo = (4, 4, 4);
    let hi = (10, 10, 10);
    let mut region = SubgridRegion::new(&p, lo, hi).expect("valid subgrid");

    // Snapshot start with all-zero parent.
    set_parent_uniform(&mut p, 0.0);
    region.snapshot_coarse_e_t(&p);

    // Snapshot end with all-one parent.
    set_parent_uniform(&mut p, 1.0);
    region.snapshot_coarse_e_t_end(&p);

    region.interpolate_coarse_e_to_fine(0.25);

    let fine = region.fine_grid();
    // Inspect every tangential component on every face — they all must
    // read back as the linear blend 0.25·(1 - 0) + 0 = 0.25.
    let expected = 0.25;
    let probes = sample_all_face_e_t(fine);
    for (label, got) in probes {
        assert!(
            (got - expected).abs() <= TOL_UNIFORM,
            "temporal blend at frac=0.25: {label} = {got}, expected {expected}"
        );
    }
}

// -----------------------------------------------------------------------
// Test 4 — temporal blend at frac = 0.75.
// -----------------------------------------------------------------------

#[test]
fn temporal_blend_at_frac_three_quarter() {
    let mut p = parent();
    let lo = (4, 4, 4);
    let hi = (10, 10, 10);
    let mut region = SubgridRegion::new(&p, lo, hi).expect("valid subgrid");

    set_parent_uniform(&mut p, 0.0);
    region.snapshot_coarse_e_t(&p);
    set_parent_uniform(&mut p, 1.0);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.75);

    let fine = region.fine_grid();
    let expected = 0.75;
    let probes = sample_all_face_e_t(fine);
    for (label, got) in probes {
        assert!(
            (got - expected).abs() <= TOL_UNIFORM,
            "temporal blend at frac=0.75: {label} = {got}, expected {expected}"
        );
    }
}

// -----------------------------------------------------------------------
// Test 5 — six-face dispatch.
// -----------------------------------------------------------------------

#[test]
fn interpolation_face_dispatch_covers_all_six_faces() {
    let mut p = parent();

    // Asymmetric parent: every coarse E component is non-zero in a way
    // that depends on (i, j, k) so no two faces share the same field
    // sample. This catches a missing face — a face that uses the wrong
    // snapshot will collapse to zero or to a value mismatching the
    // face's own coordinate dependence.
    for ((i, j, k), v) in p.ex.indexed_iter_mut() {
        *v = 0.1 + (i as f64) * 0.001 + (j as f64) * 0.01 + (k as f64) * 0.0001;
    }
    for ((i, j, k), v) in p.ey.indexed_iter_mut() {
        *v = 0.2 + (i as f64) * 0.001 + (j as f64) * 0.01 + (k as f64) * 0.0001;
    }
    for ((i, j, k), v) in p.ez.indexed_iter_mut() {
        *v = 0.3 + (i as f64) * 0.001 + (j as f64) * 0.01 + (k as f64) * 0.0001;
    }

    let lo = (4, 4, 4);
    let hi = (10, 10, 10);
    let mut region = SubgridRegion::new(&p, lo, hi).expect("valid subgrid");

    region.snapshot_coarse_e_t(&p);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.5);

    let fine = region.fine_grid();
    let probes = sample_all_face_e_t(fine);
    // No face's representative sample may be zero (which it would be if
    // the face were not written at all and the fine grid stayed at its
    // construction-time zero-initialisation).
    for (label, got) in probes {
        assert!(
            got.abs() > 1e-6,
            "face {label} appears unwritten by interpolate_coarse_e_to_fine: got {got}"
        );
    }
}

// -----------------------------------------------------------------------
// Helper: collect one representative sample per (face, tangential
// component) for the four-corner-and-one-centre probe set.
// -----------------------------------------------------------------------

fn sample_all_face_e_t(fine: &YeeGrid) -> Vec<(&'static str, f64)> {
    // Interior probe indices on each face: 1 (just inside the corner) is
    // safely inside the bracket-half domain and avoids any clamped
    // boundary cell.
    let nx = fine.nx;
    let ny = fine.ny;
    let nz = fine.nz;
    let (ix, jy, kz) = (1usize, 1usize, 1usize);

    vec![
        ("-x E_y", fine.ey[(0, jy, kz)]),
        ("-x E_z", fine.ez[(0, jy, kz)]),
        ("+x E_y", fine.ey[(nx, jy, kz)]),
        ("+x E_z", fine.ez[(nx, jy, kz)]),
        ("-y E_x", fine.ex[(ix, 0, kz)]),
        ("-y E_z", fine.ez[(ix, 0, kz)]),
        ("+y E_x", fine.ex[(ix, ny, kz)]),
        ("+y E_z", fine.ez[(ix, ny, kz)]),
        ("-z E_x", fine.ex[(ix, jy, 0)]),
        ("-z E_y", fine.ey[(ix, jy, 0)]),
        ("+z E_x", fine.ex[(ix, jy, nz)]),
        ("+z E_y", fine.ey[(ix, jy, nz)]),
    ]
}
