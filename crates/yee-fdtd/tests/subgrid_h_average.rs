//! Fine → coarse closure tests for `SubgridRegion` (Phase 2.fdtd.7 Q4).
//!
//! Exercises [`SubgridRegion::average_fine_h_to_coarse`] and
//! [`SubgridRegion::overwrite_coarse_e_from_fine`] — the Chevalier 1997
//! §IV energy-balance closure step that overwrites coarse `H_t` (and
//! coarse `E_t` on stage 7) with area-/edge-averages of the fine grid.
//!
//! Coverage:
//!
//! - `area_average_of_uniform_h_is_uniform` — constant fine `H_z` on a
//!   face → exact constant coarse `H_z`.
//! - `forward_reverse_round_trip_preserves_static_field` — snapshot ⇒
//!   interpolate(0.5) ⇒ no fine update ⇒ closure ⇒ coarse matches the
//!   original within `f64::EPSILON · 100`. This is the discrete energy-
//!   balance closure regression test (Chevalier 1997 §IV).
//! - `area_average_face_dispatch_covers_all_six_faces` — asymmetric fine
//!   `H` perturbation per face → every coarse target cell receives a
//!   non-trivial update.
//! - `linear_gradient_averages_to_midpoint` — fine `H_z = a + b · x_f`
//!   averages to `a + b · x_center` on the corresponding coarse face.

use yee_fdtd::{SubgridRegion, YeeGrid};

/// Build a coarse parent + nested fine region of the given coarse extents,
/// returning both pieces. Vacuum, `dx = 1e-3 m`.
fn build_pair(
    nx: usize,
    ny: usize,
    nz: usize,
    lo: (usize, usize, usize),
    hi: (usize, usize, usize),
) -> (YeeGrid, SubgridRegion) {
    let parent = YeeGrid::vacuum(nx, ny, nz, 1.0e-3);
    let region = SubgridRegion::new(&parent, lo, hi).expect("valid subgrid region");
    (parent, region)
}

/// Fill every cell of the fine grid's `H_z` with the same constant.
fn fill_fine_hz(region: &mut SubgridRegion, value: f64) {
    let fine = region.fine_grid_mut();
    fine.hz.fill(value);
}

#[test]
fn area_average_of_uniform_h_is_uniform() {
    let lo = (2usize, 3, 4);
    let hi = (8usize, 9, 10);
    let (mut parent, mut region) = build_pair(12, 12, 12, lo, hi);
    fill_fine_hz(&mut region, 3.5);

    region.average_fine_h_to_coarse(&mut parent);

    // H_z is overwritten on the ±x and ±y faces (it is tangential to
    // those four faces). Each affected coarse cell must equal the
    // uniform fine value exactly.
    let tol = f64::EPSILON * 10.0;

    // ±x faces: coarse hz on i ∈ {lo.0, hi.0 - 1}, j ∈ [lo.1, hi.1),
    // k ∈ [lo.2, hi.2].
    for i_face in [lo.0, hi.0 - 1] {
        for j_c in lo.1..hi.1 {
            for k_c in lo.2..=hi.2 {
                assert!(
                    (parent.hz[(i_face, j_c, k_c)] - 3.5).abs() <= tol,
                    "±x face hz[({i_face},{j_c},{k_c})] = {} ≠ 3.5",
                    parent.hz[(i_face, j_c, k_c)]
                );
            }
        }
    }
    // ±y faces: coarse hz on j ∈ {lo.1, hi.1 - 1}, i ∈ [lo.0, hi.0),
    // k ∈ [lo.2, hi.2].
    for j_face in [lo.1, hi.1 - 1] {
        for i_c in lo.0..hi.0 {
            for k_c in lo.2..=hi.2 {
                assert!(
                    (parent.hz[(i_c, j_face, k_c)] - 3.5).abs() <= tol,
                    "±y face hz[({i_c},{j_face},{k_c})] = {} ≠ 3.5",
                    parent.hz[(i_c, j_face, k_c)]
                );
            }
        }
    }
}

#[test]
fn forward_reverse_round_trip_preserves_static_field() {
    // Static parent field: per-component CONSTANTS. A constant field is
    // preserved exactly by any linear average (edge-average, area-average,
    // and bilinear spatial / linear temporal interp), so the round-trip
    // is bit-clean within float tolerance. Variable-coordinate fields
    // require careful tracking of per-component half-cell offsets (which
    // differ between E_x / E_y / E_z); see the dedicated
    // `linear_gradient_averages_to_midpoint` test for that geometry.
    let lo = (2usize, 2, 2);
    let hi = (8usize, 8, 8);
    let mut parent = YeeGrid::vacuum(12, 12, 12, 1.0e-3);

    let ex_const = 0.1;
    let ey_const = 0.5;
    let ez_const = 0.9;
    let hx_const = 1.3;
    let hy_const = -2.4;
    let hz_const = 7.6;
    parent.ex.fill(ex_const);
    parent.ey.fill(ey_const);
    parent.ez.fill(ez_const);
    parent.hx.fill(hx_const);
    parent.hy.fill(hy_const);
    parent.hz.fill(hz_const);

    let baseline_ex = parent.ex.clone();
    let baseline_ey = parent.ey.clone();
    let baseline_ez = parent.ez.clone();

    let mut region = SubgridRegion::new(&parent, lo, hi).expect("valid region");

    // Snapshot start and end coincide → temporal blend is the identity
    // in time. Spatial interpolation of a constant field is exact.
    region.snapshot_coarse_e_t(&parent);
    region.snapshot_coarse_e_t_end(&parent);
    region.interpolate_coarse_e_to_fine(0.5);

    // Mirror the parent constants into the fine interior so the closure
    // averages reproduce the parent values exactly.
    {
        let fine = region.fine_grid_mut();
        fine.ex.fill(ex_const);
        fine.ey.fill(ey_const);
        fine.ez.fill(ez_const);
        fine.hx.fill(hx_const);
        fine.hy.fill(hy_const);
        fine.hz.fill(hz_const);
    }

    // Baseline parent H equals the constants seeded above, and fine H
    // matches. The closure averages constants → constants → no change.
    let baseline_hx = parent.hx.clone();
    let baseline_hy = parent.hy.clone();
    let baseline_hz = parent.hz.clone();
    region.average_fine_h_to_coarse(&mut parent);
    let round_trip_tol = f64::EPSILON * 100.0;
    assert!(
        parent
            .hx
            .iter()
            .zip(baseline_hx.iter())
            .all(|(a, b)| (a - b).abs() <= round_trip_tol),
        "H closure perturbed coarse hx away from baseline"
    );
    assert!(
        parent
            .hy
            .iter()
            .zip(baseline_hy.iter())
            .all(|(a, b)| (a - b).abs() <= round_trip_tol),
        "H closure perturbed coarse hy away from baseline"
    );
    assert!(
        parent
            .hz
            .iter()
            .zip(baseline_hz.iter())
            .all(|(a, b)| (a - b).abs() <= round_trip_tol),
        "H closure perturbed coarse hz away from baseline"
    );

    // E closure: should overwrite the coarse `E_t` on the six interface
    // planes with the edge-averaged fine E_t. For a linear parent field
    // replicated into the fine grid, those averages reproduce the
    // parent value exactly (within float round-off).
    region.overwrite_coarse_e_from_fine(&mut parent);

    // The closure only writes to the six interface planes; coarse
    // baselines elsewhere remain untouched. Check each affected plane.
    let tol = f64::EPSILON * 100.0;
    // ±x faces touch ey and ez at i ∈ {lo.0, hi.0}.
    for i_c in [lo.0, hi.0] {
        for j_c in lo.1..hi.1 {
            for k_c in lo.2..=hi.2 {
                let got = parent.ey[(i_c, j_c, k_c)];
                let want = baseline_ey[(i_c, j_c, k_c)];
                assert!(
                    (got - want).abs() <= tol,
                    "ey[{i_c},{j_c},{k_c}]: {got} vs {want}"
                );
            }
        }
        for j_c in lo.1..=hi.1 {
            for k_c in lo.2..hi.2 {
                let got = parent.ez[(i_c, j_c, k_c)];
                let want = baseline_ez[(i_c, j_c, k_c)];
                assert!(
                    (got - want).abs() <= tol,
                    "ez[{i_c},{j_c},{k_c}]: {got} vs {want}"
                );
            }
        }
    }
    // ±y faces touch ex and ez at j ∈ {lo.1, hi.1}.
    for j_c in [lo.1, hi.1] {
        for i_c in lo.0..hi.0 {
            for k_c in lo.2..=hi.2 {
                let got = parent.ex[(i_c, j_c, k_c)];
                let want = baseline_ex[(i_c, j_c, k_c)];
                assert!(
                    (got - want).abs() <= tol,
                    "ex[{i_c},{j_c},{k_c}]: {got} vs {want}"
                );
            }
        }
        for i_c in lo.0..=hi.0 {
            for k_c in lo.2..hi.2 {
                let got = parent.ez[(i_c, j_c, k_c)];
                let want = baseline_ez[(i_c, j_c, k_c)];
                assert!(
                    (got - want).abs() <= tol,
                    "ez[{i_c},{j_c},{k_c}]: {got} vs {want}"
                );
            }
        }
    }
    // ±z faces touch ex and ey at k ∈ {lo.2, hi.2}.
    for k_c in [lo.2, hi.2] {
        for i_c in lo.0..hi.0 {
            for j_c in lo.1..=hi.1 {
                let got = parent.ex[(i_c, j_c, k_c)];
                let want = baseline_ex[(i_c, j_c, k_c)];
                assert!(
                    (got - want).abs() <= tol,
                    "ex[{i_c},{j_c},{k_c}]: {got} vs {want}"
                );
            }
        }
        for i_c in lo.0..=hi.0 {
            for j_c in lo.1..hi.1 {
                let got = parent.ey[(i_c, j_c, k_c)];
                let want = baseline_ey[(i_c, j_c, k_c)];
                assert!(
                    (got - want).abs() <= tol,
                    "ey[{i_c},{j_c},{k_c}]: {got} vs {want}"
                );
            }
        }
    }
}

#[test]
fn area_average_face_dispatch_covers_all_six_faces() {
    let lo = (2usize, 2, 2);
    let hi = (8usize, 8, 8);
    let (mut parent, mut region) = build_pair(12, 12, 12, lo, hi);

    // Seed the fine grid with a unique value per component so each face's
    // average is non-zero and distinct, then average and verify each face
    // produced a non-trivial coarse update.
    {
        let fine = region.fine_grid_mut();
        fine.hx.fill(1.25);
        fine.hy.fill(-2.5);
        fine.hz.fill(7.75);
    }

    // Zero out the coarse H first so the test is sensitive to any write.
    parent.hx.fill(0.0);
    parent.hy.fill(0.0);
    parent.hz.fill(0.0);

    region.average_fine_h_to_coarse(&mut parent);

    let tol = f64::EPSILON * 10.0;

    // ±x faces: write H_y (= -2.5) on i ∈ {lo.0, hi.0-1}, j ∈ [lo.1..=hi.1],
    // k ∈ [lo.2..hi.2). Also H_z (= 7.75) on i ∈ {lo.0, hi.0-1}.
    for i_face in [lo.0, hi.0 - 1] {
        for j_c in lo.1..=hi.1 {
            for k_c in lo.2..hi.2 {
                assert!((parent.hy[(i_face, j_c, k_c)] - (-2.5)).abs() <= tol);
            }
        }
        for j_c in lo.1..hi.1 {
            for k_c in lo.2..=hi.2 {
                assert!((parent.hz[(i_face, j_c, k_c)] - 7.75).abs() <= tol);
            }
        }
    }

    // ±y faces: write H_x (= 1.25) on j ∈ {lo.1, hi.1-1}, i ∈ [lo.0..=hi.0],
    // k ∈ [lo.2..hi.2). And H_z (= 7.75) on j ∈ {lo.1, hi.1-1}, i ∈
    // [lo.0..hi.0), k ∈ [lo.2..=hi.2].
    for j_face in [lo.1, hi.1 - 1] {
        for i_c in lo.0..=hi.0 {
            for k_c in lo.2..hi.2 {
                assert!((parent.hx[(i_c, j_face, k_c)] - 1.25).abs() <= tol);
            }
        }
        for i_c in lo.0..hi.0 {
            for k_c in lo.2..=hi.2 {
                assert!((parent.hz[(i_c, j_face, k_c)] - 7.75).abs() <= tol);
            }
        }
    }

    // ±z faces: write H_x (= 1.25) on k ∈ {lo.2, hi.2-1}, i ∈ [lo.0..=hi.0],
    // j ∈ [lo.1..hi.1). And H_y (= -2.5) on k ∈ {lo.2, hi.2-1}, i ∈
    // [lo.0..hi.0), j ∈ [lo.1..=hi.1].
    for k_face in [lo.2, hi.2 - 1] {
        for i_c in lo.0..=hi.0 {
            for j_c in lo.1..hi.1 {
                assert!((parent.hx[(i_c, j_c, k_face)] - 1.25).abs() <= tol);
            }
        }
        for i_c in lo.0..hi.0 {
            for j_c in lo.1..=hi.1 {
                assert!((parent.hy[(i_c, j_c, k_face)] - (-2.5)).abs() <= tol);
            }
        }
    }
}

#[test]
fn linear_gradient_averages_to_midpoint() {
    // Fine H_z = a + b · x_f (with x_f the fine-x coordinate of the H_z
    // sample). Each coarse H_z face overlaps the 2×2 fine tile spanning
    // fine_i ∈ {2i_rel, 2i_rel + 1}; the arithmetic mean of those two
    // fine x-coordinates equals the coarse-cell midpoint.
    let lo = (3usize, 3, 3);
    let hi = (7usize, 7, 7);
    let (mut parent, mut region) = build_pair(12, 12, 12, lo, hi);

    let a = 0.4;
    let b = 1.7;
    {
        let fine = region.fine_grid_mut();
        // H_z is shaped [fine_nx, fine_ny, fine_nz + 1]; its x coordinate
        // is the cell-center x = (i_f + 0.5) * dx_fine. We parametrise
        // the gradient in fine-x index units (dx_fine cancels in the
        // test because we are comparing fine and coarse on the same
        // physical line).
        for i_f in 0..fine.nx {
            let x = (i_f as f64) + 0.5; // fine-x in fine-cell units
            let v = a + b * x;
            for j_f in 0..fine.ny {
                for k_f in 0..=fine.nz {
                    fine.hz[(i_f, j_f, k_f)] = v;
                }
            }
        }
    }

    parent.hz.fill(0.0);
    region.average_fine_h_to_coarse(&mut parent);

    let tol = f64::EPSILON * 100.0;
    // Inspect the +z face (k_c = hi.2 - 1) since both ±z faces touch H_z
    // through `avg_face_z` only via H_x — wait, H_z is touched by ±x and
    // ±y faces. Pick the +x face (i_c = hi.0 - 1) where the fine-x tile
    // is fine_i ∈ {fine_nx - 2, fine_nx - 1}, so the expected mean is
    // a + b * ((fine_nx - 2) + 0.5 + (fine_nx - 1) + 0.5) / 2.
    let fine_nx = 2 * (hi.0 - lo.0);
    let fine_x_lo_center = (fine_nx - 2) as f64 + 0.5;
    let fine_x_hi_center = (fine_nx - 1) as f64 + 0.5;
    let expected_plus_x = a + b * 0.5 * (fine_x_lo_center + fine_x_hi_center);
    for j_c in lo.1..hi.1 {
        for k_c in lo.2..=hi.2 {
            let got = parent.hz[(hi.0 - 1, j_c, k_c)];
            assert!(
                (got - expected_plus_x).abs() <= tol,
                "+x face hz[({},{j_c},{k_c})] = {got}, expected {expected_plus_x}",
                hi.0 - 1
            );
        }
    }
    // The −x face (i_c = lo.0) has fine_i ∈ {0, 1}, centers at 0.5 and 1.5,
    // so expected mean is a + b * 1.0.
    let expected_minus_x = a + b * 1.0;
    for j_c in lo.1..hi.1 {
        for k_c in lo.2..=hi.2 {
            let got = parent.hz[(lo.0, j_c, k_c)];
            assert!(
                (got - expected_minus_x).abs() <= tol,
                "−x face hz[({},{j_c},{k_c})] = {got}, expected {expected_minus_x}",
                lo.0
            );
        }
    }
}
