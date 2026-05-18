//! Phase 2.fdtd.7.0 Q2 — `SubgridRegion` scaffold tests.
//!
//! These tests pin the type surface (sizing, halved step sizes, error
//! conditions) so downstream tracks Q3 / Q4 / Q5 can be developed against a
//! stable API. The placeholder `SubgriddedSolver::step` must currently be
//! bit-for-bit identical to the wrapped `WalkingSkeletonSolver::step` —
//! the fine grid is dormant until Q5 lands the seven-stage interleave.

use yee_fdtd::{
    FdtdSolver, SubgridContext, SubgridRegion, SubgriddedSolver, WalkingSkeletonSolver, YeeGrid,
};

const N: usize = 16;
const DX: f64 = 1.0e-3;

fn parent_grid() -> YeeGrid {
    YeeGrid::vacuum(N, N, N, DX)
}

#[test]
fn subgrid_region_new_halves_step_sizes() {
    let parent = parent_grid();
    let region = SubgridRegion::new(&parent, (2, 2, 2), (6, 6, 6))
        .expect("valid subgrid bounds should construct");
    let fine = region.fine_grid();

    // 2× refinement: every step size (spatial + temporal) is exactly half.
    assert!((fine.dx - parent.dx / 2.0).abs() <= f64::EPSILON);
    assert!((fine.dy - parent.dy / 2.0).abs() <= f64::EPSILON);
    assert!((fine.dz - parent.dz / 2.0).abs() <= f64::EPSILON);
    assert!((fine.dt - parent.dt / 2.0).abs() <= f64::EPSILON);

    // Scalar materials inherit from the parent.
    assert_eq!(fine.eps_r, parent.eps_r);
    assert_eq!(fine.mu_r, parent.mu_r);
}

#[test]
fn subgrid_region_new_cell_count_doubled_per_axis() {
    let parent = parent_grid();
    // lo=(2,2,2), hi=(5,5,5) → coarse extent 3 per axis, fine extent 6 per axis.
    let region = SubgridRegion::new(&parent, (2, 2, 2), (5, 5, 5))
        .expect("valid subgrid bounds should construct");
    let fine = region.fine_grid();

    assert_eq!(fine.nx, 6);
    assert_eq!(fine.ny, 6);
    assert_eq!(fine.nz, 6);
}

#[test]
fn subgrid_region_new_rejects_lo_ge_hi() {
    let parent = parent_grid();
    // lo == hi on x.
    assert!(SubgridRegion::new(&parent, (4, 2, 2), (4, 5, 5)).is_err());
    // lo > hi on y.
    assert!(SubgridRegion::new(&parent, (2, 5, 2), (5, 4, 5)).is_err());
}

#[test]
fn subgrid_region_new_rejects_out_of_bounds_hi() {
    let parent = parent_grid();
    // hi.0 > parent.nx (16).
    assert!(SubgridRegion::new(&parent, (2, 2, 2), (17, 5, 5)).is_err());
}

#[test]
fn subgrid_region_overlaps_cpml_errors() {
    let parent = parent_grid();
    // npml = 4 on a 16-cell parent → CPML occupies [0, 4) and [12, 16) per
    // axis. A region at lo=(2, ..) with npml=4 overlaps the low-side CPML
    // and must error.
    let ctx = SubgridContext {
        cpml_thickness: Some(4),
        ..SubgridContext::default()
    };
    let result = SubgridRegion::new_with_context(&parent, (2, 5, 5), (8, 11, 11), ctx);
    assert!(result.is_err());

    // Sanity: a region fully inside the non-CPML interior succeeds.
    let interior_ok = SubgridRegion::new_with_context(&parent, (5, 5, 5), (11, 11, 11), ctx);
    assert!(interior_ok.is_ok());
}

#[test]
fn subgrid_region_overlaps_tfsf_face_errors() {
    let parent = parent_grid();
    // TF/SF box at coarse [4, 12) per axis. Region lo=(2,5,5) hi=(8,11,11)
    // crosses the low-x TF/SF face (parent x=4).
    let tfsf = ((4, 4, 4), (12, 12, 12));
    let ctx = SubgridContext {
        tfsf_box: Some(tfsf),
        ..SubgridContext::default()
    };
    let crossing = SubgridRegion::new_with_context(&parent, (2, 5, 5), (8, 11, 11), ctx);
    assert!(crossing.is_err());

    // Region wholly inside the TF/SF box is permitted (no face crossed).
    let inside = SubgridRegion::new_with_context(&parent, (5, 5, 5), (11, 11, 11), ctx);
    assert!(inside.is_ok());

    // Region wholly outside the TF/SF box is also permitted.
    let outside = SubgridRegion::new_with_context(&parent, (12, 5, 5), (15, 11, 11), ctx);
    assert!(outside.is_ok());
}

#[test]
fn subgridded_solver_step_matches_bare_walking_skeleton_at_step_10() {
    // Drive a bare `WalkingSkeletonSolver` and a `SubgriddedSolver::new(...)`
    // (no region) with the *same* Gaussian source for 10 steps. At Q2 the
    // wrapped step is a straight delegation, so the post-step E_z fields
    // must agree bit-for-bit (max abs diff == 0.0).
    let dt = parent_grid().dt;
    let t0 = 5.0 * dt;
    let sigma = 1.5 * dt;
    let source = (N / 2, N / 2, N / 2);
    let n_steps = 10;

    let mut bare = WalkingSkeletonSolver::new(parent_grid());
    let mut wrapped = SubgriddedSolver::new(WalkingSkeletonSolver::new(parent_grid()));

    for _ in 0..n_steps {
        let t_bare = bare.current_time();
        bare.update_h_only();
        bare.apply_cpml_h();
        bare.apply_gaussian_source_ez(source.0, source.1, source.2, t_bare, t0, sigma);
        bare.update_e_only();
        bare.apply_cpml_e();
        bare.advance_clock();

        let t_wrapped = wrapped.inner().current_time();
        // Drive the wrapped solver via its `inner_mut()` escape hatch to
        // keep the source-injection sequence identical to the bare path.
        wrapped.inner_mut().update_h_only();
        wrapped.inner_mut().apply_cpml_h();
        wrapped
            .inner_mut()
            .apply_gaussian_source_ez(source.0, source.1, source.2, t_wrapped, t0, sigma);
        wrapped.inner_mut().update_e_only();
        wrapped.inner_mut().apply_cpml_e();
        wrapped.inner_mut().advance_clock();
    }

    let g_bare = bare.grid();
    let g_wrap = wrapped.inner().grid();
    let mut max_diff = 0.0f64;
    for (a, b) in g_bare.ez.iter().zip(g_wrap.ez.iter()) {
        max_diff = max_diff.max((a - b).abs());
    }
    assert_eq!(
        max_diff, 0.0,
        "SubgriddedSolver placeholder must be bit-identical to bare WalkingSkeletonSolver at step 10"
    );
}

#[test]
fn subgridded_solver_step_placeholder_matches_inner_step() {
    // Independent of the source-driven check above: with no source and no
    // region, repeated calls to `SubgriddedSolver::step` must produce a
    // grid byte-identical to a bare `WalkingSkeletonSolver::step` over the
    // same number of steps.
    let mut bare = WalkingSkeletonSolver::new(parent_grid());
    let mut wrapped = SubgriddedSolver::new(WalkingSkeletonSolver::new(parent_grid()));

    // Seed a non-zero initial E_z so the leapfrog actually has something
    // to propagate (vacuum + zero initial conditions stays zero forever).
    wrapped.inner_mut().grid_mut().ez[[N / 2, N / 2, N / 2 - 1]] = 1.0;
    bare.grid_mut().ez[[N / 2, N / 2, N / 2 - 1]] = 1.0;

    for _ in 0..10 {
        bare.step();
        wrapped.step();
    }

    let g_bare = bare.grid();
    let g_wrap = wrapped.inner().grid();
    let mut max_diff = 0.0f64;
    for (a, b) in g_bare.ez.iter().zip(g_wrap.ez.iter()) {
        max_diff = max_diff.max((a - b).abs());
    }
    assert_eq!(max_diff, 0.0);
}

#[test]
fn subgrid_region_fine_grid_mut_is_writable() {
    let parent = parent_grid();
    let mut region = SubgridRegion::new(&parent, (2, 2, 2), (6, 6, 6)).expect("valid bounds");
    // Sanity: the mut getter actually returns a writable handle. This
    // pins the API surface for Q3 / Q4 which will write material and
    // boundary state through it.
    region.fine_grid_mut().ez[[1, 1, 1]] = 42.0;
    assert_eq!(region.fine_grid().ez[[1, 1, 1]], 42.0);
}
