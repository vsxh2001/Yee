//! Integration test for the FDTD walking skeleton.
//!
//! Drives the reference solver with a Gaussian-in-time `E_z` point source at
//! the center of a 50³ vacuum grid and checks that a wavefront reaches the
//! outer cells. The grid uses the hard PEC boundary, so this is *not* a test
//! of an absorbing boundary — only of propagation.

use yee_fdtd::{WalkingSkeletonSolver, YeeGrid};

#[test]
fn pulse_propagates_outward() {
    // 50³ cells of vacuum, 1 mm cell size.
    let grid = YeeGrid::vacuum(50, 50, 50, 1.0e-3);
    let mut solver = WalkingSkeletonSolver::new(grid);

    // ~3-cycle Gaussian centered at t0 = 5·dt, sigma = 1.5·dt.
    let dt = solver.grid().dt;
    let t0 = 5.0 * dt;
    let sigma = 1.5 * dt;

    // c·dt ≈ 0.9·dx/√3 ≈ 0.52 mm/step → 120 steps covers ~62 mm, well past
    // the 25-cell ≈ 25 mm distance to the boundary.
    for _ in 0..120 {
        solver.step_with_source(25, 25, 25, t0, sigma);
    }

    let center_energy: f64 = solver.grid().ez[(25, 25, 25)].powi(2);
    // Probe one cell inside the y = 0 PEC wall; the wall itself is clamped
    // to zero each step, so we sample the layer just interior to it.
    let edge_energy: f64 = (0..50).map(|i| solver.grid().ez[(i, 1, 25)].powi(2)).sum();

    println!("center energy: {center_energy:.3e}");
    println!("edge energy:   {edge_energy:.3e}");

    assert!(
        center_energy.is_finite(),
        "center energy diverged: {center_energy}"
    );
    assert!(
        edge_energy.is_finite(),
        "edge energy diverged: {edge_energy}"
    );
    assert!(
        edge_energy > 0.0,
        "wave should have propagated outward; got edge_energy = {edge_energy}"
    );
}

#[test]
fn courant_limit_is_positive() {
    let grid = YeeGrid::vacuum(10, 10, 10, 1.0e-3);
    let cfl = grid.courant_limit();
    assert!(
        cfl > 0.0 && cfl.is_finite(),
        "Courant limit must be finite > 0, got {cfl}"
    );
    assert!(
        grid.dt < cfl,
        "dt = {} should be below CFL = {}",
        grid.dt,
        cfl
    );
}
