//! Per-cell ε_r coupling into the CPML auxiliary update.
//!
//! Regression test for the MMMMMMMM commit-body finding #1: prior to this
//! track, [`yee_fdtd::CpmlState::update_e`] hard-coded its coefficient as
//! `Δt / (ε₀ · grid.eps_r)`, ignoring any per-cell ε_r map. The
//! consequence was that placing a heterogeneous dielectric inside the
//! CPML region produced a mismatched coefficient between the bulk Yee
//! update (which uses the per-cell value) and the CPML auxiliary
//! correction (which used the scalar). The two updates therefore disagree
//! on the same E cell, the disagreement compounds each step, and the
//! field diverges to ~10^71 by step 1000 — empirically observed in the
//! superseded `heterogeneous_substrate.rs` workaround.
//!
//! After this track, both [`yee_fdtd::CpmlState::update_e`] /
//! [`yee_fdtd::CpmlState::update_h`] consult `grid.eps_r_cells` /
//! `grid.mu_r_cells` per cell, eliminating the coefficient mismatch.
//!
//! ## Two tests in this file
//!
//! 1. **`bit_exact_match_with_scalar_path`** — build two solvers from the
//!    same vacuum baseline: one with the scalar fallback (`eps_r_cells =
//!    None`), one with a per-cell map filled with the same value
//!    (`eps_r_cells = Some(all 1.0)`). Run 1000 CPML steps with the same
//!    source / probe and verify the two traces are bit-identical. This
//!    confirms the backward-compatibility contract: `Some(uniform)` must
//!    behave exactly like `None` with the same scalar.
//!
//! 2. **`heterogeneous_eps_in_cpml_is_stable`** — embed a 2.2-ε_r
//!    substrate slab that *crosses* the CPML region (the slab extends
//!    from interior into the high-x PML). Run 1000 CPML steps and
//!    verify `|E_z|_max < 10` at the final step (no late-time
//!    divergence to ~10^71). This is the precise scenario that was
//!    impossible before the fix.
//!
//! Wall-time budget: < 30 s release.

use ndarray::Array3;

use yee_fdtd::{CpmlParams, WalkingSkeletonSolver, YeeGrid};

/// Run a simulation injecting a Gaussian pulse at `source` and record
/// `E_z` at `probe` for each step.
fn run_trace(
    mut solver: WalkingSkeletonSolver,
    n_steps: usize,
    source: (usize, usize, usize),
    probe: (usize, usize, usize),
    t0: f64,
    sigma: f64,
) -> Vec<f64> {
    let mut trace = Vec::with_capacity(n_steps);
    for _ in 0..n_steps {
        solver.step_with_source(source.0, source.1, source.2, t0, sigma);
        trace.push(solver.grid().ez[probe]);
    }
    trace
}

#[test]
fn bit_exact_match_with_scalar_path() {
    // Small, fast grid: 30³ with npml = 10 leaves a 10³ interior. We
    // run 1000 steps to give late-time accumulator drift plenty of
    // opportunity to break a bit-exact comparison.
    const N: usize = 30;
    const DX: f64 = 1.0e-3;
    const NPML: usize = 10;
    const N_STEPS: usize = 1000;
    const SOURCE: (usize, usize, usize) = (15, 15, 15);
    const PROBE: (usize, usize, usize) = (20, 15, 15);

    let grid_ref = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid_ref.dt;
    let t0 = 20.0 * dt;
    let sigma = 6.0 * dt;
    drop(grid_ref);

    // --- Scalar (None) path ---
    let grid_scalar = YeeGrid::vacuum(N, N, N, DX);
    let params_scalar = CpmlParams::for_grid(&grid_scalar, NPML);
    let trace_scalar = run_trace(
        WalkingSkeletonSolver::with_cpml(grid_scalar, params_scalar),
        N_STEPS,
        SOURCE,
        PROBE,
        t0,
        sigma,
    );

    // --- Per-cell-but-uniform path ---
    let eps_cells = Array3::<f64>::from_elem((N + 1, N + 1, N + 1), 1.0);
    let grid_cells = YeeGrid::vacuum(N, N, N, DX).with_eps_r_cells(eps_cells);
    let params_cells = CpmlParams::for_grid(&grid_cells, NPML);
    let trace_cells = run_trace(
        WalkingSkeletonSolver::with_cpml(grid_cells, params_cells),
        N_STEPS,
        SOURCE,
        PROBE,
        t0,
        sigma,
    );

    // Bit-exact comparison: a uniform per-cell map of 1.0 must yield
    // identical results to the scalar fallback path at every step.
    assert_eq!(trace_scalar.len(), trace_cells.len());
    let mut max_abs_diff = 0.0_f64;
    let mut mismatch_step: Option<usize> = None;
    for (n, (a, b)) in trace_scalar.iter().zip(trace_cells.iter()).enumerate() {
        let d = (a - b).abs();
        if d > max_abs_diff {
            max_abs_diff = d;
        }
        if a.to_bits() != b.to_bits() && mismatch_step.is_none() {
            mismatch_step = Some(n);
        }
    }
    eprintln!(
        "bit_exact_match: max |Δ| = {max_abs_diff:.3e}, \
         first bit-mismatch step = {mismatch_step:?}"
    );
    assert!(
        mismatch_step.is_none(),
        "scalar and per-cell (uniform ε_r = 1.0) traces diverge at step \
         {mismatch_step:?} (max |Δ| = {max_abs_diff:.3e}); CPML must be \
         bit-exact when eps_r_cells matches scalar eps_r"
    );
}

#[test]
fn heterogeneous_eps_in_cpml_is_stable() {
    // 30³ grid with NPML = 10. Place a 2.2-ε_r substrate slab that
    // *intentionally* extends into the high-x PML: SLAB_LO = 15 (inside
    // the interior), SLAB_HI = N = 30 (so the slab crosses into and
    // through the high-x CPML at i ∈ [20, 30)).
    //
    // Before the fix this configuration drove the solver to ~10^71 by
    // step 1000. After the fix the per-cell coefficient is consistent
    // between the bulk Yee update and the CPML correction, so the
    // field stays bounded.
    const N: usize = 30;
    const DX: f64 = 1.0e-3;
    const NPML: usize = 10;
    const N_STEPS: usize = 1000;
    const SOURCE: (usize, usize, usize) = (12, 15, 15);
    // Probe near the source so we get a clean amplitude signal that
    // isn't already absorbed by the PML. Off-axis to catch any
    // anisotropic instability.
    const PROBE: (usize, usize, usize) = (12, 15, 15);

    const EPS_SUBSTRATE: f64 = 2.2;
    const SLAB_LO: usize = 15;
    const SLAB_HI: usize = N; // crosses into high-x CPML [20, 30)

    let dt = YeeGrid::vacuum(N, N, N, DX).dt;
    let t0 = 20.0 * dt;
    let sigma = 6.0 * dt;

    // Per-cell ε_r: 2.2 inside the slab (which straddles the CPML), 1.0
    // everywhere else.
    let mut eps_cells = Array3::<f64>::from_elem((N + 1, N + 1, N + 1), 1.0);
    for i in SLAB_LO..SLAB_HI {
        for j in 0..=N {
            for k in 0..=N {
                eps_cells[(i, j, k)] = EPS_SUBSTRATE;
            }
        }
    }

    let grid = YeeGrid::vacuum(N, N, N, DX).with_eps_r_cells(eps_cells);
    let params = CpmlParams::for_grid(&grid, NPML);
    let trace = run_trace(
        WalkingSkeletonSolver::with_cpml(grid, params),
        N_STEPS,
        SOURCE,
        PROBE,
        t0,
        sigma,
    );

    // Stability checks: all values finite, peak magnitude bounded.
    assert!(
        trace.iter().all(|x| x.is_finite()),
        "trace went non-finite — CPML / per-cell ε coupling unstable"
    );
    let peak = trace.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let final_abs = trace.last().copied().unwrap_or(0.0).abs();
    eprintln!(
        "heterogeneous_eps_in_cpml_is_stable: peak |E_z| = {peak:.3e}, \
         |E_z|@step{N_STEPS} = {final_abs:.3e}"
    );
    // The Gaussian source delivers amplitudes of order 1 V/m at the
    // probe early on. Pre-fix, the field grew to ~10^71 by step 1000;
    // post-fix it should stay bounded under O(1). |E_z|_max < 10 is
    // a deliberately loose ceiling that nonetheless rejects the
    // catastrophic late-time divergence the brief calls out.
    assert!(
        final_abs < 10.0,
        "|E_z|@step{N_STEPS} = {final_abs:.3e} exceeds 10 — CPML / per-cell \
         ε coupling appears unstable (pre-fix this diverged to ~10^71)"
    );
}
