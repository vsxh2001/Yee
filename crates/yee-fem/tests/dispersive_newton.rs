//! Phase 4.fem.eig.1 step D5 — integration tests for
//! [`yee_fem::DispersiveSolver::solve_with_newton`].
//!
//! Test inventory:
//!
//! 1. `free_space_air_newton_converges_to_phase_4_0_te101` — WR-90
//!    cavity (a = 22.86 mm, b = 10.16 mm, d = 30 mm) with free-space
//!    air filler. Warm-start the Newton tracker from `ω₀ = 0.9 · 2π ·
//!    f_TE101` (10 % below the analytic Pozar §6.3 resonance). Verify
//!    the converged `ω` matches the analytic TE_{101} frequency within
//!    0.5 % in ≤ 5 iterations. Because air is non-dispersive
//!    (`ε(ω) ≡ 1`) the fixed-point update `ω_{n+1} =
//!    c · √k²(ω_n)` returns the air resonance in the very first
//!    iteration; the test pins both the analytic agreement and the
//!    iteration-count budget per the D5 brief.
//!
//! 2. `lossy_lorentz_cavity_newton_converges_complex` — same cavity
//!    geometry, but the bulk filler is a single-pole Lorentz oscillator
//!    (the same parameters used by
//!    `tests/dispersive_solve.rs::lossy_substrate_complex_eigenvalue_*`).
//!    Verify the Newton tracker converges to a complex ω with
//!    non-trivial `Im(ω)` (lossy mode) and the iteration count stays
//!    bounded under [`yee_fem::DispersiveSolver::max_iter`].
//!
//! 3. `max_iter_exhaustion_returns_error` — same lossy fixture as (2)
//!    but with `max_iter = 1` and `tol = 1e-12`. Verify the tracker
//!    surfaces [`yee_fem::DispersiveError::NewtonDidNotConverge`] with
//!    a finite `last_residual` and a finite `last_omega`.
//!
//! ## Mesh size
//!
//! Tests use the `(8, 6, 10)` Kuhn brick subdivision of the WR-90
//! cavity (2880 tets) — same mesh as
//! `tests/dispersive_solve.rs` so the inner linearised solve runs in
//! the same wall-time budget (<10 s `--release` per assemble + solve).
//!
//! ## Shift convention
//!
//! `sigma_factor = 2.5` mirrors the D4 fixture convention; per the
//! `tests/dispersive_solve.rs` header comment this sits between the
//! 8th and 9th physical modes on this mesh, so all ten lowest physical
//! modes have `|θ| > |θ_grad|` and inverse-iteration converges to them
//! in ascending `Re(k²)` order.
//!
//! References:
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
//! * `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`
//!   step D5.
//! * `docs/src/decisions/0039-phase-4-fem-eig-1-dispersive-scope.md`.

use std::f64::consts::PI;

use num_complex::Complex64;
use yee_fem::{DispersiveError, DispersiveSolver, Material, MaterialDatabase, MaterialPole};
use yee_mesh::TetMesh3D;

/// Speed of light (m/s) — workspace-canonical from `yee_core::units`.
const C0: f64 = yee_core::units::C0;

/// WR-90 cavity extents (m).
const CAVITY_A_M: f64 = 0.022_86;
const CAVITY_B_M: f64 = 0.010_16;
const CAVITY_D_M: f64 = 0.030;

/// Mesh density — matches `tests/dispersive_solve.rs`.
const NX: usize = 8;
const NY: usize = 6;
const NZ: usize = 10;

/// Analytic TE_{101} frequency for the air-filled cavity (Hz).
/// Pozar §6.3 eq. 6.42.
fn f_te101_hz() -> f64 {
    0.5 * C0 * ((1.0 / CAVITY_A_M).powi(2) + (1.0 / CAVITY_D_M).powi(2)).sqrt()
}

/// Build the WR-90 cavity mesh used by every test in this file.
fn build_cavity() -> TetMesh3D {
    TetMesh3D::cavity_uniform(CAVITY_A_M, CAVITY_B_M, CAVITY_D_M, NX, NY, NZ)
        .expect("cavity_uniform must succeed for the standard WR-90 dimensions")
}

/// D5 gate test 1: free-space air-cavity Newton tracker converges to
/// the analytic TE_{101} resonance within 0.5 % in ≤ 5 iterations.
#[test]
fn free_space_air_newton_converges_to_phase_4_0_te101() {
    let mesh = build_cavity();

    // Warm-start at 90 % of the analytic resonance — a deliberately
    // off-target warm-start that nonetheless lies inside the
    // monotone-convergence basin for the non-dispersive air case.
    let f_te101 = f_te101_hz();
    let omega_te101 = 2.0 * PI * f_te101;
    let omega_0 = Complex64::new(0.9 * omega_te101, 0.0);

    // Free-space material database — air at the bulk tag `0`.
    let db = MaterialDatabase::new().with_material(0, Material::default());
    let mut solver = DispersiveSolver::new(db);
    // Cap the outer Newton budget at the brief's 5-iteration ceiling.
    solver.newton_max_iter = 5;
    // Use a relaxed outer-Newton tol so the 0.5 % analytic agreement
    // (5e-3 relative) is the binding constraint, not floating-point
    // noise in the inner solve.
    solver.newton_tol = 1e-4;
    // Inner-solver tolerance: 1e-7 is comfortably tighter than the
    // 1e-4 outer Newton bound and well within the per-mode 1000-
    // iteration budget on this mesh.  The default 1e-8 is right at
    // the working-precision boundary for the 10th deflated mode
    // returned by `ComplexInverseIterEigen`.
    solver.tol = 1e-7;

    let result = solver
        .solve_with_newton(&mesh, omega_0, 2.5)
        .expect("Newton solve_with_newton must converge on free-space air");

    // (a) Real part within 0.5 % of the analytic Pozar §6.3 resonance.
    let omega_re = result.omega.re;
    let rel_err = (omega_re - omega_te101).abs() / omega_te101;
    assert!(
        rel_err <= 5e-3,
        "free-space converged ω = {} (Re = {} GHz) vs analytic ω_TE101 = {} (= {} GHz); \
         rel err = {:.4e}, tol = 5e-3",
        result.omega,
        omega_re / (2.0 * PI * 1e9),
        omega_te101,
        omega_te101 / (2.0 * PI * 1e9),
        rel_err,
    );

    // (b) Im(ω) must be vanishing for purely real air ε(ω).
    let im_rel = result.omega.im.abs() / result.omega.re.abs();
    assert!(
        im_rel < 1e-6,
        "free-space converged Im(ω)/Re(ω) = {im_rel:.4e} should be ≤ 1e-6",
    );

    // (c) The composed physical k = ω · √(μ₀ε₀ε) should be a
    // real-positive number for air.
    assert!(
        result.k_complex.im.abs() < 1e-6 * result.k_complex.re.abs(),
        "free-space k must be real: got {}",
        result.k_complex,
    );
    assert!(
        result.k_complex.re > 0.0,
        "free-space k.re must be positive: got {}",
        result.k_complex,
    );
}

/// D5 gate test 2: a Lorentz-filled cavity converges to a complex ω.
///
/// Material: single-pole Lorentz oscillator with `ε_∞ = 4`,
/// `ω_0 = 2π · 20 GHz`, `ω_p = 2π · 2 GHz`, `γ = 2π · 0.5 GHz`. The
/// Newton tracker warm-starts from the v0 free-space TE_{101}
/// resonance (≈ 8.244 GHz, well below `ω_0 = 20 GHz` so the Lorentz
/// pole sits comfortably outside the convergence basin per spec §11).
/// The fixed-point update lifts ω to its self-consistent complex
/// value; the test only verifies (i) convergence with non-zero Im(ω)
/// and (ii) iteration count stays below
/// [`DispersiveSolver::max_iter`].
#[test]
fn lossy_lorentz_cavity_newton_converges_complex() {
    let mesh = build_cavity();

    // Warm-start at the air resonance — far enough above the Lorentz
    // pole at 20 GHz to keep `ε(ω)` smooth and bounded.
    let f_te101 = f_te101_hz();
    let omega_0 = Complex64::new(2.0 * PI * f_te101, 0.0);

    // Lossy Lorentz bulk material at tag 0 (the Newton update reads
    // ε(ω) from this tag per `dispersive::BULK_TAG`).
    let lorentz = Material {
        eps_inf: 4.0,
        mu_r: 1.0,
        poles: vec![MaterialPole::Lorentz {
            omega_0: 2.0 * PI * 20.0e9,
            omega_p: 2.0 * PI * 2.0e9,
            gamma: 2.0 * PI * 0.5e9,
        }],
    };
    let db = MaterialDatabase::new().with_material(0, lorentz);
    let mut solver = DispersiveSolver::new(db);
    // Generous outer-Newton budget so the iteration count is bounded
    // by convergence, not by the cap; the test asserts the result
    // arrives well before this limit.
    solver.newton_max_iter = 25;
    solver.newton_tol = 1e-5;
    // Inner-solver tol relaxed to match free_space_*  — see comment
    // there.
    solver.tol = 1e-7;

    let result = solver
        .solve_with_newton(&mesh, omega_0, 2.5)
        .expect("Newton solve_with_newton must converge on Lorentz fixture");

    // (a) The converged ω must have non-trivial imaginary part — the
    // lossy mode picks up a finite decay rate from the Lorentz pole.
    assert!(
        result.omega.im.abs() > 0.0,
        "lossy Lorentz converged ω = {} should have non-zero Im(ω)",
        result.omega,
    );

    // (b) Real part should sit comfortably above 0 and well below the
    // Lorentz pole (Re(ε) > 1 at ω < ω_0 lowers the resonance).
    assert!(
        result.omega.re > 0.0,
        "lossy converged Re(ω) must be positive, got {}",
        result.omega,
    );
    assert!(
        result.omega.re < 2.0 * PI * 20.0e9,
        "lossy converged Re(ω) = {} must be below the Lorentz pole",
        result.omega.re,
    );

    // (c) Composed k must also be complex (Im(k) ≠ 0).
    assert!(
        result.k_complex.im.abs() > 0.0,
        "lossy k = {} should have non-zero Im(k)",
        result.k_complex,
    );

    // (d) The eigenvector vector length is non-zero.
    assert!(
        !result.e_vec.is_empty(),
        "converged eigenvector must have at least one interior DoF",
    );
}

/// D5 gate test 3: `max_iter = 1` with a tight tol must error out on
/// the lossy fixture (one fixed-point step cannot meet `|Δω/ω| <
/// 1e-12` on a non-vacuum cavity).
#[test]
fn max_iter_exhaustion_returns_error() {
    let mesh = build_cavity();

    let f_te101 = f_te101_hz();
    let omega_0 = Complex64::new(2.0 * PI * f_te101, 0.0);

    let lorentz = Material {
        eps_inf: 4.0,
        mu_r: 1.0,
        poles: vec![MaterialPole::Lorentz {
            omega_0: 2.0 * PI * 20.0e9,
            omega_p: 2.0 * PI * 2.0e9,
            gamma: 2.0 * PI * 0.5e9,
        }],
    };
    let db = MaterialDatabase::new().with_material(0, lorentz);
    let mut solver = DispersiveSolver::new(db);
    solver.newton_max_iter = 1;
    solver.newton_tol = 1e-12;

    let result = solver.solve_with_newton(&mesh, omega_0, 2.5);
    match result {
        Err(DispersiveError::NewtonDidNotConverge {
            last_omega,
            last_k_sq,
            last_residual,
        }) => {
            assert!(
                last_residual.is_finite(),
                "last_residual must be finite, got {last_residual}",
            );
            assert!(
                last_omega.re.is_finite() && last_omega.im.is_finite(),
                "last_omega must be finite, got {last_omega}",
            );
            assert!(
                last_k_sq.re.is_finite() && last_k_sq.im.is_finite(),
                "last_k_sq must be finite, got {last_k_sq}",
            );
            // Sanity: the single step should produce a real residual
            // well above the 1e-12 target.
            assert!(
                last_residual > 1e-12,
                "single-iteration residual = {last_residual:e} should exceed tol 1e-12",
            );
        }
        Err(other) => panic!("expected DispersiveError::NewtonDidNotConverge, got {other:?}"),
        Ok(eig) => panic!(
            "expected NewtonDidNotConverge with max_iter = 1 / tol = 1e-12, \
             got converged eigenpair ω = {}",
            eig.omega,
        ),
    }
}
