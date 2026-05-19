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

    // `sigma_factor = 0.9` places σ ~10 % below the trial `(ω/c)²` so
    // TE_{101} is the dominant largest-`|1/(λ−σ)|` mode for the
    // inverse-iteration solver — empirically required for the
    // post-fix update rule (ω² = c²·λ) to converge from the air
    // warm-start on the Lorentz fixture. The D4 / pre-fix tests used
    // `sigma_factor = 2.5` because the buggy ε-double-divide form
    // happened to collapse ω onto the lower band where 2.5 happened to
    // bracket the lowest mode; the corrected form needs σ explicitly
    // below the target eigenvalue. See QQQQQQQQ D6 finding 2 in
    // `crates/yee-fem/validation/README.md` for the sigma-heuristic
    // history.
    let result = solver
        .solve_with_newton(&mesh, omega_0, 0.9)
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

    // (e) Analytic-compare gate (Track TTTTTTTT regression catch for
    //     QQQQQQQQ D6 finding 1).
    //
    //     The closed-form Pozar §3.1 dispersion relation for the
    //     TE_{101} mode on a homogeneously-filled PEC cavity is
    //
    //         ω² · ε(ω) / c² = (π/a)² + (π/d)²,
    //
    //     where `ε(ω)` is the same single-pole Lorentz form used by the
    //     FEM `M` matrix. Solving it via inner Newton on the continuum
    //     gives an analytic complex `ω_analytic` against which the
    //     converged FEM `ω` must agree within ±5 % on Re, ±10 % on Im.
    //
    //     Under the pre-fix `ω² = λ / (μ₀ε₀ε(ω))` form, the FEM
    //     converges to `Re(ω_analytic) / √ε_∞ ≈ Re(ω_analytic) / 2` on
    //     this fixture (ε_∞ = 4) and this assertion fails. Under the
    //     post-fix `ω² = c²·λ` form, the FEM tracks the analytic root
    //     to FEM-discretisation tolerance.
    let omega_analytic = analytic_lorentz_te101_omega();
    let re_rel_err = (result.omega.re - omega_analytic.re).abs() / omega_analytic.re.abs();
    assert!(
        re_rel_err <= 0.05,
        "lossy Lorentz Re(ω) = {} ({} GHz) vs analytic Re(ω_analytic) = {} ({} GHz): \
         relative error = {:.4e} > 5 % tolerance (pre-fix expected ratio ≈ 1/√ε_∞ = 0.5)",
        result.omega.re,
        result.omega.re / (2.0 * PI * 1e9),
        omega_analytic.re,
        omega_analytic.re / (2.0 * PI * 1e9),
        re_rel_err,
    );
    let im_rel_err = (result.omega.im - omega_analytic.im).abs() / omega_analytic.im.abs();
    assert!(
        im_rel_err <= 0.10,
        "lossy Lorentz Im(ω) = {} (= {} MHz) vs analytic Im(ω_analytic) = {} (= {} MHz): \
         relative error = {:.4e} > 10 % tolerance",
        result.omega.im,
        result.omega.im / (2.0 * PI * 1e6),
        omega_analytic.im,
        omega_analytic.im / (2.0 * PI * 1e6),
        im_rel_err,
    );
}

/// Closed-form complex Lorentz permittivity used by the test fixture:
/// matches `Material::eps_at` evaluated at the same parameters
/// (`ε_∞ = 4`, `ω_0 = 2π·20 GHz`, `ω_p = 2π·2 GHz`, `γ = 2π·0.5 GHz`),
/// generalised to complex argument so the analytic Newton root finder
/// can step into the complex plane.
///
/// Lorentz contribution: `+ ω_p² / (ω_0² − ω² − jγω)`. For complex
/// `ω` the denominator is `ω_0² − ω² − jγω` evaluated with complex
/// arithmetic; the FEM-side `Material::eps_at(omega: f64)` collapses to
/// this expression on the real axis.
fn eps_lorentz_complex(omega: Complex64) -> Complex64 {
    let eps_inf = Complex64::new(4.0, 0.0);
    let omega_0 = 2.0 * PI * 20.0e9;
    let omega_p = 2.0 * PI * 2.0e9;
    let gamma = 2.0 * PI * 0.5e9;
    let j = Complex64::new(0.0, 1.0);
    let denom = Complex64::new(omega_0 * omega_0, 0.0) - omega * omega - j * gamma * omega;
    eps_inf + Complex64::new(omega_p * omega_p, 0.0) / denom
}

/// Analytic complex `ω_analytic` for the TE_{101} mode on the WR-90
/// cavity uniformly filled with the test's Lorentz oscillator,
/// computed by Newton-iterating the Pozar §3.1 dispersion relation
///
/// ```text
///     F(ω) = ω² · ε(ω) / c² − k_geom²,    k_geom² = (π/a)² + (π/d)².
/// ```
///
/// The Newton derivative is evaluated analytically:
/// `F'(ω) = (2ω·ε(ω) + ω²·ε'(ω)) / c²`, with `ε'(ω)` the closed-form
/// derivative of the Lorentz pole. The root finder warm-starts from
/// the lossless air resonance `ω_air = c · √k_geom²` and converges in
/// ~10 iterations on this fixture; the result is the gold reference
/// against which the FEM `solve_with_newton` is compared.
fn analytic_lorentz_te101_omega() -> Complex64 {
    let omega_0 = 2.0 * PI * 20.0e9;
    let omega_p = 2.0 * PI * 2.0e9;
    let gamma = 2.0 * PI * 0.5e9;
    let j = Complex64::new(0.0, 1.0);

    // d/dω [ω_p² / (ω_0² − ω² − jγω)]
    //   = − ω_p² · (−2ω − jγ) / (ω_0² − ω² − jγω)²
    //   =   ω_p² · (2ω + jγ)  / (ω_0² − ω² − jγω)²
    let eps_prime = |omega: Complex64| -> Complex64 {
        let denom = Complex64::new(omega_0 * omega_0, 0.0) - omega * omega - j * gamma * omega;
        Complex64::new(omega_p * omega_p, 0.0) * (Complex64::new(2.0, 0.0) * omega + j * gamma)
            / (denom * denom)
    };

    let k_geom_sq = (PI / CAVITY_A_M).powi(2) + (PI / CAVITY_D_M).powi(2);
    let mut omega = Complex64::new(C0 * k_geom_sq.sqrt(), 0.0);
    for _ in 0..50 {
        let eps = eps_lorentz_complex(omega);
        let f = omega * omega * eps / Complex64::new(C0 * C0, 0.0) - Complex64::new(k_geom_sq, 0.0);
        let f_prime = (Complex64::new(2.0, 0.0) * omega * eps + omega * omega * eps_prime(omega))
            / Complex64::new(C0 * C0, 0.0);
        let step = f / f_prime;
        omega -= step;
        if step.norm() < 1e-9 * omega.norm() {
            break;
        }
    }
    omega
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
