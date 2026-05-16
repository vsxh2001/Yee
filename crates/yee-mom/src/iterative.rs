//! Iterative Krylov-subspace solvers for the MoM impedance matrix.
//!
//! Phase 1.6 walking-skeleton scope:
//! - Restarted GMRES(m) with default m = 30, max_restarts = 50.
//! - Diagonal Jacobi preconditioner (cheapest non-trivial preconditioner).
//!   Block-diagonal coming in Phase 1.6.1 when basis-block structure is
//!   exposed by `RwgBasis`.
//! - Convergence target: residual / b_norm <= 1e-10.
//!
//! Intended use case: n ≥ 50k, where the dense partial-pivot LU path used
//! since Phase 1.0 overflows GPU memory or wallclock budgets. For smaller
//! n the direct path remains preferable — GMRES has a higher constant
//! factor and only wins asymptotically.
//!
//! This commit lands only the public-surface skeleton ([`GmresParams`],
//! [`GmresResult`], and the [`gmres_jacobi`] entry point that returns
//! `x = 0` with `converged = false`). The Arnoldi / Givens inner loop is
//! filled in by the next commit.

use faer::Mat;
use num_complex::Complex64;

/// Tuning knobs for [`gmres_jacobi`].
#[derive(Debug, Clone, Copy)]
pub struct GmresParams {
    /// Restart parameter `m`: dimension of the Krylov subspace built between
    /// restarts. Larger `m` converges in fewer outer restarts but costs
    /// `O(m·n)` extra memory and `O(m²·n)` extra work per cycle.
    pub restart: usize,
    /// Maximum number of outer restart cycles before giving up.
    pub max_restarts: usize,
    /// Convergence threshold on the relative preconditioned residual
    /// `||M^{-1}(b - A x)||₂ / ||M^{-1} b||₂`.
    pub tolerance: f64,
}

impl Default for GmresParams {
    fn default() -> Self {
        Self {
            restart: 30,
            max_restarts: 50,
            tolerance: 1.0e-10,
        }
    }
}

/// Outcome of a [`gmres_jacobi`] call.
#[derive(Debug, Clone)]
pub struct GmresResult {
    /// Approximate solution `x ≈ A⁻¹ b`.
    pub x: Mat<Complex64>,
    /// Total number of Arnoldi iterations (summed across restart cycles).
    pub iterations: usize,
    /// Final relative preconditioned residual norm at termination.
    pub final_residual: f64,
    /// `true` iff `final_residual <= params.tolerance`.
    pub converged: bool,
}

/// Solve `A x = b` with restarted GMRES(m) and a left Jacobi (diagonal)
/// preconditioner.
///
/// Skeleton — returns the zero vector with `converged = false`. The real
/// Arnoldi + Givens implementation lands in the next commit. Kept here so
/// the public API surface (and downstream `LinearSolver` wiring) can be
/// committed first.
pub fn gmres_jacobi(
    a: faer::MatRef<Complex64>,
    _b: faer::MatRef<Complex64>,
    _params: GmresParams,
) -> GmresResult {
    let n = a.nrows();
    GmresResult {
        x: Mat::<Complex64>::zeros(n, 1),
        iterations: 0,
        final_residual: f64::INFINITY,
        converged: false,
    }
}
