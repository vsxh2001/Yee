//! Convolutional Perfectly Matched Layer (CPML) absorbing boundary.
//!
//! Implements the Roden & Gedney (2000) CPML formulation on all six outer
//! faces of a [`crate::grid::YeeGrid`]. References below are to:
//!
//! > J. A. Roden and S. D. Gedney, "Convolutional PML (CPML): An efficient
//! > FDTD implementation of the CFS-PML for arbitrary media",
//! > *Microwave Opt. Technol. Lett.* **27** (5) 334–339, 2000.
//!
//! ## Theory (R&G 2000, §III)
//!
//! The CFS-PML stretching variables are
//!
//! ```text
//! s_w = κ_w + σ_w / (α_w + jωε₀)        (R&G eq. 3)
//! ```
//!
//! In the time domain this becomes a convolution. The convolution is updated
//! recursively via the auxiliary variable ψ (R&G eq. 14):
//!
//! ```text
//! ψ_E_x(n+1) = b_x · ψ_E_x(n) + c_x · ∂H/∂x|^(n+1/2)
//! ```
//!
//! with coefficients (R&G eq. 25):
//!
//! ```text
//! b_w = exp[-(σ_w/κ_w + α_w) · Δt / ε₀]
//! c_w = (σ_w / (σ_w·κ_w + κ_w²·α_w)) · (b_w − 1)
//! ```
//!
//! The Maxwell curl update is then
//!
//! ```text
//! ∂E_z/∂t|cpml = (1/ε) · (1/κ_x · ∂H_y/∂x − 1/κ_y · ∂H_x/∂y + ψ_Ezx − ψ_Ezy)
//! ```
//!
//! ## Walking-skeleton scope (this commit)
//!
//! - [`CpmlParams`] with the standard parameter set:
//!   `σ_max = -(m+1) · ln(R_0) / (2·η₀·npml·dx)` with `R_0 = 1e-6`,
//!   `κ_max = 1`, `α_max = 0.05`, polynomial grading of order `m = 3`
//!   (R&G eq. 17 / Taflove §7.5).
//! - Helpers to sample the polynomial-grading profile (`σ(ρ)`, `κ(ρ)`,
//!   `α(ρ)`) and finalize the convolutional coefficients `(b, c)` from
//!   them (R&G eq. 25).
//!
//! The auxiliary-field state and update kernels follow in subsequent
//! commits; this module is intentionally small and unit-tested in
//! isolation first.

use yee_core::units::{EPS0, ETA0};

use crate::grid::YeeGrid;

/// CPML configuration parameters.
///
/// Defaults follow the standard Roden–Gedney 2000 / Taflove §7.5 recipe:
/// 10-cell layer, third-order polynomial grading, `R_0 = 1e-6`,
/// `κ_max = 1`, `α_max = 0.05`. `σ_max` is set so that the theoretical
/// reflection of a normally-incident wave at the inner PML edge is `R_0`.
#[derive(Debug, Clone, Copy)]
pub struct CpmlParams {
    /// Number of PML layers on each face (symmetric). Standard: 10.
    pub npml: usize,
    /// Polynomial grading order. Standard: 3.
    pub m: i32,
    /// Peak conductivity inside the PML. Populated from `R_0 = 1e-6` by
    /// [`CpmlParams::for_grid`]; this raw field is exposed for advanced
    /// callers who want to override the standard choice.
    pub sigma_max: f64,
    /// Peak coordinate-stretching factor. Standard: 1.0 (no stretching).
    pub kappa_max: f64,
    /// Peak CFS shift parameter. Standard: 0.05. Larger `α` improves
    /// low-frequency / evanescent-wave absorption at the cost of more
    /// reflection of propagating waves.
    pub alpha_max: f64,
}

impl Default for CpmlParams {
    fn default() -> Self {
        // sigma_max is grid-dependent; build a placeholder here and let
        // `for_grid` populate it. Default uses npml=10, m=3, dx=1mm.
        let npml = 10;
        let m = 3i32;
        let dx = 1.0e-3;
        Self {
            npml,
            m,
            sigma_max: sigma_max_optimal(m, npml, dx),
            kappa_max: 1.0,
            alpha_max: 0.05,
        }
    }
}

impl CpmlParams {
    /// Standard parameter set sized to the given grid.
    ///
    /// `σ_max = -(m+1) · ln(R_0) / (2·η₀·npml·dx)` with `R_0 = 1e-6`.
    pub fn for_grid(grid: &YeeGrid, npml: usize) -> Self {
        let m = 3i32;
        Self {
            npml,
            m,
            sigma_max: sigma_max_optimal(m, npml, grid.dx),
            kappa_max: 1.0,
            alpha_max: 0.05,
        }
    }
}

/// `σ_max = -(m+1) · ln(R_0) / (2·η₀·npml·dx)` with the standard
/// `R_0 = 1e-6` reflection target (Taflove eq. 7.66).
fn sigma_max_optimal(m: i32, npml: usize, dx: f64) -> f64 {
    let r0: f64 = 1.0e-6;
    -(f64::from(m) + 1.0) * r0.ln() / (2.0 * ETA0 * (npml as f64) * dx)
}

/// Sample σ, κ, α at depth `rho_over_d` ∈ (0, 1] using the standard
/// polynomial grading (R&G eq. 17 / Taflove eq. 7.79):
///
/// ```text
/// σ(ρ) = σ_max · (ρ/d)^m
/// κ(ρ) = 1 + (κ_max − 1) · (ρ/d)^m
/// α(ρ) = α_max · (1 − ρ/d)
/// ```
#[allow(dead_code)] // exercised by upcoming CpmlState commits + tests below
pub(crate) fn grading_sample(params: &CpmlParams, rho_over_d: f64) -> (f64, f64, f64) {
    let rho_m = rho_over_d.powi(params.m);
    let sigma = params.sigma_max * rho_m;
    let kappa = 1.0 + (params.kappa_max - 1.0) * rho_m;
    let alpha = params.alpha_max * (1.0 - rho_over_d);
    (sigma, kappa, alpha)
}

/// Finalize `(b, c)` from `(σ, κ, α)` and `dt` via R&G eq. 25:
///
/// ```text
/// b = exp(-(σ/κ + α) · dt / ε₀)
/// c = σ / (σ·κ + κ²·α) · (b − 1)
/// ```
#[allow(dead_code)] // exercised by upcoming CpmlState commits + tests below
pub(crate) fn finalize_coeffs(sigma: f64, kappa: f64, alpha: f64, dt: f64) -> (f64, f64) {
    let exponent = -(sigma / kappa + alpha) * dt / EPS0;
    let b = exponent.exp();
    let denom = sigma * kappa + kappa * kappa * alpha;
    let c = if denom.abs() > 1.0e-30 {
        sigma * (b - 1.0) / denom
    } else {
        0.0
    };
    (b, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_have_expected_shape() {
        let p = CpmlParams::default();
        assert_eq!(p.npml, 10);
        assert_eq!(p.m, 3);
        assert!(p.sigma_max > 0.0);
        assert_eq!(p.kappa_max, 1.0);
        assert!((p.alpha_max - 0.05).abs() < 1e-12);
    }

    #[test]
    fn sigma_max_optimal_matches_roden_gedney() {
        // For npml=10, dx=1mm, m=3, R_0=1e-6:
        // σ_max = -(4) · ln(1e-6) / (2 · 376.73 · 10 · 1e-3)
        //       ≈ -(4) · (-13.8155) / (7.5346) ≈ 7.334
        let s = sigma_max_optimal(3, 10, 1.0e-3);
        assert!(s > 7.0 && s < 8.0, "σ_max ≈ 7.33 expected, got {s}");
    }

    #[test]
    fn grading_sample_polynomial_endpoints() {
        let p = CpmlParams {
            npml: 10,
            m: 3,
            sigma_max: 10.0,
            kappa_max: 5.0,
            alpha_max: 0.1,
        };
        // At ρ/d = 0: σ = 0, κ = 1, α = α_max.
        let (s0, k0, a0) = grading_sample(&p, 0.0);
        assert!((s0 - 0.0).abs() < 1e-12);
        assert!((k0 - 1.0).abs() < 1e-12);
        assert!((a0 - 0.1).abs() < 1e-12);
        // At ρ/d = 1: σ = σ_max, κ = κ_max, α = 0.
        let (s1, k1, a1) = grading_sample(&p, 1.0);
        assert!((s1 - 10.0).abs() < 1e-12);
        assert!((k1 - 5.0).abs() < 1e-12);
        assert!((a1 - 0.0).abs() < 1e-12);
    }

    #[test]
    fn finalize_coeffs_in_neutral_limit() {
        // With σ = 0, b = exp(-α·dt/ε₀), c = 0 (degenerate denom).
        let (b, c) = finalize_coeffs(0.0, 1.0, 0.05, 1.0e-12);
        assert!(b > 0.0 && b <= 1.0);
        assert!((c - 0.0).abs() < 1e-12);
    }

    #[test]
    fn finalize_coeffs_strong_decay() {
        // With large σ, b is small and c is negative.
        let (b, c) = finalize_coeffs(10.0, 1.0, 0.0, 1.0e-12);
        assert!(b < 1.0 && b > 0.0);
        assert!(c < 0.0);
    }
}
