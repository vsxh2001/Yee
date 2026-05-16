//! Bayesian optimization with a Gaussian-process surrogate and Expected
//! Improvement acquisition.
//!
//! ## Algorithm
//!
//! For a black-box objective `f: ℝ^d → ℝ` defined on a hyper-rectangle:
//!
//! 1. Sample `n_initial` Latin-hypercube points in the bounds, evaluate the
//!    objective at each.
//! 2. Repeat `n_iters` times:
//!    a. Refit a GP via [`crate::GaussianProcess::fit_ml`] on the running
//!    history.
//!    b. Draw `n_candidates` uniform random points in the bounds.
//!    c. Score each candidate by Expected Improvement against the current
//!    `f_best`, pick the maximizer, evaluate the objective there, append
//!    the new `(x, y)` pair.
//! 3. Return the best `(x, y)` seen and the full history.
//!
//! ## Expected Improvement
//!
//! For minimization with current best `f_best`, predictive mean `μ`,
//! predictive stddev `σ`, exploration parameter `ξ`:
//!
//! ```text
//! improvement = f_best - μ - ξ
//! z           = improvement / σ                       if σ > 0
//! ei          = improvement · Φ(z) + σ · φ(z)         if σ > 0
//!             = max(improvement, 0)                   if σ == 0
//! ```
//!
//! `Φ` is the standard normal CDF and `φ` the PDF. Both are evaluated inline
//! via the Abramowitz & Stegun 7.1.26 rational approximation to `erf` (no
//! extra dependency).
//!
//! ## Randomness
//!
//! BO does not need cryptographic randomness; this module ships a small
//! xorshift64 PRNG inline. Seed via [`BoConfig::seed`] for reproducibility.

use nalgebra::DVector;

/// Configuration for [`minimize`].
#[derive(Debug, Clone)]
pub struct BoConfig {
    /// Number of Latin-hypercube initial samples drawn before BO starts.
    pub n_initial: usize,
    /// Number of BO iterations after the initial design.
    pub n_iters: usize,
    /// Number of uniform random candidates scored by Expected Improvement per
    /// iteration.
    pub n_candidates: usize,
    /// Expected Improvement exploration parameter. Larger values bias the
    /// acquisition towards higher-variance candidates.
    pub xi: f64,
    /// RNG seed for the initial design and per-iter candidate sampling.
    pub seed: u64,
}

impl Default for BoConfig {
    fn default() -> Self {
        Self {
            n_initial: 5,
            n_iters: 20,
            n_candidates: 1024,
            xi: 0.01,
            seed: 0xC0FFEE,
        }
    }
}

/// Result of a [`minimize`] run.
#[derive(Debug, Clone)]
pub struct BoResult {
    /// Best parameter vector seen.
    pub x_best: DVector<f64>,
    /// Objective value at `x_best`.
    pub y_best: f64,
    /// Full evaluation history in chronological order: initial design first,
    /// then each BO iteration's selected candidate.
    pub history: Vec<(DVector<f64>, f64)>,
}

/// Expected Improvement (minimization formulation).
///
/// - `mean` — predictive mean at the candidate point.
/// - `std` — predictive standard deviation at the candidate point. Must be
///   non-negative; values `≤ 0` are treated as the deterministic limit.
/// - `f_best` — current best (lowest) observed objective value.
/// - `xi` — exploration parameter; larger values favor uncertain candidates.
///
/// Returns a non-negative scalar; larger is better.
pub fn ei(mean: f64, std: f64, f_best: f64, xi: f64) -> f64 {
    let improvement = f_best - mean - xi;
    if std <= 0.0 {
        return improvement.max(0.0);
    }
    let z = improvement / std;
    improvement * std_normal_cdf(z) + std * std_normal_pdf(z)
}

/// Standard normal probability density function.
fn std_normal_pdf(z: f64) -> f64 {
    const INV_SQRT_2PI: f64 = 0.398_942_280_401_432_7; // 1 / sqrt(2π)
    INV_SQRT_2PI * (-0.5 * z * z).exp()
}

/// Standard normal cumulative distribution function:
/// `Φ(z) = 0.5 · (1 + erf(z / √2))`.
fn std_normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf_as(z / std::f64::consts::SQRT_2))
}

/// Abramowitz & Stegun 7.1.26 rational approximation to `erf(x)`.
///
/// Maximum absolute error ≈ 1.5·10⁻⁷, more than enough for an acquisition
/// function. Inlined so the surrogate crate keeps a clean dep tree.
fn erf_as(x: f64) -> f64 {
    // A&S 7.1.26 coefficients.
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + P * ax);
    let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * (-ax * ax).exp();
    sign * y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erf_matches_known_values() {
        // Reference values from a high-precision erf (e.g. NumPy).
        for (x, expected) in [
            (0.0, 0.0),
            (0.5, 0.520_499_877_813_046_5),
            (1.0, 0.842_700_792_949_714_9),
            (2.0, 0.995_322_265_018_952_7),
            (-1.0, -0.842_700_792_949_714_9),
        ] {
            let got = erf_as(x);
            assert!(
                (got - expected).abs() < 2e-7,
                "erf_as({x}) = {got}, expected {expected}"
            );
        }
    }

    #[test]
    fn std_normal_cdf_at_zero_is_half() {
        // A&S 7.1.26 has absolute error ≲ 1.5e-7; CDF inherits that bound.
        assert!((std_normal_cdf(0.0) - 0.5).abs() < 1e-7);
    }

    #[test]
    fn std_normal_pdf_at_zero_is_inv_sqrt_2pi() {
        let expected = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
        assert!((std_normal_pdf(0.0) - expected).abs() < 1e-12);
    }

    #[test]
    fn ei_zero_std_returns_clipped_improvement() {
        // Deterministic predictor: EI collapses to max(improvement, 0).
        assert_eq!(ei(0.5, 0.0, 1.0, 0.0), 0.5);
        assert_eq!(ei(1.5, 0.0, 1.0, 0.0), 0.0);
    }

    #[test]
    fn ei_increases_with_std_for_uncertain_candidate() {
        // For a candidate predicted slightly above the best, raising the
        // stddev increases EI (more upside through the tail).
        let f_best = 0.0;
        let mean = 0.2;
        let lo = ei(mean, 0.1, f_best, 0.0);
        let hi = ei(mean, 0.5, f_best, 0.0);
        assert!(hi > lo, "EI should grow with std: {lo} -> {hi}");
    }
}
