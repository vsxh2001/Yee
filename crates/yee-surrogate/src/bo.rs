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

use nalgebra::{DMatrix, DVector};

use crate::{GaussianProcess, MlFitConfig};

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

/// Minimize a black-box scalar objective over a hyper-rectangle via Bayesian
/// optimization with Expected Improvement. See the module docs for the
/// algorithm.
///
/// `bounds` is a `(lo, hi)` pair per dimension. The objective is called
/// `n_initial + n_iters` times in total. The returned [`BoResult::history`]
/// records every evaluation in chronological order; `x_best` / `y_best`
/// reflect the lowest `y` seen.
///
/// Panics if `bounds` is empty, any `(lo, hi)` has `hi <= lo`, or if any GP
/// refit fails. The small-n regime normally does not trip the GP, but a
/// pathological objective could; in that case the caller should widen the
/// initial `sigma_n` in [`MlFitConfig`] and retry.
pub fn minimize<F>(objective: F, bounds: Vec<(f64, f64)>, cfg: BoConfig) -> BoResult
where
    F: Fn(&DVector<f64>) -> f64,
{
    assert!(!bounds.is_empty(), "minimize: bounds must be non-empty");
    for (i, (lo, hi)) in bounds.iter().enumerate() {
        assert!(hi > lo, "minimize: bounds[{i}] has hi ({hi}) <= lo ({lo})");
    }
    assert!(
        cfg.n_initial >= 2,
        "minimize: n_initial must be ≥ 2 (got {}); GP needs ≥ 2 distinct training points",
        cfg.n_initial
    );

    let d = bounds.len();
    let mut rng = Xorshift64::new(cfg.seed);

    // Initial Latin-hypercube design.
    let initial = latin_hypercube(cfg.n_initial, &bounds, &mut rng);
    let mut history: Vec<(DVector<f64>, f64)> = Vec::with_capacity(cfg.n_initial + cfg.n_iters);
    for x in initial {
        let y = objective(&x);
        history.push((x, y));
    }

    for _ in 0..cfg.n_iters {
        // Refit the GP on the current history.
        let n = history.len();
        let mut x_train = DMatrix::<f64>::zeros(n, d);
        let mut y_train = DVector::<f64>::zeros(n);
        for (i, (x, y)) in history.iter().enumerate() {
            for j in 0..d {
                x_train[(i, j)] = x[j];
            }
            y_train[i] = *y;
        }

        let gp = GaussianProcess::fit_ml(x_train, y_train, MlFitConfig::default())
            .expect("BO: GP refit failed; consider widening initial sigma_n");

        let f_best = history
            .iter()
            .map(|(_, y)| *y)
            .fold(f64::INFINITY, f64::min);

        // Score n_candidates uniform random points, pick the EI maximizer.
        let mut best_x = uniform_in_bounds(&bounds, &mut rng);
        let (best_mean, best_var) = gp.predict(&best_x);
        let mut best_ei = ei(best_mean, best_var.sqrt(), f_best, cfg.xi);
        for _ in 1..cfg.n_candidates {
            let cand = uniform_in_bounds(&bounds, &mut rng);
            let (mean, var) = gp.predict(&cand);
            let v = ei(mean, var.sqrt(), f_best, cfg.xi);
            if v > best_ei {
                best_ei = v;
                best_x = cand;
            }
        }

        let y_new = objective(&best_x);
        history.push((best_x, y_new));
    }

    // Find best.
    let mut idx_best = 0;
    let mut y_best = history[0].1;
    for (i, (_, y)) in history.iter().enumerate().skip(1) {
        if *y < y_best {
            y_best = *y;
            idx_best = i;
        }
    }
    let x_best = history[idx_best].0.clone();
    BoResult {
        x_best,
        y_best,
        history,
    }
}

/// Latin-hypercube sample of `n` points in `bounds`.
///
/// For each dimension, splits `[0, 1]` into `n` equal strata, draws one
/// uniform point per stratum, then independently permutes the strata across
/// dimensions and scales to the requested bounds.
pub(crate) fn latin_hypercube(
    n: usize,
    bounds: &[(f64, f64)],
    rng: &mut Xorshift64,
) -> Vec<DVector<f64>> {
    let d = bounds.len();
    // Unit-cube LHS values: lhs_unit[i][j] is the i-th sample's j-th coord in [0, 1].
    let mut lhs_unit = vec![vec![0.0_f64; d]; n];
    #[allow(clippy::needless_range_loop)]
    for j in 0..d {
        // Generate one stratified value per stratum.
        let mut col: Vec<f64> = (0..n)
            .map(|k| ((k as f64) + rng.next_f64()) / (n as f64))
            .collect();
        // Fisher-Yates shuffle across samples for this dimension.
        for i in (1..n).rev() {
            let r = rng.next_u64() as usize % (i + 1);
            col.swap(i, r);
        }
        for (i, v) in col.into_iter().enumerate() {
            lhs_unit[i][j] = v;
        }
    }
    // Scale to bounds.
    lhs_unit
        .into_iter()
        .map(|row| {
            let mut x = DVector::<f64>::zeros(d);
            for (j, &(lo, hi)) in bounds.iter().enumerate() {
                x[j] = lo + (hi - lo) * row[j];
            }
            x
        })
        .collect()
}

/// One uniform sample in the hyper-rectangle `bounds`.
fn uniform_in_bounds(bounds: &[(f64, f64)], rng: &mut Xorshift64) -> DVector<f64> {
    let mut x = DVector::<f64>::zeros(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        x[j] = lo + (hi - lo) * rng.next_f64();
    }
    x
}

/// Tiny xorshift64 PRNG. Seeded reproducibly from [`BoConfig::seed`].
///
/// This is deliberately not exposed in the public API; BO callers control
/// randomness via [`BoConfig::seed`] only.
#[derive(Debug, Clone)]
pub(crate) struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    /// Construct a new PRNG. A zero seed is replaced by a non-zero constant so
    /// the recurrence does not collapse to `state = 0`.
    pub(crate) fn new(seed: u64) -> Self {
        let s = if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        };
        Self { state: s }
    }

    /// Next raw `u64` from the stream.
    pub(crate) fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Next `f64` in `[0, 1)`.
    pub(crate) fn next_f64(&mut self) -> f64 {
        // Take the top 53 bits, scale by 2^-53.
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1_u64 << 53) as f64))
    }
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

    #[test]
    fn xorshift_is_deterministic_from_seed() {
        let mut a = Xorshift64::new(42);
        let mut b = Xorshift64::new(42);
        for _ in 0..16 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn latin_hypercube_respects_bounds_and_stratification() {
        let mut rng = Xorshift64::new(7);
        let bounds = vec![(0.0, 6.0)];
        let n = 8;
        let pts = latin_hypercube(n, &bounds, &mut rng);
        assert_eq!(pts.len(), n);
        // Every stratum [k/n, (k+1)/n) (scaled to bounds) gets exactly one point.
        let mut hit = vec![false; n];
        for p in &pts {
            let u = (p[0] - 0.0) / 6.0;
            let k = (u * n as f64).floor() as usize;
            let k = k.min(n - 1);
            assert!(!hit[k], "stratum {k} hit twice");
            hit[k] = true;
            assert!(p[0] >= 0.0 && p[0] <= 6.0);
        }
        assert!(hit.iter().all(|&b| b), "all strata must be covered");
    }
}
