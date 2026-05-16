//! Active learning: pick next samples by maximum predictive variance.
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
//!    c. Score each candidate by predictive variance, pick the maximizer,
//!    evaluate the objective there, append the new `(x, y)` pair.
//! 3. Refit a final GP on the full history and return it alongside the history.
//!
//! ## Acquisition
//!
//! The variance acquisition is simply `var = GaussianProcess::predict(x).1`.
//! No exploration parameter is needed — the GP's own predictive variance is
//! the uncertainty estimate, and the candidate with the highest variance is
//! the one whose evaluation most reduces posterior uncertainty (MacKay 1992).
//!
//! This is the same iteration loop as [`crate::bo::minimize`] with Expected
//! Improvement swapped for variance; it reuses the same xorshift PRNG and
//! Latin-hypercube helpers from [`crate::bo`] for reproducible randomness.

use nalgebra::{DMatrix, DVector};

use crate::bo::{Xorshift64, latin_hypercube};
use crate::{GaussianProcess, MlFitConfig};

/// Configuration for [`active_learn`].
#[derive(Debug, Clone)]
pub struct AlConfig {
    /// Number of Latin-hypercube initial samples drawn before AL starts.
    pub n_initial: usize,
    /// Number of active-learning iterations after the initial design.
    pub n_iters: usize,
    /// Number of uniform random candidates scored by predictive variance per
    /// iteration.
    pub n_candidates: usize,
    /// RNG seed for the initial design and per-iter candidate sampling.
    pub seed: u64,
}

impl Default for AlConfig {
    fn default() -> Self {
        Self {
            n_initial: 5,
            n_iters: 20,
            n_candidates: 1024,
            seed: 0xC0FFEE,
        }
    }
}

/// Result of an [`active_learn`] run.
#[derive(Debug, Clone)]
pub struct AlResult {
    /// Full evaluation history in chronological order: initial design first,
    /// then each AL iteration's selected candidate.
    pub history: Vec<(DVector<f64>, f64)>,
    /// GP refit on the full observation history. Callers can call
    /// [`GaussianProcess::predict`] / [`GaussianProcess::predict_mean`]
    /// directly on this for downstream surrogate use.
    pub final_gp: GaussianProcess,
}

/// Run active learning: starting from a Latin-hypercube initial design,
/// iteratively pick the point of maximum predictive variance and query the
/// black-box at it. See the module docs for the algorithm.
///
/// `bounds` is a `(lo, hi)` pair per dimension. The objective is called
/// `n_initial + n_iters` times in total. The returned [`AlResult::history`]
/// records every evaluation in chronological order; [`AlResult::final_gp`]
/// is a fresh fit on the full history.
///
/// Panics if `bounds` is empty, any `(lo, hi)` has `hi <= lo`,
/// `n_initial < 2`, or if any GP refit fails. The small-n regime normally
/// does not trip the GP, but a pathological objective could; in that case
/// the caller should widen the initial `sigma_n` in [`MlFitConfig`] and
/// retry.
pub fn active_learn<F>(objective: F, bounds: Vec<(f64, f64)>, cfg: AlConfig) -> AlResult
where
    F: Fn(&DVector<f64>) -> f64,
{
    assert!(!bounds.is_empty(), "active_learn: bounds must be non-empty");
    for (i, (lo, hi)) in bounds.iter().enumerate() {
        assert!(
            hi > lo,
            "active_learn: bounds[{i}] has hi ({hi}) <= lo ({lo})"
        );
    }
    assert!(
        cfg.n_initial >= 2,
        "active_learn: n_initial must be ≥ 2 (got {}); GP needs ≥ 2 distinct training points",
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
        let gp = fit_gp_on_history(&history, d);

        // Score n_candidates uniform random points, pick the max-variance candidate.
        let mut best_x = uniform_in_bounds(&bounds, &mut rng);
        let (_, best_var) = gp.predict(&best_x);
        let mut best_v = best_var;
        for _ in 1..cfg.n_candidates {
            let cand = uniform_in_bounds(&bounds, &mut rng);
            let (_, var) = gp.predict(&cand);
            if var > best_v {
                best_v = var;
                best_x = cand;
            }
        }

        let y_new = objective(&best_x);
        history.push((best_x, y_new));
    }

    // Final GP fit on the full history; the user gets a ready-to-predict
    // surrogate alongside the raw observations.
    let final_gp = fit_gp_on_history(&history, d);

    AlResult { history, final_gp }
}

/// Build the `(X, y)` design matrix from the running history and refit a
/// GP via [`GaussianProcess::fit_ml`]. Panics on fit failure (matches the
/// BO module's behavior).
fn fit_gp_on_history(history: &[(DVector<f64>, f64)], d: usize) -> GaussianProcess {
    let n = history.len();
    let mut x_train = DMatrix::<f64>::zeros(n, d);
    let mut y_train = DVector::<f64>::zeros(n);
    for (i, (x, y)) in history.iter().enumerate() {
        for j in 0..d {
            x_train[(i, j)] = x[j];
        }
        y_train[i] = *y;
    }
    GaussianProcess::fit_ml(x_train, y_train, MlFitConfig::default())
        .expect("active_learn: GP refit failed; consider widening initial sigma_n")
}

/// One uniform sample in the hyper-rectangle `bounds`. Mirrors the
/// private helper in [`crate::bo`]; kept inline so this module needs no
/// extra `pub(crate)` exposure from `bo`.
fn uniform_in_bounds(bounds: &[(f64, f64)], rng: &mut Xorshift64) -> DVector<f64> {
    let mut x = DVector::<f64>::zeros(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        x[j] = lo + (hi - lo) * rng.next_f64();
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn al_config_defaults_match_spec() {
        let cfg = AlConfig::default();
        assert_eq!(cfg.n_initial, 5);
        assert_eq!(cfg.n_iters, 20);
        assert_eq!(cfg.n_candidates, 1024);
    }

    #[test]
    fn al_history_length_matches_budget() {
        let cfg = AlConfig {
            n_initial: 4,
            n_iters: 6,
            n_candidates: 32,
            seed: 1,
        };
        let res = active_learn(
            |x: &DVector<f64>| x[0].sin(),
            vec![(0.0, std::f64::consts::TAU)],
            cfg,
        );
        assert_eq!(res.history.len(), 10);
    }
}
