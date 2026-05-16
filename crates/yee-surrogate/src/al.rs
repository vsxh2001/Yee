//! Active learning: pick next samples by maximum predictive variance.
//!
//! This module is the variance-acquisition counterpart to [`crate::bo`]:
//! instead of Expected Improvement (which favors low-objective regions),
//! the AL loop selects the next query point by maximum GP predictive
//! variance (which favors high-uncertainty regions). The same iteration
//! loop, the same xorshift PRNG, the same Latin-hypercube initial design.
//!
//! Phase 3.al.0 ships only the public-API skeleton ([`AlConfig`],
//! [`AlResult`]) and the [`active_learn`] entry point; the iteration loop
//! body lands in a follow-up commit.

use nalgebra::DVector;

use crate::GaussianProcess;

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
/// black-box at it.
///
/// Skeleton stub: panics until the iteration loop lands in the next commit.
pub fn active_learn<F>(_objective: F, _bounds: Vec<(f64, f64)>, _cfg: AlConfig) -> AlResult
where
    F: Fn(&DVector<f64>) -> f64,
{
    unimplemented!("active_learn: iteration loop lands in the next commit")
}
