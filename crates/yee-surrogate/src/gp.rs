//! Gaussian-process regression surrogate.
//!
//! Implements scalar-valued GP regression with the squared-exponential (RBF)
//! kernel
//!
//! ```text
//! k(x, x') = sigma_f^2 · exp(-‖x - x'‖^2 / (2 · length_scale^2))
//! ```
//!
//! Training factors the Gram matrix `K + sigma_n^2 I` via a Cholesky
//! decomposition once, caches the factor for variance queries, and pre-solves
//! `α = K^{-1} y` for fast mean prediction. The cached factorization makes
//! per-query cost O(n) for the mean and O(n^2) for the variance, where `n`
//! is the number of training points.
//!
//! ## Hyperparameters
//!
//! - `length_scale` — RBF length scale `ℓ`. Controls smoothness; larger
//!   values produce a smoother posterior mean.
//! - `sigma_f` — signal standard deviation. Controls the prior amplitude.
//! - `sigma_n` — observation-noise standard deviation. Acts as a Tikhonov
//!   regularizer (`sigma_n^2 I` added to the diagonal of `K`).
//!
//! ## Hyperparameter optimization
//!
//! [`GaussianProcess::fit_ml`] maximizes the log marginal likelihood
//! `log p(y | X, θ)` over `θ = (length_scale, sigma_f, sigma_n)` via gradient
//! ascent in log-space (so positivity of each hyperparameter is preserved by
//! construction). The gradient is computed by **central differences** rather
//! than analytically: the analytic gradient requires `tr(K⁻¹ ∂K/∂θ)` which is
//! O(n³) per parameter and easy to get wrong. Central differences cost
//! 6 K-builds per iteration (two per parameter), but for the small
//! problems this surrogate is intended for (n ≲ 50) the simplicity of the
//! numerical path is the better tradeoff. Each K-build is itself O(n³), so
//! callers with larger `n` should expect optimization to dominate training
//! time.

use nalgebra::{Cholesky, DMatrix, DVector, Dyn};
use num_complex::Complex64;

use crate::{Dataset, Error, Result, Surrogate};

/// Squared-exponential Gaussian-process regressor over a real scalar output.
///
/// Construct via [`GaussianProcess::fit`].
#[derive(Debug, Clone)]
pub struct GaussianProcess {
    /// Training input matrix, shape (n, d).
    x_train: DMatrix<f64>,
    /// Training output vector, shape (n,).
    y_train: DVector<f64>,
    /// RBF kernel length scale.
    length_scale: f64,
    /// Signal standard deviation.
    sigma_f: f64,
    /// Observation noise standard deviation.
    sigma_n: f64,
    /// Pre-computed `α = K^{-1} y` for fast mean prediction.
    k_inv_y: DVector<f64>,
    /// Cached Cholesky factor of `K + sigma_n^2 I`, reused for variance.
    k_chol: Cholesky<f64, Dyn>,
}

impl GaussianProcess {
    /// Fit a GP to training data `(x, y)`.
    ///
    /// - `x` is an `(n, d)` matrix; each row is a training input.
    /// - `y` is an `(n,)` vector of corresponding scalar outputs.
    /// - `length_scale`, `sigma_f`, `sigma_n` configure the RBF kernel and
    ///   observation noise.
    ///
    /// Returns [`Error::FitFailed`] if `x` and `y` have inconsistent shapes,
    /// the dataset is empty, or the Gram matrix fails to Cholesky-factor.
    pub fn fit(
        x: DMatrix<f64>,
        y: DVector<f64>,
        length_scale: f64,
        sigma_f: f64,
        sigma_n: f64,
    ) -> Result<Self> {
        let n = x.nrows();
        if n == 0 {
            return Err(Error::EmptyDataset);
        }
        if y.len() != n {
            return Err(Error::FitFailed(format!(
                "x rows ({}) and y length ({}) disagree",
                n,
                y.len()
            )));
        }
        if length_scale <= 0.0 {
            return Err(Error::FitFailed(format!(
                "length_scale must be positive, got {length_scale}"
            )));
        }
        if sigma_f <= 0.0 {
            return Err(Error::FitFailed(format!(
                "sigma_f must be positive, got {sigma_f}"
            )));
        }
        if sigma_n < 0.0 {
            return Err(Error::FitFailed(format!(
                "sigma_n must be non-negative, got {sigma_n}"
            )));
        }

        // Build K + sigma_n^2 I.
        let mut k = DMatrix::<f64>::zeros(n, n);
        let two_l2 = 2.0 * length_scale * length_scale;
        let sf2 = sigma_f * sigma_f;
        let sn2 = sigma_n * sigma_n;
        for i in 0..n {
            for j in i..n {
                let r2 = squared_distance(&x, i, &x, j);
                let val = sf2 * (-r2 / two_l2).exp();
                k[(i, j)] = val;
                if i != j {
                    k[(j, i)] = val;
                }
            }
            k[(i, i)] += sn2;
        }

        let k_chol = Cholesky::new(k).ok_or_else(|| {
            Error::FitFailed(
                "Gram matrix is not positive-definite; \
                 try increasing sigma_n or sigma_f, or check for duplicate inputs"
                    .to_string(),
            )
        })?;

        let k_inv_y = k_chol.solve(&y);

        Ok(Self {
            x_train: x,
            y_train: y,
            length_scale,
            sigma_f,
            sigma_n,
            k_inv_y,
            k_chol,
        })
    }

    /// Number of training points.
    pub fn n_train(&self) -> usize {
        self.x_train.nrows()
    }

    /// Input dimensionality.
    pub fn input_dim(&self) -> usize {
        self.x_train.ncols()
    }

    /// Posterior mean at `x_star`. Cost is O(n · d).
    ///
    /// `x_star` must have length equal to [`Self::input_dim`].
    pub fn predict_mean(&self, x_star: &DVector<f64>) -> f64 {
        let k_star = self.kernel_vector(x_star);
        k_star.dot(&self.k_inv_y)
    }

    /// Posterior mean and variance at `x_star`. Cost is O(n^2 · d).
    ///
    /// Variance is computed as
    ///
    /// ```text
    /// k(x_star, x_star) - k_star^T · K^{-1} · k_star
    /// ```
    ///
    /// via the cached Cholesky factor. A small floor of `0.0` is applied to
    /// guard against tiny negative values that can arise from finite-precision
    /// cancellation when the query coincides with a training point.
    pub fn predict(&self, x_star: &DVector<f64>) -> (f64, f64) {
        let k_star = self.kernel_vector(x_star);
        let mean = k_star.dot(&self.k_inv_y);
        let v = self.k_chol.solve(&k_star);
        let k_ss = self.sigma_f * self.sigma_f; // k(x_star, x_star) for RBF
        let var = (k_ss - k_star.dot(&v)).max(0.0);
        (mean, var)
    }

    /// Build the `(n,)` kernel vector `k_star[i] = k(x_star, x_train[i, :])`.
    fn kernel_vector(&self, x_star: &DVector<f64>) -> DVector<f64> {
        let n = self.x_train.nrows();
        let d = self.x_train.ncols();
        assert_eq!(
            x_star.len(),
            d,
            "GaussianProcess::kernel_vector: x_star has length {} but training inputs are {}-dim",
            x_star.len(),
            d
        );
        let two_l2 = 2.0 * self.length_scale * self.length_scale;
        let sf2 = self.sigma_f * self.sigma_f;
        let mut v = DVector::<f64>::zeros(n);
        for i in 0..n {
            let mut r2 = 0.0;
            for j in 0..d {
                let dx = self.x_train[(i, j)] - x_star[j];
                r2 += dx * dx;
            }
            v[i] = sf2 * (-r2 / two_l2).exp();
        }
        v
    }

    /// Length scale this GP was fit with.
    pub fn length_scale(&self) -> f64 {
        self.length_scale
    }
    /// Signal stddev this GP was fit with.
    pub fn sigma_f(&self) -> f64 {
        self.sigma_f
    }
    /// Noise stddev this GP was fit with.
    pub fn sigma_n(&self) -> f64 {
        self.sigma_n
    }

    /// Log marginal likelihood of the training data under the fitted GP:
    ///
    /// ```text
    /// log p(y | X, θ) = -0.5 · yᵀ K⁻¹ y - 0.5 · log|K| - (n/2) · log(2π)
    /// ```
    ///
    /// Computed from cached state:
    ///
    /// - The **data-fit** term `-0.5 · yᵀ α` reuses the cached `α = K⁻¹ y`.
    /// - The **complexity** term `-0.5 · log|K| = -sum(log diag(L))` reuses
    ///   the cached lower-triangular Cholesky factor `L` of `K + σ_n² I`,
    ///   exploiting `log|K| = 2 · sum(log diag(L))`.
    /// - The **normalizer** is `-(n/2) · log(2π)`.
    ///
    /// Higher is better; useful for model selection between fitted GPs and
    /// as the objective for hyperparameter optimization via
    /// [`GaussianProcess::fit_ml`].
    pub fn log_marginal_likelihood(&self) -> f64 {
        let n = self.x_train.nrows();
        // -0.5 * y^T alpha
        let data_fit = -0.5 * self.y_train.dot(&self.k_inv_y);
        // -0.5 * log|K| = -sum(log diag(L))
        let l = self.k_chol.l();
        let mut log_diag_sum = 0.0;
        for i in 0..n {
            log_diag_sum += l[(i, i)].ln();
        }
        let complexity = -log_diag_sum;
        // -(n/2) log(2π)
        let norm = -0.5 * (n as f64) * (2.0 * std::f64::consts::PI).ln();
        data_fit + complexity + norm
    }
}

/// Configuration for marginal-likelihood hyperparameter optimization.
///
/// Passed to [`GaussianProcess::fit_ml`]. Defaults are tuned for the
/// low-dimensional, small-n surrogate problems this crate targets; expect to
/// tune `gradient_step` and `max_iters` for larger problems.
#[derive(Debug, Clone)]
pub struct MlFitConfig {
    /// Starting length scale (linear-space, not log-space).
    pub initial_length_scale: f64,
    /// Starting signal stddev.
    pub initial_sigma_f: f64,
    /// Starting noise stddev.
    pub initial_sigma_n: f64,
    /// Maximum gradient-ascent iterations.
    pub max_iters: usize,
    /// Step magnitude in log-space. Used as a learning rate when the
    /// gradient is small (‖grad‖ < 1) and as a maximum step magnitude
    /// otherwise (the gradient is normalized to unit length, then scaled by
    /// this value). The split protects against the very large gradients the
    /// log marginal likelihood produces when started far from the optimum.
    pub gradient_step: f64,
    /// Convergence threshold on the L2 norm of the log-space gradient.
    pub tol: f64,
}

impl Default for MlFitConfig {
    fn default() -> Self {
        Self {
            initial_length_scale: 1.0,
            initial_sigma_f: 1.0,
            initial_sigma_n: 1e-3,
            max_iters: 200,
            gradient_step: 0.05,
            tol: 1e-4,
        }
    }
}

impl GaussianProcess {
    /// Fit a GP by maximizing the log marginal likelihood over
    /// `(length_scale, sigma_f, sigma_n)`.
    ///
    /// Optimizes in log-space (`θ = (log ℓ, log σ_f, log σ_n)`) so the three
    /// hyperparameters stay strictly positive without an explicit constraint.
    /// The gradient is computed by central differences (step `1e-3` in
    /// log-space) on [`Self::log_marginal_likelihood`]; the optimizer is
    /// plain gradient ascent. See the module-level docs for the rationale.
    ///
    /// Iteration halts when the L2 norm of the log-space gradient drops
    /// below `cfg.tol` or `cfg.max_iters` is exhausted, whichever comes
    /// first. The returned [`GaussianProcess`] is a fresh refit with the
    /// optimized hyperparameters, so its cached `α` and Cholesky factor are
    /// consistent.
    ///
    /// Returns [`Error::FitFailed`] if any K-build along the optimization
    /// trajectory is non-PSD (i.e. Cholesky-factor fails). In that case the
    /// last successful fit is *not* returned; the caller should widen the
    /// initial `sigma_n` and retry.
    pub fn fit_ml(x: DMatrix<f64>, y: DVector<f64>, cfg: MlFitConfig) -> Result<Self> {
        if cfg.initial_length_scale <= 0.0
            || cfg.initial_sigma_f <= 0.0
            || cfg.initial_sigma_n <= 0.0
        {
            return Err(Error::FitFailed(format!(
                "fit_ml: initial hyperparameters must be strictly positive, got \
                 (length_scale={}, sigma_f={}, sigma_n={})",
                cfg.initial_length_scale, cfg.initial_sigma_f, cfg.initial_sigma_n
            )));
        }

        // Log-space parameter vector: (log ℓ, log σ_f, log σ_n).
        let mut theta = [
            cfg.initial_length_scale.ln(),
            cfg.initial_sigma_f.ln(),
            cfg.initial_sigma_n.ln(),
        ];
        let h = 1e-3_f64; // central-difference step in log-space

        // Evaluate the log marginal likelihood at a candidate θ. Wrapped so
        // we can call it repeatedly inside the central-difference loop
        // without re-deriving the (ℓ, σ_f, σ_n) decoding each time.
        let lml_at = |theta: [f64; 3]| -> Result<f64> {
            let l = theta[0].exp();
            let sf = theta[1].exp();
            let sn = theta[2].exp();
            let gp = GaussianProcess::fit(x.clone(), y.clone(), l, sf, sn)?;
            Ok(gp.log_marginal_likelihood())
        };

        // Step-scaling regime: vanilla gradient ascent overshoots violently
        // when starting far from the optimum because the log marginal
        // likelihood has very large gradients in directions where the prior
        // is mis-scaled (e.g. tiny sigma_f). We therefore treat
        // `gradient_step` as the *maximum* log-space step per iteration: if
        // ‖grad‖ exceeds 1 we scale the update by 1/‖grad‖, which preserves
        // the gradient direction but caps the move at `gradient_step` in
        // log-space. When ‖grad‖ < 1 the step is the raw gradient, so
        // convergence near the optimum still feels the curvature.
        for _ in 0..cfg.max_iters {
            // Central-difference gradient: 2 K-builds per parameter = 6 total.
            let mut grad = [0.0_f64; 3];
            for k in 0..3 {
                let mut tp = theta;
                let mut tm = theta;
                tp[k] += h;
                tm[k] -= h;
                let lp = lml_at(tp)?;
                let lm = lml_at(tm)?;
                grad[k] = (lp - lm) / (2.0 * h);
            }

            let gnorm = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
            if gnorm < cfg.tol {
                break;
            }

            let scale = if gnorm > 1.0 { 1.0 / gnorm } else { 1.0 };
            for k in 0..3 {
                theta[k] += cfg.gradient_step * scale * grad[k];
            }
        }

        let length_scale = theta[0].exp();
        let sigma_f = theta[1].exp();
        let sigma_n = theta[2].exp();
        Self::fit(x, y, length_scale, sigma_f, sigma_n)
    }
}

/// Squared Euclidean distance between row `i` of `a` and row `j` of `b`.
fn squared_distance(a: &DMatrix<f64>, i: usize, b: &DMatrix<f64>, j: usize) -> f64 {
    debug_assert_eq!(a.ncols(), b.ncols());
    let mut r2 = 0.0;
    for c in 0..a.ncols() {
        let d = a[(i, c)] - b[(j, c)];
        r2 += d * d;
    }
    r2
}

/// `Surrogate` impl that adapts the scalar GP to the existing trait shape.
///
/// Treats the **real part** of `sample.output[0]` as the regression target
/// across all training samples. Output dimensionality > 1 is rejected with
/// [`Error::InconsistentOutput`] so callers can spot when they meant to use a
/// multi-output backend. Hyperparameters are supplied at construction time;
/// callers who need multi-output GPs should construct a [`GaussianProcess`]
/// directly via [`GaussianProcess::fit`].
///
/// Stored as a wrapper rather than implemented directly on `GaussianProcess`
/// because the GP has no meaningful "untrained" state.
#[derive(Debug, Clone)]
pub struct GpSurrogate {
    inner: Option<GaussianProcess>,
    length_scale: f64,
    sigma_f: f64,
    sigma_n: f64,
}

impl Default for GpSurrogate {
    fn default() -> Self {
        Self {
            inner: None,
            length_scale: 1.0,
            sigma_f: 1.0,
            sigma_n: 1e-4,
        }
    }
}

impl GpSurrogate {
    /// Construct an untrained `GpSurrogate` with the given RBF hyperparameters.
    pub fn with_hyperparams(length_scale: f64, sigma_f: f64, sigma_n: f64) -> Self {
        Self {
            inner: None,
            length_scale,
            sigma_f,
            sigma_n,
        }
    }

    /// Access the underlying fitted [`GaussianProcess`] (if trained).
    pub fn inner(&self) -> Option<&GaussianProcess> {
        self.inner.as_ref()
    }

    /// Posterior variance at `params`. Returns `None` if not yet trained.
    pub fn predict_variance(&self, params: &[f64]) -> Option<f64> {
        let gp = self.inner.as_ref()?;
        let x = DVector::from_row_slice(params);
        Some(gp.predict(&x).1)
    }
}

impl Surrogate for GpSurrogate {
    fn train(&mut self, dataset: &Dataset) -> Result<()> {
        if dataset.is_empty() {
            return Err(Error::EmptyDataset);
        }
        let d = dataset.samples[0].params.len();
        let out_dim = dataset.samples[0].output.len();
        if out_dim != 1 {
            return Err(Error::InconsistentOutput {
                expected: 1,
                got: out_dim,
            });
        }
        for s in &dataset.samples {
            if s.params.len() != d {
                return Err(Error::DimMismatch {
                    training: d,
                    query: s.params.len(),
                });
            }
            if s.output.len() != out_dim {
                return Err(Error::InconsistentOutput {
                    expected: out_dim,
                    got: s.output.len(),
                });
            }
        }
        let n = dataset.samples.len();
        let mut x = DMatrix::<f64>::zeros(n, d);
        let mut y = DVector::<f64>::zeros(n);
        for (i, s) in dataset.samples.iter().enumerate() {
            for (j, v) in s.params.iter().enumerate() {
                x[(i, j)] = *v;
            }
            y[i] = s.output[0].re;
        }
        self.inner = Some(GaussianProcess::fit(
            x,
            y,
            self.length_scale,
            self.sigma_f,
            self.sigma_n,
        )?);
        Ok(())
    }

    fn predict(&self, params: &[f64]) -> Result<Vec<Complex64>> {
        let gp = self
            .inner
            .as_ref()
            .ok_or_else(|| Error::FitFailed("predict called before train".to_string()))?;
        if params.len() != gp.input_dim() {
            return Err(Error::DimMismatch {
                training: gp.input_dim(),
                query: params.len(),
            });
        }
        let x = DVector::from_row_slice(params);
        let mean = gp.predict_mean(&x);
        Ok(vec![Complex64::new(mean, 0.0)])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_rejects_empty() {
        let x = DMatrix::<f64>::zeros(0, 1);
        let y = DVector::<f64>::zeros(0);
        assert!(matches!(
            GaussianProcess::fit(x, y, 1.0, 1.0, 1e-3),
            Err(Error::EmptyDataset)
        ));
    }

    #[test]
    fn fit_rejects_shape_mismatch() {
        let x = DMatrix::<f64>::from_row_slice(2, 1, &[0.0, 1.0]);
        let y = DVector::<f64>::from_row_slice(&[0.0, 1.0, 2.0]);
        assert!(matches!(
            GaussianProcess::fit(x, y, 1.0, 1.0, 1e-3),
            Err(Error::FitFailed(_))
        ));
    }

    #[test]
    fn log_marginal_likelihood_matches_textbook_decomposition() {
        // Tiny n=3 problem: re-derive log p(y|X,θ) from raw K and compare.
        let x = DMatrix::<f64>::from_row_slice(3, 1, &[0.0, 1.0, 2.0]);
        let y = DVector::<f64>::from_row_slice(&[0.0, 1.0, 0.0]);
        let length_scale = 0.5;
        let sigma_f = 1.0;
        let sigma_n = 1e-2;
        let gp = GaussianProcess::fit(x.clone(), y.clone(), length_scale, sigma_f, sigma_n)
            .expect("fit");

        // Rebuild K + sigma_n^2 I directly and compute log p(y|X,θ) by the
        // textbook expression: -0.5 y^T K^-1 y - 0.5 log|K| - n/2 log(2π).
        let n = 3;
        let mut k = DMatrix::<f64>::zeros(n, n);
        let two_l2 = 2.0 * length_scale * length_scale;
        let sf2 = sigma_f * sigma_f;
        let sn2 = sigma_n * sigma_n;
        for i in 0..n {
            for j in 0..n {
                let dx = x[(i, 0)] - x[(j, 0)];
                k[(i, j)] = sf2 * (-(dx * dx) / two_l2).exp();
            }
            k[(i, i)] += sn2;
        }
        let chol = k.clone().cholesky().expect("PSD");
        let alpha = chol.solve(&y);
        let data_fit = -0.5 * y.dot(&alpha);
        // log|K| = 2 sum(log diag(L))
        let l = chol.l();
        let mut log_diag = 0.0;
        for i in 0..n {
            log_diag += l[(i, i)].ln();
        }
        let complexity = -log_diag;
        let norm = -0.5 * (n as f64) * (2.0 * std::f64::consts::PI).ln();
        let expected = data_fit + complexity + norm;

        let got = gp.log_marginal_likelihood();
        assert!(
            (got - expected).abs() < 1e-10,
            "log_marginal_likelihood = {got}, expected {expected}"
        );
    }

    #[test]
    fn variance_at_training_point_is_near_noise() {
        let x = DMatrix::<f64>::from_row_slice(3, 1, &[0.0, 1.0, 2.0]);
        let y = DVector::<f64>::from_row_slice(&[0.0, 1.0, 0.0]);
        let sigma_n = 1e-3;
        let gp = GaussianProcess::fit(x, y, 0.5, 1.0, sigma_n).unwrap();
        let x_star = DVector::<f64>::from_row_slice(&[1.0]);
        let (_, var) = gp.predict(&x_star);
        assert!(
            (var - sigma_n * sigma_n).abs() < 1e-3,
            "var at training point = {var}, expected ≈ {}",
            sigma_n * sigma_n
        );
    }
}
