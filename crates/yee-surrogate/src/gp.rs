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
//! decomposition once and pre-solves `α = K^{-1} y` for fast mean prediction.
//!
//! ## Hyperparameters
//!
//! - `length_scale` — RBF length scale `ℓ`. Controls smoothness; larger
//!   values produce a smoother posterior mean.
//! - `sigma_f` — signal standard deviation. Controls the prior amplitude.
//! - `sigma_n` — observation-noise standard deviation. Acts as a Tikhonov
//!   regularizer (`sigma_n^2 I` added to the diagonal of `K`).

use nalgebra::{Cholesky, DMatrix, DVector, Dyn};

use crate::{Error, Result};

/// Squared-exponential Gaussian-process regressor over a real scalar output.
///
/// Construct via [`GaussianProcess::fit`].
#[derive(Debug, Clone)]
pub struct GaussianProcess {
    /// Training input matrix, shape (n, d).
    x_train: DMatrix<f64>,
    /// Training output vector, shape (n,).
    #[allow(dead_code)]
    y_train: DVector<f64>,
    /// RBF kernel length scale.
    length_scale: f64,
    /// Signal standard deviation.
    sigma_f: f64,
    /// Observation noise standard deviation.
    sigma_n: f64,
    /// Pre-computed `α = K^{-1} y` for fast mean prediction.
    #[allow(dead_code)]
    k_inv_y: DVector<f64>,
    /// Cached Cholesky factor of `K + sigma_n^2 I` (reserved for variance).
    #[allow(dead_code)]
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
}
