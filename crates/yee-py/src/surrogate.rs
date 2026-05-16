//! `#[pyclass] PyGaussianProcess` — Python wrapper for
//! [`yee_surrogate::GaussianProcess`].
//!
//! Exposes the scalar-output, RBF-kernel Gaussian-process regressor with
//! numpy-friendly inputs. Training data is supplied as a `(n, d)` 2-D array
//! `x` of inputs and a `(n,)` 1-D array `y` of targets. Hyperparameters
//! `(length_scale, sigma_f, sigma_n)` are passed either explicitly to
//! [`PyGaussianProcess::fit`] or initialized for the ML optimizer
//! [`PyGaussianProcess::fit_ml`].

use nalgebra::{DMatrix, DVector};
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use yee_surrogate::Error as SurrogateError;
use yee_surrogate::gp::{GaussianProcess as RustGp, MlFitConfig};

/// Map `yee_surrogate::Error` into the appropriate Python exception.
///
/// The surrogate crate is dataset/numerical-fit oriented: every variant
/// expresses caller-supplied invalid input, so they all surface as
/// `ValueError`. This keeps the bindings consistent with how `nalgebra`
/// shape mismatches are surfaced elsewhere in `yee-py`.
fn surrogate_to_py(err: SurrogateError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// Convert a row-major numpy `(n, d)` array into an `nalgebra::DMatrix<f64>`.
///
/// `nalgebra::DMatrix` is column-major, so we cannot use `from_row_slice`
/// against the numpy buffer's row-major flat view — that would transpose
/// silently. Instead build the matrix one element at a time.
fn numpy_to_dmatrix(x: &PyReadonlyArray2<'_, f64>) -> DMatrix<f64> {
    let view = x.as_array();
    let (n, d) = (view.shape()[0], view.shape()[1]);
    let mut m = DMatrix::<f64>::zeros(n, d);
    for i in 0..n {
        for j in 0..d {
            m[(i, j)] = view[[i, j]];
        }
    }
    m
}

/// Convert a numpy `(n,)` array into an `nalgebra::DVector<f64>`.
fn numpy_to_dvector(y: &PyReadonlyArray1<'_, f64>) -> DVector<f64> {
    let view = y.as_array();
    DVector::from_iterator(view.len(), view.iter().copied())
}

/// Python class wrapping a fitted [`yee_surrogate::GaussianProcess`].
///
/// Construct only via the static methods [`PyGaussianProcess::fit`] or
/// [`PyGaussianProcess::fit_ml`]; the GP has no meaningful "untrained"
/// state, matching the underlying Rust type.
#[pyclass(name = "GaussianProcess", module = "yee._yee")]
pub struct PyGaussianProcess {
    pub(crate) inner: RustGp,
}

#[pymethods]
impl PyGaussianProcess {
    /// Fit a GP with user-supplied hyperparameters.
    ///
    /// Args:
    ///     x: (n, d) ndarray of training inputs.
    ///     y: (n,) ndarray of training targets.
    ///     length_scale: RBF kernel length scale (must be > 0).
    ///     sigma_f: signal stddev (must be > 0).
    ///     sigma_n: observation noise stddev (must be >= 0).
    ///
    /// Raises:
    ///     ValueError: if shapes are inconsistent, hyperparameters are
    ///         non-positive, or the Gram matrix is not positive-definite.
    #[staticmethod]
    fn fit<'py>(
        x: PyReadonlyArray2<'py, f64>,
        y: PyReadonlyArray1<'py, f64>,
        length_scale: f64,
        sigma_f: f64,
        sigma_n: f64,
    ) -> PyResult<Self> {
        let xm = numpy_to_dmatrix(&x);
        let yv = numpy_to_dvector(&y);
        let inner = RustGp::fit(xm, yv, length_scale, sigma_f, sigma_n).map_err(surrogate_to_py)?;
        Ok(Self { inner })
    }

    /// Maximum-likelihood fit: optimize `(length_scale, sigma_f, sigma_n)`
    /// by gradient ascent on the log marginal likelihood.
    ///
    /// All `initial_*` arguments default to the values in
    /// [`yee_surrogate::gp::MlFitConfig::default`] when omitted (`None`).
    ///
    /// Args:
    ///     x: (n, d) ndarray of training inputs.
    ///     y: (n,) ndarray of training targets.
    ///     initial_length_scale: starting length scale (must be > 0).
    ///     initial_sigma_f: starting signal stddev (must be > 0).
    ///     initial_sigma_n: starting noise stddev (must be > 0).
    ///     max_iters: maximum gradient-ascent iterations.
    ///
    /// Raises:
    ///     ValueError: as for `fit`, plus if any K-build along the
    ///         optimization trajectory is non-PSD.
    #[staticmethod]
    #[pyo3(signature = (
        x,
        y,
        initial_length_scale = None,
        initial_sigma_f = None,
        initial_sigma_n = None,
        max_iters = None,
    ))]
    fn fit_ml<'py>(
        x: PyReadonlyArray2<'py, f64>,
        y: PyReadonlyArray1<'py, f64>,
        initial_length_scale: Option<f64>,
        initial_sigma_f: Option<f64>,
        initial_sigma_n: Option<f64>,
        max_iters: Option<usize>,
    ) -> PyResult<Self> {
        let xm = numpy_to_dmatrix(&x);
        let yv = numpy_to_dvector(&y);
        let mut cfg = MlFitConfig::default();
        if let Some(v) = initial_length_scale {
            cfg.initial_length_scale = v;
        }
        if let Some(v) = initial_sigma_f {
            cfg.initial_sigma_f = v;
        }
        if let Some(v) = initial_sigma_n {
            cfg.initial_sigma_n = v;
        }
        if let Some(v) = max_iters {
            cfg.max_iters = v;
        }
        let inner = RustGp::fit_ml(xm, yv, cfg).map_err(surrogate_to_py)?;
        Ok(Self { inner })
    }

    /// Posterior mean at a single query point `x_star` (shape `(d,)`).
    fn predict_mean<'py>(&self, x_star: PyReadonlyArray1<'py, f64>) -> f64 {
        let xv = numpy_to_dvector(&x_star);
        self.inner.predict_mean(&xv)
    }

    /// Posterior `(mean, variance)` at a single query point `x_star`
    /// (shape `(d,)`).
    fn predict<'py>(&self, x_star: PyReadonlyArray1<'py, f64>) -> (f64, f64) {
        let xv = numpy_to_dvector(&x_star);
        self.inner.predict(&xv)
    }

    /// Log marginal likelihood of the training data under the fitted GP.
    fn log_marginal_likelihood(&self) -> f64 {
        self.inner.log_marginal_likelihood()
    }

    /// RBF length scale this GP was fit with.
    #[getter]
    fn length_scale(&self) -> f64 {
        self.inner.length_scale()
    }

    /// Signal stddev this GP was fit with.
    #[getter]
    fn sigma_f(&self) -> f64 {
        self.inner.sigma_f()
    }

    /// Observation-noise stddev this GP was fit with.
    #[getter]
    fn sigma_n(&self) -> f64 {
        self.inner.sigma_n()
    }
}
