//! Bayesian-optimization Python bindings.
//!
//! Exposes [`yee_surrogate::bo::minimize`] together with its
//! [`yee_surrogate::bo::BoConfig`] and [`yee_surrogate::bo::BoResult`]
//! types to Python. The objective is supplied as an arbitrary Python
//! callable taking a 1-D `float64` numpy array of length `d` and
//! returning a Python `float`.
//!
//! This commit lands only the class skeleton ([`PyBoConfig`],
//! [`PyBoResult`]) and a placeholder [`minimize`] entry point that
//! raises `NotImplementedError`. The real wrapper that drives
//! `yee_surrogate::bo::minimize` lands in the follow-up commit.

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use yee_surrogate::bo::BoConfig as RustCfg;

/// Python wrapper for [`yee_surrogate::bo::BoConfig`].
#[pyclass(name = "BoConfig", module = "yee._yee", from_py_object)]
#[derive(Clone)]
pub struct PyBoConfig {
    /// Number of Latin-hypercube initial samples drawn before BO starts.
    #[pyo3(get)]
    pub n_initial: usize,
    /// Number of BO iterations after the initial design.
    #[pyo3(get)]
    pub n_iters: usize,
    /// Number of uniform random candidates scored by Expected Improvement
    /// per iteration.
    #[pyo3(get)]
    pub n_candidates: usize,
    /// Expected Improvement exploration parameter. Larger values bias the
    /// acquisition towards higher-variance candidates.
    #[pyo3(get)]
    pub xi: f64,
    /// RNG seed for the initial design and per-iter candidate sampling.
    #[pyo3(get)]
    pub seed: u64,
}

#[pymethods]
impl PyBoConfig {
    /// Construct a config.
    ///
    /// All arguments are keyword-friendly with the same defaults as
    /// [`yee_surrogate::bo::BoConfig::default`].
    ///
    /// Args:
    ///     n_initial: number of Latin-hypercube initial samples (must
    ///         be >= 2).
    ///     n_iters: number of BO iterations after the initial design.
    ///     n_candidates: candidates scored by Expected Improvement per
    ///         iteration.
    ///     xi: Expected Improvement exploration parameter (>= 0).
    ///     seed: RNG seed for reproducibility.
    #[new]
    #[pyo3(signature = (
        n_initial = 5,
        n_iters = 20,
        n_candidates = 1024,
        xi = 0.01,
        seed = 0x00C0_FFEE,
    ))]
    fn new(n_initial: usize, n_iters: usize, n_candidates: usize, xi: f64, seed: u64) -> Self {
        Self {
            n_initial,
            n_iters,
            n_candidates,
            xi,
            seed,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "BoConfig(n_initial={}, n_iters={}, n_candidates={}, xi={}, seed={})",
            self.n_initial, self.n_iters, self.n_candidates, self.xi, self.seed,
        )
    }
}

impl PyBoConfig {
    /// Internal: convert this config into the Rust `BoConfig`.
    #[allow(dead_code)]
    fn to_rust(&self) -> RustCfg {
        RustCfg {
            n_initial: self.n_initial,
            n_iters: self.n_iters,
            n_candidates: self.n_candidates,
            xi: self.xi,
            seed: self.seed,
        }
    }
}

/// Python wrapper for [`yee_surrogate::bo::BoResult`].
///
/// Attributes are exposed as numpy arrays: `x_best` is shape `(d,)`,
/// `history_x` is shape `(n_evals, d)`, `history_y` is shape `(n_evals,)`.
/// Rows of `history_x` are ordered chronologically — first the
/// Latin-hypercube initial design, then each BO iteration's selected
/// candidate.
#[pyclass(name = "BoResult", module = "yee._yee")]
pub struct PyBoResult {
    pub(crate) x_best: Vec<f64>,
    pub(crate) y_best: f64,
    /// `n_evals × d` row-major buffer.
    pub(crate) history_x: Vec<f64>,
    pub(crate) history_y: Vec<f64>,
    pub(crate) n_evals: usize,
    pub(crate) d: usize,
}

#[pymethods]
impl PyBoResult {
    /// Best parameter vector seen, shape `(d,)`.
    #[getter]
    fn x_best<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.x_best.clone().into_pyarray(py)
    }

    /// Objective value at [`PyBoResult::x_best`].
    #[getter]
    fn y_best(&self) -> f64 {
        self.y_best
    }

    /// All evaluated points in evaluation order, shape `(n_evals, d)`.
    #[getter]
    fn history_x<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        Ok(
            numpy::ndarray::Array2::from_shape_vec((self.n_evals, self.d), self.history_x.clone())
                .expect("history_x buffer length matches (n_evals, d) by construction")
                .into_pyarray(py),
        )
    }

    /// All evaluated objective values in evaluation order, shape
    /// `(n_evals,)`.
    #[getter]
    fn history_y<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.history_y.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "BoResult(y_best={}, n_evals={}, d={})",
            self.y_best, self.n_evals, self.d,
        )
    }
}

/// Bayesian-optimization minimizer (placeholder).
///
/// The full wrapper that drives [`yee_surrogate::bo::minimize`] lands
/// in the follow-up commit. Until then this entry point exists only
/// so the module-level registration in `lib.rs` resolves; calling it
/// raises `NotImplementedError`.
#[pyfunction]
#[pyo3(signature = (objective, bounds, config = None))]
pub fn minimize<'py>(
    _py: Python<'py>,
    objective: Bound<'py, PyAny>,
    bounds: Vec<(f64, f64)>,
    config: Option<PyBoConfig>,
) -> PyResult<PyBoResult> {
    let _ = (objective, bounds, config);
    Err(PyNotImplementedError::new_err(
        "yee.bo.minimize: skeleton commit; real implementation follows",
    ))
}
