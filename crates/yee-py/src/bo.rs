//! Bayesian-optimization Python bindings.
//!
//! Exposes [`yee_surrogate::bo::minimize`] together with its
//! [`yee_surrogate::bo::BoConfig`] and [`yee_surrogate::bo::BoResult`]
//! types to Python. The objective is supplied as an arbitrary Python
//! callable taking a 1-D `float64` numpy array of length `d` and
//! returning a Python `float`.
//!
//! The Rust [`yee_surrogate::bo::minimize`] takes
//! `F: Fn(&DVector<f64>) -> f64` — i.e. an infallible closure. To bridge
//! that with a Python callable that may raise, we cache the first
//! exception in a [`std::cell::RefCell`] shared with the closure,
//! substitute `f64::INFINITY` for the failing evaluation so the Rust
//! loop can continue, and re-raise the cached `PyErr` after
//! `rust_minimize` returns.

use std::cell::RefCell;

use nalgebra::DVector;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use yee_surrogate::bo::{BoConfig as RustCfg, BoResult as RustResult, minimize as rust_minimize};

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

/// Bayesian-optimization minimizer.
///
/// Args:
///     objective: callable ``f(np.ndarray of shape (d,)) -> float``.
///         If the callable raises, the first exception is stored and
///         re-raised once the Rust loop exits; the failing evaluation
///         is treated as ``+inf`` so the inner loop can continue.
///     bounds: list of `(lo, hi)` tuples, length `d`. Each `hi` must
///         be strictly greater than its `lo`.
///     config: [`PyBoConfig`] controlling the budget and RNG seed.
///         Defaults to [`yee_surrogate::bo::BoConfig::default`] if `None`.
///
/// Returns:
///     [`PyBoResult`] with `x_best`, `y_best`, `history_x`, `history_y`.
#[pyfunction]
#[pyo3(signature = (objective, bounds, config = None))]
pub fn minimize<'py>(
    py: Python<'py>,
    objective: Bound<'py, PyAny>,
    bounds: Vec<(f64, f64)>,
    config: Option<PyBoConfig>,
) -> PyResult<PyBoResult> {
    let cfg = config.unwrap_or_else(|| PyBoConfig::new(5, 20, 1024, 0.01, 0x00C0_FFEE));
    let rust_cfg = cfg.to_rust();
    let d = bounds.len();

    // Cache the first PyErr raised by the Python callable; the Rust
    // closure must return f64, so we cannot propagate errors directly.
    let caught: RefCell<Option<PyErr>> = RefCell::new(None);

    let result: RustResult = {
        let objective = &objective;
        let caught_ref = &caught;
        let cb = move |x: &DVector<f64>| -> f64 {
            // If an error was already raised, short-circuit with +inf so
            // the Rust loop exits as quickly as possible without invoking
            // Python again.
            if caught_ref.borrow().is_some() {
                return f64::INFINITY;
            }
            // The GIL is held: rust_minimize is being called from a
            // #[pyfunction]; no detach/attach is needed.
            let arr = PyArray1::<f64>::from_slice(py, x.as_slice());
            match objective.call1((arr,)) {
                Ok(obj) => match obj.extract::<f64>() {
                    Ok(v) => v,
                    Err(e) => {
                        *caught_ref.borrow_mut() = Some(e);
                        f64::INFINITY
                    }
                },
                Err(e) => {
                    *caught_ref.borrow_mut() = Some(e);
                    f64::INFINITY
                }
            }
        };
        rust_minimize(cb, bounds, rust_cfg)
    };

    if let Some(err) = caught.into_inner() {
        return Err(err);
    }

    // Flatten the (DVector, f64) history into row-major (n_evals, d) +
    // a parallel (n_evals,) y-buffer.
    let n_evals = result.history.len();
    let mut history_x = Vec::with_capacity(n_evals * d);
    let mut history_y = Vec::with_capacity(n_evals);
    for (x, y) in &result.history {
        for j in 0..d {
            history_x.push(x[j]);
        }
        history_y.push(*y);
    }
    let x_best: Vec<f64> = result.x_best.iter().copied().collect();
    Ok(PyBoResult {
        x_best,
        y_best: result.y_best,
        history_x,
        history_y,
        n_evals,
        d,
    })
}
