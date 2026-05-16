//! Active-learning Python bindings.
//!
//! Exposes [`yee_surrogate::al::active_learn`] together with its
//! [`yee_surrogate::al::AlConfig`] and [`yee_surrogate::al::AlResult`] types
//! to Python. The objective is a Python callable taking a 1-D `float64`
//! numpy array of length `d` and returning a Python `float`.
//!
//! The Rust [`yee_surrogate::al::active_learn`] takes
//! `F: Fn(&DVector<f64>) -> f64` — i.e. an infallible closure. To bridge that
//! with a Python callable that may raise, we cache the first exception in a
//! [`std::cell::RefCell`] shared with the closure, substitute
//! `f64::INFINITY` for the failing evaluation so the Rust loop can
//! continue, and re-raise the cached `PyErr` after the Rust loop returns.

use std::cell::RefCell;

use nalgebra::DVector;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use yee_surrogate::GaussianProcess;
use yee_surrogate::al::{AlConfig as RustCfg, AlResult as RustResult, active_learn as rust_al};

use crate::surrogate::PyGaussianProcess;

/// Python wrapper for [`yee_surrogate::al::AlConfig`].
#[pyclass(name = "AlConfig", module = "yee._yee", from_py_object)]
#[derive(Clone)]
pub struct PyAlConfig {
    /// Number of Latin-hypercube initial samples drawn before AL starts.
    #[pyo3(get)]
    pub n_initial: usize,
    /// Number of active-learning iterations after the initial design.
    #[pyo3(get)]
    pub n_iters: usize,
    /// Number of uniform random candidates scored by predictive variance per
    /// iteration.
    #[pyo3(get)]
    pub n_candidates: usize,
    /// RNG seed for the initial design and per-iter candidate sampling.
    #[pyo3(get)]
    pub seed: u64,
}

#[pymethods]
impl PyAlConfig {
    /// Construct an active-learning config.
    ///
    /// Args:
    ///     n_initial: number of Latin-hypercube initial samples (must be >= 2).
    ///     n_iters: number of active-learning iterations after the initial
    ///         design.
    ///     n_candidates: candidates scored by predictive variance per
    ///         iteration.
    ///     seed: RNG seed for reproducibility.
    #[new]
    #[pyo3(signature = (
        n_initial = 5,
        n_iters = 20,
        n_candidates = 1024,
        seed = 0x00C0_FFEE,
    ))]
    fn new(n_initial: usize, n_iters: usize, n_candidates: usize, seed: u64) -> Self {
        Self {
            n_initial,
            n_iters,
            n_candidates,
            seed,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlConfig(n_initial={}, n_iters={}, n_candidates={}, seed={})",
            self.n_initial, self.n_iters, self.n_candidates, self.seed,
        )
    }
}

impl PyAlConfig {
    /// Internal: convert this config into the Rust `AlConfig`.
    fn to_rust(&self) -> RustCfg {
        RustCfg {
            n_initial: self.n_initial,
            n_iters: self.n_iters,
            n_candidates: self.n_candidates,
            seed: self.seed,
        }
    }
}

/// Python wrapper for [`yee_surrogate::al::AlResult`].
///
/// Attributes are exposed as numpy arrays: `history_x` is shape
/// `(n_initial + n_iters, d)`, `history_y` is shape `(n_initial + n_iters,)`.
/// `final_gp()` returns a [`PyGaussianProcess`] wrapping the GP refit on
/// the full history.
#[pyclass(name = "AlResult", module = "yee._yee")]
pub struct PyAlResult {
    /// Row-major `(n_evals, d)` flat buffer of evaluated inputs in
    /// chronological order.
    pub(crate) history_x: Vec<f64>,
    /// `(n_evals,)` chronological objective values.
    pub(crate) history_y: Vec<f64>,
    pub(crate) n_evals: usize,
    pub(crate) d: usize,
    /// GP refit on the full history. Cloned out of the underlying Rust
    /// `AlResult` so `final_gp()` can hand a fresh `PyGaussianProcess` to
    /// callers without consuming `self`.
    pub(crate) final_gp: GaussianProcess,
}

#[pymethods]
impl PyAlResult {
    /// All evaluated inputs in chronological order, shape `(n_evals, d)`.
    #[getter]
    fn history_x<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        Ok(
            numpy::ndarray::Array2::from_shape_vec((self.n_evals, self.d), self.history_x.clone())
                .expect("history_x buffer length matches (n_evals, d) by construction")
                .into_pyarray(py),
        )
    }

    /// All evaluated objective values in chronological order, shape
    /// `(n_evals,)`.
    #[getter]
    fn history_y<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.history_y.clone().into_pyarray(py)
    }

    /// Return a [`PyGaussianProcess`] trained on the full AL history.
    fn final_gp(&self) -> PyGaussianProcess {
        PyGaussianProcess {
            inner: self.final_gp.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!("AlResult(n_evals={}, d={})", self.n_evals, self.d)
    }
}

/// Active-learning loop: variance-acquisition sample selection on a GP
/// surrogate.
///
/// Args:
///     objective: callable ``f(np.ndarray of shape (d,)) -> float``.
///         If the callable raises, the first exception is stored and
///         re-raised once the Rust loop exits; the failing evaluation
///         is treated as ``+inf`` so the inner loop can continue.
///     bounds: list of `(lo, hi)` tuples, length `d`. Each `hi` must
///         be strictly greater than its `lo`.
///     config: [`PyAlConfig`] controlling the budget and RNG seed.
///         Defaults to [`yee_surrogate::al::AlConfig::default`] if `None`.
///
/// Returns:
///     [`PyAlResult`] with `history_x`, `history_y`, and `final_gp()`.
#[pyfunction]
#[pyo3(signature = (objective, bounds, config = None))]
pub fn active_learn<'py>(
    py: Python<'py>,
    objective: Bound<'py, PyAny>,
    bounds: Vec<(f64, f64)>,
    config: Option<PyAlConfig>,
) -> PyResult<PyAlResult> {
    if bounds.is_empty() {
        return Err(PyValueError::new_err(
            "active_learn: bounds must be non-empty",
        ));
    }
    let d = bounds.len();
    let cfg = config.unwrap_or_else(|| PyAlConfig::new(5, 20, 1024, 0x00C0_FFEE));
    let rust_cfg = cfg.to_rust();

    let caught: RefCell<Option<PyErr>> = RefCell::new(None);

    let result: RustResult = {
        let objective = &objective;
        let caught_ref = &caught;
        let cb = move |x: &DVector<f64>| -> f64 {
            if caught_ref.borrow().is_some() {
                return f64::INFINITY;
            }
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
        rust_al(cb, bounds, rust_cfg)
    };

    if let Some(err) = caught.into_inner() {
        return Err(err);
    }

    let n_evals = result.history.len();
    let mut history_x = Vec::with_capacity(n_evals * d);
    let mut history_y = Vec::with_capacity(n_evals);
    for (x, y) in &result.history {
        for j in 0..d {
            history_x.push(x[j]);
        }
        history_y.push(*y);
    }
    Ok(PyAlResult {
        history_x,
        history_y,
        n_evals,
        d,
        final_gp: result.final_gp,
    })
}
