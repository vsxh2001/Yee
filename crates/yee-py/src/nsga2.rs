//! NSGA-II Python bindings.
//!
//! Exposes [`yee_surrogate::nsga2::minimize`] together with its
//! [`yee_surrogate::nsga2::Nsga2Config`] and
//! [`yee_surrogate::nsga2::Nsga2Result`] types to Python. The objective is a
//! Python callable taking a 1-D `float64` numpy array of length `d` and
//! returning a length-`n_objectives` sequence of `float`.
//!
//! The Rust [`yee_surrogate::nsga2::minimize`] takes
//! `F: Fn(&DVector<f64>) -> Vec<f64>` — i.e. an infallible closure. To bridge
//! that with a Python callable that may raise, we cache the first exception
//! in a [`std::cell::RefCell`] shared with the closure, substitute
//! `+inf` for every component of the failing evaluation so the Rust loop
//! can continue, and re-raise the cached `PyErr` after `rust_minimize`
//! returns.

use std::cell::RefCell;

use nalgebra::DVector;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use yee_surrogate::nsga2::{
    Nsga2Config as RustCfg, Nsga2Result as RustResult, minimize as rust_minimize,
};

/// Python wrapper for [`yee_surrogate::nsga2::Nsga2Config`].
///
/// `mutation_probability == 0.0` is treated as "unset" and resolved to the
/// canonical default `1 / d` inside [`nsga2_minimize`] once `bounds.len()`
/// is known.
#[pyclass(name = "Nsga2Config", module = "yee._yee", from_py_object)]
#[derive(Clone)]
pub struct PyNsga2Config {
    /// Population size `N`.
    #[pyo3(get)]
    pub population_size: usize,
    /// Number of generations to run.
    #[pyo3(get)]
    pub n_generations: usize,
    /// SBX crossover distribution index `η_c`.
    #[pyo3(get)]
    pub crossover_eta: f64,
    /// Polynomial mutation distribution index `η_m`.
    #[pyo3(get)]
    pub mutation_eta: f64,
    /// Per-gene mutation probability. `0.0` resolves to `1 / d` inside
    /// [`nsga2_minimize`].
    #[pyo3(get)]
    pub mutation_probability: f64,
    /// RNG seed.
    #[pyo3(get)]
    pub seed: u64,
}

#[pymethods]
impl PyNsga2Config {
    /// Construct an NSGA-II config.
    ///
    /// Args:
    ///     population_size: number of individuals per generation (>= 2).
    ///     n_generations: number of generations to run.
    ///     crossover_eta: SBX distribution index (larger = closer to parents).
    ///     mutation_eta: polynomial mutation distribution index.
    ///     mutation_probability: per-gene mutation rate. ``None`` resolves to
    ///         ``1 / d`` once the problem dimension is known.
    ///     seed: RNG seed for reproducibility.
    #[new]
    #[pyo3(signature = (
        population_size = 100,
        n_generations = 100,
        crossover_eta = 20.0,
        mutation_eta = 20.0,
        mutation_probability = None,
        seed = 0x00C0_FFEE,
    ))]
    fn new(
        population_size: usize,
        n_generations: usize,
        crossover_eta: f64,
        mutation_eta: f64,
        mutation_probability: Option<f64>,
        seed: u64,
    ) -> Self {
        Self {
            population_size,
            n_generations,
            crossover_eta,
            mutation_eta,
            mutation_probability: mutation_probability.unwrap_or(0.0),
            seed,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Nsga2Config(population_size={}, n_generations={}, crossover_eta={}, \
mutation_eta={}, mutation_probability={}, seed={})",
            self.population_size,
            self.n_generations,
            self.crossover_eta,
            self.mutation_eta,
            self.mutation_probability,
            self.seed,
        )
    }
}

/// Python wrapper for [`yee_surrogate::nsga2::Nsga2Result`].
///
/// Attributes are exposed as numpy arrays: `population` is shape `(N, d)`,
/// `objectives` is shape `(N, m)`, `pareto_indices` is a 1-D array of
/// `int64` row indices into `population`/`objectives` for the non-dominated
/// front (front 0).
#[pyclass(name = "Nsga2Result", module = "yee._yee")]
pub struct PyNsga2Result {
    /// Row-major `(N, d)` flat buffer of the final population.
    pub(crate) population: Vec<f64>,
    /// Row-major `(N, m)` flat buffer of the final population's objectives.
    pub(crate) objectives: Vec<f64>,
    /// Indices into rows of `population` / `objectives` for the Pareto front.
    pub(crate) pareto_indices: Vec<i64>,
    pub(crate) n: usize,
    pub(crate) d: usize,
    pub(crate) m: usize,
}

#[pymethods]
impl PyNsga2Result {
    /// Final population, shape `(N, d)`.
    #[getter]
    fn population<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        Ok(
            numpy::ndarray::Array2::from_shape_vec((self.n, self.d), self.population.clone())
                .expect("population buffer length matches (N, d) by construction")
                .into_pyarray(py),
        )
    }

    /// Objective vectors for the final population, shape `(N, m)`.
    #[getter]
    fn objectives<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        Ok(
            numpy::ndarray::Array2::from_shape_vec((self.n, self.m), self.objectives.clone())
                .expect("objectives buffer length matches (N, m) by construction")
                .into_pyarray(py),
        )
    }

    /// 1-D array of `int64` indices into `population` / `objectives` for the
    /// non-dominated (Pareto) front.
    #[getter]
    fn pareto_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        self.pareto_indices.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "Nsga2Result(N={}, d={}, m={}, |F0|={})",
            self.n,
            self.d,
            self.m,
            self.pareto_indices.len(),
        )
    }
}

/// NSGA-II multi-objective minimizer.
///
/// Args:
///     objectives: callable ``f(np.ndarray of shape (d,)) -> sequence of
///         ``n_objectives`` floats``. If the callable raises, or returns the
///         wrong length, the first such error is stored and re-raised once
///         the Rust loop exits; failing evaluations are treated as
///         ``[+inf] * n_objectives`` so the inner loop can continue.
///     bounds: list of `(lo, hi)` tuples, length `d`. Each `hi` must
///         be strictly greater than its `lo`.
///     n_objectives: number `m >= 1` of objectives the callable returns.
///     config: [`PyNsga2Config`] controlling budget, operator parameters,
///         and RNG seed. Defaults to
///         [`yee_surrogate::nsga2::Nsga2Config::default`] if `None`.
///
/// Returns:
///     [`PyNsga2Result`] with `population`, `objectives`, `pareto_indices`.
#[pyfunction]
#[pyo3(signature = (objectives, bounds, n_objectives, config = None))]
pub fn nsga2_minimize<'py>(
    py: Python<'py>,
    objectives: Bound<'py, PyAny>,
    bounds: Vec<(f64, f64)>,
    n_objectives: usize,
    config: Option<PyNsga2Config>,
) -> PyResult<PyNsga2Result> {
    if bounds.is_empty() {
        return Err(PyValueError::new_err("nsga2: bounds must be non-empty"));
    }
    if n_objectives < 1 {
        return Err(PyValueError::new_err("nsga2: n_objectives must be >= 1"));
    }
    let d = bounds.len();
    let m = n_objectives;
    let mut cfg =
        config.unwrap_or_else(|| PyNsga2Config::new(100, 100, 20.0, 20.0, None, 0x00C0_FFEE));
    if cfg.mutation_probability == 0.0 {
        cfg.mutation_probability = 1.0 / d as f64;
    }
    let rust_cfg = RustCfg {
        population_size: cfg.population_size,
        n_generations: cfg.n_generations,
        crossover_eta: cfg.crossover_eta,
        mutation_eta: cfg.mutation_eta,
        mutation_probability: cfg.mutation_probability,
        seed: cfg.seed,
    };

    // Cache the first PyErr; the Rust closure must return Vec<f64>, so we
    // cannot propagate errors directly.
    let caught: RefCell<Option<PyErr>> = RefCell::new(None);

    let result: RustResult = {
        let objectives = &objectives;
        let caught_ref = &caught;
        let cb = move |x: &DVector<f64>| -> Vec<f64> {
            if caught_ref.borrow().is_some() {
                return vec![f64::INFINITY; m];
            }
            let arr = PyArray1::<f64>::from_slice(py, x.as_slice());
            match objectives.call1((arr,)) {
                Ok(obj) => match obj.extract::<Vec<f64>>() {
                    Ok(v) if v.len() == m => v,
                    Ok(v) => {
                        *caught_ref.borrow_mut() = Some(PyValueError::new_err(format!(
                            "objective returned {} values, expected {}",
                            v.len(),
                            m
                        )));
                        vec![f64::INFINITY; m]
                    }
                    Err(e) => {
                        *caught_ref.borrow_mut() = Some(e);
                        vec![f64::INFINITY; m]
                    }
                },
                Err(e) => {
                    *caught_ref.borrow_mut() = Some(e);
                    vec![f64::INFINITY; m]
                }
            }
        };
        rust_minimize(cb, bounds, n_objectives, rust_cfg)
    };

    if let Some(err) = caught.into_inner() {
        return Err(err);
    }

    let n = result.population.len();
    let mut population = Vec::with_capacity(n * d);
    let mut objectives_buf = Vec::with_capacity(n * m);
    for (x, y) in result.population.iter().zip(result.objectives.iter()) {
        for j in 0..d {
            population.push(x[j]);
        }
        for k in 0..m {
            objectives_buf.push(y[k]);
        }
    }
    let pareto_indices: Vec<i64> = result
        .pareto_front_indices
        .iter()
        .map(|&i| i as i64)
        .collect();

    Ok(PyNsga2Result {
        population,
        objectives: objectives_buf,
        pareto_indices,
        n,
        d,
        m,
    })
}
