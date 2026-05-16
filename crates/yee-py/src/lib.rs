//! Python bindings for Yee electromagnetic simulation.
//!
//! Wraps the Rust core via PyO3 0.28. The pymodule itself is named `_yee`
//! and is wrapped by a pure-Python package `yee` under
//! `crates/yee-py/python/yee/` for future `.pyi` stubs and convenience
//! helpers.

use pyo3::prelude::*;

mod al;
mod bo;
mod errors;
mod fdtd;
mod freq;
mod helpers;
mod nsga2;
mod solver;
mod sparams;
mod surrogate;
mod touchstone;
mod trimesh;
mod validation;

#[pymodule]
fn _yee(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<trimesh::PyTriMesh>()?;
    m.add_class::<freq::PyFreqRange>()?;
    m.add_class::<sparams::PySParameters>()?;
    m.add_class::<solver::PyPlanarMoM>()?;
    m.add_class::<surrogate::PyGaussianProcess>()?;
    m.add_class::<fdtd::PyFdtdDriverConfig>()?;
    m.add_class::<fdtd::PyFdtdDriver>()?;
    m.add_class::<fdtd::PyRadiationPattern>()?;
    m.add_class::<bo::PyBoConfig>()?;
    m.add_class::<bo::PyBoResult>()?;
    m.add_function(wrap_pyfunction!(bo::minimize, m)?)?;
    m.add_class::<nsga2::PyNsga2Config>()?;
    m.add_class::<nsga2::PyNsga2Result>()?;
    m.add_function(wrap_pyfunction!(nsga2::nsga2_minimize, m)?)?;
    m.add_class::<al::PyAlConfig>()?;
    m.add_class::<al::PyAlResult>()?;
    m.add_function(wrap_pyfunction!(al::active_learn, m)?)?;
    m.add_function(wrap_pyfunction!(helpers::s11_db, m)?)?;
    m.add_function(wrap_pyfunction!(helpers::s11_phase, m)?)?;
    m.add_function(wrap_pyfunction!(helpers::smith_xy, m)?)?;
    m.add_class::<validation::PyValidationCase>()?;
    m.add_class::<validation::PyValidationReport>()?;
    m.add_function(wrap_pyfunction!(validation::run_validation, m)?)?;
    let ts_mod = PyModule::new(py, "touchstone")?;
    touchstone::register(&ts_mod)?;
    m.add_submodule(&ts_mod)?;
    Ok(())
}
