//! Python bindings for Yee electromagnetic simulation.
//!
//! Wraps the Rust core via PyO3 0.28. The pymodule itself is named `_yee`
//! and is wrapped by a pure-Python package `yee` under
//! `crates/yee-py/python/yee/` for future `.pyi` stubs and convenience
//! helpers.

use pyo3::prelude::*;

mod al;
mod bo;
mod eigensolver;
mod errors;
mod fdtd;
mod fem;
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
    let ts_mod = PyModule::new(py, "yee.touchstone")?;
    touchstone::register(&ts_mod)?;
    m.add_submodule(&ts_mod)?;
    // PyO3's `add_submodule` only sets the child as an attribute of the parent;
    // it does NOT register the child in `sys.modules`. Without that registration,
    // `from yee.touchstone import read` raises `ModuleNotFoundError` even though
    // `yee.touchstone.read` (attribute access) works. Manually insert the child
    // so both import paths succeed.
    py.import("sys")?
        .getattr("modules")?
        .set_item("yee.touchstone", &ts_mod)?;
    let eig_mod = PyModule::new(py, "yee.eigensolver")?;
    eigensolver::register(&eig_mod)?;
    m.add_submodule(&eig_mod)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("yee.eigensolver", &eig_mod)?;
    let fem_mod = PyModule::new(py, "yee.fem")?;
    fem::register(&fem_mod)?;
    m.add_submodule(&fem_mod)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("yee.fem", &fem_mod)?;
    Ok(())
}
