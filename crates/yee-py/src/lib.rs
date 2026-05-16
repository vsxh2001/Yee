//! Python bindings for Yee electromagnetic simulation.
//!
//! Wraps the Rust core via PyO3 0.28. The pymodule itself is named `_yee`
//! and is wrapped by a pure-Python package `yee` under
//! `crates/yee-py/python/yee/` for future `.pyi` stubs and convenience
//! helpers.

use pyo3::prelude::*;

mod errors;
mod freq;
mod solver;
mod sparams;
mod touchstone;
mod trimesh;

#[pymodule]
fn _yee(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<trimesh::PyTriMesh>()?;
    m.add_class::<freq::PyFreqRange>()?;
    m.add_class::<sparams::PySParameters>()?;
    m.add_class::<solver::PyPlanarMoM>()?;
    let ts_mod = PyModule::new(py, "touchstone")?;
    touchstone::register(&ts_mod)?;
    m.add_submodule(&ts_mod)?;
    Ok(())
}
