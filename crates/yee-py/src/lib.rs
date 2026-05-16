//! Python bindings for Yee electromagnetic simulation.
//!
//! Wraps the Rust core via PyO3 0.28. The pymodule itself is named `_yee`
//! and is wrapped by a pure-Python package `yee` under
//! `crates/yee-py/python/yee/` for future `.pyi` stubs and convenience
//! helpers.

use pyo3::prelude::*;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
