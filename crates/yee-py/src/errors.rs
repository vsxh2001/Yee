//! Error mapping helpers translating Rust-side `yee_*::Error` variants into
//! the appropriate `PyErr` subclass (`ValueError`, `IOError`, `RuntimeError`).
//!
//! See spec §4 Error mapping table for the canonical mapping.

use pyo3::PyErr;
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};

/// Map `yee_core::Error` into the corresponding Python exception.
pub fn yee_to_py(err: yee_core::Error) -> PyErr {
    match err {
        yee_core::Error::Invalid(msg) => PyValueError::new_err(msg),
        yee_core::Error::Numerical(msg) => PyRuntimeError::new_err(format!("numerical: {msg}")),
        yee_core::Error::Unimplemented(msg) => {
            PyRuntimeError::new_err(format!("unimplemented: {msg}"))
        }
        yee_core::Error::Io(msg) => PyIOError::new_err(msg),
    }
}

/// Map `yee_io::Error` into the corresponding Python exception.
pub fn io_to_py(err: yee_io::Error) -> PyErr {
    match err {
        yee_io::Error::Io(msg) => PyIOError::new_err(msg),
        yee_io::Error::TouchstoneParse { line, col, msg } => {
            PyValueError::new_err(format!("touchstone parse at line {line}, col {col}: {msg}"))
        }
        yee_io::Error::NotEnabled(feature) => {
            PyRuntimeError::new_err(format!("yee-io feature `{feature}` not enabled"))
        }
        yee_io::Error::InvalidFile(msg) => PyValueError::new_err(msg),
    }
}

/// Map `yee_mesh::Error` into the corresponding Python exception.
pub fn yee_mesh_to_py(err: yee_mesh::Error) -> PyErr {
    use yee_mesh::Error as E;
    match err {
        E::Invalid(msg) => PyValueError::new_err(msg),
        E::NotEnabled => PyRuntimeError::new_err("yee-mesh `gmsh` feature not enabled"),
        E::Gmsh(code) => PyRuntimeError::new_err(format!("gmsh error code {code}")),
    }
}
