//! Notebook helper functions: numpy-friendly S-parameter conversions.
//!
//! These are intended for direct use from Jupyter / scripts via the
//! re-exported `yee.s11_db`, `yee.s11_phase`, `yee.smith_xy` symbols.
//! They operate on flat 1-D complex arrays (typically a single S-parameter
//! sweep over frequency) and return numpy arrays.

use num_complex::Complex64;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::prelude::*;

const MIN_DB: f64 = -200.0;

/// 20 · log10(|S|) clamped to `MIN_DB` at exact zero magnitude.
///
/// # Errors
///
/// Currently infallible; signature returns `PyResult` for forward-compat
/// (future shape / dtype validation may surface errors here).
#[pyfunction]
pub fn s11_db<'py>(
    py: Python<'py>,
    s: PyReadonlyArray1<'_, Complex64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let arr = s.as_array();
    let out: Vec<f64> = arr
        .iter()
        .map(|c| {
            let m = c.norm();
            if m <= 0.0 { MIN_DB } else { 20.0 * m.log10() }
        })
        .collect();
    Ok(out.into_pyarray(py))
}

/// Phase angle in degrees in (-180, 180].
///
/// # Errors
///
/// Currently infallible; signature returns `PyResult` for forward-compat.
#[pyfunction]
pub fn s11_phase<'py>(
    py: Python<'py>,
    s: PyReadonlyArray1<'_, Complex64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let arr = s.as_array();
    let out: Vec<f64> = arr.iter().map(|c| c.arg().to_degrees()).collect();
    Ok(out.into_pyarray(py))
}

/// Smith-chart (x, y) Cartesian coordinates of S₁₁ on the complex unit disk.
///
/// Returns a stacked 2-D array of shape `(N, 2)` where column 0 is the
/// real part and column 1 is the imaginary part.
///
/// # Errors
///
/// Currently infallible; signature returns `PyResult` for forward-compat.
#[pyfunction]
pub fn smith_xy<'py>(
    py: Python<'py>,
    s: PyReadonlyArray1<'_, Complex64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let arr = s.as_array();
    let n = arr.len();
    let mut buf = vec![0.0_f64; n * 2];
    for (i, c) in arr.iter().enumerate() {
        buf[i * 2] = c.re;
        buf[i * 2 + 1] = c.im;
    }
    Ok(numpy::ndarray::Array2::from_shape_vec((n, 2), buf)
        .expect("buffer length matches shape by construction")
        .into_pyarray(py))
}
