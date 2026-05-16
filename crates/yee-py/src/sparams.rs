//! `#[pyclass] PySParameters` — Python wrapper for `yee_mom::SParameters`.

use num_complex::Complex64;
use numpy::{IntoPyArray, PyArray1, PyArray3};
use pyo3::prelude::*;
use std::path::PathBuf;
use yee_mom::SParameters as RustSParameters;

#[pyclass(name = "SParameters", module = "yee._yee")]
pub struct PySParameters {
    pub(crate) inner: RustSParameters,
}

#[pymethods]
impl PySParameters {
    #[getter]
    fn n_ports(&self) -> usize {
        self.inner.n_ports
    }

    #[getter]
    fn freq_hz<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.freq_hz.clone().into_pyarray(py)
    }

    #[getter]
    fn data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray3<Complex64>> {
        let f = self.inner.freq_hz.len();
        let n = self.inner.n_ports;
        let mut buf: Vec<Complex64> = Vec::with_capacity(f * n * n);
        for row in &self.inner.data {
            buf.extend_from_slice(row);
        }
        numpy::ndarray::Array3::from_shape_vec((f, n, n), buf)
            .expect("buffer length matches shape by construction")
            .into_pyarray(py)
    }

    fn write_touchstone(&self, path: PathBuf, z0: f64) -> PyResult<()> {
        self.inner
            .write_touchstone(&path, z0)
            .map_err(crate::errors::yee_to_py)
    }
}
