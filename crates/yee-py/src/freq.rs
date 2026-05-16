//! `#[pyclass] PyFreqRange` — Python wrapper for `yee_core::FreqRange`.

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use yee_core::FreqRange as RustFreqRange;

#[pyclass(name = "FreqRange", module = "yee._yee")]
pub struct PyFreqRange {
    pub(crate) inner: RustFreqRange,
}

#[pymethods]
impl PyFreqRange {
    #[new]
    fn new(start_hz: f64, stop_hz: f64, n_points: usize) -> PyResult<Self> {
        let inner =
            RustFreqRange::new(start_hz, stop_hz, n_points).map_err(crate::errors::yee_to_py)?;
        Ok(Self { inner })
    }

    #[getter]
    fn start_hz(&self) -> f64 {
        self.inner.start_hz
    }

    #[getter]
    fn stop_hz(&self) -> f64 {
        self.inner.stop_hz
    }

    #[getter]
    fn n_points(&self) -> usize {
        self.inner.n_points
    }

    fn iter<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let values: Vec<f64> = self.inner.iter().collect();
        values.into_pyarray(py)
    }
}
