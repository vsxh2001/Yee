//! `#[pyclass] PyPlanarMoM` — Python wrapper for `yee_mom::PlanarMoM`.

use pyo3::prelude::*;
use yee_core::Solver;
use yee_mom::PlanarMoM as RustPlanarMoM;

#[pyclass(name = "PlanarMoM", module = "yee._yee")]
pub struct PyPlanarMoM {
    pub(crate) inner: RustPlanarMoM,
}

#[pymethods]
impl PyPlanarMoM {
    #[new]
    fn new() -> Self {
        Self {
            inner: RustPlanarMoM::default(),
        }
    }

    fn run(
        &self,
        mesh: &crate::trimesh::PyTriMesh,
        freq: &crate::freq::PyFreqRange,
    ) -> PyResult<crate::sparams::PySParameters> {
        let s = self
            .inner
            .run(&mesh.inner, freq.inner)
            .map_err(crate::errors::yee_to_py)?;
        Ok(crate::sparams::PySParameters { inner: s })
    }
}
