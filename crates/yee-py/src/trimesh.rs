//! `#[pyclass] PyTriMesh` — Python wrapper for `yee_mesh::TriMesh`.
//!
//! Accepts numpy arrays for vertices (`[N, 3]` f64), triangles (`[M, 3]` u32),
//! and tags (`[M]` u32). Getters return numpy arrays — a structural copy is
//! performed because the underlying Rust storage (`Vec<Vector3<f64>>`,
//! `Vec<[u32; 3]>`) is not directly viewable as a strided numpy buffer.

use nalgebra::Vector3;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use yee_mesh::TriMesh as RustTriMesh;

#[pyclass(name = "TriMesh", module = "yee._yee")]
pub struct PyTriMesh {
    pub(crate) inner: RustTriMesh,
}

#[pymethods]
impl PyTriMesh {
    #[new]
    fn new(
        vertices: PyReadonlyArray2<'_, f64>,
        triangles: PyReadonlyArray2<'_, u32>,
        tags: PyReadonlyArray1<'_, u32>,
    ) -> PyResult<Self> {
        let v = vertices.as_array();
        let t = triangles.as_array();
        let g = tags.as_array();
        if v.shape()[1] != 3 {
            return Err(PyValueError::new_err(format!(
                "vertices must have shape [N, 3]; got [{}, {}]",
                v.shape()[0],
                v.shape()[1]
            )));
        }
        if t.shape()[1] != 3 {
            return Err(PyValueError::new_err(format!(
                "triangles must have shape [M, 3]; got [{}, {}]",
                t.shape()[0],
                t.shape()[1]
            )));
        }
        if t.shape()[0] != g.shape()[0] {
            return Err(PyValueError::new_err(format!(
                "triangles and tags must have the same length; got {} and {}",
                t.shape()[0],
                g.shape()[0]
            )));
        }
        let verts: Vec<Vector3<f64>> = v
            .rows()
            .into_iter()
            .map(|row| Vector3::new(row[0], row[1], row[2]))
            .collect();
        let tris: Vec<[u32; 3]> = t
            .rows()
            .into_iter()
            .map(|row| [row[0], row[1], row[2]])
            .collect();
        let tags_vec: Vec<u32> = g.iter().copied().collect();
        let inner =
            RustTriMesh::new(verts, tris, tags_vec).map_err(crate::errors::yee_mesh_to_py)?;
        Ok(Self { inner })
    }

    fn n_tris(&self) -> usize {
        self.inner.n_tris()
    }

    #[getter]
    fn vertices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let n = self.inner.vertices.len();
        let mut buf = vec![0.0_f64; n * 3];
        for (i, v) in self.inner.vertices.iter().enumerate() {
            buf[i * 3] = v.x;
            buf[i * 3 + 1] = v.y;
            buf[i * 3 + 2] = v.z;
        }
        numpy::ndarray::Array2::from_shape_vec((n, 3), buf)
            .expect("buffer length matches shape by construction")
            .into_pyarray(py)
    }

    #[getter]
    fn triangles<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let m = self.inner.triangles.len();
        let mut buf = vec![0_u32; m * 3];
        for (i, t) in self.inner.triangles.iter().enumerate() {
            buf[i * 3] = t[0];
            buf[i * 3 + 1] = t[1];
            buf[i * 3 + 2] = t[2];
        }
        numpy::ndarray::Array2::from_shape_vec((m, 3), buf)
            .expect("buffer length matches shape by construction")
            .into_pyarray(py)
    }

    #[getter]
    fn tags<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.inner.tags.clone().into_pyarray(py)
    }
}
