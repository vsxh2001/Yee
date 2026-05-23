//! `yee.eigensolver` submodule — `PyTriMesh2D` + `PyNumericalCrossSection`.
//!
//! Python wrappers for [`yee_mesh::TriMesh2D`] and
//! [`yee_mom::ports::NumericalCrossSection`], i.e. the Phase 1.3.1.1
//! wave-port cross-section eigensolver. The shape of the binding mirrors
//! the Rust API: build a 2-D triangle mesh, hand it to a
//! `NumericalCrossSection` together with per-tag `eps_r` / `mu_r` dicts,
//! then call `.solve(freq_hz)` and read back `beta` / `z_w`.
//!
//! ## `Complex64` ↔ Python `complex`
//!
//! PyO3 0.28's automatic `num_complex::Complex64 <-> PyComplex` conversion
//! is gated behind the `num-complex` cargo feature, which this workspace
//! does **not** enable (see `pyo3` deps in the root `Cargo.toml`). All
//! complex extraction / production therefore goes through the explicit
//! [`PyComplex::from_doubles`] / `c.real(); c.imag()` path in the helpers
//! below. This is the same pattern that would be used by hand even with
//! the feature enabled — the feature only saves a few lines of plumbing.
//!
//! ## Submodule registration
//!
//! Mirrors the `yee.touchstone` pattern: we insert the module into
//! `sys.modules` from `lib.rs` so that
//! `from yee.eigensolver import NumericalCrossSection` works in addition to
//! the attribute-access form `yee.eigensolver.NumericalCrossSection`.

use num_complex::Complex64;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict};
use std::collections::HashMap;
use yee_mesh::{MaterialTag, TriMesh2D};
use yee_mom::ports::NumericalCrossSection;

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTriMesh2D>()?;
    m.add_class::<PyNumericalCrossSection>()?;
    Ok(())
}

/// Python wrapper for [`yee_mesh::TriMesh2D`].
///
/// Constructed from plain Python lists / tuples rather than numpy arrays
/// because the underlying Rust types are `Vec<[f64; 2]>` and
/// `Vec<[usize; 3]>` — list-of-tuple is the closer match and avoids
/// dtype-juggling for a small mesh-of-72-triangles use case. If a future
/// caller needs to push thousands of triangles, an `from_numpy`
/// classmethod can be added without breaking the current API.
#[pyclass(name = "TriMesh2D", module = "yee.eigensolver")]
pub struct PyTriMesh2D {
    pub(crate) inner: TriMesh2D,
}

#[pymethods]
impl PyTriMesh2D {
    #[new]
    #[pyo3(signature = (vertices, triangles, vertex_material=None, triangle_material=None))]
    fn new(
        vertices: Vec<(f64, f64)>,
        triangles: Vec<(usize, usize, usize)>,
        vertex_material: Option<Vec<MaterialTag>>,
        triangle_material: Option<Vec<MaterialTag>>,
    ) -> PyResult<Self> {
        let verts: Vec<[f64; 2]> = vertices.into_iter().map(|(x, y)| [x, y]).collect();
        let tris: Vec<[usize; 3]> = triangles.into_iter().map(|(a, b, c)| [a, b, c]).collect();
        let inner = TriMesh2D::new(verts, tris, vertex_material, triangle_material)
            .map_err(crate::errors::yee_mesh_to_py)?;
        Ok(Self { inner })
    }

    /// Number of vertices.
    fn n_verts(&self) -> usize {
        self.inner.n_verts()
    }

    /// Number of triangles.
    fn n_tris(&self) -> usize {
        self.inner.n_tris()
    }

    /// Signed area of triangle `idx` (positive for CCW winding).
    fn area(&self, idx: usize) -> PyResult<f64> {
        if idx >= self.inner.n_tris() {
            return Err(PyValueError::new_err(format!(
                "triangle index {idx} out of range (n_tris = {})",
                self.inner.n_tris()
            )));
        }
        Ok(self.inner.area(idx))
    }

    /// Centroid of triangle `idx` as `(x, y)`.
    fn centroid(&self, idx: usize) -> PyResult<(f64, f64)> {
        if idx >= self.inner.n_tris() {
            return Err(PyValueError::new_err(format!(
                "triangle index {idx} out of range (n_tris = {})",
                self.inner.n_tris()
            )));
        }
        let c = self.inner.centroid(idx);
        Ok((c[0], c[1]))
    }
}

/// Extract a `dict[int, complex]` into `HashMap<MaterialTag, Complex64>`.
///
/// See module docs — the PyO3 `num-complex` feature is not enabled, so the
/// conversion is done manually: read `PyComplex.real / .imag`, accept a
/// bare `float` as `re + 0j` for ergonomic real-only material specs.
fn extract_complex_map(
    dict: &Bound<'_, PyDict>,
    name: &str,
) -> PyResult<HashMap<MaterialTag, Complex64>> {
    let mut out = HashMap::with_capacity(dict.len());
    for (k, v) in dict.iter() {
        let key: MaterialTag = k.extract().map_err(|e| {
            PyValueError::new_err(format!("{name} key must be int (MaterialTag): {e}"))
        })?;
        let value: Complex64 = if let Ok(c) = v.cast::<PyComplex>() {
            Complex64::new(c.real(), c.imag())
        } else if let Ok(f) = v.extract::<f64>() {
            // Allow real-valued shorthand: `{0: 1.0}` means `1 + 0j`.
            Complex64::new(f, 0.0)
        } else {
            return Err(PyValueError::new_err(format!(
                "{name} value for key {key} must be complex or float"
            )));
        };
        out.insert(key, value);
    }
    Ok(out)
}

/// Python wrapper for [`yee_mom::ports::NumericalCrossSection`].
///
/// The Rust struct stores its `beta` / `z_w` caches as `Option<Complex64>`;
/// before `solve` runs they are `None`, after a successful solve they are
/// populated. The Python view exposes them as `complex | None` getters.
#[pyclass(name = "NumericalCrossSection", module = "yee.eigensolver")]
pub struct PyNumericalCrossSection {
    inner: NumericalCrossSection,
}

#[pymethods]
impl PyNumericalCrossSection {
    /// Build a cross-section descriptor from a `TriMesh2D` and per-tag
    /// material dicts.
    ///
    /// `eps_r` and `mu_r` are `dict[int, complex]` mapping a
    /// `MaterialTag` (matching `mesh.triangle_material`) to the
    /// region's complex relative permittivity / permeability.
    #[new]
    fn new(
        mesh: &PyTriMesh2D,
        eps_r: &Bound<'_, PyDict>,
        mu_r: &Bound<'_, PyDict>,
    ) -> PyResult<Self> {
        let eps_r_map = extract_complex_map(eps_r, "eps_r")?;
        let mu_r_map = extract_complex_map(mu_r, "mu_r")?;
        let inner = NumericalCrossSection::new(mesh.inner.clone(), eps_r_map, mu_r_map);
        Ok(Self { inner })
    }

    /// Run the 2-D Nedelec eigensolve at `freq_hz`.
    ///
    /// Raises `RuntimeError` if the Rust solve returns a numerical /
    /// unimplemented error, or `ValueError` if the inputs are invalid.
    fn solve(&mut self, freq_hz: f64) -> PyResult<()> {
        self.inner.solve(freq_hz).map_err(crate::errors::yee_to_py)
    }

    /// Cached propagation constant β from the most recent solve, or
    /// `None` before any solve has run.
    #[getter]
    fn beta<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyComplex>> {
        self.inner
            .beta
            .map(|b| PyComplex::from_doubles(py, b.re, b.im))
    }

    /// Cached wave impedance Z_w from the most recent solve, or
    /// `None` before any solve has run.
    #[getter]
    fn z_w<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyComplex>> {
        self.inner
            .z_w
            .map(|z| PyComplex::from_doubles(py, z.re, z.im))
    }

    /// Cached dominant-mode **longitudinal** field `E_z` from the most
    /// recent solve as a `list[complex]` (one amplitude per global mesh
    /// vertex, in `mesh` vertex order), or `None` before any solve has
    /// run.
    ///
    /// These are the nodal-Lagrange vertex-DoF amplitudes of the
    /// quasi-TEM dominant mode's longitudinal electric field, scattered
    /// out from the interior-vertex DoF set with Dirichlet (PEC) boundary
    /// vertices set to `0`. On a **homogeneous** (air-filled) guide the
    /// dominant mode is purely transverse, so this is ~zero; on an
    /// **inhomogeneous** (dielectric-loaded / microstrip) cross-section it
    /// carries the genuine longitudinal field that couples through the
    /// dielectric interface. Mirrors the
    /// [`NumericalCrossSection::mode_profile_ez`] Rust field; real-valued
    /// on the lossless path but typed `complex` to match `beta` / `z_w`.
    #[getter]
    fn mode_profile_ez<'py>(&self, py: Python<'py>) -> Option<Vec<Bound<'py, PyComplex>>> {
        self.inner.mode_profile_ez.as_ref().map(|ez| {
            ez.iter()
                .map(|c| PyComplex::from_doubles(py, c.re, c.im))
                .collect()
        })
    }
}
