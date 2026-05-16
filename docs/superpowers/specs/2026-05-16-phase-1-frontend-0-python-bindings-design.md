# Phase 1.frontend.0 — Python Bindings — Design

**Date:** 2026-05-16
**Status:** Ready for writing-plans
**Repo base SHA at design time:** `684763f` (Phase 0 + Phase 1.0 dipole plan)
**Predecessor:** Phase 0 walking-skeleton; Phase 1.0 mom-dipole sub-project (in progress in parallel on `feature/phase-1-0-mom-dipole`).

This spec is the first frontend sub-project under Phase 1. It produces a minimal Python bindings crate (`yee-py`) exposing the workspace's core types (`TriMesh`, `FreqRange`, `PlanarMoM`, `SParameters`, `touchstone`) to Python via PyO3 0.28 + maturin 1.10 with `abi3-py310` wheels. It runs in parallel with the Phase 1.0 mom-dipole physics sub-project — no shared files.

---

## 1. Scope & Success Criteria

### In scope

- New workspace member `crates/yee-py/` built with PyO3 0.28 + maturin 1.10.
- Python module `yee` exporting:
  - `yee.TriMesh(vertices: np.ndarray[N,3,f64], triangles: np.ndarray[M,3,u32], tags: np.ndarray[M,u32])`
  - `yee.FreqRange(start_hz: f64, stop_hz: f64, n_points: int)`
  - `yee.PlanarMoM().run(mesh, freq) -> SParameters`
  - `yee.SParameters` with `n_ports`, `freq_hz`, `data`, `write_touchstone(path, z0)`
  - `yee.touchstone.read(path) -> dict`, `yee.touchstone.write(path, file)`
- Numpy zero-copy where practical; copy where layout mismatch (e.g. `Vec<Vector3<f64>>` → 2D ndarray).
- Error mapping: `yee_core::Error::Invalid` → `ValueError`; `yee_io::Error::TouchstoneParse` → `ValueError` with position; `yee_*::Error::Io` → `IOError`; everything else → `RuntimeError` with the underlying message.
- `pytest` suite covering construction, validation, run, touchstone round-trip, numpy buffer protocol, error mapping.
- `abi3-py310` wheels built locally via `maturin develop --release`.
- New CI job `python-bindings` in `.github/workflows/ci.yml` running `maturin develop --release` + `pytest`.
- `crates/yee-py/README.md` documents install + Jupyter example.

### Out of scope

- PyPI publishing (deferred).
- Type-hint stubs (`.pyi`).
- Mesh helpers (`yee.shapes.dipole`) — a separate sub-project.
- Smith-chart / matplotlib helpers — separate sub-project.
- GUI (egui) — separate sub-project.
- Phase 1.0 mom-dipole physics accuracy — orthogonal track.

### Done means

1. `pip install maturin && maturin develop --release` in `crates/yee-py/` exits 0 on Linux with Rust 1.88 + Python 3.10+.
2. `python -c "import yee; print(yee.__version__)"` exits 0.
3. `pytest crates/yee-py/tests/` exits 0 with at least six tests across `test_trimesh.py`, `test_freq.py`, `test_solver.py`, `test_touchstone.py`.
4. `cargo build -p yee-py --release` exits 0 standalone (i.e. without invoking maturin).
5. CI workflow `python-bindings` exits green on every PR and push to `main`.
6. All Phase 0 + Phase 1.0 gates stay green.
7. `crates/yee-py/README.md` covers install + Jupyter example.

### Performance budget (informational)

- `maturin develop --release` cold: < 8 min (pyo3 + numpy + workspace transitive build).
- Warm rebuild after a Rust edit: < 60 s with sccache.

---

## 2. Decisions Locked During Brainstorming

| # | Decision | Choice |
|---|----------|--------|
| D1 | First frontend sub-project | Python bindings (smallest scope, highest leverage for RF engineers using notebooks) |
| D2 | Scope | Solver-only (D2 from Q2 (i)): `TriMesh`, `FreqRange`, `PlanarMoM`, `SParameters`, `touchstone`. No mesh helpers, no plotting helpers. |
| D3 | numpy ↔ Python | Numpy zero-copy (or copy where layout mismatch). The `numpy` crate is already in TECH_STACK. |
| D4 | Architecture | Approach A — separate `yee-py` crate as a workspace member, built with maturin. PyO3 isolated from library crates. |
| D5 | Parallelism with Phase 1.0 mom-dipole track | Run in parallel: different worktree, different crate, no shared files. |

---

## 3. Crate Layout

```
crates/yee-py/
├── Cargo.toml                 # cdylib + rlib; pyo3 features
├── pyproject.toml             # maturin build backend
├── README.md                  # install + Jupyter example
├── src/
│   ├── lib.rs                 # #[pymodule] yee; re-exports + version
│   ├── trimesh.rs             # #[pyclass] TriMesh wrapper
│   ├── freq.rs                # #[pyclass] FreqRange wrapper
│   ├── sparams.rs             # #[pyclass] SParameters wrapper
│   ├── solver.rs              # #[pyclass] PlanarMoM wrapper
│   ├── touchstone.rs          # #[pymodule] yee.touchstone (read/write)
│   └── errors.rs              # yee → PyErr mapping helpers
├── python/
│   └── yee/
│       └── __init__.py        # pure-Python re-export layer for type stubs later
└── tests/
    ├── test_trimesh.py
    ├── test_freq.py
    ├── test_solver.py
    ├── test_touchstone.py
    └── conftest.py            # tmp_path fixtures
```

`Cargo.toml`:

```toml
[package]
name         = "yee-py"
version      = { workspace = true }
edition      = { workspace = true }
rust-version = { workspace = true }
license      = { workspace = true }
repository   = { workspace = true }
description  = "Python bindings for Yee electromagnetic simulation."

[lib]
name       = "yee"
crate-type = ["cdylib", "rlib"]

[dependencies]
yee-core    = { workspace = true }
yee-mesh    = { workspace = true }
yee-mom     = { workspace = true }
yee-io      = { workspace = true }
pyo3        = { workspace = true }
numpy       = { workspace = true }
num-complex = { workspace = true }
nalgebra    = { workspace = true }
thiserror   = { workspace = true }
```

`pyproject.toml`:

```toml
[build-system]
requires      = ["maturin>=1.10,<2.0"]
build-backend = "maturin"

[project]
name            = "yee"
description     = "GPU-accelerated electromagnetic simulation."
requires-python = ">=3.10"
license         = { text = "GPL-3.0-or-later" }
classifiers     = [
    "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Physics",
]
dependencies    = ["numpy>=1.26"]
dynamic         = ["version"]

[project.optional-dependencies]
test = ["pytest>=7"]

[tool.maturin]
manifest-path  = "Cargo.toml"
features       = ["pyo3/extension-module"]
python-source  = "python"
module-name    = "yee._yee"
```

Pure-Python wrapper `python/yee/__init__.py`:

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import (
    TriMesh,
    FreqRange,
    PlanarMoM,
    SParameters,
    touchstone,
    __version__,
)

__all__ = [
    "TriMesh",
    "FreqRange",
    "PlanarMoM",
    "SParameters",
    "touchstone",
    "__version__",
]
```

---

## 4. Python API Surface

```python
import yee
import numpy as np

vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=np.float64)
triangles = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
tags = np.array([1, 2], dtype=np.uint32)
mesh = yee.TriMesh(vertices, triangles, tags)
mesh.n_tris()        # int
mesh.vertices        # np.ndarray[N, 3, f64]
mesh.triangles       # np.ndarray[M, 3, u32]
mesh.tags            # np.ndarray[M, u32]

freq = yee.FreqRange(start_hz=1.0e9, stop_hz=2.0e9, n_points=21)
freq.start_hz        # float
freq.stop_hz         # float
freq.n_points        # int
freq.iter()          # np.ndarray[n_points, f64]

solver = yee.PlanarMoM()
s = solver.run(mesh, freq)

s.n_ports            # int (=1 in Phase 1.0)
s.freq_hz            # np.ndarray[F, f64]
s.data               # np.ndarray[F, N, N, c128]
s.write_touchstone(path, z0=50.0)

file = yee.touchstone.read("dipole.s1p")
# file is a dict: {"z0": 50.0, "freq_unit": "Hz", "format": "RI",
#                  "n_ports": 1, "freq_hz": np.ndarray, "data": np.ndarray, "comments": list[str]}

yee.touchstone.write("out.s1p", file)
```

### Error mapping

```
yee_core::Error::Invalid(msg)        -> ValueError(msg)
yee_core::Error::Numerical(msg)      -> RuntimeError("numerical: " + msg)
yee_core::Error::Unimplemented(msg)  -> RuntimeError("unimplemented: " + msg)
yee_core::Error::Io(msg)             -> IOError(msg)
yee_io::Error::TouchstoneParse{...}  -> ValueError("touchstone parse at line N, col M: msg")
yee_io::Error::Io(msg)               -> IOError(msg)
yee_io::Error::NotEnabled(feature)   -> RuntimeError("yee-io feature `feature` not enabled")
yee_io::Error::InvalidFile(msg)      -> ValueError(msg)
```

Implemented as helper functions `yee_to_py` / `io_to_py` in `src/errors.rs`.

---

## 5. FFI Bindings Layout

### `src/lib.rs`

```rust
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
    m.add_class::<solver::PyPlanarMoM>()?;
    m.add_class::<sparams::PySParameters>()?;
    let touchstone_mod = PyModule::new_bound(py, "touchstone")?;
    touchstone::register(&touchstone_mod)?;
    m.add_submodule(&touchstone_mod)?;
    Ok(())
}
```

### `src/trimesh.rs`

```rust
use nalgebra::Vector3;
use numpy::{IntoPyArray, PyArray2, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
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
                v.shape()[0], v.shape()[1]
            )));
        }
        if t.shape()[1] != 3 {
            return Err(PyValueError::new_err(format!(
                "triangles must have shape [M, 3]; got [{}, {}]",
                t.shape()[0], t.shape()[1]
            )));
        }
        if t.shape()[0] != g.shape()[0] {
            return Err(PyValueError::new_err(format!(
                "triangles and tags must have the same length; got {} and {}",
                t.shape()[0], g.shape()[0]
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
            .unwrap()
            .into_pyarray_bound(py)
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
            .unwrap()
            .into_pyarray_bound(py)
    }

    #[getter]
    fn tags<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.inner.tags.clone().into_pyarray_bound(py)
    }
}
```

(Method signatures and numpy-crate API names are templates against the documented 0.28 surface. The implementer may need to adjust `from_shape_vec` / `into_pyarray_bound` / `rows()` calls for the exact installed minor — consult `https://docs.rs/numpy/0.28` and `https://pyo3.rs/v0.28`.)

### `src/freq.rs`

```rust
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
        let inner = RustFreqRange::new(start_hz, stop_hz, n_points)
            .map_err(crate::errors::yee_to_py)?;
        Ok(Self { inner })
    }

    #[getter] fn start_hz(&self) -> f64 { self.inner.start_hz }
    #[getter] fn stop_hz(&self) -> f64 { self.inner.stop_hz }
    #[getter] fn n_points(&self) -> usize { self.inner.n_points }

    fn iter<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let values: Vec<f64> = self.inner.iter().collect();
        values.into_pyarray_bound(py)
    }
}
```

### `src/solver.rs`

```rust
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
        Self { inner: RustPlanarMoM::default() }
    }

    fn run(
        &self,
        mesh: &crate::trimesh::PyTriMesh,
        freq: &crate::freq::PyFreqRange,
    ) -> PyResult<crate::sparams::PySParameters> {
        let s = self.inner
            .run(&mesh.inner, freq.inner)
            .map_err(crate::errors::yee_to_py)?;
        Ok(crate::sparams::PySParameters { inner: s })
    }
}
```

### `src/sparams.rs`

```rust
use numpy::{IntoPyArray, PyArray1, PyArray3};
use num_complex::Complex64;
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
    fn n_ports(&self) -> usize { self.inner.n_ports }

    #[getter]
    fn freq_hz<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.freq_hz.clone().into_pyarray_bound(py)
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
            .unwrap()
            .into_pyarray_bound(py)
    }

    fn write_touchstone(&self, path: PathBuf, z0: f64) -> PyResult<()> {
        self.inner
            .write_touchstone(&path, z0)
            .map_err(crate::errors::yee_to_py)
    }
}
```

### `src/touchstone.rs`

```rust
use numpy::{IntoPyArray, PyArray1, PyArray3};
use num_complex::Complex64;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;
use yee_io::touchstone::{self, Format, FreqUnit};

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    Ok(())
}

#[pyfunction]
fn read<'py>(py: Python<'py>, path: PathBuf) -> PyResult<Bound<'py, PyDict>> {
    let file = touchstone::read(&path).map_err(crate::errors::io_to_py)?;
    let d = PyDict::new_bound(py);
    d.set_item("z0", file.z0)?;
    d.set_item(
        "freq_unit",
        match file.freq_unit {
            FreqUnit::Hz => "Hz",
            FreqUnit::KHz => "kHz",
            FreqUnit::MHz => "MHz",
            FreqUnit::GHz => "GHz",
        },
    )?;
    d.set_item(
        "format",
        match file.format {
            Format::RealImag => "RI",
            Format::MagnitudeAngle => "MA",
            Format::DecibelAngle => "DB",
        },
    )?;
    d.set_item("n_ports", file.n_ports)?;
    d.set_item("freq_hz", file.freq_hz.into_pyarray_bound(py))?;
    let f = file.data.len();
    let n = file.n_ports;
    let mut buf: Vec<Complex64> = Vec::with_capacity(f * n * n);
    for row in &file.data {
        buf.extend_from_slice(row);
    }
    let arr = numpy::ndarray::Array3::from_shape_vec((f, n, n), buf)
        .unwrap()
        .into_pyarray_bound(py);
    d.set_item("data", arr)?;
    d.set_item("comments", file.comments.clone())?;
    Ok(d)
}

#[pyfunction]
fn write(path: PathBuf, file: &Bound<'_, PyDict>) -> PyResult<()> {
    // Convert PyDict -> yee_io::touchstone::File then call yee_io::touchstone::write.
    // Full conversion in implementation; see implementation plan.
    let _ = (path, file);
    todo!("convert PyDict to yee_io::touchstone::File")
}
```

The `write` direction (`PyDict` → `File`) requires careful field extraction with explicit type checks. The implementation plan provides the full body.

### `src/errors.rs`

```rust
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::PyErr;

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

pub fn yee_mesh_to_py(err: yee_mesh::Error) -> PyErr {
    use yee_mesh::Error as E;
    match err {
        E::Invalid(msg) => PyValueError::new_err(msg),
        E::NotEnabled => PyRuntimeError::new_err("yee-mesh `gmsh` feature not enabled"),
        E::Gmsh(code) => PyRuntimeError::new_err(format!("gmsh error code {code}")),
    }
}
```

---

## 6. Build, Test, Distribution

### Local development

```bash
cd crates/yee-py
python -m venv .venv && source .venv/bin/activate
pip install maturin pytest numpy
maturin develop --release
python -c "import yee; print(yee.__version__)"
pytest tests/
```

### CI

Add a job to `.github/workflows/ci.yml`:

```yaml
  python-bindings:
    runs-on: ubuntu-latest
    needs: lint-test
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          toolchain: "1.88"
      - uses: Swatinem/rust-cache@v2
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: maturin develop
        working-directory: crates/yee-py
        run: |
          python -m pip install --upgrade pip
          pip install maturin pytest numpy
          maturin develop --release
      - name: pytest
        working-directory: crates/yee-py
        run: pytest tests/
```

The job depends on `lint-test` so a fmt/clippy break stops the Python job from running too.

### Distribution

Phase 1.frontend.0 ships local-build only. PyPI publishing is a separate sub-project after at least one downstream user reports success.

---

## 7. Test Strategy

`tests/test_trimesh.py`:

```python
import numpy as np
import pytest
import yee


def test_trimesh_constructs_from_numpy():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([1], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    assert mesh.n_tris() == 1


def test_trimesh_rejects_bad_shape():
    v = np.array([[0, 0]], dtype=np.float64)
    t = np.array([[0, 0, 0]], dtype=np.uint32)
    g = np.array([0], dtype=np.uint32)
    with pytest.raises(ValueError, match="shape"):
        yee.TriMesh(v, t, g)


def test_trimesh_rejects_length_mismatch():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([0, 0], dtype=np.uint32)
    with pytest.raises(ValueError, match="length"):
        yee.TriMesh(v, t, g)
```

`tests/test_freq.py`:

```python
import numpy as np
import pytest
import yee


def test_freqrange_constructs():
    f = yee.FreqRange(1.0e9, 2.0e9, 21)
    assert f.start_hz == 1.0e9
    assert f.n_points == 21


def test_freqrange_iter_endpoints_exact():
    f = yee.FreqRange(1.0e9, 2.0e9, 3)
    pts = f.iter()
    assert pts.shape == (3,)
    assert pts[0] == 1.0e9
    assert pts[-1] == 2.0e9


def test_freqrange_rejects_invalid():
    with pytest.raises(ValueError):
        yee.FreqRange(2.0e9, 1.0e9, 5)
```

`tests/test_solver.py`:

```python
import numpy as np
import yee


def test_planar_mom_runs_and_returns_sparams():
    v = np.array([[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([1, 2], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    s = yee.PlanarMoM().run(mesh, freq)
    assert s.n_ports == 1
    assert s.freq_hz.shape == (3,)
    assert s.data.shape == (3, 1, 1)
    assert s.data.dtype == np.complex128
```

`tests/test_touchstone.py`:

```python
import numpy as np
import yee


def test_sparams_write_and_read_roundtrip(tmp_path):
    v = np.array([[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([1, 2], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    s = yee.PlanarMoM().run(mesh, freq)

    path = tmp_path / "test.s1p"
    s.write_touchstone(str(path), 50.0)

    file = yee.touchstone.read(str(path))
    assert file["n_ports"] == 1
    assert len(file["freq_hz"]) == 3
```

Rust-side smoke test in `crates/yee-py/src/lib.rs`:

```rust
#[cfg(test)]
mod tests {
    #[test]
    fn crate_compiles() {
        // PyO3 modules are tested via pytest; this smoke test ensures the
        // crate compiles standalone without invoking maturin.
    }
}
```

---

## 8. Risk Register

| # | Risk | Likelihood | Impact | Mitigation |
|---|------|-----------|--------|------------|
| R1 | `numpy` crate API changes between PyO3 minors | Med | Med | Pin exact `numpy = "0.28"` matching `pyo3 = "0.28"`. |
| R2 | `maturin develop` requires a Python venv | Confirmed | Low | Document venv setup in `crates/yee-py/README.md`; CI uses `actions/setup-python@v5`. |
| R3 | `abi3-py310` excludes free-threaded ABI (3.13t/3.14t) | Confirmed | Low | TECH_STACK already notes this. Revisit when PEP 803 lands. |
| R4 | Workspace `Cargo.lock` grows with pyo3 + numpy + maturin transitive deps | Confirmed | Low | One-shot cost; sccache absorbs subsequent rebuilds. |
| R5 | `manylinux_2_28` glibc baseline incompatible with older Linux distros | Low | Low | Phase 1.frontend.0 ships local-build only; PyPI deferred. |
| R6 | `crate-type = ["cdylib", "rlib"]` warns downstream Rust consumers | Low | Low | Acceptable — `yee-py` is a Python-bindings crate, not intended as a Rust dependency. |
| R7 | mom-001 19% accuracy gap propagates to Python users | Med | Low | Python bindings expose whatever `yee-mom` produces. Track A debugging continues independently. Documented in README. |
| R8 | First-time PyO3 friction (`Bound<'py, T>` 0.21+ API patterns) | Med | Low | Implementer brief links `https://pyo3.rs/v0.28`; the spec includes templates. |
| R9 | numpy crate's `IntoPyArray` / `into_pyarray_bound` naming varies across minors | Med | Med | Spec notes the templates may need 0.28-specific adjustment; implementer consults docs. |

---

## 9. Next Step

After approval, invoke the `superpowers:writing-plans` skill to produce a task-by-task implementation plan with TDD-shaped steps for each module, the pytest fixtures, the maturin build configuration, and the new CI job.
