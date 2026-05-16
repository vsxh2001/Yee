# Phase 1.frontend.0 — Python Bindings — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a `yee-py` workspace crate that exposes `TriMesh`, `FreqRange`, `PlanarMoM`, `SParameters`, and the `touchstone` module to Python via PyO3 0.28 + maturin 1.10, with `abi3-py310` wheels, numpy zero-copy / copy ergonomics, and a pytest suite.

**Architecture:** Approach A — separate `crates/yee-py/` workspace member, `crate-type = ["cdylib", "rlib"]`. PyO3 attributes isolated from library crates. `maturin develop --release` builds + installs the wheel into the current Python venv. CI runs the bindings in a dedicated `python-bindings` job.

**Tech Stack:** Rust 1.88, PyO3 0.28, numpy 0.28 crate, maturin 1.10, Python 3.10+ (abi3-py310), pytest 7+, numpy>=1.26.

**Companion spec:** `docs/superpowers/specs/2026-05-16-phase-1-frontend-0-python-bindings-design.md`

**Parallelism:** Runs on a separate branch `feature/phase-1-frontend-0-py` from a separate worktree. Does NOT depend on Phase 1.0 mom-dipole physics completion — exposes whatever `yee-mom` currently produces.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `Cargo.toml` (root) | Modify | Add `yee-py` workspace member; add `pyo3`, `numpy` to `[workspace.dependencies]` if missing |
| `crates/yee-py/Cargo.toml` | Create | `cdylib` + `rlib`, PyO3 + numpy deps |
| `crates/yee-py/pyproject.toml` | Create | maturin build backend, project metadata |
| `crates/yee-py/README.md` | Create | install + Jupyter example |
| `crates/yee-py/src/lib.rs` | Create | `#[pymodule]` registration |
| `crates/yee-py/src/errors.rs` | Create | error-mapping helpers |
| `crates/yee-py/src/trimesh.rs` | Create | `PyTriMesh` |
| `crates/yee-py/src/freq.rs` | Create | `PyFreqRange` |
| `crates/yee-py/src/solver.rs` | Create | `PyPlanarMoM` |
| `crates/yee-py/src/sparams.rs` | Create | `PySParameters` |
| `crates/yee-py/src/touchstone.rs` | Create | `read`/`write` Python functions |
| `crates/yee-py/python/yee/__init__.py` | Create | Pure-Python re-export shim |
| `crates/yee-py/tests/conftest.py` | Create | pytest fixtures (none needed at v0.1; placeholder) |
| `crates/yee-py/tests/test_trimesh.py` | Create | TriMesh construction + validation |
| `crates/yee-py/tests/test_freq.py` | Create | FreqRange construction + iter |
| `crates/yee-py/tests/test_solver.py` | Create | PlanarMoM run produces SParameters |
| `crates/yee-py/tests/test_touchstone.py` | Create | Round-trip through `touchstone.read`/`write` |
| `.github/workflows/ci.yml` | Modify | Add `python-bindings` job |
| `crates/yee-py/.gitignore` | Create | Ignore `.venv/`, `target/`, `*.so`, `__pycache__/` |

---

## Conventions

- Single worktree for this plan: `/home/hadassi/Code/Yee/worktrees/py-bindings` on branch `feature/phase-1-frontend-0-py`. The orchestrator creates the worktree before Task 1.
- Each task is one logical change with TDD shape where possible. Tests live in `tests/test_*.py` and run via `pytest`. Rust-side unit tests are not used for FFI (pytest is the source of truth).
- Every task ends with a commit. Frequent small commits.
- `numpy`-crate API (`0.28`) function names may differ in minor revisions; the plan provides the documented form. Each agent dispatch must verify against `https://docs.rs/numpy/0.28` and `https://pyo3.rs/v0.28` and adjust if names have drifted.

---

## Task 0: Worktree setup (orchestrator-only)

**Files:** none (git operations only).

- [ ] **Step 1: Create the worktree**

```bash
cd /home/hadassi/Code/Yee
BASE_SHA=$(git rev-parse main)
echo "BASE_SHA=$BASE_SHA"   # record this base
git worktree add -b feature/phase-1-frontend-0-py worktrees/py-bindings "$BASE_SHA"
git worktree list
```

Expected: a new worktree at `worktrees/py-bindings` on branch `feature/phase-1-frontend-0-py`.

- [ ] **Step 2: Pre-flight Python tooling**

```bash
which python3 || sudo apt-get install -y python3 python3-pip python3-venv
python3 --version    # expect 3.10+
python3 -m venv worktrees/py-bindings/.venv
source worktrees/py-bindings/.venv/bin/activate
pip install --upgrade pip
pip install maturin pytest numpy
maturin --version    # expect 1.10+
```

If the system Python is < 3.10, install `python3.11` via the host package manager before continuing.

---

## Task 1: Workspace wiring + `yee-py` crate skeleton

**Files:**
- Modify: `Cargo.toml` (root)
- Create: `crates/yee-py/Cargo.toml`
- Create: `crates/yee-py/src/lib.rs`
- Create: `crates/yee-py/.gitignore`

- [ ] **Step 1: Add `yee-py` to the workspace members**

In root `Cargo.toml`:

```toml
[workspace]
resolver = "2"
members = [
    "crates/yee-core",
    "crates/yee-cuda",
    "crates/yee-mesh",
    "crates/yee-mom",
    "crates/yee-fdtd",
    "crates/yee-io",
    "crates/yee-cli",
    "crates/yee-py",
]
```

Verify `pyo3` and `numpy` are present in `[workspace.dependencies]` (they should be, from Phase 0 pre-flight). If not, add them:

```toml
pyo3  = { version = "0.28", features = ["abi3-py310", "extension-module"] }
numpy = "0.28"
```

- [ ] **Step 2: Create `crates/yee-py/Cargo.toml`**

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

- [ ] **Step 3: Create `crates/yee-py/src/lib.rs` with the empty pymodule**

```rust
//! Python bindings for Yee electromagnetic simulation.
//!
//! See `crates/yee-py/README.md` for install + usage. The pymodule itself is
//! named `_yee` and is wrapped by a pure-Python package `yee` under
//! `crates/yee-py/python/yee/` so that future `.pyi` stubs and convenience
//! helpers can sit alongside the compiled extension.

use pyo3::prelude::*;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
```

- [ ] **Step 4: Create `crates/yee-py/.gitignore`**

```
.venv/
__pycache__/
*.so
*.dylib
*.pyd
target/
build/
dist/
*.egg-info/
.pytest_cache/
```

- [ ] **Step 5: Verify the crate compiles under the workspace**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings
. "$HOME/.cargo/env"
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache
cargo build -p yee-py --release 2>&1 | tail -5
```

Expected: `yee-py` builds. The cdylib output appears under `target/release/libyee.so` (Linux) or `target/release/yee.dll` (Windows).

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml crates/yee-py/
git commit -m "yee-py: add workspace crate skeleton with empty pymodule"
```

---

## Task 2: maturin build configuration + pure-Python shim

**Files:**
- Create: `crates/yee-py/pyproject.toml`
- Create: `crates/yee-py/python/yee/__init__.py`
- Create: `crates/yee-py/python/yee/py.typed`

- [ ] **Step 1: Create `crates/yee-py/pyproject.toml`**

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

- [ ] **Step 2: Create `crates/yee-py/python/yee/__init__.py`**

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import __version__

__all__ = ["__version__"]
```

(This minimal version is expanded in later tasks as classes land.)

- [ ] **Step 3: Create `crates/yee-py/python/yee/py.typed`** (empty file marker for PEP 561 type-hint support, even though we ship no stubs in v0.1):

```
```

(File is intentionally empty.)

- [ ] **Step 4: First `maturin develop` smoke build**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings/crates/yee-py
source /home/hadassi/Code/Yee/worktrees/py-bindings/.venv/bin/activate
. "$HOME/.cargo/env"
export CARGO_TARGET_DIR=/home/hadassi/Code/Yee/worktrees/py-bindings/target
export RUSTC_WRAPPER=sccache
maturin develop --release 2>&1 | tail -10
python -c "import yee; print(yee.__version__)"
```

Expected: maturin reports success; the Python `import yee` succeeds and prints `0.0.0`.

- [ ] **Step 5: Commit**

```bash
git add crates/yee-py/pyproject.toml crates/yee-py/python/
git commit -m "yee-py: maturin pyproject + Python re-export shim"
```

---

## Task 3: Error mapping module

**Files:**
- Create: `crates/yee-py/src/errors.rs`
- Modify: `crates/yee-py/src/lib.rs` (declare `mod errors;`)

- [ ] **Step 1: Add `mod errors;` to `crates/yee-py/src/lib.rs`**

Update `src/lib.rs` to:

```rust
//! Python bindings for Yee electromagnetic simulation.

use pyo3::prelude::*;

mod errors;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
```

- [ ] **Step 2: Create `crates/yee-py/src/errors.rs`**

```rust
//! Conversion from yee-* `Error` types to PyO3 `PyErr` instances.

use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::PyErr;

/// Map a `yee_core::Error` to a Python exception.
///
/// - `Invalid` → `ValueError`
/// - `Numerical` → `RuntimeError` prefixed with `numerical:`
/// - `Unimplemented` → `RuntimeError` prefixed with `unimplemented:`
/// - `Io` → `IOError`
pub fn yee_to_py(err: yee_core::Error) -> PyErr {
    match err {
        yee_core::Error::Invalid(msg) => PyValueError::new_err(msg),
        yee_core::Error::Numerical(msg) => {
            PyRuntimeError::new_err(format!("numerical: {msg}"))
        }
        yee_core::Error::Unimplemented(msg) => {
            PyRuntimeError::new_err(format!("unimplemented: {msg}"))
        }
        yee_core::Error::Io(msg) => PyIOError::new_err(msg),
    }
}

/// Map a `yee_io::Error` to a Python exception.
///
/// - `TouchstoneParse { line, col, msg }` → `ValueError` with position
/// - `Io` → `IOError`
/// - `NotEnabled` → `RuntimeError`
/// - `InvalidFile` → `ValueError`
pub fn io_to_py(err: yee_io::Error) -> PyErr {
    match err {
        yee_io::Error::Io(msg) => PyIOError::new_err(msg),
        yee_io::Error::TouchstoneParse { line, col, msg } => PyValueError::new_err(format!(
            "touchstone parse at line {line}, col {col}: {msg}"
        )),
        yee_io::Error::NotEnabled(feature) => {
            PyRuntimeError::new_err(format!("yee-io feature `{feature}` not enabled"))
        }
        yee_io::Error::InvalidFile(msg) => PyValueError::new_err(msg),
    }
}

/// Map a `yee_mesh::Error` to a Python exception.
///
/// - `Invalid` → `ValueError`
/// - `NotEnabled` → `RuntimeError`
/// - `Gmsh(code)` → `RuntimeError` with the error code
pub fn yee_mesh_to_py(err: yee_mesh::Error) -> PyErr {
    use yee_mesh::Error as E;
    match err {
        E::Invalid(msg) => PyValueError::new_err(msg),
        E::NotEnabled => PyRuntimeError::new_err("yee-mesh `gmsh` feature not enabled"),
        E::Gmsh(code) => PyRuntimeError::new_err(format!("gmsh error code {code}")),
    }
}
```

- [ ] **Step 3: Verify the crate still builds + the empty pymodule still works**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings/crates/yee-py
source /home/hadassi/Code/Yee/worktrees/py-bindings/.venv/bin/activate
. "$HOME/.cargo/env"
export CARGO_TARGET_DIR=/home/hadassi/Code/Yee/worktrees/py-bindings/target
maturin develop --release 2>&1 | tail -5
python -c "import yee; print(yee.__version__)"
```

Expected: still succeeds.

- [ ] **Step 4: Commit**

```bash
git add crates/yee-py/src/errors.rs crates/yee-py/src/lib.rs
git commit -m "yee-py: error mapping helpers (yee_core/yee_io/yee_mesh → PyErr)"
```

---

## Task 4: `PyTriMesh`

**Files:**
- Create: `crates/yee-py/src/trimesh.rs`
- Modify: `crates/yee-py/src/lib.rs` (add `mod trimesh;` and `add_class`)
- Modify: `crates/yee-py/python/yee/__init__.py` (re-export `TriMesh`)
- Create: `crates/yee-py/tests/test_trimesh.py`

- [ ] **Step 1: Write the failing pytest**

`crates/yee-py/tests/test_trimesh.py`:

```python
"""Tests for yee.TriMesh."""

import numpy as np
import pytest

import yee


def test_trimesh_constructs_from_numpy():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([1], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    assert mesh.n_tris() == 1


def test_trimesh_rejects_bad_vertex_shape():
    v = np.array([[0, 0]], dtype=np.float64)  # missing z column
    t = np.array([[0, 0, 0]], dtype=np.uint32)
    g = np.array([0], dtype=np.uint32)
    with pytest.raises(ValueError, match="vertices must have shape"):
        yee.TriMesh(v, t, g)


def test_trimesh_rejects_bad_triangle_shape():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    t = np.array([[0, 1]], dtype=np.uint32)  # missing third index
    g = np.array([0], dtype=np.uint32)
    with pytest.raises(ValueError, match="triangles must have shape"):
        yee.TriMesh(v, t, g)


def test_trimesh_rejects_length_mismatch():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([0, 0], dtype=np.uint32)  # length 2 vs 1 triangle
    with pytest.raises(ValueError, match="length"):
        yee.TriMesh(v, t, g)


def test_trimesh_getters_return_arrays_with_correct_shape():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64)
    t = np.array([[0, 1, 2]], dtype=np.uint32)
    g = np.array([42], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    assert mesh.vertices.shape == (3, 3)
    assert mesh.triangles.shape == (1, 3)
    assert mesh.tags.shape == (1,)
    assert int(mesh.tags[0]) == 42
```

- [ ] **Step 2: Run; expect ImportError because `TriMesh` not yet exposed**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings/crates/yee-py
source /home/hadassi/Code/Yee/worktrees/py-bindings/.venv/bin/activate
pytest tests/test_trimesh.py 2>&1 | tail -15
```

Expected: failure / `AttributeError: module 'yee' has no attribute 'TriMesh'`.

- [ ] **Step 3: Create `crates/yee-py/src/trimesh.rs`**

```rust
//! `PyTriMesh` — Python wrapper around `yee_mesh::TriMesh`.

use nalgebra::Vector3;
use numpy::ndarray::Array2;
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
        let inner = RustTriMesh::new(verts, tris, tags_vec)
            .map_err(crate::errors::yee_mesh_to_py)?;
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
        Array2::from_shape_vec((n, 3), buf)
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
        Array2::from_shape_vec((m, 3), buf)
            .unwrap()
            .into_pyarray_bound(py)
    }

    #[getter]
    fn tags<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.inner.tags.clone().into_pyarray_bound(py)
    }
}
```

If the `numpy` 0.28 minor on your host names its conversions differently (e.g. `IntoPyArray::into_pyarray` instead of `into_pyarray_bound`), adjust per `https://docs.rs/numpy/0.28`. The semantics — "convert an owned ndarray into a Python-managed PyArray" — are stable across the minor.

- [ ] **Step 4: Wire into `crates/yee-py/src/lib.rs`**

```rust
//! Python bindings for Yee electromagnetic simulation.

use pyo3::prelude::*;

mod errors;
mod trimesh;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<trimesh::PyTriMesh>()?;
    Ok(())
}
```

- [ ] **Step 5: Update `crates/yee-py/python/yee/__init__.py`**

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import TriMesh, __version__

__all__ = ["TriMesh", "__version__"]
```

- [ ] **Step 6: Rebuild + run tests**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings/crates/yee-py
source /home/hadassi/Code/Yee/worktrees/py-bindings/.venv/bin/activate
maturin develop --release 2>&1 | tail -10
pytest tests/test_trimesh.py -v 2>&1 | tail -20
```

Expected: 5 tests pass.

- [ ] **Step 7: Commit**

```bash
git add crates/yee-py/src/trimesh.rs crates/yee-py/src/lib.rs crates/yee-py/python/yee/__init__.py crates/yee-py/tests/test_trimesh.py
git commit -m "yee-py: PyTriMesh with numpy zero-copy constructors + accessors"
```

---

## Task 5: `PyFreqRange`

**Files:**
- Create: `crates/yee-py/src/freq.rs`
- Modify: `crates/yee-py/src/lib.rs`
- Modify: `crates/yee-py/python/yee/__init__.py`
- Create: `crates/yee-py/tests/test_freq.py`

- [ ] **Step 1: Write the failing tests**

`crates/yee-py/tests/test_freq.py`:

```python
"""Tests for yee.FreqRange."""

import numpy as np
import pytest

import yee


def test_freqrange_constructs():
    f = yee.FreqRange(1.0e9, 2.0e9, 21)
    assert f.start_hz == 1.0e9
    assert f.stop_hz == 2.0e9
    assert f.n_points == 21


def test_freqrange_iter_endpoints_exact():
    f = yee.FreqRange(1.0e9, 2.0e9, 3)
    pts = f.iter()
    assert pts.shape == (3,)
    assert pts[0] == 1.0e9
    assert pts[-1] == 2.0e9


def test_freqrange_rejects_inverted_band():
    with pytest.raises(ValueError):
        yee.FreqRange(2.0e9, 1.0e9, 5)


def test_freqrange_rejects_zero_points():
    with pytest.raises(ValueError):
        yee.FreqRange(1.0e9, 2.0e9, 0)


def test_freqrange_rejects_non_finite():
    with pytest.raises(ValueError):
        yee.FreqRange(float("inf"), 2.0e9, 5)
```

- [ ] **Step 2: Run; expect failure (no FreqRange)**

```bash
pytest tests/test_freq.py -v 2>&1 | tail -10
```

- [ ] **Step 3: Create `crates/yee-py/src/freq.rs`**

```rust
//! `PyFreqRange` — Python wrapper around `yee_core::FreqRange`.

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
        values.into_pyarray_bound(py)
    }
}
```

- [ ] **Step 4: Wire into `crates/yee-py/src/lib.rs`**

```rust
//! Python bindings for Yee electromagnetic simulation.

use pyo3::prelude::*;

mod errors;
mod freq;
mod trimesh;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<trimesh::PyTriMesh>()?;
    m.add_class::<freq::PyFreqRange>()?;
    Ok(())
}
```

- [ ] **Step 5: Update `crates/yee-py/python/yee/__init__.py`**

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import FreqRange, TriMesh, __version__

__all__ = ["FreqRange", "TriMesh", "__version__"]
```

- [ ] **Step 6: Rebuild + run**

```bash
maturin develop --release
pytest tests/test_freq.py -v 2>&1 | tail -15
```

Expected: 5 tests pass.

- [ ] **Step 7: Commit**

```bash
git add crates/yee-py/src/freq.rs crates/yee-py/src/lib.rs crates/yee-py/python/yee/__init__.py crates/yee-py/tests/test_freq.py
git commit -m "yee-py: PyFreqRange with validation + iter() returning np.ndarray"
```

---

## Task 6: `PySParameters`

**Files:**
- Create: `crates/yee-py/src/sparams.rs`
- Modify: `crates/yee-py/src/lib.rs`
- Modify: `crates/yee-py/python/yee/__init__.py`

(No tests added here — `SParameters` is exercised by `test_solver.py` and `test_touchstone.py` in later tasks.)

- [ ] **Step 1: Create `crates/yee-py/src/sparams.rs`**

```rust
//! `PySParameters` — Python wrapper around `yee_mom::SParameters`.

use numpy::ndarray::Array3;
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
    fn n_ports(&self) -> usize {
        self.inner.n_ports
    }

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
        Array3::from_shape_vec((f, n, n), buf)
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

- [ ] **Step 2: Wire into `crates/yee-py/src/lib.rs`**

```rust
//! Python bindings for Yee electromagnetic simulation.

use pyo3::prelude::*;

mod errors;
mod freq;
mod sparams;
mod trimesh;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<trimesh::PyTriMesh>()?;
    m.add_class::<freq::PyFreqRange>()?;
    m.add_class::<sparams::PySParameters>()?;
    Ok(())
}
```

- [ ] **Step 3: Update `crates/yee-py/python/yee/__init__.py`**

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import FreqRange, SParameters, TriMesh, __version__

__all__ = ["FreqRange", "SParameters", "TriMesh", "__version__"]
```

- [ ] **Step 4: Verify build**

```bash
maturin develop --release 2>&1 | tail -5
python -c "import yee; print(yee.SParameters)"
```

Expected: `<class 'builtins.SParameters'>` or similar — class is reachable, no instance test here.

- [ ] **Step 5: Commit**

```bash
git add crates/yee-py/src/sparams.rs crates/yee-py/src/lib.rs crates/yee-py/python/yee/__init__.py
git commit -m "yee-py: PySParameters with freq_hz/data accessors + write_touchstone"
```

---

## Task 7: `PyPlanarMoM` and end-to-end run test

**Files:**
- Create: `crates/yee-py/src/solver.rs`
- Modify: `crates/yee-py/src/lib.rs`
- Modify: `crates/yee-py/python/yee/__init__.py`
- Create: `crates/yee-py/tests/test_solver.py`

- [ ] **Step 1: Write the failing pytest**

`crates/yee-py/tests/test_solver.py`:

```python
"""Tests for yee.PlanarMoM."""

import numpy as np

import yee


def test_planar_mom_runs_and_returns_sparams():
    v = np.array(
        [[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0]], dtype=np.float64
    )
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([1, 2], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    s = yee.PlanarMoM().run(mesh, freq)
    assert s.n_ports == 1
    assert s.freq_hz.shape == (3,)
    assert s.data.shape == (3, 1, 1)
    assert s.data.dtype == np.complex128


def test_planar_mom_returns_numerical_error_without_port():
    import pytest

    v = np.array(
        [[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0]], dtype=np.float64
    )
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([0, 0], dtype=np.uint32)  # no port tags
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    with pytest.raises(RuntimeError, match="numerical"):
        yee.PlanarMoM().run(mesh, freq)
```

- [ ] **Step 2: Run; expect failure**

```bash
pytest tests/test_solver.py -v 2>&1 | tail -10
```

Expected: `AttributeError: module 'yee' has no attribute 'PlanarMoM'`.

- [ ] **Step 3: Create `crates/yee-py/src/solver.rs`**

```rust
//! `PyPlanarMoM` — Python wrapper around `yee_mom::PlanarMoM`.

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
```

- [ ] **Step 4: Wire into `crates/yee-py/src/lib.rs`**

```rust
//! Python bindings for Yee electromagnetic simulation.

use pyo3::prelude::*;

mod errors;
mod freq;
mod solver;
mod sparams;
mod trimesh;

#[pymodule]
fn _yee(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<trimesh::PyTriMesh>()?;
    m.add_class::<freq::PyFreqRange>()?;
    m.add_class::<sparams::PySParameters>()?;
    m.add_class::<solver::PyPlanarMoM>()?;
    Ok(())
}
```

- [ ] **Step 5: Update `crates/yee-py/python/yee/__init__.py`**

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import FreqRange, PlanarMoM, SParameters, TriMesh, __version__

__all__ = ["FreqRange", "PlanarMoM", "SParameters", "TriMesh", "__version__"]
```

- [ ] **Step 6: Rebuild + run**

```bash
maturin develop --release 2>&1 | tail -5
pytest tests/test_solver.py -v 2>&1 | tail -15
```

Expected: 2 tests pass. Note: this test depends on `yee-mom`'s current state (whatever physics accuracy `PlanarMoM` produces is what Python sees). It does NOT assert mom-001's ±5% gate — that's Track A's problem.

- [ ] **Step 7: Commit**

```bash
git add crates/yee-py/src/solver.rs crates/yee-py/src/lib.rs crates/yee-py/python/yee/__init__.py crates/yee-py/tests/test_solver.py
git commit -m "yee-py: PyPlanarMoM run wrapper + end-to-end pytest"
```

---

## Task 8: `touchstone` submodule (read + write)

**Files:**
- Create: `crates/yee-py/src/touchstone.rs`
- Modify: `crates/yee-py/src/lib.rs`
- Modify: `crates/yee-py/python/yee/__init__.py`
- Create: `crates/yee-py/tests/test_touchstone.py`

- [ ] **Step 1: Write the failing pytest**

`crates/yee-py/tests/test_touchstone.py`:

```python
"""Tests for yee.touchstone read/write."""

import numpy as np

import yee


def test_sparams_write_and_read_roundtrip(tmp_path):
    v = np.array(
        [[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0]], dtype=np.float64
    )
    t = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    g = np.array([1, 2], dtype=np.uint32)
    mesh = yee.TriMesh(v, t, g)
    freq = yee.FreqRange(1.0e9, 1.5e9, 3)
    s = yee.PlanarMoM().run(mesh, freq)

    path = tmp_path / "test.s1p"
    s.write_touchstone(str(path), 50.0)

    file = yee.touchstone.read(str(path))
    assert file["n_ports"] == 1
    assert file["z0"] == 50.0
    assert file["freq_unit"] == "Hz"
    assert file["format"] == "RI"
    assert len(file["freq_hz"]) == 3
    assert file["data"].shape == (3, 1, 1)


def test_touchstone_read_rejects_missing_file(tmp_path):
    import pytest

    path = tmp_path / "nope.s1p"
    with pytest.raises(IOError):
        yee.touchstone.read(str(path))
```

- [ ] **Step 2: Run; expect failure (no touchstone submodule)**

```bash
pytest tests/test_touchstone.py -v 2>&1 | tail -10
```

- [ ] **Step 3: Create `crates/yee-py/src/touchstone.rs`**

```rust
//! `yee.touchstone` submodule — Python wrappers for read/write.

use numpy::ndarray::Array3;
use numpy::{IntoPyArray, PyArray1, PyArray3, PyReadonlyArray1, PyReadonlyArray3};
use num_complex::Complex64;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;
use yee_io::touchstone::{self, File as TsFile, Format, FreqUnit};

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
    d.set_item("freq_unit", freq_unit_to_str(file.freq_unit))?;
    d.set_item("format", format_to_str(file.format))?;
    d.set_item("n_ports", file.n_ports)?;
    let freq = file.freq_hz.clone().into_pyarray_bound(py);
    d.set_item("freq_hz", freq)?;
    let f = file.data.len();
    let n = file.n_ports;
    let mut buf: Vec<Complex64> = Vec::with_capacity(f * n * n);
    for row in &file.data {
        buf.extend_from_slice(row);
    }
    let arr = Array3::from_shape_vec((f, n, n), buf)
        .unwrap()
        .into_pyarray_bound(py);
    d.set_item("data", arr)?;
    d.set_item("comments", file.comments.clone())?;
    Ok(d)
}

#[pyfunction]
fn write(path: PathBuf, file: &Bound<'_, PyDict>) -> PyResult<()> {
    let z0: f64 = file
        .get_item("z0")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing key: z0"))?
        .extract()?;
    let freq_unit_str: String = file
        .get_item("freq_unit")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing key: freq_unit"))?
        .extract()?;
    let format_str: String = file
        .get_item("format")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing key: format"))?
        .extract()?;
    let n_ports: usize = file
        .get_item("n_ports")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing key: n_ports"))?
        .extract()?;
    let freq_hz_obj = file
        .get_item("freq_hz")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing key: freq_hz"))?;
    let freq_hz: PyReadonlyArray1<f64> = freq_hz_obj.extract()?;
    let data_obj = file
        .get_item("data")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing key: data"))?;
    let data: PyReadonlyArray3<Complex64> = data_obj.extract()?;
    let comments: Vec<String> = match file.get_item("comments")? {
        Some(c) => c.extract()?,
        None => Vec::new(),
    };

    let freq_hz_vec: Vec<f64> = freq_hz.as_array().iter().copied().collect();
    let data_arr = data.as_array();
    let f = data_arr.shape()[0];
    if data_arr.shape()[1] != n_ports || data_arr.shape()[2] != n_ports {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "data shape {:?} inconsistent with n_ports {}",
            data_arr.shape(),
            n_ports
        )));
    }
    let mut data_vec: Vec<Vec<Complex64>> = Vec::with_capacity(f);
    for k in 0..f {
        let mut row = Vec::with_capacity(n_ports * n_ports);
        for i in 0..n_ports {
            for j in 0..n_ports {
                row.push(data_arr[(k, i, j)]);
            }
        }
        data_vec.push(row);
    }

    let ts = TsFile {
        z0,
        freq_unit: str_to_freq_unit(&freq_unit_str)?,
        format: str_to_format(&format_str)?,
        n_ports,
        freq_hz: freq_hz_vec,
        data: data_vec,
        comments,
    };
    touchstone::write(&path, &ts).map_err(crate::errors::io_to_py)
}

fn freq_unit_to_str(u: FreqUnit) -> &'static str {
    match u {
        FreqUnit::Hz => "Hz",
        FreqUnit::KHz => "kHz",
        FreqUnit::MHz => "MHz",
        FreqUnit::GHz => "GHz",
    }
}

fn format_to_str(f: Format) -> &'static str {
    match f {
        Format::RealImag => "RI",
        Format::MagnitudeAngle => "MA",
        Format::DecibelAngle => "DB",
    }
}

fn str_to_freq_unit(s: &str) -> PyResult<FreqUnit> {
    match s.to_ascii_lowercase().as_str() {
        "hz" => Ok(FreqUnit::Hz),
        "khz" => Ok(FreqUnit::KHz),
        "mhz" => Ok(FreqUnit::MHz),
        "ghz" => Ok(FreqUnit::GHz),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown freq_unit: {other}"
        ))),
    }
}

fn str_to_format(s: &str) -> PyResult<Format> {
    match s.to_ascii_uppercase().as_str() {
        "RI" => Ok(Format::RealImag),
        "MA" => Ok(Format::MagnitudeAngle),
        "DB" => Ok(Format::DecibelAngle),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown format: {other}"
        ))),
    }
}
```

The exact `yee_io::touchstone::Format` and `FreqUnit` variant names must be verified against `crates/yee-io/src/touchstone.rs` from the merge base — if Agent D's Phase 0 work named them differently (e.g. `Format::Db` not `Format::DecibelAngle`), adjust the match arms here. The spec uses the names from the Phase 0 cap that ships on `main`.

- [ ] **Step 4: Wire into `crates/yee-py/src/lib.rs`**

```rust
//! Python bindings for Yee electromagnetic simulation.

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
    let ts_mod = PyModule::new_bound(py, "touchstone")?;
    touchstone::register(&ts_mod)?;
    m.add_submodule(&ts_mod)?;
    Ok(())
}
```

- [ ] **Step 5: Update `crates/yee-py/python/yee/__init__.py`**

```python
"""Yee electromagnetic simulation — Python bindings."""

from yee._yee import (
    FreqRange,
    PlanarMoM,
    SParameters,
    TriMesh,
    __version__,
    touchstone,
)

__all__ = [
    "FreqRange",
    "PlanarMoM",
    "SParameters",
    "TriMesh",
    "__version__",
    "touchstone",
]
```

- [ ] **Step 6: Rebuild + run**

```bash
maturin develop --release 2>&1 | tail -5
pytest tests/test_touchstone.py -v 2>&1 | tail -15
```

Expected: 2 tests pass.

- [ ] **Step 7: Commit**

```bash
git add crates/yee-py/src/touchstone.rs crates/yee-py/src/lib.rs crates/yee-py/python/yee/__init__.py crates/yee-py/tests/test_touchstone.py
git commit -m "yee-py: yee.touchstone.read/write with PyDict ↔ yee_io::File mapping"
```

---

## Task 9: README + Jupyter example

**Files:**
- Create: `crates/yee-py/README.md`

- [ ] **Step 1: Write `crates/yee-py/README.md`**

````markdown
# yee — Python bindings for Yee electromagnetic simulation

Minimal Python interface to the [Yee](https://github.com/yee-em/yee) workspace.
Wraps the Rust core via PyO3 0.28 and ships as `abi3-py310` wheels.

## Install (local development)

```bash
# from a fresh checkout
python -m venv .venv
source .venv/bin/activate
pip install maturin numpy pytest
cd crates/yee-py
maturin develop --release
python -c "import yee; print(yee.__version__)"
```

## Jupyter example

```python
import numpy as np
import yee

# Build a two-triangle test mesh with a delta-gap port between the
# differently-tagged triangles.
vertices = np.array(
    [[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0]], dtype=np.float64
)
triangles = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
tags = np.array([1, 2], dtype=np.uint32)
mesh = yee.TriMesh(vertices, triangles, tags)

# Sweep from 1 to 2 GHz, 21 points
freq = yee.FreqRange(1.0e9, 2.0e9, 21)

# Run
solver = yee.PlanarMoM()
s = solver.run(mesh, freq)
print(s.n_ports, s.freq_hz.shape, s.data.shape)

# Export Touchstone
s.write_touchstone("toy.s1p", 50.0)

# Read it back
file = yee.touchstone.read("toy.s1p")
print(file["z0"], file["format"], file["data"].shape)
```

## Accuracy

The Phase 1.0 mom-001 (half-wave dipole) integration in `yee-mom` currently
reports `Z_in` within roughly 20% of the Balanis reference (target: ±5%).
Track A development continues to close that gap. Python users see whatever
the underlying solver produces — no calibration is applied at this layer.

## License

GPL v3.0 or later.
````

- [ ] **Step 2: Commit**

```bash
git add crates/yee-py/README.md
git commit -m "yee-py: README with install + Jupyter example"
```

---

## Task 10: CI job

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Append a `python-bindings` job to `.github/workflows/ci.yml`**

Add this job (after the existing `lint-test` job):

```yaml
  python-bindings:
    runs-on: ubuntu-latest
    needs: lint-test
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          toolchain: "1.88"
          components: rustfmt, clippy
      - uses: Swatinem/rust-cache@v2
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: install build deps
        run: |
          python -m pip install --upgrade pip
          pip install maturin pytest numpy
      - name: maturin develop
        working-directory: crates/yee-py
        run: maturin develop --release
      - name: pytest
        working-directory: crates/yee-py
        run: pytest tests/
```

- [ ] **Step 2: Verify the YAML is well-formed locally**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings
which yamllint && yamllint .github/workflows/ci.yml || echo "yamllint absent; skipping"
```

(Skipped if yamllint not installed — the workflow runner will catch syntax errors on push.)

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: add python-bindings job (maturin develop + pytest)"
```

---

## Task 11: Full workspace regression + merge + push

**Files:** none (verification + git operations).

- [ ] **Step 1: Run the full workspace gates from the worktree**

```bash
cd /home/hadassi/Code/Yee/worktrees/py-bindings
. "$HOME/.cargo/env"
source .venv/bin/activate
export CARGO_TARGET_DIR="$PWD/target"
export RUSTC_WRAPPER=sccache

cargo check --workspace --no-default-features
cargo test --workspace --no-default-features
cargo clippy --workspace --all-targets --no-default-features -- -D warnings
cargo fmt --check --all
cargo doc --workspace --no-default-features --no-deps

# Python side
cd crates/yee-py
maturin develop --release
pytest tests/
```

Expected: every step exits 0; pytest reports at least the test counts: trimesh (5) + freq (5) + solver (2) + touchstone (2) = **14 tests**.

- [ ] **Step 2: Lane check before merge**

```bash
cd /home/hadassi/Code/Yee
git -C worktrees/py-bindings diff --stat main..HEAD
```

Expected: every path is in one of:
- `Cargo.toml` (root, workspace member addition)
- `Cargo.lock`
- `crates/yee-py/**`
- `.github/workflows/ci.yml`

Anything else → reject before merge.

- [ ] **Step 3: Merge into main**

```bash
cd /home/hadassi/Code/Yee
git merge --no-ff feature/phase-1-frontend-0-py \
  -m "Merge Phase 1.frontend.0: Python bindings (yee-py)"
git log --oneline -10
```

- [ ] **Step 4: Cleanup**

```bash
git worktree remove worktrees/py-bindings
git branch -d feature/phase-1-frontend-0-py
git worktree list
```

- [ ] **Step 5: Push + tag**

```bash
git push origin main
git tag -a phase-1-frontend-0 -m "Phase 1.frontend.0 — Python bindings v0.1"
git push origin phase-1-frontend-0
```

---

## Self-Review

**1. Spec coverage:**

- Spec §1 in scope items: TriMesh (Task 4), FreqRange (Task 5), PlanarMoM (Task 7), SParameters (Task 6), touchstone (Task 8). Error mapping (Task 3). README (Task 9). CI (Task 10). All covered.
- Spec §3 crate layout: Tasks 1+2 lay the workspace + maturin config. Tasks 3–8 populate `src/*.rs` files. All files in the spec layout are created by some task.
- Spec §4 Python API surface: every signature in the spec appears in code blocks in Tasks 4–8.
- Spec §5 FFI bindings: source files match spec one-for-one. Task 8 explicitly notes the `yee_io::Format` / `FreqUnit` variant verification step.
- Spec §6 build/test/distribution: Task 2 establishes maturin; Tasks 4/5/7/8 add pytest; Task 10 adds CI. Distribution (PyPI) is explicitly deferred per spec.
- Spec §7 test strategy: all four pytest files exist and contain the assertions from the spec.
- Spec §8 risks: R1 (pin numpy=0.28) addressed by Task 1 dep block. R2 (venv) addressed in Task 0 + README. R3 (free-threaded ABI) is a deferred concern, documented in README via tech-stack reference. R7 (mom-001 accuracy passthrough) is documented in the README. R8/R9 (PyO3 / numpy crate API variability) is called out in the relevant tasks.

**2. Placeholder scan:**

- Tasks 4 and 8 mention "verify the exact numpy/yee_io variant names against installed minor versions" — this is forward-looking validation guidance, not a TBD. Each task gives the canonical names AND tells the implementer how to adjust if drift is detected. Acceptable.
- No `// TODO` / `fill in details` / "appropriate validation" patterns.

**3. Type consistency:**

- `PyTriMesh` (Task 4), `PyFreqRange` (Task 5), `PySParameters` (Task 6), `PyPlanarMoM` (Task 7) — all `#[pyclass(name = "...", module = "yee._yee")]`, all exposed via `m.add_class::<...>()?` in `lib.rs`, all re-exported in `python/yee/__init__.py`. Names match.
- `crate::errors::yee_to_py` (Task 3) consumed by Tasks 5, 6, 7, 8. Type checks.
- `crate::errors::yee_mesh_to_py` (Task 3) consumed by Task 4.
- `crate::errors::io_to_py` (Task 3) consumed by Task 8.

No issues that require revision.
