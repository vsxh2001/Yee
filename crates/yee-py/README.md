# yee-py

Python bindings for the Yee electromagnetic-simulation workspace. `yee-py`
wraps `yee-core`, `yee-mesh`, `yee-mom`, and `yee-io` via
[PyO3 0.28](https://pyo3.rs) + [maturin 1.10](https://www.maturin.rs) and
ships as an `abi3-py310` wheel, giving a single build that is
forward-compatible across Python 3.10 through 3.14. The Rust extension module
is `yee._yee`; the public surface is re-exported from the top-level `yee`
package.

---

## Install (local development)

### pip + maturin

```bash
# from the repo root
python -m venv .venv
source .venv/bin/activate
pip install maturin numpy pytest
cd crates/yee-py
maturin develop --release
python -c "import yee; print(yee.__version__)"
```

### uv

```bash
uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install maturin numpy pytest
cd crates/yee-py
maturin develop --release
```

The wheel is `abi3-py310`: one build, forward-compatible across Python 3.10–3.14.
No recompile is needed when upgrading the interpreter within that range.

---

## Jupyter example

```python
import numpy as np
import yee

# ── 1. Build a 2-triangle planar mesh ─────────────────────────────────────────
# Tag 1 = port triangle, tag 2 = ground triangle (port edge is where the two
# tags meet, per the yee-mom basis convention).
vertices = np.array(
    [[0.0, 0.0, 0.0],
     [0.1, 0.0, 0.0],
     [0.1, 0.1, 0.0],
     [0.0, 0.1, 0.0]],
    dtype=np.float64,
)
triangles = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
tags      = np.array([1, 2],                 dtype=np.uint32)

mesh = yee.TriMesh(vertices, triangles, tags)
print(f"mesh: {mesh.n_tris()} triangles")

# ── 2. Define a linear frequency sweep ────────────────────────────────────────
freq = yee.FreqRange(1.0e9, 1.5e9, 3)
print(f"freq: {freq.start_hz/1e9:.2f}–{freq.stop_hz/1e9:.2f} GHz, "
      f"{freq.n_points} points")

# ── 3. Run the solver ─────────────────────────────────────────────────────────
# PlanarMoM.run() is the Phase 0 stub: it always returns
#   RuntimeError("unimplemented: PlanarMoM::run not implemented in phase 0")
# When Phase 1.0 mom-dipole physics land, run() returns a real SParameters
# and the try/except simply stops firing.
solver = yee.PlanarMoM()

try:
    s = solver.run(mesh, freq)
except RuntimeError as exc:
    print(f"[phase 0 stub] solver raised: {exc}")
    s = None

# ── 4. Touchstone write + read round-trip ─────────────────────────────────────
n_ports = 1
freq_hz = np.array([1.0e9, 1.25e9, 1.5e9], dtype=np.float64)
data    = np.array(
    [[[0.10 + 0.20j]],
     [[0.15 + 0.25j]],
     [[0.20 + 0.30j]]],
    dtype=np.complex128,
)
file_dict = {
    "z0":        50.0,
    "freq_unit": "Hz",
    "format":    "RI",
    "n_ports":   n_ports,
    "freq_hz":   freq_hz,
    "data":      data,
    "comments":  ["yee-py Jupyter example"],
}

path = "/tmp/example.s1p"
yee.touchstone.write(path, file_dict)
parsed = yee.touchstone.read(path)

print(f"touchstone round-trip: n_ports={parsed['n_ports']}, "
      f"z0={parsed['z0']} Ω, data shape={parsed['data'].shape}, "
      f"dtype={parsed['data'].dtype}")

if s is not None:
    s.write_touchstone("/tmp/solver_output.s1p", z0=50.0)
    print(f"solver output: n_ports={s.n_ports}, "
          f"freq shape={s.freq_hz.shape}, data shape={s.data.shape}")
```

Expected output under the Phase 0 stub:

```
mesh: 2 triangles
freq: 1.00–1.50 GHz, 3 points
[phase 0 stub] solver raised: unimplemented: PlanarMoM::run not implemented in phase 0
touchstone round-trip: n_ports=1, z0=50.0 Ω, data shape=(3, 1, 1), dtype=complex128
```

---

## Notebook helpers

Three numpy-friendly convenience functions for plotting S-parameters
straight from a notebook. Each accepts a 1-D `complex128` array (typically
one S₁₁ trace over frequency) and returns a numpy array:

| Function | Returns |
|----------|---------|
| `yee.s11_db(s)` | `float64[N]` — 20·log₁₀(\|S\|), clamped to −200 dB at exact zero. |
| `yee.s11_phase(s)` | `float64[N]` — phase angle in degrees, in (−180, 180]. |
| `yee.smith_xy(s)` | `float64[N, 2]` — Cartesian (Re, Im) for Smith-chart plotting. |

```python
import numpy as np
import matplotlib.pyplot as plt
import yee

freq_hz = np.linspace(1.0e9, 2.0e9, 101)
# `s11` would normally come from `solver.run(...).data[:, 0, 0]`
s11 = 0.3 * np.exp(-1j * 2 * np.pi * freq_hz / 1.5e9)

fig, (ax_db, ax_smith) = plt.subplots(1, 2, figsize=(10, 4))
ax_db.plot(freq_hz / 1e9, yee.s11_db(s11))
ax_db.set_xlabel("Frequency (GHz)"); ax_db.set_ylabel("|S₁₁| (dB)")

xy = yee.smith_xy(s11)
ax_smith.plot(xy[:, 0], xy[:, 1])
ax_smith.set_aspect("equal"); ax_smith.set_title("Smith")
```

---

## Public API at a glance

| Symbol | Description |
|--------|-------------|
| `yee.TriMesh(vertices, triangles, tags)` | Build a surface mesh from numpy arrays. `vertices`: `float64[N, 3]`; `triangles`: `uint32[M, 3]`; `tags`: `uint32[M]`. |
| `yee.FreqRange(start_hz, stop_hz, n_points)` | Validate and iterate a linear frequency sweep. Raises `ValueError` on bad inputs. |
| `yee.PlanarMoM().run(mesh, freq)` | Solve planar MoM; returns `SParameters`. Currently raises `RuntimeError("unimplemented: ...")` until Phase 1.0 physics land. |
| `yee.SParameters` | Container: `.n_ports` (int), `.freq_hz` (`float64[F]`), `.data` (`complex128[F, N, N]`), `.write_touchstone(path, z0)`. |
| `yee.touchstone.read(path)` | Read a Touchstone v1.1 file; returns a dict with keys `z0`, `freq_unit`, `format`, `n_ports`, `freq_hz`, `data`, `comments`. |
| `yee.touchstone.write(path, file_dict)` | Write a Touchstone v1.1 file from the same dict schema. |

---

## Error model

All Rust errors map to Python exceptions in `crates/yee-py/src/errors.rs`.

### `yee_core::Error` → Python

| Rust variant | Python exception |
|---|---|
| `Error::Invalid(msg)` | `ValueError` |
| `Error::Numerical(msg)` | `RuntimeError` (prefixed `"numerical: "`) |
| `Error::Unimplemented(msg)` | `RuntimeError` (prefixed `"unimplemented: "`) |
| `Error::Io(msg)` | `IOError` |

### `yee_io::Error` → Python

| Rust variant | Python exception |
|---|---|
| `Error::Io(msg)` | `IOError` |
| `Error::TouchstoneParse { line, col, msg }` | `ValueError` with `"touchstone parse at line N, col M: msg"` |
| `Error::InvalidFile(msg)` | `ValueError` |
| `Error::NotEnabled(feature)` | `RuntimeError` |

### `yee_mesh::Error` → Python

| Rust variant | Python exception |
|---|---|
| `Error::Invalid(msg)` | `ValueError` |
| `Error::NotEnabled` | `RuntimeError` |
| `Error::Gmsh(code)` | `RuntimeError` with `"gmsh error code N"` |

---

## Status

**Phase 1.frontend.0** (this crate) is complete:

- PyO3 0.28 + maturin 1.10 extension module, `abi3-py310` wheel.
- `TriMesh`, `FreqRange`, `PlanarMoM`, `SParameters` bindings.
- Touchstone v1.1 read/write submodule.
- 14 pytests passing; `cargo clippy -p yee-py --all-targets -- -D warnings` clean.
- CI `python-bindings` job wired.

**Phase 1.0 mom-dipole physics** is a parallel track. The NEC-4 finite-radius
reference for a half-wave dipole at `a = 5 mm` is **87 + j41 Ω**. Iteration to
within tolerance is in progress. Whatever the upstream `yee-mom` solver
produces, `yee-py` exposes without transformation or calibration — accuracy
is owned by the physics layer, not the bindings layer.

**`PlanarMoM.run` today** returns
`RuntimeError("unimplemented: PlanarMoM::run not implemented in phase 0")` for
every input. The test suite asserts this contract explicitly. When Phase 1.0
lands and the stub is replaced, the exception stops firing and the example
above falls through to the success path without changes.

---

## License

GPL-3.0-or-later. See [`LICENSE`](../../LICENSE).
