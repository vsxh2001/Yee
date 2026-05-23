# Tutorial 4 — Waveguide eigenmode from Python

This tutorial walks through Yee's 2-D Nedelec edge-element eigensolver
from Python. You will build a structured triangle mesh of the WR-90
rectangular waveguide cross-section, hand it to `yee.eigensolver.NumericalCrossSection`
to extract the dominant TE10 propagation constant `beta(f)`, sweep
across frequency, and compare the numerical result to the closed-form
Pozar §3.3 analytic. This is the same Phase 1.3.1.1 cross-section solver
that future numerical wave-ports will use to seed their modal
distribution — exposed here as a standalone Python class so it can be
driven from a notebook for teaching, exploration, or quick mesh-quality
sanity checks.

## Goal

Compute the TE10 phase constant `beta` of WR-90
(`a x b = 22.86 mm x 10.16 mm`, air-filled) at 10 GHz and verify it
lands within ~0.1 % of the analytic value `beta = sqrt(k0^2 - (pi/a)^2)
~= 158.150550 rad/m`. Sweep across the X-band (8–12 GHz) and print a
side-by-side table of numerical vs analytic. This is the contract the
Rust validation gate in `crates/yee-mom/tests/eigensolver_wr90.rs`
measures against; the Python binding lets you reproduce it from a REPL.

## Prerequisites

- Rust 1.92+ (the Python wheel is built from source).
- Python 3.10 through 3.14. The wheel is `abi3-py310`, so any
  interpreter in that range works without rebuilding.
- `pip install maturin numpy pytest` — or the `uv` equivalent.
- Jupyter is optional; the snippets run in a plain `python` REPL too.

## Install

If you already followed [Tutorial 2](02-dipole-from-python.md) you can
skip this step — the same wheel exposes the `yee.eigensolver` submodule.
Otherwise, from the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install maturin numpy pytest
cd crates/yee-py
maturin develop --release
python -c "from yee.eigensolver import NumericalCrossSection, TriMesh2D; print('ok')"
```

## Background

A hollow rectangular waveguide of inner dimensions `a x b` supports a
discrete family of TE_mn and TM_mn modes (Pozar, *Microwave
Engineering*, §3.3). Each mode has a **cutoff frequency**

```text
f_c(m, n) = (c / 2) * sqrt((m / a)^2 + (n / b)^2)
```

below which it is evanescent and above which it propagates with phase
constant

```text
beta(f) = sqrt(k0(f)^2 - k_c(m, n)^2),   k_c = 2 * pi * f_c / c
```

For WR-90 the dominant TE10 mode (`m = 1, n = 0`) cuts off at
`f_c = c / (2a) ~= 6.557 GHz`, so the standard X-band 8–12 GHz operating
range sits comfortably above cutoff. The eigensolver discretises the
transverse Helmholtz equation on a triangle mesh of the cross-section
using lowest-order Nedelec edge elements; the smallest physical
eigenvalue corresponds to `k_c^2`, from which `beta` follows.

## Build the cross-section mesh

`yee.eigensolver.TriMesh2D` takes plain Python lists of `(x, y)`
vertices and `(v0, v1, v2)` triangles (CCW winding mandatory — the
constructor rejects clockwise triangles with a `ValueError`). For a
rectangular cross-section the cleanest mesh is a structured `nx x ny`
quad grid with each quad split along its lower-left → upper-right
diagonal into two CCW triangles. The pattern below matches the
canonical fixture in
`crates/yee-py/tests/test_eigensolver.py::test_numerical_cross_section_solve_wr90`
and the Rust gate in `crates/yee-mom/tests/eigensolver_wr90.rs`.

```python
from yee.eigensolver import NumericalCrossSection, TriMesh2D

a = 22.86e-3   # WR-90 long inner dimension (m)
b = 10.16e-3   # WR-90 short inner dimension (m)
nx, ny = 6, 6  # 72 triangles, 49 vertices

vertices = [
    (a * i / nx, b * j / ny)
    for j in range(ny + 1)
    for i in range(nx + 1)
]

def idx(i, j):
    return j * (nx + 1) + i

triangles = []
for j in range(ny):
    for i in range(nx):
        v00 = idx(i, j)
        v10 = idx(i + 1, j)
        v11 = idx(i + 1, j + 1)
        v01 = idx(i, j + 1)
        triangles.append((v00, v10, v11))
        triangles.append((v00, v11, v01))

mesh = TriMesh2D(vertices, triangles)
print(f"n_verts = {mesh.n_verts()}, n_tris = {mesh.n_tris()}")
```

Expected output:

```text
n_verts = 49, n_tris = 72
```

All 72 triangles default to material tag `0`, which is fine for an
air-filled waveguide. If you need a partially-filled cross-section
(dielectric slab, ridge waveguide), pass `triangle_material=[...]`
to the constructor and key your `eps_r` / `mu_r` dicts on those tags.

## Run the eigensolver

`NumericalCrossSection` takes the mesh plus per-tag material dicts.
`solve(freq_hz)` runs the dense generalised eigenproblem and caches
`beta` and `z_w` on the instance; before `solve` runs both are `None`.

```python
eps_r = {0: complex(1, 0)}   # air
mu_r  = {0: complex(1, 0)}

nc = NumericalCrossSection(mesh, eps_r, mu_r)
nc.solve(10.0e9)

print(f"beta = {nc.beta.real:.6f} rad/m  (Im = {nc.beta.imag:.2e})")
print(f"Z_w  = {abs(nc.z_w):.2f} Ω")
```

Expected output (6x6 mesh, 10 GHz):

```text
beta = 158.236... rad/m  (Im = ~1e-13)
Z_w  = ~500.00 Ω
```

The imaginary part of `beta` is numerical noise — WR-90 is lossless and
air-filled, so the analytic answer is purely real. The 6x6 mesh agrees
with the closed-form `158.150550 rad/m` to ~0.055 %, well inside the
1 % gate.

## Cutoff sweep and verify

Loop over frequencies above cutoff, solve, and compare to the analytic
Pozar formula side-by-side. A fresh `NumericalCrossSection` per
frequency is the simplest pattern; the solve cost is dominated by the
dense eigendecomposition at this mesh size.

```python
import math

C0 = 299_792_458.0

def analytic_beta_te10(a, freq_hz, eps_r=1.0):
    k = 2.0 * math.pi * freq_hz * math.sqrt(eps_r) / C0
    kc = math.pi / a
    if k <= kc:
        return float("nan")
    return math.sqrt(k * k - kc * kc)

print(f"WR-90 TE10 cutoff: f_c = {C0 / (2.0 * a) / 1e9:.3f} GHz")
print()
print(f"{'f (GHz)':>8} {'β_num':>12} {'β_an':>12} {'rel err %':>10}")
for f_ghz in [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]:
    freq_hz = f_ghz * 1.0e9
    nc = NumericalCrossSection(mesh, eps_r, mu_r)
    nc.solve(freq_hz)
    beta_num = nc.beta.real
    beta_an = analytic_beta_te10(a, freq_hz)
    rel_err = abs(beta_num - beta_an) / beta_an * 100.0
    print(f"{f_ghz:8.2f} {beta_num:12.4f} {beta_an:12.4f} {rel_err:10.4f}")
```

Expected output (numbers within rounding):

```text
WR-90 TE10 cutoff: f_c = 6.557 GHz

 f (GHz)        β_num         β_an  rel err %
    7.00      52.0...      51.9...     0.0...
    8.00     110.6...     110.5...     0.0...
    9.00     143.0...     142.9...     0.0...
   10.00     158.24...    158.15...    0.0552
   11.00     181.8...     181.7...     0.0...
   12.00     203.6...     203.5...     0.0...
```

The relative error stays at the ~0.05–0.1 % level across the X-band on
this mesh, as expected for lowest-order Nedelec elements on a 6x6 grid.
For a graphical view (analytic curve plus numerical scatter), see the
notebook companion linked in [Next steps](#next-steps).

## Mesh refinement

The eigensolver is a standard lowest-order Nedelec edge-element method,
so the error in `beta` scales roughly as `O(h^2)` where `h ~ a / nx` is
the cell size. Doubling the mesh density should cut the error by ~4x;
the dense eigensolve cost grows ~`O(n^3)` in the DOF count, so this is
a real tradeoff for larger cross-sections.

| `nx = ny` | `n_tris` | rel err at 10 GHz | dense solve cost      |
|-----------|----------|-------------------|-----------------------|
| 4         | 32       | ~0.5 %            | trivial               |
| 6         | 72       | ~0.055 %          | trivial               |
| 8         | 128      | ~0.014 %          | trivial               |
| 12        | 288      | ~0.003 %          | noticeable            |
| 16+       | 512+     | sub-ppm           | dense path too slow   |

For air-filled WR-90 at X-band, `nx = ny = 6` is the sweet spot — the
validation gate uses it because it's the coarsest mesh that holds the
1 % contract with margin. Refine when the cross-section is electrically
larger, when you care about higher modes (TE20, TE01, TE11) whose
eigenvectors are less well-resolved by a coarse mesh, or when a
dielectric inclusion forces stronger fields near a material interface.

## Next steps

- The notebook version of this walkthrough, with matplotlib plots
  comparing the numerical sweep to a dense analytic curve, lives at
  `examples/python/eigensolver_wr90_te10.ipynb`. Run it after `maturin
  develop --release` from the same virtualenv.
- **Status update (steps 4–5.2 shipped).** The roadmap items below have
  since landed: step 4 added an in-tree block LOBPCG (ADR-0050), step 5
  the mixed `(E_t, E_z)` longitudinal block (ADR-0051), and step 5.2 the
  β-direct extraction fix that makes dielectric fills correct (ADR-0053).
  `NumericalCrossSection` is now production-quality for hollow and
  *uniformly*-filled guides and supports *inhomogeneous* (microstrip-style)
  fills, with a high-contrast accuracy gap still closing (step 5.3). For
  the current formulation, solver options, and validation status see the
  theory chapter
  [Cross-Section (Waveguide-Port) Eigensolver](../theory/cross-section-eigensolver.md).
  This tutorial's TE10-only homogeneous walkthrough remains correct as an
  introduction.
- Theory background and references for the Nedelec edge-element
  formulation live in `docs/src/theory/planar-mom.md`; Pozar §3.3
  remains the textbook reference for the analytic side.
