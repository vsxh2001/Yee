# Tutorial 2 — Half-wave dipole from Python

This tutorial walks through Yee from a Jupyter notebook. You will install
the `yee-py` extension module via maturin, build a thin-cylinder dipole
TriMesh from numpy, hand it to `PlanarMoM().run(...)` for an impedance
solve, write the result to a Touchstone `.s1p`, and plot the reflection
coefficient on a Smith chart with the `yee plot` CLI subcommand. This is
the studio-tool path: Python in front, Rust solver underneath,
Touchstone files as the lingua franca.

## Goal

Compute the input impedance of a half-wave dipole at its first resonance
(L = 1 m, radius = 5 mm, lateral cylinder surface, delta-gap feed) and
verify it lands near the NEC-4 finite-radius reference of
`Z = 87 + j41 ohms`. Write the result to disk. Plot it.

## Prerequisites

- Rust 1.88+ (the Python wheel is built from source).
- Python 3.10 through 3.14. The wheel is `abi3-py310`, so any
  interpreter in that range works without rebuilding.
- `pip install maturin numpy pytest` — or the `uv` equivalent.
- Jupyter is optional; the snippets run in a plain `python` REPL too.

## Install

From the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install maturin numpy pytest
cd crates/yee-py
maturin develop --release
python -c "import yee; print(yee.__version__)"
```

`maturin develop --release` compiles the Rust extension and drops the
wheel into your active virtualenv. After this, `import yee` works
anywhere in that environment.

## Build the mesh

The NEC-4 reference is for a *finite-radius* dipole, so a wire-thin
representation won't do — we need the lateral surface of a cylinder.
The Rust test fixture lives in
`crates/yee-mom/tests/fixtures/cylinder.rs`; we reproduce its formula in
Python for clarity. The tag convention is what tells `yee-mom` where to
place the delta-gap port: the central edge ring is the boundary between
two differently-tagged triangle bands (tags `1` and `2`), and the basis
builder picks that up automatically.

```python
import math
import numpy as np
import yee

def thin_cylinder(length_m, radius_m, n_axial, n_around):
    assert n_axial >= 2 and n_axial % 2 == 0
    assert n_around >= 3
    dz = length_m / n_axial
    z0 = -length_m / 2.0
    dtheta = (2.0 * math.pi) / n_around

    verts = []
    for i in range(n_axial + 1):
        z = z0 + i * dz
        for j in range(n_around):
            theta = j * dtheta
            verts.append([radius_m * math.cos(theta),
                          radius_m * math.sin(theta),
                          z])

    tris, tags = [], []
    central = n_axial // 2
    for i in range(n_axial):
        for j in range(n_around):
            jn = (j + 1) % n_around
            a = i * n_around + j
            b = i * n_around + jn
            c = (i + 1) * n_around + jn
            d = (i + 1) * n_around + j
            tris.append([a, b, c]); tris.append([a, c, d])
            tag = 1 if i == central - 1 else (2 if i == central else 0)
            tags.append(tag); tags.append(tag)

    verts = np.array(verts, dtype=np.float64)
    tris  = np.array(tris,  dtype=np.uint32)
    tags  = np.array(tags,  dtype=np.uint32)
    return yee.TriMesh(verts, tris, tags)

mesh = thin_cylinder(length_m=1.0, radius_m=0.005,
                     n_axial=12, n_around=32)
print(f"mesh: {mesh.n_tris()} triangles")
```

We use `n_axial=12, n_around=32` here (768 triangles) — coarser than
the validation gate, which uses `(24, 176)` for ~8.4k triangles. Coarse
meshes converge to within roughly 20 % of NEC-4 on this geometry; the
shipped gate refines until the error drops below 5 % on the real part
and 10 % on the imaginary part.

## Run the solver

```python
f0 = 299_792_458.0 / 2.0           # half-wave resonance for L = 1 m
freq = yee.FreqRange(f0, f0 + 1.0, 1)

solver = yee.PlanarMoM()
s = solver.run(mesh, freq)

s11 = s.data[0, 0, 0]
z0 = 50.0
z_in = z0 * (1 + s11) / (1 - s11)
print(f"S11 = {s11:.4f}")
print(f"Z_in = {z_in.real:.2f} + j{z_in.imag:.2f} ohms")
```

Honest expectations:

- The `mom-001` validation gate **passes** at the NEC-4 reference of
  `87 + j41 ohms` with `(24, 176)` triangles; see
  `crates/yee-mom/tests/dipole.rs`. That is the contract the project
  measures against.
- At the coarser `(12, 32)` mesh shown above you will see something in
  the *ballpark* — within roughly 20 % of the reference on a good
  laptop. This is fine for a tutorial; refine to match the gate.
- If you are reading this before Phase 1.0 mom-dipole physics merge,
  `PlanarMoM().run(...)` raises
  `RuntimeError("unimplemented: PlanarMoM::run not implemented in phase 0")`.
  Wrap the call in a `try / except RuntimeError` if you want a clean
  stub-mode run; see `crates/yee-py/README.md`.

## Write Touchstone

```python
s.write_touchstone("dipole.s1p", z0=50.0)
```

The file is Touchstone v1.1, real-imaginary format. Any commercial RF
tool will read it.

## Plot

The Yee CLI ships a `plot` subcommand (Track O) that reads Touchstone
and emits PNG or SVG via the `yee-plotters` crate. For a single
frequency point the dB plot is uninteresting, so a Smith chart is the
right view:

```bash
yee plot dipole.s1p --kind smith --output dipole_smith.png
```

The output format is inferred from the extension; pass
`--output dipole.svg` for SVG instead. Other kinds are `db` (the
classic `|S11|` in dB versus frequency) and `phase`.

For a real sweep, replace the single-point `FreqRange` with something
like `yee.FreqRange(120e6, 180e6, 61)` and re-run. The `db` plot will
show a resonant null near 150 MHz.

## Next

Continue to [Tutorial 3 — FDTD cavity resonance](03-fdtd-cavity.md),
which steps off the planar-MoM rail and runs Yee's time-domain solver
in two boundary configurations.
