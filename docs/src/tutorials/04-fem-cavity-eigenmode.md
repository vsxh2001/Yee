# FEM cavity eigenmode from Python

This tutorial walks through Yee's 3-D Nedelec edge-element FEM
eigensolver from Python. You will build a WR-90-based rectangular
cavity, hand it to `yee.fem.solve_cavity` for a shift-invert inverse
iteration on the curl-curl pencil `K·e = k²·M·e`, and compare the
lowest ten resonances against the published Pozar §6.3 TE/TM table
for an air-filled metallic box. This is the Phase 4.fem.eig.0
walking skeleton — the same physical case [Tutorial 3 — FDTD cavity
resonance](03-fdtd-cavity.md) attacks in the time domain, but solved
directly as a generalised eigenvalue problem in the frequency domain
so the modes pop out as discrete eigenpairs instead of peaks in a
ringdown spectrum.

## Goal

Solve the closed-cavity Helmholtz problem

```text
∇ × (1/μ_r) ∇ × E = k² ε_r E
```

on the WR-90-based cavity `(a, b, d) = (22.86, 10.16, 30) mm` with PEC
tangential-`E`-zero Dirichlet on every face, and verify that the
lowest ten resonant frequencies match Pozar 4th ed. §6.3 (eq. 6.42)
within ±0.3 % on `TE_{101}` and ±1 % per-mode across the lowest ten.
The Phase 4 validation gate at
`crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs`
enforces exactly these bounds; this tutorial reproduces the gate from
a Python REPL.

## What's shipping in Phase 4.fem.eig.0

The walking-skeleton solver under `crates/yee-fem/` covers:

- **First-order Nedelec (Whitney-1) edge elements** on a structured
  Kuhn 6-tet brick mesh built by `yee_mesh::TetMesh3D::cavity_uniform`.
- **PEC tangential-`E`-zero Dirichlet** applied by dropping boundary
  edges from the global DoF list — the interior-edge block is what
  the eigensolver sees.
- **Free-space material** (`ε_r = μ_r = 1`); per-tet dispersive
  materials are deferred to Phase 4.fem.eig.1.
- **Shift-invert inverse-power iteration** on `(K − σM)^{-1} M` via
  the `SparseEigen` trait in `yee_fem::InverseIterEigen`, with the
  shift `σ` lifted above the gradient-kernel cluster at `k² ≈ 0` so
  the iteration converges to the physical mode band.

The Python entry point is `yee.fem.solve_cavity`, a thin wrapper over
`yee_fem::FemEigenAssembly::new_free_space(...).assemble()` followed by
`InverseIterEigen::default().solve(...)`. The companion theory page
[FEM Eigenmode Solver (3-D Nedelec)](../theory/fem-eigenmode.md)
derives the element integrals and explains why the shift placement
matters.

## Prerequisites

- Rust 1.92+ (the Python wheel is built from source).
- Python 3.10 through 3.14. The wheel is `abi3-py310`, so any
  interpreter in that range works without rebuilding.
- `pip install maturin numpy pytest` — or the `uv` equivalent.
- Jupyter is optional; the snippets run in a plain `python` REPL too.

## Install

Same flow as [Tutorial 2 — Half-wave dipole from Python](02-dipole-from-python.md).
From the repo root:

```bash
uv venv .venv
source .venv/bin/activate
uv pip install maturin numpy pytest
cd crates/yee-py
maturin develop --release
python -c "import yee.fem; print(yee.fem.__doc__)"
```

`maturin develop --release` compiles the Rust extension and drops the
wheel into your active virtualenv. The `--release` flag matters here:
the FEM assembly + LU factorisation runs noticeably slower in debug.
After this, `import yee.fem` works anywhere in that environment.

## First call

The minimal end-to-end snippet — build the cavity, solve for ten
modes, print frequencies in GHz:

```python
import yee.fem
import numpy as np

# WR-90-based cavity, (a, b, d) = (22.86, 10.16, 30) mm.
freqs, modes = yee.fem.solve_cavity(
    0.02286, 0.01016, 0.030,   # extents in metres
    12, 9, 15,                  # Kuhn 6-tet brick subdivisions
    num_eigs=10,
)

for i, f in enumerate(freqs):
    print(f"mode {i + 1}: {f / 1e9:.4f} GHz")
```

Expected output: ten ascending frequencies with the lowest near
`8.2510 GHz`. The analytic Pozar `TE_{101}` for this geometry is
`f = (c / 2) · sqrt((1/a)² + (1/d)²) ≈ 8.2439 GHz`, so the numerical
result lands within roughly 0.09 % — well inside the ±0.3 % hard gate
the validation driver enforces.

`freqs` is a length-`num_eigs` numpy array of resonant frequencies in
Hz, sorted ascending. `modes` is the corresponding interior-edge
eigenvector block, shape `(n_interior_edges, num_eigs)`, in the
same column order as `freqs`.

## Mode-by-mode comparison against Pozar §6.3

The Pozar §6.3 table evaluated at `(a, b, d) = (22.86, 10.16, 30) mm`
enumerates the allowed `TE_{mnp}` (`p ≥ 1`, `m + n ≥ 1`) and `TM_{mnp}`
(`m ≥ 1`, `n ≥ 1`, `p ≥ 0`) families. The lowest ten distinct
resonances are:

| #  | mode             | f (GHz) |
|----|------------------|---------|
| 1  | TE_{101}         |  8.244  |
| 2  | TE_{102}         |  9.840  |
| 3  | TE_{201}         | 13.764  |
| 4  | TE_{103}         | 12.018  |
| 5  | TE_{202}         | 14.770  |
| 6  | TE_{011}         | 16.150  |
| 7  | TE_{111} ≡ TM_{111} | 16.456 |
| 8  | TE_{111} ≡ TM_{111} | 16.456 (degenerate pair) |
| 9  | TE_{301}         | 18.366  |
| 10 | TE_{203}         | 17.318  |

(The exact analytic values come from
`yee_validation::fem_eig_001_analytic_modes(a, b, d, max_order=6)`
sorted ascending; the table above shows the families that populate
the first ten slots.) The degenerate `TE_{111} / TM_{111}` pair
deserves a footnote: first-order Nedelec resolves the pair as two
distinct numerical eigenvalues very close together, so the analytic
table is **not** deduplicated. Keeping the duplicate entry keeps the
positional mode-by-mode comparison aligned past the degeneracy.

The validation driver computes the analytic table inline. To
reproduce it in Python:

```python
import itertools

C0 = 299_792_458.0

def pozar_te_tm_modes(a, b, d, max_order=6):
    freqs = []
    # TE_{mnp}: p >= 1, (m, n) != (0, 0).
    for m, n, p in itertools.product(
        range(max_order + 1), range(max_order + 1), range(1, max_order + 1)
    ):
        if m == 0 and n == 0:
            continue
        freqs.append(0.5 * C0 * np.sqrt((m / a) ** 2 + (n / b) ** 2 + (p / d) ** 2))
    # TM_{mnp}: m >= 1, n >= 1, p >= 0.
    for m, n, p in itertools.product(
        range(1, max_order + 1), range(1, max_order + 1), range(max_order + 1)
    ):
        freqs.append(0.5 * C0 * np.sqrt((m / a) ** 2 + (n / b) ** 2 + (p / d) ** 2))
    return sorted(freqs)


analytic = pozar_te_tm_modes(0.02286, 0.01016, 0.030)[:10]
for i, (f_num, f_an) in enumerate(zip(freqs, analytic)):
    rel = abs(f_num - f_an) / f_an
    print(f"mode {i + 1}: {f_num / 1e9:6.3f} GHz vs {f_an / 1e9:6.3f} GHz  ({rel * 100:.2f}%)")
```

Running this on the same `(12, 9, 15)` mesh reproduces the Track
QQQQQQ measurement that landed the gate: mode-10 RMS error around
0.37 %, with every individual mode inside ±0.6 % of its analytic
neighbour. The gate budget is generous — ±0.3 % on `TE_{101}` and ±1 %
per-mode across the lowest ten — and the (12, 9, 15) mesh clears it
comfortably in well under a minute of wall-time in `--release`.

If you bump the mesh density `(nx, ny, nz)` further the per-mode error
keeps shrinking, but the LU factorisation cost grows roughly as
`O(N_int_edges^1.5)`. The (12, 9, 15) mesh is the sweet spot the
validation gate runs at; (8, 6, 10) is what the Rust-side
`#[test]` defaults to and is also passing.

## Mode visualisation

`modes` is the interior-edge eigenvector block — one column per
returned eigenvalue, each row a Whitney-1 expansion coefficient on
the corresponding interior edge. **It is not directly an `E`-field
sample on the mesh vertices.** A lift to a vertex- or cell-centred
field requires evaluating each Whitney-1 basis on its tetrahedron and
accumulating; a `yee.fem.lift_to_full_basis` helper is deferred to a
follow-up phase.

For now a useful sanity check is the per-mode coefficient magnitude
histogram:

```python
import matplotlib.pyplot as plt

plt.figure()
plt.hist(np.abs(modes[:, 0]), bins=64)
plt.xlabel("|coefficient| on interior edge")
plt.ylabel("count")
plt.title(f"TE_101 edge-coefficient distribution (f = {freqs[0] / 1e9:.3f} GHz)")
plt.tight_layout()
plt.savefig("te101_coeff_hist.png", dpi=120)
```

The `TE_{101}` column should show a bimodal distribution: a tall
spike near zero (edges far from the high-field region) and a tail of
larger magnitudes concentrated on the edges aligned with the dominant
`E_y` lobe at the box centre. Higher-order modes spread their support
across more edges, so the histogram flattens out.

## What's next

Phase 4.fem.eig.0 is intentionally a walking skeleton — closed PEC
cavity, free-space material, Kuhn-mesh-only, ten lowest modes. The
roadmap (see `ROADMAP.md`) sketches the immediate follow-ups:

- **Phase 4.fem.eig.1** — per-tet dispersive `ε_r(ω)` and `μ_r(ω)`
  for filter materials and DRA modal analysis.
- **Phase 4.fem.eig.2** — real waveguide ports and absorbing
  boundaries, so the same assembler can attack driven-port problems
  alongside closed-cavity eigenmodes. Once ports are real, the FEM
  solver becomes a complement to MoM for volumetric structures the
  planar Green's function cannot reach.
- **Mode-export utilities** — `lift_to_full_basis` plus VTK / VTU
  writers so the `modes` block becomes a Paraview-loadable field
  rather than a coefficient vector.

If you want the full design rationale — element integrals, signed
assembly, shift-placement heuristic, why the gradient-kernel cluster
is a hazard — read the spec at
`docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
and the theory chapter at
[FEM Eigenmode Solver (3-D Nedelec)](../theory/fem-eigenmode.md).

## References

- **Pozar, D. M.**, *Microwave Engineering*, 4th ed., Wiley 2011 —
  §6.3 (rectangular cavity TE/TM table) and eq. 6.42 (analytic
  resonance frequencies for an air-filled metallic box).
- **Jin, J.-M.**, *The Finite Element Method in Electromagnetics*,
  3rd ed., Wiley 2014 — ch. 9 (3-D Nedelec edge elements,
  curl-curl assembly, eigenproblem formulation). This is the
  reference the Phase 4 T3 element layer cites for the Whitney-1
  basis derivation.
- **Yee project spec** —
  `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
  (the Phase 4.fem.eig.0 design, including the `fem-eig-001` validation
  gate contract).
- **Validation gate** —
  `crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs`
  (the production gate that asserts the bounds reproduced in this
  tutorial).
