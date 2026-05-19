# Dispersive FEM cavity eigenmode from Python

This tutorial walks through Yee's **Phase 4.fem.eig.1** dispersive
eigensolver from Python. It is the lossy-material follow-on to
[Tutorial 4 — FEM cavity eigenmode from Python](04-fem-cavity-eigenmode.md):
where the v0 walking skeleton returned real resonant frequencies on an
air-filled metallic box, the v1 solver returns **complex eigenfrequencies**
`f = f' − j f''` on the same geometry once you fill it with a
single-pole Drude / Lorentz / Debye dielectric. The imaginary part
encodes loss, the implied Q factor falls out as
`Q = −Re(f) / (2 Im(f))`, and a Newton-Raphson `ω`-tracker drives the
nonlinear `K(ω) e = (ω/c)² M(ω) e` eigenproblem to its fixed point
in a handful of outer iterations.

## Goal

Solve the dispersive closed-cavity Helmholtz problem

```text
∇ × (1/μ_r(ω)) ∇ × E = (ω/c)² ε_r(ω) E
```

on a small rectangular PEC cavity with a frequency-dependent filler,
and recover a complex `f_TE101`. The Phase 4.fem.eig.1 spec
([`docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`](../../superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md))
scopes the v1 deliverable to:

- **Single-pole dispersive `ε_r(ω)`** — Drude, Lorentz, Debye —
  reusing the Phase 2.fdtd.3 ADE `Material` enum verbatim (now lifted
  into `yee-core::material`).
- **Newton-Raphson `ω`-tracker** with analytic Hellmann–Feynman
  derivative; Beyn 2012 / contour-integral nonlinear eigensolve is
  deferred to Phase 4.fem.eig.1.5.
- **Validation gate `fem-eig-002`**: a lossy Drude-loaded
  `(10, 5, 20) mm` cavity must hit TE_{101} within **±0.5 % on
  Re(f)** and **±5 % on Im(f)** against a hand-derived analytic root.

The v1 scope, deferrals, and tolerances are recorded in
[ADR-0039](../decisions/0039-phase-4-fem-eig-1-dispersive-scope.md).

## Prerequisites

- Rust 1.92+ (`rust-toolchain.toml` pins the toolchain; the Python
  wheel is built from source).
- Python 3.10 through 3.14 — the wheel is `abi3-py310`, so any
  interpreter in that range works without rebuilding.
- `uv pip install maturin numpy pytest` — or the equivalent
  `pip install` invocation.

Install matches [Tutorial 4](04-fem-cavity-eigenmode.md) exactly:

```bash
uv venv .venv
source .venv/bin/activate
uv pip install maturin numpy pytest
cargo build --release -p yee-fem -p yee-py
cd crates/yee-py
maturin develop --release
python -c "import yee.fem; help(yee.fem.solve_cavity_dispersive)"
```

`maturin develop --release` is non-negotiable on the dispersive path:
each Newton step re-assembles a complex sparse pencil and factors a
`Complex64` LU, both of which are an order of magnitude slower in
debug than in release.

> **`max_iter` defaults to 8.** This is the spec §9 ceiling and is
> sufficient for the published `fem-eig-002` Drude case (which
> converges in 3–5 outer iterations). For unusually lossy or
> warm-start-far cases, raise the cap via `with_tuning` on the Rust
> side (`DispersiveSolver::with_tuning(db, inner_max_iter, inner_tol)`)
> and pass a larger `max_iter` to `yee.fem.solve_cavity_dispersive`.
> If your run hits the cap without converging, the result dict's
> `"converged"` field will be `False`; iterating beyond ~12 outer steps
> is almost always a Newton-basin problem (spec §11 risk #6) rather
> than a tolerance one.

## Free-space sanity

The Phase 4 plan step D4 includes a regression test
(`assemble_complex_at_real_eps_matches_real_assemble`) that the
complex assembly path **reproduces** the Phase 4.fem.eig.0 real
assembler bit-for-bit when every material is `eps_inf = 1, mu_r = 1`
with no poles. The Python binding inherits this guarantee. The check:

```python
import math
import yee.fem

C0 = 299_792_458.0
A, B, D = 0.022_86, 0.010_16, 0.030  # WR-90-based cavity
NX, NY, NZ = 8, 6, 10

f_te101 = 0.5 * C0 * math.sqrt((1.0 / A) ** 2 + (1.0 / D) ** 2)

materials = [{"tag": 0, "eps_inf": 1.0, "mu_r": 1.0, "poles": []}]
result = yee.fem.solve_cavity_dispersive(
    A, B, D, NX, NY, NZ,
    materials,
    omega_warm_start_hz=0.9 * f_te101,
    max_iter=8,
    tol=1e-6,
)

print(f"converged   : {result['converged']}")
print(f"iterations  : {result['iterations']}")
print(f"f (GHz)     : {result['frequency_hz'] / 1e9:.6f}")
print(f"|Im(f)| (Hz): {abs(result['frequency_hz'].imag):.3e}")
```

Expected output: `converged: True`, `Re(f) ≈ 8.2439 GHz` (matches the
[Tutorial 4](04-fem-cavity-eigenmode.md) walking-skeleton value to
±0.3 %), and `|Im(f)| < 1 MHz`. The composed wavenumber
`k = ω · √(μ₀ ε₀ ε(ω))` returned under `"k_complex"` is real-positive
to working precision. This is the lossless free-space fixed point of
the dispersive solver, and it is the same answer the v0
`yee.fem.solve_cavity` returns. (See
`crates/yee-py/tests/test_fem_dispersive.py::test_free_space_air_matches_solve_cavity`
for the production-test version.)

## Single Drude pole — measurable loss

Filling the same WR-90-based cavity with a single-pole Drude
oscillator shifts Re(f) downward (because `ε_∞ > 1`) and introduces a
non-zero Im(f). The Drude parameters below are tuned for a `tan δ`
that the gate can resolve at ±5 %:

```python
import math
import yee.fem

materials = [
    {
        "tag": 0,
        "eps_inf": 3.78,
        "mu_r": 1.0,
        "poles": [
            {
                "kind": "drude",
                "omega_p": 2.0 * math.pi * 0.4e9,   # plasma frequency
                "gamma":   2.0 * math.pi * 2.0e9,   # collision rate
            },
        ],
    },
]

# Smaller `fem-eig-002` cavity — (10, 5, 20) mm.
r = yee.fem.solve_cavity_dispersive(
    a=0.010, b=0.005, d=0.020,
    nx=8, ny=4, nz=16,
    materials=materials,
    omega_warm_start_hz=8.62e9,
    max_iter=8,
    tol=1e-6,
)

f = r["frequency_hz"]
q = -f.real / (2.0 * f.imag)
print(f"f = {f.real / 1e9:.4f} + ({f.imag / 1e6:+.2f} MHz) j")
print(f"Q = {q:.1f}")
```

Expected output (modulo per-mesh quadrature noise inside the
`fem-eig-002` ±0.5 % / ±5 % envelope):

```text
f ≈ 8.62 + (-9.50 MHz) j GHz
Q ≈ 450
```

The negative `Im(f)` is the physical decay direction (positive-time
convention `e^{-jωt}` with `Im(ω) < 0`); the Q factor follows the
spec §9 closed-form. The warm-start
`omega_warm_start_hz = 8.62e9` is the analytic root computed by hand
in spec §9.1, but the solver also converges from the v0 free-space
air resonance (~16.77 GHz) — see ADR-0039 consequences. The Drude
pole at `ω_p = 2π · 0.4 GHz` sits two decades below the resonance,
which keeps the Newton iteration comfortably inside its
quadratic-convergence basin.

## Three pole models compared

The Phase 4.fem.eig.1 binding accepts any of the three Phase 2.fdtd.3
single-pole models. The pole-spec dictionary keys match the Rust-side
`Material::Drude`, `Material::Lorentz`, `Material::Debye` variants
one-to-one:

| Kind      | Closed form `ε(ω) − ε_∞`              | Required pole-dict keys                     |
|-----------|---------------------------------------|---------------------------------------------|
| `drude`   | `−ω_p² / (ω² − jγω)`                  | `omega_p`, `gamma`                          |
| `lorentz` | `Δε ω₀² / (ω₀² − ω² + 2jδω)`          | `omega_0`, `delta_eps`, `delta`             |
| `debye`   | `Δε / (1 + jωτ)`                      | `delta_eps`, `tau`                          |

A side-by-side comparison fixture (drop into the same WR-90 cavity
geometry as the free-space sanity above, leaving `omega_warm_start_hz`
at the air resonance):

```python
import math

WP   = 2.0 * math.pi * 0.4e9   # rad/s, plasma frequency
GAM  = 2.0 * math.pi * 2.0e9   # rad/s, Drude collision rate
W0   = 2.0 * math.pi * 5.0e9   # rad/s, Lorentz centre
DELTA = 2.0 * math.pi * 0.5e9  # rad/s, Lorentz damping
DE   = 2.0                     # Δε
TAU  = 1.0e-10                 # s, Debye relaxation time

cases = {
    "drude":   [{"kind": "drude",   "omega_p": WP, "gamma": GAM}],
    "lorentz": [{"kind": "lorentz", "omega_0": W0, "delta_eps": DE, "delta": DELTA}],
    "debye":   [{"kind": "debye",   "delta_eps": DE, "tau": TAU}],
}
```

Each entry plugged into a materials list `[{"tag": 0, "eps_inf": 3.78,
"mu_r": 1.0, "poles": <case>}]` and passed to
`yee.fem.solve_cavity_dispersive` returns a complex `frequency_hz`
whose real part lies below the air resonance and whose imaginary part
carries the per-model loss. Drude is the lossy-conductor /
free-electron model; Lorentz is the bound-oscillator dielectric model
with a resonance at `ω₀`; Debye is the relaxation-only polar-liquid
model. Plan step D3 of the Phase 4.fem.eig.1 implementation plan
(`docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`)
contains the closed-form `dε/dω` per variant that the Newton tracker
relies on.

## Validation

The corresponding validation gate, **fem-eig-002**, lives at
`crates/yee-validation/tests/fem_eig_002_lossy_sio2_cavity.rs` once
plan step D6 / Track QQQQQQQQ lands. The gate is the same Drude
fixture above (`a, b, d = 10, 5, 20 mm`, `ε_∞ = 3.78`,
`ω_p = 2π · 0.4 GHz`, `γ = 2π · 2.0 GHz`) and enforces:

1. `|Re(f_FEM) − Re(f_analytic)| / Re(f_analytic) ≤ 0.5 %`
2. `|Im(f_FEM) − Im(f_analytic)| / |Im(f_analytic)| ≤ 5 %`
3. Newton converges in **≤ 8** outer iterations.
4. The bisection fallback (spec §11) does **not** trigger on the
   gate.

The row will land in `crates/yee-fem/validation/README.md` alongside
the v0 `fem-eig-001` entry. The published-reference rationale — Pozar
§3.1 closed-form complex propagation constant, Bucur et al. (1996)
fused-silica permittivity at ~10 GHz, deliberately exaggerated Drude
loss for measurable Im(f) — is in spec §9.

## What's next

Phase 4.fem.eig.1 is a strict extension of the v0 walking skeleton:
the Phase 4.fem.eig.0 `fem-eig-001` gate stays green unmodified
(verified by D1's `complex_matches_real_for_pure_real_eps_mu` test
and D4's `assemble_complex_at_real_eps_matches_real_assemble`
regression). The roadmap lays out the immediate follow-ons:

- **Phase 4.fem.eig.1.1** — multi-pole expansions (Drude + multiple
  Lorentz oscillators) for materials whose published model needs more
  than one pole.
- **Phase 4.fem.eig.1.2** — magnetic dispersion `μ_r(ω)` for ferrite
  cases.
- **Phase 4.fem.eig.1.5** — Beyn 2012 / Sakurai–Sugiura
  contour-integral nonlinear eigensolve, reserved for the case where
  Newton-with-bisection-fallback proves insufficient.
- **Phase 4.fem.eig.3** — dielectric-resonator-antenna validation
  (Petosa Ch. 3) with the puck modelled as a single-pole Drude.

## References

- **Pozar, D. M.**, *Microwave Engineering*, 4th ed., Wiley 2012 —
  §3.1 (lossy-waveguide propagation constants, closed-form complex
  `γ(ω)`) and §6.3 (rectangular cavity resonator with material
  loss; the analytic baseline that fem-eig-002 extends to complex ε).
- **Jin, J.-M.**, *The Finite Element Method in Electromagnetics*,
  3rd ed., Wiley 2014 — §9.5 (lossy-material FEM eigenvalue problems;
  Hellmann–Feynman differentiation for nonlinear eigenproblems —
  the analytic derivative the Newton tracker uses).
- **ADR-0039** —
  [`docs/src/decisions/0039-phase-4-fem-eig-1-dispersive-scope.md`](../decisions/0039-phase-4-fem-eig-1-dispersive-scope.md):
  v1 scope (Newton-only, Beyn 2012 deferred), Material relocation
  to `yee-core`, complex-symmetric inner-product convention, and the
  fem-eig-002 ±0.5 % / ±5 % gate envelope.
- **Yee project spec** —
  [`docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`](../../superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md).
- **Yee project plan** —
  [`docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`](../../superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md);
  D3 (closed-form `dε/dω`), D5 (Newton tracker pseudocode), D6
  (fem-eig-002 gate construction).
- **Phase 2.fdtd.3** — `crates/yee-fdtd/src/material.rs`: the
  single-pole Drude / Lorentz / Debye `Material` enum reused
  verbatim on the FEM side.
- **Phase 4.fem.eig.0** —
  [Tutorial 4 — FEM cavity eigenmode from Python](04-fem-cavity-eigenmode.md);
  the lossless walking skeleton this tutorial extends.
- **Validation driver** —
  `crates/yee-validation/tests/fem_eig_002_lossy_sio2_cavity.rs`
  (lands with Track QQQQQQQQ).
