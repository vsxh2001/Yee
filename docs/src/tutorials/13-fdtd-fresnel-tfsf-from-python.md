# FDTD TF/SF Fresnel Transmission from Python

This tutorial runs the FDTD TF/SF Fresnel-transmission validation gate
(fdtd-204) from Python using `yee.run_fresnel_tfsf`.

## Background

### Total-Field / Scattered-Field (TF/SF)

`PlaneWaveSource` implements the TF/SF technique (Taflove §5.7): a continuous
sinusoidal plane wave is injected across a rectangular surface separating the
**total-field** region (incident + scattered) from the **scattered-field**
region (scattered only). CPML boundaries on all six faces absorb the outgoing
fields.

### Fresnel transmission through a dielectric slab

When a plane wave hits a lossless dielectric slab of thickness *d* and
permittivity ε_r at normal incidence, the amplitude transmission coefficient
is given by the Born & Wolf §1.6.2 transfer-matrix formula:

```
n₂ = √ε_r

δ  = 2π · f · n₂ · d / c          (electrical phase depth in slab)

r₁₂ = (1 − n₂) / (1 + n₂)         (vacuum → slab reflection)
r₂₃ = (n₂ − 1) / (n₂ + 1)         (slab → vacuum reflection)
t₁₂ = 2 / (1 + n₂)
t₂₃ = 2·n₂ / (n₂ + 1)

t_slab = t₁₂ · t₂₃ · e^{jδ} / (1 + r₁₂ · r₂₃ · e^{j2δ})

|t_analytic| = |t_slab|
```

For the default scenario (ε_r = 2.2, d = 5 mm, f = 10 GHz):

```
n₂ ≈ 1.4832,  δ ≈ 1.553 rad  (near quarter-wave)
|t_analytic| ≈ 0.927
```

The slab is near-quarter-wave at 10 GHz, giving a sub-unity transmission
(about 27% field-amplitude reduction) that is well measurable while avoiding
the deep transmission minimum of a half-wave transformer.

## Quick start

```python
from yee import run_fresnel_tfsf

r = run_fresnel_tfsf()          # 80³ grid, 600 steps, 10 GHz
print(r)
# FresnelTfsfResult(t_measured=..., t_analytic=0.9270, rel_err=..., passed=True)
assert r.passed, f"fdtd-204 gate failed: rel_err={r.rel_err:.3f}"
print("fdtd-204: PASS")
```

## Geometry

| Parameter            | Value                                                      |
|---------------------|------------------------------------------------------------|
| Grid                | 80 × 80 × 80 cells, dx = 1 mm                              |
| CPML                | npml = 10, all six faces                                   |
| TF box              | i₀=12, i₁=69, j₀=1, j₁=78, k₀=1, k₁=78                  |
| Dielectric slab     | i ∈ [50, 55), ε_r = 2.2; elsewhere ε_r = 1.0              |
| TF/SF frequency     | 10 GHz (CW sinusoid, ramp_steps = 50)                      |
| Probe (incident)    | (25, 40, 40) — vacuum, inside TF box, before slab          |
| Probe (transmitted) | (62, 40, 40) — vacuum, inside TF box, after slab           |
| n_steps             | 600                                                        |
| Settling window     | steps 0–199 discarded; measurement over [200, 600)         |

### Measurement protocol

Two independent simulations are run with identical `PlaneWaveSource` setups:

1. **Vacuum run** (no slab): records E_z at `probe_inc = (25, 40, 40)` for
   the incident amplitude baseline.
2. **Slab run** (`with_eps_r_cells` set): records E_z at
   `probe_trans = (62, 40, 40)` for the transmitted amplitude.

```
A_inc   = max |E_z_vacuum(probe_inc)|  over steps [200, 600)
A_trans = max |E_z_slab(probe_trans)| over steps [200, 600)
t_measured = A_trans / A_inc
```

## Gate criteria (fdtd-204)

```
|t_measured / t_analytic − 1| < 0.05   (5%)
```

The `passed` attribute of `FresnelTfsfResult` encodes this check.

## Full gate run (~5–15 min in release mode)

```python
from yee import run_fresnel_tfsf

r = run_fresnel_tfsf()   # default: 80³, 600 steps
print(f"t_measured  = {r.t_measured:.4f}")
print(f"t_analytic  = {r.t_analytic:.4f}  (Born-Wolf §1.6.2, ≈0.927)")
print(f"rel_err     = {r.rel_err:.4f}  (gate < 0.05)")
print(f"passed      = {r.passed}")
```

## Customising parameters

All parameters are keyword-only:

```python
from yee import run_fresnel_tfsf

# Quick API smoke check (5 steps — not a physics gate).
r = run_fresnel_tfsf(n_steps=5)
print("t_analytic:", r.t_analytic)   # analytic value always valid

# Thicker slab: 10 cells × 1 mm = 10 mm.
r2 = run_fresnel_tfsf(slab_i0=50, slab_i1=60, n_steps=600)
print(r2)

# Different permittivity (e.g. PTFE, ε_r ≈ 2.1).
r3 = run_fresnel_tfsf(eps_r=2.1, n_steps=600)
print(r3)
```

| Parameter  | Default | Description                                     |
|------------|---------|-------------------------------------------------|
| `nx`       | `80`    | Grid cells in the x-direction (propagation)     |
| `ny`       | `80`    | Grid cells in the y-direction (transverse)      |
| `nz`       | `80`    | Grid cells in the z-direction (transverse)      |
| `dx`       | `1e-3`  | Cell size in metres (1 mm)                      |
| `eps_r`    | `2.2`   | Slab relative permittivity                      |
| `slab_i0`  | `50`    | Slab start index (inclusive)                    |
| `slab_i1`  | `55`    | Slab end index (exclusive) — 5-cell slab        |
| `freq_hz`  | `10e9`  | TF/SF frequency in Hz (10 GHz)                  |
| `n_steps`  | `600`   | Number of FDTD time steps                       |
| `settle`   | `200`   | Steps discarded before measuring peak amplitude |

## `FresnelTfsfResult` fields

| Field        | Description                                                  |
|--------------|--------------------------------------------------------------|
| `t_measured` | Measured amplitude transmission coefficient A_trans / A_inc  |
| `t_analytic` | Analytic prediction from Born & Wolf §1.6.2 transfer-matrix  |
| `rel_err`    | \|t_measured − t_analytic\| / t_analytic                    |
| `passed`     | `True` iff `rel_err < 0.05` (fdtd-204 gate ±5%)             |

## Connection to the Rust gate

The underlying Rust unit tests are in `crates/yee-validation/src/lib.rs`:

```bash
# Fast analytic smoke (< 1 ms, no FDTD):
cargo test -p yee-validation --lib -- fdtd_204 --release

# Full physics gate (~5-15 min in release):
cargo test -p yee-validation -- fdtd_204_live_gate --release --ignored
```

## See also

- [FDTD: CPML, NTFF, TF/SF, lumped, dispersive](../theory/fdtd-details.md)
- [FDTD dipole radiation pattern from Python](12-fdtd-dipole-pattern-from-python.md)
- [FDTD cavity resonance from Python](11-fdtd-cavity-resonance-from-python.md)
- [FDTD lossy cavity from Python](10-fdtd-lossy-cavity-from-python.md)
