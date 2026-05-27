# FDTD Dipole Radiation Pattern from Python

This tutorial runs the FDTD short-dipole radiation-pattern validation gate
(fdtd-203) from Python using `yee.run_dipole_pattern`.

## Background

A z-polarised short dipole radiates in a doughnut pattern described by
Balanis *Antenna Theory* §4.2:

```
E_θ ∝ sin θ
```

Maximum radiation is broadside (θ = 90°); nulls occur at endfire
(θ = 0° and 180°). `run_dipole_pattern` drives a sinusoidal J_z
current source through a 5-cell FDTD dipole at 1 GHz on a 60³ grid,
accumulates the NTFF surface currents for 800 steps, and sweeps the
far field over θ ∈ [0°, 180°] in 5° steps.

## Quick start

```python
from yee import run_dipole_pattern

# Smoke check: n_steps=5 only verifies the API, not the physics.
r = run_dipole_pattern(n_steps=5)
print(r)  # DipolePatternResult(...)
print("theta_deg shape:", r.theta_deg_array().shape)  # (37,)
```

## Full gate run (~30 s in release mode)

```python
from yee import run_dipole_pattern

r = run_dipole_pattern()  # 60³ grid, 800 steps, 1 GHz
print(f"θ=  0°: {r.e_theta_0:.4f}  (expected ~0, endfire null)")
print(f"θ= 45°: {r.e_theta_45:.4f}  (expected ~0.707 = sin 45°)")
print(f"θ= 90°: {r.e_theta_90:.4f}  (expected ~1.0, broadside peak)")
print(f"θ=135°: {r.e_theta_135:.4f}  (expected ~0.707 = sin 135°)")
print(f"θ=180°: {r.e_theta_180:.4f}  (expected ~0, endfire null)")
assert r.passed, f"fdtd-203 gate failed: {r}"
print("fdtd-203: PASS")
```

## Plotting the pattern

```python
import numpy as np
import matplotlib.pyplot as plt
from yee import run_dipole_pattern

r = run_dipole_pattern()
theta = r.theta_deg_array()
e_fdtd = r.e_theta_array()
e_analytic = np.sin(np.deg2rad(theta))
e_analytic /= e_analytic.max()

fig, ax = plt.subplots()
ax.plot(theta, e_fdtd, label="FDTD (60³, 800 steps)")
ax.plot(theta, e_analytic, "--", label="Analytic sin θ")
ax.set_xlabel("θ (degrees)")
ax.set_ylabel("|E_θ| (normalized)")
ax.set_title("Short-Dipole Radiation Pattern")
ax.legend()
plt.show()
```

## `DipolePatternResult` fields

| Field | Description |
|-------|-------------|
| `e_theta_0` | \|E_θ\|(0°) — endfire null |
| `e_theta_45` | \|E_θ\|(45°) ≈ sin 45° ≈ 0.707 |
| `e_theta_90` | \|E_θ\|(90°) — broadside peak, normalized to 1.0 |
| `e_theta_135` | \|E_θ\|(135°) ≈ sin 135° ≈ 0.707 |
| `e_theta_180` | \|E_θ\|(180°) — endfire null |
| `passed` | `True` iff all gate criteria hold |
| `theta_deg_array()` | 37-point θ array (numpy float64) |
| `e_theta_array()` | 37-point \|E_θ\| array (numpy float64), max=1.0 |

## Connection to the Rust gate

The underlying Rust integration test is
`crates/yee-fdtd/tests/dipole_pattern.rs` (Phase 2.fdtd.4).
Run it directly with:

```bash
cargo test -p yee-fdtd --test dipole_pattern --release -- --include-ignored
```

## See also

- [FDTD: CPML, NTFF, TF/SF, lumped, dispersive](../theory/fdtd-details.md)
- [FDTD cavity resonance from Python](11-fdtd-cavity-resonance-from-python.md)
- [FDTD lossy cavity from Python](10-fdtd-lossy-cavity-from-python.md)
