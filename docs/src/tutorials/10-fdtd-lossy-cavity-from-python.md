# FDTD Lossy Cavity from Python

This tutorial shows how to use `yee.run_cavity_q()` to simulate a
lossy rectangular PEC cavity, extract the Q-factor from the TE₁₀₁
ring-down, and verify the result against the analytic formula.

## Why Q-factor matters

The quality factor Q of a resonator measures how many radians of
oscillation it takes to dissipate the stored energy:

```text
Q = ω₀ · W / P_loss
```

High-Q resonators are used in microwave filters (low insertion loss
within the passband), oscillators (low phase noise), and cavity
perturbation measurements (non-destructive material characterisation).
A first-principles FDTD simulation can predict Q before any hardware
is built, provided the E-update kernel correctly models the ohmic loss
σ · E.

## Physics background

A rectangular PEC cavity with dimensions `a × b × d` filled uniformly
with a medium of electric conductivity σ supports the TE₁₀₁ dominant
mode at resonant frequency (Pozar §6.3):

```text
f₁₀₁ = (c / 2) · √((1/a)² + (1/d)²)
```

When the medium has conductivity σ the electric field decays
exponentially after the source is removed (Taflove §3.7):

```text
E_y(t) ∝ exp(−t / τ),   τ = 2ε₀ / σ
```

The Q-factor follows directly:

```text
Q = ω₁₀₁ · ε₀ / σ = 2π · f₁₀₁ · ε₀ / σ
```

For σ₀ = 2.96 × 10⁻³ S/m and `a = d = 0.20 m` (`f₁₀₁ ≈ 1.06 GHz`):

```text
Q_analytic = 2π × 1.06 × 10⁹ × 8.85 × 10⁻¹² / 2.96 × 10⁻³ ≈ 20
```

The FDTD solver implements the CA/CB E-update (Taflove §3.7 eq. 3.63):

```text
E_y^{n+1} = CA · E_y^n + CB · (curl H)^{n+1/2}

CA = (2ε₀ − σ Δt) / (2ε₀ + σ Δt)
CB = 2Δt / ((2ε₀ + σ Δt) · dx)
```

For σ₀ and Δt ≈ 17.3 ps the stability coefficient CA ≈ 0.9971 < 1,
confirming the update is energy-dissipating and stable.

## Running the simulation

The `run_cavity_q()` function encapsulates the full simulation pipeline:
build a lossy vacuum grid, inject a Gaussian pulse, let the TE₁₀₁ mode
ring down, fit the log-linear decay, and return the extracted Q.

```python
from yee import run_cavity_q

# Default: 20×10×20 grid, dx=10 mm, σ=2.96e-3 S/m → Q≈20
result = run_cavity_q()

print(f"f₁₀₁       = {result.f101_hz * 1e-9:.4f} GHz")
print(f"Q_analytic = {result.q_analytic:.4f}")
print(f"Q_measured = {result.q_measured:.4f}")
print(f"rel_err    = {result.rel_err:.4e}")
print(f"passed     = {result.passed}")
```

Expected output (tolerances ±5 %):

```
f₁₀₁       = 1.0607 GHz
Q_analytic = 20.0001
Q_measured = 19.9924
rel_err    = 3.8e-04
passed     = True
```

## Interpreting the result

| Field | Meaning |
|---|---|
| `f101_hz` | Analytic TE₁₀₁ frequency from Pozar §6.3 |
| `q_analytic` | `2π · f₁₀₁ · ε₀ / σ` |
| `q_measured` | Extracted from log-linear ring-down fit on the last 2/3 of `n_ring` steps |
| `rel_err` | `|q_measured − q_analytic| / q_analytic` |
| `passed` | `True` iff `rel_err < 0.05` (fdtd-202 gate: ±5 %) |

## Visualising the ring-down

The probe time series (E_y at a strong TE₁₀₁ antinode) is available
via `result.probe_array()`:

```python
import numpy as np
import matplotlib.pyplot as plt

from yee import run_cavity_q

result = run_cavity_q(n_ring=6000)
probe = result.probe_array()

# Grid time step at dx=10 mm, safety factor 0.9
dt = 0.9 / (299_792_458 * (3 / 0.01**2) ** 0.5)
t_ns = np.arange(len(probe)) * dt * 1e9  # nanoseconds

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# Linear scale — shows the decaying sinusoid
ax1.plot(t_ns, probe)
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("E_y (arb.)")
ax1.set_title("TE₁₀₁ ring-down — linear scale")

# Log scale — straight line confirms single-mode exponential decay
ax2.semilogy(t_ns, np.abs(probe).clip(1e-30))
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("|E_y|")
ax2.set_title("TE₁₀₁ ring-down — log scale (slope = −1/τ)")

plt.tight_layout()
plt.savefig("cavity_q_ringdown.png", dpi=150)
print("Saved cavity_q_ringdown.png")
```

The log-scale plot should show a straight-line decay after the initial
multi-mode transient (~1-2 ns), confirming that the TE₁₀₁ mode
dominates the ring-down. The slope of the line equals −1/τ, from which
Q = π · f₁₀₁ · τ.

## Varying the conductivity

Q scales inversely with σ:

```python
from yee import run_cavity_q

for sigma in [2.96e-2, 2.96e-3, 2.96e-4]:
    r = run_cavity_q(sigma=sigma, n_ring=6000 if sigma > 1e-3 else 60_000)
    print(f"σ = {sigma:.2e} S/m  →  Q_analytic = {r.q_analytic:.1f},  "
          f"Q_measured = {r.q_measured:.1f},  passed = {r.passed}")
```

Note: for very high Q (low σ) the ring-down takes many more time steps
to decay sufficiently for a reliable fit — increase `n_ring`
proportionally (Q ≈ 200 needs `n_ring ≈ 60_000`).

## Next steps

- **General material builder:** a follow-on `PyYeeGrid` API (Phase
  2.fdtd.py.1) will expose per-cell ε_r and σ from Python, enabling
  inhomogeneous cavities.
- **Frequency domain:** the ring-down DFT scan used by fdtd-201 can be
  combined with the Q extraction here to produce a full S-parameter
  response for a driven resonator.
- **Validation:** the same physics is exercised by `yee validate fdtd-202`
  (or `yee validate all`) — run it to confirm the gate passes on any
  machine.
