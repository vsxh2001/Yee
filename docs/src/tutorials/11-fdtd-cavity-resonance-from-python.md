# FDTD Cavity Resonance from Python

This tutorial shows how to use `yee.run_cavity_resonance()` to simulate a
rectangular PEC cavity, extract the TE₁₀₁ resonant frequency via a
broadband DFT scan, and verify the result against the analytic formula.

## Why cavity resonance matters

Resonant frequency accuracy is a fundamental test of any time-domain
electromagnetic solver. A systematic frequency shift compared to the
analytic formula is the direct fingerprint of **numerical dispersion** —
the artefact by which the discretised Yee grid propagates plane waves
at a speed slightly below *c*. Measuring this shift on a simple PEC
cavity with a known closed-form answer lets you:

- Quantify the grid-dispersion error before running a more complex
  design (the error scales as `O(Δx²)`, so halving the cell size
  quarter-reduces the error);
- Confirm that the PEC boundary implementation correctly places an
  anti-node at the cavity walls;
- Validate that the DFT post-processing pipeline resolves peaks
  correctly at a target frequency resolution.

## Physics background

A rectangular PEC cavity with dimensions `a × b × d` (width, height,
depth) supports the dominant TE₁₀₁ mode at resonant frequency
(Pozar §6.3):

```text
f₁₀₁ = (c/2) · √((1/a)² + (1/d)²)
```

The mode has no dependence on the `b` (height) dimension and a single
half-wavelength along both `a` and `d`. For a square cross-section
`a = d = 0.20 m` and `b = 0.10 m` (`dx = 10 mm`, 20×10×20 cells):

```text
f₁₀₁ = (299 792 458 / 2) · √((1/0.20)² + (1/0.20)²) ≈ 1.0607 GHz
```

The FDTD measurement works by injecting a broadband Gaussian pulse at an
off-centre source cell and recording the E_y probe at an opposite
off-centre cell. The pulse bandwidth is wide enough to excite the TE₁₀₁
mode and several higher modes. After `n_steps` time steps the DFT of
the probe series is scanned over `n_freq_bins` candidate frequencies in
`[0.65·f_ref, 1.50·f_ref]`; the candidate with the highest DFT magnitude
is reported as the extracted frequency.

On a 20×10×20 grid with `dx = 10 mm` the expected numerical dispersion
is approximately 0.5–1%, well within the ±2.5% gate tolerance.

## Running the simulation

The `run_cavity_resonance()` function encapsulates the full pipeline:
build a vacuum grid, inject a Gaussian pulse, apply PEC boundaries, run
the time loop, scan the DFT, and return the extracted frequency.

```python
from yee import run_cavity_resonance

# Default: 20×10×20 grid, dx=10 mm, 30 000 steps
result = run_cavity_resonance()

print(f"f_analytic  = {result.f_analytic_hz * 1e-9:.6f} GHz")
print(f"f_extracted = {result.f_extracted_hz * 1e-9:.6f} GHz")
print(f"rel_err     = {result.rel_err:.4e}")
print(f"passed      = {result.passed}")
```

Expected output (exact values depend on floating-point constants and
grid dispersion; typical run):

```
f_analytic  = 1.060660 GHz
f_extracted = 1.061140 GHz
rel_err     = 4.5e-04
passed      = True
```

## Interpreting the result

| Field | Meaning |
|---|---|
| `f_analytic_hz` | Analytic TE₁₀₁ frequency from Pozar §6.3 |
| `f_extracted_hz` | Peak frequency from the DFT magnitude scan |
| `rel_err` | `|f_extracted_hz − f_analytic_hz| / f_analytic_hz` |
| `passed` | `True` iff `rel_err < 0.025` (fdtd-201 gate: ±2.5 %) |

The typical relative error on the default 20×10×20 grid is under 1%,
which is consistent with second-order numerical dispersion at
`Δx = 10 mm` (roughly 14 cells per wavelength at 2.1 GHz).

## Plotting the probe time series

The full probe time series is available via `result.probe_array()`:

```python
import numpy as np
import matplotlib.pyplot as plt

from yee import run_cavity_resonance

result = run_cavity_resonance(n_steps=30_000)
probe = result.probe_array()

# Grid time step at dx=10 mm, Courant safety factor 0.9
dx = 0.01  # metres
c = 299_792_458.0
dt = 0.9 / (c * (3.0 / dx**2) ** 0.5)  # seconds
t_ns = np.arange(len(probe)) * dt * 1e9  # nanoseconds

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# Time-domain probe signal
ax1.plot(t_ns, probe)
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("E_y (arb.)")
ax1.set_title("PEC cavity probe — time domain")

# DFT magnitude (brute-force scan over the same range used by the driver)
f_ref = result.f_analytic_hz
freqs = np.linspace(0.65 * f_ref, 1.50 * f_ref, 400)
dft_mag = np.array([
    abs(sum(probe[n] * np.exp(-2j * np.pi * f * n * dt)
            for n in range(len(probe))))
    for f in freqs
])
ax2.plot(freqs * 1e-9, dft_mag)
ax2.axvline(result.f_extracted_hz * 1e-9, color="r", linestyle="--",
            label=f"extracted {result.f_extracted_hz * 1e-9:.4f} GHz")
ax2.axvline(f_ref * 1e-9, color="g", linestyle=":",
            label=f"analytic {f_ref * 1e-9:.4f} GHz")
ax2.set_xlabel("Frequency (GHz)")
ax2.set_ylabel("DFT magnitude (arb.)")
ax2.set_title("TE₁₀₁ DFT scan")
ax2.legend()

plt.tight_layout()
plt.savefig("cavity_resonance.png", dpi=150)
print("Saved cavity_resonance.png")
```

The DFT magnitude plot should show a sharp peak very close to the
analytic TE₁₀₁ frequency, with a small shift attributable to numerical
dispersion. The time-domain signal shows a multi-mode transient during
the first ~5 ns followed by a quasi-sinusoidal steady-state as the
higher modes damp out via the PEC boundaries (all energy is trapped, so
the envelope stays flat in the lossless case).

## Gate and validation

The gate tolerance is **±2.5%**, set conservatively to accommodate
coarse grids in CI. The underlying FDTD Rust integration test (in
`crates/yee-fdtd/tests/cavity_resonance.rs`) is `#[ignore]`-gated
because it runs 30 000 steps on a 20×10×20 grid (~5–15 s in release
mode). Run it directly with:

```bash
cargo test -p yee-fdtd --test cavity_resonance --release -- --ignored --nocapture
```

The same physics case is registered in the validation aggregator as
`fdtd-201` (skipped by default to keep `yee validate all` fast):

```bash
yee validate all   # shows fdtd-201 as Skipped
```

A tighter **±0.5% refinement** path is available by increasing the grid
resolution to `dx = 5 mm` (40×20×40 cells) and `n_steps = 60_000`:

```python
result = run_cavity_resonance(nx=40, ny=20, nz=40, dx=0.005, n_steps=60_000)
print(f"rel_err at dx=5mm: {result.rel_err:.4e}")
```

At half the cell size the dispersion error quarters, giving
`rel_err < 0.2%` in typical runs.

## Next steps

- **Lossy cavity Q-factor:** combine the resonance frequency extraction
  here with the ring-down fit from Tutorial 10 to measure Q from a
  single broadband simulation.
- **Higher-order modes:** the `fdtd-201-x` gate (registered as Skipped)
  targets the TE₂₀₁ mode on a 24×4×16 grid — see ADR-0066.
- **CPML open-boundary variant:** replace the PEC walls with a CPML
  absorber and drive the cavity from a wave-port; see the CPML
  reflection gate (`yee validate fdtd-cpml`).
