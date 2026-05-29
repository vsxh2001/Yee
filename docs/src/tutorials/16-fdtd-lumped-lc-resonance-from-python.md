# Tutorial 16: FDTD lumped series-LC resonance from Python

This tutorial demonstrates the `run_lc_resonance()` function, which runs the
**fdtd-206** validation gate: a lumped series-LC port oscillating at its
natural resonant frequency in a small PEC box.

## Physics

A series R-L-C circuit (R small, underdamped) has a natural resonant frequency

```
f₀ = 1 / (2π √(LC))
```

(Pozar §2.4; Hayt & Kemmerly §14.1). With L = 1 nH and
C = 1/(4π² f₀² L) ≈ 25.330 pF, the analytic f₀ = 1 GHz exactly.

The FDTD lumped-port model integrates the series-RLC ODE alongside the Yee
leapfrog update (Taflove & Hagness §15.10). The gate verifies the ODE
integration accurately reproduces the analytic resonant frequency.

## Running from Python

```python
from yee import run_lc_resonance

result = run_lc_resonance()
print(result)
# LcResonanceResult(f_measured_hz=1.0000e+09, f_analytic_hz=1.0000e+09, rel_err=<1e-3, passed=True)

assert result.passed, f"Gate failed: rel_err={result.rel_err:.4e}"
print(f"f_measured = {result.f_measured_hz/1e9:.4f} GHz")
print(f"f_analytic = {result.f_analytic_hz/1e9:.4f} GHz")
print(f"rel_err    = {result.rel_err*100:.3f}%")
```

## Geometry and parameters

| Parameter | Value |
|-----------|-------|
| Grid | 5 × 5 × 40 cells, dx = 1 mm |
| Port cell | (2, 2, 20) |
| Inductance L | 1 nH |
| Capacitance C | ≈ 25.330 pF |
| Resistance R | 1 Ω |
| Q factor | ≈ 6.28 |
| Analytic f₀ | 1.000 GHz |
| Gate tolerance | ±2 % |

## Simulation procedure

1. **Kick phase** (30 steps): a narrow Gaussian pulse on E_z excites the LC
   with a broadband impulse.
2. **Ring-down phase** (5 000 steps): the LC port oscillates at f₀ while
   the resistor R damps the envelope (τ ≈ Q/πf₀ ≈ 2 ns).
3. **DFT scan**: 1 000 bins from 0.5 GHz to 1.5 GHz (2 MHz spacing) find
   the peak of the inductor-current spectrum.
4. **Gate**: |f_measured − f₀| / f₀ < 2 %.

## Result fields

| Field | Description |
|-------|-------------|
| `f_measured_hz` | DFT peak frequency (Hz) |
| `f_analytic_hz` | Analytic f₀ = 1/(2π√LC) (Hz) |
| `rel_err` | |f_measured − f₀| / f₀ |
| `passed` | `True` iff rel_err < 2 % |

## References

- Pozar, "Microwave Engineering," 4th ed., §2.4 (series resonator).
- Hayt & Kemmerly, "Engineering Circuit Analysis," §14.1.
- Taflove & Hagness, "Computational Electrodynamics," 3rd ed., §15.10.
- ADR-0080: `docs/src/decisions/0080-fdtd-206-lc-resonance.md`
