# FDTD lumped series-LC resonance from Python

This tutorial validates the FDTD lumped-element series-RLC port by checking
that the simulated ring-down frequency matches the analytic LC resonance
`f₀ = 1/(2π√LC)`.

## Physics

A lumped series-LC circuit (R small, underdamped) has a natural resonant
frequency given by circuit theory (Pozar §2.4):

```
f₀ = 1 / (2π √(LC))
```

With L = 1 nH and C ≈ 25.33 pF, the analytic value is f₀ = 1 GHz exactly.

The FDTD simulator integrates the series-RLC ODE alongside the Yee field
update (Taflove & Hagness §15.10). The gate validates this path for the
series-RLC case — the pure-resistor path was validated in Phase 2.fdtd.6.

## Running the gate

```python
from yee import run_lc_resonance

result = run_lc_resonance()
print(result)
# LcResonanceResult(f_measured_hz=1.0000e+09, f_analytic_hz=1.0000e+09, rel_err=..., passed=True)

assert result.passed
print(f"f_measured = {result.f_measured_hz:.4e} Hz")
print(f"f_analytic = {result.f_analytic_hz:.4e} Hz")
print(f"rel_err    = {result.rel_err * 100:.3f} %")
```

## Gate tolerance

| Field | Value |
|---|---|
| L | 1 nH |
| C | ≈25.33 pF |
| R | 1 Ω |
| Q | ≈6.28 |
| f₀ analytic | 1 GHz |
| Gate | `abs(f_measured − f₀) / f₀ < 2 %` |

## Result fields

- `f_measured_hz` — peak DFT frequency from the inductor-current ring-down
- `f_analytic_hz` — analytic 1/(2π√LC) = 1 GHz
- `rel_err` — `|f_measured − f_analytic| / f_analytic`
- `passed` — `True` iff `rel_err < 2 %`

## References

- Pozar, *Microwave Engineering*, 4th ed., §2.4
- Taflove & Hagness, *Computational Electrodynamics*, 3rd ed., §15.10
- ADR-0080 (Phase 2.fdtd.6.1 design decision)
