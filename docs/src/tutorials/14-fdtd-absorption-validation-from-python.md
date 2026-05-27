# FDTD Absorption Validation from Python

This tutorial runs three FDTD physics-validation gates — CPML boundary
absorption (cpml-001), near-to-far-field pattern accuracy (ntff-001), and
dispersive Drude-material reflection (dispersive-001) — from Python using the
`yee` package. Together they verify that the Yee FDTD core handles
wave-absorbing boundaries, far-field projection, and frequency-dispersive
materials correctly before any application-level simulation is attempted.

## 1  CPML Reflection Gate (`run_cpml_reflection`)

The Convolutional Perfectly-Matched Layer (CPML) is the primary absorbing
boundary used by all FDTD simulations in Yee. To quantify how well it
suppresses back-reflections, `run_cpml_reflection` runs two 50³-cell vacuum
simulations (dx = 1 mm, 300 steps): one with PEC (hard) boundaries and one
with 10-cell CPML boundaries, both driven by the same Gaussian pulse at the
centre. The peak reflected amplitude at a probe near the boundary is measured
in each run, and the reduction in dB is reported.

```python
from yee import run_cpml_reflection

r = run_cpml_reflection()
print(r.reduction_db, r.passed)
# e.g. 69.3 True
print(r)
# CpmlReflectionResult(reduction_db=69.30, passed=True)
assert r.passed, f"cpml-001 gate failed: {r.reduction_db:.2f} dB < 30 dB"
```

Gate criterion (cpml-001, Roden–Gedney 2000): `reduction_db ≥ 30 dB`. In
practice the CPML achieves ≥ 69 dB on this grid, so the gate is conservatively
set at 30 dB to allow for grid-parameter variations.

## 2  NTFF Broadside/Endfire Gate (`run_ntff_broadside`)

The near-to-far-field (NTFF) transform projects surface-current equivalents
from a Huygens box just inside the CPML to the far field. For a z-directed
point dipole, the far-field pattern is `E_θ ∝ sin θ` (Balanis §4.2), giving a
maximum at broadside (θ = 90°) and a null at endfire (θ = 0°). The
`run_ntff_broadside` function builds a 50³ vacuum grid (dx = 1 mm, 10-cell
CPML), drives an E_z source at the centre, runs 2000 steps at 15 GHz, and
evaluates the NTFF far field at both angles. A large broadside-to-endfire ratio
(in dB) confirms that the NTFF projection is working correctly.

```python
from yee import run_ntff_broadside

r = run_ntff_broadside()
print(r.ratio_db, r.passed)
# e.g. 38.5 True
print(r)
# NtffResult(ratio_db=38.50, passed=True)
assert r.passed, f"ntff-001 gate failed: {r.ratio_db:.2f} dB < 20 dB"
```

Gate criterion (ntff-001, Balanis §4.2): `ratio_db ≥ 20 dB`. A perfect null
at endfire would give infinite dB; the 20 dB threshold accommodates the finite
grid resolution that fills the null partially.

## 3  Dispersive Drude Material Gate (`run_dispersive_drude`)

Frequency-dispersive materials are handled by the Auxiliary Differential
Equation (ADE) method. The Drude model,

```
ε_r(ω) = ε_∞ − ω_p² / (ω² + jγω)
```

describes metals and electron-plasma materials. `run_dispersive_drude` builds
an 80³ vacuum grid (dx = 1 mm, 10-cell CPML), fills cells i ∈ [50, 70) with a
Drude material (ε_∞ = 1, ω_p = 2π·20 GHz, γ = 2π·5 GHz), injects a broadband
Gaussian E_z pulse, and measures the reflected-wave DFT amplitude at 10 GHz.
The analytic Fresnel coefficient Γ = (1 − n) / (1 + n) with n = √ε_r(ω_probe)
provides the reference.

```python
from yee import run_dispersive_drude

r = run_dispersive_drude()
print(r.gamma_measured, r.gamma_analytic, r.rel_err, r.passed)
# e.g. 0.5832 0.5941 0.0183 True
print(r)
# DispersiveDrudeResult(gamma_measured=0.5832, gamma_analytic=0.5941, rel_err=1.83%, passed=True)
assert r.passed, (
    f"dispersive-001 gate failed: rel_err={r.rel_err:.2%} > 20% "
    f"(measured={r.gamma_measured:.4f}, analytic={r.gamma_analytic:.4f})"
)
```

Gate criterion (dispersive-001, Taflove §9): `rel_err ≤ 20%`. The ADE Drude
update on a 1 mm cell at 10 GHz typically achieves well under 5% relative
error.

## 4  Running All Three Together

The three gates are independent and can be run in sequence in a single script:

```python
import yee

gates = {
    "cpml-001":        yee.run_cpml_reflection(),
    "ntff-001":        yee.run_ntff_broadside(),
    "dispersive-001":  yee.run_dispersive_drude(),
}

for name, r in gates.items():
    status = "PASS" if r.passed else "FAIL"
    print(f"{name}: {status}  ({r})")

all_passed = all(r.passed for r in gates.values())
assert all_passed, "One or more FDTD absorption gates failed"
print("All FDTD absorption gates passed.")
```

## Result Field Reference

### `CpmlReflectionResult`

| Field          | Description                                                          |
|----------------|----------------------------------------------------------------------|
| `reduction_db` | CPML reflection reduction vs PEC in dB (positive = CPML is better)  |
| `passed`       | `True` iff `reduction_db ≥ 30.0`                                     |

### `NtffResult`

| Field      | Description                                                             |
|------------|-------------------------------------------------------------------------|
| `ratio_db` | Broadside-to-endfire amplitude ratio in dB (θ=90° vs θ=0°)             |
| `passed`   | `True` iff `ratio_db ≥ 20.0`                                            |

### `DispersiveDrudeResult`

| Field            | Description                                                         |
|------------------|---------------------------------------------------------------------|
| `gamma_measured` | Measured \|Γ\| from FDTD DFT of reflected vs incident amplitude     |
| `gamma_analytic` | Analytic \|Γ\| from Fresnel (1−n)/(1+n) with Drude ε_r(ω_probe)    |
| `rel_err`        | \|gamma_measured − gamma_analytic\| / gamma_analytic                |
| `passed`         | `True` iff `rel_err ≤ 0.20`                                         |

## Connection to the Rust Gates

The underlying Rust unit tests live in `crates/yee-validation/src/lib.rs`:

```bash
# cpml-001 gate (fast, ~300 steps × 2):
cargo test -p yee-validation --lib -- cpml001 --release

# ntff-001 gate (2000 steps):
cargo test -p yee-validation --lib -- ntff001 --release

# dispersive-001 gate (800 steps × 2):
cargo test -p yee-validation --lib -- dispersive001 --release
```

## See also

- [FDTD TF/SF Fresnel transmission from Python](13-fdtd-fresnel-tfsf-from-python.md)
- [FDTD dipole radiation pattern from Python](12-fdtd-dipole-pattern-from-python.md)
- [FDTD cavity resonance from Python](11-fdtd-cavity-resonance-from-python.md)
- [FDTD lossy cavity from Python](10-fdtd-lossy-cavity-from-python.md)
- [FDTD theory details](../theory/fdtd-details.md)
- [ADR-0074: Phase 1.validation.2 FDTD gate integration (cpml-001, ntff-001, dispersive-001)](../decisions/0074-phase-1-validation-2-fdtd-aggregator-gates.md)
- [ADR-0077: Phase 2.fdtd.py.4 FDTD absorption Python drivers](../decisions/0077-phase-2-fdtd-py-4-absorption-python.md)
