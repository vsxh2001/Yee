# Phase 2.fdtd.py.2 — Python Dipole-Pattern Driver + fdtd-203 Gate

**Date:** 2026-05-27  
**Status:** Proposed  
**Related ADR:** ADR-0075  
**Phase:** 2.fdtd.py.2  

---

## 1. Context

`FdtdDriver` (in `crates/yee-fdtd/src/driver.rs`) is a complete end-to-end
FDTD radiation-pattern solver that wires a vacuum `YeeGrid` + CPML +
Hann-ramped sinusoidal J_z dipole source + `NtffState` into a single `run()`
call that returns a 37-point θ-cut of `|E_θ|` normalised to its maximum.

The driver is already partially wrapped for Python in `crates/yee-py/src/fdtd.rs`
(`PyFdtdDriver`, `PyFdtdDriverConfig`, `PyRadiationPattern`), giving advanced
users low-level control.  What is **missing** is:

1. A high-level convenience function `run_dipole_pattern()` (mirroring
   `run_cavity_q` / `run_cavity_resonance`) that runs the gate scenario with
   sensible defaults and returns a `DipolePatternResult` with a `passed` flag.
2. Validation gate **fdtd-203** registered in `Report::run_all()` (as `Skipped`
   due to wall-time, like fdtd-201).
3. Tutorial `12-fdtd-dipole-pattern-from-python.md`.

The corresponding `crates/yee-fdtd/tests/dipole_pattern.rs` integration test
already exists (`#[ignore]`-gated) and defines the gate criteria against the
Balanis §4.2 short-dipole `E_θ ∝ sin θ` reference.

This delivers a Phase 2 FDTD **radiation-pattern validation milestone** and
extends the Python API into the NTFF / far-field domain.

---

## 2. Reference

- Balanis, *Antenna Theory: Analysis and Design*, 4th ed., Wiley, 2016,
  §4.2 — infinitesimal dipole far-field `E_θ = j η₀ k I₀ ℓ e^{−jkr}/(4πr) · sin θ`.
  Pattern is axially symmetric; normalized `|E_θ|(θ) = sin θ` (maximum at
  θ = 90°, nulls at θ = 0° and 180°).

The gate criteria (lifted from `dipole_pattern.rs`) are:

| Point | θ | Expected | Tolerance |
|-------|---|----------|-----------|
| Endfire null  | 0°   | 0     | < 0.15    |
| 45° lobe      | 45°  | 0.707 | ±0.15     |
| Broadside peak| 90°  | 1.0   | ±0.05     |
| 135° lobe     | 135° | 0.707 | ±0.15     |
| Endfire null  | 180° | 0     | < 0.15    |

These are loose; the 60³ grid at 1 GHz (λ/dx = 60) gives only a coarse
staircase approximation. Tighter tolerances belong to a finer-grid gate in
a future increment.

---

## 3. Decision

### 3.1 `run_dipole_pattern()` pyfunction

Added to `crates/yee-py/src/fdtd.rs`, following the `run_cavity_q` /
`run_cavity_resonance` pattern:

```rust
#[pyfunction]
#[pyo3(signature = (
    nx = 60usize,
    ny = 60usize,
    nz = 60usize,
    dx = 5.0e-3_f64,
    n_steps = 800usize,
    source_freq_hz = 1.0e9_f64,
))]
pub fn run_dipole_pattern(
    nx, ny, nz, dx, n_steps, source_freq_hz,
) -> PyDipolePatternResult
```

Default scenario matches `dipole_pattern.rs`:
- 60³ grid, dx = 5 mm
- 800 steps, 1 GHz source
- Dipole at centre (30,30,30), 5 cells long
- ntff_surface_pad = 4, cpml_thickness = 10

The function builds a `YeeGrid::vacuum`, wraps it in `FdtdDriver`, calls
`driver.run()`, extracts the five sample points by finding the nearest θ value
in the returned `theta_deg` array, applies the gate criteria, and returns
`PyDipolePatternResult`.

The GIL is released during the FDTD time loop (same as `PyFdtdDriver::run`).

### 3.2 `PyDipolePatternResult` struct

```rust
#[pyclass(name = "DipolePatternResult", module = "yee._yee")]
pub struct PyDipolePatternResult {
    #[pyo3(get)] pub e_theta_0: f64,      // |E_θ|(0°)
    #[pyo3(get)] pub e_theta_45: f64,     // |E_θ|(45°)
    #[pyo3(get)] pub e_theta_90: f64,     // |E_θ|(90°) — normalized peak
    #[pyo3(get)] pub e_theta_135: f64,    // |E_θ|(135°)
    #[pyo3(get)] pub e_theta_180: f64,    // |E_θ|(180°)
    #[pyo3(get)] pub passed: bool,
    pub theta_deg_vec: Vec<f64>,
    pub e_theta_vec: Vec<f64>,
}
```

With `theta_deg_array()` and `e_theta_array()` methods returning numpy arrays.

### 3.3 Module registration

In `crates/yee-py/src/lib.rs`:
```rust
m.add_class::<fdtd::PyDipolePatternResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_dipole_pattern, m)?)?;
```

### 3.4 fdtd-203 validation gate

In `crates/yee-validation/src/lib.rs`, a new
`run_fdtd_203_dipole_pattern()` function returns:

```rust
CaseResult {
    id: "fdtd-203".into(),
    description: "FDTD short-dipole sin θ NTFF radiation pattern (Balanis §4.2)".into(),
    status: CaseStatus::Skipped,
    notes: "wall-time ~30 s release — run via cargo test \
            -p yee-fdtd --test dipole_pattern --release -- --include-ignored \
            or from Python: from yee import run_dipole_pattern; \
            assert run_dipole_pattern().passed".into(),
    ..Default::default()
}
```

Registered in `Report::run_all()` alongside fdtd-201 and fdtd-201-x.

### 3.5 Tutorial 12

`docs/src/tutorials/12-fdtd-dipole-pattern-from-python.md` — short
worked example:

```python
from yee import run_dipole_pattern
import numpy as np, matplotlib.pyplot as plt

r = run_dipole_pattern()
print(f"θ=0°:  {r.e_theta_0:.4f}  (expected ~0, null)")
print(f"θ=90°: {r.e_theta_90:.4f}  (expected ~1, broadside peak)")
assert r.passed
```

Registered in `docs/src/SUMMARY.md` under Tutorials.

---

## 4. DoD (machine-checkable)

1. `cargo build -p yee-py` exits 0 — new Rust code compiles.
2. `cargo test -p yee-py` exits 0 — existing pytest compatibility.
3. `cargo clippy --workspace --all-targets -- -D warnings` exits 0.
4. `cargo fmt --check --all` exits 0.
5. `python -c "from yee import run_dipole_pattern, DipolePatternResult; \
   r = run_dipole_pattern(n_steps=5); \
   assert hasattr(r, 'passed'); print('smoke OK')"` exits 0 (quick
   smoke, n_steps=5 is enough to check plumbing; pattern not meaningful
   but struct round-trips).
6. `python -c "from yee import run_dipole_pattern; \
   r = run_dipole_pattern(); assert r.passed, f'gate failed: {r}'"` exits
   0 in release (`maturin develop --release`).
7. `cargo test -p yee-validation` exits 0 — fdtd-203 registered as Skipped.

---

## 5. Lane

Allowed paths:
- `crates/yee-py/**`
- `crates/yee-validation/src/lib.rs`
- `docs/src/tutorials/12-fdtd-dipole-pattern-from-python.md`
- `docs/src/SUMMARY.md`
- `docs/superpowers/specs/2026-05-27-phase-2-fdtd-py-2-dipole-pattern-python-design.md`
- `docs/superpowers/plans/2026-05-27-phase-2-fdtd-py-2-dipole-pattern-python.md`
- `docs/src/decisions/0075-phase-2-fdtd-py-2-dipole-pattern-python.md`

**Out of lane (surface as findings):**
- `crates/yee-fdtd/src/**` — FdtdDriver already complete; no changes needed.
- `crates/yee-fdtd/tests/**` — dipole_pattern.rs already has the `#[ignore]`-gated test.
- Any other crate.

---

## 6. Escape hatch

Blocked > 15 min → surface finding and stop. Do not weaken existing tests.
