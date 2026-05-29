# Phase 2.fdtd.py.4 — Implementation Plan

**Spec:** `docs/superpowers/specs/2026-05-27-phase-2-fdtd-py-4-absorption-python-design.md`
**ADR:** 0077
**Branch:** `feature/phase-2-fdtd-py-4-absorption-python`

---

## Step 1 — Expose physics helpers from yee-validation

**File:** `crates/yee-validation/src/lib.rs`

Change `fn cpml001_run()`, `fn ntff001_run()`, `fn dispersive001_run()`,
and `fn dispersive001_fresnel_gamma()` from private `fn` to `pub fn`.

Note: `cpml001_run_trace()` stays private (it's an implementation detail of
`cpml001_run()`).

---

## Step 2 — Add Python result types and driver functions

**File:** `crates/yee-py/src/fdtd.rs`

Append below the existing Phase 2.fdtd.py.3 block:

```
// ---------------------------------------------------------------------------
// Phase 2.fdtd.py.4 — cpml-001 / ntff-001 / dispersive-001 Python drivers
// ---------------------------------------------------------------------------
```

### 2a — cpml-001

```rust
#[pyclass(name = "CpmlReflectionResult", module = "yee._yee")]
pub struct PyCpmlReflectionResult {
    #[pyo3(get)] pub reduction_db: f64,
    #[pyo3(get)] pub passed: bool,
}
#[pymethods] impl PyCpmlReflectionResult { fn __repr__ ... }

#[pyfunction]
pub fn run_cpml_reflection() -> PyCpmlReflectionResult {
    let reduction_db = yee_validation::cpml001_run();
    PyCpmlReflectionResult {
        reduction_db,
        passed: reduction_db >= 30.0,
    }
}
```

### 2b — ntff-001

```rust
#[pyclass(name = "NtffResult", module = "yee._yee")]
pub struct PyNtffResult {
    #[pyo3(get)] pub ratio_db: f64,
    #[pyo3(get)] pub passed: bool,
}
#[pymethods] impl PyNtffResult { fn __repr__ ... }

#[pyfunction]
pub fn run_ntff_broadside() -> PyNtffResult {
    let ratio_db = yee_validation::ntff001_run();
    PyNtffResult {
        ratio_db,
        passed: ratio_db >= 20.0,
    }
}
```

### 2c — dispersive-001

```rust
#[pyclass(name = "DispersiveDrudeResult", module = "yee._yee")]
pub struct PyDispersiveDrudeResult {
    #[pyo3(get)] pub gamma_measured: f64,
    #[pyo3(get)] pub gamma_analytic: f64,
    #[pyo3(get)] pub rel_err: f64,
    #[pyo3(get)] pub passed: bool,
}
#[pymethods] impl PyDispersiveDrudeResult { fn __repr__ ... }

#[pyfunction]
pub fn run_dispersive_drude() -> PyDispersiveDrudeResult {
    let (gamma_measured, gamma_analytic) = yee_validation::dispersive001_run();
    let rel_err = (gamma_measured - gamma_analytic).abs() / gamma_analytic;
    PyDispersiveDrudeResult { gamma_measured, gamma_analytic, rel_err, passed: rel_err <= 0.20 }
}
```

---

## Step 3 — Register in yee-py lib.rs

**File:** `crates/yee-py/src/lib.rs`

Add after the existing `run_fresnel_tfsf` registrations:
```rust
m.add_class::<fdtd::PyCpmlReflectionResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_cpml_reflection, m)?)?;
m.add_class::<fdtd::PyNtffResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_ntff_broadside, m)?)?;
m.add_class::<fdtd::PyDispersiveDrudeResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_dispersive_drude, m)?)?;
```

---

## Step 4 — Update __init__.py

**File:** `crates/yee-py/python/yee/__init__.py`

Add to imports and `__all__`:
```python
CpmlReflectionResult,
DispersiveDrudeResult,
NtffResult,
run_cpml_reflection,
run_dispersive_drude,
run_ntff_broadside,
```

---

## Step 5 — Add pytest cases

**File:** `crates/yee-py/tests/test_fdtd.py`

Add a new section "Phase 2.fdtd.py.4" with at least 6 test functions:

```python
# ── cpml-001 ────────────────────────────────────────────────────────────────
def test_run_cpml_reflection_returns_result():
    r = yee.run_cpml_reflection()
    assert hasattr(r, "reduction_db")
    assert hasattr(r, "passed")

def test_run_cpml_reflection_passes_gate():
    r = yee.run_cpml_reflection()
    assert r.passed, f"cpml-001 gate failed: {r.reduction_db:.2f} dB < 30 dB"

# ── ntff-001 ─────────────────────────────────────────────────────────────────
def test_run_ntff_broadside_returns_result():
    r = yee.run_ntff_broadside()
    assert hasattr(r, "ratio_db")
    assert hasattr(r, "passed")

def test_run_ntff_broadside_passes_gate():
    r = yee.run_ntff_broadside()
    assert r.passed, f"ntff-001 gate failed: {r.ratio_db:.2f} dB < 20 dB"

# ── dispersive-001 ───────────────────────────────────────────────────────────
def test_run_dispersive_drude_returns_result():
    r = yee.run_dispersive_drude()
    assert hasattr(r, "gamma_measured")
    assert hasattr(r, "gamma_analytic")
    assert hasattr(r, "rel_err")
    assert hasattr(r, "passed")

def test_run_dispersive_drude_passes_gate():
    r = yee.run_dispersive_drude()
    assert r.passed, (
        f"dispersive-001 gate failed: rel_err={r.rel_err:.2%} > 20% "
        f"(measured={r.gamma_measured:.4f}, analytic={r.gamma_analytic:.4f})"
    )
```

---

## Step 6 — Write tutorial 14

**File:** `docs/src/tutorials/14-fdtd-absorption-validation-from-python.md`

Title: "FDTD absorption validation from Python"

Sections:
1. Introduction (what these three gates test and why they matter)
2. CPML reflection gate (`run_cpml_reflection`)
3. NTFF broadside/endfire gate (`run_ntff_broadside`)
4. Dispersive Drude material gate (`run_dispersive_drude`)
5. What the gates prove together (absorbing boundary + far-field transform + dispersive material)

---

## Step 7 — Update SUMMARY.md

**File:** `docs/src/SUMMARY.md`

Add after tutorial 13 entry:
```
- [FDTD absorption validation from Python](tutorials/14-fdtd-absorption-validation-from-python.md)
```

Add after ADR-0076 entry:
```
- [ADR-0077: Phase 2.fdtd.py.4 FDTD absorption Python drivers](decisions/0077-phase-2-fdtd-py-4-absorption-python.md)
```

---

## Step 8 — Write ADR-0077

**File:** `docs/src/decisions/0077-phase-2-fdtd-py-4-absorption-python.md`

Standard ADR format: context, decision, consequences.

---

## Verification command

```bash
cargo test -p yee-validation --test integration 2>&1 | grep "test result: ok"
cargo clippy --workspace --all-targets -- -D warnings 2>&1 | grep -v "^warning"
cargo fmt --check --all 2>&1
```

Expected: all exit 0.

The pytest cases in `test_fdtd.py` are verified by the CI `python-bindings` job
(maturin develop + pytest). Local verification via:
```bash
cd crates/yee-py && maturin develop --release && pytest tests/test_fdtd.py -k "cpml or ntff or dispersive" -v
```

---

## Pattern file

Imitate: `crates/yee-py/src/fdtd.rs` (Phase 2.fdtd.py.3 block near EOF) and
`crates/yee-py/tests/test_fdtd.py` (Phase 2.fdtd.py.3 section).
