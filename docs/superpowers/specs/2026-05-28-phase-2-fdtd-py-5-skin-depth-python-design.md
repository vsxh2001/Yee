# Spec: Phase 2.fdtd.py.5 — Ohmic Skin-Depth Python Driver

**Date:** 2026-05-28
**Status:** Proposed
**ADR:** [ADR-0079](../../../docs/src/decisions/0079-phase-2-fdtd-py-5-skin-depth-python.md)

---

## 1  Background

Phase 2.fdtd.9 shipped the `fdtd-205` Ohmic skin-depth spatial penetration
gate (ADR-0078) in `yee-validation`: `fdtd205_run()` / `SkinDepthResult`.
The gate is registered in `Report::run_all()` and passes (Gate A 1.05 %,
Gate B 2.22 %).

The established `fdtd.py.*` pattern exposes each Rust validation gate to Python
via a lightweight wrapper in `yee-py/src/fdtd.rs`. Phases 2.fdtd.py.0–4
shipped analogous bindings for fdtd-202, 201, 203, 204, cpml-001, ntff-001, and
dispersive-001. fdtd-205 is the remaining gate without a Python callable.

---

## 2  Scope

**In-lane:** `crates/yee-py/**`, `docs/**`

**Out-of-scope (not in this phase):**

- No changes to `yee-validation`, `yee-fdtd`, or `run_all()` — fdtd-205 is
  already registered.
- No new Rust gate logic.

---

## 3  Deliverables

| Artifact | Path |
|---|---|
| `PySkinDepthResult` pyclass | `crates/yee-py/src/fdtd.rs` |
| `run_skin_depth()` pyfunction | `crates/yee-py/src/fdtd.rs` |
| Module registration | `crates/yee-py/src/lib.rs` |
| Python export | `crates/yee-py/python/yee/__init__.py` |
| 3 pytest cases | `crates/yee-py/tests/test_fdtd.py` |
| Tutorial | `docs/src/tutorials/15-fdtd-ohmic-skin-depth-from-python.md` |
| ADR-0079 | `docs/src/decisions/0079-phase-2-fdtd-py-5-skin-depth-python.md` |
| SUMMARY.md updates | `docs/src/SUMMARY.md` |

---

## 4  API Design

```python
from yee import run_skin_depth, SkinDepthResult

r = run_skin_depth()
assert isinstance(r, SkinDepthResult)
assert r.passed  # Gate A (<10%) and Gate B (<15%) both pass
print(r.delta_analytic_m)   # 0.01 m (10 mm)
print(r.rel_err_1delta)     # ≈ 0.0105  (1.05%)
print(r.rel_err_2delta)     # ≈ 0.0222  (2.22%)
print(r)
# SkinDepthResult(delta_analytic_m=1.0000e-02, rel_err_1delta=1.0500e-02,
#                 rel_err_2delta=2.2200e-02, passed=True)
```

### `SkinDepthResult` fields

| Field | Type | Description |
|---|---|---|
| `delta_analytic_m` | `float` | Analytic δ = √(2/(ω μ₀ σ)) in metres |
| `amp_surface` | `float` | Peak \|E_x\| at the conductor surface |
| `amp_1delta` | `float` | Peak \|E_x\| one skin depth in |
| `amp_2delta` | `float` | Peak \|E_x\| two skin depths in |
| `ratio_1delta` | `float` | amp_1delta / amp_surface |
| `ratio_2delta` | `float` | amp_2delta / amp_surface |
| `rel_err_1delta` | `float` | \|ratio_1δ − e⁻¹\| / e⁻¹ (Gate A) |
| `rel_err_2delta` | `float` | \|ratio_2δ − e⁻²\| / e⁻² (Gate B) |
| `passed` | `bool` | Gate A (< 10 %) AND Gate B (< 15 %) |

---

## 5  Gate Criteria

Same as the Rust fdtd-205 gate (ADR-0078):

- **Gate A:** `rel_err_1delta < 0.10` (10 %)
- **Gate B:** `rel_err_2delta < 0.15` (15 %)

---

## 6  DoD (machine-checkable)

1. `cargo check -p yee-py` → exit 0
2. `cargo clippy -p yee-py -- -D warnings` → exit 0
3. `cargo fmt --check -p yee-py` → exit 0
4. After `maturin develop`:  
   `python -c "from yee import run_skin_depth; r = run_skin_depth(); assert r.passed"` → exit 0
5. `pytest crates/yee-py/tests/test_fdtd.py -k skin_depth -v` → 3/3 pass

---

## 7  Risk

Minimal — pure extension of the established pattern. `fdtd205_run()` is already
`pub` and all fields of `SkinDepthResult` are `pub`. The Rust call is fast (~8 s
debug) and is not `#[ignore]`-gated, so pytest runs are acceptable.
