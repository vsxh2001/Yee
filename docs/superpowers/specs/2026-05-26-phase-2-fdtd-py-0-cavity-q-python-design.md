# Phase 2.fdtd.py.0 — Python Cavity-Q Driver + fdtd-202 Aggregator Registration

**Phase:** 2.fdtd.py.0  
**Date:** 2026-05-26  
**Status:** proposed  
**ADR:** [ADR-0072](../../../docs/src/decisions/0072-phase-2-fdtd-py-0-cavity-q-python.md)  
**Plan:** [2026-05-26-phase-2-fdtd-py-0-cavity-q-python.md](../plans/2026-05-26-phase-2-fdtd-py-0-cavity-q-python.md)

---

## 1  Motivation

Phase 2.fdtd.8 shipped per-cell electric conductivity (σ) in `YeeGrid` and
validated it via the fdtd-202 integration test (Q-factor ring-down, 0.04 %
error, 0.38 s wall-time).  However, two gaps remain:

1. **No Python surface.** The `yee-py` `FdtdDriver` only exposes the
   dipole+NTFF radiation-pattern path (Phase 1.frontend.2).  Python users
   have no way to run a lossy cavity, set per-cell conductivity, or extract a
   Q-factor.  The entire fdtd-202 physics is Rust-only.

2. **Not in the validation aggregator.** `yee validate all` / `Report::run_all()`
   does not include fdtd-202.  The three FDTD slots in `run_all()` (cpml-001,
   ntff-001, dispersive-001) are all `CaseStatus::Skipped` (long-running
   integration tests).  fdtd-202 runs in 0.38 s and is not `#[ignore]`-gated —
   it qualifies for the fast path.

This increment closes both gaps in one well-scoped track.

---

## 2  Scope

### 2.1  Python driver (walking skeleton)

Add a high-level `run_cavity_q()` free function to `yee-py/src/fdtd.rs` (same
module as `PyFdtdDriver`) that mirrors the fdtd-202 Rust test, returning a
`PyCavityQResult` PyO3 class:

```python
from yee import run_cavity_q

result = run_cavity_q(
    nx=20, ny=10, nz=20,
    dx=0.01,          # metres
    sigma=2.96e-3,    # S/m
    n_src=200,        # source-injection steps
    n_ring=6000,      # ring-down steps to record
)
print(result.q_measured)   # ≈ 20.0
print(result.rel_err)      # ≈ 0.0004
print(result.passed)       # True (rel_err < 0.05)
# result.probe_array — numpy f64 array, the ring-down time series
```

`PyCavityQResult` fields (all `#[pyo3(get)]`):
| Field | Type | Description |
|---|---|---|
| `q_measured` | `f64` | Q extracted from log-linear ring-down fit |
| `q_analytic` | `f64` | Analytic Q = 2π·f₁₀₁·ε₀/σ |
| `f101_hz` | `f64` | Analytic TE₁₀₁ resonant frequency (Hz) |
| `rel_err` | `f64` | `|q_measured − q_analytic| / q_analytic` |
| `passed` | `bool` | `rel_err < 0.05` |
| `probe_array` | `PyArray1<f64>` | ring-down time series (n_ring samples) |

### 2.2  Validation aggregator registration

Add `run_fdtd_202_lossy_cavity_q() -> CaseResult` to
`crates/yee-validation/src/lib.rs` using the same physics logic
(duplicated from cavity_q.rs — no cross-crate pub leak needed).  Wire into
`Report::run_all()` between `run_dispersive_001()` and `run_fem_eig_001()`.

The gate is `|Q_measured − Q_analytic| / Q_analytic < 5 %` (same as the
integration test).  Case ID: `"fdtd-202"`.

### 2.3  Tutorial

`docs/src/tutorials/10-fdtd-lossy-cavity-from-python.md` — a worked example
that:
1. Imports `yee.run_cavity_q`.
2. Calls it with the canonical σ₀ = 2.96e-3 S/m.
3. Explains the physics (TE₁₀₁ ring-down, Q = ε₀ω/σ, Taflove §3.7).
4. Shows how to plot the ring-down with matplotlib.

Register in `docs/src/SUMMARY.md`.

---

## 3  Design decisions

### D1 — High-level driver, not raw YeeGrid API

Exposing `YeeGrid` as a general-purpose Python object is a larger scope
(general material builders, stepping loop, field access arrays).  Per the
walking-skeleton principle, ship the minimal useful surface first.
`run_cavity_q()` is a single-call API that covers the canonical fdtd-202 use
case.  A general `PyYeeGrid` can be a follow-on (Phase 2.fdtd.py.1).

### D2 — Duplicate physics logic in yee-validation (no pub leak)

The fdtd-202 helper functions (`run_lossy_cavity`, `fit_log_decay`,
`analytic_f101`, `analytic_q`) live in the private integration test.  Exporting
them as `pub` from `yee-fdtd` would change the crate's public API contract.
Instead, the validation crate duplicates the ~80 lines needed — the logic is
simple and self-contained.  If a `yee-fdtd::analysis` module is added later,
the duplicated code moves there.

### D3 — `probe_array` exposed for downstream plotting

Returning the full ring-down time series lets Python callers plot the decay,
inspect transient behavior, and verify the fit.  It is not needed for the gate,
but it is the kind of "show me the data" affordance that makes the Python API
useful in a Jupyter notebook.

---

## 4  Non-goals (deferred)

- General `PyYeeGrid` with full material builder API — Phase 2.fdtd.py.1.
- Per-cell ε_r from Python — Phase 2.fdtd.py.1.
- fdtd-201 / fdtd-201x in the aggregator — those are `#[ignore]`-gated
  (minutes-scale) and do not qualify for the fast aggregator path.
- Plotters integration (PNG artifact from the aggregator) — follow-on.

---

## 5  Validation (DoD)

All items must be machine-checkable before the track declares done:

| # | Check | Command | Exit |
|---|---|---|---|
| V1 | Rust workspace builds | `cargo check --workspace` | 0 |
| V2 | Clippy clean | `cargo clippy --workspace --all-targets -- -D warnings` | 0 |
| V3 | Format clean | `cargo fmt --check --all` | 0 |
| V4 | Validation aggregator includes fdtd-202 PASS | `cargo test -p yee-validation --lib -- run_fdtd_202 --nocapture` | 0 |
| V5 | Python wheel builds | `maturin develop --manifest-path crates/yee-py/Cargo.toml` | 0 |
| V6 | Cavity-Q pytest passes | `pytest crates/yee-py/tests/test_fdtd.py -k cavity_q -v` | 0 |
| V7 | run_all reports fdtd-202 | `cargo test -p yee-validation --lib -- run_all --nocapture 2>&1 \| grep fdtd-202` | 0 (line present) |
