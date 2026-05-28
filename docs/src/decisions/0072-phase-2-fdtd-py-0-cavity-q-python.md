# ADR-0072 — Phase 2.fdtd.py.0: Python Cavity-Q Driver + fdtd-202 Aggregator Gate

**Status:** Accepted (2026-05-26)  
**Phase:** 2.fdtd.py.0

---

## Context

Phase 2.fdtd.8 (ADR-0071) shipped per-cell electric conductivity (σ) and the
fdtd-202 Q-factor integration test (Q ≈ 20, rel err 0.04 %, 0.38 s wall-time,
NOT `#[ignore]`-gated).  Two gaps remained:

1. **No Python surface.** The `yee-py` `FdtdDriver` / `FdtdDriverConfig` /
   `RadiationPattern` API (Phase 1.frontend.2) covers only the dipole + NTFF
   radiation-pattern path.  Lossy cavity simulation requires direct use of
   `YeeGrid::set_sigma_box`, a stepping loop, and the exponential-decay fit —
   all Rust-only.  Python users (notebooks, surrogate sweeps) cannot access
   this capability.

2. **Not in `yee validate all`.** `Report::run_all()` does not include
   fdtd-202.  The three FDTD slots in `run_all()` (cpml-001, ntff-001,
   dispersive-001) are all `CaseStatus::Skipped` (long-running integration
   tests).  fdtd-202 is fast (0.38 s) and non-ignored — it qualifies for the
   fast aggregator path, unlike fdtd-201 / fdtd-201x.

---

## Decision

### Python driver (walking skeleton)

Add a `run_cavity_q()` free function and `CavityQResult` class to
`crates/yee-py/src/fdtd.rs`.  The function mirrors the fdtd-202 integration
test exactly — same geometry defaults (20×10×20 cells, dx=10 mm), same source
(off-centre Gaussian, 200 steps), same ring-down extraction (6000 steps, last
2/3 window, log-linear fit).  All parameters are keyword arguments with the
fdtd-202 defaults.

**Rationale for high-level API:** A general `PyYeeGrid` (raw grid + stepping)
would require exposing the full `YeeGrid` builder, field arrays, and stepping
API to Python — a significant PyO3 surface with many edge cases.  Per the
walking-skeleton principle, a single-call `run_cavity_q()` delivers the
fdtd-202 use case immediately.  A general `PyYeeGrid` is queued for Phase
2.fdtd.py.1.

`CavityQResult` fields: `q_measured`, `q_analytic`, `f101_hz`, `rel_err`,
`passed`, `probe_array()` method (returns numpy f64 array).

### Validation aggregator registration

Add `run_fdtd_202_lossy_cavity_q() -> CaseResult` to
`crates/yee-validation/src/lib.rs`.  Physics helpers duplicated from
`cavity_q.rs` (not exported from `yee-fdtd` to avoid API churn).  Wire into
`Report::run_all()` between `run_dispersive_001()` and `run_fem_eig_001()`.

### Tutorial

`docs/src/tutorials/10-fdtd-lossy-cavity-from-python.md` — worked example
calling `run_cavity_q()`, explaining the ring-down physics (Taflove §3.7),
and showing a matplotlib plot of the probe decay.

---

## Consequences

- **`yee validate all` will show fdtd-202 as PASS** for the first time — the
  first Phase-2 FDTD gate visible in the aggregator output.
- **Python users** can run a lossy FDTD cavity in two lines and get a
  validated Q-factor result.
- **fdtd-201 / fdtd-201x** remain `#[ignore]`-gated integration tests (minutes
  wall-time); they do not enter `run_all()`.
- **General `PyYeeGrid`** (per-cell ε_r from Python, field array access, etc.)
  is deferred to Phase 2.fdtd.py.1.
- **fdtd-202 physics is unchanged** — this increment adds a Python surface
  and a validation aggregator entry; it does not touch `yee-fdtd` core.

---

## Alternatives considered

**Alt A — Export `run_lossy_cavity` from `yee-fdtd/src/`.**  
Rejected: would add a pub function to the `yee-fdtd` public API that exists
only for validation.  The validation crate duplicating 80 lines is simpler.

**Alt B — General `PyYeeGrid` first.**  
Rejected: larger scope, more edge cases (array ownership, step ordering,
boundary conditions), not needed for the fdtd-202 use case.  Walking skeleton
first.
