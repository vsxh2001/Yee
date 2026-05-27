# ADR-0077: Phase 2.fdtd.py.4 — FDTD absorption validation Python drivers

**Date:** 2026-05-27
**Status:** Proposed → Shipped (post-merge)
**Supersedes:** N/A
**Relates to:** ADR-0074 (cpml/ntff/dispersive gate implementation), ADR-0072–0076 (Phase 2.fdtd.py.0–3)

---

## Context

ADR-0074 (Phase 1.validation.2) added real physics implementations for the
`cpml-001` (CPML reflection), `ntff-001` (NTFF broadside/endfire), and
`dispersive-001` (Drude-slab Fresnel) gates in `yee-validation`. These run
inline in `Report::run_all()` and their `#[ignore]`-gated Rust unit tests
exercise them from Rust.

However, unlike the fdtd-202/201/203/204 gates that were each exposed
individually to Python via Phase 2.fdtd.py.0–3, these three gates have no
Python-callable equivalents in `yee-py`. A Python user wanting to invoke
just the CPML reflection test or the NTFF broadside check must run
`run_validation()` (the full aggregator) or drop into Rust.

---

## Decision

Add three Python-callable driver functions in `yee-py`:

| Python function | Rust gate | Gate threshold |
|-----------------|-----------|----------------|
| `run_cpml_reflection()` | cpml-001 | reduction_db ≥ 30 dB |
| `run_ntff_broadside()`  | ntff-001 | ratio_db ≥ 20 dB |
| `run_dispersive_drude()` | dispersive-001 | rel_err ≤ 20% |

Each returns a structured Python result object with `.passed`.

**Implementation path:**
1. Export `cpml001_run()`, `ntff001_run()`, `dispersive001_run()`, and
   `dispersive001_fresnel_gamma()` as `pub fn` from `yee-validation`.
2. Add thin wrappers in `crates/yee-py/src/fdtd.rs` that call these pub
   functions and return PyO3 result objects.
3. Register in `crates/yee-py/src/lib.rs` and `crates/yee-py/python/yee/__init__.py`.
4. Add 6 pytest cases (2 per driver: API plumbing + gate assertion).
5. Add tutorial 14 combining all three gates.

The `run_all()` aggregator registration and gate status (cpml/ntff/dispersive
already run inline, not Skipped) are **not changed**.

---

## Consequences

**Positive:**
- Python users can call `run_cpml_reflection()`, `run_ntff_broadside()`,
  `run_dispersive_drude()` individually, matching the pattern of py.0–3.
- Tutorial 14 gives a single narrative combining absorbing-boundary, far-field,
  and dispersive-material validation.
- Closes the gap between the three ADR-0074 Rust gates and the Python API.

**Neutral:**
- The pub `cpml001_run()` etc. are implementation details exposed for binding
  purposes; they are not primary public API (documented with `/// pub for
  Python bindings only`).
- No change to gate tolerances or aggregator behavior.

**Risks:**
- The three physics functions are fast (< 10 s each), so no wall-time concern.
- The `pub` promotion of private functions increases the surface area of
  `yee-validation`'s API slightly; acceptable given the precedent of
  `fdtd204_t_analytic` which is `pub` for the same reason.
