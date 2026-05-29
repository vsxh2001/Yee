# ADR-0079: Phase 2.fdtd.py.5 — Ohmic Skin-Depth Python Driver

**Date:** 2026-05-28
**Status:** Accepted
**Supersedes:** —
**Superseded by:** —

---

## Context

Phase 2.fdtd.9 (ADR-0078) shipped the `fdtd-205` Ohmic skin-depth spatial
penetration gate in `yee-validation`:

- `fdtd205_run() → SkinDepthResult` (public)
- `run_fdtd_205() → CaseResult` (private, registered in `Report::run_all()`)
- Gate A: `rel_err_1δ < 10 %` (measured 1.05 %)
- Gate B: `rel_err_2δ < 15 %` (measured 2.22 %)
- Scenario: 5×5×130 grid, σ = 2.5331 S/m, f = 1 GHz, δ = 10 mm = 10 cells

The Python API pattern (Phases 2.fdtd.py.0–4) exposes each gate via a
`run_*()` function and a matching result class in `yee-py`. fdtd-205 is the
only gate not yet exposed.

---

## Decision

Expose `fdtd205_run()` to Python as:

- `yee.run_skin_depth() → SkinDepthResult` — delegates directly to
  `yee_validation::fdtd205_run()` with no parameter overrides (the scenario
  is fully determined by the validation constants). No GIL release needed:
  the call takes ~8 s debug, but is short enough that the overhead of
  `py.detach()` is not justified for a one-shot validation function.
- `yee.SkinDepthResult` — `#[pyclass]` with all 9 public scalar fields
  exposed via `#[pyo3(get)]` and a `__repr__` matching the established
  fdtd.py.* style.

`run_all()` registration is unchanged (already done in Phase 2.fdtd.9).

---

## Alternatives Considered

1. **Add σ / f / grid kwargs to `run_skin_depth`** — rejected: the scenario
   is frozen by the validation gate; parameterisation would require duplicating
   the PMC-boundary logic and risks a test diverging from the canonical gate.

2. **GIL release via `py.detach()`** — not needed for an 8 s call invoked once
   in a validation context (unlike `run_dipole_pattern` which is a long-running
   simulation). Omitting it keeps the code simpler.

---

## Consequences

- `from yee import run_skin_depth, SkinDepthResult` works after `maturin develop`.
- `pytest crates/yee-py/tests/test_fdtd.py -k skin_depth` → 3/3 pass.
- fdtd-205 gate status is consistent between Rust (`yee validate all`) and
  Python (`run_skin_depth().passed`).
- No new dependency. No yee-fdtd or yee-validation changes.
