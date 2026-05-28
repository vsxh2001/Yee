# Phase 2.fdtd.py.4 — Python drivers for cpml-001 / ntff-001 / dispersive-001

**Date:** 2026-05-27
**ADR:** 0077
**Phase:** 2.fdtd.py.4

---

## 1. Context

ADR-0074 (Phase 1.validation.2) wired the three FDTD physics gates —
`cpml-001` (CPML reflection attenuation), `ntff-001` (NTFF broadside/endfire
ratio), and `dispersive-001` (Drude-slab Fresnel reflection) — into
`yee-validation::Report::run_all()` using real physics. Their `#[ignore]`-gated
unit tests exercise them from Rust, but no Python-callable counterparts exist
in `yee-py`.

The Phase 2.fdtd.py.0–3 series established a clear pattern:

| Phase       | Gate        | Python function               | Tutorial |
|-------------|-------------|-------------------------------|----------|
| 2.fdtd.py.0 | fdtd-202    | `run_cavity_q`                | 10       |
| 2.fdtd.py.1 | fdtd-201    | `run_cavity_resonance`        | 11       |
| 2.fdtd.py.2 | fdtd-203    | `run_dipole_pattern`          | 12       |
| 2.fdtd.py.3 | fdtd-204    | `run_fresnel_tfsf`            | 13       |
| **2.fdtd.py.4** | cpml-001/ntff-001/dispersive-001 | **`run_cpml_reflection`, `run_ntff_broadside`, `run_dispersive_drude`** | **14** |

This phase fills the gap by adding individual Python callables for all three
ADR-0074 gates in one track.

---

## 2. Gate physics recap

### cpml-001 — CPML reflection attenuation
- Grid: 50³ vacuum at dx=1 mm, CPML thickness 10 cells.
- Source: Gaussian E_z pulse at (25,25,25), t0=20·dt, σ=6·dt.
- Probe: (38,25,25) — captures the outgoing wave and any reflected tail.
- Two runs: PEC (no ABC) and CPML.
- Metric: `reduction_db = PEC_reflection_db − CPML_reflection_db` (positive = CPML beats PEC).
- **Gate: ≥ 30 dB** (Roden–Gedney 2000 target; currently measured ≥ 69 dB).
- Wall-time: < 1 s.

### ntff-001 — NTFF broadside/endfire ratio
- Grid: 50³ vacuum at dx=1 mm, CPML=10, box_margin=15.
- Source: Gaussian E_z at (25,25,25), drive frequency 15 GHz.
- Run: 2000 steps.
- NTFF evaluated at broadside (θ=π/2, φ=0) and endfire (θ=0, φ=0).
- Metric: `20·log10(|E_broadside|/|E_endfire|)` in dB.
- **Gate: ≥ 20 dB** (analytic E_z dipole broadside/endfire → ∞ at the null).
- Wall-time: < 5 s.

### dispersive-001 — Drude-slab Fresnel reflection
- Grid: 80³ at dx=1 mm, CPML=10, 800 steps.
- Drude material: ε_inf=1, ω_p=2π·20 GHz, γ=2π·5 GHz, slab i∈[50,70).
- Source and probe at x=(20,30), DFT at 10 GHz.
- Metric: `|Γ_measured / Γ_analytic − 1|` (relative error).
- Analytic: Fresnel `Γ = (1−n)/(1+n)`, n=√ε_r(ω_probe).
- **Gate: rel_err ≤ 20%** (ADE Drude tolerance per Taflove §9).
- Wall-time: < 10 s.

---

## 3. Design decisions

### D1 — Make physics functions public in yee-validation

`cpml001_run() -> f64`, `ntff001_run() -> f64`,
`dispersive001_run() -> (f64, f64)`, and
`dispersive001_fresnel_gamma(eps_r: Complex64) -> Complex64`
are currently private. Making them `pub` is the minimal change needed and
matches the precedent of `fdtd204_t_analytic` (also a `pub fn` in
`yee-validation`).

The full `run_cpml_001()` / `run_ntff_001()` / `run_dispersive_001()` wrappers
(which return `CaseResult`) stay private — Python callers don't need the
aggregator machinery.

### D2 — Inline physics in yee-py (no direct delegation)

Following py.0–3 precedent, `yee-py` duplicates the default-parameter physics
calls inline (no cross-crate delegation), but uses the public helpers from
`yee-validation` for analytic references (like `fdtd204_t_analytic`). For
these three gates the analytic reference is only needed for `dispersive-001`
(`dispersive001_fresnel_gamma`), so we import just that.

Actually, the simplest approach for cpml-001 and ntff-001 is to call the pub
`cpml001_run()` / `ntff001_run()` directly from `yee-py`, since the physics
encapsulated there is stable and non-trivially complex (two-run setup for CPML,
NTFF state machinery for ntff). For `dispersive-001` we similarly call
`dispersive001_run()` and `dispersive001_fresnel_gamma()`.

### D3 — Result types

Three new Python result classes:
- `PyCpmlReflectionResult`: `reduction_db`, `passed` (≥30 dB).
- `PyNtffResult`: `ratio_db`, `e_broadside`, `e_endfire`, `passed` (≥20 dB).
- `PyDispersiveDrudeResult`: `gamma_measured`, `gamma_analytic`, `rel_err`, `passed` (≤20%).

### D4 — Registration in run_all()

These three gates already run inline in `run_all()` (not Skipped) as of
ADR-0074. No changes needed to `run_all()` registration.

### D5 — Tutorial 14

`docs/src/tutorials/14-fdtd-absorption-validation-from-python.md` — combines
all three in one tutorial ("FDTD absorption validation from Python").

---

## 4. Definition of Done

1. `cargo test -p yee-py --test test_fdtd -- cpml ntff dispersive 2>&1 | grep "test result: ok"` exits 0.
2. `cargo test -p yee-validation --test integration 2>&1 | grep "test result: ok"` exits 0 (no regression).
3. `cargo clippy --workspace --all-targets -- -D warnings` exits 0.
4. `cargo fmt --check --all` exits 0.
5. `crates/yee-py/tests/test_fdtd.py` has at least 2 test functions per driver (API plumbing + gate check).
6. `docs/src/tutorials/14-fdtd-absorption-validation-from-python.md` exists and is linked from `docs/src/SUMMARY.md`.
7. `docs/src/decisions/0077-phase-2-fdtd-py-4-absorption-python.md` exists.

---

## 5. Out of scope

- Changing cpml-001/ntff-001/dispersive-001 gate tolerances.
- Changing the `run_all()` registration or Skipped/Passed status of these gates.
- Adding parameters to the physics (e.g., variable grid size) beyond the minimum needed for smoke tests.
- Tutorial translation or internationalization.
