# ADR-0076: Phase 2.fdtd.py.3 FDTD TF/SF Fresnel-Transmission Python Driver + fdtd-204 Gate

**Date:** 2026-05-27  
**Status:** Accepted  
**Phase:** 2.fdtd.py.3  

---

## Context

Phases 2.fdtd.py.0–2 established the pattern of exposing FDTD validation
gates to Python with a `run_*()` function and a `*Result` class. The
TF/SF plane-wave source (Phase 2.fdtd.5.3.2) and the per-cell ε_r
infrastructure (MMMMMMMM commit) have both shipped but have never been
tested together in a quantitative validation case.

The `heterogeneous_substrate.rs` test in `yee-fdtd` validates per-cell ε_r
Fresnel reflection using a **point source**; there is no existing test
combining TF/SF with per-cell ε_r. This ADR captures the decision to add
`fdtd-204` as that combined test.

---

## Decision

**Ship `run_fresnel_tfsf()` / `FresnelTfsfResult` in `yee-py` and register
`fdtd-204` in the `yee-validation` aggregator.**

### Physical scenario

Normal-incidence TF/SF plane wave at 10 GHz through a 5-cell (5 mm) lossless
dielectric slab (ε_r = 2.2) in an 80³ grid with CPML.

The slab is near-quarter-wave (δ ≈ π/2), giving |t_analytic| ≈ 0.927 — well
measurable without Fabry-Perot cancellation effects.

### Measurement

Two runs (vacuum + slab) with the same `PlaneWaveSource`. The transmitted
amplitude is the peak |E_z| at a probe in the vacuum region after the slab;
the incident amplitude is the peak |E_z| at a probe in the vacuum region
before the slab. The ratio is compared to the analytic transfer-matrix
prediction.

### Gate criterion

`|t_measured / t_analytic − 1| < 0.05` (5%)

### Wall-time

~5–15 min release (two 80³ × 600-step runs). Registered Skipped in
`run_all()`. A `#[ignore]`-gated unit test in `yee-validation` enables
manual gate verification.

---

## Alternatives Considered

1. **Reflect instead of transmit**: The reflection coefficient requires a
   two-run SUBTRACTION to isolate the reflected wave, which complicates the
   measurement. The transmission approach uses a direct amplitude ratio.

2. **Half-space slab extending to CPML**: Would create a dielectric-CPML
   impedance mismatch (r ≈ 0.195 amplitude) generating a spurious backward
   wave inside the slab. The 5-cell slab with vacuum buffer avoids this.

3. **Thicker slab (half-wave at 10 GHz)**: A 10-cell slab at 10 GHz is near
   a half-wave transformer, giving near-zero overall reflection and
   |t_analytic| ≈ 1. This is poorly measured and provides little sensitivity.
   The 5-cell quarter-wave slab at |t| ≈ 0.927 is both measurable and
   physically meaningful.

4. **Use DFT instead of peak amplitude**: For CW TF/SF, both approaches give
   the same result in steady state. Peak amplitude is simpler to implement.

---

## Consequences

- **New capability demonstrated**: TF/SF + per-cell ε_r validated together
  for the first time in a quantitative gate.
- **SUMMARY.md gap fixed**: ADR-0075 (missing from prior session) is added.
- **Tutorial 13** (`13-fdtd-fresnel-tfsf-from-python.md`) extends the
  Python FDTD tutorial series.
- **Pattern consistency**: follows py.0/1/2 structure exactly.

---

## References

- Born & Wolf, *Principles of Optics*, §1.6.2 (slab transfer matrix)  
- Pozar, *Microwave Engineering*, 4th ed., §2.3 (dielectric slab)  
- `crates/yee-fdtd/tests/heterogeneous_substrate.rs` — point-source Fresnel reference  
- `crates/yee-fdtd/tests/plane_wave_propagation.rs` — TF/SF quiet-zone reference  
- ADRs 0071–0075: fdtd-202/201/201-x/cpml-203 gates (py.0/1/2 pattern)
