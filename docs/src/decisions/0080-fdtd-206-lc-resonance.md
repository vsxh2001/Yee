# ADR-0080: fdtd-206 lumped series-LC resonant frequency gate (Phase 2.fdtd.6.1)

**Status:** Accepted  
**Date:** 2026-05-29  
**Supersedes:** none  
**Related:** ADR-0017 (Phase 2.fdtd.6 lumped RLC port + explicit 2.fdtd.6.1 deferral)

---

## Context

ADR-0017 shipped the `LumpedRlcPort` in Phase 2.fdtd.6 with an energy-dissipation
validation gate for the pure-resistor path. The series-RLC path ("compiles and
self-tests") was explicitly deferred to Phase 2.fdtd.6.1 with the note:

> "Phase 2.fdtd.6.1 is non-optional. The energy-dissipation gate is sufficient
> for a placeholder but not sufficient to call the lumped-port subsystem
> 'validated' in the sense the rest of `yee-fdtd` is validated."

ADR-0017 framed Phase 2.fdtd.6.1 as: calibrated TEM stripline + Γ gate + multi-
axis ports + parallel-RLC topology. The Γ gate requires an upstream calibrated
TEM stripline with a known Z₀, which `yee-fdtd` does not yet have.

However, the series-RLC ODE accuracy can be validated directly without a
stripline: the **natural resonant frequency f₀ = 1/(2π√LC)** of a series-LC
circuit is an exact analytic reference that depends only on the lumped-element
ODE integration, not on field propagation or Z₀ calibration. This is a complete
§4-compliant published-benchmark gate against the exact circuit-theory formula
(Pozar §2.4; Hayt & Kemmerly §14.1).

---

## Decision

Ship **fdtd-206** as Phase 2.fdtd.6.1-a (ODE-accuracy sub-step of the full
2.fdtd.6.1):

- **Geometry:** 5×5×40 cell PEC box at dx=1 mm; series-LC port at the centre
  cell (2,2,20) with L=1 nH, C=25.330 pF → f₀=1 GHz analytic, R=1 Ω → Q≈6.28.
- **Gate:** DFT scan of the inductor-current ring-down; peak frequency within
  ±2 % of the analytic 1/(2π√LC)=1 GHz.
- **Wall time:** < 0.1 s (5 000 steps, 5×5×40 grid). NOT `#[ignore]`-gated.
  Registered in `run_all()` as Passed.
- **Python wrapper:** `run_lc_resonance()` / `LcResonanceResult` in `yee-py`
  (Phase 2.fdtd.py.6). Tutorial 16.

### Why not the full Γ gate?

The Γ-against-analytic gate (Γ = (Z_L − Z₀)/(Z_L + Z₀)) requires a calibrated
TEM stripline, which is a separate infrastructure increment. The f₀ gate closes
the ODE-accuracy gap without that infrastructure. The Γ gate is the next
increment in 2.fdtd.6.1 and is explicitly noted below as the follow-on.

---

## Consequences

**What ships (fdtd-206):**
- `LcResonanceResult` + `fdtd206_run()` in `yee-validation`.
- `run_fdtd_206_lumped_lc_resonance()` in `run_all()` → Passed.
- Self-contained integration test `yee-fdtd/tests/lumped_lc_resonance.rs`.
- Python `run_lc_resonance()` + `LcResonanceResult`.
- Tutorial 16.

**What remains deferred (Phase 2.fdtd.6.1-b and beyond):**
- Γ gate (requires calibrated TEM stripline).
- Multi-axis / oriented-arbitrary lumped ports.
- Parallel-RLC topology.

These are the Phase 2.fdtd.6.1-b follow-ons; the ODE accuracy is now
established so the follow-on only needs to add the stripline geometry.

---

## References

- Pozar, "Microwave Engineering," 4th ed., §2.4 (series resonator).
- Hayt & Kemmerly, "Engineering Circuit Analysis," §14.1.
- Taflove & Hagness, "Computational Electrodynamics," 3rd ed., §15.10.
- ADR-0017 (Phase 2.fdtd.6 scope and 2.fdtd.6.1 deferral).
- `docs/superpowers/specs/2026-05-29-fdtd-206-lc-resonance-design.md`
- `docs/superpowers/plans/2026-05-29-fdtd-206-lc-resonance.md`
