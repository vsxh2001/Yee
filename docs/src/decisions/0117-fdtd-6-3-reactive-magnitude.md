# ADR-0117: Phase 2.fdtd.6.3 — reactive-magnitude correctness of the two-way lumped port

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0116 (2.fdtd.6.2 two-way lumped port — this completes its punted
reactive DoD), ADR-0115 (F2.3 lumped FDTD — still blocked on this), the
lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0116 shipped the stable two-way semi-implicit lumped R-L-C port. Its gate
`lumped_rlc_twoway_001` asserts the **resistive** Γ exactly (Δ|Γ| ≤ 0.07 at
4/6/9 GHz) and stability (no NaN), but only **prints** the reactive |Γ| — the
reactive magnitude was deferred as a "measurement follow-on."

Wiring F2.3 (ADR-0115) onto that primitive surfaced the deferred piece as the
real blocker. The lumped-LC filter's FDTD |S21| comes out **flat ≈ 1.0 across
the whole band** (zero selectivity). Isolating the primitive with the gate's own
single-load sweep shows the cause is **not** the F2.3 driver geometry — it is the
reactive update coefficients:

| load | |Γ|_fdtd (4/6/9 GHz) | |Γ|_analytic | verdict |
|------|---------------------|--------------|---------|
| resistor (R = Z₀_eff) | 0.370 / 0.333 / 0.322 | 0.333 | **exact** |
| inductor (L = 13.2 nH) | 0.013 / 0.015 / 0.019 | 0.60 / 0.45 / 0.32 | **transparent** (under-couples) |
| capacitor (C = 53.4 fF) | 1.080 / 0.971 / 0.936 | 0.32 / 0.45 / 0.60 | **near-open** (over-couples) |
| series R-L-C | 0.023 / 0.022 / 0.024 | ~0.29–0.33 | transparent |

The inductor reflects almost nothing; the capacitor reflects almost everything.
Because the **resistor limit is exact** and the update is stable, this is a
**coefficient/scaling bug in the L and C terms**, not the ill-posedness that
defers the MoM microstrip port (ADR-0064) — the field↔lumped geometric mapping
(`β = dt·dz/(2ε₀·dA)`, `V = E·dz`) is correct for R, so the bug is the `dt`/`dz`/
`dA`/`ε₀` grouping on the reactive (`L/dt`, `dt/(2C)`, `dt/C`) coefficients.

## Decision

Fix the reactive coefficients of `LumpedRlcPort::correct_e` (two-way path) so the
discrete impedance `Z_d(ω)` the update realizes matches `R + jωL + 1/(jωC)`. The
disciplined route is to **derive `Z_d(ω)`** from the semi-implicit branch update
(z-transform), compare term-by-term to the continuous impedance, and correct the
mis-scaled L/C factor. Keep the exact resistor limit and the public API.

**Strengthen the gate, do not weaken it:** `lumped_rlc_twoway_001` must now
**assert** the reactive |Γ| (pure-L, pure-C, series-RLC) within a loose tol at the
sweep frequencies — converting the deliberately-punted print into a real check.
The resistor-exact and stability assertions stay.

## Consequences

**Ships:** a quantitatively-correct reactive two-way lumped-element FDTD port —
completing ADR-0116's deferred DoD. **Unblocks F2.3**: with reactive |Γ| correct,
the L‖C resonators load the line and the filter |S21| acquires its band-pass
shape, so `fdtd_lumped_001` passes (the goal's "EM simulation" component).

**Gate:** `lumped_rlc_twoway_001` GREEN in CI with the new reactive-magnitude
assertions; resistor + stability non-regressed.

**Escape hatch (recorded honestly):** if, after deriving `Z_d(ω)`, the reactive
|Γ| cannot be brought within the loose tol without destabilising the update, the
precise partial (the derived `Z_d(ω)`, the attempted fix, the residual table) is
surfaced and this becomes a deferred research item — the gate is **not** weakened
back to a print, and F2.3 stays unmerged. Faking a pass is forbidden.

**Not in scope:** the F2.3 board sim itself (rides on this); SRF/ESR parasitics;
multi-element parasitic coupling.

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-3-reactive-magnitude-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-3-reactive-magnitude.md`.
- Taflove & Hagness, *Computational Electrodynamics*, lumped-element FDTD
  (Piket-May 1994) — the canonical R/L/C `E_z` update and its discrete impedance.
