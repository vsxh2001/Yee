# ADR-0129: Filter Phase F2.3-e — finer-grid lumped EM sim (toward the 20 dB gate)

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0128 (F2.3-d — coarse grid saturates the notch at ~5 dB; the
maintainer chose to invest in a finer grid, keeping the strict gate), ADR-0125
(the dx-stable aperture port), ADR-0124 (the dx-sweep that found the *single-cell*
port got WORSE with finer dx — no longer the case with the aperture port), the
lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0128 proved the coarse-grid (dx = 0.4 mm) lumped-FDTD saturates at ~5 dB
stopband rejection (gate needs 20) — a port-accuracy / discretisation floor, not a
measurement bug. The maintainer chose (AskUserQuestion, 2026-05-31) to **invest in
a finer grid** to chase the original 20 dB gate (gate kept strict).

Crucial enabler: the **aperture port is dx-stable** (ADR-0125 — its reactance no
longer collapses as O(dx²)). So — unlike the *single-cell* port, where ADR-0124's
dx-sweep found finer dx made the response WORSE — finer dx should now genuinely
**sharpen** the resonance (better trace/substrate/tank geometry resolution) and
shrink the residual port-accuracy error. This is the lever.

## Decision

Refine F2.3 (aperture port + CW drive) toward the 20 dB gate, in bounded steps:

1. **dx-refinement decisive test (first increment):** re-run F2.3 (aperture + CW)
   at dx = 0.2 mm (and 0.1 mm if container runtime allows) at the gate frequencies
   (2.0 / 2.4 GHz). Does the 2.4 GHz rejection **climb toward 20 dB** as dx shrinks,
   and does the over-unity passband (the short-board de-embed artifact) resolve?
   - **climbs** → finer dx is the path: find the dx that meets the gate (bounded by
     container runtime; the maintainer sanctioned multi-week compute) → F2.3 merges
     → EM-sim ships.
   - **still capped** → the residual is the port-accuracy + short-board de-embed →
     the next sub-increment is a higher-accuracy aperture port (sub-cell reactance
     correction) + a matched/longer-board de-embed (resolving the over-unity).
2. Use the bounded container; if dx = 0.1 mm is too slow for the full DUT/thru pair,
   run dx = 0.2 mm + extrapolate the trend, and `log` the runtime/cost.

Keep `fdtd_lumped_001`'s strict 20 dB bar (maintainer's choice). Never weaken;
never fake.

## Consequences

**Ships (if finer dx reaches 20 dB):** the goal's EM-simulation component at the
strict cross-validation bar — full-wave FDTD of the lumped board reproducing the
analytic band-pass to ≥ 20 dB. With F2.0/F2.1/F2.2/F2.4 → lumped engine complete.

**Gate:** `fdtd_lumped_001` GREEN at the strict 20 dB tol on the F2.3 branch before
merge; lumped/CPML/aperture gates non-regressed.

**Cost:** finer dx is N⁴; the maintainer sanctioned the multi-week investment. Each
increment runs bounded (a dx step at the gate freqs) and reports the trend +
runtime, so the path/cost stays visible.

**Not in scope (this increment):** the higher-accuracy port (only if dx-refinement
alone is insufficient); tight beyond 20 dB; the studio UI (Track B).

---

## References
- ADR-0128 (the coarse-grid saturation + the maintainer decision); ADR-0125 (the
  dx-stable aperture port); ADR-0124 (the single-cell dx-sweep — finer-was-worse,
  now fixed by the aperture port).
- `docs/superpowers/specs/2026-05-31-f2-3-e-finer-grid-emsim-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-e-finer-grid-emsim.md`.
