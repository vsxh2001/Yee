# ADR-0200: R.2b — measured complex Γ of a known one-port, differentially

**Status:** Accepted
**Date:** 2026-07-07
**Related:** ADR-0195 (R.2 — the THRU-Γ negative result this answers),
ADR-0194 (R.1 vias), ADR-0189 (the 3-probe split + its residual flag).
**Spec:** `docs/superpowers/specs/2026-07-07-r2b-shorted-line-gamma-design.md`

## Decision

Gate `engine-sparams-003` (`crates/yee-engine/tests/board_short_gamma.rs`,
two release solves) measures the complex reflection of a via-shorted
microstrip line — Γ = −e^{−2jβd} in TL theory — **differentially over two
short distances**, asserting the phase-slope *difference* against
`−4π·Δd·√ε_eff/c` so the termination's own reactance cancels exactly.

**Measured: |Γ| mean 0.998 / 0.995 (±15 % gate); differential round-trip
phase slope error 2.32 % (±5 % gate); single-run slope ratios 1.143 / 1.082
vs the ideal short.**

## The instrumented iterations that forced each design element

Four runs, each converting a wrong number into a scenario fix (all preserved
in the gate's module docs):

1. **Single centre-cell via**: |Γ| ≈ 0.98 but a 3× phase slope — one Ez
   column under a 10-cell-wide trace is a partial shunt, and the passed wave
   reflected off the open line end (composite reflector). → full-width
   **via fence**.
2. **Fence with line continuing beyond**: still ~1.8× — the fence leaks. →
   the trace now **ends at the fence plane** (one reflection plane).
3. **Single-plane, single-distance**: +14 % slope excess — the fence is an
   *inductive* short (~0.3 nH at these dimensions), and an inductive short's
   dφ/df adds apparent depth. Real termination physics, not an error. → the
   assert went **differential**, which cancels it; the single-run ratios
   stay as a bounded sanity check ([0.95, 1.35]), and the measured excess
   behaves exactly as a fixed reactance should (constant ~1.7–1.9 mm
   absolute, shrinking relative share with distance).
4. **Quality-bin selection**: with a total reflector the standing-wave nulls
   sweep across the probe triple, and the 3-probe fit degenerates where a
   null sits near the middle probe (measured |Γ| → 0 fallback bins). Bins
   are selected by the fit's own flags — the ADR-0189 residual |Im cos βd|
   and the fitted β against Hammerstad–Jensen — 41/41 pass in the final
   scenario.

## What this closes and what stays open

Closed: the engine demonstrably measures a **complex Γ** — magnitude and
phase delay — against TL theory on a known one-port; the ADR-0195 question
("can complex_reflection be trusted?") is answered yes, with the fit-quality
selectors as the usage contract. Open (unchanged): fixture de-embedding for
a *small* THRU S11 — a calibration-standards problem deferred until a use
case forces it. The gate runs under the blanket `yee-engine gates` CI step
(include-ignored) automatically.
