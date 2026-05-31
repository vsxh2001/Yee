# ADR-0133: Filter Phase F2.3-h — clean forward-wave launch for the lumped EM sim

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0132 (F2.3-g — the 2-point de-embed made S21 PHYSICAL, but the
PEC-box soft source barely launches forward power → launch/probe-floor-limited;
the notch-at-f0 is likely artifact), ADR-0108 (`run_line_eeff` time-gated
incident-wave on a PEC line), ADR-0014/0021/0026 (TF/SF source), ADR-0125/0127
(port correct in isolation), ADR-0115 (the gate), the lumped-LC → PCB goal
(maintainer chose "keep investing"), [[project-lumped-lc-and-studio-redesign]]

---

## Context

F2.3-g (ADR-0132) achieved the first **physical** F2.3 de-embed (no over-unity) via
2-point forward/backward separation — but exposed a new limiter: the PEC-box soft
`E_z` source **reflects almost entirely** (input ≈ pure standing wave) and the bare
thru **barely couples forward power to the output region** (β_out=0 at 1.6/1.8 GHz,
|b₂|~0.02 vs |a₁|~7–14). So `S21 = (b₂/a₁)_dut/(b₂/a₁)_thru` divides small,
partly-degenerate readings — the "deep notch at f0" is **likely a floor artifact,
not a real result**. The de-embed math is sound; the **forward-wave launch** is now
the wall.

## Decision

Give F2.3 a **clean forward-wave launch + a well-resolved output probe** so the
travelling-wave amplitudes `a₁` (incident at input) and `b₂` (transmitted at
output) are trustworthy (β>0 at all gate freqs, `b₂` well above the floor), then
re-measure S21 and **disambiguate** the notch-at-f0:

- Use the `run_line_eeff` **time-gated incident-wave** pattern (ADR-0108) for the
  forward reference: a pulse launched into a long-enough line gives a clean
  incident `a₁` (time-gated before the first reflection); and/or a directional /
  TF-SF-style launch that injects predominantly forward (less source reflection).
- Lengthen the line + place the output reference region where a propagating forward
  wave is clean (clear of the PEC end wall / evanescent zones — fixing β_out=0).
- Keep the CW steady-state for the DUT *response* (the tanks must ring up), but
  reference it to a trustworthy forward `a₁` (hybrid: time-gated incident `a₁` +
  CW-settled `b₂`, or a high-amplitude directional CW launch).

**Outcome gate (disambiguation):**
- a clean band-pass emerges (peak @2.0 GHz, notch @2.4 GHz ≥20 dB) → EM-sim **ships**.
- a clean **inverted** response (notch @f0) persists with a trustworthy launch →
  a **real topology inversion** (shunt tanks shorting at f0) → a cheap F2.3
  placement fix (next).
- still floor-degenerate / inconclusive → the FDTD S21 of a high-Q microstrip
  filter is a genuine multi-layer measurement-research wall → surface the
  cumulative picture to the maintainer.
- Keep `fdtd_lumped_001`'s strict bar. Never weaken; never fake.

## Consequences

**Ships (if a trustworthy launch reveals a ≥20 dB band-pass):** the goal's EM-sim
component → lumped-LC 6/6. Otherwise it **definitively** classifies the residual
(topology bug → cheap fix; or a real research wall → maintainer decision), instead
of the current floor-ambiguous state.

**Gate:** `fdtd_lumped_001` GREEN at 20 dB before merge; gates non-regressed.

**Not in scope:** the topology-inversion fix (next, if the launch reveals a real
inverted response); the sub-cell port correction; the studio Verify stage.

---

## References
- ADR-0132 (the physical de-embed + the launch-floor limiter); ADR-0108
  (`run_line_eeff` time-gated incident-wave, PEC); ADR-0014/0021/0026 (TF/SF);
  ADR-0125/0127 (port correct in isolation); ADR-0115 (the gate).
- `docs/superpowers/specs/2026-05-31-f2-3-h-clean-forward-launch-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-h-clean-forward-launch.md`.
