# ADR-0132: Filter Phase F2.3-g — PEC-box 2-point standing-wave CW de-embed

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0131 (F2.3-f — matched-CPML failed: monotone + over-unity, hit the
ADR-0108 CPML-into-substrate instability), ADR-0108 (`run_line_eeff` uses PEC +
forward/standing-wave on a PEC line, NOT CPML, for this microstrip geometry),
ADR-0125/0127 (the aperture port is correct in isolation), ADR-0115 (the gate),
the lumped-LC → PCB goal (maintainer chose "keep investing — non-CPML de-embed"),
[[project-lumped-lc-and-studio-redesign]]

---

## Context

Across every de-embed tried — short-board DUT/thru (F2.3-c/d), finer grid
(F2.3-e), matched-CPML (F2.3-f) — the F2.3 board response is **monotone + over-unity,
never a band-pass**, despite the aperture port being proven correct in isolation.
The matched-CPML also hit the **documented** CPML-into-substrate instability
(ADR-0108). The maintainer chose (AskUserQuestion, 2026-05-31) to **keep investing
with a non-CPML de-embed**.

Key insight: **over-unity is a bad-de-embed signature** — the prior measurements
take a single load-cell voltage and never separate the forward (incident/
transmitted) from the backward (reflected) travelling wave, so the standing wave on
the short board corrupts the ratio. A proper **2-point standing-wave de-embed** on a
**PEC box** (stable — `run_line_eeff`'s pattern, no CPML) extracts the true
travelling-wave amplitudes and may reveal the band-pass the bad de-embeds masked.

## Decision

De-embed F2.3's S21 with a **PEC-box 2-point standing-wave CW** method:

- **PEC box** (microstrip in a PEC-bounded box, NO CPML — stable, per ADR-0108),
  long enough that the elements clear the ends and a standing wave develops.
- Drive a Hann-ramped **CW** sinusoid at the input, settle to steady state.
- At each port, sample the line voltage at **≥2 points** (known spacing). The
  standing wave `V(x) = a·e^{−jβx} + b·e^{+jβx}` → solve the 2-point system for the
  forward `a` and backward `b` travelling-wave amplitudes (β from a thru-line ε_eff
  calibration or a 3-point fit). `S21 = b₂/a₁` (transmitted-forward at port 2 over
  incident-forward at port 1), thru-normalized.
- Verify the result is **physical** (no over-unity) and re-assess the shape: does a
  **band-pass** now form (peak at 2.0 GHz, notch at 2.4 GHz)?

**Outcome gate:**
- band-pass forms + notch ≥ 20 dB + passband ≈ 0 dB → EM-sim **ships** (merge F2.3).
- physical band-pass but the notch is shallow (< 20 dB) → the residual is now
  cleanly the **aperture-port accuracy** → the sub-cell reactance correction (the
  "higher-accuracy port" half of the maintainer's choice) is the next sub-increment.
- still monotone / no band-pass even with a clean 2-point de-embed → the board
  integration genuinely doesn't resonate (a deeper finding to surface).
- Keep `fdtd_lumped_001`'s strict 20 dB bar. Never weaken; never fake.

## Consequences

**Ships (if a clean de-embed reveals a ≥20 dB band-pass):** the goal's EM-sim
component at the strict gate → lumped-LC 6/6.

**Gate:** `fdtd_lumped_001` GREEN at 20 dB before merge; lumped/CPML/aperture gates
non-regressed.

**De-risks:** isolates whether the persistent monotone/over-unity is a de-embed
artifact (fixed here) vs a real board-integration failure vs port accuracy — three
distinct, actionable outcomes.

**Not in scope (this increment):** the sub-cell port correction (next, only if the
clean de-embed shows a shallow band-pass); the studio Verify stage.

---

## References
- ADR-0131 (matched-CPML failed); ADR-0108 (`run_line_eeff` PEC + standing-wave, no
  CPML); ADR-0125/0127 (port correct in isolation); ADR-0115 (the gate).
- `docs/superpowers/specs/2026-05-31-f2-3-g-pec-2point-deembed-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-g-pec-2point-deembed.md`.
