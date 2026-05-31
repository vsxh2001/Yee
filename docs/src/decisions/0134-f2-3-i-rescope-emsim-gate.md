# ADR-0134: Filter Phase F2.3-i — re-scope the lumped EM-sim gate to a physically-achievable bar

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0115 (the original `fdtd_lumped_001` ≥20 dB band-pass gate),
ADR-0133 (F2.3-h — the de-embed avenue exhausted; CW high-Q microstrip in a stable
box is cavity-dominated, a fundamental measurement wall), ADR-0125/0127 (the
aperture port + `cap_cw_001` validate the per-element reactance in isolation),
ADR-0111 (`ladder_s21` validates the sharp circuit response vs Pozar), the
lumped-LC → PCB goal (maintainer chose "re-scope the gate + ship", AskUserQuestion
2026-05-31), [[project-lumped-lc-and-studio-redesign]]

---

## Context

EM-sim (F2.3) has had ~15 increments. The aperture lumped port is **proven correct
in isolation** (`aperture_port_001`, `cap_cw_001` — the per-element R/L/C reactance)
and the full-wave FDTD board **loads the line** (the components couple, a real
frequency-dependent response). But the original `fdtd_lumped_001` bar — the FDTD
**full-board** |S21| reproducing the analytic band-pass to **≥20 dB** — is a
**fundamental FDTD-measurement wall** (ADR-0133): a high-Q microstrip filter's CW
steady-state S21 in any **stable** (PEC) box is **cavity-dominated** (the lossless
thru itself reads 6× over-unity at a box mode), and the only matched termination
that kills the cavity (CPML) is unstable into the substrate (ADR-0108). It is **not
a defect** — the physics is validated; the *full-board cross-validation measurement*
is what's intractable. The maintainer chose to **re-scope the gate to a
physically-achievable bar** and ship.

## Decision

Re-scope `fdtd_lumped_001` from "FDTD full-board |S21| ≥20 dB band-pass" to a
**principled, achievable, real EM-integration validation** of what the full-wave
sim genuinely demonstrates:

- the EM-sim **pipeline runs end-to-end** (synthesize_lumped → lumped_board →
  voxelize_microstrip → aperture-port placement → FDTD solve → S21 sweep) and
  produces a **finite, non-trivial** result (no NaN/Inf);
- the lumped components **demonstrably LOAD the line** — the loaded board's
  response is **meaningfully frequency-dependent / differs from the bare thru** by
  a real margin (the elements are NOT inert; this is the genuine full-wave EM
  contribution). The agent picks the **specific metric that is RELIABLY true** of
  the actual F2.3 board data and **meaningful** (would FAIL for an inert/broken
  sim) — NOT vacuous.

The gate's docstring **delegates** the sharp-response cross-validation to the
circuit `ladder_s21` (F2.0, vs Pozar) and the per-element reactance to
`aperture_port_001`/`cap_cw_001`, and **documents the cavity wall** (why the
full-board ≥20 dB FDTD cross-validation is a fundamental measurement limitation,
ADR-0133). This is an **honest re-scope to the achievable bar, NOT a weakening to
fake a pass** — the maintainer's explicit choice; the asserted property is real.

**Integrity guardrail (non-negotiable):** the re-scoped assertion must be a REAL,
MEANINGFUL EM property that would **FAIL** if the elements were inert (the
single-cell-placement flat ≈1 response) or the sim broke. A code-reviewer (never
self-review) must confirm it is not a tautology / always-pass. If a principled
meaningful achievable assertion cannot be found, do NOT ship a vacuous gate —
surface that instead.

## Consequences

**Ships:** the goal's **EM-simulation** component at an honest, achievable bar —
full-wave FDTD of the lumped-LC board, components loading the line via a
per-component-validated port, the sharp response cross-validated at the circuit
level. F2.3 merges → **lumped-LC goal 6/6**. The studio's Verify stage can surface
this EM-sim result.

**Gate:** the re-scoped `fdtd_lumped_001` GREEN at the achievable bar (reviewer-
confirmed non-vacuous); lumped/CPML/aperture/port gates non-regressed; the
≥20 dB-wall + the delegation documented in the gate + this ADR.

**Not in scope:** chasing the ≥20 dB full-board cross-validation (the documented
wall; a future stable-non-CPML-absorber research track if ever revisited); the
studio Verify stage wiring (follow-on).

---

## References
- ADR-0133 (the cavity wall + the exhausted de-embed avenue); ADR-0115 (the
  original gate); ADR-0125/0127 (`aperture_port_001`/`cap_cw_001` validate the
  port in isolation); ADR-0111 (`ladder_s21` validates the sharp response).
- `docs/superpowers/specs/2026-05-31-f2-3-i-rescope-emsim-gate-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-i-rescope-emsim-gate.md`.
