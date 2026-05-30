# ADR-0127: Phase 2.fdtd.6.10 — CW capacitor steady-state diagnostic

**Status:** Accepted — shipped (merge `f053164`). VERDICT: **the aperture cap arm
is CORRECT**; F2.3's failure is a measurement (pulse-drive) limit → F2.3 needs a CW
per-frequency drive (the next increment).
**Date:** 2026-05-31
**Related:** ADR-0126 (F2.3-c — the "longer record → worse" symptom this settles),
ADR-0125 (the aperture port), ADR-0115 (the F2.3 gate), the lumped-LC → PCB goal,
[[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0126 found that F2.3 on the aperture port loads the line but no band-pass
forms, and a **longer record made the in-band loss worse** — the shunt-tank
capacitor looked like a deepening near-short. Two hypotheses: a *measurement*
limit (the pulse DFT of an unsettled transient on a short standing-wave line) vs a
*cap-update windup bug*. This is the deciding diagnostic.

## Decision

Ship a CW single-frequency steady-state diagnostic (`cap_cw_001`, `#[ignore]`'d,
release) that drives the aperture capacitor — and a shunt L‖C tank — with a CW
sinusoid at f0 to steady state, and **asserts** the isolated arm's behaviour
(line-de-embed-noise-free): the cap reactance is `−jX` in every settled window;
`|Z|` does **not** drift over the record (`< 0.10`); `V_C` is a bounded
oscillation (no windup, `< 0.10`); the realized `|Z| ≈ β + 1/(jωC)` analytic
(`< 0.10`); and the L‖C tank **resonates** (`|Z| > single-arm |X|`). An on-guide
field-mediated probe is *recorded* (not asserted) — it scatters, confirming the
short-line de-embed is the noise source.

## Consequences — VERDICT: MEASUREMENT-LIMIT, cap arm correct

Measured under clean CW: `Z_cap = 97.9 − j100.1 Ω` (Im exactly `−1/(ωC)`,
capacitive every window); `|Z|` drift **0.003**; `V_C` envelope growth **−0.0000**
(bounded); realized vs analytic accuracy **0.001**; the L‖C tank **resonates**
(mean `|Z| = 157.8 Ω` > single-arm `100 Ω`; parallel resonance `+149.6 + j49.8 Ω`).

So **the aperture port — inductor (6.9) and capacitor — is proven correct**; there
is no cap-update bug (`lumped.rs` unchanged). F2.3's "longer → worse" is the
**pulse drive + short-line standing-wave de-embed** (a windowed `V_T/I` of a 90°
branch current beats against the carrier). **The remaining EM-sim fix is purely
downstream: F2.3 needs a CW per-frequency steady-state drive** (drive each sweep
frequency as a CW sinusoid to steady state, measure the steady-state DUT/thru S21)
— the next increment (F2.3-d). The physics is settled; the band-pass should form
once F2.3 measures at steady state.

**Gate:** `cap_cw_001` GREEN (genuine, tight assertions); no regression
(`aperture_port_001`, `lumped_rlc_twoway_001`, the lumped/cpml gates green); diff =
the test only.

**Not in scope:** the F2.3 CW driver (F2.3-d, next); a CI gate job for `cap_cw_001`
(it is an on-demand diagnostic like `aperture_port_001`; a CI-lane follow-up if
wanted).

---

## References
- ADR-0126 (the symptom); ADR-0125 (the aperture port + the CW caveat).
- `docs/superpowers/specs/2026-05-31-fdtd-6-10-cw-capacitor-diagnostic-design.md`;
  `docs/superpowers/plans/2026-05-31-fdtd-6-10-cw-capacitor-diagnostic.md`.
