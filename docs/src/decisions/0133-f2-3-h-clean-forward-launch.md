# ADR-0133: Filter Phase F2.3-h — clean forward-wave launch for the lumped EM sim

**Status:** Investigated — the launch fix RESOLVED `a₁` (time-gated incident,
β>0, well-resolved), but `b₂` is **cavity-degenerate**: a PEC box under CW steady
state is an intrinsic high-Q cavity whose own resonances dominate `b₂` + the thru
(thru over-unity 6.06× at a box mode). The de-embed avenue is **exhausted**; the
residual is a **fundamental FDTD-measurement wall** (high-Q-filter-CW-steady-state
⊥ stable-cavity-free-termination). Re-surfaced for a maintainer decision. Branch
`23a52e0` (unmerged). See Outcome.
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

## Outcome (2026-05-31) — `a₁` resolved; `b₂` cavity-degenerate → fundamental wall

Built (branch `23a52e0`): `calibrate_launch` (time-gated Gaussian-pulse pre-pass
per freq, `run_line_eeff` pattern → a pure forward incident `a₁_gated` + numerical
β), `inject_directional_source` (two `E_z` sheets one cell apart, downstream
retarded by `β·dx` — poor-man's Huygens/TF-SF, forward-biased), CW `b₂` referenced
to `a₁_gated`.

- **`a₁` RESOLVED** (the deliverable that succeeded): β = 69–89 rad/m (positive,
  ~15% of the FR-4 quasi-TEM guess at every gate freq), `|a₁_gated|` ~0.49–0.66
  (above floor) — the F2.3-g `a₁` ambiguity is fixed.
- **`b₂` STILL cavity-floor-limited:** in CW steady state the PEC box is an
  **intrinsic high-Q cavity** — the input stays a near-pure standing wave
  (refl/fwd ≈ 0.93–1.00 in DUT *and* thru; the directional source fixes the
  transient, not the steady state), β_out=0 at 1.6/1.8 GHz, and the **bare thru
  `b₂` is dominated by box cavity resonances** (`|b₂/a₁|_thru` swings 0.05 → 0.45 →
  **6.06 (2.2 GHz, over-unity on a lossless thru)** → 0.18). |S21| sweep:
  1.6→31, 1.8→34, 2.0→56, 2.2→78, 2.4→43, 2.6→29 dB — NOT a band-pass; the deepest
  point is at **2.2 GHz, tracking a box cavity mode, not f0**.

**Classification: a fundamental FDTD-measurement wall.** A CW steady-state S21 of a
high-Q microstrip filter in a PEC box is intrinsically a **cavity** measurement —
the box's resonances dominate `b₂` + the thru over both the filter response and the
launch directionality. This is the **third outcome, definitively**: not a cheap
topology fix. `fdtd_lumped_001` RED (56 dB in-band IL ≫ 6 dB bar), **NOT weakened**.

**The de-embed avenue is exhausted** (short-board → over-unity; finer-grid →
collapse; matched-CPML → unstable; PEC 2-point → physical but launch-floor; clean
launch → `a₁` fixed but `b₂` cavity-bound). The **fundamental tension**: a high-Q
filter S21 needs **CW steady state** (the tanks must ring up) → in any **stable**
(PEC) box that is **cavity-dominated** → and the only matched termination that
kills the cavity is **CPML, which is unstable into the substrate** (ADR-0108/0131).
Every approach hits one horn. The aperture port is proven correct **in isolation**
(6.9/6.10); the **circuit `ladder_s21`** already validates the sharp response
(F2.0). 

**→ Maintainer decision (AskUserQuestion, 2026-05-31), after ~15 increments + the
de-embed avenue exhausted:** (a) the one remaining technique-class — a **stable
non-CPML absorbing termination** (a graded lossy-material / tapered-resistive
microstrip terminator: pure loss = unconditionally stable, absorbs the line wave =
kills the cavity, matched = no over-unity) — uncertain multi-week; (b) **re-scope**
`fdtd_lumped_001` to a physically-achievable bar (the FDTD board loads the line +
the circuit `ladder_s21` validates the sharp response); (c) **accept / defer**
EM-sim (goal 5/6 with the UI shipped; documented FDTD-measurement limitation).

---

## References
- ADR-0132 (the physical de-embed + the launch-floor limiter); ADR-0108
  (`run_line_eeff` time-gated incident-wave, PEC; CPML-into-substrate unstable);
  ADR-0131 (matched-CPML unstable); ADR-0014/0021/0026 (TF/SF);
  ADR-0125/0127 (port correct in isolation); ADR-0111 (`ladder_s21`); ADR-0115.
- `docs/superpowers/specs/2026-05-31-f2-3-h-clean-forward-launch-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-h-clean-forward-launch.md`.
