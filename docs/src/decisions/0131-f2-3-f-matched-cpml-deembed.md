# ADR-0131: Filter Phase F2.3-f — matched-CPML board de-embed for the lumped EM sim

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0129 (F2.3-e — finer dx DISPROVEN; the residual is the short-board
DUT/thru de-embed, not grid), ADR-0122 (per-axis CPML — enables an x-matched board),
ADR-0123 (the matched-line CPML≠matched-at-DC dead end — sidestepped here by CW),
ADR-0125 (the aperture port), ADR-0115 (the gate), the lumped-LC → PCB goal
(maintainer chose "higher-accuracy port + matched de-embed"),
[[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0129 disproved finer-grid: the F2.3 notch collapses (5.1→0.19 dB at dx
0.4→0.2 mm) and the passband flips over-unity→loss — the **short-board DUT/thru
de-embed does not converge**, and that (not grid resolution) is the dominant wall.
The aperture port is proven correct *in isolation* (ADR-0125/0127). The maintainer
chose (AskUserQuestion, 2026-05-31) the **higher-accuracy port + matched de-embed**
path. The de-embed is the dominant symptom and a prerequisite to cleanly measure
the port accuracy — attack it first.

**Why it's now tractable:** per-axis CPML (ADR-0122, `with_axes`) lets the F2.3
microstrip be terminated with **x-only CPML matched ends** (a matched output that
absorbs the transmitted wave) while keeping PEC transverse walls. Under the **CW**
drive (ADR-0128, single-frequency — NOT the broadband/DC regime where ADR-0123
found CPML≠matched-at-DC), the CPML absorbs cleanly at the carrier, so the output
reference plane reads a **clean transmitted wave**, not a standing-wave artifact.

## Decision

Replace F2.3's short-board lumped-resistor DUT/thru measurement with a
**matched-CPML board de-embed**:

- Terminate the microstrip with **x-only CPML** at the output (and input) end
  (`CpmlParams::for_grid(..).with_axes([true,false,false])`), PEC transverse walls,
  so the transmitted (and source-side reflected) waves are absorbed after one pass.
- Under the CW per-frequency drive, measure the **transmitted-wave amplitude** at an
  output reference plane (de-embedded; with the matched end there is no backward
  reflection to corrupt it), thru-normalized: `S21(f) = |V_out,ss(f)| /
  |V_thru,ss(f)|`. This should yield a **physical** |S21| (no over-unity) that is
  **dx-stable** (the de-embed converges).
- Lengthen the board if needed so the reference planes sit clear of the
  element/port discontinuities.

**Outcome gate:**
- |S21| now physical + dx-stable AND the notch reaches the strict 20 dB → EM-sim
  **ships** (merge F2.3).
- physical + dx-stable but the notch is still shallow → the residual is now cleanly
  the **aperture-port accuracy** → the next sub-increment is the sub-cell reactance
  correction (the "higher-accuracy port" half of the maintainer's choice).
- Keep `fdtd_lumped_001`'s strict 20 dB bar. Never weaken; never fake.

## Consequences

**Ships (if the matched de-embed + the proven port reach 20 dB):** the goal's
EM-sim component at the strict gate → lumped-LC 6/6.

**Gate:** `fdtd_lumped_001` GREEN at 20 dB on the F2.3 branch before merge; the
lumped/CPML/aperture gates non-regressed.

**De-risks the maintainer's path:** isolates whether the residual is the de-embed
(fixed here) or the port accuracy (the sub-cell correction next), instead of
conflating them.

**Not in scope (this increment):** the sub-cell port reactance correction (next,
only if the clean de-embed still shows a shallow notch); the studio Verify stage.

---

## References
- ADR-0129 (finer-dx disproven; de-embed is the wall); ADR-0122 (per-axis CPML);
  ADR-0123 (CPML≠matched-at-DC — sidestepped by CW); ADR-0125/0127 (the port is
  correct in isolation).
- `docs/superpowers/specs/2026-05-31-f2-3-f-matched-cpml-deembed-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-f-matched-cpml-deembed.md`.
