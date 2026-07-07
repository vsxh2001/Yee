# ADR-0191: A.1 inset-fed patch — machinery shipped, the match is measurably a model gap

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0190 (A.0), ADR-0189 (the directional fit reused here), ADR-0192 (A.2
per-face boundary, whose first measurement happened here).
**Spec:** `docs/superpowers/specs/2026-07-06-a1-a3-antenna-path-design.md`

## What shipped

1. **`yee_layout::inset_fed_patch`** — closed-form matched-feed synthesis: Balanis
   slot conductance `G₁` (piecewise in `W/λ₀`; the mutual term `G₁₂` neglected,
   documented) → `R_edge = 1/(2G₁)` → inset depth `x₀ = (L/π)·acos(√(Z₀/R_edge))`;
   geometry as a 4-rect metal union (feed through the notch). Plus
   **`inset_fed_patch_with_depth`** — the explicit-depth variant that is the A.3
   design-loop knob. Unit-gated vs hand arithmetic (G₁ = 1.03 mS, x₀ = 11.4 mm).
2. **`sparams::directional_reflection_db`** — single-run |Γ| = |bwd|/|fwd| from the
   three-probe standing-wave fit (the slotted-line measurement): no reference run, no
   subtraction artifacts, half the solve cost of the A.0 method. Unit-gated (a
   synthetic Γ = 0.3 reads −10.46 dB).
3. Gate `engine-antenna-002`: dip position ±10 % of the designed 2.45 GHz + ≥ 1 dB
   depth tripwire.

## Measured — and what it teaches

| variant | dip | depth |
|---|---|---|
| PEC lid (A.0 boundary) | 2.475 GHz (+1.0 %) | −1.1 dB |
| **Open top (A.2 per-face CPML)** | **2.425 GHz (−1.0 %)** | **−1.2 dB** |

The resonance is exactly where the closed forms put it, under both boundaries — but the
**match is poor and the lid was not the cause**. Working backward from |Γ| = 0.87, the
input resistance at the inset is ~5 Ω: the G₁-only model's `R_edge = 485 Ω` is far too
high for this thick / high-ε_r substrate, so the depth formula pushed the feed toward
the patch centre where `R → 0`. That is a *closed-form model gap*, precisely the class
of error the project's EM-in-the-loop machinery exists to close — **A.3 tunes
`x_inset` against the measured return loss** instead of trusting `G₁`.

## Consequences

The single-run directional |S11| becomes the antenna track's standard observable. The
A.1 gate honestly certifies resonance placement and coupling; the matched-antenna
assert (≤ −10 dB) moves to A.3 where the loop earns it.
