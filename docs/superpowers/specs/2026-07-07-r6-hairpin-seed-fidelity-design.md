# R.6 — hairpin seed fidelity: corner correction + resonator impedance

**Date:** 2026-07-07
**Track:** RF-TOOL-ROADMAP R.6 (feeds R.4c)
**Related:** ADR-0197 (R.4 — the measured +17 % detune and the
tap-realizability wall this addresses), ADR-0109 (F1.2.2 hairpin dims).

## Problem — two quantified seed defects from the R.4 instrumented runs

1. **Corner detune**: the fold-corrected midline model still measured the
   seed resonance ~+17 % high (5.95 GHz vs designed 5.0 on h = 0.8 mm FR-4).
   The U's two 90° corners cut the current path relative to the midline —
   the classic bend-length effect the midline model cannot see.
2. **Tap-realizability wall on thick stacks**: at 50 Ω resonator lines on
   h = 1.6 mm FR-4, the fat (3 mm) line's fold consumes the half-wave and
   `TapNotRealizable` correctly rejects. The classic fix is **thinner,
   higher-impedance resonator lines** (Zr ≈ 70 Ω) — the tap formula already
   carries `(Z0/Zr)`; the dims layer doesn't.

## Scope

`yee_filter::HairpinOptions { fold_widths, resonator_z_ohm, corner_widths }`
+ `dimension_hairpin_opts`:

- **Corner correction** (default ON, `corner_widths = 0.85`): each corner
  shortens the electrical length by ≈ `κ·w`, so
  `arm = (λ_g/2 − fold)/2 + κ·w`. κ = 0.85 is **calibrated from the single
  R.4 measured data point** (resonance 5.95 vs 5.0 GHz → 2.62 mm deficit
  over two corners of a 1.529 mm line, net of open-end lengthening) and
  documented as such — refined when more stacks are measured, and
  independently checked by the R.6 gate run.
- **Resonator impedance** (`resonator_z_ohm`, default `None` → spec Z0):
  resonator width/ε_eff/λ_g/fold/gaps all computed at Zr; the tap solve uses
  `(Z0/Zr)`; the feed stays a Z0 line — `HairpinDimensions` gains
  `feed_width_m` so layouts stop conflating the two widths.

`dimension_hairpin`/`_with_fold` delegate with defaults; the corner fix
(being a fix, like the fold correction before it) applies to them too —
`hairpin_dim_001` evolves, and dependent scenarios re-measure.

## Gates

- Unit: corrected arm formula pinned; Zr = 70 Ω on the previously-rejected
  h = 1.6 mm stack now dimensions (tap on the arm, thinner resonator line,
  Z0-width feed); defaults reproduce `dimension_hairpin_with_fold`'s widths.
- Full-wave: re-run `engine-bpf-bo-001` — the seed's passband location is
  the direct check of κ. Assert numbers evolve to whatever the corrected
  seed honestly supports (if a real passband forms at dx = 0.2 mm, the
  original centre/peak asserts return; if the coupling floor still buries
  it, the machinery asserts stay and the ADR records the new seed numbers).
