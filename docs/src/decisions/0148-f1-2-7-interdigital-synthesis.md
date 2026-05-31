# ADR-0148: Filter F1.2.7 — interdigital dimensional synthesis

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0144 (F1.2.5 combline synthesis — the pattern this mirrors), ADR-0145
(F1.2.6 combline layout — the F1.2.8 interdigital-layout follow-on mirrors it), ADR-0138
(hairpin λg/4 arm — the length formula), [[lumped-lc-and-studio-redesign]]

---

## Context

Interdigital is the last greyed studio gallery technique. Like edge-coupled / hairpin /
combline it is a coupled-resonator **band-pass**, so it shares the coupling-matrix synthesis
(Qe/Mij from the lowpass prototype + FBW) and the inter-resonator gap solve. What is
interdigital-specific is the resonator: a **straight λg/4 line short-circuited at one end,
adjacent resonators grounded at alternating ends**, with **no loading capacitor**.

## Decision

Ship the interdigital **engine** (`dimension_interdigital` in `yee-filter`), gated, as the
first of the usual three steps (engine → layout → studio-lighting), mirroring combline.

**Key realization — interdigital is the θ = π/2 limit of combline.** A short-circuited stub
of electrical length θ has input susceptance `B(f) = −(1/Z0)·cot(θ·f/f0)`. Combline shortens
the line to θ0 < π/2 and adds `C_L = cot(θ0)/(2π·f0·Z0)` to reach `B(f0)=0`. Interdigital
takes **θ = π/2** (full λg/4): `cot(π/2)=0` → `B(f0)=0` with **no cap**. `dimension_combline`
deliberately errors at θ0 = π/2, so `dimension_interdigital` is a distinct function — same
`solve_gap` coupling + synthesis, different resonator (λg/4, no cap).

`dimension_interdigital(project, substrate) -> Result<InterdigitalDimensions, DimError>`
(no θ0 parameter — θ fixed): spec-`Z0` Hammerstad-Jensen width; `resonator_length_m =
(π/2)/β(f0) = λg/4`; inter-resonator gaps via the shared `solve_gap` (`target_k = FBW·m`).
`InterdigitalDimensions { line_width_m, resonator_length_m, gaps_m, target_k }` — the combline
struct minus `loading_cap_f`/`theta0_rad`.

## Consequences

**Ships:** the interdigital dimensional engine, gated by `dim_interdigital_001` (three parts):
(1) the **published H&L Qe/M benchmark** (5-pole 0.1 dB Cheb, FBW 0.10/0.15 → Qe/M₁₂/M₂₃ vs
H&L's published 11.468/0.07975/0.06077 and 7.645/0.11962/0.09115, < 1% and < 1e-3) — the
non-tautological synthesis core, compared to the book not the synthesizer; (2) the
**interdigital-distinct λg/4 resonance** — independent short-circuited-stub susceptance
`B(f) = −(1/Z0)·cot((π/2)·f/f0)` root-finds to f0 with **no cap** (vs combline needing C_L);
(3) dims solved/bracketed (no clamping), `resonator_length_m == (π/2)/β(f0)`, all positive,
**no loading-cap field**, error paths. Completes the coupled-resonator engine family.

**Gate honesty:** the H&L Qe/M (gate 1) is legitimately shared with combline (the coupling
matrix is prototype-derived, technique-independent) and is non-tautological (vs the book);
gate 2 is the interdigital-specific physics (λg/4 resonance, no cap), distinguishing it from
combline. Not a self-consistency tautology.

**Not in scope:** the grounded-alternating comb board (`dimension_interdigital_layout`,
F1.2.8) and studio lighting (later App increment) — the next two steps. The alternating-ground
even/odd coupling refinement (deferred EM follow-on; first-order reuses the shared coupled-
microstrip coupling, as combline did around the cap). Tap/Qe→feed (F1.2.1).

---

## References
- `crates/yee-filter/src/dimension.rs` (`dimension_interdigital` / `InterdigitalDimensions`,
  mirroring `dimension_combline`); `crates/yee-filter/tests/dim_interdigital_001.rs`;
  `crates/yee-layout` (`microstrip_width`, `eps_eff`, `coupling_coefficient`).
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, §5 (interdigital
  λg/4 short-circuited-at-alternating-ends resonators; the §5.2.5 worked Chebyshev example).
- `docs/superpowers/specs/2026-05-31-f1-2-7-interdigital-synthesis-design.md`;
  `docs/superpowers/plans/2026-05-31-f1-2-7-interdigital-synthesis.md`.
