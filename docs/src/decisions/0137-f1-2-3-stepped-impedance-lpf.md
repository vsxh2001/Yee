# ADR-0137: Filter F1.2.3 — Stepped-impedance low-pass filter synthesis + dimensions

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0136 (App.2.0 recommender — recommends `SteppedImpedance` for low-pass;
this makes it real), the product vision (`docs/superpowers/specs/2026-05-31-ideal-filter-design-app-vision.md` §5),
the edge-coupled / hairpin dimensional synthesis (`yee-filter::dimension`),
[[project-filter-design-final-goal]], [[lumped-lc-and-studio-redesign]]

---

## Context

Yee's distributed synthesis is band-pass-only (edge-coupled, hairpin). The 2026-05-31
product-vision pass ranked the **stepped-impedance low-pass filter** as the highest
value-per-effort breadth move: the simplest distributed topology (alternating high-Z /
low-Z microstrip line sections), a textbook closed-form design (Pozar §8.6), the first
**low-pass** capability, and the synthesizer the App.2.0 recommender already points at
(it recommends `SteppedImpedance` for low-pass ≥ 500 MHz but nothing backs it).

## Decision

Add closed-form stepped-impedance low-pass synthesis + dimensions to
`yee-filter::dimension`, mirroring the edge-coupled / hairpin pattern:

- From the low-pass prototype g-values (`yee_synth::prototype`), map each reactive
  element `g_k` to a short microstrip line section, alternating shunt-capacitor (low-Z)
  / series-inductor (high-Z) starting with a shunt capacitor: electrical length
  `βl = g_k·Z_low/Z₀` (low-Z) or `βl = g_k·Z₀/Z_high` (high-Z) — Pozar §8.6.
- Physical width via `yee_layout::microstrip_width`; physical length from the guided
  wavelength at the section width (`yee_layout::eps_eff`). Plus a placeable layout
  (`dimension_stepped_impedance_layout`).

**Scope: synthesis + dimensions + gate only.** Lighting the studio's `SteppedImpedance`
gallery card with a live flow is deferred — the studio Spec→Synthesis is band-pass-only,
and a low-pass response path through it is a distinct follow-on increment.

## Consequences

**Ships:** Yee's first **low-pass** distributed capability; the recommender's
`SteppedImpedance` recommendation is now backed by a real, validated synthesizer +
dimensions + layout. Pure closed-form, no new dependency, WASM-safe.

**Gate (`dim_stepped_001`, published-benchmark, non-vacuous):** Pozar Example 8.6 —
Butterworth N=6, f_c=2.5 GHz, Z₀=50 Ω, Z_high=120 Ω, Z_low=20 Ω → the six section
electrical lengths (source→load) match Pozar's table within ±1.0°:
`[11.85°, 33.76°, 44.28°, 46.12°, 32.41°, 12.34°]`. Six distinct values from real
g-values — a constant recommender/synthesizer fails. The test derives βl from the
formula (does not hardcode the computed value), asserts the low-Z-first alternation,
and that physical lengths are positive and finite.

**Not in scope:** studio low-pass UI (follow-on); elliptic stepped-Z; stub low-pass; EM
verification (the cavity wall ADR-0133 is untouched).

---

## References
- Pozar, *Microwave Engineering* §8.6 (Stepped-Impedance Low-Pass Filters), Example 8.6.
- `crates/yee-filter/src/dimension.rs` (edge-coupled / hairpin pattern);
  `crates/yee-synth/src/lib.rs` (`prototype`, `Prototype`).
- `docs/superpowers/specs/2026-05-31-f1-2-3-stepped-impedance-lpf-design.md`;
  `docs/superpowers/plans/2026-05-31-f1-2-3-stepped-impedance-lpf.md`.
