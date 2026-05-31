# ADR-0136: App.2.0 — Guided technique-recommender (dual-UI entry)

**Status:** Accepted
**Date:** 2026-05-31
**Related:** the product vision (`docs/superpowers/specs/2026-05-31-ideal-filter-design-app-vision.md`),
ADR-0130 (the Dioxus studio), ADR-0135 (studio live on Pages),
[[project-filter-design-final-goal]], [[lumped-lc-and-studio-redesign]]

---

## Context

A product pass (maintainer brief: "think what the ideal filter-design web app would
look like, browse competitors, do the product work") browsed the landscape: Tier-A
commercial suites (Nuhertz FilterSolutions, Keysight Genesys, AWR iFilter) offer a
**guided** entry (Nuhertz "FilterQuick") that *recommends* a topology from a
requirement; every Tier-B free web calculator (Marki, rftools.io, …) is an **expert**
tool — you must already know the topology and just turn a crank. Yee's studio Technique
stage is likewise expert-only (a gallery of cards). The most product-distinctive gap
is a **guided "recommend-a-technique" entry**, and the maintainer chose it as the next
increment.

## Decision

Add a **guided technique-recommender** as a dual-UI entry on the Technique stage,
split into a validatable pure-domain engine + a thin studio consumer:

- **Engine (`yee-filter`):** `recommend_technique(&FilterSpec) -> TechniqueRecommendation`
  + a `RealizationTechnique` enum (LumpedLc, EdgeCoupled, Hairpin, Combline,
  Interdigital, SteppedImpedance). A deterministic decision tree keyed on response /
  centre-or-cutoff frequency / fractional bandwidth, returning a primary technique, a
  plain-language **rationale** naming the deciding factor, and **ranked alternatives**.
  The thresholds (≈500 MHz distributed-feasibility floor; FBW bands 5% / 20%) are
  documented engineering judgment (Pozar Ch. 8; Hong & Lancaster; Matthaei-Young-Jones)
  **pinned by a gate** so they cannot silently drift.
- **UI (`yee-studio-web`):** a Guided panel atop the expert gallery — a small form →
  "Recommend" → the highlighted primary + rationale + alternatives, routing **live**
  techniques (edge-coupled, lumped) into the flow and honestly labeling **Soon** ones
  (with the nearest live alternative offered). The expert gallery stays = dual-UI.

## Consequences

**Ships:** the studio gains a novice entry — "I want a band-pass at 2.4 GHz, 5% wide"
→ "edge-coupled, because …" → into the flow. Matches the Nuhertz FilterQuick pattern;
differentiates Yee from every free calculator. Pure validatable logic, WASM-safe.

**Gate (non-vacuous):** `cargo test -p yee-filter` asserts canonical cases —
(BP,100 MHz,5%)→LumpedLc, (BP,2.4 GHz,5%)→EdgeCoupled, (BP,2.4 GHz,25%)→EdgeCoupled,
(BP,5 GHz,2%)→Interdigital, (LP,1 GHz)→SteppedImpedance, (LP,50 MHz)→LumpedLc,
(HP,1 GHz)→LumpedLc — plus non-empty rationale and primary-∉-alternatives. A constant
recommender fails. Studio gate: `dx build --platform web --release` EXIT 0 with the
Guided panel rendering.

**Honest scope:** the recommender can recommend techniques not yet built (the four
"Soon" gallery cards); the UI labels those honestly and offers the nearest live
alternative. Building those techniques is the *breadth* track (vision §5), separate
increments. The EM-verify wall (ADR-0133) is untouched — out of scope here.

**Not in scope:** building the Soon techniques; learned/ML recommendation; auto-order
tuning. Deterministic decision tree only.

---

## References
- The product vision doc (§4 gap analysis, §5 the chosen increment).
- `crates/yee-filter/src/lib.rs` (FilterSpec/Response/Topology); `crates/yee-studio-web/src/stages.rs` (technique_stage gallery).
- `docs/superpowers/specs/2026-05-31-app-2-0-guided-technique-recommender-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-0-guided-technique-recommender.md`.
- Pozar, *Microwave Engineering* Ch. 8; Hong & Lancaster, *Microstrip Filters for
  RF/Microwave Applications*; Matthaei, Young & Jones.
