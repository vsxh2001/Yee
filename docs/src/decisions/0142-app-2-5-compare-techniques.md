# ADR-0142: App.2.5 — Compare techniques side-by-side

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0136 (App.2.0 recommender — recommend one), ADR-0141 (App.2.4
`verify_view` — the metric-extraction pattern this mirrors), ADR-0138/0139 (the
hairpin + low-pass flows being compared), the product vision (§1 dual entry, §5 P2),
[[lumped-lc-and-studio-redesign]]

---

## Context

At the complete-app boundary, the maintainer chose "deepen the flows — optimize/compare"
as the next direction. The Technique stage already offers a guided recommendation
(App.2.0) and an expert gallery, but no **side-by-side comparison**: to weigh
edge-coupled vs hairpin vs lumped for a given spec, the user must select each, walk its
flow, and remember the board size / verdict. The per-flow engines + the App.2.4
metric-extraction already provide everything needed to compare.

## Decision

Add a **Compare** view (Technique stage) built on a pure host-testable helper:

- `compare_techniques(spec) -> Vec<TechniqueComparison>` synthesizes every **live**
  technique that realizes the spec's response class (band-pass → edge-coupled / hairpin /
  lumped; low-pass → stepped-impedance; high-pass → none yet) and collects a comparable
  row per technique — board size, PASS/FAIL, order, worst ripple / return loss / stopband
  rejection — pulled directly from each design's existing graded structs (real engine
  output; `realizable=false` when a design fails to dimension).
- A Compare table on the Technique stage marks the recommended technique
  (`recommend_technique`) and offers a "Use this" routing into each realizable row,
  completing the dual entry (recommend → compare → pick).

## Consequences

**Ships:** the user can see, for their spec, how every applicable technique compares —
board size, verdict, key metrics — side by side, with the recommendation marked and
one-click routing. Real engine output for every cell; built entirely on the existing
per-flow engines + the App.2.4 metric extraction. No new physics.

**Gate (non-vacuous):** `compare_techniques` on a band-pass spec returns the three
band-pass techniques with metrics equal to each design's graded fields and rows that are
not all identical (e.g. hairpin vs edge-coupled board sizes differ); a low-pass spec
returns exactly the stepped-impedance row; a high-pass spec returns `[]`. A
constant/empty helper fails. Plus `dx build` EXIT 0 + no regression.

**Not in scope:** a multi-technique response-overlay plot (follow-on); tune-and-watch
sliders (the spec form already re-derives live); cross-response comparison. The EM wall
(ADR-0133) is untouched.

---

## References
- `crates/yee-studio-web/src/{engine.rs (design_*_from, the graded structs, verify_view),
  stages.rs (technique_stage)}`; `crates/yee-filter/src/recommend.rs`.
- `docs/superpowers/specs/2026-05-31-app-2-5-compare-techniques-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-5-compare-techniques.md`.
