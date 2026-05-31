# ADR-0150: App.2.8 — Light the Interdigital technique in the studio

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0148 (F1.2.7 interdigital engine), ADR-0149 (F1.2.8 interdigital layout — the
two prerequisites), ADR-0146 (App.2.7 combline lighting — the pattern this mirrors), ADR-0136
(the recommender), ADR-0142/0143 (compare/overlay — interdigital joins them),
[[lumped-lc-and-studio-redesign]]

---

## Context

Interdigital's synthesis (F1.2.7) and board layout (F1.2.8) shipped, but the studio gallery
still greys it (`technique_status(Interdigital) = Soon(EdgeCoupled)` stand-in) and there is no
`Topology::Interdigital`. With both engine + layout in place, lighting it is a clean mirror of
App.2.7 combline. Interdigital is a coupled-resonator **band-pass** — same synthesis as
edge-coupled / hairpin / combline; only the realization differs (full λg/4 lines
short-circuited at alternating ends, **no loading cap**).

## Decision

Add `Topology::Interdigital` and light it as a live studio technique (`crates/yee-studio-web`),
mirroring every site combline occupies:

- `derive_geometry` → a `Topology::Interdigital` arm calling `dimension_interdigital` +
  `dimension_interdigital_layout` (no `theta0`), `loading_cap_f: None`; the shared
  synthesis/response/verdict computed once.
- `topology_name` = "interdigital (λ/4, alt. short)", `length_label` = "resonator length (mm)";
  the topbar / verify / layout distributed-flow match groups gain `| Topology::Interdigital`.
- `technique_status(Interdigital)` → `Live(Topology::Interdigital)`;
  `topology_response(Interdigital)` → `Bandpass`; `technique_topology` / `technique_label`
  arms; the gallery card → `selects: Some(Topology::Interdigital)`.
- `compare_techniques` gains an interdigital row (the band-pass list is hardcoded — 5 rows:
  edge-coupled / hairpin / combline / interdigital / lumped); `overlay_curves` joins it to the
  shared coupled-resonator ideal curve (interdigital differs only physically, like combline).
- The combline-distinct loading cap is simply **absent** (`combline_loading_cap_f = None`); the
  resonator table surfaces the λg/4 length instead.

## Consequences

**Ships:** interdigital is a live, routable studio technique end-to-end (synthesis → board →
Gerber/KiCad, in compare + overlay), honestly surfacing the λg/4 alternating-short resonator
with no cap. **Completes the coupled-resonator gallery:** edge-coupled / hairpin / lumped /
combline / interdigital (band-pass) + stepped-impedance (low-pass), all live. Clean mirror of
combline; no new physics (the engine + layout already shipped + are gated).

**Gate:** a non-vacuous host test (`design_demo_from(.., Interdigital)` yields a layout that
differs from edge-coupled AND hairpin AND combline for the same spec, shared coupling/verdict,
`combline_loading_cap_f` None, resonator length > 0) + `dx build` EXIT 0 + no regression.

**Not in scope:** the via / 3-D short render (the alternating grounding is in the comb geometry,
not extra UI); precise tap/Qe→feed (F1.2.1). EM-verify wall (ADR-0133/0147) untouched.

**Milestone:** with the gallery complete, the easy gallery wins are exhausted — the next
direction (high-pass response class, the maintainer-funded EM fork, or a new product
direction) is surfaced to the maintainer.

---

## References
- `crates/yee-studio-web/src/{engine.rs (derive_geometry Combline arm), stages.rs (Topology /
  technique_status / topology_response / gallery)}`; `crates/yee-filter/src/dimension.rs`
  (`dimension_interdigital` / `dimension_interdigital_layout`).
- `docs/superpowers/specs/2026-05-31-app-2-8-studio-interdigital-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-8-studio-interdigital.md`.
