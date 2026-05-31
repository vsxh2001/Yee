# ADR-0146: App.2.7 — Light the Combline technique in the studio

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0144 (F1.2.5 combline engine), ADR-0145 (F1.2.6 combline layout — the
two prerequisites), ADR-0138 (App.2.1 hairpin lighting — the pattern this mirrors),
ADR-0136 (the recommender recommends combline), ADR-0142/0143 (compare/overlay — combline
joins them automatically), [[lumped-lc-and-studio-redesign]]

---

## Context

Combline's synthesis (F1.2.5) and board layout (F1.2.6) shipped, but the studio gallery
still greys it "Soon" and the recommender can only route the edge-coupled stand-in. With
both the engine and the layout in place, lighting combline is a clean mirror of the
App.2.1 hairpin lighting (which itself became trivial once `dimension_hairpin_layout`
existed). Combline is a coupled-resonator **band-pass** technique — same synthesis as
edge-coupled/hairpin; only the realization (short-circuited θ0 resonators + loading caps)
and the board differ.

## Decision

Light combline as a live studio technique (`crates/yee-studio-web`):

- `Topology::Combline` (rail = distributed); a `derive_geometry` arm →
  `dimension_combline` + `dimension_combline_layout` (θ0 = π/4 = λg/8, the compact
  default); the shared synthesis/response/verdict computed once.
- The Combline gallery card becomes selectable; `technique_status(Combline)` →
  `Live`; `topology_response(Combline)` = `Bandpass` (so the recommender routes it).
- The layout stage renders the combline board (generic `board_svg`) + the resonator
  table + the combline-distinct **loading cap C_L** (a single value — uniform θ0/Z0 →
  same `C_L = cot(θ0)/(2π·f0·Z0)` per resonator). Compare + overlay (App.2.5/2.6) iterate
  live techniques, so combline joins both automatically; export emits Gerber/KiCad.

## Consequences

**Ships:** combline is a live, routable studio technique end-to-end (synthesis → board →
Gerber/KiCad, in compare + overlay), honestly surfacing the loading cap. Completes the
maintainer's combline pick. The gallery now has four live band-pass techniques
(edge-coupled / hairpin / lumped / combline) + low-pass (stepped-Z). Clean mirror of
hairpin; no new physics (the engine + layout already shipped + are gated).

**Gate:** a non-vacuous host test (`design_demo_from(.., Combline)` yields a layout that
differs from edge-coupled AND hairpin for the same spec, shared coupling/verdict, and a
real `C_L > 0`) + `dx build` EXIT 0 + no regression.

**Not in scope:** the SMD-cap hybrid board render (caps surfaced as a value/table line,
not drawn as footprints — a polish follow-on); discrete E-series `C_L` selection;
interdigital (a separate technique, now gateable via the same H&L Qe/M approach).

---

## References
- `crates/yee-studio-web/src/{engine.rs (derive_geometry Hairpin arm), stages.rs
  (Topology / technique_status / topology_response / gallery)}`;
  `crates/yee-filter/src/dimension.rs` (`dimension_combline` / `dimension_combline_layout`).
- `docs/superpowers/specs/2026-05-31-app-2-7-studio-combline-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-7-studio-combline.md`.
