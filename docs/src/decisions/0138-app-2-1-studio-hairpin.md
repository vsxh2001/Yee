# ADR-0138: App.2.1 — Light the Hairpin technique in the studio

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0136 (App.2.0 recommender — recommends Hairpin; this makes it routable),
F1.2.2 (the `dimension_hairpin` engine, shipped but un-surfaced), ADR-0130 (the studio),
the product vision (`docs/superpowers/specs/2026-05-31-ideal-filter-design-app-vision.md` §5),
[[lumped-lc-and-studio-redesign]]

---

## Context

The hairpin band-pass dimensioner (`dimension_hairpin` / `dimension_hairpin_layout`)
shipped in F1.2.2 but the studio gallery still greys **Hairpin** as "Soon" — the engine
is un-surfaced, and the App.2.0 recommender can only route a Hairpin recommendation to
the edge-coupled stand-in. The product vision (§5) calls for filling the gallery.

A hairpin filter is the **same coupled-resonator band-pass synthesis** as edge-coupled
(identical prototype / coupling matrix / swept response / mask verdict); only the
physical realization differs (U-folded λ/4 arms vs straight λ/2 lines). So surfacing it
is low-risk: everything in the studio's `Designed` except the geometry-derived fields is
already correct for hairpin, and the board renders generically from `yee_layout::Layout`.

## Decision

Light Hairpin as a live studio technique (`crates/yee-studio-web` only):

- Add `Topology::Hairpin` (rail = the distributed six). Thread the topology into the
  engine's geometry derivation: `design_demo_from(spec, topology)` →
  `derive_geometry(project, topology)` branches edge-coupled vs the existing
  `dimension_hairpin*`; the shared synthesis/response/verdict fields are computed once.
  The `designed` memo recomputes on (spec, topology).
- The Hairpin gallery card becomes selectable; the App.2.0 recommender's
  `technique_status` maps `Hairpin → Live(Topology::Hairpin)`. Layout/Export render the
  hairpin board (generic `Layout`) + Gerber/KiCad from the real hairpin layout.

## Consequences

**Ships:** the studio fills a gallery card — Hairpin is a real, routable band-pass
technique driven by the already-validated `dimension_hairpin` engine; the recommender's
Hairpin recommendation now routes to its real dimensioner. Visible end-to-end app
breadth, low risk (shared synthesis, generic board render, no engine change).

**Gate (studio):** `dx build --platform web --release` EXIT 0; a NEW non-vacuous host
test asserts `design_demo_from(demo_spec(), Topology::Hairpin)` yields a real layout that
**differs** from the edge-coupled layout for the same spec (the card routes to the real
dimensioner, not a stub); existing edge-coupled + lumped engine tests unregressed;
clippy/fmt/check clean.

**Not in scope:** Combline / Interdigital (no engine yet); the stepped-Z low-pass studio
path (needs a low-pass response path — separate increment); the hairpin `qe`→tap-offset
refinement (deferred in the engine since F1.2.1). The EM-verify wall (ADR-0133) is
untouched.

---

## References
- `crates/yee-studio-web/src/{engine.rs, stages.rs, main.rs}`;
  `crates/yee-filter/src/dimension.rs` (`dimension_hairpin*`, `HairpinDimensions`).
- `docs/superpowers/specs/2026-05-31-app-2-1-studio-hairpin-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-1-studio-hairpin.md`.
