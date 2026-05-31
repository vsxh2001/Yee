# ADR-0143: App.2.6 — Multi-technique response overlay

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0142 (App.2.5 compare table — this is its visual companion), ADR-0136
(recommender), ADR-0138/0139 (the hairpin + low-pass flows), ADR-0133 (the EM wall —
why only *ideal/realized-circuit* responses are overlaid, not a realized distributed EM
response), [[lumped-lc-and-studio-redesign]]

---

## Context

The App.2.5 Compare panel tabulates per-technique board size + verdict + metrics but
shows no response shape. The maintainer's "deepen the flows — optimize/compare"
direction calls for the user to *see* the responses, not just the numbers. A single
chart overlaying the techniques' swept `|S21|` for the current spec is the visual
companion to the compare table.

**Honesty constraint:** edge-coupled and hairpin share the *same* coupled-resonator
synthesis (identical coupling matrix → identical ideal `|S21|`); they differ only
physically (board layout/size — shown in the compare table), not in the ideal response.
A realized *distributed* EM response would differ, but that needs full-wave EM (the
deferred ADR-0133 wall). So the genuinely distinct curves are the coupled-resonator
ideal (shared by edge-coupled + hairpin) and the lumped realized-ladder; for low-pass,
the stepped ideal.

## Decision

Add a response overlay to the Compare panel via a pure helper:

- `overlay_curves(spec) -> Vec<OverlayCurve>` returns the **distinct** swept responses
  for the spec's response class — band-pass: the coupled-resonator ideal (labelled as
  edge-coupled / hairpin) + the lumped realized-ladder; low-pass: the stepped ideal;
  high-pass: none. Real engine sweeps on the shared frequency grid, **labelled
  truthfully** (the distributed techniques are one shared ideal curve, not two).
- `svg::response_overlay(curves, bands)` renders one `|S21|` polyline per curve in a
  distinct colour + a legend + the shaded mask bands, mirroring `response_plot`.
- The Compare panel renders it below the table.

## Consequences

**Ships:** the Compare view is complete — table (board + metrics) **and** chart
(response shape vs mask) — honestly showing the distributed techniques' shared ideal
response and the lumped realized-ladder as distinct curves. Built on the existing
per-flow sweeps + the `response_plot` renderer; no new physics.

**Gate (non-vacuous):** `overlay_curves` on a band-pass spec returns 2 curves on the
same grid whose sweeps are **not identical** (lumped realized ≠ coupled-resonator
ideal — the real difference the overlay exists to show), each equal to the corresponding
design's `.sweep`; low-pass → 1; high-pass → `[]`; `response_overlay` renders one
polyline per curve + a legend; `dx build` EXIT 0; no regression.

**Not in scope:** a realized *distributed* response (the EM wall, ADR-0133); interactive
zoom; `|S11|` overlay. The labels stay honest — no faked per-technique distributed
curves.

---

## References
- `crates/yee-studio-web/src/{engine.rs (sweeps, design_*_from), svg.rs (response_plot),
  stages.rs (compare_panel)}`.
- `docs/superpowers/specs/2026-05-31-app-2-6-response-overlay-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-6-response-overlay.md`.
