# ADR-0158: F2.2-cli — wire lumped-board manufacturing export into `yee filter synth`

**Status:** Accepted
**Date:** 2026-06-04
**Related:** ADR-0102/0106 (the planar `--gerber`/`--kicad-pcb` CLI wiring this mirrors for lumped),
the F2.0/F2.1/F2.2 lumped track (`synthesize_lumped`, `lumped_board`), ADR-0134 (lumped-LC goal
6/6), `FILTER-DESIGN-ROADMAP.md` §8 manufacturing export, [[project-filter-design-final-goal]].

---

## Context

The lumped-LC track is library-complete: `yee_filter::synthesize_lumped` → `LumpedLadder`, and
`yee_filter::lumped_board(ladder, substrate, footprint)` → a `LumpedBoard` whose `.layout` already
carries the signal line, ground rail, **and every component pad as copper polygons**, feedable
straight into `yee_export::{layout_to_gerber, layout_to_kicad_pcb}` (board.rs docs). But the CLI
`yee filter synth` only ever builds and exports the **planar edge-coupled** layout
(`dimension_edge_coupled`); its `--gerber`/`--kicad-pcb`/`--layout-svg` flags have no lumped path.
So a user can synthesize a lumped filter but cannot emit its board's manufacturing files from the
CLI — a **planar↔lumped parity gap** in the otherwise-complete lumped track. (Planar gets CLI export
via ADR-0102/0106; lumped does not.)

## Decision

Add a **`--lumped`** path to `yee filter synth` (with `--footprint <0402|0603|0805>`, default 0603):
when `--lumped` is set, the export flags build the geometry from `synthesize_lumped` →
`lumped_board(...)` → `board.layout` instead of the planar edge-coupled layout, then reuse the
**existing** `--gerber`/`--kicad-pcb`/`--layout-svg` writers verbatim (they take `&Layout`). The
SVG/Gerber/KiCad all derive from the **one** `board.layout`, so — exactly as the planar path's
single-`Layout` invariant — they can never diverge.

**Gate `cli-lumped-export`** (fast, pure-compute — NO FEM): `yee filter synth <bandpass spec>
--lumped --gerber <out>` produces a non-empty, structurally-valid RS-274X Gerber whose copper-region
count covers the board's components (≥ 2·N pad regions for an N-component board) + the signal
line/rail; assert it differs from the planar Gerber for the same spec (proves the lumped path is
actually taken, not silently falling through to planar). Mirror the shipped planar `cli_gerber` /
`cli_kicad_pcb` gate style.

## Consequences

**Ships lumped manufacturing-file CLI reachability** — planar↔lumped CLI parity; the lumped track's
"manufacturing export" stage is now user-reachable end-to-end (spec → `synthesize_lumped` →
`lumped_board` → Gerber/KiCad from one command). Pure-compute + deterministic → a fast non-`#[ignore]`
gate (no boxed heavy run). Small, self-contained CLI increment.

**Not in scope:** drill/multi-layer/silkscreen refinements (the planar export is also single
copper-layer + outline today), the studio-UI lumped surfacing (gated on the eframe-vs-Dioxus
framework verdict — out of scope here; this is CLI-only), per-component value annotation on the
board (values come from the F2.1 BOM keyed by ref-des, not the geometry).

---

## References
- Library path (complete): `yee_filter::lumped_board` (`crates/yee-filter/src/board.rs:181`) →
  `LumpedBoard.layout`; `yee_export::{layout_to_gerber, layout_to_kicad_pcb}`
  (`crates/yee-export/src/lib.rs:142/357`).
- The planar CLI export this mirrors: `yee-cli::filter::run_synth`
  (`crates/yee-cli/src/filter.rs:47`, the `--gerber`/`--kicad-pcb` branch ~line 189).
- Spec: `docs/superpowers/specs/2026-06-04-f2-2-cli-lumped-board-export-design.md`;
  plan: `docs/superpowers/plans/2026-06-04-f2-2-cli-lumped-board-export.md`.
