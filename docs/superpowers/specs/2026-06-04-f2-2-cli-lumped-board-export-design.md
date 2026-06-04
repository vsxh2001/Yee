# F2.2-cli — lumped-board manufacturing export in `yee filter synth` — design

**Date:** 2026-06-04
**ADR:** [ADR-0158](../../src/decisions/0158-f2-2-cli-lumped-board-export.md)

## Problem

The lumped track is library-complete (`synthesize_lumped` → `lumped_board` → `board.layout` →
`layout_to_gerber/kicad`) but `yee filter synth` only exports the planar edge-coupled layout. Lumped
manufacturing files aren't CLI-reachable — a planar↔lumped parity gap.

## Goal

A `--lumped` path in `yee filter synth` that routes the existing `--gerber`/`--kicad-pcb`/
`--layout-svg` writers to the lumped `board.layout`, gated by `cli-lumped-export`.

## Non-goals

Studio-UI lumped surfacing (framework-verdict-gated), drill/multi-layer/silkscreen, on-board value
annotation. CLI-only, single copper layer + outline (matching the planar export today).

## Architecture

- `crates/yee-cli/src/main.rs`: add `--lumped` (bool) + `--footprint <0402|0603|0805>` (default
  `0603`) flags to the `filter synth` subcommand; thread them into `run_synth`.
- `crates/yee-cli/src/filter.rs::run_synth`: add a `lumped: bool` + `footprint` param. In the export
  block (currently planar-only, ~line 189): if `lumped`, build `let board =
  lumped_board(&ladder, &substrate, footprint); let layout = board.layout;` (the `ladder` from
  `synthesize_lumped(&proj)`; the `substrate` from the existing `--eps-r`/`--h-mm`), and feed THAT
  `layout` to the same `layout_to_gerber`/`layout_to_kicad_pcb`/SVG writers. The single-`Layout`
  invariant holds (one `board.layout` for all three). Non-lumped path unchanged.
- The non-export console summary (synthesis printout) stays; optionally print the lumped BOM/
  placements if `--lumped` (nice-to-have, not required for the gate).

## Data flow (lumped path)

`FilterSpec` → `synthesize` (proto/coupling) → `synthesize_lumped` → `LumpedLadder` →
`lumped_board(ladder, substrate, footprint)` → `board.layout` → `layout_to_gerber`/`_kicad_pcb`/SVG.

## Gate `cli-lumped-export` (fast, non-ignored)

In the yee-cli tests (mirror the planar `cli_gerber`/`cli_kicad_pcb` gate): run the synth handler (or
the CLI binary) on a small bandpass spec with `--lumped --gerber <tmp>`; assert the Gerber is
non-empty, structurally valid (has the RS-274X `%FS`/`G36`/`G37`/`M02` markers the planar gate
checks), its copper-region (`G36`) count ≥ 2·N (N components → ≥2 pads each) + the signal line/rail,
AND it differs from the planar Gerber for the same spec (proves the `--lumped` branch is taken, not
a silent planar fall-through). Pure-compute; runs in the normal `cargo test` (no boxed/heavy run).

## Risks

- `lumped_board` needs a `Substrate` + `Footprint`; confirm the CLI already builds a `Substrate`
  from `--eps-r`/`--h-mm` (the planar dims path does) — reuse it. (low)
- `synthesize_lumped` may error for some specs (returns `LumpedError`); handle it (the gate uses a
  spec known to synthesize, like the `lumped_001` fixture). (low)
