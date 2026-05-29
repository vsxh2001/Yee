# ADR-0102: Filter — `yee filter synth --gerber` export wiring

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0098 (`yee filter synth` dims + `--layout-svg`), ADR-0100
(`yee-export::layout_to_gerber`), `FILTER-DESIGN-ROADMAP.md` (F1.4)

---

## Context

F1.4.0 (ADR-0100) shipped `yee_export::layout_to_gerber`, but no user command
calls it. `yee filter synth` already builds the substrate + the `Layout` for its
`--layout-svg` output (ADR-0098); writing the Gerber is the same plumbing with a
different emitter. Wiring it completes the "synthesize → fab file" path on the
CLI.

## Decision

Add a `--gerber <PATH>` flag to `yee filter synth`. When given, write
`yee_export::layout_to_gerber(&dimension_edge_coupled_layout(&project,
&substrate)?, &GerberOptions::default())` to the path (reusing the substrate +
layout already built for `--layout-svg`). `yee-cli` gains a `yee-export`
dependency.

## Consequences

**Ships:** `yee filter synth --gerber out.gbr`. Single copper layer (F1.4.0
skeleton). On a dimensioning `Err` (unrealizable coupling) the existing non-zero
exit + diagnostic already cover it (the Gerber write sits after the dims success
path).

**Gate (`yee-cli` test):** `cli_gerber` — for the committed Chebyshev N=5 fixture
(FR-4 default), `--gerber <tmp>` writes a file whose contents contain
`%FSLAX46Y46*%` and `M02*` (a structurally-valid single-layer Gerber).

**Not in scope:** board outline / drill / multi-layer / KiCad / STEP (F1.4.1+); a
studio export button (a follow-on). Pure CLI wiring over the shipped emitter.

---

## References
- ADR-0100 (`layout_to_gerber`); ADR-0098 (the `--layout-svg` plumbing it mirrors).
- `docs/superpowers/specs/2026-05-30-filter-cli-gerber-export-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-cli-gerber-export.md`.
