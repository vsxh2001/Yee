# ADR-0105: Filter Phase F1.4.1b — `yee-export` KiCad `.kicad_pcb` S-expr export

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0100 (F1.4.0 Gerber walking skeleton), ADR-0103 (F1.4.1a Gerber
Edge.Cuts outline), ADR-0089 (desktop+web app — WASM-safe light flow),
`FILTER-DESIGN-ROADMAP.md`

---

## Context

F1.4.0 (ADR-0100) and F1.4.1a (ADR-0103) gave `yee-export` two RS-274X Gerber
emitters: copper (`layout_to_gerber`) and board outline
(`layout_to_gerber_outline`). Gerber is the fab hand-off format, but it is
awkward to *open and edit*. The product goal (re-confirmed by the user on
2026-05-30: "web ui for filter design, incl topologies, full filter simulation,
**kicad export**") names KiCad explicitly — and "KiCad export" most directly
means a `.kicad_pcb` board file the user can open in the KiCad PCB editor, not
just Gerbers.

## Decision

Add a third pure-text emitter to `yee-export`:
`layout_to_kicad_pcb(&Layout, &KicadPcbOptions) -> String`, producing a KiCad 7
S-expression board (`(kicad_pcb (version 20221018) … )`) with:

- a `(layers …)` table declaring `F.Cu` / `B.Cu` / `Edge.Cuts`;
- one filled `(gr_poly … (layer "F.Cu") (fill solid))` per `Layout` trace
  polygon;
- the board outline as a `(gr_poly … (layer "Edge.Cuts") (width 0.1) (fill
  none))` rectangle = `bbox ± outline_margin_mm` (the same geometry as the
  F1.4.1a Gerber outline).

Coordinates are KiCad mm **floats** (metres × 1e3, ~6 decimals) via a new private
`xy_mm` helper — deliberately distinct from the Gerber `mm_to_fixed46` integer
fixed-point, because the two formats use different coordinate encodings.

Walking-skeleton scope only: copper + outline. Footprints, pads, vias, drill,
zones, net classes, 3-D models, B.Cu routing, and KiCad-display Y-axis
reconciliation are explicitly deferred (F1.4.1c+). The emitter stays pure
`String` / WASM-safe so the studio can offer client-side KiCad export.

## Consequences

**Ships:** `yee-export::layout_to_kicad_pcb` — the pipeline can now emit a
`.kicad_pcb` the user opens directly in KiCad, alongside the Gerbers.

**Gate (structural, like gerber-001):** `cargo test -p yee-export` passes;
`kicad-001` asserts the S-expr starts with `(kicad_pcb`, is paren-balanced, names
`F.Cu` + `Edge.Cuts` in its layer table, and has one `F.Cu` `gr_poly` per trace
plus one `Edge.Cuts` `gr_poly`; `kicad-002` parses the `(xy …)` pairs back and
confirms they equal the trace vertices (mm) and the bbox±margin rectangle. A
structural gate (not a "KiCad opens it" gate) is the bar — the CI box has no
KiCad — consistent with the gerber-001 structural-gate precedent for export
writers.

**Not in scope:** CLI `--kicad-pcb` flag / studio export button (crosses into
`yee-cli` / `yee-studio` lanes — a follow-on); footprints/pads/vias/zones; Y-axis
convention reconciliation.

---

## References
- ADR-0100 / ADR-0103 (the Gerber emitters this sits beside).
- `docs/superpowers/specs/2026-05-30-filter-f1-4-1b-kicad-pcb-export-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-f1-4-1b-kicad-pcb-export.md`.
- KiCad 7 board file format (`.kicad_pcb`, S-expression, `version 20221018`).
