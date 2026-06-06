# ADR-0164: JLCPCB production-ready output track — LCSC part autopick + BOM + CPL

**Status:** Accepted (track kickoff; maintainer goal-refinement 2026-06-06)
**Date:** 2026-06-06
**Related:** ADR-0100 (`layout_to_gerber`), ADR-0103 (board-outline Gerber), ADR-0105/0106 (KiCad export),
F2.1 `parts.rs` (E-series BOM `select_components`), F2.2 `board.rs` (`lumped_board` → `Layout` +
`Placement`), [[project-filter-design-final-goal]].

---

## Context

The maintainer refined the final goal (2026-06-06): **a full pipeline from spec & design to a JLCPCB
production-ready board + BOM, including auto-picking parts.** The lumped track already produces, from a
`FilterSpec`: synthesis → dimensioning → `lumped_board` (`Layout` of copper pads + a `Placement` list,
ref-des → footprint → position) → `layout_to_gerber` (single copper layer) + board-outline Gerber + KiCad
`.kicad_pcb`, and a value-selecting BOM (`select_components` → `Bom` of E-series-nearest L/C values). What
is MISSING for a JLCPCB-orderable result:

1. **Real part picking.** The BOM carries *values* (E-series), not orderable parts. JLCPCB assembly needs
   each line mapped to a real **LCSC part number** (ideally a JLCPCB "Basic" part — free assembly —
   matching value + footprint + a sane voltage/tolerance). This is the maintainer's "auto-picking parts."
2. **JLCPCB BOM CSV** — JLCPCB's assembly BOM schema: `Comment, Designator, Footprint, LCSC Part #`.
3. **JLCPCB CPL / centroid CSV** — `Designator, Mid X, Mid Y, Layer, Rotation` (from the `Placement` list).
4. **Gerber completeness** (polish) — JLCPCB fab wants copper + outline (have) and ideally mask/silk/drill.

## Decision

Build a **JLCPCB production-ready output** path in ordered bricks. Constraint: **WASM-safe, offline,
deterministic** (the studio runs client-side; no network/API key) — so the LCSC catalog is a **bundled,
curated, real-part table** (seeded from JLCPCB's published Basic Parts list; extensible), NOT a live query.
Research-first: source real LCSC C-numbers, don't invent them.

| # | Brick | Gate |
|---|-------|------|
| **J1** | **LCSC part autopick + bundled parts table.** A new `yee-filter` module (`jlcpcb.rs` or `parts_lcsc.rs`): a curated static table of real LCSC parts (common 0402/0603/0805 R/L/C values → LCSC C-number, value, footprint, basic/extended) seeded from JLCPCB's published Basic Parts list; `autopick(BomLine) -> Option<LcscPart>` picks the nearest in-table part of matching kind+footprint within the E-series tolerance, preferring Basic. | NON-circular: a synthesized 3-pole filter's every BOM line resolves to a real LCSC part whose value is within the E-series tolerance of the chosen value + footprint matches; the picked C-numbers are well-formed (`C\d+`) + present in the bundled table; basic-preferred. Document the table as a curated seed (coverage + how to extend). |
| J2 | **JLCPCB BOM CSV** (`Comment, Designator, Footprint, LCSC Part #`) from the autopicked `Bom` + the `Placement` ref-des. | Valid CSV, one row per distinct part (grouped designators), the JLCPCB column header exact; every row has an autopicked LCSC #. |
| J3 | **JLCPCB CPL/centroid CSV** (`Designator, Mid X, Mid Y, Layer, Rotation`) from the `Placement` list (mm, top layer). | Valid CSV, one row per component, designators match the BOM, coordinates within the board outline. |
| J4 | **CLI + studio wiring.** `yee filter synth <spec> --jlcpcb <dir>` writes Gerbers + `bom.csv` + `cpl.csv` (the JLCPCB upload set); the studio Export stage adds JLCPCB BOM/CPL downloads. | A CLI gate: the `--jlcpcb` dir contains a valid Gerber + bom.csv + cpl.csv; the BOM LCSC #s are non-empty; designators consistent across BOM/CPL. |
| J5 | **Gerber completeness** (if needed for fab): mask/silk/drill layers. | JLCPCB-loadable Gerber set (deferred / as-needed). |

**Start J1** — the autopick + parts table is the heart of the new ask and the prerequisite for J2/J4. The
LCSC table is the data risk; research JLCPCB's published Basic Parts list for real C-numbers (a curated
seed covering the common RF L/C decades), and gate that picks are real + value/footprint-correct.

## Consequences

- Closes spec→JLCPCB-orderable for the **lumped** track (the manufacturable path; the planar full-wave is
  a separate fidelity track per ADR-0162). The studio (the deliverable) can emit a JLCPCB upload set.
- The bundled LCSC table is a curated seed — honest about coverage (it will not cover every value); the
  picker returns `None` (surfaced, not faked) when no in-table part matches, and the table is documented
  as extensible from JLCPCB's published list. NO invented part numbers.
- Scope J1: `crates/yee-filter/src/` (new module + lib re-export) + `crates/yee-filter/tests/`. Pure
  data/`f64`, WASM-safe. This ADR is the design record; J1 may get a short spec if the table sourcing is
  non-trivial.
- **Not in scope:** a live LCSC/JLCPCB API (offline constraint), exotic footprints, the planar-filter
  JLCPCB path (lumped first), J5 Gerber layers until fab needs them.

## References
- BOM: `yee_filter::parts` (`select_components`, `Bom`, `BomLine`, `ESeries`).
- Board/placements: `yee_filter::board` (`lumped_board`, `Layout`, `Placement`, `Footprint`).
- Export: `yee_layout::layout_to_gerber`, KiCad (ADR-0105); CLI `yee filter synth` (`yee_cli::filter`).
- JLCPCB assembly BOM/CPL schema + Basic Parts list (researched in J1; cite the source in the table doc).
