# ADR-0168: `yee filter synth --jlcpcb` auto-routes the orderable topology (T4)

**Status:** Accepted (T4 — the ADR-0167 verdict'd follow-on)
**Date:** 2026-06-06
**Related:** ADR-0167 (T3 — `synthesize_orderable` selector), ADR-0166 (T2 — `top_c_board`), ADR-0164 (J4 —
`yee filter synth --jlcpcb` writes Gerber+BOM+CPL for the fixed alternating ladder),
[[project-filter-design-final-goal]] (spec → JLCPCB-orderable board + BOM, one command).

---

## Context

ADR-0167 (T3) shipped `synthesize_orderable(&FilterProject, Footprint) -> OrderableBoard` — the routing brain
that returns whichever lumped topology (alternating ladder / top-C) yields a fully-orderable board, or an
honest `fully_orderable=false`. But the user-facing CLI `yee filter synth --jlcpcb <dir>` (J4, ADR-0164) still
hardcodes the **alternating ladder** (`write_jlcpcb_set`: `synthesize_lumped` → `lumped_board` →
`join_placed_parts`). So a narrow-band spec the ladder can't make orderable emits a half-blank BOM even though
top-C would resolve it — the selector's value is not reaching the user.

T4 wires the selector into the CLI: `--jlcpcb` auto-routes the topology and reports which one was chosen +
whether the board is fully orderable. **Wrinkle:** `synthesize_orderable` places on an internal FR-4 reference
substrate (for its pure-compute gate), but the CLI honors the user's `--eps-r`/`--h-mm` (the board geometry's
connecting-trace width is `microstrip_width(z0, eps_r, h)`). T4 must NOT regress that.

## Decision

1. **`yee-filter`: add `synthesize_orderable_on(project, substrate: &Substrate, footprint) -> Result<OrderableBoard, _>`** —
   identical policy to `synthesize_orderable` but places both candidate boards on the GIVEN substrate.
   `synthesize_orderable(project, footprint)` becomes a thin delegate (`reference_substrate()`), so the T3 API
   and the `topology-select-001` gate are unchanged. (The topology DECISION + BOM orderability are
   substrate-independent — component values come from the LC synthesis — so routing is identical; only the
   board geometry/CPL coords differ.)
2. **`yee-cli`: `write_jlcpcb_set` calls `synthesize_orderable_on(proj, substrate, footprint)`** instead of the
   fixed-ladder path. Emit the copper + outline Gerber from `ob.board.layout`, the CPL from
   `ob.board.placements`, the BOM from `ob.parts` (reusing `jlcpcb_bom_csv`/`jlcpcb_cpl_csv` unchanged).
   Report the chosen topology + orderability:
   - `println!("  topology: {alternating ladder | top-C-coupled} (auto-selected)")`.
   - If `ob.fully_orderable` → `"  all N parts matched a JLCPCB Basic part"` (as today).
   - Else → the honest `"  NOTE: M of N parts have no JLCPCB Basic match …"` PLUS, when BOTH topologies
     blanked, a `"  (neither lumped topology is fully orderable for this spec; consider the distributed/planar
     track)"` pointer. Never fabricate orderability.

**Gate `crates/yee-cli/tests/`** (`cli-jlcpcb-autoroute`), invoking `run_synth` with `--jlcpcb <tmpdir>`:
- **0.5 GHz/20 %/0402** (the T3 discriminating spec) → the run reports **top-C-coupled** + the written
  `bom.csv` has **zero blank LCSC #s** (every row carries a `C\d+`); CPL designators == BOM designators.
- **1 GHz/70 %/0402** (wideband) → reports **alternating ladder** + zero-blank BOM.
- **2 GHz/5 %/0402** (GHz-narrow) → reports `fully_orderable=false`, the note fires, and `bom.csv` carries
  the real blank rows (honest, not dropped).
Non-circular: the gate reads the actually-written CSV files + the run's reported topology, not the in-process
selector return.

## Consequences

- **The deliverable's headline path works end-to-end:** `yee filter synth <spec> --jlcpcb <dir>` returns an
  orderable upload set across the broadest spec range either lumped topology covers, naming the topology — the
  user runs one command, no topology knowledge required.
- **T5 follow-on (noted, not in T4):** the studio Export stage's JLCPCB download buttons call
  `synthesize_orderable`/`_on` + surface the chosen topology + orderability badge (the WASM/UI lane).
- Scope T4: `crates/yee-filter/src/{topology.rs, lib.rs}` (the `_on` variant + re-export) + `crates/yee-cli/src/filter.rs`
  + `crates/yee-cli/tests/`. The `_on` addition keeps WASM-safety (pure data). This ADR is the design record.
- **Not in scope:** the studio wiring (T5); a `--topology` override flag (the auto-route is the default; a
  manual override is a trivial later add); J5 Gerber-completeness (ADR-0164).

## References
- Selector: `yee_filter::{synthesize_orderable, OrderableBoard, BoardTopology}` (ADR-0167).
- CLI: `yee_cli::filter::{run_synth, write_jlcpcb_set, JLCPCB_SERIES}` (ADR-0164 J4).
- Boards/CSVs: `yee_filter::{lumped_board, top_c_board, jlcpcb_bom_csv, jlcpcb_cpl_csv}`,
  `yee_export::{layout_to_gerber, layout_to_gerber_outline}`.
