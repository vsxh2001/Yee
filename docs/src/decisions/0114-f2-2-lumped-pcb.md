# ADR-0114: Filter Phase F2.2 — lumped-LC PCB board generator

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0111 (F2.0 LC ladder), ADR-0112 (F2.1 BOM), ADR-0097/0109
(distributed `*_layout` → `yee_layout::Layout` pattern), ADR-0100/0105 (Gerber/
KiCad export), the lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

The lumped-LC goal names "to the pcb level." F2.0 gives the LC ladder and F2.1
the BOM, but nothing places the components on a board. The distributed track
already maps a design to a `yee_layout::Layout` (`dimension_edge_coupled_layout`)
that `layout_to_gerber`/`layout_to_kicad_pcb` render — the lumped track needs the
analogous board generator.

## Decision

Add `crates/yee-filter/src/board.rs`: `lumped_board(&LumpedLadder, &Substrate,
Footprint) -> LumpedBoard` placing each resonator's L/C as **SMD footprints (2
copper pads each)** along a 50 Ω signal microstrip with a ground rail — series
branches in-line, shunt branches on a stub to ground — and returning a
`yee_layout::Layout` (so existing Gerber/KiCad export works) plus a `Placement`
list (ref-des → footprint → position) for BOM cross-ref. `Footprint` = SMD
0402/0603/0805 (IPC land patterns). Pure-geometry, WASM-safe; yee-filter already
deps yee-layout (no cycle).

Gate `lumped_pcb_001` (cheb N=5): `2·N` placements with unique ref-des, the right
pad count, **no pad overlap** (rect-disjoint), a finite positive-area bbox
containing all pads, and correct series/shunt y-placement.

## Consequences

**Ships:** a manufacturable lumped-LC board (`Layout` + placements) — the goal's
"to the pcb level." Feeds Gerber/KiCad export (already shipped) and F2.3 (the
voxelized board for FDTD).

**Gate:** `cargo test -p yee-filter` green incl. `lumped_pcb_001`. Pure-geometry.

**Not in scope:** KiCad-native `(footprint)` objects + courtyards/3D (F2.2b);
soldermask/silkscreen; auto-routing/meander matching; FDTD (F2.3); UI.

---

## References
- `docs/superpowers/specs/2026-05-30-f2-2-lumped-pcb-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-2-lumped-pcb.md`. IPC-7351 land patterns.
