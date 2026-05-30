# Filter Phase F2.2 — lumped-LC PCB board generator — Design Spec

**ADR:** ADR-0114 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

The lumped-LC goal names **"to the pcb level."** F2.2 = place the LC ladder's
components as **SMD footprints + pads + connecting traces on a board**, producing
a `yee_layout::Layout` (so the existing Gerber / KiCad export works) plus a
placement list (ref-des → footprint → position) for BOM cross-ref. Pure-geometry,
WASM-safe, no FDTD.

## Placement (walking skeleton)

A horizontal **signal microstrip** runs left→right at 50 Ω width; a **ground rail**
(copper) runs along the bottom. For each `LcResonator` (left to right):

- **Series branch** (series L–C): break the signal line with a gap and place the
  component footprints **in-line** bridging the gap (L then C pads in series along
  the signal path).
- **Shunt branch** (parallel L–C): place the component footprints on a short
  **stub** dropping from the signal line to the ground rail (L and C side by side,
  each from line to ground).

Each component is one SMD footprint = **two rectangular copper pads** at the
footprint's pad pitch. `Footprint` enum: `Smd0402 | Smd0603 | Smd0805` with
pad length/width/pitch from standard IPC land patterns (0603 default). Pads are
added to `Layout.traces` as copper `Polygon` rects; the signal line + ground rail
are copper rects too.

## Changes (`crates/yee-filter/**` ONLY — yee-filter already deps yee-layout, no cycle)

- New `crates/yee-filter/src/board.rs`:
  - `pub enum Footprint { Smd0402, Smd0603, Smd0805 }` with `fn pad() -> PadSpec`
    (pad_len_m, pad_w_m, pitch_m) from IPC land patterns.
  - `pub enum BranchKind { Series, Shunt }`
  - `pub struct Placement { ref_des: String, footprint: Footprint, kind: BranchKind,
    center_m: (f64, f64) }`
  - `pub struct LumpedBoard { layout: yee_layout::Layout, placements: Vec<Placement> }`
    (+ serde where types allow).
  - `pub fn lumped_board(ladder: &LumpedLadder, substrate: &yee_layout::Substrate,
    footprint: Footprint) -> LumpedBoard`: lay out the signal line + ground rail +
    per-resonator L/C footprints per the rules above; ref-des `L1,C1,L2,C2,…`;
    record placements; assemble the `Layout` (substrate, all copper rects as
    `traces`, `bbox` from the polygons, ports at the two line ends).
- Re-export from `lib.rs`.

## DoD (machine-checkable; pure-geometry, NO FDTD)

1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-filter` exit 0 — incl. gate `lumped_pcb_001` (cheb N=5):
   - `placements.len() == 2·N` (an L + a C footprint per resonator); ref-des unique.
   - Pad count in `layout.traces` ≥ `2·placements.len()` (2 pads/footprint) plus
     the signal line + ground rail rects.
   - **No pad overlap:** every pair of pad rectangles is disjoint (axis-aligned
     rect intersection test) — placement spacing is valid.
   - `layout.bbox` is finite, positive-area, and contains all pads.
   - Each shunt footprint reaches the ground rail; each series footprint sits on
     the signal line (assert y-coordinates by branch).

## Out of scope

KiCad-native `(footprint)`/pad objects (this skeleton emits pads as copper
polygons that the existing `layout_to_gerber`/`layout_to_kicad_pcb` already
render; proper footprint objects + courtyards + 3D = F2.2b); solder mask /
silkscreen designators; auto-routing / impedance-matched meander; the FDTD sim
(F2.3); UI. Component VALUES come from F2.1's BOM (cross-referenced by ref-des).

## Why now

Unblocked, dispatchable, mostly pure-geometry, and directly delivers the goal's
"to the pcb level." Depends only on shipped F2.0 (`LumpedLadder`) +
`yee_layout` (already a yee-filter dep). The existing export turns the result
into Gerber/KiCad immediately.
