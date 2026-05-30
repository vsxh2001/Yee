# Filter Phase F1.4.1b — `yee-export` KiCad `.kicad_pcb` S-expr export — Design Spec

**ADR:** ADR-0105 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Emit a KiCad-native board file (`.kicad_pcb`, KiCad 7 S-expression format) for a
`yee_layout::Layout`, so the filter-design pipeline produces a file the user can
open directly in KiCad's PCB editor — not just the intermediate Gerber. This is
the literal "KiCad export" endpoint of the product goal. Walking-skeleton scope:
top-copper trace polygons on `F.Cu` + the board outline on `Edge.Cuts`. Pure
text, WASM-safe (same constraint as the Gerber emitters).

## Why KiCad-native in addition to Gerber
Gerber (F1.4.0/F1.4.1a) is the fab hand-off format but is awkward to *edit*. A
`.kicad_pcb` opens in KiCad for inspection, tweaking, and re-routing before the
user re-plots Gerbers. Both are pure-text writers over the same `Layout`.

## Changes (`crates/yee-export/**` ONLY)
- `src/lib.rs`: add
  - `KicadPcbOptions { copper_layer: String (default "F.Cu"), outline_margin_mm: f64 (default 1.0), generator: String (default "yee-export") }` + `Default`.
  - `pub fn layout_to_kicad_pcb(layout: &Layout, opts: &KicadPcbOptions) -> String`.
  - a private `xy_mm(x_m, y_m) -> String` (metres → mm **float**, 6 significant
    decimals, e.g. `3.059`) — KiCad S-expr coordinates are mm floats, NOT the
    Gerber 4.6 fixed-point integers, so this is a separate helper (do NOT reuse
    `mm_to_fixed46`). Document why.
- Document every public item (`#![warn(missing_docs)]`).

### Emitted structure (KiCad 7, `version 20221018`)
```
(kicad_pcb
  (version 20221018)
  (generator "<opts.generator>")
  (general (thickness <substrate.height_m*1e3>))
  (paper "A4")
  (layers
    (0 "F.Cu" signal)
    (31 "B.Cu" signal)
    (44 "Edge.Cuts" user)
  )
  (setup)
  ; one filled gr_poly per trace polygon, on the copper layer
  (gr_poly (pts (xy X Y) (xy X Y) ...) (layer "<opts.copper_layer>") (width 0) (fill solid))
  ...
  ; the board outline as a gr_poly on Edge.Cuts (bbox ± outline_margin_mm), width 0.1, fill none
  (gr_poly (pts (xy X Y) (xy X Y) (xy X Y) (xy X Y)) (layer "Edge.Cuts") (width 0.1) (fill none))
)
```
The board-outline rectangle reuses the SAME bbox±margin geometry as
`layout_to_gerber_outline` (lower-left CCW); the Edge.Cuts `gr_poly` is closed by
KiCad implicitly (no explicit repeat of the first point needed for `gr_poly`).

Coordinates: KiCad mm, no Y-flip in the skeleton (the studio/Gerber and KiCad
share the layout's metre frame scaled to mm; a Y-axis convention pass is a later
brick if a round-trip mismatch ever shows).

## DoD (machine-checkable; pure text, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-export --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-export` exit 0 (sub-second; pure text).
4. New gate `tests/kicad_001_structure.rs`: build a `Layout` with ≥2 rect
   traces; assert the output (a) starts with `(kicad_pcb`, (b) has balanced
   parentheses, (c) contains a `(layers` block naming `F.Cu` and `Edge.Cuts`,
   (d) has exactly one `(gr_poly` with `(layer "F.Cu")` per trace polygon, (e)
   has exactly one `(gr_poly` with `(layer "Edge.Cuts")`.
5. New gate `tests/kicad_002_geometry.rs`: parse the `(xy ...)` pairs out of the
   `F.Cu` polygons and confirm they equal the trace vertices in mm (metres×1e3)
   within 1e-6, and that the Edge.Cuts rectangle equals `bbox ± outline_margin_mm`.

## Out of scope
Footprints, pads, vias, drill, net classes, zones, 3-D models, B.Cu routing,
Y-axis-convention reconciliation with KiCad's display frame. Those are F1.4.1c+.
CLI/studio wiring (an `--kicad-pcb` flag / export button) is a separate follow-on
(it crosses into `yee-cli`/`yee-studio` lanes) — define the API here only.
