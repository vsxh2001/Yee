# Filter Phase F1.4.1b — `yee-export` KiCad `.kicad_pcb` export — Plan

**Spec:** `2026-05-30-filter-f1-4-1b-kicad-pcb-export-design.md` · **ADR:** ADR-0105

## Lane
`crates/yee-export/**` ONLY (`src/lib.rs`, `tests/`). Do NOT edit `yee-layout`/
other crates — consume `yee_layout::{Layout, BBox, Point2, Polygon}`. Out of lane
→ finding. Keep `yee-export` WASM-safe (pure `String`; no fs / native dep).

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-f1-4-1b-kicad-pcb`.

## Pattern files
- `crates/yee-export/src/lib.rs` — the F1.4.0 `layout_to_gerber` + the F1.4.1a
  `layout_to_gerber_outline` + `OutlineOptions` + the private `coord_word_xy`
  helper + module coordinate-model docs to MIRROR (house style: options struct
  with `Default`, a `pub fn layout_to_X(&Layout, &XOptions) -> String`, a private
  coordinate helper, every public item documented).
- `crates/yee-export/tests/gerber_003_outline_structure.rs` +
  `gerber_004_outline_geometry.rs` — the test idiom (build a `Layout`, assert
  structure; parse coords back and compare).
- `crates/yee-layout/src/lib.rs` — `BBox { min, max }`, `Point2 { x, y }`,
  `Polygon { verts }`, `Layout { substrate, bbox, traces, ports }`, `Substrate
  { height_m, .. }`, `Polygon::rect`, `BBox::from_polygons`.

## Steps
1. `src/lib.rs`: add `KicadPcbOptions` (+ `Default`) and `layout_to_kicad_pcb`
   per the spec (kicad_pcb header + layers + one `gr_poly` per trace on the
   copper layer + an Edge.Cuts outline `gr_poly` = bbox±`outline_margin_mm`).
   Add a private `xy_mm(x_m, y_m)` (metres → mm float, ~6 decimals). Reuse the
   bbox±margin corner math from the F1.4.1a outline emitter (lower-left CCW).
   Document every public item; explain why `xy_mm` is a float helper distinct
   from the Gerber `mm_to_fixed46`.
2. `tests/kicad_001_structure.rs` + `tests/kicad_002_geometry.rs` per DoD 4–5.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-export --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-export --jobs 2
```
Pure text — sub-second. Do NOT run `cargo test --workspace`, FDTD, mom-001. This
box is MEMORY-CONSTRAINED: keep to `-p yee-export` with `--jobs 2`; do NOT build
the whole workspace (it can OOM-kill).

## Escape hatch
Blocked > 15 min — uncertainty on the KiCad 7 S-expr header tokens (version /
layers table), or whether `gr_poly` needs an explicit closing point (it does NOT
— KiCad closes `gr_poly` implicitly), or a `bbox` accessor is missing → STOP and
surface. Do NOT add footprints/pads/vias/zones; do NOT add a file writer; do NOT
edit yee-layout. A structurally-valid skeleton (opens conceptually; paren-balanced;
correct layers + gr_poly per trace) is the bar — do NOT chase a
pixel-perfect-in-KiCad result you cannot verify without KiCad installed.

## Done when
DoD 1–5 pass; `git diff --stat <base>..HEAD` = only `crates/yee-export/**` + the
3 committed docs; `yee-export` still pure/WASM-safe (no fs/native dep).
