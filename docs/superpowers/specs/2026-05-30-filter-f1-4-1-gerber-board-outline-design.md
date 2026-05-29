# Filter Phase F1.4.1a — `yee-export` Gerber board outline — Design Spec

**ADR:** ADR-0103 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Add a board-outline (Edge.Cuts) RS-274X emitter to `yee-export`: a single closed
rectangular contour around the layout bbox + margin, stroked. Pure text,
WASM-safe, no FDTD.

## API (`crates/yee-export/src/lib.rs`)
```rust
pub struct OutlineOptions { pub layer_name: String, pub margin_mm: f64 }
impl Default for OutlineOptions { fn default() -> Self { Self { layer_name: "Edge.Cuts".into(), margin_mm: 1.0 } } }
pub fn layout_to_gerber_outline(layout: &yee_layout::Layout, opts: &OutlineOptions) -> String;
```

## Emission
- Header `%FSLAX46Y46*%` then `%MOMM*%`; `G04 <layer_name>*`; aperture
  `%ADD10C,0.100*%` (0.1 mm stroke) + `D10*`.
- Rectangle corners from `layout.bbox` (`BBox { min: Point2, max: Point2 }`,
  metres) expanded by `margin_m = margin_mm * 1e-3` on each side:
  `(min.x-m, min.y-m)`, `(max.x+m, min.y-m)`, `(max.x+m, max.y+m)`,
  `(min.x-m, max.y+m)`.
- `D02*` move to corner 0, `D01*` draw to corners 1, 2, 3, then `D01*` back to
  corner 0 (explicit close). NO `G36`/`G37` (stroked, not filled).
- Footer `M02*`.
- Reuse the existing private `mm_to_fixed46` / `coord_word` helpers (metres → mm
  → 4.6). If a helper needs a raw `(x_m, y_m)` rather than a `Point2`, add a thin
  private variant; do not change `coord_word`'s existing callers.

## DoD (machine-checkable; pure text, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-export --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-export` exit 0 (existing gerber-001/002 still pass).
4. **`gerber-003` (structure, `tests/`):** for a known `Layout`, the outline output
   starts `%FSLAX46Y46*%` / `%MOMM*%`; contains exactly one `%ADD` aperture; has
   exactly one `D02*` and ≥4 `D01*`; contains NO `G36*`/`G37*`; ends `M02*`.
5. **`gerber-004` (geometry, `tests/`):** parse the `X<int>Y<int>` words and assert
   the four distinct corners equal `bbox.min/max ± margin` (metres→mm→4.6, within
   1e-6 mm). Use a `Layout` whose bbox is known (e.g. `Polygon::rect` traces).

## Out of scope
Drill, multi-copper, mask/silk, KiCad/STEP, non-rectangular outline, any file or
bundle writer (the crate stays pure — callers write the string). No change to
`layout_to_gerber` (copper).
