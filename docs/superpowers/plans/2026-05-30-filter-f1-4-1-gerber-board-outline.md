# Filter Phase F1.4.1a — `yee-export` Gerber board outline — Plan

**Spec:** `2026-05-30-filter-f1-4-1-gerber-board-outline-design.md` · **ADR:** ADR-0103

## Lane
`crates/yee-export/**` ONLY (`src/lib.rs`, `tests/`). Do NOT edit `yee-layout`/
other crates — consume `yee_layout::{Layout, BBox, Point2}`. Out of lane →
finding. Keep `yee-export` WASM-safe (pure `String`; no fs / native dep).

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-f1-4-1-gerber-outline`.

## Pattern files
- `crates/yee-export/src/lib.rs` — the F1.4.0 `layout_to_gerber` + `GerberOptions`
  + the private `mm_to_fixed46`/`coord_word` helpers + module coordinate-model
  docs to mirror. `tests/gerber_001_structure.rs` + `tests/gerber_002_roundtrip.rs`
  — the test idiom (build a `Layout`, assert structure + parse coords back).
- `crates/yee-layout/src/lib.rs` — `BBox { min: Point2, max: Point2 }`,
  `Point2 { x, y }`, `Layout { bbox, traces, .. }`, `BBox::from_polygons`,
  `Polygon::rect`.

## Steps
1. `src/lib.rs`: add `OutlineOptions` (+ `Default`) and `layout_to_gerber_outline`
   per the spec (stroked rectangular contour from `bbox ± margin`). Reuse
   `mm_to_fixed46`; add a private `coord_word_xy(x_m, y_m)` if needed (don't break
   `coord_word`). Document every public item.
2. `tests/gerber_003_outline_structure.rs` + `tests/gerber_004_outline_geometry.rs`
   per DoD 4–5.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-export --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-export --jobs 2
```
Pure text — sub-second. Do NOT run `cargo test --workspace`, FDTD, mom-001.

## Escape hatch
Blocked > 15 min — uncertainty whether Edge.Cuts should be stroked vs region (it
is STROKED — a contour, not a `G36` fill), or a `bbox` accessor is missing →
STOP and surface. Do NOT region-fill the outline; do NOT add a file writer; do
NOT edit yee-layout.

## Done when
DoD 1–5 pass; `git diff --stat <base>..HEAD` = only `crates/yee-export/**` + the
3 committed docs; `yee-export` still pure/WASM-safe (no fs/native dep).
