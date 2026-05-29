# Filter Phase F1.4.0 — `yee-export` Gerber walking skeleton — Design Spec

**ADR:** ADR-0100 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
A new `yee-export` crate emitting a single-copper-layer RS-274X Gerber from a
`yee_layout::Layout` (polygons → G36/G37 filled regions). Walking skeleton —
minimal end-to-end pipe. Pure text, WASM-safe, no FDTD, no native dep.

## Crate (`crates/yee-export/`)
- `Cargo.toml`: workspace-inherited package fields; deps `yee-layout = { workspace
  = true }` + `serde` (only if needed for `GerberOptions`; otherwise omit).
  Lints: `unsafe_code = "forbid"`, `missing_docs = "warn"` (match the repo).
- `src/lib.rs`:
  ```rust
  pub struct GerberOptions { pub layer_name: String }
  impl Default for GerberOptions { fn default() -> Self { Self { layer_name: "F.Cu".into() } } }
  pub fn layout_to_gerber(layout: &yee_layout::Layout, opts: &GerberOptions) -> String;
  ```

## RS-274X emission (single layer)
Header, in order:
- `%FSLAX46Y46*%`  (format: absolute, leading-zero omission, X/Y = 4 integer + 6
  decimal digits)
- `%MOMM*%`        (units = millimetres)
- a layer/function comment is OK (`G04 <layer_name>*`)
- one aperture: `%ADD10C,0.010*%` then `D10*` selected (regions need a current
  aperture selected even though the fill ignores its size)

Per polygon (read the `Layout` polygon vertex list + units from `yee-layout` —
`Point2 { x, y }` in **metres**):
- `G36*`
- move to the first vertex: `X<ix>Y<iy>D02*`
- `X<ix>Y<iy>D01*` for each subsequent vertex (and a closing draw back to the
  first vertex if the polygon isn't already closed)
- `G37*`

Coordinate conversion: metres → mm (`*1e3`) → 4.6 fixed-point integer
(`round(mm * 1e6)`), emitted as a plain signed integer (e.g. `3.0590 mm` →
`3059000` → `X3059000`). Document the conversion inline.

Footer: `M02*`.

## DoD (machine-checkable; pure text, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-export --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-export` exit 0 (sub-second).
4. **`gerber-001` (structure, `tests/`):** build a known small `Layout` (use a
   `yee-layout` generator — e.g. `edge_coupled_bpf` with tiny params, or two
   `Polygon::rect`s if a `Layout` can be hand-built). Assert the output: starts
   with `%FSLAX46Y46*%` then `%MOMM*%`; contains exactly one `G36*` and one
   `G37*` per polygon (count them, equal to `layout` polygon count); defines ≥1
   aperture (`%ADD`); ends with `M02*`.
5. **`gerber-002` (coordinate round-trip, `tests/`):** for one polygon, parse the
   `X<int>Y<int>` integers out of its `G36…G37` region and assert they reproduce
   that polygon's vertices (apply the inverse 4.6/mm/metres conversion) within
   the 1e-6 mm quantisation. Confirms the coordinate emission is correct.

## Workspace wiring (sanctioned cross-lane for a new crate)
Add `"crates/yee-export"` to the root `Cargo.toml` `[workspace] members`. If a
workspace dep entry is wanted for downstream crates, add `yee-export` to
`[workspace.dependencies]` — but do NOT wire any consumer (yee-cli/studio) in
this increment (that is a follow-on).

## Out of scope
KiCad footprint/PCB, STEP/3-D, drill, board outline, multi-layer, mask/silk; any
consumer wiring. WASM-safety: keep it pure (no native/FDTD dep).
