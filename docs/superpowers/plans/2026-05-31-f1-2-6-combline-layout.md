# Filter F1.2.6 — Combline board layout — Plan

**Spec:** `2026-05-31-f1-2-6-combline-layout-design.md` · **ADR:** ADR-0145

## Lane
`crates/yee-filter/**` ONLY (compose `yee_layout` primitives — do NOT edit yee-layout).
Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/f1-2-6-combline-layout`.

## Pattern files (READ FIRST — all under the WORKTREE `crates/...`, edit there, not main)
- `crates/yee-filter/src/dimension.rs` — `dimension_stepped_impedance_layout` (THE
  primitive-composition pattern to mirror: `Polygon::rect(x0,y0,w,h)`, `PortRef{at,
  width_m, ref_impedance_ohm}`, `Point2::new`, `BBox::from_polygons`, `Layout{substrate,
  traces, ports, bbox}`); `dimension_combline` + `ComblineDimensions` (line_width_m,
  theta0_rad, resonator_length_m, loading_cap_f, gaps_m, target_k) — call it for the dims;
  `dimension_hairpin_layout` (the feed/port idiom + the uniform-gap doc note pattern).
- `crates/yee-layout/src/lib.rs` — `Polygon` (`rect`), `PortRef`, `Point2`, `BBox`
  (`from_polygons`), `Layout`, `Substrate`. (Reference only — do NOT edit.)
- The spec §Method (the comb geometry) + ADR-0145.

## Steps
1. `dimension_combline_layout(project, theta0_rad, substrate) -> Result<Layout, DimError>`:
   `let dims = dimension_combline(project, theta0_rad, substrate)?;` then compose:
   - N resonator-line rects (vertical: x = running offset, y from 0 to resonator_length_m,
     w = line_width_m, h = resonator_length_m). x_i = Σ_{j<i}(line_width_m + gaps_m[j]).
   - a ground-spine rect at the short-circuit end spanning the x-range.
   - N cap-pad rects at the open ends (≈ line_width_m squares).
   - input/output feed rects + 2 `PortRef`s (ref impedance = project.spec.z0_ohm), neutral
     feed length/width defaults (mirror stepped-Z / hairpin).
   - `bbox = BBox::from_polygons(&traces)`. Re-export from the crate root. Document.
2. Gate `crates/yee-filter/tests/dim_combline_layout_001.rs` (spec §DoD): build the demo
   order-5 combline layout (θ0=π/4); assert N resonator lines with dims == dimension_combline
   (line_width × resonator_length), consecutive resonator x-pitch = line_width + gaps_m[i]
   (monotone), exactly 2 ports (ref z0), bbox positive/finite, all rects positive extent.
   Pattern the test on an existing geometry/layout test (geo-003 / a hairpin layout test).

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-filter --test dim_combline_layout_001` — quote "test result: ok" + the
  resonator count + the x-pitches vs (line_width + gaps).
- `cargo test -p yee-filter` (full crate, no regressions).
- `cargo clippy -p yee-filter --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
(yee-filter is light/pure — host fine; NO Docker box, NO FDTD.)

Commit: `yee-filter: combline board layout generator + geometry gate (F1.2.6, ADR-0145)`
+ the Co-Authored-By trailer.

## Escape hatch
Compose the comb HONESTLY (aligned grounded resonator lines + cap pads) — do NOT reuse
`edge_coupled_bpf` (staggered open lines, wrong for combline). If composing the spine/pads
is fiddly, the resonator lines + feeds + ports + a spine are the must-haves (cap pads can
be minimal). Do NOT edit yee-layout (compose primitives in yee-filter). Do NOT fake the
geometry or weaken the per-section-gap assertion. Blocked > 30 min → surface.

## Done when
`dimension_combline_layout` + the non-vacuous geometry gate are green; clippy/fmt/check
clean; diff = `crates/yee-filter/**`. Then I (dispatcher) verify + adversarial-review +
merge. Studio lighting (App.2.7) is the next increment.
