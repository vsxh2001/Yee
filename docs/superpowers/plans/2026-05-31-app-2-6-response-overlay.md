# App.2.6 — Multi-technique response overlay — Plan

**Spec:** `2026-05-31-app-2-6-response-overlay-design.md` · **ADR:** ADR-0143

## Lane
`crates/yee-studio-web/src/{engine.rs, svg.rs, stages.rs}` ONLY. NO yee-filter/physics
edits. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-6-response-overlay`.

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/svg.rs` — `response_plot(sweep, bands) -> String` (MIRROR
  it for the multi-curve `response_overlay`: same viewBox/axis/scaling/mask-band code,
  one polyline per curve + a legend). The `MUTED`/`ACCENT` colour consts.
- `crates/yee-studio-web/src/engine.rs` — `SweepPoint`, each design's `.sweep`
  (`Designed`, `LumpedDesigned`, `SteppedLowpassDesigned`), `design_demo_from` /
  `design_lumped_from` / `design_stepped_from`, `mask_bands`, the App.2.5
  `compare_techniques` + its test (mirror the test shape), the App.2.4 `verify_view`.
- `crates/yee-studio-web/src/stages.rs` — `compare_panel` (where the overlay renders,
  below the table), `synthesis_stage` (the `div { class: "plot", dangerous_inner_html }`
  idiom for an SVG string).
- The spec §Honesty-constraint + §Method (the distinct curves) + ADR-0143.

## Steps
1. `engine.rs`: `OverlayCurve` + `overlay_curves(spec)` per the spec — band-pass → the
   coupled-resonator ideal (label names edge-coupled + hairpin) + the lumped realized
   (`Err`→ not realizable/empty); low-pass → stepped ideal; high-pass → `[]`. Pure, doc.
2. `engine.rs`: the non-vacuous test (spec DoD §1) — 2 curves for band-pass, same grid,
   lumped ≠ coupled-resonator at ≥1 point, each `.sweep` equals the design's; low-pass →
   1; high-pass → `[]`.
3. `svg.rs`: `response_overlay(curves, bands)` — mask bands + one polyline per curve in a
   distinct colour + a legend. (A small palette of ~3 distinct colours.)
4. `stages.rs`: render the overlay in `compare_panel` below the table (the `plot` /
   `dangerous_inner_html` idiom); honest legend; high-pass → no empty chart.

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — the new non-vacuous `overlay_curves` test passes;
  existing tests unregressed. Quote "test result: ok".
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0 (the overlay
  SVG renders in the Compare panel).

Commit on the branch: `yee-studio-web: multi-technique response overlay (App.2.6,
ADR-0143)` + the Co-Authored-By trailer.

## Escape hatch
Keep the labels HONEST — the distributed techniques share one ideal curve; do NOT render
two identical-but-relabelled distributed curves to fake "3 techniques". `overlay_curves`
pure (component reads `spec()` + passes it). If `response_overlay`'s legend/scaling
fights the existing `response_plot` code, surface it rather than duplicating axis logic
incorrectly. NEVER edit yee-filter/physics. Blocked > 30 min → surface.

## Done when
`overlay_curves` + its non-vacuous test green; the Compare panel shows the response
overlay (coupled-resonator ideal + lumped realized for band-pass; stepped for low-pass)
vs the mask with honest legend; dx build EXIT 0; existing flows unregressed;
clippy/fmt/check clean; diff = `crates/yee-studio-web/src/{engine.rs, svg.rs, stages.rs}`.
Then I verify + review + merge.
