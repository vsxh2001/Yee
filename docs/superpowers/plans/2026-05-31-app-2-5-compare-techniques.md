# App.2.5 — Compare techniques side-by-side — Plan

**Spec:** `2026-05-31-app-2-5-compare-techniques-design.md` · **ADR:** ADR-0142

## Lane
`crates/yee-studio-web/src/{engine.rs, stages.rs}` (+ `main.rs` only if a signal must
be threaded — not expected). NO yee-filter/physics edits. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-5-compare`.

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/engine.rs` — `design_demo_from(spec, topology)`,
  `design_lumped_from(spec) -> Result<LumpedDesigned, _>`, `design_stepped_from(spec)`,
  `Designed` (`.report`, `.board_size_mm`, `.layout`, `.order()`), `LumpedDesigned`
  (`.verdict`, `.board_size_mm`, `.order()`), `SteppedLowpassDesigned` (`.pass`,
  `worst_*`, `.stopband`, `.board_size_mm`, `.order`); the App.2.4 `verify_view` (the
  same metric-extraction shape — mirror its field reads) + its test (mirror for the new
  test).
- `crates/yee-studio-web/src/stages.rs` — `technique_stage` (where the panel lands;
  it already has the `spec` signal from the guided panel), `guided_panel` (panel/table
  RSX idiom), `technique_status` (RealizationTechnique → live/Soon), `topology_label`,
  `route_into`, `render_recommendation` (the "Use this" button idiom).
- `crates/yee-filter/src/recommend.rs` — `recommend_technique` (to mark the recommended row).
- The spec §Method (the `TechniqueComparison` shape + the per-response technique sets) + ADR-0142.

## Steps
1. `engine.rs`: `TechniqueComparison` + `compare_techniques(spec)` per the spec — key on
   `spec.response`; build each technique's design; pull metrics directly from its graded
   struct (ripple/RL/rejection/pass/order/board size); `realizable=false` on a failed
   dimension (lumped `Err`/`None`, or distributed `layout.is_none()`). Pure; documented.
2. `engine.rs`: the non-vacuous test (spec DoD §1) — band-pass → 3 techniques, metrics
   equal each design's fields, rows not all identical; low-pass → `[SteppedImpedance]`;
   high-pass → `[]`.
3. `stages.rs`: a Compare panel in `technique_stage` — a table over
   `compare_techniques(spec())`, recommended row marked (`recommend_technique`), a "Use
   this" per realizable row (`route_into`), honest single-row / empty-row notes.

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — the new non-vacuous `compare_techniques` test passes;
  existing tests unregressed. Quote "test result: ok".
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0 (the Compare
  panel renders).

Commit on the branch: `yee-studio-web: compare-techniques panel (App.2.5, ADR-0142)` +
the Co-Authored-By trailer.

## Escape hatch
If a metric field isn't where expected, surface it — do NOT fabricate. Keep
`compare_techniques` pure (the component reads `spec()` and passes it). Building three
designs per render is fine (the engines are light/synchronous); if it's a perf concern,
memoize — do NOT cache incorrectly. NEVER edit yee-filter/physics. Blocked > 30 min →
surface.

## Done when
`compare_techniques` + its non-vacuous test are green; the Technique stage shows a real
side-by-side comparison (recommended marked, "Use this" routes); dx build EXIT 0;
existing flows unregressed; clippy/fmt/check clean; diff = `crates/yee-studio-web/src/
{engine.rs, stages.rs}`. Then I verify + review + merge.
