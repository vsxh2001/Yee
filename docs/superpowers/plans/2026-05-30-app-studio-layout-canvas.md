# App — yee-studio layout preview canvas — Plan

**Spec:** `2026-05-30-app-studio-layout-canvas-design.md` · **ADR:** ADR-0101

## Lane
`crates/yee-studio/**` ONLY (`src/lib.rs`, `src/app.rs`, tests; no `Cargo.toml`
change needed — `yee-layout`/`egui_plot` already deps). Do NOT edit
`yee-filter`/`yee-layout`/other crates. Out of lane → finding. Keep `StudioState`
egui-free + WASM-safe (the new `layout` is `Result<yee_layout::Layout, String>`).

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/app-studio-layout-canvas`.

## Pattern files
- `crates/yee-studio/src/lib.rs` — `StudioState`, `apply_derived` (where `dims`
  is computed — add `layout` right beside it, same idiom), `from_spec` (init).
- `crates/yee-studio/src/app.rs` — `show_response_plot` (the existing
  `egui_plot::Plot` + `Polygon::new(name, PlotPoints::from(box_pts))` mask-box
  pattern to mirror) and `show_dimensions` (panel-invocation idiom).
- `crates/yee-filter/src/dimension.rs` — `dimension_edge_coupled_layout` signature
  + return type.
- `crates/yee-layout/src/lib.rs` — `Layout.traces: Vec<Polygon>`, `Polygon.verts:
  Vec<Point2>`, `Point2 { x, y }` (metres).

## Steps
1. `lib.rs`: add `layout` derived field; init in `from_spec`; compute in
   `apply_derived` via `dimension_edge_coupled_layout`.
2. `app.rs`: `show_layout` fn (Plot, data_aspect 1.0, one Polygon per trace, mm);
   wire into the central panel.
3. test: `studio_state_layout` per DoD 4.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-studio --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-studio --jobs 2
nice -n 19 cargo check -p yee-studio --no-default-features --target wasm32-unknown-unknown --jobs 2
nice -n 19 cargo tree -p yee-studio --no-default-features --target wasm32-unknown-unknown -i egui  # egui MUST be ABSENT
```
The native eframe build is cached from App.1.2a; incremental. Do NOT run
`cargo test --workspace`, FDTD, mom-001.

## Escape hatch
Blocked > 15 min — the `egui_plot` 0.35 `Polygon`/`Plot` API for drawing the
layout fights you (fill/aspect), OR `StudioState` can't hold `layout` without an
egui/native type → STOP and surface (the API issue + what you tried, or the
WASM-safety blocker). Fallback within budget: draw `Line` outlines instead of
filled `Polygon`. Do NOT weaken DoD 5 (egui must stay absent from the headless
wasm tree). Do NOT edit yee-filter/yee-layout.

## Done when
DoD 1–5 pass; `git diff --stat <base>..HEAD` = only `crates/yee-studio/**`
(+ docs); `StudioState` still compiles egui-free for `--no-default-features` wasm32.
