# App.D.1L — lumped-LC studio stages on the Dioxus shell — Plan

**Spec:** `2026-05-30-app-d1-lumped-studio-stages-design.md` · **ADR:** ADR-0120

## Lane
`crates/yee-studio-web/**` ONLY (`src/{main,stages,engine,svg}.rs`,
`assets/studio.css`, the crate's `Cargo.toml` to add `yee-filter` dep if not
present). Do NOT edit `yee-filter`/`yee-synth`/`yee-layout`/`yee-export`/
`yee-studio` or `StudioState` — call into them read-only. Out of lane → finding.

## Base
Work in the existing worktree `worktrees/studio-dx` on branch
`feature/app-d0-dioxus-poc` (the POC, tip `b8cfb90`). Re-fetch and rebase the
branch on current `main` if it lags (it predates the F2.x merges — ensure
`yee-filter`'s lumped API is available). DO NOT merge to main (maintainer-gated).

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/main.rs` — Shell A (top bar + stage rail + canvas),
  stage routing. Add the lumped stages to the rail/router.
- `crates/yee-studio-web/src/stages.rs` — the existing Synthesis + Layout stage
  renderers (RSX + design-system classes). MIRROR their style for the lumped
  stages; do not invent a new look.
- `crates/yee-studio-web/src/engine.rs` — the existing engine adapter (drives
  `yee_synth`/`yee_filter` on the fixture). ADD the lumped adapter calls here.
- `crates/yee-studio-web/src/svg.rs` — the |S21|-vs-mask plot + board SVG helpers.
  Reuse the plot for `ladder_s21`; add a lumped-board SVG.
- `crates/yee-studio-web/assets/studio.css` — the design-system tokens. Reuse;
  extend only if a new component needs it.
- Engine APIs (read their signatures/docs): `yee_filter::{synthesize_lumped,
  LumpedLadder, LcResonator, select_components, ESeries, Bom, monte_carlo_yield,
  YieldResult, lumped_board, LumpedBoard, ladder_s21}` and the spec mask type.

## Steps
1. `engine.rs`: lumped adapter — from the fixture/spec, build `LumpedLadder`, `Bom`
   (E24 + E96), `YieldResult` (seed fixed, n≈500), `LumpedBoard`. Return plain
   view structs.
2. Technique stage: Lumped LC live + selectable → sets the topology, routes
   downstream stages to the lumped renderers.
3. `stages.rs` + `svg.rs`: the four lumped renderers (Synthesis table + ladder_s21
   SVG/mask/PASS-FAIL; Components+BOM table E24/E96; Tolerance yield card; Layout
   board SVG + table) — mirror the existing stage style.
4. Build for web + serve; eyeball the polish.

## Verify (this crate is WASM — DO NOT run heavy cargo; light checks only)
- `cargo fmt --check -p yee-studio-web` → exit 0.
- `cargo clippy -p yee-studio-web --target wasm32-unknown-unknown -- -D warnings`
  (or the POC's documented check command) → exit 0.
- Web build: `dx build --platform web --release` (dx 0.6.3, as the POC used) →
  succeeds; static bundle produced.
- Serve the bundle locally (python http.server on a port) and confirm the lumped
  stages render REAL values (the actual ladder L/C, the real BOM rows for E24 and
  E96, the real yield %, the real board SVG) — not placeholders. Capture the served
  URL + a note of the rendered numbers in the report.

## Escape hatch
If a lumped engine API is missing/awkward for the view (e.g. a field not exposed),
surface it as a finding (do NOT edit `yee-filter` — that's out of lane; define a
thin view-side adapter instead). If the wasm build breaks on a dep that is not
WASM-safe, STOP and surface which dep (the lumped-LC crates are pure-math/WASM-safe
by design — a break means an accidental non-WASM pull-in). Blocked > 60 min →
surface. Do NOT merge to main; do NOT touch engine crates or eframe.

## Done when
DoD 1–4: wasm build green, fmt/clippy clean, the four lumped stages render real
engine output, served locally; diff = `crates/yee-studio-web/**` only; branch stays
`feature/app-d0-dioxus-poc` (unmerged). Report the served URL + rendered numbers.
