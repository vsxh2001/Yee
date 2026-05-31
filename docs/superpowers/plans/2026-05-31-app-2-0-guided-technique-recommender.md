# App.2.0 — Guided Technique-Recommender — Plan

**Spec:** `2026-05-31-app-2-0-guided-technique-recommender-design.md` · **ADR:** ADR-0136

## Lane
`crates/yee-filter/**` (engine + gate) + `crates/yee-studio-web/**` (guided UI).
Out of lane → finding, do NOT fix in place.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-0-recommender`.

## Pattern files (READ FIRST)
- `crates/yee-filter/src/lib.rs` — `FilterSpec` / `Response` / `Topology` / `synthesize`.
  Match its doc style (`#![warn(missing_docs)]`, every public item documented) and its
  inline `#[cfg(test)] mod tests` style for the gate (or a `tests/` integration file).
- `crates/yee-studio-web/src/stages.rs` — `technique_stage` (the gallery card pattern,
  RSX style, the `Topology` enum at line ~24) + `main.rs` routing (`topology` Signal,
  `Stage::rail`, the `Stage::Technique => technique_stage(...)` arm).
- The vision doc §4–5 and ADR-0136 (the decision tree + the gate table — implement
  exactly those thresholds).

## Steps
1. **Engine (yee-filter):** add `RealizationTechnique`, `TechniqueRecommendation`,
   `recommend_technique(&FilterSpec)` per the spec's decision tree. Re-export from the
   crate root. Document every public item. Rationale strings name the deciding factor.
2. **Gate (yee-filter):** the canonical spec→technique table (spec DoD §1) as a test —
   each case asserts the exact expected technique; plus non-empty-rationale and
   primary-∉-alternatives invariants. Must be NON-VACUOUS (a constant recommender fails).
3. **UI (yee-studio-web):** a Guided panel atop `technique_stage`: a form (response,
   f0/cutoff, fbw, optional stopband target) → "Recommend" → calls `recommend_technique`,
   renders the primary (highlighted) + rationale + ranked alternatives. Map
   `RealizationTechnique`→`Topology`; "Use this" routes live techniques (set the topology
   Signal + jump to Spec); Soon techniques show honestly + offer the nearest live one.
   Seed the editable `FilterSpec` from the form. Keep the expert gallery below.

## Verify (run these; expected EXIT 0)
- `cargo test -p yee-filter` — the gate passes (yee-filter is light/pure — host is fine,
  no Docker box needed). Quote the "test result: ok" line for the gate test.
- `cargo clippy -p yee-filter --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- Studio: install `cargo install dioxus-cli --version 0.6.3 --locked` + `rustup target
  add wasm32-unknown-unknown`, then `dx build --platform web --release` in
  `crates/yee-studio-web` → EXIT 0; confirm the Guided panel renders (grep the built
  `target/dx/yee-studio-web/release/web/public/index.html` for the recommend-form text,
  or describe the RSX).

Commit on the branch: `yee-filter: recommend_technique guided dual-UI engine + gate
(App.2.0, ADR-0136)` and a second `yee-studio-web: guided technique-recommender panel
(App.2.0)` + the Co-Authored-By trailer.

## Escape hatch
If the Dioxus UI blocks > 30 min (RSX/signal wiring, dx build), **land the engine + gate
(cargo test -p yee-filter green) and surface the UI as a follow-on** — the validatable
core ships; do NOT fake the UI and do NOT weaken the gate. If a canonical gate case
seems wrong, surface it (don't silently change the expected technique to pass).

## Done when
`recommend_technique` + the non-vacuous gate are green; the studio dx-builds with the
Guided panel rendering the recommendation; workspace check/clippy/fmt clean; diff =
`crates/yee-filter/**` + `crates/yee-studio-web/**`. Then I (dispatcher) verify + review.
