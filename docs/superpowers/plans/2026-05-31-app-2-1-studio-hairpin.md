# App.2.1 — Light the Hairpin technique in the studio — Plan

**Spec:** `2026-05-31-app-2-1-studio-hairpin-design.md` · **ADR:** ADR-0138

## Lane
`crates/yee-studio-web/**` ONLY. NO yee-filter / yee-layout edits (the
`dimension_hairpin` / `dimension_hairpin_layout` engine already exists). Out of lane →
finding, do NOT fix.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-1-studio-hairpin`.

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/engine.rs` — `Designed`, `design_demo_from`,
  `derive_geometry`, the `Geometry` bundle, the edge-coupled `dimension_edge_coupled*`
  calls, `ResonatorRow`, and the in-crate `#[cfg(test)]` engine tests. **Mirror the
  edge-coupled geometry path for hairpin.**
- `crates/yee-studio-web/src/stages.rs` — `Topology` (line ~24), `Stage::rail`,
  `technique_stage` (the card list + `selects`), `technique_status` (App.2.0
  recommender map), `topology_label`, `layout_stage` (the board + resonator table).
- `crates/yee-studio-web/src/main.rs` — the `designed` memo / `use_effect` that calls
  `design_demo_from` (thread `topology()` in as a dependency), and the
  `Stage::Layout`/`Synthesis` routing arms (`lumped_flow` branch).
- `crates/yee-filter/src/dimension.rs` — `dimension_hairpin`, `dimension_hairpin_layout`,
  `HairpinDimensions` (the engine you are surfacing — read its fields: `line_width_m`,
  `arm_length_m`, `gaps_m`, `fold_spacing_m`).
- The spec §Method (the 6 steps) + ADR-0138.

## Steps (spec §Method order)
1. `Topology::Hairpin` + `Stage::rail(Hairpin) => &Stage::DISTRIBUTED`.
2. `design_demo_from(spec, topology)` + `design_demo()` updated; `derive_geometry(project,
   topology)` branches edge-coupled vs `dimension_hairpin*`. Shared fields computed once.
   Hairpin populates `layout`, `board_size_mm`, `line_eps_eff`, `dim_error`, and a
   hairpin `resonators` table from `HairpinDimensions`.
3. `main.rs` memo depends on `(spec, topology)`.
4. Hairpin gallery card `selects: Some(Topology::Hairpin)`.
5. `technique_status`: `Hairpin => Live(Topology::Hairpin)`; `topology_label` Hairpin arm.
6. Layout/Export render the hairpin board (generic `Layout` — should need no change) +
   a topology-aware table label.

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — incl. the NEW non-vacuous host test (spec DoD §3:
  `design_demo_from(demo_spec(), Topology::Hairpin)` → `Some(layout)` AND its layout
  differs from the edge-coupled layout for the same spec). Existing tests still pass.
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cargo install dioxus-cli --version 0.6.3 --locked` (if needed) + `rustup target add
  wasm32-unknown-unknown`; `cd crates/yee-studio-web && dx build --platform web --release`
  → EXIT 0. (Studio is light to build; no Docker box, no FDTD.)

Commit on the branch: `yee-studio-web: light the Hairpin technique (App.2.1, ADR-0138)`
+ the Co-Authored-By trailer.

## Escape hatch
If the hairpin `Layout` does NOT render generically (the SVG path assumes edge-coupled
geometry) or the `ResonatorRow` table is too edge-coupled-specific to reuse, SURFACE the
specific blocker — do NOT fake a hairpin render with edge-coupled geometry, do NOT stub
the layout. A lit Hairpin card MUST route to the real `dimension_hairpin` output. If
threading topology into the memo causes a Dioxus borrow/stale-signal issue you can't
resolve in 30 min, surface it. Blocked > 30 min → stop and surface.

## Done when
Hairpin is a live, routable technique driven by the real `dimension_hairpin` engine; the
non-vacuous host test (hairpin layout ≠ edge-coupled layout) passes; dx build EXIT 0;
clippy/fmt/check clean; edge-coupled + lumped flows unregressed; diff =
`crates/yee-studio-web/**`. Then I (dispatcher) verify + review + merge.
