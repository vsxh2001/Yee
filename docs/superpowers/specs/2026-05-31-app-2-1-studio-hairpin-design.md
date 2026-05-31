# App.2.1 — Light the Hairpin technique in the studio — Design Spec

**ADR:** ADR-0138 · **Date:** 2026-05-31 · **Status:** Accepted
**Vision:** `2026-05-31-ideal-filter-design-app-vision.md` §5 (band-pass breadth — fill
the gallery). Follows ADR-0136 (the guided recommender recommends Hairpin) and F1.2.2
(the `dimension_hairpin` engine already exists, un-surfaced).

## Problem

The studio's `dimension_hairpin` / `dimension_hairpin_layout` engine shipped in F1.2.2,
but the studio gallery still greys **Hairpin** as "Soon" — the engine is un-surfaced.
The guided recommender (App.2.0) recommends Hairpin for moderate-bandwidth band-pass but
can only route it to the edge-coupled stand-in. Lighting Hairpin fills a gallery card and
makes the recommendation routable to its real dimensioner.

## Key insight (why this is low-risk)

A hairpin filter is the **same coupled-resonator band-pass synthesis** as edge-coupled
(identical prototype, coupling matrix, swept ideal response, PASS/FAIL mask verdict) —
**only the physical realization differs** (U-folded λ/4 arms vs straight λ/2 lines). So
everything in `Designed` except the *geometry-derived* fields (layout, resonator table,
board size, line ε_eff, dim_error) is topology-independent and already correct for
hairpin. The board renders generically from `yee_layout::Layout` (SVG polygons), so the
existing layout/export rendering already handles a hairpin `Layout`. The only new work is
**branching the geometry derivation on topology** and the **routing/UI**.

## Method

All in `crates/yee-studio-web` (lane). Mirror the edge-coupled path.

1. **`Topology::Hairpin`** added to the studio enum (stages.rs). `Stage::rail(Hairpin)`
   = `DISTRIBUTED` (same 6 stages as edge-coupled — Spec/Technique/Synthesis/Layout/
   Verify/Export).
2. **Engine (engine.rs):** thread the topology into the geometry derivation. Change
   `design_demo_from(spec)` → `design_demo_from(spec, topology)` (and `design_demo()` →
   `design_demo_from(demo_spec(), Topology::EdgeCoupled)`); `derive_geometry(project,
   topology)` branches: `EdgeCoupled`/`LumpedLc`-distributed-fallback → the existing
   `dimension_edge_coupled*`; `Hairpin` → `dimension_hairpin` + `dimension_hairpin_layout`.
   The shared fields (project, g_values, coupling, sweep, mask_bands, report) are computed
   once, topology-independent. The hairpin geometry populates `layout` (the generic
   `Layout`), `board_size_mm`, `line_eps_eff` (at the hairpin line width), `dim_error`
   (propagating `DimError`), and a hairpin-appropriate `resonators` table (arm length
   `≈λg/4`, gaps, line width — from `HairpinDimensions`). Keep `LumpedLc` routing to the
   lumped engine unchanged.
3. **Reactive memo (main.rs):** the `designed` memo recomputes on **(spec, topology)** —
   add `topology()` as a dependency so switching technique re-derives the geometry.
4. **Technique gallery (stages.rs):** the Hairpin card gets `selects: Some(Topology::Hairpin)`.
5. **Recommender panel (stages.rs `technique_status`):** `RealizationTechnique::Hairpin`
   → `TechStatus::Live(Topology::Hairpin)` (was `Soon(EdgeCoupled)`); `topology_label`
   gains a Hairpin arm.
6. **Layout / Export stages:** render the hairpin board (generic `Layout` SVG — should
   need no change) + a topology-aware resonator table label (arm-length vs resonator-
   length). Export emits the hairpin Gerber/KiCad from the real hairpin `Layout`.

## Changes

- `crates/yee-studio-web/src/{stages.rs, engine.rs, main.rs}` (+ `svg.rs` only if the
  hairpin `Layout` needs a render tweak — it should not). Update the in-crate engine
  tests for the `design_demo_from` signature.
- NO yee-filter / yee-layout change (the hairpin engine already exists).

## DoD (machine-checkable)

1. `dx build --platform web --release` (working-dir `crates/yee-studio-web`) EXIT 0.
2. The built wasm (`target/dx/yee-studio-web/release/web/public/assets/*.wasm`) is the
   served bundle; the Hairpin technique is **live + routable** — verify the card no
   longer renders as "Soon" for Hairpin (the routing path exists) by confirming the
   built bundle + a host unit test (below) exercise the hairpin geometry.
3. `cargo test -p yee-studio-web` green — including a NEW host test asserting
   `design_demo_from(demo_spec(), Topology::Hairpin)` returns a `Some(layout)` (the demo
   spec dimensions as a hairpin) **and** that the hairpin layout differs from the
   edge-coupled layout for the same spec (the card routes to the REAL hairpin
   dimensioner, not a stub — non-vacuous). The existing edge-coupled / lumped engine
   tests still pass (no regression).
4. `cargo clippy -p yee-studio-web --all-targets -- -D warnings` + `cargo fmt --check`
   clean; `cargo check --workspace` green.

## Out of scope

Combline / Interdigital (no engine yet); stepped-Z low-pass studio path (needs a
low-pass response path — separate); the hairpin `qe`→tap-offset refinement (F1.2.1,
deferred in the engine). EM verify (ADR-0133 wall) untouched.

## Why

Fills a gallery card with a real, already-validated dimensioner; makes the recommender's
Hairpin recommendation routable; advances the *visible* end-to-end app — low-risk because
the synthesis is shared and the board renders generically.
