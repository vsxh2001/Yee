# App.D.1L — lumped-LC studio stages on the Dioxus shell — Design Spec

**ADR:** ADR-0120 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

The lumped-LC engine is shipped (F2.0 synth, F2.1 BOM, F2.4 tolerance, F2.2 board)
but invisible — no UI. The Dioxus POC renders only the distributed Synthesis +
Layout. The goal wants a **polished UI** surfacing **component-choosing, BOM,
tolerance**. Build the lumped-LC stage set on the POC shell.

## Goal

A slick, real-engine lumped-LC flow in the Dioxus studio (POC branch): Technique →
Synthesis → Components+BOM → Tolerance → Layout, in the existing design system,
SVG graphics, honest ideal-vs-realized labelling. Engine crates untouched.

## Architecture

Reuse the POC's structure (`crates/yee-studio-web/src/{main,stages,engine,svg}.rs`
+ `assets/studio.css`). The view layer calls the shipped engine:

- `engine.rs` gains a lumped-LC adapter: from the committed Chebyshev N=5 fixture
  (or the studio spec), call `yee_filter::synthesize_lumped` → `LumpedLadder`,
  `select_components(&ladder, ESeries)` → `Bom`, `monte_carlo_yield(...)` →
  `YieldResult`, `lumped_board(...)` → `LumpedBoard`. Pure data structs handed to
  the view (WASM-safe — these are the light-flow crates, no FDTD).
- `stages.rs` gains the lumped stage renderers; `svg.rs` gains the lumped board
  SVG (footprints/pads/traces) and reuses the |S21|-vs-mask plot for `ladder_s21`.
- Technique stage: a topology gallery where **Lumped LC** is live (selecting it
  sets a `Topology::LumpedLc` in the studio state and routes downstream stages).

### Stage content
1. **Technique** — Lumped LC selectable (live); distributed entries stay greyed
   "Soon". Selection drives the flow.
2. **Synthesis (lumped)** — resonator table (index, series/shunt, L [nH], C [pF]) +
   ideal `ladder_s21` |S21| SVG vs mask + PASS/FAIL chip.
3. **Components + BOM** — E24/E96 toggle; BOM table (ref, kind L/C, ideal value,
   chosen E-series value, deviation %, tolerance %, qty). Deviation colour-coded.
4. **Tolerance / yield** — run `monte_carlo_yield` (fixed seed, n≈500); show yield
   % (E24 vs E96), worst-case RL + rejection, and the honest "narrowband lumped
   needs tight parts / looser spec" note when yield is low.
5. **Layout (lumped)** — `lumped_board` SVG (footprints to scale, pads, signal
   trace, ground) + dimension callouts + placement table.

## DoD (machine-checkable)

1. `crates/yee-studio-web` builds for `wasm32-unknown-unknown`; `cargo fmt --check`
   + `cargo clippy -p yee-studio-web --target wasm32-unknown-unknown -- -D warnings`
   (or the crate's existing check) clean.
2. The lumped stages render **real** engine output: the actual `LumpedLadder`
   values, the real `ladder_s21` curve + PASS/FAIL, the real `Bom` (E24 *and* E96),
   the real `YieldResult`, the real `lumped_board` SVG — not placeholders.
3. Builds a static web bundle served locally (a URL the maintainer opens); looks
   slick (design system, SVG, refined type).
4. Engine crates (`yee-filter`, `yee-synth`, `yee-layout`, `yee-export`) and
   `StudioState`/`yee-studio` (eframe) **unchanged**. Work stays on the POC branch.

## Out of scope

Merge to main / eframe retirement (maintainer-gated); EM-Verify stage (Track A);
desktop/webview; distributed-flow stage build-out; CI wiring; mobile.

## Why

It makes the goal's "polished UI + component-choosing + BOM + tolerance" real and
reviewable from the shipped engine, completing the lumped journey in the browser
(minus EM-Verify, which Track A unblocks) — and gives the maintainer a concrete,
fuller build to judge per ADR-0110.
