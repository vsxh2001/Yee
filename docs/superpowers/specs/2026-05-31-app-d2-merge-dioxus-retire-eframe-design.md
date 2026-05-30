# App.D.2 — merge the Dioxus studio + retire the eframe view — Design Spec

**ADR:** ADR-0130 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

The maintainer approved merging the Dioxus studio (the polished-UI component) +
retiring the eframe `yee-studio`. The Dioxus studio (`yee-studio-web`, branch
`feature/app-d0-dioxus-poc`) has the full real lumped-LC flow + POC distributed
stages; it must land on `main` as THE studio, with CI wired and eframe retired.

## Goal

`yee-studio-web` on `main` as the studio (real lumped flow + the distributed
Synthesis/Layout; Spec/Technique/Export real-enough), CI green (wasm build +
fmt/clippy), the eframe view retired, `StudioState`/engine intact. Reviewer before
merge.

## Architecture / scope

- **Stage build-out (yee-studio-web):** Spec → a real spec input form (drives
  `synthesize`); Technique → topology gallery (Lumped LC live; distributed
  selectable or honestly "Soon"); Export → param-sheet + download affordances.
  Match the existing design system + SVG idiom. (Lumped Synthesis/Components/BOM/
  Tolerance/Layout already real — ADR-0120.) A *fine* distributed polish pass is
  out of scope (incremental after merge); the bar is "shippable, real-where-it
  claims-to-be, honest stubs labelled."
- **CI:** add/extend a job so `yee-studio-web` builds for `wasm32-unknown-unknown`
  + `cargo fmt --check`/`clippy -D warnings`. The existing `wasm-build` job targets
  `yee-studio` (eframe) — repoint or add `yee-studio-web`. Keep the workspace +
  docs jobs green.
- **Retire eframe:** remove `crates/yee-studio`'s eframe app (`app.rs` + the eframe
  `main` / `desktop`-feature render path), keeping `StudioState` + the egui-free
  core + the engine crates (`yee-studio-web` + the validation paths reuse them). If
  a clean delete risks orphaning refs in one pass, **deprecate** (gate off) +
  delete in a tight follow-up. No engine/physics change.
- **Merge** `--no-ff` to `main` after review.

## DoD (machine-checkable)

1. `cargo fmt --check --all` + `cargo clippy --workspace --all-targets -- -D
   warnings` (or the no-default-features CI variant) exit 0.
2. `yee-studio-web` builds for `wasm32-unknown-unknown`; the CI job is wired.
3. The workspace builds + the docs (mdBook) build green; no orphaned `yee-studio`
   refs after retirement (`cargo check --workspace`).
4. The Dioxus studio serves locally with the real lumped flow rendering real engine
   data (ladder/BOM/yield/board) — a smoke check.
5. Code-review APPROVED (no orphaned refs, workspace green, lumped flow real);
   merged `--no-ff` to `main`.

## Out of scope

EM-Verify stage (Track A); desktop/webview packaging; mobile; a fine distributed
polish pass.

## Why

It ships the goal's polished-UI component: a pure-Rust web-first studio holding the
real lumped-LC journey on `main`, replacing the chunky eframe tool — the maintainer's
approved direction.
