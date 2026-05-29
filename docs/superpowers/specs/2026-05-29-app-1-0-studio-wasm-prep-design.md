# App.1.0 — `yee-studio` WASM-prep — Design Spec

**Phase:** App.1.0 · **ADR:** ADR-0092 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
Gate `yee-studio`'s eframe shell behind a default `desktop` Cargo feature so
`StudioState` (the WASM-safe flow logic) compiles with NO eframe/egui in the dep
graph. Natively verifiable App.1 readiness (the `wasm32` target is not installed
yet; the actual WASM build is App.1). `yee-studio` only; no new dependency.

## Changes (`crates/yee-studio` only)
- `Cargo.toml`:
  - Make `eframe`, `egui`, `egui_plot` `optional = true`.
  - `[features]` → `default = ["desktop"]`; `desktop = ["dep:eframe", "dep:egui", "dep:egui_plot"]`.
  - Keep `[[bin]] yee-studio` (its body is cfg-gated; a non-desktop build yields
    a bin with an empty/`#[cfg]`-stubbed `main`, or gate the whole bin — see plan).
- `src/lib.rs`: `#[cfg(feature = "desktop")] pub mod app;`. The `StudioState`,
  `MaskRegionView`, and all flow logic stay un-gated (feature-independent).
- `src/main.rs`: gate the eframe entry behind `#[cfg(feature = "desktop")]`;
  provide a `#[cfg(not(feature = "desktop"))] fn main()` that prints a one-line
  "built without the desktop feature" notice (so the bin target still links).
- Remove the now-resolved `TODO(App.1)` comment on `mod app` (or update it to
  point at App.1's remaining wasm32/trunk work).

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-studio --all-targets -- -D warnings` exit 0 (default features).
3. `cargo clippy -p yee-studio --no-default-features --all-targets -- -D warnings` exit 0.
4. `cargo build -p yee-studio` (default) exit 0 — desktop app builds as in App.0.
5. `cargo build -p yee-studio --no-default-features` exit 0 — and `cargo tree -p
   yee-studio --no-default-features` shows **no `eframe`/`egui`/`egui_plot`** in
   the dep graph (grep the tree output; assert absent).
6. `cargo test -p yee-studio` (default) green; `cargo test -p yee-studio
   --no-default-features` green (the `studio_state_recompute` tests are
   egui-free, so they run under both).

## Out of scope
The `wasm32-unknown-unknown` build, `trunk`/`index.html`, static deploy (App.1).
No behaviour change to the desktop app.
