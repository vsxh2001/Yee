# ADR-0092: App.1.0 — `yee-studio` WASM-prep (gate the eframe shell behind a feature)

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0090 (App.0 yee-studio), ADR-0089 (app architecture — WASM-safety)

---

## Context

App.0 shipped `yee-studio` with `StudioState` (egui-free flow logic) + an `app`
module (eframe shell). ADR-0089 requires the light flow to compile to WASM for
the web app (App.1). App.0 left a `TODO(App.1)`: `pub mod app` is unconditional,
so the lib pulls `eframe`/`wgpu` even for a consumer that only wants
`StudioState`. Before the full WASM build (trunk/wasm-pack — and the `wasm32`
target is not yet installed here), the cheap, natively-verifiable step is to
**gate the eframe shell behind a Cargo feature** so `StudioState` builds without it.

## Decision

In `yee-studio`:
- Make `eframe`, `egui`, `egui_plot` **optional** deps; add
  `[features] default = ["desktop"]`, `desktop = ["dep:eframe", "dep:egui",
  "dep:egui_plot"]`.
- Gate the shell: `#[cfg(feature = "desktop")] pub mod app;` and gate the
  `[[bin]] yee-studio` / `src/main.rs` body behind `#[cfg(feature = "desktop")]`
  (a non-`desktop` build has no windowed binary).
- `StudioState` + `MaskRegionView` + `voxelize`-free flow logic stay
  **always-compiled** (feature-independent), eframe-free.

This makes `cargo build -p yee-studio --no-default-features` compile only the
WASM-safe flow logic — provable WASM-readiness without the wasm32 toolchain.
The full WASM target build + `trunk` deploy remain App.1 (when the toolchain is
installed).

## Consequences

**Ships:** the `desktop` feature + cfg-gating; default build is byte-identical to
App.0 (desktop on by default). Gate: `cargo build -p yee-studio` (default) builds
the app as before; `cargo build -p yee-studio --no-default-features` builds
`StudioState` with NO eframe/egui in the dep graph; `cargo test -p yee-studio`
(default) still green; the `studio_state_recompute` tests run under
`--no-default-features` too (they don't touch egui). clippy `-D warnings` clean
on both feature sets.

**Not in scope:** the actual `wasm32-unknown-unknown` build, `trunk`/`index.html`,
and static deploy (App.1 — needs the wasm toolchain installed). No new dependency.

---

## References
- ADR-0090 (`TODO(App.1)` in `yee-studio/src/lib.rs`), ADR-0089.
- `docs/superpowers/specs/2026-05-29-app-1-0-studio-wasm-prep-design.md`;
  `docs/superpowers/plans/2026-05-29-app-1-0-studio-wasm-prep.md`.
