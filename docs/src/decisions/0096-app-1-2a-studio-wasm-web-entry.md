# ADR-0096: App.1.2a — yee-studio wasm32 web entry (eframe WebRunner)

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0089 (filter-design app architecture: one eframe codebase →
desktop + web), ADR-0092 (App.1.0 — `desktop` feature gate), ADR-0095 (App.1.1 —
`--no-default-features` wasm32 build check + CI gate), `FILTER-DESIGN-ROADMAP.md`
(§5a App track)

---

## Context

App.1.0 (ADR-0092) put the `eframe` shell behind a default `desktop` feature so
the egui-free `StudioState` flow compiles for `wasm32-unknown-unknown`, and
App.1.1 (ADR-0095) proved exactly that (`--no-default-features` → wasm32, egui
absent, CI `wasm-build` job). What is *not* yet proven: the **full eframe UI**
(`src/app.rs`'s `StudioApp`) compiling to wasm32 with a browser entry point —
the actual web app, not just the headless logic.

A 2026-05-30 recon recommended splitting App.1.2 in two so the heavy/uncertain
parts are isolated:

- **App.1.2a (this ADR)** — make the full eframe app *compile* for
  `wasm32-unknown-unknown` behind a new `web` feature, with a
  `#[cfg(target_arch = "wasm32")]` `eframe::WebRunner` entry. Gate is a pure
  `cargo check` against the already-installed wasm32 target — **no toolchain
  install**, so it is the light, low-risk half.
- **App.1.2b (deferred)** — the `trunk` bundle (`index.html` / `Trunk.toml`),
  `cargo install trunk`, and static deploy. Heavier (toolchain install) and
  carries the wgpu-WebGL2 runtime risk; it builds on 1.2a.

## Decision

Add a **`web`** Cargo feature to `yee-studio` and a wasm32 browser entry, so the
*same* eframe `StudioApp` runs natively (desktop) and in the browser (web) — one
codebase, per ADR-0089. No `StudioState`/logic changes; this is wiring only.

```toml
[features]
default = ["desktop"]
desktop = ["dep:eframe", "dep:egui", "dep:egui_plot"]
web     = ["dep:eframe", "dep:egui", "dep:egui_plot",
          "dep:wasm-bindgen", "dep:web-sys", "dep:console_error_panic_hook"]
```

`pub mod app` and the eframe/egui/egui_plot deps become gated on
`any(feature = "desktop", feature = "web")`. `src/main.rs` gains a three-way
`cfg` split:

- `all(feature = "desktop", not(target_arch = "wasm32"))` → `eframe::run_native`
  (unchanged native path);
- `all(feature = "web", target_arch = "wasm32")` → a `#[wasm_bindgen(start)]`
  entry running `eframe::WebRunner::new().start("the_canvas_id", WebOptions, …)`
  via `wasm_bindgen_futures::spawn_local`, seeded from the same `default_spec()`;
- otherwise → the existing no-GUI `println!` stub.

eframe's `wgpu` feature is reused for web (wgpu 29 supports
`wasm32-unknown-unknown` with WebGPU → WebGL2 fallback). The implementer is free
to split the eframe dependency per-target (`[target.'cfg(...)'.dependencies]`)
and/or set `default-features = false` with an explicit web-appropriate feature
set if the native default features (x11/wayland/glow) break the wasm32 check —
whatever makes the gate pass, within the `yee-studio` lane.

## Consequences

**Ships:** a `web` feature + wasm32 `WebRunner` entry for `yee-studio`. The
existing native desktop build and the `--no-default-features` headless build are
unchanged (App.1.0/1.1 gates still hold).

**Gate (App.1.2a):** `cargo check -p yee-studio --target wasm32-unknown-unknown
--features web` exits 0 (the full eframe `StudioApp` compiles for the browser).
The wasm32 target is already installed locally and in the App.1.1 CI job. If the
WebGPU bindings require `RUSTFLAGS='--cfg=web_sys_unstable_apis'`, that is
recorded for App.1.2b's `Trunk.toml`. **No new runtime gate** — this is a
compile-only milestone; actually loading the app in a browser is App.1.2b.

**Not in scope:** `trunk`/`index.html`/`Trunk.toml`, `cargo install trunk`, the
static deploy, and any in-browser runtime verification — all App.1.2b. No
`StudioState` or flow-logic change.

**Constraint (ADR-0089):** the light-flow crates (`yee-synth`/`yee-filter`/
`yee-layout`/`yee-plotters`) stay WASM-safe; this ADR touches only `yee-studio`.

---

## References

- eframe 0.34 web support: `eframe::WebRunner` + `#[wasm_bindgen(start)]`
  (eframe_template `src/main.rs`); wgpu 29 wasm32 WebGPU/WebGL2.
- `docs/superpowers/specs/2026-05-30-app-1-2a-studio-wasm-web-entry-design.md`;
  `docs/superpowers/plans/2026-05-30-app-1-2a-studio-wasm-web-entry.md`.
