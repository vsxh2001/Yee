# App.1.2a — yee-studio wasm32 web entry (eframe WebRunner) — Design Spec

**Phase:** App.1.2a · **ADR:** ADR-0096 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

Make the **full** `yee-studio` eframe UI (`StudioApp`, `src/app.rs`) *compile*
for `wasm32-unknown-unknown` behind a new `web` feature, with a
`#[cfg(target_arch = "wasm32")]` `eframe::WebRunner` browser entry. Compile-only:
the `trunk` bundle + deploy + in-browser run are App.1.2b. Wiring only — no
`StudioState`/flow-logic change. One eframe codebase → desktop + web (ADR-0089).

## Changes (`crates/yee-studio/**` ONLY)

### `Cargo.toml`
- Add a `web` feature mirroring `desktop`'s eframe/egui/egui_plot plus the web
  deps:
  ```toml
  web = ["dep:eframe", "dep:egui", "dep:egui_plot",
         "dep:wasm-bindgen", "dep:web-sys", "dep:console_error_panic_hook"]
  ```
- Add the new optional deps: `wasm-bindgen`, `web-sys`,
  `console_error_panic_hook` (all `optional = true`). Add
  `wasm-bindgen-futures` under `[target.'cfg(target_arch = "wasm32")'.dependencies]`
  (only pulled on wasm; not optional needed there, or gate via the `web` path —
  implementer's call as long as native builds don't pull it).
- If eframe's native default features (x11/wayland/glow) break the wasm32 check,
  split the eframe dep per-target and/or `default-features = false` + an explicit
  web feature set. Whatever makes the gate pass within this lane.
- Pick versions that match the workspace/`Cargo.lock` already-present transitives
  where possible (eframe 0.34 already pulls `wasm-bindgen`/`web-sys`/`js-sys`),
  to minimize new `Cargo.lock` churn.

### `src/lib.rs`
- Change `#[cfg(feature = "desktop")] pub mod app;` →
  `#[cfg(any(feature = "desktop", feature = "web"))] pub mod app;`
- Update the `app` module doc-comment so it no longer says only `desktop` pulls
  eframe (now `web` does too).

### `src/main.rs` — three-way cfg split
- Native desktop entry: gate the existing `run_native` `main` and its
  `default_spec`/imports on `all(feature = "desktop", not(target_arch = "wasm32"))`
  (so a `--features web` *native* check doesn't try to run_native on wasm, and a
  wasm build never compiles the native path).
- Web entry (`all(feature = "web", target_arch = "wasm32")`):
  ```rust
  #[cfg(all(feature = "web", target_arch = "wasm32"))]
  #[wasm_bindgen::prelude::wasm_bindgen(start)]
  pub fn start_web() {
      console_error_panic_hook::set_once();
      let web_options = eframe::WebOptions::default();
      wasm_bindgen_futures::spawn_local(async {
          eframe::WebRunner::new()
              .start("the_canvas_id", web_options,
                     Box::new(|_cc| Ok(Box::new(StudioApp::new(
                         StudioState::from_spec(default_spec()))))))
              .await
              .expect("yee-studio web start failed");
      });
  }
  ```
  Use a wasm-shared `default_spec()` (the same satisfiable Chebyshev N=5 BPF) —
  factor it so both the native and web paths see it (e.g. gate `default_spec` on
  `any(feature="desktop", feature="web")`). Confirm the `WebRunner::start`
  closure return type matches eframe 0.34 (`Result<Box<dyn App>, _>` vs
  `Box<dyn App>`); the native path here uses `Ok(Box::new(...))`, so match it.
- Stub `main`: gate on
  `not(any(all(feature="desktop", not(target_arch="wasm32")), all(feature="web", target_arch="wasm32")))`
  so exactly one `main`/entry compiles for every feature×target combination.
  Verify there is no "multiple/zero `main`" error in any of: default (desktop,
  native), `--no-default-features` (stub), `--features web --target wasm32`.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. **Primary gate:** `cargo check -p yee-studio --target wasm32-unknown-unknown
   --features web` exit 0 (the full eframe `StudioApp` compiles for the browser).
   If WebGPU needs `RUSTFLAGS='--cfg=web_sys_unstable_apis'` to pass, use it and
   record the exact RUSTFLAGS in the ADR + a `// ` note for App.1.2b.
3. **No regression — native default still builds:**
   `cargo check -p yee-studio` exit 0 (desktop path unchanged).
4. **No regression — headless wasm still builds (App.1.1 gate):**
   `cargo check -p yee-studio --no-default-features --target wasm32-unknown-unknown`
   exit 0, and egui/eframe absent from that dep tree
   (`cargo tree -p yee-studio --no-default-features --target wasm32-unknown-unknown -i egui`
   reports egui not in the tree / errors as "not found").
5. `cargo clippy -p yee-studio --features web --target wasm32-unknown-unknown -- -D warnings`
   exit 0 (best-effort; if clippy-on-wasm is flaky in this env, fall back to the
   native `cargo clippy -p yee-studio -- -D warnings` and note it).

## Out of scope
`trunk`/`index.html`/`Trunk.toml`; `cargo install trunk`; the static deploy; any
in-browser runtime test (all App.1.2b). Any `StudioState`/flow change. A CI job
for the web build (a follow-on; App.1.1's `wasm-build` job already proves the
headless path — extending it to `--features web` can be App.1.2b's CI step).
