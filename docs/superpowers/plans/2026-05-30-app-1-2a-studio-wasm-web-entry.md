# App.1.2a — yee-studio wasm32 web entry — Plan

**Spec:** `2026-05-30-app-1-2a-studio-wasm-web-entry-design.md` · **ADR:** ADR-0096

## Lane
`crates/yee-studio/**` ONLY. Out of lane → finding, do NOT fix. Do NOT touch
`yee-synth`/`yee-filter`/`yee-layout`/`yee-plotters` (they stay WASM-safe; this is
their downstream consumer). No `StudioState`/flow-logic change — wiring only.

## Base
New worktree off current `main` (base SHA pinned in the brief). Branch
`feature/app-1-2a-studio-web-entry`.

## Pattern files
- `crates/yee-studio/Cargo.toml` — the existing `desktop` feature + optional
  eframe/egui/egui_plot; mirror its shape for `web`.
- `crates/yee-studio/src/main.rs` — the existing two-arm
  `#[cfg(feature="desktop")]` / `#[cfg(not(feature="desktop"))]` `main`; extend to
  the three-way native/web/stub split.
- `crates/yee-studio/src/lib.rs` — `#[cfg(feature="desktop")] pub mod app;` (widen
  to `any(desktop, web)`).
- eframe 0.34 web entry: eframe_template `src/main.rs` (`WebRunner::new().start(...)`
  + `#[wasm_bindgen(start)]`). Confirm the exact 0.34 `start` signature / closure
  return type against the installed eframe (`cargo doc -p eframe --open` is
  unavailable headless — read the source under `~/.cargo` or docs.rs for 0.34).

## Steps
1. `Cargo.toml`: add the `web` feature + optional `wasm-bindgen`/`web-sys`/
   `console_error_panic_hook` deps + `wasm-bindgen-futures` as a `cfg(wasm32)`
   target dep. Prefer versions already in `Cargo.lock` (eframe 0.34 pulls
   wasm-bindgen 0.2 / web-sys 0.3 transitively) to minimize churn.
2. `lib.rs`: widen the `app` module gate to `any(feature="desktop", feature="web")`;
   fix its doc-comment.
3. `main.rs`: factor `default_spec()` to `any(desktop, web)`; gate the native
   `run_native` `main` on `all(desktop, not(wasm32))`; add the
   `#[wasm_bindgen(start)]` `WebRunner` web entry on `all(web, wasm32)`; re-gate the
   stub `main` so exactly one entry compiles per feature×target.
4. Run the DoD gates (below). Iterate the eframe feature wiring until gate 2 passes
   (try per-target eframe dep / `default-features=false` / the `web_sys_unstable_apis`
   RUSTFLAGS as needed — all within this lane). Record any RUSTFLAGS in the ADR.

## Verify (exit 0; nice -n 19, --jobs 2 — these are wasm `check`s, one-time compile)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo check -p yee-studio --target wasm32-unknown-unknown --features web --jobs 2          # PRIMARY gate
nice -n 19 cargo check -p yee-studio --jobs 2                                                          # native unchanged
nice -n 19 cargo check -p yee-studio --no-default-features --target wasm32-unknown-unknown --jobs 2    # App.1.1 gate holds
nice -n 19 cargo tree -p yee-studio --no-default-features --target wasm32-unknown-unknown -i egui      # egui must be ABSENT
nice -n 19 cargo clippy -p yee-studio --features web --target wasm32-unknown-unknown --jobs 2 -- -D warnings
```
Do NOT run `cargo test --workspace`, FDTD, mom-001, or any release build. The
wasm dep-tree compile (eframe+wgpu for wasm32) is a one-time cost; niced + --jobs 2.

## Escape hatch
Blocked > 15 min — the `--features web --target wasm32` check will not pass
without `trunk`/extra tooling, OR wgpu-on-wasm fundamentally fails to compile
(not just a feature-flag/RUSTFLAGS fix), OR the eframe 0.34 `WebRunner::start`
signature differs materially from the spec snippet → STOP, surface the exact
compiler error + what feature/RUSTFLAGS combinations you tried + the eframe 0.34
`start` signature you found. Do NOT install `trunk` (that is App.1.2b). Do NOT
weaken the App.1.1 headless gate (gate 4: egui must stay absent from the
`--no-default-features` wasm tree).

## Done when
DoD 1–5 pass; `git diff --stat <base>..HEAD` shows only `crates/yee-studio/**`
(+ the 3 committed docs already on the base if rebased in); the native desktop
build and the `--no-default-features` headless wasm build are both unregressed.
