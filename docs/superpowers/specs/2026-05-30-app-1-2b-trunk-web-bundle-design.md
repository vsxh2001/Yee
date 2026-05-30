# App.1.2b — `trunk` web bundle + CI build for `yee-studio` — Design Spec

**ADR:** ADR-0107 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Make the `yee-studio` browser app **buildable and runnable in a browser**, not
just compilable. App.1.2a (ADR-0096) proved the full eframe `StudioApp` compiles
to `wasm32-unknown-unknown` behind the `web` feature with a
`#[wasm_bindgen(start)] start_web()` entry that mounts on
`<canvas id="the_canvas_id">`. What is missing is the bundling layer: an
`index.html`, a `Trunk.toml`, and a CI job that runs `trunk build` to emit the
deployable static bundle (HTML + JS glue + `.wasm`). After this, `trunk serve`
gives a working in-browser filter-design UI, and CI proves the bundle builds.

This is the "working web UI" component of the product goal — reached without
running the heavy wasm bundle on the local (memory-constrained) box: the full
`trunk build` runs on a **GitHub CI runner**.

## Changes
Two lanes (one focused agent, both scoped):

### `crates/yee-studio/**`
- `index.html` — minimal trunk entry:
  - a `<canvas id="the_canvas_id"></canvas>` (the id `start_web` looks up);
  - a trunk rust directive building THIS crate's `yee-studio` bin for wasm with
    the `web` feature and `--no-default-features` (so the native `desktop`
    feature/ winit stack does not leak in):
    `<link data-trunk rel="rust" data-bin="yee-studio" data-cargo-no-default-features data-cargo-features="web" data-wasm-opt="2" />`
  - a little CSS so the canvas fills the viewport; a `<title>Yee Filter Studio</title>`.
- `Trunk.toml` — `[build] target = "index.html"` + `dist = "dist"`; pin nothing
  exotic. Keep it minimal; feature selection lives in the `data-trunk` link.

### `.github/workflows/ci.yml`
- Extend the existing `wasm-build` job (or add a sibling `web-bundle` job) to,
  AFTER the existing `--no-default-features` compile gate:
  - install trunk (`cargo binstall trunk` if available, else
    `cargo install --locked trunk`; or the `jetli/trunk-action`);
  - install the `wasm-bindgen-cli` matching the locked `wasm-bindgen` if trunk
    needs it (trunk usually fetches it itself — prefer letting trunk manage);
  - run `trunk build crates/yee-studio/index.html` (release);
  - upload `crates/yee-studio/dist` as a build artifact (`actions/upload-artifact`).
- Do NOT add a GitHub Pages deploy step — `docs.yml` already owns the Pages
  deployment for the mdBook; a second Pages deploy would clash. Pages/static
  hosting of the app is a follow-on (App.1.2c).

## DoD (machine-checkable)
1. `crates/yee-studio/index.html` exists, contains `id="the_canvas_id"`, a
   `data-trunk rel="rust"` link selecting `data-cargo-features="web"` +
   `data-cargo-no-default-features`, and is well-formed HTML.
2. `crates/yee-studio/Trunk.toml` exists and parses as TOML (`[build]` table with
   a `target`).
3. `.github/workflows/ci.yml` is valid YAML (parses) and has a step/job that runs
   `trunk build` on `crates/yee-studio/index.html` and uploads the `dist`
   artifact; the existing `--no-default-features` wasm gate is unchanged.
4. The native build is unaffected: `cargo build -p yee-studio` (default
   `desktop`) still compiles. (Light; the agent MAY skip running it and instead
   confirm no Cargo.toml change — this PR adds no deps.)

**Local-verification note (memory-constrained box):** do NOT run `trunk build` or
`cargo build --features web --target wasm32` locally — both are heavy and the box
OOM-kills on large builds. App.1.2a already proved the wasm `--features web`
compile is green; the `trunk build` proof is delegated to CI. Local checks = file
existence + HTML/TOML/YAML well-formedness + a grep for the required tokens.

## Out of scope
GitHub Pages / public hosting of the app (App.1.2c — must not clash with the
mdBook Pages deploy); a JS framework; server calls (App.2 `yee-server`); any UI
change (the StudioApp already renders synth + dims + layout panels).
