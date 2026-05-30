# ADR-0107: App.1.2b — `trunk` web bundle + CI build for `yee-studio`

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0096 (App.1.2a — wasm `web` feature + `WebRunner` entry),
ADR-0095 (App.1.1 — CI `wasm-build` gate), ADR-0089 (desktop+web app),
`FILTER-DESIGN-ROADMAP.md`

---

## Context

App.1.2a (ADR-0096) made the full eframe `StudioApp` compile to
`wasm32-unknown-unknown` behind the `web` feature, with a
`#[wasm_bindgen(start)] start_web()` entry that mounts on
`<canvas id="the_canvas_id">`. But "compiles to wasm" is not "runs in a browser":
there is no `index.html`, no bundler config, and nothing emits the deployable
HTML+JS+`.wasm` triple. The product goal names a **working web UI**; that needs
the bundling layer.

The local development box is memory-constrained (large Rust builds OOM-kill it),
and the standing "use less CPU" constraint rules out running the heavy full-eframe
wasm bundle locally. So the heavy `trunk build` must run on a **CI runner**, not
this box.

## Decision

Add the trunk bundling layer for `yee-studio`:

- `crates/yee-studio/index.html` — a trunk entry with `<canvas
  id="the_canvas_id">` and a `data-trunk rel="rust"` directive building the
  `yee-studio` bin for wasm with `--no-default-features --features web` (so the
  native winit/x11 stack does not leak into the wasm build).
- `crates/yee-studio/Trunk.toml` — minimal `[build]` config.
- `.github/workflows/ci.yml` — extend the existing `wasm-build` job to install
  trunk, run `trunk build --release crates/yee-studio/index.html`, and upload the
  `dist/` bundle as a CI artifact. The pre-existing `--no-default-features` wasm
  compile gate is kept unchanged.

No GitHub Pages deploy: `docs.yml` already owns the Pages deployment (mdBook), and
a second would clash. Public hosting is a follow-on (App.1.2c).

## Consequences

**Ships:** `trunk serve` (locally, on any non-constrained machine) gives a
working in-browser Yee filter-design UI (synth + dims + layout panels), and CI
produces a downloadable static `dist/` artifact on every push — proving the bundle
builds. This is the "working web UI" goal component, reached without running the
heavy bundle on the constrained local box.

**Gate:** structural/CI — `index.html` (canvas id + web-feature directive),
`Trunk.toml` (valid `[build]`), `ci.yml` (valid YAML; a `trunk build` + artifact
step). The real build proof is the green CI job; the local checks are
well-formedness only (per the memory constraint).

**Not in scope:** Pages/public hosting (App.1.2c); `yee-server` calls (App.2); any
UI change.

---

## References
- ADR-0096 (the wasm `web` entry this bundles); ADR-0095 (the CI wasm gate this
  extends).
- `docs/superpowers/specs/2026-05-30-app-1-2b-trunk-web-bundle-design.md`;
  `docs/superpowers/plans/2026-05-30-app-1-2b-trunk-web-bundle.md`.
- Trunk (<https://trunkrs.dev>) — `data-trunk rel="rust"` directive,
  `data-cargo-features` / `data-cargo-no-default-features`.
