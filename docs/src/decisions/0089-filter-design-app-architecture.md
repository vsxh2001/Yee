# ADR-0089: Filter-design delivered as a desktop + web app (egui/eframe native + WASM)

**Status:** Accepted
**Date:** 2026-05-29
**Related:** `FILTER-DESIGN-ROADMAP.md` (§0 Vision, §5a App/Studio track),
ADR-0011 (egui 0.34 / wgpu 29), ADR-0018 (yee-gui validation panel)

---

## Context

The stated **final goal** is an *app / web app* for end-to-end RF filter design
— not just a CLI/library flow. The F-series (F0–F4) builds the design *engine*
(synthesis → coupling matrix → layout → EM verification → export); this ADR
fixes how that engine is delivered as a usable application on desktop **and** in
the browser, without a second UI codebase.

Constraints from the existing stack:
- The GUI is `egui` 0.34 + `eframe` + `wgpu` 29 (ADR-0011) — `eframe` already
  targets **native** and **WASM** (browser, WebGL/WebGPU) from one codebase.
- The light flow (spec → `yee-synth` → `yee-filter` → `yee-layout` →
  `yee-plotters` spec-mask) is pure Rust with no native-only deps → **WASM-safe**
  and already shipped (F0/F0.1/F0.2/F1.0).
- The heavy flow (FDTD/FEM verification, `yee-surrogate` dimensional synthesis,
  Gmsh meshing, CAD export) is native-compute-bound and not practical to run in
  the browser.

## Decision

Ship the filter designer as **one `egui`/`eframe` application, built for two
targets — native (desktop) and WASM (web)** — rather than a separate JS/TS
frontend. Split compute by weight:

- **Client (egui app, native or WASM):** the light flow — spec entry, synthesis,
  coupling-matrix display, `yee-layout` geometry preview, `draw_sparam_with_mask`
  ideal-response view. Runs fully in-process / in-browser; no server needed for
  the design-front-end stages.
- **`yee-server` (new crate, axum, native):** the heavy steps — FDTD/FEM
  verification, surrogate dimensional synthesis, mesh, and KiCad/Gerber/STEP
  export — behind JSON/artifact endpoints. The **web** client calls it over
  HTTP; the **desktop** app links the engine directly (same trait-boundary, two
  transports).
- **`yee-studio` (new crate):** the `eframe` app itself (stage-gated panels),
  seeded from the existing `yee-gui` panels.

### Why egui/eframe WASM rather than a JS frontend

The team's UI investment, plot code (`yee-plotters`/egui Smith/S-param/VSWR/
spec-mask), and domain types are all Rust. `eframe` compiles the same app to
WASM, so a web app costs a build target + a thin server, not a rewrite. A JS SPA
would duplicate the plotting + domain model and add a language boundary.

### Why split light/heavy rather than all-WASM or all-server

All-WASM can't run the heavy native EM at usable speed; all-server makes the
(already-shipped, instant) synthesis/design stages needlessly round-trip. The
split lets the design front-end work offline/instantly in the browser and only
reaches the server for the genuinely heavy F1.1+ steps.

## Consequences

**Adds (App/Studio track, §5a):** `yee-studio` (eframe app, native first then
WASM) and `yee-server` (axum). The light flow is reused as-is (no changes to the
shipped F0–F1.0 crates needed beyond keeping them WASM-safe — no native-only deps
in `yee-synth`/`yee-filter`/`yee-layout`/`yee-plotters`).

**Constraint going forward:** keep the light-flow crates WASM-compatible (no
`std::fs`/threads-only/native-only deps on their default path); push anything
native-only behind the `yee-server` boundary.

**Deferred:** WebGPU vs WebGL backend choice for the WASM viewport (egui handles
both); auth/multi-user for a hosted instance; packaging/signing the desktop app.

---

## References
- `FILTER-DESIGN-ROADMAP.md` §0 / §5a; ADR-0011 (egui/wgpu), ADR-0018 (yee-gui).
- eframe native+web targets; the shipped F0/F0.1/F0.2/F1.0 light-flow crates.
