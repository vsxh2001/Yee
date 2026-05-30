# Yee Filter Studio — product redesign (Dioxus) — Design Spec

**Date:** 2026-05-30 · **Status:** Accepted (POC-first) · **ADR:** ADR-0110

## Problem

The shipped `yee-studio` (egui/eframe) works but reads as a "chunky traditional tool":
no filter-topology selection/visualization, no real board/component view, weak
styling. The maintainer wants a **truly polished, pure-Rust, web-first** product
covering the full filter-design journey.

## Competitive landscape (research)

- **Integrated commercial** (Ansys Nuhertz FilterSolutions→HFSS; Cadence AWR
  Microwave Office + iFilter→AXIEM): full spec→synth→layout→EM→optimize→parts-DB,
  but **$$$$, desktop, chunky multi-window, fragmented** across tools.
- **Free web calculators** (Marki, Hokua, Wcalc, Pasternack): web/mobile/free but
  **calculators, not designers** — one quantity at a time, **no built-in EM
  verify**, no manufacturing output; all say *"first approximation, verify in 3D
  EM elsewhere."*
- **Open desktop** (QucsStudio+openEMS, Qucs-RFlayout, scikit-rf): the pieces
  exist but are **script-glued, desktop, fragmented**.

**Wedge:** one **slick, web-first** app holding the *entire* journey — spec →
technique → synthesis → dimensioned layout *with real materials* → **built-in
FDTD-verified response** → manufacturable files — and **honest about
ideal-vs-realized** (the thing calculators hide). Exactly Yee's stack.

## Product flow (approved)

Shell **A** — slim **stage rail** (left) + dominant central canvas. Six stages:

1. **Spec** — f0, bandwidth, order/auto, ripple, return loss, stopband mask, Z0;
   live realizability check.
2. **Technique** — topology gallery (edge-coupled + hairpin available; combline /
   interdigital / lumped / stepped-impedance greyed "Soon") + medium
   (microstrip/stripline) + substrate library (FR-4, Rogers, alumina…).
3. **Synthesis** — g-values, Qe, coupling matrix M, ideal |S21| vs mask, PASS/FAIL.
4. **Layout + Materials** — top-view board (copper traces, ground, ports,
   dimension callouts, layer toggles) + **material stackup** cross-section
   (εr / h / Cu-thickness / loss, substrate library) + **components table**
   (per-resonator W/L/gap → Z0e/Z0o/εeff/realized-k). Editable; recomputes live.
5. **Verify (EM)** — FDTD-simulated *realized* response overlaid on ideal + mask
   (closes the ideal-vs-realized gap); tune / auto-optimize. *(rides on
   F1.1b.2/F1.2.1; later.)*
6. **Export** — Gerber · KiCad `.kicad_pcb` · Touchstone · STEP + a **final
   parameter sheet** (design summary, per-section geometry, verified IL/RL/
   rejection, downloads).

## Technical decision: view layer → **Dioxus**

Pure Rust, RSX. **Web target renders to real DOM + CSS** → SaaS-class polish
(animations/transitions/typography) — the ceiling egui can't reach. **Desktop**
= the same UI via webview. One codebase.

- **Reused untouched:** `StudioState` (egui-free core) + engine crates
  (`yee-synth`, `yee-filter`, `yee-layout`, `yee-export`, `yee-plotters`). The
  view calls into them; physics/synthesis/export unchanged.
- **Retired:** the eframe view (`crates/yee-studio` `app.rs` + the eframe `main`).
  `StudioState` and the proven light-flow survive; only egui rendering goes.
- **Board / plots:** rendered as **SVG** in the DOM (crisp, scalable,
  CSS-styleable) — not a raster canvas.
- **Design system:** a theme-tokens module (palette, type scale, spacing, base
  components) is the "slick" foundation, applied app-wide. Dark, flat, teal
  accent (`#2dd4bf`), copper `#e6b24d`, refined type — per the approved mockups.

## Decomposition (each its own spec→plan→build; fan out after the shell)

- **App.D.0 — Dioxus shell + design system + StudioState bridge** *(this POC).*
- **App.D.1** Spec + Technique (topology gallery) · **D.2** Synthesis (matrix +
  response/mask SVG) · **D.3** Layout+Materials+Components · **D.4** Export
  (param sheet + downloads) — largely disjoint → **parallel fan-out**.
- **App.D.5** EM-Verify — later (FDTD-in-loop dependency).

## THIS deliverable — the POC (App.D.0)

A **proof-of-concept** to validate (a) Dioxus delivers the polish, (b) the
`StudioState`/engine bridge works in Rust→WASM, (c) it builds + runs on web.
Maintainer will judge the POC before committing to the full build ("do a poc,
then will see").

**Scope:**
- New crate `crates/yee-studio-web` (Dioxus app; does NOT disturb `yee-studio`
  yet — eframe retirement happens when the full build replaces it, not in the POC).
- The **design-system CSS** (tokens + base components) from the approved aesthetic.
- **Shell A**: top bar + stage rail (all six stages shown; rail navigable).
- **Real engine data** (not lorem): the POC drives `yee_synth`/`yee_filter` on the
  committed Chebyshev N=5 fixture and renders **two real stages slickly** —
  **Synthesis** (coupling-matrix grid + ideal |S21|/|S11| **SVG** plot vs mask +
  PASS/FAIL) and **Layout+Materials** (board top-view SVG from
  `dimension_edge_coupled_layout` + stackup + resonator table). Spec/Technique/
  Export rendered as styled-but-static stage stubs (prove the shell, not full
  interactivity).
- **Web build**: compiles to `wasm32-unknown-unknown` and serves locally for the
  maintainer to open in a browser and judge polish. (Desktop/webview deferred.)

**DoD (POC):**
1. `crates/yee-studio-web` builds for web (wasm32); `cargo fmt`/`clippy -D warnings`
   clean on the crate.
2. Drives the real engine: the Synthesis stage shows the actual computed coupling
   matrix + a real |S21| SVG vs the mask + the real PASS/FAIL; the Layout stage
   shows the real dimensioned board SVG + the real per-resonator table.
3. Serves locally (a URL the maintainer opens); looks slick (the approved dark/
   teal/copper design system, SVG graphics, refined type/spacing).
4. `StudioState` + engine crates **unchanged**; `yee-studio` (eframe) untouched
   (parallel, retired only when the full build lands).

## Out of scope (POC)

Desktop/webview target; full interactivity on every stage; EM-verify stage;
topology beyond a static gallery; the eframe deletion; CI wiring; mobile.

## Engine honesty (the maintainer's material/strips point)

The Synthesis response is the **ideal closed-form prototype**; the Layout
dimensions DO use material (εr, h) + strips (W, S) via Hammerstad-Jensen +
Kirschning-Jansen. Metal thickness, loss, dispersion, and the **full-wave
realized response** are NOT yet in the response — the POC labels the response
"ideal (prototype)" honestly; closing the loop is stage 5 / F1.2.1 (later).
