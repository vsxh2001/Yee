# ADR-0110: Yee Filter Studio redesign — Dioxus view layer (POC-first)

**Status:** Accepted (POC-first)
**Date:** 2026-05-30
**Related:** ADR-0089 (desktop+web app), ADR-0090/0092/0096/0107 (eframe studio +
wasm), ADR-0097/0109 (dimensional synthesis), `FILTER-DESIGN-ROADMAP.md`,
[[project-filter-design-final-goal]]

---

## Context

The shipped eframe `yee-studio` reads as a chunky traditional EDA tool: no
topology selection, no real board/material view, weak styling. The maintainer
wants a **truly polished, pure-Rust, web-first** product for the full
filter-design journey. Competitive research (Ansys Nuhertz, Cadence/AWR iFilter,
Marki/Hokua web calculators, QucsStudio/openEMS) shows everyone either fragments
the flow across tools or is chunky/desktop/$$$; the open wedge is one slick
web-first app holding the whole journey with built-in EM verification.

## Decision

Rebuild the studio **view layer in Dioxus** (pure Rust, RSX). Web renders to real
DOM+CSS (SaaS-class polish egui can't reach); desktop via webview; one codebase.
`StudioState` (egui-free core) + the engine crates (`yee-synth`, `yee-filter`,
`yee-layout`, `yee-export`, `yee-plotters`) are **reused untouched**; board/plots
render as **SVG** in the DOM; a **design-system** (theme tokens + base components)
is the polish foundation. The eframe view is **retired** (when the full build
lands — not in the POC). Product flow = Shell A stage rail with six stages
(Spec → Technique → Synthesis → Layout+Materials → Verify → Export); decomposed
into App.D.0 shell + D.1–D.4 stages (parallel) + D.5 EM-verify (later). Full
design + flow + competitive analysis: the design spec.

**POC-first (this increment):** a `crates/yee-studio-web` Dioxus proof-of-concept
— design system + Shell A + `StudioState`/engine bridge, rendering **two real
stages** (Synthesis: coupling matrix + ideal |S21| SVG vs mask + PASS/FAIL;
Layout+Materials: real dimensioned board SVG + stackup + resonator table) from the
live engine, built for web and served locally for the maintainer to judge polish
before committing to the full build.

## Consequences

**Ships (POC):** a slick, real-data Dioxus web POC; validates polish + the
Rust→WASM engine bridge + the web build. `StudioState`/engine/eframe `yee-studio`
all untouched (the POC is additive).

**Then:** on POC approval, the full App.D.1–D.4 stages fan out in parallel atop
the shell; the eframe view is retired; App.D.5 (EM-verify) follows the FDTD-in-loop
work.

**Not in scope (POC):** desktop/webview, full per-stage interactivity, EM-verify,
eframe deletion, CI wiring, mobile.

---

## References
- `docs/superpowers/specs/2026-05-30-studio-redesign-dioxus-design.md`.
- Research: Ansys Nuhertz FilterSolutions; Cadence AWR iFilter; Marki web Microstrip
  Filter Tool; Tidy3D coupled-line FDTD; the 2025 Rust GUI survey (Dioxus/Slint/egui).
