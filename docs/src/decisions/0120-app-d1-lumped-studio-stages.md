# ADR-0120: App.D.1L — lumped-LC studio stages on the Dioxus shell

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0110 (Dioxus redesign, POC `b8cfb90`), ADR-0111/0112/0113/0114
(F2.0 lumped synth / F2.1 BOM / F2.4 tolerance / F2.2 lumped board — the engine
this surfaces), the lumped-LC → PCB goal (maintainer chose "approve direction,
build out stages"), [[project-lumped-lc-and-studio-redesign]]

---

## Context

The maintainer approved building out the Dioxus studio stages. The active goal is
the **lumped-LC** journey, whose engine is shipped (F2.0 ladder synth, F2.1
component-choosing + BOM, F2.4 tolerance/yield, F2.2 lumped board) but **not
surfaced in any UI**. The POC (App.D.0) renders only the *distributed* Synthesis +
Layout stages. The goal explicitly names "polished UI, component choosing, BOM,
tolerance" — so the highest-value build-out is the **lumped-LC stage set**.

## Decision

On the Dioxus POC branch (`feature/app-d0-dioxus-poc`), add a **Lumped-LC flow**
driven by the real shipped engine, rendered in the existing dark/teal/copper
design system with SVG graphics:

1. **Technique**: make **Lumped LC** a selectable, live topology (alongside the
   greyed distributed entries) — selecting it routes the downstream stages to the
   lumped engine.
2. **Synthesis (lumped)**: `yee_filter::synthesize_lumped` → the LC ladder
   (series/shunt resonators, L/C values) as a styled table + the ideal `ladder_s21`
   |S21| **SVG** vs the spec mask + PASS/FAIL.
3. **Components + BOM**: `yee_filter::select_components` (E24/E96 toggle) → the
   **BOM table** (per-element L/C, chosen E-series value, deviation %, qty) — the
   goal's "component choosing" + "BOM" made visible.
4. **Tolerance / yield**: `yee_filter::monte_carlo_yield` → the yield % + worst-case
   RL/rejection (E24 vs E96), with the honest narrowband-yield insight surfaced —
   the goal's "tolerance consideration".
5. **Layout (lumped)**: `yee_filter::lumped_board` → the dimensioned lumped board
   **SVG** (footprints, pads, traces) + the placement/footprint table.

All physics/synthesis/export crates are **reused untouched** (the view calls in).
Stays on the user-gated POC branch; the maintainer reviews the fuller build before
any merge / eframe retirement (ADR-0110's gate).

## Consequences

**Ships (to the POC branch, for maintainer review):** the goal's "polished UI +
component-choosing + BOM + tolerance" components, rendered slickly from real engine
data — the lumped-LC journey end-to-end in the browser (EM-Verify still pending
Track A). Closes 2 of the goal's open components into a reviewable artifact.

**Gate:** the crate builds for wasm32 + `cargo clippy -D warnings`/`fmt` clean; the
lumped stages render the *real* computed ladder / BOM / yield / board (not
placeholder); served locally for the maintainer. Engine crates unchanged.

**Not in scope:** merging to `main` / retiring eframe (maintainer-gated, ADR-0110);
the EM-Verify stage (Track A dependency); desktop/webview; the distributed-flow
stage build-out (App.D.1/D.2/D.3 proper — follow-on).

---

## References
- `docs/superpowers/specs/2026-05-30-app-d1-lumped-studio-stages-design.md`;
  `docs/superpowers/plans/2026-05-30-app-d1-lumped-studio-stages.md`.
- ADR-0110 + `2026-05-30-studio-redesign-dioxus-design.md` (the shell + flow).
