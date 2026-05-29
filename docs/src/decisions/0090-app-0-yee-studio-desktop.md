# ADR-0090: App.0 — `yee-studio` filter-design desktop app skeleton

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0089 (filter-design app architecture), `FILTER-DESIGN-ROADMAP.md`
§5a (App track), ADR-0084/0086/0087 (the shipped light-flow crates)

---

## Context

ADR-0089 fixed the final deliverable as a desktop + web app (one `egui`/`eframe`
codebase, native + WASM; light flow client-side, heavy EM on a native server).
**App.0** is the first product increment: a native `eframe` desktop app that
wires the *already-shipped light flow* (F0/F0.1/F0.2/F1.0) into stage-gated
panels. It needs no EM and no server — it consumes only `yee-synth`/`yee-filter`
(synthesis + ideal response + spec-mask) — so it can ship now, in parallel with
the F1.1+ engine work.

## Decision

New crate **`yee-studio`** (lib + bin), an `eframe` app seeded from `yee-gui`'s
shape (lib.rs/main.rs/app.rs, `eframe` `["wgpu"]` + `egui` + `egui_plot`):

- **`StudioState`** (lib, pure + testable, NO egui): holds an editable
  `FilterSpec`; `recompute()` derives the `FilterProject` (`yee_filter::
  synthesize`), the ideal-response sweep, the |S21| dB trace, the spec-mask
  regions, and the `check_mask` verdict. All headless-testable.
- **`StudioApp`** (`impl eframe::App`): a left **spec-editor** panel (f0, FBW,
  order, ripple, return-loss, stopband points, approximation), a central
  **synthesis** panel (g-values, coupling matrix, Qe, mask PASS/FAIL), and an
  **`egui_plot`** |S21|-vs-mask view (live, recomputed on edit).
- **`main.rs`**: thin `eframe::run_native(StudioApp::default())`.

The **layout** preview panel is deferred to a later App increment: turning a
spec into the *right* physical dimensions needs the F1.2 coupling-matrix→dims
mapping, which is not built yet. App.0 covers spec → synthesis → spec-mask plot.

## Consequences

**Ships:** `yee-studio` crate (workspace member). Gate (per §4, crate test):
`StudioState::recompute()` for a satisfiable Chebyshev BPF yields the published
F0 g-values + a mask **PASS**, and a too-low order yields **FAIL** — headless,
no rendering. The `eframe` bin is build-only (a windowed app is not run in CI,
matching `yee-gui`). `cargo build -p yee-studio` + `clippy -D warnings` green.

**Constraint (ADR-0089):** `StudioState` (the logic) stays WASM-safe (no
native-only deps) so App.1 can compile it to WASM; only the `eframe`/windowing
shell is native-specific.

**Not in scope:** the layout preview (needs F1.2); EM/server calls (App.2);
WASM build (App.1); project save/load + export (App.3). No new external
dependency beyond the egui stack already in the workspace.

---

## References
- ADR-0089 (app architecture), `FILTER-DESIGN-ROADMAP.md` §5a;
  `crates/yee-gui/` (eframe app pattern); the shipped `yee-synth`/`yee-filter`.
- `docs/superpowers/specs/2026-05-29-app-0-yee-studio-desktop-design.md`;
  `docs/superpowers/plans/2026-05-29-app-0-yee-studio-desktop.md`.
