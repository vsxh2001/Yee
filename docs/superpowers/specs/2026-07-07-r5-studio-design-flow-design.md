# R.5 — Studio spec→design→export flow

**Date:** 2026-07-07
**Track:** RF-TOOL-ROADMAP R.5
**Related:** ADR-0179 (Tauri studio scaffold), ADR-0180 (S.3 views), ADR-0197
(R.4 — the design machinery this surfaces), Phase 0 Touchstone gates,
`yee-export` Gerber gates.

## Problem

The Tauri studio runs raw FDTD jobs (`run_job`) and plots fields, but a
designer cannot *design* in it: no spec entry, no S-parameter response view,
no export. All the machinery exists in the libraries (synthesis, dimensions,
layouts, `coupling_matrix_s_params`, `yee_io::touchstone`,
`yee_export::layout_to_gerber*`) and is only reachable from the CLI.

## Scope — walking skeleton

One new Tauri command, `design_filter` (pure core `design_filter_impl` in
`studio/src-tauri/src/design.rs`), runs the **closed-form design flow**:
spec → `synthesize` → `dimension_hairpin_with_fold` → `hairpin_bpf_sections`
→ `coupling_matrix_s_params` over `f0·(1 ± 2·FBW)` → response curves +
ready-to-save artifacts (`.s2p` via `yee_io::touchstone::to_string`, copper +
outline Gerber via `yee_export`). Instant (no FDTD), so the studio stays
interactive; dimensioning errors (e.g. `TapNotRealizable`) surface as
designer-grade error strings.

React side: `FilterDesignPanel` (spec form: f0, FBW, order, ripple/
Butterworth, ε_r, h) + `SparamPlot` (dependency-free SVG, |S21| solid /
|S11| dashed, −60 dB floor) + export buttons (Blob downloads).

**Full-wave verify/refine in the studio is the follow-on** (R.5b): stream the
R.4 BO loop's per-solve progress over the existing `JobEvent` protocol. The
walking skeleton deliberately ships the instant design side first — the same
decomposition the engine track used (S.2 shell before S.3 views).

## Gates

- **`studio-design-e2e-001`** (`studio/src-tauri/tests/design_e2e.rs`,
  headless command layer, in the `studio-build` CI job): a scripted design
  flow whose artifacts are **byte-checked** — the `.s2p` round-trips through
  `yee_io::touchstone` (read → re-render → byte-identical) and the Gerbers
  are byte-identical to `yee_export`'s own output for the same synthesized
  layout; the response is a real band-pass design (passband ≈ 0 dB at f0,
  stopband < −15 dB, in-band S11 < −6 dB). Plus the unrealizable-spec path
  (narrow FBW on a thick board) returns the dims' qe-range error.
- **vitest DOM gates** (`studio/src/sparam.test.tsx`): `SparamPlot` renders
  both traces inside the viewBox (floor clamping) with the band caption;
  `FilterDesignPanel` renders the six-input spec form.

## Consequences

The studio covers spec → design → response → manufacturable artifacts
end-to-end with byte-checked fidelity to the validated libraries. R.5b
(full-wave loop streaming) and antenna specs are follow-ons on the same
panel pattern.
