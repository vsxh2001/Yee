# ADR-0198: R.5 — the studio designs: spec → response → byte-checked exports

**Status:** Accepted
**Date:** 2026-07-07
**Related:** ADR-0179 (Tauri scaffold), ADR-0180 (S.3 views), ADR-0197 (R.4
design machinery), RF-TOOL-ROADMAP R.5.
**Spec:** `docs/superpowers/specs/2026-07-07-r5-studio-design-flow-design.md`

## Decision

The Tauri studio gains the **filter design flow** as a second command beside
`run_job`:

- `design_filter` (`studio/src-tauri/src/design.rs`, pure core
  `design_filter_impl`): spec → `synthesize` → `dimension_hairpin_with_fold`
  (the R.4a fold-corrected, qe→tapped, per-section dims) →
  `hairpin_bpf_sections` → `coupling_matrix_s_params` over `f0·(1 ± 2·FBW)` →
  response curves, synthesized dimensions, and ready-to-save artifacts: a
  `.s2p` (via `yee_io::touchstone::to_string`; the coupling-matrix model is
  passive by construction) and copper + outline Gerbers (via
  `yee_export::layout_to_gerber`/`_outline`). Closed form and instant — the
  studio stays interactive; dims errors (`TapNotRealizable` with its
  realizable qe range) surface verbatim as designer-grade messages.
- React: `FilterDesignPanel` (f0/FBW/order/ripple/ε_r/h spec form, dims
  readout, export buttons via Blob downloads) + `SparamPlot`
  (dependency-free SVG, |S21| solid / |S11| dashed, −60 dB floor) — the
  studio's first response-vs-frequency view (everything before R.5 plotted
  fields).

**Deliberately out of the walking skeleton:** running the full-wave verify /
R.4 BO loop from the studio (R.5b — streams per-solve progress over the
existing `JobEvent` protocol) and antenna specs. Same decomposition as the
engine track: shell first, loops next.

## Gates

- **`studio-design-e2e-001`** (headless command-layer e2e, `studio-build` CI
  job): a scripted design flow with **byte-checked exports** — the `.s2p`
  round-trips `yee_io::touchstone` read → re-render **byte-identical**, and
  both Gerbers are **byte-identical** to `yee_export`'s output for the same
  synthesized layout (the studio adds no bytes of its own); the response is a
  real band-pass design (S21(f0) > −1 dB, stopband < −15 dB, in-band
  S11 < −6 dB); the unrealizable path (FBW = 0.10 on the 1.6 mm board — tap
  beyond the fold-shortened arm) returns the qe-range error instead of
  panicking.
- **vitest DOM gates** (`sparam.test.tsx`): `SparamPlot` keeps every vertex
  inside the viewBox under floor clamping and captions the band;
  `FilterDesignPanel` renders the six-input form. 14 vitest tests total.

## Consequences

A designer can enter a spec, see the designed response, and leave with a
valid `.s2p` and manufacturable Gerbers — without touching the CLI. The
export fidelity is pinned byte-for-byte to the validated libraries, so studio
artifacts can never silently diverge from `yee export`'s. R.5b (full-wave
loop streaming) queued.
