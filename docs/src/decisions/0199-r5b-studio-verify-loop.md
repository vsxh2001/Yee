# ADR-0199: R.5b — the studio verifies: full-wave loop + the shared board fixture

**Status:** Accepted
**Date:** 2026-07-07
**Related:** ADR-0198 (R.5 design flow), ADR-0197 (R.4 — the measurement
fixture this promotes to a library), ADR-0186/0187/0189 (the S.9/S.10/S.12
fixture itself).
**Spec:** `docs/superpowers/specs/2026-07-07-r5-studio-design-flow-design.md`
(R.5b section)

## Decision

Two things land together, deliberately:

1. **`yee_engine::board`** — the voxelize → `JobSpec` measurement fixture
   that every board gate (S.8–S.12, R.0–R.4) carried as copy-pasted test
   code, promoted to a library API: `two_port_board_job(layout, opts)` builds
   the S.9/S.10-certified fixture (CPML-xy + PEC ground/lid, aperture ports
   on the feed cross-sections, two ordered 3-probe triples for the S.12
   directional observables) from any two-port `yee_layout::Layout`;
   `reference_through_line` builds the matching straight-line reference on
   the shared bbox/grid. `TwoPortBoardOptions::for_band` carries the R.4
   gate's defaults. Unit-gated (fixture shape, shared-grid property,
   short-feed rejection). yee-engine gains yee-voxel/yee-layout deps (no
   cycle; the crate is native-only).
2. **`verify_filter`** (studio, `verify.rs`) — the full-wave verify loop:
   rebuild the designed hairpin with verify-length feeds, run the
   straight-line reference + DUT through the engine (phases streamed to the
   webview as `verify://progress {phase, step, total}`), and return the
   measured directional |S21| next to the design's coupling-matrix curve.
   React: a "Verify (full-wave)" button on the design panel, phase-tagged
   progress, and a measured-vs-designed `SparamPlot` overlay (the plot
   gained a `labels` prop).

Because the studio measures through the **same builder the gates use**, the
studio's fixture cannot silently drift from the certified one — the R.5
byte-check philosophy applied to measurement instead of export.

## Gate `studio-verify-e2e-001`

Headless command-layer e2e at reduced fidelity (coarse dx = 0.4 mm, 700
steps, ~1–2 min debug — a **pipe** gate; the fixture's physics is gated in
yee-engine/yee-filter):

1. progress streams for both phases, reference strictly before DUT, monotone
   within each phase;
2. **the identity case measures exactly**: verifying a straight line whose
   reference is the same line runs two bit-identical solves, so the
   directional |S21| must read **0 dB to numerical identity** (< 1e−9 dB) —
   a strong end-to-end determinism + post-processing check;
3. short feeds surface the shared builder's error, not a panic.

Runs in the `studio-build` CI job beside the R.5 design e2e.

## Consequences

The studio now covers the full R.5 ambition: spec → instant design + exports
→ full-wave verify with live progress — every stage riding validated,
shared library code. The gates' fixture duplication can now be retired
opportunistically (migrate test-by-test onto `yee_engine::board` when each
file is next touched). Queued: R.4c (fine-grid BO), R.2b, R.0b; a studio BO
loop button becomes trivial once R.4c makes it worth minutes of user time.
