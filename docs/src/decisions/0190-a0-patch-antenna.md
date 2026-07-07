# ADR-0190: A.0 — the antenna track opens; patch resonance verified at 0.0 %

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0182..0189 (the measurement stack this inherits), the project goal
"design an antenna or filter with the engine".
**Spec:** `docs/superpowers/specs/2026-07-06-a0-patch-antenna-design.md`

## Decision

The antenna track (A.*) opens as Part 3 of ENGINE-STUDIO-ROADMAP, walking-skeleton
first, mirroring how the filter verify chain was built:

1. **Synthesis** (`yee-layout`, pure closed forms, WASM-safe): `patch_antenna_dims`
   (Balanis §14.2 transmission-line model — `W = c/(2f₀)·√(2/(ε_r+1))`, Hammerstad
   `ε_eff`, open-end ΔL, `L = c/(2f₀√ε_eff) − 2ΔL`) and `edge_fed_patch` (feed line +
   patch `Layout`, one port). `open_end_delta_l` is now a public `yee-layout` function
   (the stub gate previously inlined it). Unit-gated against hand-computed values
   (2.45 GHz FR-4: W = 37.26 mm, ε_eff = 4.09, L ≈ 28.8 mm).
2. **Verify** (`engine-antenna-001`, `yee-engine/tests`, release CI via the existing
   `--include-ignored` step): voxelize → S.5 materials job → CPML-xy walls → S.10
   aperture drive → S.7 two-run |S11|.

## Measured

**The |S11| dip lands at 2.450 GHz — 0.0 % from the Balanis design frequency**
(25 MHz raster; gate ±10 %), −5.4 dB deep vs a +2.0 dB band median (7.4 dB prominence,
gate ≥ 2 dB). Two ~200 s release solves.

Honest observable caveat (documented in the gate): away from resonance the two-run
subtraction |S11| reads large positive values (+34 dB at 1.9 GHz) — the reference
feed's **open end** reflects differently from the patch-loaded end, so the difference
signal is large everywhere; only the **localized dip** is physical, and that is the
only thing asserted. A directional S11 (the S.12 three-probe fit at the feed) is the
clean upgrade when A.1 needs return-loss *magnitude*.

## Consequences and the A.* queue

A patch designed by closed forms resonates exactly where the engine says it should —
the antenna design chain is open. Queued: **A.1** inset-fed matching (+ directional
S11 for real return-loss numbers), **A.2** radiation pattern over the protocol (NTFF
exposure + per-face z-boundary: open top, PEC ground), **A.3** the design loop (S.11/12
machinery) on patch length vs measured resonance.
