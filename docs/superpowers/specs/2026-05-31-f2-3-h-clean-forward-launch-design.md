# Filter Phase F2.3-h — clean forward-wave launch — Design Spec

**ADR:** ADR-0133 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

F2.3-g (ADR-0132) made the de-embed physical (no over-unity) but is
launch/probe-floor-limited: the PEC-box soft source reflects (weak forward launch),
and the output reference region sees no clean propagating forward wave (β_out=0 at
some gate freqs, |b₂| at the floor). So the S21 — and the "notch at f0" — is not
yet trustworthy. Fix the forward-wave launch + output probe.

## Goal

Well-resolved travelling-wave amplitudes `a₁` (incident at input) and `b₂`
(transmitted at output) — β>0 at all gate freqs, `b₂` well above the floor — so the
F2.3 S21 is trustworthy and the notch-at-f0 is disambiguated (band-pass / real
inversion / research wall). Keep the strict gate.

## Method

In the F2.3 driver (`crates/yee-voxel/src/lumped_sim.rs`), building on the F2.3-g
PEC-box 2-point de-embed:

- **Clean forward reference `a₁`:** adopt the `run_line_eeff` time-gated
  incident-wave pattern (ADR-0108) — launch the drive into a long-enough lead-in so
  the incident forward wave is cleanly time-gated *before* the first reflection,
  giving a trustworthy `a₁`. (And/or a directional source: a soft source with a
  PEC/absorbing backing, or a TF-SF-style launch, so it injects predominantly
  forward, reducing the near-pure-standing-wave problem.)
- **Output probe:** lengthen the line + place the output reference region clear of
  the PEC end wall and any evanescent zone so a propagating forward wave is clean
  there (β_out>0 at all gate freqs, `b₂` above the floor).
- **DUT response:** keep the CW steady-state (the high-Q tanks must ring up), but
  reference it to the trustworthy forward `a₁`. A hybrid (time-gated incident `a₁`
  reference + CW-settled DUT/thru `b₂`) is acceptable; document the scheme.
- Keep the aperture-port placement. dx = 0.4 mm; bounded freq set; keep runs minutes
  (do NOT reintroduce finer-dx / multi-hour runs).

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `run_board_solve`: the clean forward launch (time-gated incident reference and/or
  directional source) + longer line + output-probe placement so `a₁`/`b₂` are
  trustworthy. Keep the 2-point de-embed + aperture ports. Module doc updated. Bring
  the branch up to `main` first. Do NOT edit yee-fdtd/yee-filter.

## DoD (machine-checkable; container-iterated)

1. fmt + clippy -D warnings on yee-voxel → exit 0.
2. The de-embed run + REPORTED: are `a₁`/`b₂` now well-resolved (β>0 at all gate
   freqs, `b₂` above the floor, thru not degenerate)? The full |S21| sweep.
3. Disambiguation verdict (all honest, no gate weakening):
   - clean band-pass ≥20 dB → `fdtd_lumped_001` GREEN → ships.
   - clean inverted response (notch @f0) persists → real topology inversion → next.
   - still floor-degenerate / inconclusive → surface the cumulative measurement
     wall.
4. Runs bounded (minutes, dx=0.4 mm); no regression to other yee-voxel gates.

## Out of scope

The topology-inversion fix (next); the sub-cell port correction; the studio Verify
stage; finer-dx (proven worse, F2.3-e).

## Why

It converts F2.3-g's physical-but-floor-limited measurement into a trustworthy S21,
turning the ambiguous "notch at f0" into a definitive classification — either
shipping EM-sim, exposing a cheap topology bug, or proving the research wall for a
maintainer decision.
