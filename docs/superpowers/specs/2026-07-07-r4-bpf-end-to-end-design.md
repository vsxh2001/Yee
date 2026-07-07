# R.4 — BPF end-to-end (F1.2.1 core) + multi-knob surrogate BO

**Date:** 2026-07-07
**Track:** RF-TOOL-ROADMAP R.4 (= the deferred FILTER-DESIGN-ROADMAP F1.2.1 core)
**Related:** ADR-0109 (hairpin dims, F1.2.2 — documented the uniform-gap layout
limitation and the deferred qe→tap), ADR-0188/0189 (the refine loop + the
directional observable it needed), ADR-0196 (R.3), `coupling_matrix_s_001`
(the validated coupling-matrix→S-parameter reference).

## Problem

The hairpin BPF chain stops one step short of end-to-end: synthesis produces
N−1 *distinct* gaps and an external Q, but the convenience layout collapses
the gaps to their mean (`HairpinParams` carries a single `coupling_gap_m`) and
taps the feed at a *placeholder* `arm_length/3` with an explicit "qe→tap is
deferred to F1.2.1" note. No full-wave gate verifies a synthesized BPF against
its coupling-matrix response, and no multi-knob refinement exists (the S.11/12
loop is single-knob secant).

## Scope — two shippable increments

### R.4a — per-section geometry + qe→tap + full-wave verify

1. **Per-section hairpin layout** (`yee-layout`): `HairpinSectionParams` with
   `gaps_m: Vec<f64>` (length N−1) + `hairpin_bpf_sections`, accumulating each
   resonator's x-base with its own gap. The existing `hairpin_bpf` (and its
   committed `geo-003` gate) stays untouched.
2. **qe→tap** (`yee-filter::dimension`): the tapped half-wave-resonator
   external-Q relation (Hong & Lancaster ch. 6 tapped-line coupling; energy
   derivation): for a λ/2 resonator of impedance `Zr`, length `L = 2·arm`,
   tapped by a `Z0` feed at distance `t` from the **open end**,

   `Qe(t) = (π/2)·(Z0/Zr) / cos²(π·t/L)`

   (voltage antinode at the open ends, null at the fold; W = C_l·V0²·L/4,
   P = V_t²/2Z0, Qe = ω0·W/P with ω0 = πv/L). Inverted:

   `t = (L/π)·acos( sqrt( (π/2)·(Z0/Zr) / qe ) )`,   `qe = g0·g1/FBW`

   with a hard error when `qe < (π/2)(Z0/Zr)` (tap cannot couple that
   strongly) or when `t` falls off the arm. `HairpinDimensions` gains
   `tap_offset_m`; `dimension_hairpin_layout` uses the per-section generator
   and the real tap.
3. **Gate `engine-bpf-verify-001`** (`yee-filter/tests/engine_bpf_verify.rs`,
   `#[ignore]`, release, ~2 solves): synthesize a low-order Chebyshev/
   Butterworth hairpin BPF (FBW sized so every gap ≥ 2·dx — coupled gaps are
   the resolution driver), voxelize, run DUT + straight-reference through the
   job protocol (the `engine_lpf_verify` skeleton), extract directional |S21|
   (3-probe standing-wave fit at both ports — a resonant DUT reflects hard),
   and hold it against **`coupling_matrix_s_params`**: measured centre
   frequency and −3 dB bandwidth vs the coupling-matrix curve within
   walking-skeleton tolerances (set from the first honest measurement, S.8
   style), plus stopband rejection. Tolerances tighten in R.4b, not here.

### R.4b — multi-knob surrogate BO

`yee_surrogate::bo::minimize` (closure objective, d-dim bounds) over ≥2 knobs
(synthesis f0 + a global gap scale, the two dominant error axes) driving the
R.4a pipeline; objective = weighted |f0_meas − f0_spec| + |BW_meas − BW_spec|.
Gate `engine-bpf-bo-001`: BO closes both to spec tolerance within a fixed
solve budget, beating the synthesis seed. Budget: each objective call is one
release solve, so `n_initial + n_iters ≲ 10`.

## Out of scope

Fold-spacing self-coupling refinement, tap-point reactance compensation,
conductor loss (R.0b), full-wave k/Qe extraction (F1.1b.2) — the BO closes
the *measured response* to spec, which is the tool-level contract.

## Consequences

The filter path becomes what the antenna path already is: synthesize →
measure → close the loop on the engine, now with the industry-standard BPF
topology and multi-knob optimization. Combined with R.2's Touchstone export,
a designed BPF leaves the tool as a verified `.s2p`.
