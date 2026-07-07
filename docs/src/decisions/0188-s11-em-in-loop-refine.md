# ADR-0188: S.11 / F1.2.1.0 EM-in-the-loop refinement — loop landed, convergence blocked on a named artifact

**Status:** Accepted (diagnostic landed; convergence gate deferred to S.12)
**Date:** 2026-07-06
**Related:** ADR-0187 (the measurement this consumes), ADR-0129/0131 (the SAME artifact
class in F2.3 and its fix), FILTER-DESIGN-ROADMAP F1.2.1.
**Spec:** `docs/superpowers/specs/2026-07-06-s11-em-in-loop-refine-design.md`

## What landed

The first **closed design loop** on the engine, mechanically complete:
synthesize (`f_c` knob) → voxelize → verify over the job protocol → fit the measured
|S21| curve → correct the synthesis frequency (secant on the measured map) → repeat.
`engine_lpf_refine.rs` runs it end-to-end (manual diagnostic, deliberately not in CI).
The observable is a **whole-curve fitted Butterworth cutoff** (least-squares in dB over
the 69-point band, deep stopband excluded) — itself the survivor of two rejected
observables (threshold crossings are metric-dependent on these skirts: two defensible
detectors read the same board 1.7 and 2.9 GHz).

## Why the convergence gate does NOT ship (measured)

Secant iteration on the N = 5 scenario:

| iter | synth f_c | fitted cutoff | fit rms |
|---|---|---|---|
| 0 | 2.000 GHz | 1.460 GHz | **5.51 dB** |
| 1 | 2.740 GHz | 2.530 GHz | 1.30 dB |
| 2 | 2.373 GHz | 1.740 GHz | **5.45 dB** |
| 3 | 2.494 GHz | 2.330 GHz | 1.22 dB |

The loop oscillates (−27 % → +26.5 % → −13 % → +16.5 %) because the map is non-smooth:
the fit rms alternates between ~1.3 dB (genuinely Butterworth-shaped curves) and
~5.5 dB (curves carrying bumps at ~1 GHz spacing — matching the **port-to-port round
trip** over the ~75 mm board). The S.10 aperture ports brought passband S11 to −9 dB,
but two −9 dB ports still leave multi-reflection ripple riding on the skirt, and a
scalar-cutoff observable cannot converge on it. Also recorded: an N = 3 mini-board was
rejected outright (short body leaks an over-substrate air path; +3..+5 dB stopband
bumps).

## The named fix (S.12)

This is the artifact class F2.3 already solved (ADR-0129/0131, "standing-wave /
over-unity"): fit the **standing wave** on the line with three probes,
`V(x) = a·e^{−jβx} + b·e^{+jβx}`, and take the forward amplitude — a directional |S21|
immune to the reflected wave. `yee-voxel` has `fit_standing_wave` already; S.12 ports
that observable to the S-parameter path (pure post-processing on three probe series —
no engine change), re-runs this loop, and restores the convergence asserts
(final error ≤ half seed and ≤ 10 %).

## Consequences

Loop plumbing is proven and stays; the refinement result is honestly *not claimed*.
S.12 (directional S-parameters) is the next increment and unblocks both this gate and
higher-fidelity absolute S-parameters generally.
