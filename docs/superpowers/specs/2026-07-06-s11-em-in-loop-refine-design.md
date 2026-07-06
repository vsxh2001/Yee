# S.11 / F1.2.1.0 — EM-in-the-loop dimension refinement (walking skeleton)

**Date:** 2026-07-06
**Phase:** S.11 (ENGINE-STUDIO-ROADMAP) = F1.2.1.0 (the first slice of the filter
roadmap's F1.2.1 "EM-in-the-loop refinement"). Consumes the S.10 measurement fidelity.
**Plan:** `docs/superpowers/plans/2026-07-06-s11-em-in-loop-refine.md`

## Problem

S.10 left a clean, attributable design-side error: the synthesized LPF's measured
cutoff sits 15 % below design because the closed-form dimensions are seeds (staircased
high-Z widths, un-de-embedded junctions). "Designing with the engine" means closing
that error with the EM measurement in the loop. The full F1.2.1 program (surrogate-BO
over k/Qe for BPFs) is large; the walking skeleton is **one scalar knob, one correction
step** — the loop structure, not the optimizer.

## Design

The knob: the **design frequency handed to the dimensioner**. A stepped-impedance
LPF's electrical lengths scale as `1/f_c`, so to first order the measured cutoff scales
with the synthesis frequency. One multiplicative correction:

```
f_c'  =  f_c_design · (f_c_design / f_measured)
```

re-synthesize the layout at `f_c'`, re-verify, and the measured cutoff should land near
the design target. This is a secant step on a monotone scalar map — no optimizer crate,
no surrogate; those arrive with the BPF (multi-knob) slice.

Scenario sized for a 4-solve budget: **N = 3 Butterworth, f_c = 2 GHz, FR-4**, the
S.9/S.10-certified measurement stack (CPML-xy walls, aperture ports), smaller margins,
~8500 steps. Self-contained gate in `yee-filter/tests` (the S.8 gate stays untouched —
its certified numbers must remain reproducible).

## Validation gate — engine-refine-001

1. The seed verify measures a cutoff meaningfully off-design (it is; that is the
   premise — assert |err_seed| ≥ 5 % so the gate is honest about having something to
   fix).
2. The refined verify's cutoff error is **at most half** the seed error.
3. The refined error is within the walking-skeleton band |err| ≤ 10 %.
4. Both measured tables printed for the record; final numbers recorded in ADR-0188.

## Non-goals

Multi-knob refinement (per-section lengths/widths), BO/surrogate integration
(`yee-surrogate` exists — wired in when the BPF slice needs it), BPF k/Qe extraction,
convergence iteration beyond one correction step (add when a scenario needs it).
