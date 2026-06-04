# F1.2.1.1 — EM coupling design-curve gap correction — design

**Date:** 2026-06-04
**ADR:** [ADR-0159](../../src/decisions/0159-f1-2-1-1-em-coupling-design-curve-gap-correct.md)

## Problem

The full-wave 3-pole filter S21 floors ~−27 dB because the analytic dimensioner sized the gaps with
the impedance-k, which diverges ~37 % from the realized resonant-k (k_imp≠k_eps, K2). Finer mesh was
refuted. The fix: correct each gap using the EM-measured coupling. The FEM K(gap) was shown smooth +
monotone, so a per-gap root-find (Hong-Lancaster full-wave coupling-design) converges cleanly.

## Goal

B1: `correct_gap_fem_k` — root-find the gap realizing a target K via `coupled_resonator_k`. B2: apply
to the filter + re-grade the corrected S21. B3 (fallback): ASM-Broyden if gap interactions block B2.

## Non-goals

The fixed-box_w change (not needed — K(gap) smooth), Qe-feed dimensioning (ADR-0156), finer-mesh as
the lever (refuted), full ASM unless B2 falls short.

## Architecture (B1 — this spec's brick)

- **Home:** `crates/yee-fem/src/coupled_resonator_k.rs` (add `correct_gap_fem_k`) or a sibling module;
  gate in `crates/yee-fem/tests/`.
- **`correct_gap_fem_k(geom_template, k_target, gap_lo, gap_hi, tol, max_evals) -> GapCorrection`**:
  a 1-D root-find of `f(gap) = coupled_resonator_k({geom_template, gap_s: gap}).k_fem − k_target` on
  `[gap_lo, gap_hi]`. K(gap) is monotone-decreasing (confirmed), so use **bisection** (robust to the
  occasional extraction outlier) or secant; return the gap, the achieved K_fem, the eval count, and a
  `converged` flag. Each eval = one `coupled_resonator_k` FEM sweep (~280 s).
- **Robustness:** guard non-finite/`!peaks_resolvable` K (skip/penalize, don't crash); bisection
  tolerates a single noisy eval better than secant. Bracket `[gap_lo,gap_hi]` from the smooth-curve
  range (e.g. [1.0, 4.0] mm, clamped to GAP_MIN/MAX).

## Gate `fem-coupling-correct-001` (`#[ignore]`'d, `--release`, boxed)

Non-circular: pick a target K (a synthesis `FBW·m`-style value reachable on the K(gap) curve, e.g.
0.040) and a deliberately-mis-dimensioned SEED gap (the analytic-k gap, off-target ≥10 %). Run
`correct_gap_fem_k` → assert: `converged`, `|K_fem(g*) − k_target|/k_target < 0.08` (reachable on the
smooth curve, vs F1.2.1.0's 13.4 %), in ≤ ~5-6 FEM evals, and the corrected gap differs from the seed
by a real margin. Print the root-find trajectory (gap, K_fem per eval). Heavy → orchestrator runs it.

## Testing

- Fast unit test (non-ignored, no FEM): the root-find logic on a synthetic monotone f(gap) (e.g. a
  closed-form decreasing function) converges to the root in ≤ N steps — validates the bisection/secant
  + bracketing without the FEM cost.
- The heavy gate `#[ignore]`'d + a `fem-coupling-correct-001` `--release` CI step; the orchestrator
  runs the heavy ~5-eval sweep (misfire-split: agent writes code, I run).

## Risks

- Extraction outliers in K(gap) (F1.2.1.0 saw one) → bisection (not secant) for robustness; the gate's
  ≤6-eval bound catches a non-converging case honestly. (med-low)
- B2 (next brick) gap-interaction: the isolated-pair K vs the in-filter coupling may differ → B2 grades
  honestly + B3 (ASM) escalates. (research-open, B2's risk, not B1's)
