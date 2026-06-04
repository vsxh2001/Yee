# F1.2.1.1 — EM coupling design-curve gap correction — plan

**Spec:** [2026-06-04-f1-2-1-1-em-coupling-design-curve-design.md](../specs/2026-06-04-f1-2-1-1-em-coupling-design-curve-design.md)
**ADR:** [ADR-0159](../../src/decisions/0159-f1-2-1-1-em-coupling-design-curve-gap-correct.md)
**Fork from:** main (`d97e5a9` or later).

## Brick B1 — `correct_gap_fem_k` + `fem-coupling-correct-001` gate

1. `crates/yee-fem/src/coupled_resonator_k.rs` (or a new `coupling_correct.rs` re-exported from lib):
   `pub struct GapCorrection { gap_m, k_fem, k_target, n_evals, converged }`; `pub fn
   correct_gap_fem_k(base: &CoupledResonatorGeom, k_target: f64, gap_lo: f64, gap_hi: f64, tol_frac:
   f64, max_evals: usize, n_pts: usize) -> GapCorrection`. BISECTION on `g ↦
   coupled_resonator_k(&{..base, gap_s: g}, n_pts).k_fem − k_target` (K monotone-DEC, so k>target ⇒
   gap too small ⇒ move right). Skip non-finite / !peaks_resolvable evals defensively. Re-export.
2. Fast unit test (non-ignored, NO FEM): bisection on a synthetic monotone-decreasing closure (e.g.
   `f(g)=0.09−0.02·g·1e3` style) converges to the root within tol in ≤ ceil(log2(range/tol)) steps.
3. Gate `crates/yee-fem/tests/fem_coupling_correct_001.rs` (`#[ignore]`'d + `--release`): base =
   the K1/K2 probe geom; `k_target = 0.040`; seed gap = the analytic-k gap (off ≥10 %); bracket
   [1.0e-3, 4.0e-3]; tol_frac 0.08; max_evals 6; n_pts 61. Assert `converged`, `|k_fem−k_target|/
   k_target < 0.08`, `n_evals ≤ 6`, corrected gap ≠ seed by a real margin. Print the trajectory.
4. CI: add a `fem-coupling-correct-001` step to the `fem-eigen-gate` `--release` job.

## B2 (next brick, after B1 merges) — apply to the filter + re-grade

5. In `microstrip_filter_s21.rs`: replace the analytic dimensioned gaps with B1-corrected gaps (root-
   find each target_k[i]), re-build the filter, re-run the S21 sweep, re-grade vs the Cheb mask.
   HONEST gate: record corrected in-band peak / mask margin vs the −27 dB baseline; assert the
   measured improvement; if it clears the mask, assert it.

## B3 (fallback, only if B2 short) — ASM-Broyden multi-D over all gaps (per [[reference-em-in-loop-space-mapping]]).

## Dispatch (B1 now)
- SPLIT: a code-only agent writes B1 (`correct_gap_fem_k` + the fast unit test + the gate + CI step),
  verifies it COMPILES (clippy/fmt + the fast unit test, all fast) + commits; does NOT run the heavy
  gate. The orchestrator runs the heavy `fem-coupling-correct-001` (~5-6 FEM evals, ~25 min) boxed,
  reads the trajectory, verifies converged <8%.
- Then: code-reviewer (NEVER self-review — gate honesty: non-circular, real convergence, bisection
  correct), fix P0/P1, merge `--no-ff` from `/home/hadassi/Code/Yee`, push, ROADMAP/memory, cleanup.
