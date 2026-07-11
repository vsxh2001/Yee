# ADR-0211: FS.5a — Monte-Carlo yield analysis (deterministic, model-agnostic)

**Date:** 2026-07-11 · **Status:** accepted · **Track:** FS.5 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-11-fs5a-yield-mc-design.md`

## Context

Commercial parity (HFSS Optimetrics, ADS yield) requires answering "what
fraction of manufactured devices meets spec under fabrication tolerances?"
Yee had the models (closed forms, GP surrogate, full engine) but no yield
primitive. FS.5a is the walking skeleton.

## Decision

`yee_surrogate::yield_mc` (module name sidesteps the `yield` keyword):

1. **Model-agnostic estimator.** `yield_estimate(pass, &ToleranceSpec,
   n_samples, seed)` draws independent Gaussian parameter perturbations and
   counts passes; `pass: FnMut(&[f64]) -> bool` may wrap a closed form, a
   trained `GaussianProcess`, or a full engine run. The estimator never
   learns what model it is sampling.
2. **Deterministic in-crate RNG, no new dependency.** `splitmix64`
   (public-domain; verified in-module against the C reference vector for
   seed 1234567) + Box-Muller. Pure wrapping-integer + f64 arithmetic ⇒
   the same seed reproduces the estimate bit-for-bit on every platform,
   so gate numbers are pinnable and CI-stable. `rand` was rejected: one
   consumer, and its float streams are not guaranteed stable across major
   versions.
3. **Wilson 95 % score interval**, not Wald: yield analysis operates near
   yield → 1 where the Wald interval collapses to zero width (unit test
   pins non-collapse at 0/1000 and 1000/1000 passes).

## Measured gates (all instant, non-ignored, `crates/yee-surrogate/tests/yield_mc.rs`)

- **`yield-mc-001`** — pass iff x < z, x ~ N(0,1) ⇒ yield = Φ(z)
  (Abramowitz–Stegun 7.1.26 reference). MC at n = 100 000 brackets Φ(z)
  within 1.5× its own CI for z ∈ {−1, 0, 0.5, 1, 2}.
- **`yield-mc-002`** — same seed ⇒ bit-identical `YieldEstimate`;
  different seeds diverge.
- **`surrogate-yield-001`** (the roadmap FS.5 gate) — patch-resonance
  closed form f = c/(2L√ε_eff), FR-4 nominal L = 29 mm ± 0.1 mm,
  ε_r = 4.4 ± 0.05, spec ±40 MHz. GP trained in σ-normalized coordinates
  on a 9×9 ±3σ grid. **Measured: brute-force 0.9721 vs surrogate 0.9720
  (Δ = 1e-4)**, both on the analytic linearization 2Φ(2.2)−1 ≈ 0.972.
  Same seed for both runs, so the sample streams are identical and the
  delta isolates pure surrogate error at the spec boundary. Assert
  |Δ| ≤ 0.02.

## Consequences

- Yield over *engine-verified* responses is now a composition, not new
  machinery: train the GP on engine samples (the R.4 BO loop already
  produces them), wrap `predict_mean` in a `pass` closure, call
  `yield_estimate`. FS.5c exposes exactly that in the studio.
- Follow-ons behind the same API: correlated tolerances (Cholesky of a
  covariance in `ToleranceSpec`), non-Gaussian marginals, importance
  sampling for rare-failure yields.
