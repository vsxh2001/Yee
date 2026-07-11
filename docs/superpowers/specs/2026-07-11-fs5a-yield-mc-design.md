# FS.5a — Monte-Carlo yield analysis (design)

**Date:** 2026-07-11 · **Track:** FS.5 (optimization maturity, `FULL-SUITE-ROADMAP.md`) · **ADR:** 0211

## Problem

Commercial suites (HFSS Optimetrics, ADS yield analysis) answer "what fraction
of manufactured boards meets spec, given dimension tolerances?" Yee has the
ingredients — closed-form models in `yee-filter`/`yee-layout`, a validated GP
surrogate in `yee-surrogate` — but no yield primitive. FS.5a is the walking
skeleton: a deterministic Monte-Carlo yield estimator, gate-certified against
the analytic normal CDF and against brute-force MC on a closed-form testcase
(the roadmap's FS.5 gate).

## Non-goals

- Space mapping (FS.5b), studio exposure (FS.5c), non-Gaussian tolerance
  distributions, correlated tolerances, importance sampling / low-discrepancy
  sequences. All are follow-ons behind the same API.

## Design

New module `yee_surrogate::yield_mc` (module name avoids the `yield` keyword):

- **Deterministic RNG, no new dependency.** `splitmix64` (public-domain
  algorithm, 5 lines) seeds and advances a 64-bit stream; standard normals via
  Box-Muller on pairs of uniforms. Same seed ⇒ bit-identical estimate on every
  platform (pure integer + f64 arithmetic).
- **`ToleranceSpec { nominal: Vec<f64>, sigma: Vec<f64> }`** — independent
  Gaussian tolerances per design parameter (the universal PCB-fab datum:
  ±3σ etch/thickness tolerance).
- **`yield_estimate(pass, &spec, n_samples, seed) -> YieldEstimate`** where
  `pass: impl FnMut(&[f64]) -> bool` is the spec test. Works identically for a
  closed-form pass function and for a GP-surrogate-backed one — the estimator
  does not know or care what model sits behind the closure.
- **`YieldEstimate { yield_frac, ci95_half_width, n_pass, n_samples }`** with
  the Wilson score interval (well-behaved at yield → 0/1 where the Wald
  interval collapses).

## Gates (all instant, non-ignored)

1. **`yield-mc-001` (analytic):** pass iff `x₀ < nominal₀ + z·σ₀` in 1-D ⇒
   yield = Φ(z) exactly. Assert the MC estimate brackets Φ(z) within its own
   95 % CI (with a 1.5× guard factor) for z ∈ {−1, 0, 0.5, 1, 2} at
   n = 100 000. Φ via Abramowitz–Stegun 7.1.26 (|ε| < 1.5·10⁻⁷, far below the
   MC noise floor).
2. **`yield-mc-002` (determinism):** same seed ⇒ bit-identical
   `YieldEstimate`; different seed ⇒ different sample stream (n_pass differs
   for at least one of several seeds).
3. **`surrogate-yield-001` (roadmap FS.5 gate):** closed-form testcase — patch
   resonance `f(L, εr) = c / (2 L √εeff)` with spec `|f − f₀| ≤ band`.
   Brute-force MC on the closed form vs the same yield run through a
   `GaussianProcess` trained on the closed form; the two estimates agree
   within combined MC noise + a surrogate-error allowance (assert
   |Δyield| ≤ 0.05 with measured value pinned in the test comment).

## Lane

`crates/yee-surrogate/**`, `docs/**`, `FULL-SUITE-ROADMAP.md`.
