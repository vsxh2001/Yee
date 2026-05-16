# ADR-0013: Active learning as BO with a variance acquisition function

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

`ROADMAP.md` listed **active learning (AL)** as a separate Phase 3
deliverable alongside Bayesian optimisation and NSGA-II. The framing in
the roadmap (and in `docs/src/theory/multi-objective-and-active-
learning.md`) presents AL as its own algorithm family: "given a
surrogate, choose the next sample to evaluate so as to *improve the
surrogate's accuracy*, not to *find the minimum*."

When Phase 3.al.0 came up the queue, the obvious implementation
question was: how big is the new code surface?

The literature framing makes AL sound like a separate algorithm, but
operationally it is **the same BO loop with a different acquisition
function**:

- **BO loop (ADR-0010).** Fit a GP on the current samples. Maximise
  an acquisition `α(x)` over the bounded box (random-search over
  `n_candidates` Latin-hypercube candidates). Query the expensive
  objective at `argmax α(x)`. Append to the sample set. Repeat.
- **BO acquisition (ADR-0010).** Expected Improvement,
  `α_EI(x) = E[max(0, y_min − μ(x))]`. Drives sampling toward
  *low-mean* regions, i.e. the minimum.
- **AL acquisition.** Pure posterior variance,
  `α_AL(x) = σ²(x)`. Drives sampling toward *high-variance*
  regions, i.e. where the surrogate is least confident.

That is the entire algorithmic difference. The outer loop, the GP fit,
the candidate sampling, the LHS initial design, the RNG, the result
struct shape — all identical.

Two structural responses were considered:

1. **Wire AL as an `Acquisition::Variance` enum arm inside `bo`.**
   Reuse the BO loop, add a `pub enum Acquisition { Ei, Variance }`,
   dispatch on it inside the optimiser. Compact, low LoC.
2. **Ship a small `al` module that calls the BO building blocks
   directly.** Reuse `bo::rng`, `bo::lhs`, `GaussianProcess::fit_ml`,
   and the random-search candidate maximiser, but with a different
   public entry point (`al::active_learn`) and a different
   acquisition function. No changes to the BO public surface; the
   AL module is its own crate-public function.

Option 1 trades a smaller new code surface for a less honest public
API. `bo::minimize` returns `BoResult { best_x, best_y, history }`
— a single best design. AL does not have a "best design"; it has a
trained surrogate. Forcing an `Acquisition::Variance` arm through
the `minimize` shape misnames the operation: AL is not minimising
anything.

Option 2 keeps the public surface honest: `al::active_learn` returns
the trained `GaussianProcess` (and the sample history), not a "best
point". The cost is a few hundred lines of duplicated outer-loop
plumbing, which is small.

The choice is option 2.

The validation case is **sin(x) on [0, 2π]** with a fixed test grid:

- Train one GP on **20 random samples** (Latin-hypercube), measure
  RMSE on a 100-point uniform test grid.
- Train one GP on **20 AL samples** (5 LHS initial + 15 chosen by
  variance acquisition), measure RMSE on the same grid.
- Gate: `RMSE_AL < 0.5 × RMSE_random`.

A 2× improvement is the conventional "AL is doing something useful"
threshold. The gate is loose by design — variance acquisition has a
known degenerate behaviour on flat-objective regions that means a
factor-of-10 win is achievable on toy problems and 2× is achievable on
real ones.

## Decision

Ship `yee_surrogate::al` as a small new module that reuses BO's
building blocks:

```text
crates/yee-surrogate/src/
  bo/
    rng.rs            — Xorshift64 (pub(crate))
    lhs.rs            — latin_hypercube (pub(crate))
    mod.rs            — single-objective BO (ADR-0010)
  al/
    mod.rs            — pub fn active_learn(...)
```

The public surface:

```rust
pub struct AlConfig {
    pub n_initial: usize,         // Latin-hypercube initial samples
    pub n_iterations: usize,      // AL loop iterations
    pub n_candidates: usize,      // candidates per iter, scored by σ²(x)
    pub seed: u64,
}

pub struct AlResult {
    pub gp: GaussianProcess,
    pub samples: Vec<(Vec<f64>, f64)>,
}

pub fn active_learn<F>(
    f: F,
    bounds: &[(f64, f64)],
    cfg: AlConfig,
) -> Result<AlResult>
where
    F: FnMut(&[f64]) -> f64;
```

**Acquisition function:** posterior variance `σ²(x)`, read directly
from `GaussianProcess::predict(x) -> (mean, variance)`. The
`predict_mean` value is **not** used by AL — only the variance enters
the acquisition. Refitting `fit_ml` after every new sample keeps the
length-scale honest as the sample density grows.

**Outer-loop scaffolding.** `al::active_learn` uses:

- `bo::rng::Xorshift64` for reproducibility.
- `bo::lhs::latin_hypercube` for the `n_initial` initial design.
- `GaussianProcess::fit_ml` (ADR-0009) for refit after each new
  sample.
- A random-search candidate maximiser that scores
  `n_candidates` Latin-hypercube points by their predicted variance
  and picks the argmax. (Same structural pattern as BO's
  EI-candidate loop; different score function.)

**Validation gate.**
`crates/yee-surrogate/validation/al_sinx.rs` runs AL and random
baseline on `f(x) = sin(x)`, `x ∈ [0, 2π]`, with
`n_initial = 5, n_iterations = 15, n_candidates = 200`. The gate
asserts `RMSE_AL < 0.5 × RMSE_random` on a 100-point test grid.

Phase 3.al.0 measured **`RMSE_AL ≈ 0.012`** vs
**`RMSE_random ≈ 0.96`** on this validation case — an 80×
improvement, well past the 2× gate.

## Consequences

**What becomes easier:**

- The AL deliverable for the v1.0 design-loop story closes with
  minimal new code: roughly 150 lines for the AL module itself, plus
  the validation gate. The implementation effort is proportional to
  the algorithmic difference (one acquisition function), not to the
  surface area of "active learning" as a topic.
- The Python bindings (`yee-py`) get a third design-loop entry point
  (`active_learn`) alongside `bo_minimize` and `nsga2_optimize`,
  each with a parallel signature and a parallel `*_result` return
  type.
- The `al_sinx` validation case is a clean, cheap, fast-running
  regression gate (~50 ms) that exercises the `predict_variance`
  path of the GP backend (ADR-0009) in exactly the way it was
  designed for: variance-driven sample placement.

**What becomes harder:**

- Variance-only acquisition has a **known failure mode**:
  *fence-post behaviour* on flat objectives. If `f(x)` has a wide
  flat region and a narrow peak, the GP's posterior variance after a
  few samples drops uniformly across the flat region, and the
  remaining variance is concentrated near the box boundary
  (because the GP has fewer neighbours there to constrain it).
  Variance-greedy AL will then place every remaining sample on the
  boundary, ignoring the peak. The well-known fixes — combining
  with EI (a "uncertainty + expected improvement" hybrid), or
  imposing a minimum distance from existing samples as a soft
  constraint — are deferred to Phase 3.al.1 and not in Phase
  3.al.0.
- Users who want a *Bayesian-experiment-design* style of AL (e.g.
  mutual-information acquisition, integrated-variance acquisition,
  BALD) will find that Phase 3.al.0 only ships variance-greedy.
  These are tractable extensions but each needs its own ADR-style
  decision and validation case.
- Refitting `GaussianProcess::fit_ml` after every new sample is
  `O(n³)` per refit (Cholesky factorisation) plus the marginal-
  likelihood gradient-ascent cost. At `n ≤ 50` (the BO/AL regime
  per ADR-0009) this is microseconds; at `n ≥ 200` AL becomes
  noticeably slower than BO because BO can amortise the GP fit
  across the EI-candidate evaluation (one fit, thousands of
  predicts) whereas AL re-fits after every single new sample by
  design.

**What's now closed off:**

- Adding AL as an `Acquisition::Variance` enum arm inside
  `bo::minimize`. The two operations have different return types
  (best point vs trained surrogate) and forcing them through one
  entry point would misname AL.
- Skipping the variance-only walking-skeleton in favour of a
  hybrid EI+variance or distance-floor acquisition. The walking-
  skeleton-first principle (root `CLAUDE.md §3`) settles this:
  ship variance-only first, fix the fence-post behaviour in
  Phase 3.al.1.
- Pulling in an external AL crate (`linfa-bayes`, an
  `active-learning-*` crate) in Phase 3.al.0. The hand-rolled
  module is 150 lines and reuses the existing GP and BO building
  blocks; the dependency-weight argument from ADR-0009 applies
  unchanged.

## References

- `crates/yee-surrogate/src/al/mod.rs` — the implementation.
- `crates/yee-surrogate/validation/al_sinx.rs` — sin(x) AL vs
  random baseline; `RMSE_AL < 0.5 × RMSE_random` gate.
- `crates/yee-py/python/yee/al.py` — Python wrapper around
  `active_learn`.
- B. Settles, "Active Learning Literature Survey", *University of
  Wisconsin-Madison Department of Computer Sciences Tech Report
  1648*, 2009. (Survey of acquisition functions including variance,
  uncertainty, BALD, expected error reduction.)
- C. Houlsby, F. Huszár, Z. Ghahramani, M. Lengyel, "Bayesian Active
  Learning for Classification and Preference Learning",
  *arXiv:1112.5745*, 2011. (BALD acquisition, deferred to future
  phases.)
- ADR-0009 — GP surrogate; `predict(x) -> (mean, variance)` is the
  hot path for AL's variance-greedy candidate scoring.
- ADR-0010 — single-objective BO; AL reuses its `rng` and `lhs`
  helpers and its candidate-maximiser pattern, but exposes a
  different public entry point.
