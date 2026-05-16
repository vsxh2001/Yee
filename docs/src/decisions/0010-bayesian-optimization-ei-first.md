# ADR-0010: Ship single-objective Expected-Improvement BO before multi-objective NSGA-II

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

`ROADMAP.md` lists Bayesian optimisation (BO) and NSGA-II
multi-objective evolutionary search in the same bullet under the
"design-loop" theme. Phase 3.bo.0 is the first concrete piece of that
theme to land in `yee-surrogate`, sitting on top of the Gaussian-
process backend from Phase 3.gp.0/1 (see ADR-0009).

The two algorithms address overlapping but distinct user needs:

- **Bayesian optimisation, single objective.** "Find the design that
  minimises a single expensive scalar." Returns one best point.
  Builds a surrogate, picks the next evaluation by maximising an
  acquisition function (Expected Improvement, Probability of
  Improvement, UCB). The textbook small-data optimiser when each
  evaluation costs minutes of wall-time, which is exactly the regime
  a Yee user is in when sweeping a patch antenna or matching
  network.
- **NSGA-II multi-objective.** "Find the Pareto front trading three
  competing objectives — size vs bandwidth vs gain." Returns a set
  of non-dominated points. Built on population-based evolution; no
  surrogate; tolerates discrete / mixed parameters more
  comfortably than vanilla BO.

The decision is **which one ships first**. Three factors drove the
choice:

1. **What does the existing stack already give us?** Phase 3.gp.0/1
   shipped the GP plus log-marginal-likelihood hyperparameter
   optimisation. A single-objective BO loop is *exactly* "wrap the
   GP in an acquisition-function maximiser and a loop". The marginal
   work is a few hundred lines: a `BoConfig`, an `ei` function, a
   candidate-sampling step, an outer loop. Nothing in the GP stack
   directly accelerates NSGA-II, which has no surrogate
   dependency.
2. **What's the minimal validation case?** Single-objective BO has
   well-established synthetic benchmarks (deceptive 1-D functions,
   2-D Rosenbrock, the Branin / Hartmann families) with known
   minima. Multi-objective benchmarks (ZDT, DTLZ) involve Pareto-
   front-quality metrics (hypervolume, generational distance) that
   themselves need implementation and validation. The shorter
   walking-skeleton path is single-objective.
3. **Who has been asking?** The user pull on the design-loop theme
   has been "automate the patch-antenna sweep", which is a single-
   objective problem in practice (minimise S₁₁ at the design
   frequency). The multi-objective ask comes up but trails behind.

Within the single-objective lane, two sub-questions remained:

- **Acquisition function.** Expected Improvement (EI) is the
  textbook default and is provably consistent under mild
  assumptions. UCB requires a confidence-tradeoff parameter `β`
  that depends on the problem and is not free to set; PI is more
  myopic than EI in practice. Phase 3.bo.0 ships EI. UCB / PI are
  trivial follow-ups on the same `Acquisition` trait if a user case
  motivates them.
- **Candidate-maximisation strategy.** EI is non-convex over the
  search space, and finding its maximum is itself an optimisation
  problem. The conventional choices are:
  - **Random search over k candidates.** Sample `k = 1000` Latin-
    hypercube points, evaluate EI at each, take the argmax. Trivial
    to implement; `O(k)` per BO iteration; effective at low
    dimension (`d ≤ 5`).
  - **L-BFGS on EI from random restarts.** Sample 10 starts, run
    L-BFGS-B from each, take the best optimum. Standard practice
    in mature BO libraries (`BoTorch`, `Optuna`); requires an EI
    gradient. Higher implementation cost; better in higher
    dimension.
  - **CMA-ES on EI.** Population-based; no gradient. Reasonable
    middle ground; another full dependency to take on.

  Phase 3.bo.0 ships **random-search candidates** as the walking
  skeleton. The `bo_synthetic` validation case
  (deceptive 1-D function) demonstrates BO beats random search on
  the *outer* loop in 20 iterations with 1000 random EI candidates
  per iteration. At `d ≤ 5` this is empirically adequate; at higher
  dimension the EI argmax becomes the dominant source of error and
  L-BFGS becomes worth its weight.

The walking-skeleton-first principle (root `CLAUDE.md` §3) settles
the meta-question: ship the smallest end-to-end loop that
demonstrates the surrogate stack works against a published synthetic
benchmark. EI + GP + random-search EI candidates is that minimum.
NSGA-II ships when there's a real user case pulling for it.

## Decision

`yee-surrogate::bo` ships a single-objective Bayesian optimiser:

```rust
pub struct BoConfig {
    pub n_initial: usize,       // Latin-hypercube initial samples
    pub n_iterations: usize,    // BO loop iterations
    pub n_ei_candidates: usize, // random-search EI candidates / iter
    pub seed: u64,
}

pub struct BoResult {
    pub best_x: Vec<f64>,
    pub best_y: f64,
    pub history: Vec<(Vec<f64>, f64)>,
}

pub fn minimize<F>(
    f: F,
    bounds: &[(f64, f64)],
    cfg: BoConfig,
) -> Result<BoResult>
where
    F: FnMut(&[f64]) -> f64;
```

**Acquisition function:** Expected Improvement, computed from the
GP's `(predict_mean, predict_variance)` posterior. The minimum so
far is tracked in the outer loop; EI is maximised (equivalently,
negative EI is minimised) over the bounded box.

**Candidate maximiser:** random search over `n_ei_candidates` Latin-
hypercube samples per iteration, with a fixed-seed `xorshift` RNG so
runs are reproducible.

**Initial design:** Latin-hypercube sampling of `n_initial` points
before the BO loop starts.

NSGA-II and other multi-objective methods are **deferred**. Single-
objective EI is the walking skeleton; the `Acquisition` and
`CandidateMaximizer` traits in the implementation leave room for UCB,
PI, and (eventually) L-BFGS-on-EI to land as in-tree alternatives
without an API break.

## Consequences

**What becomes easier:**

- The smallest user-visible ML loop ships now. A `yee-py` user can
  call `bo_minimize` against a Python callback that runs a Yee
  solve, get a best design back in 20-50 evaluations, and validate
  the surrogate stack end-to-end against their own problem.
- The implementation exercises the GP backend (ADR-0009) in
  exactly the way it was designed to be used:
  `(predict_mean, predict_variance)` is read at thousands of
  candidate points per iteration, and the cached Cholesky factor in
  `GaussianProcess` is the hot path.
- Adding UCB or PI later is a `match` arm on a future
  `Acquisition` enum, not an API break.
- The deceptive-1-D / 2-D Rosenbrock validation cases give a clean
  regression gate: BO converges to the known minimum within budget,
  and beats random search by a wide margin on the deceptive
  problem.

**What becomes harder:**

- Users with multi-objective designs (size vs bandwidth vs gain;
  S₁₁ vs gain vs cost) cannot get a Pareto front from Yee in this
  phase. Their workaround is to scalarise — combine objectives into
  a weighted sum and let single-objective BO chase the weighted
  scalar — which loses Pareto-front information but recovers a
  point optimum.
- Random-search EI maximisation degrades as dimension rises. At
  `d ≤ 5` the synthetic-benchmark performance is good; at `d ≥ 10`
  the EI argmax becomes increasingly noisy and the BO convergence
  rate suffers. Users hitting this will see slow progress and
  should be told (in `yee-surrogate`'s README) that the candidate
  count must rise super-linearly with dimension, or that an L-BFGS
  candidate maximiser is needed before BO is competitive at scale.
- L-BFGS-on-EI requires an analytic EI gradient (the closed-form
  expression in terms of GP posterior gradients is standard but
  non-trivial). This is a real implementation cost that is deferred
  to a future phase.

**What's now closed off:**

- Reaching for an external BO crate (`bayesopt`, `egobox-ego`) in
  Phase 3.bo.0. The hand-rolled loop sits directly on top of the
  in-tree GP from ADR-0009 and avoids the same dependency-weight
  argument that ADR justified.
- Bundling NSGA-II into the same phase. Multi-objective evolution
  is its own walking-skeleton project with its own validation
  matrix.
- Exposing the acquisition function as a generic numeric
  programming surface (e.g. arbitrary user callbacks). The trait
  is internal in Phase 3.bo.0; `bo_minimize` is the only public
  entry point.

## References

- `crates/yee-surrogate/src/bo/` — the implementation.
- `crates/yee-surrogate/src/bo/ei.rs` — Expected Improvement
  acquisition.
- `crates/yee-surrogate/src/bo/lhs.rs` — Latin-hypercube initial
  sampling.
- `crates/yee-surrogate/validation/bo_synthetic.rs` — BO beats
  random search on a deceptive 1-D function within 20 iterations.
- `crates/yee-py/python/yee/bo.py` — Python wrapper around
  `bo_minimize` (Phase 3.frontend.0).
- D. R. Jones, M. Schonlau, W. J. Welch, "Efficient Global
  Optimization of Expensive Black-Box Functions", *Journal of Global
  Optimization*, vol. 13, no. 4, pp. 455-492, 1998. (The original
  EI / EGO paper.)
- K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, "A Fast and Elitist
  Multiobjective Genetic Algorithm: NSGA-II", *IEEE Trans.
  Evolutionary Computation*, vol. 6, no. 2, pp. 182-197, 2002. (The
  NSGA-II paper, for the deferred follow-up.)
- ADR-0009 — the GP backend this BO consumes.
