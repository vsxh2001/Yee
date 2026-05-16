# Phase 3.bo.0 — Bayesian optimization with Expected Improvement

**Status:** Draft  
**Owner:** TBD  
**Phase:** 3.bo.0  
**Depends on:** Phase 3.gp.0 (shipped, `GaussianProcess`), 3.gp.1 (shipped, `fit_ml`)  
**Blocks:** Phase 3.bo.1 (multi-objective NSGA-II), Phase 3.al.0 (active learning)

## Assumption being challenged

ROADMAP Phase 3 lists "Bayesian optimization, NSGA-II for multi-objective" as one bullet, implying both ship together. Wrong assumption: single-objective BO is a much smaller deliverable than NSGA-II, and gates the multi-objective work anyway because the acquisition machinery is shared. Ship single-objective BO with Expected Improvement first; NSGA-II is a follow-up phase.

The GP surrogate is already in main. It predicts mean **and variance**. EI consumes both. The minimum-viable BO loop is:

1. Sample N₀ random initial points.
2. Fit GP via `fit_ml`.
3. Maximize EI over a candidate set (random or grid) to pick next point.
4. Evaluate true objective at next point, append to training set.
5. Refit GP, repeat until budget exhausted.

This is ~150 LOC plus a test.

## Scope

In:
- New module `crates/yee-surrogate/src/bo.rs`
- Expected Improvement acquisition: `ei(mean, std, f_best, xi)` where `xi` is the exploration parameter (default 0.01)
- `BayesianOptimizer` struct: takes a black-box `Fn(DVector<f64>) -> f64` + bounds + budget
- Random-search candidate proposal (1024 candidates per iter; pick the one with max EI)
- Validation test: minimize `f(x) = (x - 3)² + sin(5x)` on `x ∈ [0, 6]`, budget 20 evaluations, must beat random search by ≥2× on regret

Out:
- L-BFGS for EI maximization (random-search candidates are enough at d ≤ 5 for the validation case)
- Multi-objective Pareto fronts
- Constrained optimization
- Batch BO

## Public API

```rust
//! Bayesian optimization with Gaussian process surrogate + Expected Improvement.

pub struct BoConfig {
    pub n_initial: usize,        // random initial samples (default 5)
    pub n_iters: usize,          // BO iterations after initial (default 20)
    pub n_candidates: usize,     // candidate set size per iter (default 1024)
    pub xi: f64,                 // EI exploration parameter (default 0.01)
    pub seed: u64,               // RNG seed for reproducibility
}

pub struct BoResult {
    pub x_best: DVector<f64>,
    pub y_best: f64,
    pub history: Vec<(DVector<f64>, f64)>,  // every evaluation in order
}

pub fn minimize<F>(
    objective: F,
    bounds: Vec<(f64, f64)>,      // (lo, hi) per dimension
    cfg: BoConfig,
) -> BoResult
where
    F: Fn(&DVector<f64>) -> f64;
```

Internally:
1. Sample `n_initial` Latin-hypercube points within bounds, evaluate objective.
2. Loop `n_iters` times:
   a. `gp = GaussianProcess::fit_ml(x_history, y_history, MlFitConfig::default())`
   b. Sample `n_candidates` uniform random points within bounds.
   c. For each candidate, compute `(mean, var) = gp.predict(...)`, then `ei = expected_improvement(mean, var.sqrt(), y_best, xi)`.
   d. Pick the candidate with max EI; evaluate objective; append.
3. Return best `(x, y)` seen + full history.

### Expected Improvement formula

For minimization with current best `f_best`, predictive mean `μ`, predictive stddev `σ`, exploration `ξ`:

```
improvement = f_best - μ - ξ
z           = improvement / σ
ei          = improvement * Φ(z) + σ * φ(z)    if σ > 0
            = max(improvement, 0)              if σ == 0
```

where `Φ` is the standard normal CDF and `φ` the PDF. Implement both inline (no `statrs` dep) using `libm::erf` (already a transitive dep via `nalgebra`).

## Definition of done

1. `crates/yee-surrogate/src/bo.rs` exists; public API as above.
2. `crates/yee-surrogate/src/lib.rs` re-exports `BoConfig`, `BoResult`, `minimize`.
3. `crates/yee-surrogate/tests/bo_synthetic.rs` ships and passes:
   - On `f(x) = (x - 3)² + sin(5x)` on `[0, 6]`, BO with budget=20 reaches `y_best < 0.0` (the function has a global minimum near `x ≈ 2.79, y ≈ -1.05`).
   - Random search (same budget) on the same function reaches `y_best > -0.5` on average — pure-random must be markedly worse.
4. README updated with a "Bayesian optimization" section.
5. Verification chain (`cargo build / clippy / test / fmt --check` on `-p yee-surrogate`) all exit 0.

## Lane

`crates/yee-surrogate/**` only. No new workspace deps. Use `rand` only if already in surrogate's Cargo.toml or workspace deps; otherwise use a simple linear-congruential generator inline (BO doesn't need cryptographic randomness).

## Escape hatch

If `fit_ml` underflows on the small initial training set (n=5), surface and revisit — the GP needs ≥ 3 distinct points for the Cholesky to be well-conditioned. Either expand `n_initial` default to 8 or fall back to `GaussianProcess::fit` with fixed hyperparameters for the first iteration.
