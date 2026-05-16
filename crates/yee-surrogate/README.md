# yee-surrogate

> ML surrogate models for parameterized EM simulation outputs. **Phase 3 deliverable.**

This crate provides a `Surrogate` trait abstraction over fast approximators for
parameter -> S-parameter (or other simulator output) maps. The Phase 3 walking
skeleton ships a trivial nearest-neighbor baseline so the dataset / training /
prediction plumbing exists end-to-end before any heavy ML dependency is pulled
into the workspace.

## Why a surrogate layer?

A single full-wave run (MoM or FDTD) is expensive. Optimization, tolerance
analysis, and interactive design exploration all repeat-evaluate the same
geometry family with shifted parameters. A surrogate trained on a modest sweep
of full-solver runs lets the GUI scrub design parameters and watch S11 / Smith
update at interactive rates while the high-fidelity solver runs in the
background to backfill the dataset.

## Scope (Phase 3 walking skeleton)

- `Sample` — `(params: Vec<f64>, output: Vec<Complex64>)` pair.
- `Dataset` — append-only collection of samples.
- `Surrogate` trait — `train(&Dataset)` + `predict(&[f64])`.

## Backends

| Backend           | Output | Uncertainty | When to use                                                              |
|-------------------|--------|-------------|--------------------------------------------------------------------------|
| `NearestNeighbor` | any    | none        | Sanity check, undertrained-fallback, or coarse last-resort interpolation |
| `GaussianProcess` | scalar | posterior variance via RBF kernel | Low-dimensional designs (≲ 20 params), small N, calibrated uncertainty for active learning |

### `NearestNeighbor`

Returns the output of the closest training sample by Euclidean L2 distance in
parameter space. Always defined, no hyperparameters, no failure modes other
than empty/inconsistent input. Useful as a sanity check and as a fallback when
the model is undertrained.

### `GaussianProcess`

Squared-exponential (RBF) kernel GP regressor, scalar output:

```text
k(x, x') = sigma_f² · exp(-‖x - x'‖² / (2 · length_scale²))
```

Training Cholesky-factors `K + sigma_n² I` once and caches both `α = K⁻¹ y`
(for fast mean queries) and the Cholesky factor (for variance queries via a
single triangular solve per query). Per-query cost: O(n·d) mean, O(n²·d)
variance.

Surrogate-trait coverage is provided by `GpSurrogate`, which adapts the
scalar GP to the existing `Dataset` shape by treating the real part of
`sample.output[0]` as the regression target. For multi-output use cases,
call `GaussianProcess::fit` directly per output channel.

Usage:

```rust,ignore
use nalgebra::{DMatrix, DVector};
use yee_surrogate::GaussianProcess;

let x = DMatrix::from_column_slice(n, 1, &x_train);
let y = DVector::from_row_slice(&y_train);
let gp = GaussianProcess::fit(x, y, /*length_scale=*/ 0.5, /*sigma_f=*/ 1.0, /*sigma_n=*/ 1e-4)?;
let (mean, var) = gp.predict(&DVector::from_row_slice(&[x_star]));
```

## Hyperparameter optimization

`GaussianProcess::fit_ml` maximizes the log marginal likelihood

```text
log p(y | X, θ) = -0.5 · yᵀ K⁻¹ y - 0.5 · log|K| - (n/2) · log(2π)
```

over `θ = (length_scale, sigma_f, sigma_n)` via gradient ascent in log-space
(so each hyperparameter stays strictly positive by construction). Gradients
are computed by **central differences** on `log_marginal_likelihood`, not
analytically: the analytic gradient requires `tr(K⁻¹ ∂K/∂θ)`, which is
O(n³) per parameter and easy to get wrong. Central differences cost 6
K-builds per iteration, which is the simpler tradeoff for the small-n (≲ 50)
problems this surrogate targets. The optimizer caps per-iteration step
magnitude in log-space so the very large gradients near a poorly-scaled
starting point can't overshoot into underflow.

Usage:

```rust,ignore
use nalgebra::{DMatrix, DVector};
use yee_surrogate::{GaussianProcess, MlFitConfig};

let x = DMatrix::from_column_slice(n, 1, &x_train);
let y = DVector::from_row_slice(&y_train);
let cfg = MlFitConfig {
    initial_length_scale: 1.0,
    initial_sigma_f: 1.0,
    initial_sigma_n: 1e-3,
    ..Default::default()
};
let gp = GaussianProcess::fit_ml(x, y, cfg)?;
println!("optimized log marginal likelihood = {}", gp.log_marginal_likelihood());
let (mean, var) = gp.predict(&DVector::from_row_slice(&[x_star]));
```

The returned `GaussianProcess` is a fresh refit with the optimized
hyperparameters, so its cached `α` and Cholesky factor are consistent with
the returned `(length_scale, sigma_f, sigma_n)` accessors.

## Bayesian optimization

`bo::minimize` runs a single-objective Bayesian-optimization loop on top of
`GaussianProcess`. Each iteration refits a GP via `fit_ml` on the running
evaluation history, then maximizes Expected Improvement over a uniform
random candidate set to pick the next point. The implementation is
~150 LOC plus tests and adds no new crate dependencies — the standard
normal CDF/PDF used by EI are inlined via the Abramowitz & Stegun 7.1.26
rational approximation to `erf`, and randomness comes from a small inline
`xorshift64` seeded by `BoConfig::seed` for reproducibility.

### Acquisition: Expected Improvement (minimization)

For current best `f_best`, predictive mean `μ`, stddev `σ`, exploration `ξ`:

```text
improvement = f_best - μ - ξ
z           = improvement / σ                        if σ > 0
ei          = improvement · Φ(z) + σ · φ(z)          if σ > 0
            = max(improvement, 0)                    if σ == 0
```

The exploration parameter `xi` (default `0.01`) biases toward higher-variance
candidates; raise it if BO converges too eagerly to a local minimum.

### Initial design

`n_initial` Latin-hypercube points (default 5). Each dimension is split into
`n_initial` equal strata, one stratified value is drawn per stratum, then
strata are permuted independently across dimensions. This produces a
space-filling initial design without the clustering that pure-random
sampling can exhibit at small `n`.

### Usage

```rust,ignore
use nalgebra::DVector;
use yee_surrogate::{minimize, BoConfig};

let objective = |x: &DVector<f64>| (x[0] - 3.0).powi(2) + (5.0 * x[0]).sin();
let bounds = vec![(0.0, 6.0)];
let cfg = BoConfig {
    n_initial: 5,
    n_iters: 20,
    n_candidates: 1024,
    xi: 0.01,
    seed: 0xC0FFEE,
};
let res = minimize(objective, bounds, cfg);
println!("best x = {:?}, y = {}", res.x_best, res.y_best);
// res.history holds every (x, y) evaluation in order.
```

### Validation

`tests/bo_synthetic.rs` runs the deceptive 1-D objective
`f(x) = (x - 3)² + sin(5x)` on `[0, 6]` with budget 5 + 20 and asserts
BO `y_best < 0.0` (a fine sweep places the global minimum at `x ≈ 3.422`,
`y ≈ -0.8077`). A 25-call pure-random baseline run from the same seed set
is required to lose head-to-head against BO on the best `y` across seeds.

### Scope and out-of-scope

In:
- Single-objective minimization with continuous bounded parameters.
- EI acquisition on a uniform random candidate set.

Out (Phase 3.bo.1+):
- L-BFGS / multi-start gradient acquisition optimization.
- Constrained optimization.
- Batch BO.

## Active learning

`al::active_learn` (re-exported as `active_learn`) runs a variance-acquisition
active-learning loop on top of `GaussianProcess`. It is the same iteration
loop as `bo::minimize` with Expected Improvement swapped for predictive
variance: instead of picking the next point that most likely beats the
current best, AL picks the point of maximum posterior uncertainty so that
querying it most reduces the GP's overall variance (MacKay 1992).

Use this when the goal is **an accurate surrogate over the whole parameter
space**, not finding a minimum. The two acquisitions are dual: BO concentrates
on low-objective regions, AL concentrates on high-uncertainty regions. For a
fixed simulator budget, AL is what you want before BO when you have no
specific objective in mind yet — e.g. building a 25-sample surrogate of an
S-parameter response so the GUI can scrub geometry parameters interactively.

### Acquisition: predictive variance

```text
score(x) = var(x)   where (_, var) = GaussianProcess::predict(x)
```

No exploration parameter is needed; the GP's predictive variance is already
the right uncertainty estimate.

### Usage

```rust,ignore
use nalgebra::DVector;
use yee_surrogate::{active_learn, AlConfig};

let objective = |x: &DVector<f64>| x[0].sin();
let bounds = vec![(0.0, std::f64::consts::TAU)];
let cfg = AlConfig {
    n_initial: 5,
    n_iters: 20,
    n_candidates: 1024,
    seed: 0xC0FFEE,
};
let res = active_learn(objective, bounds, cfg);
// res.history holds every (x, y) evaluation in order.
// res.final_gp is a GP refit on the full history — call predict on it directly.
let (mean, var) = res.final_gp.predict(&DVector::from_row_slice(&[1.0]));
```

### Validation

`tests/al_synthetic.rs` runs the smooth 1-D objective `sin(x)` on `[0, 2π]`
with budget 5 + 20 and asserts that the AL-trained GP's test-grid RMSE is
strictly less than half of a GP trained on 25 uniform-random points from
the same seed family. Observed ratios are typically ~0.01 (AL is two
orders of magnitude more accurate than random on this problem); the 0.5
threshold leaves a wide margin so the test is not seed-fragile.

### Scope and out-of-scope

In:
- Single-output variance-acquisition AL with continuous bounded parameters.
- Uniform random candidate set per iteration; pick by max variance.
- Returns the full history plus a final GP refit so callers can call
  `predict` / `predict_mean` directly.

Out (Phase 3.al.1+):
- Multi-acquisition AL (variance + objective trade-off).
- Multi-fidelity / cost-aware AL (Lam 2015).
- Look-ahead / non-myopic AL (Gonzalez 2016).
- Batch / parallel AL queries.

## Multi-objective optimization (NSGA-II)

`nsga2::minimize` (re-exported as `nsga2_minimize`) runs a textbook
NSGA-II loop (Deb et al. 2002) on top of a Latin-hypercube initial design,
binary tournament selection over `(rank, crowding distance)`, simulated
binary crossover (SBX), and polynomial mutation. The return contains the
final population, its objective vectors, and the indices that lie on the
non-dominated front.

```rust,ignore
use nalgebra::DVector;
use yee_surrogate::{Nsga2Config, nsga2_minimize};

let d = 30;
let zdt1 = |x: &DVector<f64>| -> Vec<f64> {
    let f1 = x[0];
    let g = 1.0 + 9.0 / (d as f64 - 1.0) * x.iter().skip(1).sum::<f64>();
    vec![f1, g * (1.0 - (f1 / g).sqrt())]
};
let result = nsga2_minimize(zdt1, vec![(0.0, 1.0); d], 2, Nsga2Config::default());
```

`tests/nsga2_zdt1.rs` validates the implementation on the canonical ZDT1
30-variable two-objective benchmark and asserts the Inverted Generational
Distance against the analytic Pareto front `f2 = 1 - sqrt(f1)` is below
`0.05`; the canonical NSGA-II paper reports IGD ≈ 0.005 at 250 generations.

## Future direction (Phase 3.1+)

- Anisotropic RBF / Matérn kernels for the GP backend.
- Analytic-gradient + L-BFGS hyperparameter optimization for n ≫ 50 where
  the numerical-gradient cost dominates.
- MLP / residual-MLP backend for medium-dimensional spaces and amortized
  inference.
- Fourier neural operator (FNO) / DeepONet for field-level outputs, not just
  scalar S-parameters.
- On-disk dataset format (Arrow / Parquet) shared with `yee-io`.
- Active-learning loop driving the full-wave solver: pick the next sample by
  maximizing predicted uncertainty * cost-adjusted information gain.

All future backends sit behind the same `Surrogate` trait and the same
`Dataset` storage so the GUI and CLI never have to know which model is loaded.
