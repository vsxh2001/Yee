# ADR-0009: Hand-rolled RBF-kernel Gaussian process in `yee-surrogate`, no external GP crate

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

Phase 3.gp.0 introduced a Gaussian-process (GP) surrogate as the
first non-trivial backend for the `Surrogate` trait in
`yee-surrogate`. Phase 3.gp.1 added marginal-likelihood
hyperparameter optimisation (`fit_ml`). The GP is the workhorse
behind Phase 3.bo.0's Bayesian-optimisation loop (see ADR-0010), and
its design choice is the single largest dependency call in the
`yee-surrogate` crate.

The candidate implementations in May 2026 fall in three buckets:

1. **Use a Rust GP crate.** The ecosystem options are:
   - **`friedrich`** (BSD-3, last published 0.5.0 in 2022). RBF /
     Matérn kernels, dense Cholesky, no hyperparameter optimisation
     out of the box, light maintenance.
   - **`linfa-gaussian-process`** (under the `linfa` ML toolkit).
     Promising but only the `linfa-clustering` GP-mixture-model
     side is mature; the regression GP is partial as of `linfa
     0.7`.
   - **`egobox-gp`** (Apache-2.0, derived from ONERA's `egobox`).
     Full-featured but a heavyweight dependency that drags in
     `egobox-doe`, `egobox-moe`, and an `ndarray-linalg` /
     OpenBLAS dependency to link against system BLAS.
2. **Push the GP to Python.** Expose only the data-collection side
   from Rust and let users plug in `scikit-learn`, `GPy`, or
   `GPyTorch` on the Python side. Familiar to the ML audience;
   imposes a runtime Python dependency on what is supposed to be a
   pure-Rust simulator.
3. **Hand-roll the GP in Rust.** Implement the kernel evaluation,
   Cholesky factorisation, predict-mean / predict-variance, and
   log-marginal-likelihood gradient ascent directly in
   `yee-surrogate`, on top of `nalgebra`'s already-present `DMatrix`
   and `Cholesky` types.

The constraints that narrow the field:

- **Yee is a pure-Rust simulator.** Option 2 introduces a runtime
  Python dependency on the surrogate path, which inverts the
  layering: `yee-py` is supposed to wrap `yee-surrogate`, not the
  other way around. Rejected.
- **The `Surrogate` trait is small** — `fit(&mut, &[Sample])`,
  `predict(&self, x) -> (mean, variance)`. The GP backend needs to
  fit cleanly behind this trait, with no escape hatch into kernel
  internals or hyperprior configuration leaking out. A heavyweight
  crate whose public API is broader than the trait is harder to wrap
  cleanly.
- **`nalgebra` is already a transitive dependency** of `yee-mom`
  (small-matrix utilities) and `yee-core` (units). Adding a GP that
  uses `nalgebra::Cholesky` does not enlarge the workspace's
  build-time footprint.
- **`ndarray-linalg`** (the LA backend that `egobox-gp` and
  `linfa-gaussian-process` ultimately reach for) requires a system
  BLAS, which is a deployment concern for the wheel-build CI and a
  source of friction on Windows. `nalgebra`'s pure-Rust Cholesky
  has no system dependency.
- **The expected operating regime is small N.** The GP is used as a
  surrogate over expensive solver evaluations: a typical Bayesian-
  optimisation loop runs 20-50 evaluations of a few-minute MoM
  solve. Cholesky factorisation at `n ≤ 50` is ~125 µs on a single
  CPU core; the constant-factor advantage of a BLAS-backed crate is
  irrelevant at this scale.

Within option 3, two sub-questions:

- **Kernel choice.** The RBF (squared-exponential) kernel is the
  textbook default and is what Phase 3.bo.0's Expected-Improvement
  acquisition function assumes implicitly (it assumes the posterior
  variance is smooth). Matérn-5/2 is the alternative most often
  cited. The decision is to ship RBF in Phase 3.gp.0 and leave
  Matérn behind a future configuration knob; the validation case
  (`gp_baseline` on sin(x)) is the same shape for both kernels and
  not a discriminator.
- **Hyperparameter optimisation.** Phase 3.gp.0 ships hand-tuned
  hyperparameters (length scale, signal variance, noise variance)
  embedded in `GpConfig`. Phase 3.gp.1 adds a `fit_ml` method that
  optimises the log marginal likelihood by numerical-gradient
  ascent in log-space (`log_lengthscale`, `log_sigma_f`,
  `log_sigma_n`). The numerical-gradient choice trades CPU for
  implementation simplicity: an analytic gradient is
  `O(n³)` per parameter (one matrix solve per partial), the
  numerical version is `2 · O(n³)` per parameter (two
  forward log-likelihood evaluations per partial via central
  differences). At `n ≤ 50` and 3 hyperparameters the whole
  optimisation runs in milliseconds, so the analytic-gradient
  speedup buys nothing observable. At `n > 200` the picture flips.

## Decision

Implement the GP surrogate in-tree in `yee-surrogate`:

```rust
pub struct GaussianProcess {
    samples: Vec<Sample>,
    lengthscale: f64,
    sigma_f: f64,      // signal std-dev
    sigma_n: f64,      // observation noise std-dev
    factor: Option<Cholesky<f64, Dyn>>,  // cached K + sigma_n^2 I
    alpha: Option<DVector<f64>>,         // cached K^{-1} y
}

impl GaussianProcess {
    pub fn fit(&mut self, samples: &[Sample]) -> Result<()>;
    pub fn predict_mean(&self, x: &[f64]) -> f64;
    pub fn predict(&self, x: &[f64]) -> (f64, f64);  // (mean, variance)
    pub fn log_marginal_likelihood(&self) -> f64;
    pub fn fit_ml(&mut self, samples: &[Sample], cfg: MlFitConfig) -> Result<()>;
}
```

Kernel: **RBF / squared-exponential** with isotropic length scale.
Linear algebra: `nalgebra::Cholesky<f64, Dyn>` on the
`K + σ_n² I` Gram matrix.

`fit_ml` optimises three hyperparameters
(`log_lengthscale`, `log_sigma_f`, `log_sigma_n`) by
gradient-ascent on the log marginal likelihood, using **central-
difference numerical gradients**:

```text
∂L/∂θ_i ≈ ( L(θ + ε e_i) − L(θ − ε e_i) ) / (2ε)
```

with a small `ε` (1e-4 in log-space), Armijo-style backtracking line
search, and a configurable iteration cap from `MlFitConfig`.

No external GP crate is taken. No system BLAS / LAPACK dependency
is introduced. The implementation lives entirely under
`crates/yee-surrogate/src/gp/`.

## Consequences

**What becomes easier:**

- Pure-Rust, no system BLAS, no Python dependency, no transitive
  dependency on a heavyweight ML toolkit. The wheel-build CI
  (`publish-wheels.yml` under `manylinux_2_28`) does not need any
  additional native packaging.
- The `Surrogate` trait surface stays small. The GP fits behind
  `fit` / `predict` exactly the same way the
  `NearestNeighborBaseline` does, and the trait does not leak GP-
  specific concepts (kernel choice, hyperpriors) to its consumers.
- Bayesian optimisation (ADR-0010) consumes the GP through its
  `predict(x) -> (mean, variance)` method without coupling to any
  particular implementation, which keeps the door open to swapping
  backends later.
- Numerical-gradient `fit_ml` is short, auditable, and easy to
  validate against a textbook log-likelihood for sin(x) — the
  `gp_ml` validation test does exactly this.

**What becomes harder:**

- Numerical gradient is `O(n³)` per hyperparameter at the cost of two
  Cholesky factorisations per finite-difference step. At `n ≤ 50`
  (the regime the trait is used in) this is microseconds; at
  `n > 200` it becomes noticeable, and at `n > 500` it becomes the
  dominant cost of a BO iteration. Swapping to an analytic gradient
  is a future bounded change (touch `fit_ml` and add a
  `dlog_likelihood/dθ` helper) but it is real work.
- No advanced GP features ship in Phase 3.gp.0/1:
  - No multi-output GP (one independent GP per output is the user-
    visible workaround).
  - No structured kernels (Matérn, ARD with per-dimension length
    scales, periodic).
  - No exact MCMC for hyperparameters; the marginal-likelihood
    optimum is a point estimate, not a posterior.
  Each of these is recoverable by feature-gating an alternative
  backend (option 3 → option 1 swap), but they are deferred.
- Maintenance burden lives in the workspace. A bug in the Cholesky
  cache or the predict-variance formula is ours to fix; we cannot
  upstream it.

**What's now closed off:**

- Taking a runtime Python dependency on `scikit-learn` or `GPyTorch`
  for the surrogate path. The Rust trait is the source of truth.
- Linking against a system BLAS (`OpenBLAS`, `Accelerate`, MKL) for
  the surrogate path. The Cholesky is `nalgebra`'s pure-Rust
  implementation.
- Adding `egobox-gp` or `linfa-gaussian-process` as a direct
  dependency in Phase 3.gp.x. A future ADR can revisit this when
  `n > 200` cases land.

## References

- `crates/yee-surrogate/src/gp/` — the implementation.
- `crates/yee-surrogate/src/gp/kernel.rs` — RBF kernel.
- `crates/yee-surrogate/src/gp/fit_ml.rs` — central-difference
  gradient-ascent log-marginal-likelihood optimiser.
- `crates/yee-surrogate/validation/gp_baseline.rs` — sin(x)
  validation; GP RMS beats nearest-neighbour baseline.
- `crates/yee-surrogate/validation/gp_ml.rs` — `fit_ml` improves on
  deliberately-bad hand-tuned hyperparameters.
- C. E. Rasmussen and C. K. I. Williams, *Gaussian Processes for
  Machine Learning*, MIT Press, 2006. Algorithm 2.1 (predictive
  mean / variance via Cholesky) and §5.4.1 (log marginal
  likelihood).
- `nalgebra` crate documentation, `Cholesky<f64, Dyn>`:
  <https://docs.rs/nalgebra/>.
- ADR-0010 — Bayesian optimisation consumes this GP through the
  `Surrogate` trait.
