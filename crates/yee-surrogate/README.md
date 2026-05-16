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
