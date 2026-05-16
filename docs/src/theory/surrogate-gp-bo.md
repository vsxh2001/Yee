# Gaussian Process Surrogates and Bayesian Optimization — Theory of Operation

This page is the theory-of-operation reference for the Gaussian-process
(GP) surrogate and Bayesian-optimization (BO) loop shipped in the
`yee-surrogate` crate. Same audience as the planar-MoM and FDTD pages:
an engineer reading source with a textbook open. Equations are written
in KaTeX so the inline notation can stay close to the Rust source.

## 1. Introduction

Yee's solvers are deterministic but expensive. A single MoM matrix fill
for mom-001 burns several CPU-minutes; a Phase-3 microstrip filter sweep
with twenty design variables and a hundred candidate geometries would
take days of raw `mom::solve` time before any optimizer could think.
`ROADMAP.md` Phase 3 makes that affordable by routing the design loop
through a *surrogate* — a cheap statistical model fitted to a handful of
real solver evaluations — and then asking the surrogate, not the solver,
which candidate to try next.

The two building blocks for that loop are:

- **Phase 3.gp.0 / 3.gp.1** — a Gaussian-process regressor over the
  parameter-to-objective map, with hand-tuned and likelihood-fit
  hyperparameters respectively.
- **Phase 3.bo.0** — a single-objective Bayesian-optimization driver
  that wraps the GP with an Expected-Improvement acquisition.

The pair gives us a Phase-3 walking-skeleton optimizer: it never beats a
hand-tuned gradient method on smooth convex problems, but it shines on
the regime that actually matters here — expensive, noisy, derivative-
free objectives in modest (≲ 20) dimensions. That's exactly the EM-
sweep shape.

## 2. Gaussian process regression

A Gaussian process is a distribution over functions, fully specified by
a mean function $m(x)$ and a covariance kernel $k(x, x')$. For
regression we take $m \equiv 0$ (after centering $y$) and let the
kernel carry all structural assumptions:

$$
f \sim \mathcal{GP}\bigl(0,\; k(x, x')\bigr).
$$

Yee uses the *radial-basis-function* (RBF, a.k.a. squared-exponential)
kernel,

$$
k(x, x') = \sigma_f^2 \exp\!\left(-\frac{\lVert x - x' \rVert^2}{2 \ell^2}\right),
$$

with signal variance $\sigma_f^2$ and a single isotropic length scale
$\ell$. The RBF kernel encodes "the function is infinitely
differentiable and varies on a scale of $\ell$"; longer $\ell$ means
smoother, slower-varying functions. We additionally model observation
noise as i.i.d. Gaussian with variance $\sigma_n^2$, which adds
$\sigma_n^2 I$ to the diagonal of the training-data Gram matrix.

Given training inputs stacked into a matrix $X \in \mathbb{R}^{n \times d}$
and outputs $\mathbf{y} \in \mathbb{R}^n$, the Gram matrix is
$K_{ij} = k(x_i, x_j)$. The conditioning step that turns a prior into a
posterior solves a single dense linear system,

$$
\boldsymbol{\alpha} = \bigl(K + \sigma_n^2 I\bigr)^{-1} \mathbf{y},
$$

which `GaussianProcess::fit` performs by Cholesky-factoring
$K + \sigma_n^2 I$ once and caching both the factor and the solve
result. The Cholesky is the workhorse: it is numerically stable (the
matrix is symmetric positive-definite by construction whenever
$\sigma_n^2 > 0$ and the training inputs are distinct), costs $O(n^3)$
once, and amortizes to $O(n^2)$ per query.

At a test point $x_*$, write $\mathbf{k}_* = [k(x_*, x_1), \ldots, k(x_*, x_n)]^\top$.
The posterior mean and variance are

$$
\mu_* = \mathbf{k}_*^\top \boldsymbol{\alpha},
\qquad
\sigma_*^2 = k(x_*, x_*) - \mathbf{k}_*^\top K^{-1} \mathbf{k}_*.
$$

The mean is one dot product (cost $O(n \cdot d)$), so streaming queries
are cheap. The variance reuses the cached Cholesky factor with one
triangular solve, costing $O(n^2 \cdot d)$ per query — still trivial at
the $n \lesssim 50$ scale we target.

A typical `fit` call from a caller's perspective is:

```rust,ignore
use nalgebra::{DMatrix, DVector};
use yee_surrogate::GaussianProcess;

let x = DMatrix::from_column_slice(n, 1, &x_train);
let y = DVector::from_row_slice(&y_train);
let gp = GaussianProcess::fit(x, y, /*length_scale=*/ 0.5,
                              /*sigma_f=*/ 1.0, /*sigma_n=*/ 1e-4)?;
let (mean, var) = gp.predict(&DVector::from_row_slice(&[x_star]));
```

## 3. Log marginal likelihood

Hand-picking $(\ell, \sigma_f, \sigma_n)$ works for toys; for real
sweeps we want the data to choose. The natural objective is the
*log marginal likelihood* (LML), the log-probability of the observed
$\mathbf{y}$ under the GP prior after marginalising out $f$:

$$
\log p(\mathbf{y} \mid X, \theta)
= \underbrace{-\tfrac{1}{2} \mathbf{y}^\top K_\theta^{-1} \mathbf{y}}_{\text{data fit}}
\;\underbrace{- \tfrac{1}{2} \log\lvert K_\theta \rvert}_{\text{complexity penalty}}
\;\underbrace{- \tfrac{n}{2} \log 2\pi}_{\text{constant}}.
$$

The decomposition is the standard one (Rasmussen & Williams 2006, §5.4).
The data-fit term rewards kernels that explain $\mathbf{y}$ — small
quadratic form $\mathbf{y}^\top K^{-1} \mathbf{y}$ means the prior
covariance lines up with the observed correlations. The
log-determinant term penalises models with too much explanatory power:
a kernel that can fit anything has a large $\lvert K \rvert$, raising
the penalty. The remaining $-\tfrac{n}{2} \log 2\pi$ is independent of
$\theta$ and exists only to make the LML a proper log-density; it can
be dropped during optimization with no change to the argmax.

Two properties make the LML the right objective for hyperparameter
fitting, not just one among many:

- **Automatic Occam's razor.** The data-fit / complexity tradeoff is
  intrinsic to the marginalisation, not a separately tuned
  regularisation knob. A model with too short a length scale fits noise
  and pays a complexity bill; a model with too long a length scale
  smooths over signal and pays a data-fit bill. The LML balances them.
- **It's the type-II likelihood.** Maximizing $\log p(\mathbf{y} \mid X, \theta)$
  is empirical Bayes over the kernel hyperparameters, which converges
  to the true generative parameters under mild conditions as
  $n \to \infty$. No held-out validation set is needed.

`GaussianProcess::log_marginal_likelihood` computes this term-by-term
on the cached Cholesky factor: $\log\lvert K \rvert = 2 \sum_i \log L_{ii}$
where $L$ is the lower-triangular Cholesky factor, and the data-fit
term is one inner product against the cached
$\boldsymbol{\alpha} = K^{-1} \mathbf{y}$.

## 4. Hyperparameter optimization

`GaussianProcess::fit_ml` maximizes the LML over
$\theta = (\ell, \sigma_f, \sigma_n)$ by gradient ascent. Two
implementation choices are worth calling out.

**Log-space parameterization.** All three hyperparameters must stay
strictly positive (a negative length scale is meaningless, a
non-positive noise breaks the Cholesky). We optimize over
$(\log \ell,\, \log \sigma_f,\, \log \sigma_n)$ so positivity is built
into the parameterization and no constraint or projection is needed.
This is the standard trick — Rasmussen & Williams §5.4 use the same
substitution.

**Numerical gradient.** The analytic LML gradient is

$$
\frac{\partial}{\partial \theta_k} \log p(\mathbf{y} \mid X, \theta)
= \tfrac{1}{2} \operatorname{tr}\!\left[
\bigl(\boldsymbol{\alpha}\boldsymbol{\alpha}^\top - K^{-1}\bigr)
\frac{\partial K}{\partial \theta_k}
\right],
$$

which is $O(n^3)$ per parameter (the $K^{-1} \partial K/\partial\theta_k$
trace) and easy to get wrong on the kernel-derivative side. For the
small-$n$ ($\lesssim 50$) regime this surrogate targets, central
differences in log-space at step $10^{-3}$ cost six fresh K-builds per
iteration — six $O(n^3)$ Choleskys instead of one — but eliminate an
entire class of analytic-derivative bugs. The tradeoff would invert at
$n \gtrsim 100$ or for many-parameter kernels (ARD with one length
scale per input dimension): there the analytic gradient is cheaper *per
step*, and the derivative bookkeeping pays for itself.

The optimizer is plain gradient ascent with one safety knob: when the
log-space gradient norm exceeds 1, the update is rescaled to unit
length and multiplied by `MlFitConfig::gradient_step` (default `0.05`).
The very large gradients the LML produces near poorly-scaled starting
points would otherwise overshoot $\sigma_n$ into floating-point
underflow before the second iteration. Near the optimum the gradient is
small and the raw step is used, preserving curvature-sensitive
convergence.

```rust,ignore
use nalgebra::{DMatrix, DVector};
use yee_surrogate::{GaussianProcess, MlFitConfig};

let cfg = MlFitConfig {
    initial_length_scale: 1.0,
    initial_sigma_f: 1.0,
    initial_sigma_n: 1e-3,
    ..Default::default()
};
let gp = GaussianProcess::fit_ml(x, y, cfg)?;
let (mean, var) = gp.predict(&DVector::from_row_slice(&[x_star]));
```

