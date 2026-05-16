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
- `NearestNeighbor` — returns the output of the closest training sample by
  Euclidean L2 distance in parameter space. Useful as a sanity check and as a
  fallback when the model is undertrained.

## Future direction (Phase 3.1+)

- Gaussian-process regression with anisotropic RBF / Matérn kernels for
  low-dimensional design spaces (<= ~20 parameters) with calibrated
  uncertainty for active-learning sample acquisition.
- MLP / residual-MLP backend for medium-dimensional spaces and amortized
  inference.
- Fourier neural operator (FNO) / DeepONet for field-level outputs, not just
  scalar S-parameters.
- On-disk dataset format (Arrow / Parquet) shared with `yee-io`.
- Active-learning loop driving the full-wave solver: pick the next sample by
  maximizing predicted uncertainty * cost-adjusted information gain.

All future backends sit behind the same `Surrogate` trait and the same
`Dataset` storage so the GUI and CLI never have to know which model is loaded.
