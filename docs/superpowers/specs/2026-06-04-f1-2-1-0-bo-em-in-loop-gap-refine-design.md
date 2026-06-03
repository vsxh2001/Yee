# F1.2.1.0 — single-gap surrogate-BO EM-in-loop dimensional refinement — design

**Date:** 2026-06-04
**ADR:** [ADR-0157](../../src/decisions/0157-f1-2-1-0-bo-em-in-loop-gap-refine.md)

## Problem

The analytic dimensioner (`dimension_edge_coupled`) picks gaps by bisecting the impedance-k
(`coupling_coefficient`); ADR-0155/K2 showed the EM-realized resonant-k diverges from it by ~17-26 %.
F1.2.1.0 closes that gap for ONE inter-resonator gap via EM-in-the-loop Bayesian optimization — the
design engine's walking skeleton.

## Goal

A gated `bo-coupling-001`: starting from the off-target analytic-seed gap, a 1-D EI/GP BO loop with
`yee_fem::coupled_resonator_k` in the objective drives the EM-measured resonant-k to the synthesis
`target_k`, converging measurably tighter than the seed.

## Non-goals (deferred to later F1.2.1 bricks)

Multi-gap refinement, multi-D BO, surrogate reuse across gaps, Qe/feed dimensioning (ADR-0156),
`yee-server` API. ONE gap only.

## Architecture

- **Home:** `crates/yee-validation/tests/bo_coupling_001.rs` (new). `yee-validation` already deps
  `yee-filter` + `yee-fem`. Do NOT add `yee-fem` to `yee-filter` (WASM-safety). A thin private
  helper in the test (or a small non-WASM module under `yee-validation/src`) holds the refine loop.
- **Seed:** `dimension_edge_coupled(project, substrate)` → `EdgeCoupledDimensions { gaps_m, target_k }`
  for a small fixture (a 0.5 dB Chebyshev edge-coupled filter; pick one inter-resonator gap `i`).
- **Objective** (passed to `minimize`): given a normalized `x ∈ [0,1]`, unscale to a gap
  `g = g_lo + x·(g_hi − g_lo)` within `[GAP_MIN_M, GAP_MAX_M]` bracketing the seed; build the
  coupled-pair geometry at `g` (reuse `coupled_resonator_k`'s `CoupledResonatorGeom`, with the
  fixture's W/h/ε_r and `gap_s = g`); run `coupled_resonator_k(geom, n_pts)` → `k_fem`; return
  `|k_fem − target_k[i]|`.
- **Optimizer:** `yee_surrogate::minimize(objective, vec![(0.0, 1.0)], BoConfig { n_initial: 3,
  n_iters: 9, .. })`. Unscale `result.x_best` → refined gap. (Normalize to [0,1] — mandatory for the
  GP length-scale; see ADR Risk 1.)
- **Bracket:** `[g_lo, g_hi]` around the seed wide enough to contain the gap where `k_fem = target_k`
  (e.g. `[seed·0.5, seed·1.5]` clamped to `[GAP_MIN_M, GAP_MAX_M]`); pick the fixture so `target_k` is
  reachable inside it (the seed off-target but the target attainable by moving the gap).

## Data flow

`dimension_edge_coupled` → seed gap + `target_k` → `minimize` over normalized gap, each eval =
`coupled_resonator_k` FEM sweep → `k_fem` → `|k_fem − target_k|` → `BoResult.x_best` → unscale →
refined gap → measure `k_fem(refined)` → grade.

## Gate `bo-coupling-001` (`#[ignore]`'d, `--release`, boxed)

1. `|k_fem(seed) − target_k| / target_k ≥ 0.10` (seed genuinely off — not vacuous).
2. `|k_fem(refined) − target_k| < |k_fem(seed) − target_k|` (BO strictly improved).
3. `|k_fem(refined) − target_k| / target_k < 0.08` (converged, tighter than the seed).
Print the seed gap / k_fem(seed) / target_k / refined gap / k_fem(refined) / the BO history (the
`BoResult.history` of (gap, |Δk|)) — auditable.

## Testing

- A fast unit test (non-ignored, NO FEM solve): the objective's gap-normalization round-trips
  ([0,1]↔metres), the fixture's seed + target_k are finite/positive, the bracket contains the seed.
- The heavy gate `#[ignore]`'d + a new `bo-coupling-001` `--release` CI step (mirror `fem-coupling`).
- Heavy run = ~12 sequential FEM evals ≈ ~57 min, boxed (`scripts/yee-box.sh`). **The orchestrator
  runs it** (the agent writes code only — the misfire-lesson split).

## Risks

- GP unit-scale → [0,1] normalization (mandatory).
- `target_k` reachability → fixture choice; if BO can't reach `target_k` in-bracket, surface as an
  honest finding (do NOT widen the tolerance).
- noiseless objective (FEM deterministic) → GP `sigma_n` small; monotone 1-D → easy for BO.
