# Phase 3.al.0 — Active learning for GP surrogates

**Status:** Draft  
**Owner:** TBD  
**Phase:** 3.al.0  
**Depends on:** Phase 3.gp.0 (shipped), 3.gp.1 (shipped), 3.bo.0 (shipped)  
**Blocks:** Phase 3.al.1 (multi-acquisition active learning), Phase 3.nl.0 (NL design surface — uses AL internally to bootstrap)

## Assumption being challenged

ROADMAP Phase 3 frames active learning as "Solver picks the next simulation points to maximize surrogate accuracy." This implicitly assumes a single acquisition strategy.

The minimal AL loop is just BO with a different acquisition: instead of Expected Improvement (which favors low-objective regions), use **predictive variance** (which favors high-uncertainty regions). Same GP infrastructure, same iteration loop, different acquisition.

The Phase 3.bo.0 module already has 90% of the machinery. Adding AL is ~80 lines plus a test.

## Scope

In:
- New module `crates/yee-surrogate/src/al.rs` (or reuse `bo.rs` with a flag)
- Public API: `active_learn` — same shape as `bo::minimize` but with variance acquisition
- Validation: 1-D `sin(x)` on `[0, 2π]`. Run AL with 5 initial + 20 iters. After 25 total samples, GP fit on the AL-chosen set has lower RMSE on a uniform 100-point test grid than GP fit on 25 random samples.

Out:
- Multi-fidelity AL (Lam 2015)
- Look-ahead / non-myopic AL (Gonzalez 2016)
- Cost-aware AL with non-uniform query cost
- Batch AL (parallel queries per iter)

## Public API

```rust
//! Active learning: pick next samples by maximum predictive variance.

pub struct AlConfig {
    pub n_initial: usize,        // (default 5)
    pub n_iters: usize,          // (default 20)
    pub n_candidates: usize,     // (default 1024)
    pub seed: u64,
}

pub struct AlResult {
    pub history: Vec<(DVector<f64>, f64)>,   // every (x, y) observation in order
    pub final_gp: GaussianProcess,            // GP fit on the full history
}

/// Run active learning: starting from a Latin-hypercube initial design,
/// iteratively pick the point of maximum predictive variance and query the
/// black-box at it.
pub fn active_learn<F>(
    objective: F,
    bounds: Vec<(f64, f64)>,
    cfg: AlConfig,
) -> AlResult
where
    F: Fn(&DVector<f64>) -> f64;
```

Implementation is BO with EI swapped for variance acquisition:

```rust
fn variance_acquisition(_mean: f64, var: f64) -> f64 { var }
```

Pick the candidate with **maximum** variance (vs. EI which picks max improvement). The rest of the loop is identical to `bo::minimize`.

`final_gp` is the GP trained on the full observation history; users can call `predict` / `predict_mean` directly on it.

## Definition of done

1. `crates/yee-surrogate/src/al.rs` exists; `AlConfig`, `AlResult`, `active_learn` public.
2. `lib.rs` re-exports.
3. `crates/yee-surrogate/tests/al_synthetic.rs`:
   - Objective: `sin(x)` on `[0, 2π]`.
   - Seed AL config: `n_initial=5, n_iters=20`.
   - Compare AL-trained GP RMSE vs. random-trained GP RMSE on a 100-point test grid.
   - Assert `al_rmse < 0.5 * random_rmse` (AL should be at least 2× better than random on this smooth-but-non-trivial 1-D problem).
4. README "Active learning" section added.
5. Full verification chain (`cargo build / clippy / test --release / fmt --check` on `-p yee-surrogate`) green.

## Lane

`crates/yee-surrogate/**` only. Reuse `bo`'s xorshift RNG + Latin-hypercube helpers — refactor into a tiny private `_internal` module if needed, but stay in-crate.

## Why this matters

This closes ROADMAP Phase 3's "Active learning loops" deliverable for v1.0. Once a user has a single full-fidelity simulator (an FdtdDriver-driven dipole run, say, that takes 10s), AL lets them produce a usable surrogate over the parameter space in 25 simulations instead of 500.

The same module is the foundation for Phase 3.nl.0: the NL design surface will, in its inner loop, run AL to bootstrap a surrogate before BO refines.

## Escape hatch

If the variance-acquisition loop converges to a single corner (degenerate "fence-post" behavior), surface and revisit. The standard fix is to combine variance with distance-from-existing-samples — Latin-hypercube the initial set wider, or add a small `xi`-equivalent floor that discourages re-sampling near existing points. Defer to Phase 3.al.1.

## References

- Settles, "Active Learning Literature Survey", 2010
- MacKay, "Information-Based Objective Functions for Active Data Selection", Neural Computation 1992 (variance acquisition's original derivation)
- Cohn / Ghahramani / Jordan, "Active Learning with Statistical Models", JAIR 1996
