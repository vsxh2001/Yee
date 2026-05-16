//! bo_synthetic — Bayesian optimization beats pure-random search on a
//! deceptive 1-D objective.
//!
//! Objective: `f(x) = (x - 3)^2 + sin(5x)` on `x ∈ [0, 6]`.
//!
//! The quadratic `(x - 3)^2` has its bowl at `x = 3`, but the `sin(5x)`
//! ripple pulls the true minimum to `x ≈ 2.79` with `f ≈ -1.05`. Pure
//! random search wastes evaluations far from the bowl; BO with an EI
//! acquisition concentrates evaluations in the promising region after a
//! few iterations.
//!
//! Spec budget: 20 BO evaluations after 5 LHS initial samples = 25 total
//! objective calls. Random search uses the same 25-call budget for an
//! apples-to-apples comparison.
//!
//! Note on the global minimum: a 1-D fine sweep of `f` on `[0, 6]` puts the
//! global minimum at `x ≈ 3.422`, `y ≈ -0.8077`. The spec text mentions
//! `y ≈ -1.05` as the minimum, but that does not occur on this interval —
//! `f(2.79) ≈ +1.03`. The spec's `y_best < 0.0` threshold remains the right
//! signal that BO actually reached the valley and is not just exploring.
//!
//! Assertions:
//!   1. BO `y_best < 0.0` (the deep valley is reachable in 20 iters).
//!   2. BO beats pure-random search head-to-head: across the seed set, BO's
//!      best is strictly less than random search's best.

use nalgebra::DVector;
use yee_surrogate::{BoConfig, minimize};

/// The deceptive 1-D objective: bowl + ripple.
fn objective(x: &DVector<f64>) -> f64 {
    let v = x[0];
    (v - 3.0).powi(2) + (5.0 * v).sin()
}

/// Pure random search baseline. Uses the same xorshift seed → reproducible.
fn random_search(seed: u64, budget: usize) -> f64 {
    let mut state = if seed == 0 {
        0x9E37_79B9_7F4A_7C15
    } else {
        seed
    };
    let mut best = f64::INFINITY;
    for _ in 0..budget {
        // Inline xorshift64 — same recurrence the BO module uses.
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let u = ((state >> 11) as f64) * (1.0 / ((1_u64 << 53) as f64));
        let x = 6.0 * u; // map [0,1) → [0,6)
        let y = objective(&DVector::from_row_slice(&[x]));
        if y < best {
            best = y;
        }
    }
    best
}

#[test]
fn bo_beats_random_on_deceptive_1d() {
    let bounds = vec![(0.0, 6.0)];
    let total_budget = 25; // 5 initial + 20 BO iters
    let seeds: [u64; 3] = [0xC0FFEE, 0xDEADBEEF, 0xFEEDFACE];

    let mut bo_results: Vec<f64> = Vec::new();
    let mut rs_results: Vec<f64> = Vec::new();

    for &seed in &seeds {
        let cfg = BoConfig {
            n_initial: 5,
            n_iters: 20,
            n_candidates: 1024,
            xi: 0.01,
            seed,
        };
        let res = minimize(objective, bounds.clone(), cfg);
        let rs = random_search(seed, total_budget);
        println!(
            "bo_synthetic[seed=0x{seed:X}]: BO y_best = {:.6} at x = {:.6}; \
             random y_best = {:.6}",
            res.y_best, res.x_best[0], rs
        );
        bo_results.push(res.y_best);
        rs_results.push(rs);
    }

    let bo_best_across_seeds = bo_results.iter().cloned().fold(f64::INFINITY, f64::min);
    let rs_best_across_seeds = rs_results.iter().cloned().fold(f64::INFINITY, f64::min);
    let bo_mean = bo_results.iter().sum::<f64>() / (bo_results.len() as f64);
    let rs_mean = rs_results.iter().sum::<f64>() / (rs_results.len() as f64);
    println!(
        "bo_synthetic summary: BO best-across-seeds = {bo_best_across_seeds:.6} (mean {bo_mean:.6}); \
         random best-across-seeds = {rs_best_across_seeds:.6} (mean {rs_mean:.6})"
    );

    // Spec assertion 1: BO reaches below 0 (the valley around x ≈ 3.42).
    assert!(
        bo_best_across_seeds < 0.0,
        "BO best across seeds = {bo_best_across_seeds}, expected < 0.0 \
         (global minimum on [0,6] is y ≈ -0.8077 near x ≈ 3.42)"
    );

    // Spec assertion 2: BO beats random head-to-head across seeds. The
    // budget here (25 calls on a 1-D `[0,6]` interval) is large enough that
    // random search will usually land *somewhere* in the valley, but BO
    // concentrates evaluations and tightens onto the true minimum where
    // random misses by a measurable margin.
    assert!(
        bo_best_across_seeds < rs_best_across_seeds,
        "BO ({bo_best_across_seeds}) should beat random ({rs_best_across_seeds}) \
         across seeds (BO mean {bo_mean:.6}, random mean {rs_mean:.6})"
    );
}
