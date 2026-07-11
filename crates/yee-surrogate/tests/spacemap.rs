//! Gate `surrogate-sm-001` (FS.5b.0, ADR-0213): aggressive space mapping
//! beats direct BO on the same fine-evaluation budget — the roadmap FS.5
//! gate in its closed-form edition (the EM-fine R.4-scenario gate is
//! FS.5b.1).
//!
//! Testcase: patch two-mode design. Response `[f10, f01]` GHz with
//! `f10 = c/(2 L √εe)`, `f01 = c/(2 W √εe)`. The coarse model uses the
//! zeroth-order `εe = (εr + 1)/2`; the fine model adds a
//! Hammerstad–Jensen-style width dependence AND a fringing length
//! extension ΔL per side — a physically shaped space warp, not a toy
//! offset. Spec: (2.45, 3.10) GHz. Instant, non-ignored.

use nalgebra::DVector;
use yee_surrogate::{BoConfig, ExtractConfig, SpaceMapConfig, space_map};

const C0: f64 = 299_792_458.0;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;

fn coarse(z: &[f64]) -> Vec<f64> {
    let ee = (EPS_R + 1.0) / 2.0;
    vec![
        C0 / (2.0 * z[0] * ee.sqrt()) / 1e9,
        C0 / (2.0 * z[1] * ee.sqrt()) / 1e9,
    ]
}

/// Hammerstad–Jensen-style ε_eff(W/h) + Hammerstad fringing ΔL — the
/// "expensive full-wave solver" stand-in. Same physics family as the
/// coarse model, warped the way real patches are.
fn fine(x: &[f64]) -> Vec<f64> {
    let (l, w) = (x[0], x[1]);
    let ee = |width: f64| -> f64 {
        (EPS_R + 1.0) / 2.0 + (EPS_R - 1.0) / 2.0 / (1.0 + 12.0 * H_M / width).sqrt()
    };
    let dl = |width: f64, e: f64| -> f64 {
        0.412 * H_M * ((e + 0.3) * (width / H_M + 0.264)) / ((e - 0.258) * (width / H_M + 0.8))
    };
    // The 10-mode resonates along L (W is the radiating width) and vice
    // versa for the 01-mode.
    let (e_l, e_w) = (ee(w), ee(l));
    vec![
        C0 / (2.0 * (l + 2.0 * dl(w, e_l)) * e_l.sqrt()) / 1e9,
        C0 / (2.0 * (w + 2.0 * dl(l, e_w)) * e_w.sqrt()) / 1e9,
    ]
}

/// Coarse model inverted exactly for the spec — z_star in closed form.
fn coarse_optimum(spec: &[f64]) -> Vec<f64> {
    let ee = (EPS_R + 1.0) / 2.0;
    vec![
        C0 / (2.0 * spec[0] * 1e9 * ee.sqrt()),
        C0 / (2.0 * spec[1] * 1e9 * ee.sqrt()),
    ]
}

fn spec_err_pct(y: &[f64], spec: &[f64]) -> f64 {
    y.iter()
        .zip(spec)
        .map(|(a, b)| ((a - b) / b * 100.0) * ((a - b) / b * 100.0))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn surrogate_sm_001_asm_beats_direct_bo_on_equal_budget() {
    let spec = [2.45, 3.10];
    let z_star = coarse_optimum(&spec);
    let budget = 5;

    // --- Aggressive space mapping.
    let cfg = SpaceMapConfig {
        max_fine_evals: budget,
        tol: 1e-4,
        scale: z_star.clone(),
        extract: ExtractConfig::default(),
    };
    let sm = space_map(&fine, &coarse, &z_star, &z_star, &cfg);
    let sm_err = spec_err_pct(&fine(&sm.x), &spec);

    // --- Direct BO on the fine mismatch, same fine-evaluation budget.
    // (n_initial + n_iters = budget; bounds bracket the coarse optimum
    // generously — BO gets a fair search box.)
    let obj = |x: &DVector<f64>| -> f64 {
        let y = fine(x.as_slice());
        y.iter()
            .zip(&spec)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
    };
    let bounds = vec![
        (0.5 * z_star[0], 1.5 * z_star[0]),
        (0.5 * z_star[1], 1.5 * z_star[1]),
    ];
    let bo = yee_surrogate::minimize(
        obj,
        bounds,
        BoConfig {
            n_initial: 3,
            n_iters: budget - 3,
            ..BoConfig::default()
        },
    );
    let bo_err = spec_err_pct(&fine(bo.x_best.as_slice()), &spec);

    // Measured 2026-07-11: ASM 0.00143 % spec error in 4 fine evals
    // (converged); direct BO at the same 5-eval budget 44.8 % — the
    // coarse-model alignment is worth ~4 orders of magnitude here.
    // Asserts are set far inside those margins.
    assert!(
        sm_err <= 0.1,
        "ASM must meet spec to 0.1%: got {sm_err:.4}% in {} fine evals (converged: {})",
        sm.n_fine_evals,
        sm.converged
    );
    assert!(
        sm.n_fine_evals <= budget,
        "budget exceeded: {}",
        sm.n_fine_evals
    );
    assert!(
        bo_err >= 5.0 * sm_err,
        "space mapping must beat BO >=5x on equal budget: ASM {sm_err:.4}% vs BO {bo_err:.4}%"
    );
}

#[test]
fn surrogate_sm_001_is_deterministic() {
    let spec = [2.45, 3.10];
    let z_star = coarse_optimum(&spec);
    let cfg = SpaceMapConfig {
        max_fine_evals: 5,
        tol: 1e-4,
        scale: z_star.clone(),
        extract: ExtractConfig::default(),
    };
    let a = space_map(&fine, &coarse, &z_star, &z_star, &cfg);
    let b = space_map(&fine, &coarse, &z_star, &z_star, &cfg);
    assert_eq!(a.x, b.x);
    assert_eq!(a.n_fine_evals, b.n_fine_evals);
    assert_eq!(a.misalignment, b.misalignment);
}
