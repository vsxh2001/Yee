//! Aggressive space mapping (FS.5b.0, ADR-0213): design with a cheap
//! **coarse** model, correct with a handful of expensive **fine**-model
//! evaluations — Bandler's ASM, the textbook fit for Yee's closed-form +
//! full-wave-EM pair.
//!
//! The loop: `z_star` is the coarse-optimal design (the caller solves the
//! cheap problem however it likes). At each iterate the fine model is
//! evaluated once, the coarse model is **aligned** to that fine response
//! by parameter extraction (`z_k = argmin_z ‖coarse(z) − fine(x_k)‖`),
//! and the misalignment `e_k = z_k − z_star` drives a Broyden
//! quasi-Newton step. When the extraction of the fine response *is* the
//! coarse optimum, the fine design meets the spec as well as the coarse
//! one did — typically in 3–6 fine evaluations instead of the tens a
//! black-box optimizer spends.
//!
//! Everything here is deterministic (Gauss–Newton extraction, Broyden
//! updates, no randomness): identical inputs reproduce bit-for-bit.

use nalgebra::{DMatrix, DVector};

/// Configuration for [`extract`] (Gauss–Newton parameter extraction).
#[derive(Debug, Clone)]
pub struct ExtractConfig {
    /// Maximum Gauss–Newton iterations.
    pub max_iters: usize,
    /// Stop when the residual norm improves by less than this relative
    /// amount between iterations.
    pub rel_tol: f64,
    /// Relative finite-difference step for the Jacobian (scaled by
    /// `max(|z_i|, fd_floor)`).
    pub fd_step: f64,
    /// Absolute floor for the finite-difference step scale.
    pub fd_floor: f64,
}

impl Default for ExtractConfig {
    fn default() -> Self {
        Self {
            max_iters: 30,
            rel_tol: 1e-12,
            fd_step: 1e-6,
            fd_floor: 1e-9,
        }
    }
}

/// Align a coarse model to an observed response: Gauss–Newton on
/// `‖coarse(z) − target‖²` from `z0`, central-difference Jacobian (the
/// coarse model is cheap by definition), step-halving line search.
/// Returns the extracted parameters.
pub fn extract(
    coarse: &dyn Fn(&[f64]) -> Vec<f64>,
    target: &[f64],
    z0: &[f64],
    cfg: &ExtractConfig,
) -> Vec<f64> {
    let n = z0.len();
    let mut z = z0.to_vec();
    let mut r = residual(coarse, &z, target);
    let mut r_norm = r.norm();
    for _ in 0..cfg.max_iters {
        // Central-difference Jacobian of the residual.
        let mut jac = DMatrix::<f64>::zeros(target.len(), n);
        for j in 0..n {
            let h = cfg.fd_step * z[j].abs().max(cfg.fd_floor);
            let (mut zp, mut zm) = (z.clone(), z.clone());
            zp[j] += h;
            zm[j] -= h;
            let (rp, rm) = (residual(coarse, &zp, target), residual(coarse, &zm, target));
            for i in 0..target.len() {
                jac[(i, j)] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        // Normal equations; bail out on a singular (flat) Jacobian.
        let jtj = jac.transpose() * &jac;
        let jtr = jac.transpose() * &r;
        let Some(step) = jtj.lu().solve(&jtr) else {
            break;
        };
        // Step-halving line search on the residual norm.
        let mut alpha = 1.0;
        let mut improved = false;
        for _ in 0..20 {
            let z_try: Vec<f64> = z
                .iter()
                .zip(step.iter())
                .map(|(zi, s)| zi - alpha * s)
                .collect();
            let r_try = residual(coarse, &z_try, target);
            if r_try.norm() < r_norm {
                let gain = (r_norm - r_try.norm()) / r_norm.max(1e-300);
                z = z_try;
                r = r_try;
                r_norm = r.norm();
                improved = true;
                if gain < cfg.rel_tol {
                    return z;
                }
                break;
            }
            alpha *= 0.5;
        }
        if !improved {
            break;
        }
    }
    z
}

fn residual(model: &dyn Fn(&[f64]) -> Vec<f64>, z: &[f64], target: &[f64]) -> DVector<f64> {
    let y = model(z);
    DVector::from_iterator(target.len(), y.iter().zip(target).map(|(a, b)| a - b))
}

/// Configuration for [`space_map`].
#[derive(Debug, Clone)]
pub struct SpaceMapConfig {
    /// Maximum number of fine-model evaluations.
    pub max_fine_evals: usize,
    /// Convergence tolerance on `‖z_k − z_star‖`, measured relative to
    /// `scale` (per-component divisors, e.g. the nominal magnitudes).
    pub tol: f64,
    /// Per-component scale for the misalignment norm. Must match the
    /// parameter dimension.
    pub scale: Vec<f64>,
    /// Extraction settings.
    pub extract: ExtractConfig,
}

/// Result of an aggressive-space-mapping run.
#[derive(Debug, Clone)]
pub struct SpaceMapResult {
    /// The fine-space design at termination.
    pub x: Vec<f64>,
    /// Fine evaluations spent (the cost that matters).
    pub n_fine_evals: usize,
    /// Final scaled misalignment `‖(z_k − z_star)/scale‖`.
    pub misalignment: f64,
    /// Whether the tolerance was met within the budget.
    pub converged: bool,
}

/// Aggressive space mapping (Broyden variant): drive the fine design so
/// that the coarse-model extraction of its response equals the
/// coarse-optimal design `z_star`. Starts at `x0` (classically
/// `x0 = z_star`) with `B₀ = I`.
///
/// Panics if `scale`/`z_star`/`x0` dimensions disagree — caller bugs.
pub fn space_map(
    fine: &dyn Fn(&[f64]) -> Vec<f64>,
    coarse: &dyn Fn(&[f64]) -> Vec<f64>,
    z_star: &[f64],
    x0: &[f64],
    cfg: &SpaceMapConfig,
) -> SpaceMapResult {
    let n = z_star.len();
    assert_eq!(x0.len(), n, "x0 and z_star must have the same length");
    assert_eq!(cfg.scale.len(), n, "scale must match the parameter dim");
    assert!(cfg.max_fine_evals > 0, "need at least one fine evaluation");

    let mut x = DVector::from_row_slice(x0);
    let mut b = DMatrix::<f64>::identity(n, n);
    let mut e_prev: Option<DVector<f64>> = None;
    let mut n_fine = 0usize;

    loop {
        let y = fine(x.as_slice());
        n_fine += 1;
        // Extract from the coarse optimum — the best-informed start.
        let z = extract(coarse, &y, z_star, &cfg.extract);
        let e = DVector::from_iterator(n, z.iter().zip(z_star).map(|(a, b)| a - b));
        let mis = e
            .iter()
            .zip(&cfg.scale)
            .map(|(ei, s)| (ei / s) * (ei / s))
            .sum::<f64>()
            .sqrt();
        if mis < cfg.tol {
            return SpaceMapResult {
                x: x.as_slice().to_vec(),
                n_fine_evals: n_fine,
                misalignment: mis,
                converged: true,
            };
        }
        if n_fine >= cfg.max_fine_evals {
            return SpaceMapResult {
                x: x.as_slice().to_vec(),
                n_fine_evals: n_fine,
                misalignment: mis,
                converged: false,
            };
        }
        // Broyden update from the previous step (skip on the first).
        if let Some(ep) = &e_prev {
            // h = the step just taken; B += ((Δe − B h) hᵀ) / (hᵀ h).
            let h = -&b.clone().lu().solve(ep).expect("B stays invertible");
            let de = &e - ep;
            let denom = h.dot(&h);
            if denom > 0.0 {
                b += (de - &b * &h) * h.transpose() / denom;
            }
        }
        let step = b
            .clone()
            .lu()
            .solve(&e)
            .expect("Broyden matrix became singular");
        x -= step;
        e_prev = Some(e);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Coarse two-output model with an exactly known inverse.
    fn coarse(z: &[f64]) -> Vec<f64> {
        vec![10.0 / z[0], 6.0 / z[1]]
    }

    #[test]
    fn extraction_recovers_known_design() {
        let target = coarse(&[2.0, 3.0]);
        let z = extract(&coarse, &target, &[1.5, 2.0], &ExtractConfig::default());
        assert!(
            (z[0] - 2.0).abs() < 1e-9 && (z[1] - 3.0).abs() < 1e-9,
            "{z:?}"
        );
    }

    #[test]
    fn identical_models_converge_in_one_fine_eval() {
        let cfg = SpaceMapConfig {
            max_fine_evals: 5,
            tol: 1e-9,
            scale: vec![1.0, 1.0],
            extract: ExtractConfig::default(),
        };
        let z_star = [2.0, 3.0];
        let r = space_map(&coarse, &coarse, &z_star, &z_star, &cfg);
        assert!(r.converged);
        assert_eq!(r.n_fine_evals, 1, "fine = coarse must converge immediately");
    }

    #[test]
    fn shifted_fine_model_is_found() {
        // fine(x) = coarse(x − 0.1): the classic aligned-spaces case ASM
        // solves in a couple of Broyden steps.
        let fine = |x: &[f64]| coarse(&[x[0] - 0.1, x[1] - 0.1]);
        let cfg = SpaceMapConfig {
            max_fine_evals: 8,
            tol: 1e-8,
            scale: vec![1.0, 1.0],
            extract: ExtractConfig::default(),
        };
        let z_star = [2.0, 3.0];
        let r = space_map(&fine, &coarse, &z_star, &z_star, &cfg);
        assert!(r.converged, "misalignment {}", r.misalignment);
        assert!(
            (r.x[0] - 2.1).abs() < 1e-6 && (r.x[1] - 3.1).abs() < 1e-6,
            "{:?}",
            r.x
        );
        assert!(r.n_fine_evals <= 4, "took {} fine evals", r.n_fine_evals);
    }
}
