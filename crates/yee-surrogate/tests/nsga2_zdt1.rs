//! NSGA-II on the canonical ZDT1 two-objective benchmark.
//!
//! ZDT1 (Zitzler-Deb-Thiele 2000, problem 1):
//!
//! - d = 30 variables, all in [0, 1].
//! - f1(x) = x[0]
//! - g(x)  = 1 + 9/(d-1) · Σ_{i=1..d-1} x[i]
//! - f2(x) = g(x) · (1 - sqrt(f1(x) / g(x)))
//!
//! Analytic Pareto front: f2 = 1 - sqrt(f1), f1 ∈ [0, 1].
//!
//! Validation metric: Inverted Generational Distance (IGD). For each point
//! on a uniformly sampled reference front, take the L2 distance to the
//! closest point in the optimizer's reported Pareto front, average over
//! the reference. Smaller is better; the canonical paper reaches
//! IGD ≈ 0.005 on ZDT1 at 250 generations. This test uses a 200-generation
//! budget at pop=100 and asserts IGD < 0.05 with a front of at least
//! 50 distinct individuals.

use nalgebra::DVector;
use yee_surrogate::{Nsga2Config, nsga2_minimize};

#[test]
fn nsga2_recovers_zdt1_pareto_front() {
    let d = 30;
    let zdt1 = |x: &DVector<f64>| -> Vec<f64> {
        let f1 = x[0];
        let g = 1.0 + 9.0 / (d as f64 - 1.0) * x.iter().skip(1).sum::<f64>();
        let f2 = g * (1.0 - (f1 / g).sqrt());
        vec![f1, f2]
    };

    let bounds = vec![(0.0, 1.0); d];
    let cfg = Nsga2Config {
        population_size: 100,
        n_generations: 200,
        crossover_eta: 20.0,
        mutation_eta: 20.0,
        mutation_probability: 1.0 / d as f64,
        seed: 0xC0FFEE,
    };

    let result = nsga2_minimize(zdt1, bounds, 2, cfg);

    // Extract the optimizer's reported Pareto front in objective space.
    let front: Vec<(f64, f64)> = result
        .pareto_front_indices
        .iter()
        .map(|&i| (result.objectives[i][0], result.objectives[i][1]))
        .collect();

    // 100 evenly spaced reference points on the analytic front.
    let reference: Vec<(f64, f64)> = (0..100)
        .map(|i| {
            let f1 = i as f64 / 99.0;
            (f1, 1.0 - f1.sqrt())
        })
        .collect();

    // Inverted Generational Distance.
    let igd: f64 = reference
        .iter()
        .map(|(rf1, rf2)| {
            front
                .iter()
                .map(|(f1, f2)| ((f1 - rf1).powi(2) + (f2 - rf2).powi(2)).sqrt())
                .fold(f64::INFINITY, f64::min)
        })
        .sum::<f64>()
        / reference.len() as f64;

    println!("ZDT1 IGD       = {igd:.6}");
    println!("front size     = {}", front.len());

    // The NSGA-II paper reaches IGD ≈ 0.005 at 250 generations; the 200-gen
    // budget here should clear < 0.05 comfortably.
    assert!(igd < 0.05, "ZDT1 IGD too large: {igd}");
    assert!(front.len() >= 50, "front too sparse: {}", front.len());
}
