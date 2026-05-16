//! NSGA-II multi-objective optimizer (Deb 2002).
//!
//! Non-dominated Sorting Genetic Algorithm II for minimizing `m` objectives
//! `f: ℝ^d → ℝ^m` over a hyper-rectangle, returning the final population and
//! the indices of the Pareto-optimal set.
//!
//! ## Algorithm
//!
//! 1. Initialize a Latin-hypercube population `P` of size `N`.
//! 2. Evaluate all objectives.
//! 3. For each generation:
//!    a. Binary tournament selection on `(rank, crowding)` lexicographic order.
//!    b. SBX (Simulated Binary Crossover) + polynomial mutation → offspring `Q` of size `N`.
//!    c. Combine `P ∪ Q` → `R` of size `2N`.
//!    d. Fast non-dominated sort `R` into fronts `F1, F2, …`.
//!    e. Build new `P` from fronts in order; when adding the next front would
//!       exceed `N`, partial-select from the boundary front by descending
//!       crowding distance.
//! 4. Return final `P` plus the indices of `F1`.
//!
//! Reference: Deb, Pratap, Agarwal, Meyarivan, "A fast and elitist
//! multiobjective genetic algorithm: NSGA-II", IEEE Trans. Evol. Comput. 6(2),
//! 2002.
//!
//! ## Randomness
//!
//! All randomness flows through the seedable [`Xorshift64`] PRNG borrowed from
//! [`crate::bo`]. Seed via [`Nsga2Config::seed`] for reproducibility.

use nalgebra::DVector;

use crate::bo::Xorshift64;

/// Configuration for [`minimize`].
#[derive(Debug, Clone)]
pub struct Nsga2Config {
    /// Population size `N`. Offspring is sampled at the same size each
    /// generation, so each generation evaluates the objective `N` times after
    /// the initial design.
    pub population_size: usize,
    /// Number of generations to run.
    pub n_generations: usize,
    /// SBX (Simulated Binary Crossover) distribution index `η_c`. Larger
    /// values produce children closer to the parents.
    pub crossover_eta: f64,
    /// Polynomial mutation distribution index `η_m`. Larger values produce
    /// mutated values closer to the original.
    pub mutation_eta: f64,
    /// Per-gene mutation probability. The canonical default is `1 / d`.
    pub mutation_probability: f64,
    /// RNG seed for the initial design, selection, crossover, and mutation.
    pub seed: u64,
}

impl Default for Nsga2Config {
    fn default() -> Self {
        Self {
            population_size: 100,
            n_generations: 100,
            crossover_eta: 20.0,
            mutation_eta: 20.0,
            mutation_probability: 0.01,
            seed: 0xC0FFEE,
        }
    }
}

/// Result of an NSGA-II [`minimize`] run.
#[derive(Debug, Clone)]
pub struct Nsga2Result {
    /// Final population, in arbitrary order.
    pub population: Vec<DVector<f64>>,
    /// Objective vectors `f(x)` for each individual in [`Self::population`].
    /// Each entry has length `m = n_objectives`.
    pub objectives: Vec<DVector<f64>>,
    /// Indices into [`Self::population`] that lie on the non-dominated front
    /// (front F1).
    pub pareto_front_indices: Vec<usize>,
}

/// Stub entry point. The full NSGA-II pipeline lands in follow-up commits;
/// for now this only validates the inputs so the public API compiles.
pub fn minimize<F>(
    _objectives: F,
    bounds: Vec<(f64, f64)>,
    n_objectives: usize,
    cfg: Nsga2Config,
) -> Nsga2Result
where
    F: Fn(&DVector<f64>) -> Vec<f64>,
{
    assert!(!bounds.is_empty(), "nsga2: bounds must be non-empty");
    assert!(n_objectives >= 1, "nsga2: n_objectives must be ≥ 1");
    assert!(
        cfg.population_size >= 2,
        "nsga2: population_size must be ≥ 2"
    );
    // Touch the RNG so `Xorshift64` is exercised in tree from this module.
    let _ = Xorshift64::new(cfg.seed);
    Nsga2Result {
        population: Vec::new(),
        objectives: Vec::new(),
        pareto_front_indices: Vec::new(),
    }
}
