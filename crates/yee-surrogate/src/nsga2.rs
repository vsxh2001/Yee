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

/// Pareto-dominance test for minimization.
///
/// Returns `true` iff `a` dominates `b`: every component of `a` is `≤` the
/// corresponding component of `b`, and at least one is strictly less.
fn dominates(a: &[f64], b: &[f64]) -> bool {
    debug_assert_eq!(a.len(), b.len());
    let mut strictly_better = false;
    for (ai, bi) in a.iter().zip(b.iter()) {
        if ai > bi {
            return false;
        }
        if ai < bi {
            strictly_better = true;
        }
    }
    strictly_better
}

/// Fast non-dominated sort (Deb 2002, §III.A).
///
/// Input: a slice of objective vectors (`m`-dim each, minimization).
///
/// Output: a list of fronts `[F0, F1, …]` where `F0` is the Pareto-optimal
/// set. Each front is a list of indices into the input slice.
///
/// Runtime: `O(M · N²)` time, `O(N²)` space, where `N` is the input length and
/// `M` is the objective count.
pub(crate) fn fast_non_dominated_sort(objectives: &[DVector<f64>]) -> Vec<Vec<usize>> {
    let n = objectives.len();
    if n == 0 {
        return Vec::new();
    }

    // For each i: indices of solutions that i dominates, and a count of how
    // many dominate i.
    let mut dominated: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut domination_count: Vec<usize> = vec![0; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let oi = objectives[i].as_slice();
            let oj = objectives[j].as_slice();
            if dominates(oi, oj) {
                dominated[i].push(j);
                domination_count[j] += 1;
            } else if dominates(oj, oi) {
                dominated[j].push(i);
                domination_count[i] += 1;
            }
        }
    }

    let mut fronts: Vec<Vec<usize>> = Vec::new();
    let mut current: Vec<usize> = (0..n).filter(|&i| domination_count[i] == 0).collect();

    while !current.is_empty() {
        let mut next: Vec<usize> = Vec::new();
        for &i in &current {
            // Iterate by index to avoid a simultaneous mutable + immutable
            // borrow of `dominated`.
            for k in 0..dominated[i].len() {
                let j = dominated[i][k];
                domination_count[j] -= 1;
                if domination_count[j] == 0 {
                    next.push(j);
                }
            }
        }
        fronts.push(current);
        current = next;
    }

    fronts
}

/// Crowding distance assignment for a single front (Deb 2002, §III.B).
///
/// For each objective, sort the front by that objective, set the boundary
/// points' contribution to `+∞`, and add the normalized neighbor gap for
/// every interior point. The returned vector has one entry per element of
/// `front`, in the same order.
///
/// If the front has `≤ 2` points all distances are `+∞`. If an objective's
/// range is zero across the front, that objective contributes nothing
/// (avoids `0 / 0`).
pub(crate) fn crowding_distance(front: &[usize], objectives: &[DVector<f64>]) -> Vec<f64> {
    let len = front.len();
    let mut dist = vec![0.0_f64; len];
    if len == 0 {
        return dist;
    }
    if len <= 2 {
        for d in &mut dist {
            *d = f64::INFINITY;
        }
        return dist;
    }
    let m = objectives[front[0]].len();
    for k in 0..m {
        // Indices into `front` sorted by objective k.
        let mut order: Vec<usize> = (0..len).collect();
        order.sort_by(|&a, &b| {
            objectives[front[a]][k]
                .partial_cmp(&objectives[front[b]][k])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let fmin = objectives[front[order[0]]][k];
        let fmax = objectives[front[order[len - 1]]][k];
        dist[order[0]] = f64::INFINITY;
        dist[order[len - 1]] = f64::INFINITY;
        let range = fmax - fmin;
        if range <= 0.0 {
            continue;
        }
        for i in 1..(len - 1) {
            let prev = objectives[front[order[i - 1]]][k];
            let next = objectives[front[order[i + 1]]][k];
            dist[order[i]] += (next - prev) / range;
        }
    }
    dist
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

#[cfg(test)]
mod tests {
    use super::*;

    fn dv(values: &[f64]) -> DVector<f64> {
        DVector::from_row_slice(values)
    }

    #[test]
    fn dominates_basic() {
        assert!(dominates(&[0.0, 0.0], &[1.0, 1.0]));
        assert!(dominates(&[0.0, 1.0], &[1.0, 1.0]));
        assert!(!dominates(&[1.0, 1.0], &[1.0, 1.0]));
        assert!(!dominates(&[0.0, 2.0], &[1.0, 1.0]));
    }

    #[test]
    fn fast_non_dominated_sort_two_fronts() {
        // F0 = {(0,2), (1,1), (2,0)} mutually non-dominated.
        // F1 = {(1,2), (2,1)} each dominated by exactly one F0 point.
        let objs = vec![
            dv(&[0.0, 2.0]),
            dv(&[1.0, 1.0]),
            dv(&[2.0, 0.0]),
            dv(&[1.0, 2.0]),
            dv(&[2.0, 1.0]),
        ];
        let fronts = fast_non_dominated_sort(&objs);
        assert_eq!(fronts.len(), 2);
        let mut f0 = fronts[0].clone();
        f0.sort();
        assert_eq!(f0, vec![0, 1, 2]);
        let mut f1 = fronts[1].clone();
        f1.sort();
        assert_eq!(f1, vec![3, 4]);
    }

    #[test]
    fn fast_non_dominated_sort_empty() {
        let fronts = fast_non_dominated_sort(&[]);
        assert!(fronts.is_empty());
    }

    #[test]
    fn crowding_distance_boundary_infinite_interior_finite() {
        // Front of 4 points along a line in objective space: boundaries get
        // ∞, interior gets a finite positive number.
        let objs = vec![
            dv(&[0.0, 3.0]),
            dv(&[1.0, 2.0]),
            dv(&[2.0, 1.0]),
            dv(&[3.0, 0.0]),
        ];
        let front: Vec<usize> = (0..4).collect();
        let cd = crowding_distance(&front, &objs);
        assert!(cd[0].is_infinite());
        assert!(cd[3].is_infinite());
        assert!(cd[1].is_finite() && cd[1] > 0.0);
        assert!(cd[2].is_finite() && cd[2] > 0.0);
    }

    #[test]
    fn crowding_distance_small_front_is_infinite() {
        let objs = vec![dv(&[0.0, 1.0]), dv(&[1.0, 0.0])];
        let cd = crowding_distance(&[0, 1], &objs);
        assert!(cd.iter().all(|d| d.is_infinite()));
    }
}
