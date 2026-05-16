# ADR-0012: Implement NSGA-II as a separate `nsga2` module, not under `bo`

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

`ROADMAP.md` originally listed Bayesian optimisation (BO) and NSGA-II
multi-objective evolutionary search **in the same bullet** under the
Phase 3 "design-loop" theme. ADR-0010 already explained why BO shipped
first (Phase 3.bo.0, single-objective Expected Improvement): the
walking-skeleton path was shorter and the user pull on the patch-
antenna sweep was single-objective in practice.

Phase 3.bo.1 is now in the queue, and the question is **where** NSGA-II
should land in the `yee-surrogate` crate layout. Two structural options
present themselves:

1. **Under `bo`.** Treat NSGA-II as a sibling of the
   Expected-Improvement loop inside the existing `bo` module. The
   public surface would be `bo::minimize` (single-objective EI) and
   `bo::pareto` (multi-objective NSGA-II), both presented to users as
   "design-optimisation entry points".
2. **As a separate `nsga2` module.** Treat NSGA-II as algorithmically
   distinct from BO and ship it as `yee_surrogate::nsga2::optimize`
   alongside `yee_surrogate::bo::minimize`, with both crates' public
   surfaces re-exported at the crate root.

The structural argument turns on whether NSGA-II shares enough with BO
to live next to it:

- **BO uses a surrogate.** It fits a Gaussian Process (ADR-0009),
  evaluates an acquisition function (Expected Improvement) over
  thousands of candidate points per iteration, and queries the
  expensive objective only at the argmax of the acquisition. The
  surrogate **is** the algorithm.
- **NSGA-II uses no surrogate.** It is a pure evolutionary
  algorithm: maintain a population of `N` candidates, generate
  `N` offspring by simulated binary crossover (SBX) plus polynomial
  mutation, rank parents + offspring by fast non-dominated sort,
  break ties by crowding distance, take the top `N` as the next
  generation. The expensive objective is evaluated at every offspring
  on every generation — there is no notion of "the surrogate told us
  which design to try". A surrogate-assisted NSGA-II variant exists
  in the literature (SAEA-NSGA-II, K-RVEA), but it is a future
  research direction, not the textbook algorithm.

So the algorithmic overlap between BO and NSGA-II is **near-zero**.
What they share is plumbing:

- A reproducible RNG (xorshift64 with a fixed seed).
- Latin-hypercube sampling for an initial design.
- The `f: FnMut(&[f64]) -> f64` user-callback shape (extended to
  `Fn(&[f64]) -> Vec<f64>` for the multi-objective case).
- The `bounds: &[(f64, f64)]` parameter-box convention.

These are small shared helpers, not a shared algorithm. Putting them
in `bo` and asking NSGA-II to import them upward creates a
module-layering inversion (`bo` is supposed to be a leaf, not a
utility module). Putting NSGA-II into `bo` because they share helpers
miscommunicates the relationship between the two algorithms to users
reading the crate.

The validation case for NSGA-II is **ZDT1** (Zitzler-Deb-Thiele
benchmark 1), a 30-dimensional bi-objective problem with a known
analytical Pareto front (`f₁ ∈ [0, 1]`, `f₂ = 1 − √f₁`). The
quality of an approximated front is measured by **Inverted
Generational Distance (IGD)**: for each point on the true front, find
its closest point in the approximation; take the mean of those
distances. Lower is better; IGD < 0.05 is the conventional "good
approximation" threshold cited across the NSGA-II / MOEA/D / SMS-EMOA
benchmarking literature.

## Decision

Ship NSGA-II as a separate top-level module **`yee_surrogate::nsga2`**,
not nested under `bo`:

```text
crates/yee-surrogate/src/
  bo/
    mod.rs            — single-objective BO (ADR-0010)
    ei.rs             — Expected Improvement acquisition
    lhs.rs            — Latin-hypercube initial design
    rng.rs            — xorshift64 reproducible RNG
  nsga2/
    mod.rs            — pub fn optimize(...)
    sort.rs           — fast non-dominated sort + crowding distance
    sbx.rs            — simulated binary crossover
    mutation.rs       — polynomial mutation
    replacement.rs    — elitist (μ + λ) survivor selection
```

The public surface:

```rust
pub struct Nsga2Config {
    pub population_size: usize,    // N (parents) = N (offspring)
    pub n_generations: usize,
    pub eta_c: f64,                // SBX distribution index (default 20)
    pub eta_m: f64,                // polynomial mutation eta (default 20)
    pub p_crossover: f64,          // SBX probability (default 0.9)
    pub p_mutation: f64,           // mutation prob per gene (default 1/d)
    pub seed: u64,
}

pub struct Nsga2Result {
    pub pareto_front_x: Vec<Vec<f64>>,
    pub pareto_front_y: Vec<Vec<f64>>,
    pub history: Vec<Generation>,
}

pub fn optimize<F>(
    f: F,
    bounds: &[(f64, f64)],
    n_objectives: usize,
    cfg: Nsga2Config,
) -> Result<Nsga2Result>
where
    F: FnMut(&[f64]) -> Vec<f64>;
```

**Algorithmic choices follow Deb et al. 2002 verbatim:**

- **Initial population.** Latin-hypercube of `N` points in the
  bounded box (reuses `bo::lhs`).
- **Selection.** Binary tournament on (rank, crowding distance).
- **Crossover.** Simulated Binary Crossover (SBX) with distribution
  index `η_c = 20` and per-pair probability 0.9.
- **Mutation.** Polynomial mutation with `η_m = 20` and per-gene
  probability `1/d` (one expected mutation per offspring).
- **Survivor selection.** Combine parents + offspring (2N
  candidates), fast non-dominated sort, fill the next generation
  rank-by-rank, break the boundary rank by descending crowding
  distance until exactly `N` survive. This is the canonical
  elitist `(μ + λ) → μ` step.
- **Termination.** Fixed generation budget. No convergence detector
  in Phase 3.bo.1; users size `n_generations` from their budget.

**Shared helpers.** `nsga2` re-uses `bo::rng::Xorshift64` and
`bo::lhs::latin_hypercube` directly (these are `pub(crate)` already
under `bo`). Nothing else is shared.

**Validation gate.** `crates/yee-surrogate/validation/nsga2_zdt1.rs`
runs NSGA-II on ZDT1 with `N = 100, n_generations = 250` and
asserts **IGD < 0.05** against an analytical Pareto front of 500
sampled points.

## Consequences

**What becomes easier:**

- Multi-objective design problems are now solvable in pure Rust from
  the `yee-surrogate` crate without reaching for an external
  evolutionary-algorithm library. The intended user case is
  "trade size vs bandwidth vs gain on a patch antenna" — three
  objectives, ~5 design variables, a few hundred MoM solves of
  budget.
- The module layout matches the algorithmic reality: BO and NSGA-II
  are sibling top-level modules, both consumed through their
  `optimize` / `minimize` entry points. Users reading the crate root
  see two design-loop tools, not one with a confusing sub-mode.
- The shared `bo::lhs` and `bo::rng` helpers stay in `bo` (where they
  were already, for the single-objective LHS-init use). NSGA-II
  imports them with explicit `use bo::lhs::latin_hypercube;`, which
  is a clear and visible coupling.

**What becomes harder:**

- Users have to **pick** between `bo::minimize` and
  `nsga2::optimize`. Random-search candidate maximisation of EI
  (in `bo`) and SBX offspring generation (in `nsga2`) are
  algorithmically different beasts, and the crate documentation has
  to make the "single objective → BO, multi-objective → NSGA-II"
  decision rule explicit. The README and `yee-py` Python wrapper
  both grow a "which optimiser do I want?" decision-tree section.
- Surrogate-assisted multi-objective optimisation
  (SAEA-NSGA-II, K-RVEA) is not in scope for Phase 3.bo.1. A user
  who wants to use the GP from ADR-0009 to filter offspring before
  evaluation has to wait for a future ADR (Phase 3.bo.2 candidate).
- The decision to keep `bo::lhs` and `bo::rng` as shared helpers
  inside `bo` is a small layering compromise: a hypothetical future
  `cma_es` module would also want to import them, at which point
  they should be promoted to a top-level `common` module. That
  refactor is deferred until a second importer actually exists.

**What's now closed off:**

- Combining BO and NSGA-II under a single `bo` module umbrella in
  the public API. The crate root re-exports two distinct entry
  points; future tooling (CLI subcommands, Python bindings) must
  reflect the two-tool shape.
- Bringing in an external evolutionary-algorithm crate
  (`oxigen`, `revolution`, `metaheuristics-nature`) for Phase
  3.bo.1. The hand-rolled NSGA-II keeps the dependency footprint
  flat and the pure-Rust / no-system-libs property intact for the
  surrogate path (same argument as ADR-0009).
- Re-using BO's random-search EI candidate maximiser as NSGA-II's
  offspring generator. They are different patterns and the code
  reads more clearly with two implementations than one shared one.

## References

- `crates/yee-surrogate/src/nsga2/` — the implementation.
- `crates/yee-surrogate/src/nsga2/sort.rs` — fast non-dominated
  sort + crowding distance, Deb 2002 Algorithm 1 and §III.B.
- `crates/yee-surrogate/src/nsga2/sbx.rs` — simulated binary
  crossover, Deb & Agrawal 1995.
- `crates/yee-surrogate/src/nsga2/mutation.rs` — polynomial
  mutation, Deb & Goyal 1996.
- `crates/yee-surrogate/validation/nsga2_zdt1.rs` — ZDT1 with
  IGD < 0.05 gate.
- K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, "A Fast and Elitist
  Multiobjective Genetic Algorithm: NSGA-II", *IEEE Trans.
  Evolutionary Computation*, vol. 6, no. 2, pp. 182-197, 2002.
- E. Zitzler, K. Deb, L. Thiele, "Comparison of Multiobjective
  Evolutionary Algorithms: Empirical Results", *Evolutionary
  Computation*, vol. 8, no. 2, pp. 173-195, 2000. (ZDT benchmarks.)
- ADR-0009 — GP surrogate (not used by NSGA-II in Phase 3.bo.1;
  may be used by a future SAEA variant).
- ADR-0010 — single-objective BO; the sibling module.
