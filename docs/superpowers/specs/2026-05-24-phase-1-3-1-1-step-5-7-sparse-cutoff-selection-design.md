# Phase 1.3.1.1 step 5.7 — sparse cutoff-mode selection (unlock mesh scaling)

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.7 (foundational: mesh scaling + the last ε_r=10.2 lever).
**Depends on:** step 5.6 (p=2-robust converged-eigenvector selection; merge `f75fc87`),
step 5.3 (faer sparse-LU β-direct shift-invert; `6eca76a`), step 4
(`SparseEigen`/`LobpcgEigen` pattern in yee-fem; `4c2f4e1`).
**Blocks:** ε_r=10.2 ≤5% closure (needs finer meshes than the dense cap
affords) AND realistic-size cross-section solves generally.

## 1. Goal

Replace the **dense `O(n³)` `complex_eigenvalues` cutoff-pencil
eigendecomposition** inside `cutoff_candidates` (`solve.rs`) with a
**sparse shift-and-invert** that returns the low-cutoff candidate modes
without forming the full dense spectrum. This lifts the ~6×6 (n≈457)
mesh cap that limits the *entire* cross-section eigensolver, enabling
the finer meshes where (a) the p=2 higher-order convergence manifests
(closing ε_r=10.2 ≤5%) and (b) realistic cross-sections become solvable.

No gate weakened; FR-4 / uniform / homogeneous bit-identical (the sparse
path must agree with the dense path at the existing meshes).

## 2. Background

step 5.6 fixed mode *selection* but found p=2 ≈ p1 (~16.6% at ε_r=10.2)
because both plateau at the **dense-cap mesh** (~6×6) — the dense
`complex_eigenvalues` (`cutoff_candidates`, solve.rs) is `O(n³)` and caps
practical meshes there. The p=2 advantage + the high-contrast
convergence need finer meshes, which need a sparse selection.

**The challenge — the gradient null cluster.** The cutoff pencil
`A x = k_c² B x` (A = curl-curl stiffness, B = ε_r mass) has a large
curl-free gradient null space at `k_c² ≈ 0` (every nodal gradient;
*enlarged* at p=2 by `∇(λ_aλ_b)`). A naive sparse solve near 0 is
dominated by that cluster. The physical dominant quasi-TEM mode sits at
a small but **strictly positive** `k_c²` *above* the gradient floor.

## 3. Approach

`crates/yee-mom/src/eigensolver/solve.rs` — `cutoff_candidates`:

1. **Sparse shift-invert above the gradient floor.** Build the cutoff
   pencil as faer sparse matrices; factor `(A − σ B)` once via `sp_lu`
   at a **small positive shift** `σ` placed above the gradient null
   floor and below the physical dominant cutoff (e.g.
   `σ = small_fraction · (k₀²·ε_r,max)`, or a fraction of the analytic
   air-cutoff — the physical low-order modes' k_c² are O((π/a)²)). Find
   the few eigenpairs nearest σ via inverse iteration / a block LOBPCG
   (reuse the step-4 `LobpcgEigen` pattern + the step-5.3 sparse LU),
   yielding the low-cutoff *physical* candidates while the gradient
   cluster at `k_c² < σ` is excluded by the shift.
2. **Hand the candidates to the existing step-5.6 selection unchanged:**
   shift-invert each candidate's β-direct RQ, screen the converged
   β-direct eigenvector for transverse-dominance, keep highest-β². The
   selection logic from step 5.6 is order- and sparsity-agnostic.
3. **Keep the dense `cutoff_candidates` as a small-`n` fallback /
   reference** (selected by a DoF threshold, or behind the dense
   `solve_dense_mixed_rq`); the sparse path is the production selection
   for large `n`. The sparse path must reproduce the dense candidates
   (same dominant mode) at the validation meshes.

## 4. Validation / DoD

- DoD-1. Sparse `cutoff_candidates` returns the same dominant mode as the
  dense path at the existing meshes (FR-4, homogeneous, uniform,
  vertical-slab) — those gates bit-identical or within a tight tol.
- DoD-2. A finer-mesh convergence study (now affordable): ε_r=10.2 at
  e.g. 12×12 / 16×16 / 24×24 → reference 582.95. **§4 high-contrast
  closure if ≤5% is reached** (promote to a failing gate); else document
  the improved p=2 convergence trend + queue step-5.8 (the residual
  would then be a genuinely deeper modeling question).
- DoD-3. Mesh-scaling demonstrated: a cross-section solve at n well
  beyond the old ~457 dense cap completes in reasonable time (a
  `--release` example or a timed test).
- DoD-4. No regression (FR-4, uniform anchor, ε_r=1 canary, wr90,
  coupling guards, Z_w). No new `Cargo.toml` dependency (faer + the
  yee-fem LOBPCG pattern; if a cross-crate reuse is awkward, re-implement
  the small block-LOBPCG in yee-mom — no new dep). Lint clean.
  `reference.rs` untouched.

## 5. Risks

(a) **Shift σ placement** — too low hits the gradient cluster, too high
misses the dominant mode. Mitigation: σ as a fraction of the analytic
air-cutoff (k_c² ~ (π/a)²), validated against the dense path at the
existing meshes (DoD-1); make σ adaptive if needed (raise until the
returned modes are transverse-dominated).
(b) **Gradient-cluster contamination of the sparse iteration** even with
a positive shift (the cluster is large at p=2). Mitigation: gradient-
deflation (project against the discrete-gradient range) is the standard
Nedelec remedy if the shift alone is insufficient — but try the shift
first (simpler). Escape-hatch if neither suffices in budget.
(c) **ε_r=10.2 still short even at finer mesh** — then the residual is
deeper than discretization (re-examine the reference / formulation);
DoD-2's document-and-queue branch. Low probability (step 5.4/5.6 point
squarely at mesh resolution).

## 6. References

* Saad, *Numerical Methods for Large Eigenvalue Problems* (shift-invert,
  spectral transformations). Knyazev 2001 (LOBPCG).
* Jin §9 (Nedelec gradient null space, tree-cotree / gradient deflation).
* `crates/yee-fem/src/solve.rs` (`LobpcgEigen`/`build_shifted` pattern),
  `crates/yee-mom/src/eigensolver/solve.rs` (`cutoff_candidates`,
  `beta_direct_shift_invert`), `reference.rs`.
* ADR-0054/0056/0057 (the β-direct solve + selection chain).
