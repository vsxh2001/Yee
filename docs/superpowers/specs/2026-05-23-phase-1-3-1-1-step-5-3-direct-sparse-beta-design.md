# Phase 1.3.1.1 step 5.3 — direct β-direct-pencil sparse shift-and-invert solve

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.3 (close the inhomogeneous §4 gap).
**Depends on:** step 5.2 (β-direct extraction + hybrid; `3b6d899`),
step 5.1 (verified reference; `981926b`).
**Blocks:** CLAUDE.md §4 published-benchmark closure of the
*inhomogeneous* cross-section eigensolver.

## 1. Goal

Replace the step-5.2 **hybrid** mixed solve (cutoff-pencil select +
β-direct Rayleigh quotient on the *cutoff*-pencil eigenvector) with a
**direct** solve of the β-direct pencil `(k₀²B − A) x = β² B₁ x` via a
**faer sparse shift-and-invert**, recovering the *true* β-direct
eigenvector and enabling finer meshes. This resolves **both** residual
sources step 5.2 left open:

- **(b) RQ eigenvector mismatch** — using the true β-direct eigenvector
  makes β² exact for that mode (no cutoff-vs-β-direct angle error);
- **(a) discretization** — sparse LU makes finer meshes tractable
  (`solve_dense_mixed` is dense `O(n³)`, ~18 s at 12×12).

Target: numerical β within **≤5%** of the step-5.1 reference at a
representative contrast (ε_r=4.4 FR-4 primary; ε_r=10.2 RT/duroid 6010
stretch). The uniform-fill closed-form anchor stays exact; no gate
weakened.

## 2. Background

step 5.2 fixed the β-*extraction* (uniform-fill anchor exact, rel
1.5e-4) but the *inhomogeneous* mixed path ships a hybrid because the
naive direct β-direct solve (spec-5.2 "option A") **drifts**: the
β-direct pencil's spurious curl-free gradient modes land at
`β² ≈ k₀²⟨ε_r⟩` (top of spectrum), interleaved with the physical mode,
and `(K − σB₁) ≈ −A` is near-singular among them — a global eigen-sweep
thrashes on the cluster and a naive selection grabs a spurious
`E_z ≈ 0` branch. The hybrid sidesteps this by selecting on the cutoff
pencil, but then its β² is a Rayleigh quotient on the *wrong*
(cutoff-pencil) eigenvector — a mesh-stable ~17% bias at ε_r=10.2
(reviewed: ~34° eigenvector angle).

The fix is a **targeted** (not global) solve: shift-and-invert with a
shift `σ` placed near the *physical* β² isolates the physical eigenpair
and amplifies it away from the spurious cluster.

## 3. Approach

`crates/yee-mom/src/eigensolver/solve.rs` (+ faer plumbing; faer is
already a yee-mom dependency):

1. **Physics-informed shift `σ`.** Estimate the physical β² before the
   sparse solve. Cheapest: run the existing cutoff-pencil selection
   (Stage 1 of the step-5.2 hybrid) to get `k_c²` and the
   field-weighted `⟨ε_r⟩`, giving `σ₀ = (k₀² − k_c²)⟨ε_r⟩` (the hybrid's
   own β² estimate — known to be within ~17% of the true value, ample
   for a shift). The true physical β² is the eigenvalue of
   `(k₀²B − A) x = β² B₁ x` *nearest* `σ₀` whose eigenvector is
   transverse-energy-dominated.
2. **Sparse shift-and-invert.** Build `K − σ B₁` (with `K = k₀²B − A`)
   as a faer `SparseColMat`, factor once via `sp_lu`, and run
   inverse iteration `z ← (K − σB₁)⁻¹ B₁ z` to converge the eigenpair
   nearest `σ` — the **true β-direct eigenvector** `x_true`. Then
   `β² = R(x_true)` is exact for that mode (RQ stationary *at* its own
   eigenvector ⇒ no mismatch bias). Mirror the existing dense
   `inverse_iterate` (`solve.rs:546`) on the sparse operator.
3. **Spurious-mode guard.** If the converged mode is `E_z`-dominated or
   curl-free (gradient cluster), reject and re-shift / deflate. The
   transverse-energy filter already in the hybrid is the screen.
4. **Keep the dense path** as the small-`n` reference / fallback (its
   own tests + the uniform-fill anchor). Sparse is selected for the
   production / finer-mesh path.
5. **Mesh refinement.** With sparse LU, run the inhomogeneous case at
   finer meshes (e.g. 16×16, 24×24, 32×32) and report the β
   convergence toward the reference — this is the (a)-vs-(b)
   discriminator the step-5.2 review wanted: if direct-eigenvector β
   *still* plateaus short, (a) discretization dominates and even finer
   / higher-order is needed; if it converges to ≤5%, (b) was dominant
   and is now resolved.

## 4. Validation / DoD

- DoD-1. Direct sparse β-direct solve lands the physical mode (not the
  spurious branch) on the horizontal slab; reports β at ≥3 mesh
  densities.
- DoD-2. **§4 closure (primary, FR-4 ε_r=4.4):** numerical β within ≤5%
  of `reference.rs` `slab_loaded_beta(ε_r=4.4)`; the inhomogeneous
  reconciliation becomes a **failing gate** at this contrast.
- DoD-3. **Stretch (ε_r=10.2):** report the converged β vs reference
  583; if ≤5%, tighten that gate too; if not, document the residual +
  whether it is now discretization-limited (per §3.5) and queue
  step-5.4 (higher-order elements).
- DoD-4. **No regression:** uniform-fill anchor still ≤1%; ε_r=1 canary
  4e-14-class; `eigensolver_wr90` green; coupling guards green; Z_w
  reduction green.
- DoD-5. No new `Cargo.toml` dependency (faer already present). Lint
  floor clean. `reference.rs` untouched.

## 5. Risks

(a) **Shift `σ₀` lands in the spurious cluster.** Mitigation: `σ₀` is
the hybrid β² estimate (within ~17% of physical, well below the
`k₀²⟨ε_r⟩` cluster for the dominant mode); the transverse-energy filter
rejects a spurious capture; re-shift if needed.
(b) **faer sparse LU on the indefinite `K − σB₁`.** It is non-symmetric-
indefinite but nonsingular at a generic shift; `sp_lu` handles it (same
as the yee-fem `build_shifted` path). If faer's sparse LU rejects it,
fall back to a dense shifted LU at validation `n` and document the
sparse path as deferred.
(c) **Discretization still dominates at ε_r=10.2** even with the true
eigenvector. Then DoD-3 documents it + queues step-5.4 (higher-order /
curl-conforming p-refinement); DoD-2 (FR-4) still closes the §4 gap at a
representative contrast.

## 6. References

* Jin, *FEM in EM* 3rd ed. §8.4 (mixed waveguide eigenproblem).
* Saad, *Numerical Methods for Large Eigenvalue Problems* (shift-invert).
* `crates/yee-fem/src/solve.rs` — the `build_shifted` + faer `sp_lu`
  shift-invert pattern to mirror (different crate/pencil; pattern only).
* `crates/yee-mom/src/eigensolver/solve.rs` — `solve_dense_mixed`,
  `inverse_iterate`.
* `crates/yee-mom/src/eigensolver/reference.rs` — the validation oracle.
* ADR-0053 (step-5.2 as-built + the residual decomposition).
