# Phase 1.3.1.1 step 5.4 — interface-graded mesh for high-contrast inhomogeneous accuracy

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.4 (close the high-contrast inhomogeneous residual).
**Depends on:** step 5.3 (direct β-direct sparse solve; merge `6eca76a`),
which made finer meshes tractable and proved the residual is
discretization-dominated.
**Blocks:** tightening the §4 inhomogeneous gate from FR-4 (closed) to
high-contrast (ε_r=10.2).

## 1. Goal

Close the ε_r=10.2 high-contrast inhomogeneous residual (β=489 vs
reference 583, ~16%) that step 5.3 proved is **first-order-element
discretization** of the dielectric-interface field peak — **cheaply,
via interface-graded `h`-refinement**, before resorting to (expensive)
higher-order `p`-refinement.

Target: ε_r=10.2 horizontal-slab β within **≤5%** of the verified
`reference.rs` at a feasible DoF count, leveraging the step-5.3 sparse
solve. No gate weakened; FR-4 + uniform + homogeneous stay green.

## 2. Background — why h-refinement first

step 5.3's mesh-convergence study showed β plateaus under *uniform*
refinement (489→487→486) far short of 583. The diagnosed cause (ADR-0054
as-built) is that first-order Nedelec/nodal elements under-resolve the
**field peak at the high-contrast dielectric interface** — the dominant
mode concentrates energy in the ε_r=10.2 layer with a sharp gradient at
the interface that uniform coarse elements smear.

The textbook responses are `h`-refinement (smaller elements where the
field varies fastest — i.e. *graded* toward the interface) and
`p`-refinement (higher polynomial order). `p`-refinement on curl-
conforming elements is a large implementation (second-order Nedelec, the
curl no longer constant per triangle → embedded Gauss quadrature, new
DoF bookkeeping). **`h`-refinement graded to the interface is far
cheaper** — it reuses the existing first-order element matrices and the
step-5.3 sparse solve, needing only a non-uniform mesh generator that
concentrates rows near `y = d₁`. Walking-skeleton discipline: try the
cheap, high-leverage fix first; only escalate to `p` if `h` plateaus.

## 3. Approach

`crates/yee-mom/src/eigensolver/mesh.rs` (+ the test mesh builders):

1. **Graded horizontal-slab mesh.** Add a `TriMesh2D` builder that
   places `y`-grid lines geometrically clustered toward the interface
   `y = d₁` (e.g. a symmetric geometric grading with ratio `r`, finest
   cell at the interface), keeping the material tag assignment
   (dielectric below, air above) exact at the interface. Keep `x`
   uniform (the field varies slowly in `x` for the dominant mode).
2. **Convergence study.** With the step-5.3 sparse solve, sweep grading
   strength + total DoF; report β → reference. Find the DoF/grading that
   reaches ≤5% (if achievable with first-order elements).
3. **Disposition:**
   - **If graded `h` reaches ≤5%:** make the ε_r=10.2 reconciliation a
     failing ≤5% gate (alongside the FR-4 gate). §4 closed at high
     contrast too. Document the grading + DoF needed.
   - **If graded `h` plateaus short** (first-order convergence rate too
     slow even graded): document the residual + the achieved
     improvement, keep the ε_r=10.2 reconciliation a non-failing
     diagnostic, and queue **step-5.5** (`p`-refinement / second-order
     Nedelec) with the evidence that `h` alone is insufficient.
4. **Optional (timing, P2 from step 5.3):** the cutoff-pencil mode
   selection still does a dense `O(n³)` `complex_eigenvalues`, capping
   practical `cargo test` meshes at ~12×12. If the graded mesh needs
   more DoF than that allows under the dense selection, make the
   selection sparse too (shift-invert on the cutoff pencil) — or run the
   high-DoF convergence as a `--release` example rather than an inline
   test. Prefer the latter (cheaper) unless the gate needs the high-DoF
   point inline.

## 4. Validation / DoD

- DoD-1. Graded-mesh builder + a convergence study (β vs reference at
  ≥3 grading/DoF points) for ε_r=10.2.
- DoD-2. Either: ε_r=10.2 β ≤5% vs reference → failing gate (§4 closed
  at high contrast); OR a documented plateau + step-5.5 queued, with the
  best achieved β recorded.
- DoD-3. **No regression:** FR-4 gate (≤5%), uniform-fill anchor (≤1%),
  ε_r=1 canary (4e-14-class), `eigensolver_wr90`, coupling guards, Z_w
  all stay green. The graded mesh is *additive* — the uniform-mesh
  cases keep their builders/values.
- DoD-4. No new `Cargo.toml` dependency. Lint floor clean.
  `reference.rs` untouched.

## 5. Risks

(a) **First-order convergence rate too slow even graded.** The field
singularity at a high-contrast dielectric corner/edge can need `p` or
very aggressive grading. Mitigation: this is the explicit DoD-2 "plateau
→ step-5.5" branch; `h`-refinement is the cheap *attempt*, not a
guaranteed close. Even a partial improvement (e.g. 16%→8%) is progress +
evidence for the `p` decision.
(b) **DoF blow-up under the dense cutoff selection.** §3.4 — prefer a
release-example convergence study over an inline high-DoF test; only
make the selection sparse if the inline gate demands it.
(c) **Grading distorts the interface tag boundary.** Mitigation: place a
grid line *exactly* at `y = d₁` so the material partition stays sharp;
grade the cells on either side, not across the interface.

## 6. References

* Jin, *FEM in EM* 3rd ed. §9 (h- vs p-refinement, edge elements).
* Babuška & Suri, "The p and h-p versions of the FEM" (convergence rates).
* ADR-0054 as-built (the discretization diagnosis + the mesh-convergence
  plateau evidence).
* `crates/yee-mom/src/eigensolver/{mesh,solve,reference}.rs`.
