# Phase 1.3.1.1 step 5.5 — second-order Nedelec/nodal elements for high-contrast accuracy

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.5 (close the high-contrast inhomogeneous gate).
**Depends on:** step 5.3 (direct sparse β-direct solve; `6eca76a`),
step 5.4 (graded-h plateau — proved h cannot fix it; `95ac64f`).
**Blocks:** tightening the §4 inhomogeneous gate from FR-4 (closed) to
high-contrast (ε_r=10.2).

## 0. Why this, and why now (decision note)

ADR-0055 deprioritised p-refinement as "marginal-case." On
re-evaluation at the cross-section-chain milestone, every *breadth*
alternative (mom-002 numerical-port, FDTD Q6/Q7, FEM real-waveguide-port)
is **grind-risky and poorly dispatchable** — mom-002 is a 10-track
forensic quagmire entangled with the `MultilayerGreens` placeholder and
a subtle cross-formulation coupling; FDTD Q6 is an open energy-balance
problem; the real-port is entangled with the deprioritised fem-eig-006
modal-projection grind. step-5.5 is the one **high-confidence,
cleanly-dispatchable, won't-grind** option: the fix is textbook (Jin §9,
Webb 1993), the DoD is crisp (ε_r=10.2 ≤5% vs the verified reference),
step 5.4 *proved* it is the right lever (h ruled out), and the
interface-aligned mesh guarantees fast p-convergence (no corner
singularity — the slab interface is a straight, element-aligned line, so
each element is within one smooth material). A marginal-but-certain win
beats a high-variance quagmire gamble in an autonomous loop.

## 1. Goal

Add **second-order (p=2) curl-conforming Nedelec edge elements** (for
`E_t`) and **quadratic nodal-Lagrange elements** (for `E_z`) to the
cross-section eigensolver, and close the ε_r=10.2 high-contrast
inhomogeneous gate to ≤5% vs `reference.rs`. The first-order path stays
the default; p=2 is selected for high-contrast inhomogeneous solves.

No gate weakened; FR-4 / uniform / homogeneous stay green.

## 2. Background

step 5.4 proved the ε_r=10.2 residual (~16%, β=489 vs reference 583) is
the **first-order-element convergence rate** at the interface field peak
— graded `h` plateaus at ~15.6%. Second-order Nedelec elements raise the
energy convergence from `O(h)` to `O(h²)` (eigenvalue `O(h⁴)`), which —
on an interface-aligned mesh (node exactly at `y = d₁`, so each element
is within one material and the solution is smooth per-element) —
converges to the reference far faster. This is the standard `p`-fix; step
5.4's plateau is the evidence it is warranted.

## 3. Approach

`crates/yee-mom/src/eigensolver/assembly.rs` (+ `mesh.rs` DoF map):

1. **Second-order Nedelec (E_t).** Whitney-1 (first-order) has 1 DoF/edge;
   second-order Nedelec on a triangle has 2 DoF/edge + 2 interior DoF
   (8 total/triangle, Jin §9.4 / Webb 1993 "hierarchal" edge elements).
   The curl is **no longer constant per triangle** → the
   `∫(1/μ_r)(∇×N)(∇×N)` and `∫ε_r N·N` element integrals need a Gauss
   quadrature rule exact for the integrand degree (the existing
   first-order code uses closed-form constants; add a 2-D triangle Gauss
   rule, e.g. 6-point degree-4).
2. **Quadratic nodal-Lagrange (E_z).** 6 nodes/triangle (3 vertex + 3
   edge-midpoint). Standard `∫(1/μ_r)∇L·∇L`, `∫ε_r L·L`, and the
   `∫(1/μ_r)∇L·N` coupling at quadratic order, via the same quadrature.
3. **DoF bookkeeping** (`mesh.rs` / a higher-order DoF map): edge DoFs
   (2/edge for E_t, 1 midpoint-node/edge for E_z) + interior. PEC
   Dirichlet elimination extends to the new DoFs.
4. **Solve.** Reuse the step-5.3 direct sparse β-direct shift-invert
   (`solve_dense_mixed`) on the larger pencil — the assembly produces the
   same `(A, B, B₁)` block structure, just higher-dimensional. Mode
   selection + the β-direct extraction are order-agnostic.
5. **Order selection.** A `p: u8` (or `ElementOrder` enum) parameter on
   the assembly; first-order stays the default (homogeneous/FR-4/uniform
   gates unchanged); p=2 selected for the ε_r=10.2 case.

## 4. Validation / DoD

- DoD-1. p=2 assembly + a convergence study: ε_r=10.2 β → reference
  582.95 at ≥3 mesh densities, showing the faster (p=2) convergence vs
  the step-5.4 first-order plateau.
- DoD-2. **§4 high-contrast closure:** ε_r=10.2 horizontal-slab β within
  **≤5%** of `reference.rs` → a failing gate (alongside FR-4). If p=2
  still falls short (unexpected — would indicate a remaining issue, e.g.
  interface mis-alignment or a quadrature/DoF bug), document + surface as
  a finding, do NOT weaken.
- DoD-3. **No regression:** first-order path bit-identical (FR-4 gate,
  uniform anchor, ε_r=1 canary, `eigensolver_wr90`, coupling guards, Z_w
  all green — p=2 is additive, selected only for the high-contrast case).
- DoD-4. **p=2 sanity on a known case:** p=2 on the homogeneous WR-90
  reproduces the analytic TE10 β at least as accurately as p=1 (a clean
  no-singularity check that the p=2 element matrices are correct).
- DoD-5. No new `Cargo.toml` dependency. Lint floor clean.
  `reference.rs` untouched.

## 5. Risks

(a) **Element-matrix bug** (the curl is non-constant — quadrature must be
exact for the integrand). Mitigation: DoD-4 (p=2 reproduces analytic
TE10 on the homogeneous guide — a sharp correctness check independent of
the high-contrast case); a unit test comparing a p=2 element matrix to a
hand/independent-quadrature value (mirror the step-5 `local_b_ze`
independent-quadrature test).
(b) **DoF bookkeeping** (orientation/sign of the 2nd edge DoF). Mitigation:
the homogeneous canary + DoD-4; a small-mesh symmetry check.
(c) **Cost.** p=2 ~doubles-to-triples DoF; the dense cutoff-pencil
selection (step-5.3) is `O(n³)`. Run the high-DoF convergence as a
`--release` example, keep the inline gate at a feasible mesh.
(d) **Still short of ≤5%.** If the interface field genuinely has a
weak singularity (it should not for a straight element-aligned
interface), p may also converge slowly — then it is a deeper modeling
question (step-5.6). Low probability given the geometry.

## 6. References

* Jin, *FEM in EM* 3rd ed. §9.4 (higher-order edge elements).
* Webb, "Hierarchal vector basis functions of arbitrary order for
  triangular and tetrahedral finite elements", IEEE TAP 1999 (and Webb
  1993). Ingerson/Savage hierarchal Nedelec.
* ADR-0054/0055 (the discretization diagnosis + the h-plateau evidence).
* `crates/yee-mom/src/eigensolver/{assembly,mesh,solve,reference}.rs`.
