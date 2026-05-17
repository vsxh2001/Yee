# ADR-0029: Phase 4.fem.eig.0 3-D FEM eigenmode walking-skeleton scope

## Status

Accepted — 2026-05-18 (spec only; implementation deferred to
follow-up tracks — see ADR-0032).

## Context

The shipped solver portfolio leaves a real gap on resonant 3-D
problems. FDTD ringdown fails on lossless cavities; planar MoM
has no volumetric interior; the 2-D eigensolver from Phase
1.3.1.1 finds guided-wave modes, not 3-D resonant modes.
`ROADMAP.md` Phase 4 calls out cavity-filter design, Q-factor
extraction, DRAs, and accelerator cavities as targets that need a
3-D FEM eigensolver outright. The physics is textbook (Jin
Ch. 9–10); the engineering risk is sparse generalized eigensolve
plumbing in pure Rust. Track OOOOO (merge `6f75b5a`) lands the
spec.

## Decision

Phase 4.fem.eig.0 ships the minimum end-to-end pipeline:

1. **First-order Nedelec (Whitney-1) edge elements on tetrahedra,
   PEC cavity, lossless, real `ε_r` / `μ_r`.** Nodal Lagrange is
   rejected — its `∇×` kernel admits spurious modes; Nedelec's
   kernel is exactly the discrete gradient subspace, clusters at
   `k² = 0`, and shift-invert skips it. Higher-order (`p ≥ 2`)
   is 4.fem.eig.1.
2. **Generalized sparse eigensolve via shift-invert LOBPCG
   behind a `SparseEigen` trait.** Default `LobpcgEigen` over
   `lobpcg` + `faer` sparse LU preconditioner; `arpack-rs`
   non-default behind `eig-arpack` (same CI-lint-policy escape
   hatch as Phase 1.3.1.1); hand-rolled deflated inverse-power
   iteration is the documented fallback. Trait is the load-
   bearing decision; the library is the swap point.
3. **`TetMesh3D` in `yee-mesh`, paralleling `TriMesh2D`.** No
   FEM logic in `yee-mesh`; sources are hand-rolled Kuhn-
   decomposed cuboid fixtures and Gmsh `.msh`.
4. **New `crates/yee-fem/` peer crate**, not a child of
   `yee-mom::eigensolver` — non-MoM users should not pay FEM
   compile cost.

Validation gate **fem-eig-001**: WR-90-based rectangular cavity
(`a = 22.86 mm`, `b = 10.16 mm`, `d = 30 mm`). TE₁₀₁ at 9.660 GHz
within ±0.3%; lowest ten modes match Pozar §6.3 within ±1%
pairwise; no spurious mode below TE₁₀₁.

CPU-only, single-threaded, FP64. No GPU, no losses (4.fem.eig.2),
no DRA (4.fem.eig.3), no periodic BCs, no anisotropic /
dispersive media. Driven 3-D FEM and FEM-BEM hybrid are separate
Phase 4 sub-projects sharing the `K` assembly path.

## Consequences

- **Closes the resonant 3-D gap end-to-end.** Cavity filters,
  Q-factor (after 4.fem.eig.2's complex `ε_r`), DRA analysis
  (after 4.fem.eig.3) all reachable on the same surface.
- **Pure-Rust sparse generalized eigensolve is the load-bearing
  risk.** `lobpcg` is not battle-tested at FEM scale; the trait
  swap-point isolates the choice to one PR.
- **Hand-rolled Kuhn cavity is the v0 mesh source.** Gmsh tet
  ingestion is a 4.fem.eig.0.2 follow-up.
- **First-order accuracy floor is real.** Marginal gate → `p`-
  refinement (4.fem.eig.1), not `h`-refinement to absurd counts.

## References

- `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
- Track OOOOO merge commit `6f75b5a`.
- ADR-0032 — Phase 4.fem.eig.0 implementation plan (companion).
- J.-M. Jin, *The Finite Element Method in Electromagnetics*,
  3rd ed., Wiley 2014, Ch. 9–10.
- D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley 2012,
  §6.3 — fem-eig-001 reference.
- ADR-0022 — Phase 1.3.1.1 2-D Nedelec cross-section eigensolver
  spec; direct analog one dimension down.
- ADR-0023 — `TriMesh2D` precedent.
- CLAUDE.md §3, §4.
