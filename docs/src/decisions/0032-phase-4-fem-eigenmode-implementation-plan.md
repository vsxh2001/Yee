# ADR-0032: Phase 4.fem.eig.0 FEM eigenmode implementation plan

## Status

Accepted — 2026-05-18 (plan only; track execution deferred to
follow-up agents).

## Context

ADR-0029 locked the Phase 4.fem.eig.0 scope. The next step is
~2 400 LOC of plan: new `TetMesh3D`, the `yee-fem` element /
assembly / solve modules, the Kuhn-decomposed cavity mesh, the
gate. Track SSSSS (merge `0d596ab`) lands the plan. The sparse
generalized eigensolve library is the load-bearing risk; this
ADR records the `SparseEigen` trait abstraction that isolates it.

## Decision

The plan splits into **nine tracks** (seven core + two optional):

- **T1 ‖ T2 parallel start.** T1: `yee-fem` workspace member,
  stub surface returning `Error::Unimplemented`. T2: `TetMesh3D`
  in `yee-mesh` with positive-volume validation and a boundary-
  edge classifier.
- **T3 (on T1) — local 6-edge Nedelec element matrices.**
  Barycentric gradients, `N_{ij} = λ_i ∇λ_j − λ_j ∇λ_i`,
  constant curl per tet, 4-point Gauss quadrature; unit-tested
  against Jin Ch. 9 tabulated reference.
- **T4 (on T2+T3) — global sparse assembly + PEC elimination.**
  Orientation-aware scatter into `nalgebra_sparse::CsrMatrix<f64>`.
  Pattern: `crates/yee-mom/src/eigensolver/assembly.rs`.
- **T5 (on T1, sibling of T3/T4) — `SparseEigen` trait +
  `LobpcgEigen` shift-invert** over `lobpcg` + `faer` sparse LU;
  fallback is hand-rolled deflated inverse-iter behind the trait.
- **T6 (on T2) — `cavity_uniform` Kuhn-decomposed mesh.** Six
  tets per brick, orientation-preserving, no slivers.
- **T7 (on T4+T5+T6) — fem-eig-001 production gate.** WR-90
  cavity on (8,6,10) mesh; TE₁₀₁ ±0.3%, ten-mode RMS ±1%, no
  spurious below TE₁₀₁.
- **T8 ‖ T9 (optional) post-T7.** Python binding + tutorial.

Four load-bearing decisions locked here:

- **`SparseEigen` trait isolates the library to one PR.** Same
  swap-point discipline as `yee_cuda::backend::Backend`.
- **fem-eig-001's ±0.3% gate is **not** weakened.** Escape hatch
  refines (8,6,10) → (12,9,15); further failure is investigated
  upstream (gradient capture, T4 sign, Jin reference).
- **`arpack-rs` stays non-default.** Same CI-lint-policy escape
  hatch as Phase 1.3.1.1.
- **Hand-rolled Kuhn cavity is the only v0 mesh source.** Gmsh
  `.msh` tet ingestion is 4.fem.eig.0.2.

Critical path: `T1 → T3 → T4 → T7` and `T2 → T6 → T7` and
`T1 → T5 → T7`. Peak active set is three — within CLAUDE.md §5's
envelope.

## Consequences

- **Sparse-eigen library risk is isolated.** Trait is the seam;
  a library swap doesn't perturb T1/T3/T4/T6/T7.
- **mom-001 and shipped solvers stay green.** `yee-fem` is a new
  peer crate, not a refactor.
- **Hand-rolled Kuhn mesh keeps the gate accuracy floor
  independent of Gmsh.**
- **T8 / T9 are optional.** PyO3 abi3-py310 ABI mismatch with
  `lobpcg` defers T8 to 4.fem.eig.0.1 as a finding.
- **fem-eig-002 (lossy Q) and fem-eig-003 (DRA) remain open.**

## References

- `docs/superpowers/plans/2026-05-18-phase-4-fem-eigenmode.md`
- `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
- Track SSSSS merge commit `0d596ab`.
- ADR-0029 — Phase 4.fem.eig.0 scope lock (this plan's parent).
- ADR-0023 — `TriMesh2D` precedent; `TetMesh3D` mirrors invariants.
- J.-M. Jin, *FEM in Electromagnetics*, 3rd ed., Wiley 2014, §9.4.
- D. M. Pozar, *Microwave Engineering*, 4th ed., 2012, §6.3.
- H. W. Kuhn, *IBM J. Res. Dev.*, 1960 — T6 brick decomposition.
- CLAUDE.md §3, §4, §5, §6.
