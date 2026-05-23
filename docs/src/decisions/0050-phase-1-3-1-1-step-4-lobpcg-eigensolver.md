# ADR-0050 — Phase 1.3.1.1 step 4: in-tree block LOBPCG over arpack-rs

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 4 (sparse eigensolver upgrade)

## Context

ROADMAP "Pending (high priority)" names step 4 as "sparse arpack-rs /
LOBPCG eigensolver". Track NNNNNN (`fb6be04`) already shipped the
`SparseEigen` / `SparseEigenComplex` traits plus `InverseIterEigen` /
`ComplexInverseIterEigen` — deflated shift-invert inverse-power
iteration over a faer sparse LU of `(K − σM)`, computing eigenpairs
**one at a time**. That escape-hatch meets `fem-eig-001` but is weak on
the clustered / degenerate spectra that real waveguide cross-sections
exhibit (TE/TM degeneracies, `TE_{mn}`/`TE_{nm}` pairs): sequential
deflation converges slowly and accumulates M-orthogonality error across
a cluster.

The question is **how to get a robust block eigensolver**:

* **Option A** — bind system ARPACK via `arpack-rs`. The canonical
  Krylov–Schur shift-invert sparse eigensolver.
* **Option B** — adopt a published Rust `lobpcg` crate.
* **Option C** — implement an **in-tree pure-Rust block LOBPCG**
  (Knyazev 2001) on top of faer + the existing sparse LU, behind the
  same `SparseEigen` trait.

## Decision

**Option C**: in-tree pure-Rust block LOBPCG (`LobpcgEigen:
SparseEigen`). Existing `InverseIterEigen` is retained; LOBPCG is an
additive impl the caller selects.

## Rationale

(1) **No external toolchain.** `arpack-rs` (Option A) binds
Fortran/LAPACK. Making it a default dependency violates CLAUDE.md §3
("feature flags default OFF for anything requiring an external
toolchain") and risks the CI lint-clean / `manylinux` wheel build.
Option C adds **zero** new `Cargo.toml` dependencies — it reuses faer
(already linked) for the small dense Rayleigh-Ritz subproblem and the
sparse LU already built by `build_shifted`.

(2) **The escape-hatch already proved the trait seam.** Track NNNNNN's
`solve.rs:147` note ("the eventual LOBPCG / ARPACK swap is one PR") and
the cross-section eigensolver spec §fallback both anticipated this.
The `lobpcg` crate (Option B) was unavailable/unusable at NNNNNN time;
nothing has changed that. An in-tree implementation is ~200 lines on
primitives the crate already has.

(3) **Pure-Rust LA ethos.** Consistent with faer-as-reference and the
deliberate avoidance of `nalgebra-lapack` / `ndarray-linalg` system
bindings elsewhere in the workspace.

(4) **The trait stays the swap point.** If a future >10⁵-DoF
cross-section ever needs ARPACK's Krylov–Schur, it can be added as an
*optional `arpack` feature* implementing the same `SparseEigen` trait —
this decision does not foreclose that, it defers it until a workload
demands it.

## Consequences

* New `LobpcgEigen` type in `crates/yee-fem/src/solve.rs`; no new
  crate dependency.
* `InverseIterEigen` remains the default for existing consumers
  (`fem-eig-001`, `eigensolver_wr90`); pointing those at `LobpcgEigen`
  is a consumer-lane follow-on, out of scope for step 4 (this phase
  proves parity + cluster robustness via `yee-fem` self-tests).
* Complex arm (`ComplexLobpcgEigen`) deferred to step 4.1; lossy
  dispersive cavities (`fem-eig-002`) keep `ComplexInverseIterEigen`.
* A new degenerate-cluster unit test pins the capability LOBPCG adds
  over inverse iteration.

## References

* Knyazev, "Toward the Optimal Preconditioned Eigensolver: LOBPCG",
  SIAM J. Sci. Comput. 23(2), 2001.
* Phase 1.3.1.1 step 4 spec
  `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-4-lobpcg-eigensolver-design.md`
  + plan `docs/superpowers/plans/2026-05-23-phase-1-3-1-1-step-4-lobpcg-eigensolver.md`.
* ADR-0022 — Phase 1.3.1.1 eigensolver spec deferred-impl.
* `crates/yee-fem/src/solve.rs` — `SparseEigen`, `InverseIterEigen`.
