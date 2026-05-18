//! # yee-fem
//!
//! Yee Finite Element Method (FEM) solver — the **Phase 4 eigenmode
//! beachhead**.
//!
//! ## Scope
//!
//! Phase 4.fem.eig.0 ships a single end-to-end pipeline: a hand-rolled
//! tetrahedral mesh of a closed rectangular metallic air-filled cavity is
//! consumed by this crate, which assembles first-order Nedelec (Whitney-1)
//! curl-curl stiffness `K` and vector mass `M` sparse matrices, applies the
//! PEC tangential-`E`-zero Dirichlet condition by row/column elimination,
//! and solves the generalised eigenproblem `K e = k² M e` via shift-invert
//! LOBPCG for the ten smallest positive eigenvalues. Validation gate
//! `fem-eig-001` enforces TE_{101} resonance at 9.660 GHz within ±0.3% on
//! a WR-90-based cavity (a = 22.86 mm, b = 10.16 mm, d = 30 mm) against
//! the Pozar §6.3 analytic table.
//!
//! ## Status
//!
//! **Phase 4.fem.eig.0 step T1 scaffold.** This crate currently exports only
//! the module skeleton and no public solver API. The element matrices,
//! global assembly, and sparse eigen-solver land in steps T3, T4, and T5
//! respectively per the SSSSS implementation plan
//! (`docs/superpowers/plans/2026-05-18-phase-4-fem-eigenmode.md`).
//!
//! See the companion design spec at
//! `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
//! and ADRs 0029 / 0032 in `docs/src/decisions/` for the scope decisions
//! that gate every module below.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod assembly;
pub mod element;
pub mod solve;

pub use assembly::{AssembledMatrices, FemEigenAssembly};
pub use solve::{EigenpairList, InverseIterEigen, SparseEigen};
