//! 2-D FEM eigensolver for wave-port cross-section modal analysis.
//!
//! **Phase 1.3.1.1 step 2 (this commit):** Nedelec edge-element +
//! nodal-Lagrange `E_z` element-matrix assembly on a
//! [`yee_mesh::TriMesh2D`]. The dense generalized-eigensolve fallback
//! that consumes the assembled matrices lands as a follow-up commit
//! (step 3) under a `solve` submodule. Used internally by
//! [`crate::ports::NumericalCrossSection::solve`] to extract the
//! dominant propagation constant `β` and wave impedance `Z_w` of an
//! arbitrary cross-section at a single frequency.
//!
//! The assembly follows Jin, *The Finite Element Method in
//! Electromagnetics*, 3rd ed., §8.5 ("Modes of a Waveguide"):
//!
//! * `E_t` is expanded in first-order curl-conforming Nedelec (Whitney-1)
//!   edge basis functions, one DoF per interior edge.
//! * `E_z` is expanded in linear nodal-Lagrange (Whitney-0) basis
//!   functions, one DoF per interior vertex.
//! * PEC walls impose homogeneous Dirichlet on tangential `E_t` (drop
//!   boundary-edge DoFs) and on `E_z` (drop boundary-vertex DoFs).
//!
//! The result is a real-symmetric (lossless case) generalized
//! eigenproblem `S x = k_c² T x` whose smallest strictly-positive
//! eigenvalue gives the dominant mode's cutoff; the propagation
//! constant follows from `β² = k_0² − k_c²`. Sparse shift-and-invert
//! is Phase 1.3.1.1 step 4 (escape-hatched away from `arpack-rs`).

// The assembly path lands in this commit but is consumed by
// `NumericalCrossSection::solve` only after the Phase 1.3.1.1 step 3
// (dense eigensolve) commit follows. Until then the assembly items
// are unreferenced outside their own unit tests, so we allow dead-code
// at the module root rather than tagging every helper individually.
#![allow(dead_code)]

pub(crate) mod assembly;
pub(crate) mod mesh;
