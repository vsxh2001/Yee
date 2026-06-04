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
//!
//! ## Eigensolver options (Phase 1.3.1.1 step 4)
//!
//! [`solve`] exposes the [`solve::SparseEigen`] trait with two
//! real-coefficient implementations: the default
//! [`solve::InverseIterEigen`] (deflated shift-invert inverse-power
//! iteration, one mode at a time) and [`solve::LobpcgEigen`] (in-tree
//! block LOBPCG, Knyazev 2001), which resolves clustered / degenerate
//! cross-section spectra a sequential deflation handles poorly.
//! `LobpcgEigen` adds no new dependency — its small dense Rayleigh-Ritz
//! step reuses `nalgebra` and its preconditioner is the same faer
//! sparse LU of `(K − σM)` (ADR-0050; the `arpack` Krylov–Schur path is
//! deferred behind the same trait). The complex-coefficient
//! [`solve::ComplexInverseIterEigen`] keeps the lossy dispersive
//! `fem-eig-002` path; a `ComplexLobpcgEigen` peer is a step-4.1
//! follow-on.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod assembly;
pub mod coupled_resonator_k;
pub mod dispersive;
pub mod element;
pub mod material;
pub mod microstrip_mesh;
pub mod microstrip_port;
pub mod microstrip_port_numerical;
pub mod open_boundary;
pub mod pml_mesh;
pub mod solve;

pub use assembly::{AssembledMatrices, AssembledMatricesComplex, FemEigenAssembly};
pub use coupled_resonator_k::{
    CoupledKResult, CoupledResonatorGeom, GapCorrection, correct_gap_fem_k, coupled_resonator_k,
};
pub use dispersive::{BULK_TAG, DispersiveEigenpair, DispersiveError, DispersiveSolver};
pub use element::{
    LOCAL_EDGES, NedelecTetElement, NedelecTetElementComplex, assemble_abc_face_block,
    assemble_abc2_face_block, assemble_port_face_block, assemble_port_face_block_gauss_pts,
    assemble_port_face_rhs_gauss_pts, assemble_port_modal_rhs, assemble_tet_element,
    assemble_tet_element_complex, assemble_tet_element_complex_anisotropic,
};
pub use material::{Material, MaterialDatabase, MaterialPole};
pub use microstrip_mesh::{
    AIR_TAG, FR4_EPS_R, FR4_TAG, TraceRect, layered_microstrip_filter_mesh, layered_microstrip_mesh,
};
pub use microstrip_port::{
    beta_microstrip, microstrip_port, microstrip_port_windowed, modal_e_t_microstrip,
    modal_e_t_microstrip_windowed,
};
pub use microstrip_port_numerical::{
    MicrostripPortGeom, microstrip_port_numerical, microstrip_port_numerical_at,
};
pub use open_boundary::{
    AbcOrder, DrivenSystem, FaceKind, OpenBoundarySolver, PmlConfig, PmlMeshMeta, PmlRegion,
    PortDefinition, PortId, PortMode, SParameters, SParametersMatrix,
};
pub use pml_mesh::{FaceIndexMap, PmlAxis, PmlClass, extend_mesh_with_pml};
pub use solve::{
    ComplexInverseIterEigen, ComplexLobpcgEigen, EigenpairList, EigenpairListComplex,
    InverseIterEigen, LobpcgEigen, SparseEigen, SparseEigenComplex,
};
