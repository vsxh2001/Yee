//! Global FEM assembly — placeholder.
//!
//! This module will own:
//!
//! * Local 6-edge Nedelec tet element matrices (`K_local`, `M_local`) —
//!   delivered by step **T3** of the SSSSS plan
//!   (`docs/superpowers/plans/2026-05-18-phase-4-fem-eigenmode.md`).
//! * Global edge enumeration via [`yee_mesh::TetMesh3D`], orientation-aware
//!   scatter of the local blocks into `nalgebra_sparse::CsrMatrix<f64>`, and
//!   PEC tangential-`E`-zero Dirichlet elimination by row/column drop —
//!   delivered by step **T4** of the same plan.
//!
//! The Phase 4 design spec at
//! `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
//! locks the public surface (`FemEigenAssembly<'m>`, `AssembledMatrices`).
//! Nothing is exposed here yet — this scaffold module exists so the crate
//! compiles end-to-end while T3/T4 are in flight.
