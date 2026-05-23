//! 2-D FEM eigensolver for wave-port cross-section modal analysis.
//!
//! **Phase 1.3.1.1 (this lane).** Nedelec edge-element + nodal-Lagrange
//! `E_z` element-matrix assembly on a [`yee_mesh::TriMesh2D`] (step 2),
//! feeding a mixed `(E_t, E_z)` β-direct generalized eigensolve (steps
//! 5–5.3). Used internally by
//! [`crate::ports::NumericalCrossSection::solve`] to extract the dominant
//! propagation constant `β` and wave impedance `Z_w` of an arbitrary
//! cross-section at a single frequency.
//!
//! The **production** solve ([`solve::solve_dense_mixed`], step 5.3) is a
//! **faer sparse shift-and-invert** of the β-direct pencil
//! `(k_0² B − A) x = β² B_1 x` at a physics-informed shift, recovering the
//! true β-direct eigenvector (so β² is exact for the mode). The earlier
//! dense [`nalgebra::SymmetricEigen`] transverse-only path
//! ([`solve::solve_dense`]) and the step-5.2 dense Rayleigh-quotient hybrid
//! ([`solve::solve_dense_mixed_rq`]) are retained as small-`n` references /
//! regression anchors.
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

pub(crate) mod assembly;
pub(crate) mod mesh;
pub(crate) mod reference;
pub(crate) mod solve;

pub(crate) use solve::solve_dense_mixed;

/// Finite-element polynomial order for the cross-section mixed
/// `(E_t, E_z)` assembly.
///
/// **First-order is the default** (Whitney-1 Nedelec for `E_t`, linear
/// nodal-Lagrange for `E_z`) — the path validated by the WR-90 TE10 gate,
/// the FR-4 inhomogeneous closure, the uniform-fill analytic anchor, and
/// the ε_r=1 homogeneous canary. It stays bit-for-bit unchanged.
///
/// [`ElementOrder::Second`] selects the **second-order** family (Phase
/// 1.3.1.1 step 5.5): curl-conforming Nedelec-first-kind order-2 for `E_t`
/// (2 DoF/edge + 2 interior; Jin §9.4 / Webb hierarchal vector bases) and
/// quadratic nodal-Lagrange for `E_z` (6 nodes/triangle). It is selected
/// only for the high-contrast inhomogeneous case (ε_r=10.2), where the
/// first-order convergence rate plateaus short of the verified reference
/// (step 5.4). The curl is non-constant per triangle at second order, so
/// every element integral goes through a triangle Gauss rule (see
/// [`assembly::tri_gauss_deg4`]).
///
/// The selector is exercised by the lib-internal step-5.5 J3/J4 anchors
/// (which call the order-specific assemblers directly). The production
/// [`crate::ports::NumericalCrossSection::solve`] path (out of the step-5.5
/// lane) stays on the first-order [`assembly::assemble_mixed`]; wiring the
/// selector through that public boundary is a follow-up — hence the
/// `dead_code` allow on the non-test lib build.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub(crate) enum ElementOrder {
    /// Whitney-1 Nedelec + linear nodal-Lagrange (1 DoF/edge, 1 DoF/vertex).
    /// The default, validated path.
    #[default]
    First,
    /// Nedelec-first-kind order-2 + quadratic nodal-Lagrange (2 DoF/edge +
    /// 2 interior for `E_t`; vertex + edge-midpoint nodes for `E_z`).
    Second,
}
