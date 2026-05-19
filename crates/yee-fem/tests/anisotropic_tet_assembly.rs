//! Phase 4.fem.eig.3.5 P3 — unit tests for
//! [`yee_fem::assemble_tet_element_complex_anisotropic`].
//!
//! The anisotropic per-tet path drives the CFS-PML volumetric
//! absorber (Roden–Gedney 2000); it must:
//!
//! 1. Reduce bit-for-bit to the scalar
//!    [`yee_fem::assemble_tet_element_complex`] entry point when the
//!    input tensor is a scalar multiple of the identity (gate
//!    `scalar_equivalence_when_tensor_is_scalar_times_identity`).
//! 2. Produce complex-symmetric mass blocks for a diagonal complex
//!    `ε_tensor` (gate `diagonal_anisotropic_block_is_complex_symmetric`).
//! 3. Reject off-diagonal tensors with `Error::Unimplemented` —
//!    rotated PML is deferred to Phase 4.fem.eig.3.5.1 per ADR-0043
//!    §4 (gate `off_diagonal_tensor_rejected_until_v3_5_1`).
//!
//! References:
//!
//! * Phase 4.fem.eig.3.5 spec
//!   `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
//!   §4.3.
//! * Phase 4.fem.eig.3.5 plan
//!   `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-cfs-pml.md`
//!   step P3.
//! * ADR-0043 — Cartesian-aligned-only scope decision.

#![allow(non_snake_case)]

use nalgebra::{SMatrix, Vector3};
use num_complex::Complex64;
use yee_fem::{assemble_tet_element_complex, assemble_tet_element_complex_anisotropic};

/// A canonical unit-corner tet, vertices ordered for positive signed
/// volume (`yee_mesh::TetMesh3D::new` enforces this for production
/// meshes; the test uses the analytic vertex layout directly).
fn canonical_unit_tet() -> [Vector3<f64>; 4] {
    [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ]
}

/// Build a diagonal `3 × 3` complex tensor.
fn diag(a: Complex64, b: Complex64, c: Complex64) -> SMatrix<Complex64, 3, 3> {
    let mut m = SMatrix::<Complex64, 3, 3>::zeros();
    m[(0, 0)] = a;
    m[(1, 1)] = b;
    m[(2, 2)] = c;
    m
}

#[test]
fn scalar_equivalence_when_tensor_is_scalar_times_identity() {
    let vertices = canonical_unit_tet();
    let eps = Complex64::new(2.5, 0.3);
    let mu = Complex64::new(1.1, -0.07);

    let scalar = assemble_tet_element_complex(vertices, eps, mu);
    let aniso = assemble_tet_element_complex_anisotropic(
        vertices,
        diag(eps, eps, eps),
        diag(
            Complex64::new(1.0, 0.0) / mu,
            Complex64::new(1.0, 0.0) / mu,
            Complex64::new(1.0, 0.0) / mu,
        ),
    )
    .expect("anisotropic helper rejected diag tensor");

    let k_diff = (scalar.k_local - aniso.k_local).norm();
    let m_diff = (scalar.m_local - aniso.m_local).norm();
    assert!(
        k_diff < 1.0e-12,
        "scalar / anisotropic stiffness blocks must agree to round-off; got |ΔK| = {k_diff:e}"
    );
    assert!(
        m_diff < 1.0e-12,
        "scalar / anisotropic mass blocks must agree to round-off; got |ΔM| = {m_diff:e}"
    );
}

#[test]
fn diagonal_anisotropic_block_is_complex_symmetric() {
    let vertices = canonical_unit_tet();
    let eps_tensor = diag(
        Complex64::new(2.0, 0.1),
        Complex64::new(1.5, 0.05),
        Complex64::new(3.0, 0.2),
    );
    let mu_inv_tensor = diag(
        Complex64::new(0.9, -0.01),
        Complex64::new(1.0, 0.0),
        Complex64::new(1.05, 0.02),
    );

    let aniso = assemble_tet_element_complex_anisotropic(vertices, eps_tensor, mu_inv_tensor)
        .expect("anisotropic helper rejected diag tensor");

    // Complex-symmetric, NOT Hermitian — same convention as the scalar
    // entry point (ADR-0039 §6).
    let k_asym = (aniso.k_local - aniso.k_local.transpose()).norm();
    let m_asym = (aniso.m_local - aniso.m_local.transpose()).norm();
    assert!(
        k_asym < 1.0e-12,
        "diagonal anisotropic K must be complex-symmetric (K == K^T); got |K - K^T| = {k_asym:e}"
    );
    assert!(
        m_asym < 1.0e-12,
        "diagonal anisotropic M must be complex-symmetric (M == M^T); got |M - M^T| = {m_asym:e}"
    );
}

#[test]
fn off_diagonal_tensor_rejected_until_v3_5_1() {
    let vertices = canonical_unit_tet();
    let mut eps_tensor = diag(
        Complex64::new(2.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(2.0, 0.0),
    );
    // Inject an off-diagonal entry — represents a rotated PML axis
    // outside the Cartesian-aligned v3.5 envelope (ADR-0043 §4).
    eps_tensor[(0, 1)] = Complex64::new(0.1, 0.0);
    eps_tensor[(1, 0)] = Complex64::new(0.1, 0.0);

    let mu_inv_tensor = diag(
        Complex64::new(1.0, 0.0),
        Complex64::new(1.0, 0.0),
        Complex64::new(1.0, 0.0),
    );

    let result = assemble_tet_element_complex_anisotropic(vertices, eps_tensor, mu_inv_tensor);
    assert!(
        result.is_err(),
        "off-diagonal ε_tensor must produce Error::Unimplemented; got {result:?}"
    );
    let err_str = format!("{}", result.unwrap_err());
    assert!(
        err_str.contains("off-diagonal") || err_str.contains("unimplemented"),
        "error message should call out the off-diagonal restriction; got: {err_str}"
    );
}
