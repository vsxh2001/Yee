//! Phase 4.fem.eig.2 step E2 — unit tests for
//! [`yee_fem::element::assemble_port_face_block`] and
//! [`yee_fem::element::assemble_port_modal_rhs`].
//!
//! Gate test inventory:
//!
//! 1. `port_face_block_complex_symmetric` — on the unit right triangle in
//!    the xy-plane with `n̂ = ẑ`, the returned 3×3 wave-port face block
//!    satisfies `B == B^T` to `1e-12` (complex-symmetric, NOT Hermitian).
//!    Mirrors the E1 ABC face-block invariant with `β_mode` replacing
//!    `k₀`.
//! 2. `port_face_block_imaginary_coefficient` — every entry of `B` is
//!    purely imaginary; `Re(B[i][j]) = 0` up to floating-point round-off
//!    because the prefactor `j · β_mode · A / μ_r` is purely imaginary and
//!    the `(n̂ × t_i) · (n̂ × t_j)` Gram form is real.
//! 3. `port_face_block_beta_zero_returns_zero` — when `β_mode = 0`
//!    (modal cutoff), the entire 3×3 block is identically zero
//!    (no special-case branch needed; the prefactor is zero).
//! 4. `port_modal_rhs_zero_e_t_returns_zero` — when the incident modal
//!    tangential E-field is zero, every entry of the RHS vector is
//!    identically zero.
//! 5. `port_modal_rhs_uniform_e_t_dot_normal` — with a uniform incident
//!    `E_t = x̂`, the RHS entry on each edge is proportional to the
//!    closed-form `t_i · E_t` projection times the face-centroid
//!    quadrature weight `2 j β · A / 3`.
//! 6. `port_modal_rhs_complex_symmetric_with_block` — at the modal-
//!    eigenvalue self-consistent point, the RHS vector lies in the range
//!    of the wave-port face block (equivalently, is orthogonal to the
//!    null space `span{(1, 1, 1)}` of the per-face Gram matrix — the
//!    closure of the three edge tangents `Σ t_i = 0`). This is the
//!    face-level necessary condition for `B · e_mode = b_mode / 2` to be
//!    solvable.
//!
//! References:
//! * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
//!   3rd ed., Wiley 2014, §10.5 (wave-port modal decomposition).
//! * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §3.3.
//! * Phase 4.fem.eig.2 spec
//!   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
//!   §4.3.

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::element::{assemble_port_face_block, assemble_port_modal_rhs};

/// Canonical unit right-triangle face in the xy-plane with outward
/// normal `+ẑ`. Vertices in CCW order seen from `+ẑ`.
fn unit_right_triangle() -> [Vector3<f64>; 3] {
    [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ]
}

/// Representative TE_{10} `β_mode` at 10 GHz on WR-90
/// (`a = 22.86 mm`): `β = sqrt(k₀² − (π/a)²)`.
/// Numerical value `≈ 158.0 rad/m`. The exact value does not matter for
/// these tests; what matters is that it is a non-zero real number of
/// O(100).
fn beta_te10_wr90_at_10ghz() -> f64 {
    let k0 = 2.0 * PI / 0.03; // λ ≈ 3 cm at 10 GHz in vacuum
    let kc = PI / 0.022_86; // cutoff wavenumber for TE_{10} on WR-90
    (k0 * k0 - kc * kc).sqrt()
}

/// E2 DoD criterion 1: the wave-port face block is complex-symmetric
/// (`B == B^T`), NOT Hermitian. Same invariant as the ABC face block
/// (E1) with `β_mode` replacing `k₀`.
#[test]
fn port_face_block_complex_symmetric() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let block = assemble_port_face_block(
        face_vertices,
        outward_normal,
        beta_te10_wr90_at_10ghz(),
        1.0,
    );

    let asymmetry = (block - block.transpose()).norm();
    assert!(
        asymmetry < 1e-12,
        "wave-port face block must be complex-symmetric (B == B^T); \
         got ||B - B^T|| = {asymmetry:e}"
    );
}

/// E2 DoD criterion 2: every entry of the wave-port face block is
/// purely imaginary. The prefactor `j · β_mode · A / μ_r` is purely
/// imaginary and the Gram form `(n̂ × t_i) · (n̂ × t_j)` is real.
#[test]
fn port_face_block_imaginary_coefficient() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let beta = beta_te10_wr90_at_10ghz();
    let block = assemble_port_face_block(face_vertices, outward_normal, beta, 1.0);

    // Tolerance scaled by the magnitude of the imaginary prefactor.
    let tol = f64::EPSILON * 100.0 * beta;
    for i in 0..3 {
        for j in 0..3 {
            let re = block[(i, j)].re;
            assert!(
                re.abs() < tol,
                "B[{i}][{j}] must be purely imaginary; got Re = {re:e} (tol = {tol:e})"
            );
        }
    }
}

/// E2 DoD criterion 3: at modal cutoff `β_mode = 0`, the wave-port
/// face block is identically zero. The prefactor `j · β_mode · A / μ_r`
/// vanishes; no special-case branch is needed.
#[test]
fn port_face_block_beta_zero_returns_zero() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let block = assemble_port_face_block(face_vertices, outward_normal, 0.0, 1.0);

    for i in 0..3 {
        for j in 0..3 {
            let value = block[(i, j)];
            assert_eq!(
                value.re, 0.0,
                "B[{i}][{j}].re must be exactly zero at cutoff; got {}",
                value.re
            );
            assert_eq!(
                value.im, 0.0,
                "B[{i}][{j}].im must be exactly zero at cutoff; got {}",
                value.im
            );
        }
    }
}

/// E2 DoD criterion 4: when the incident modal tangential E-field is
/// zero, the wave-port modal RHS vector is identically zero. The
/// per-face quadrature `(A / 3) · (t_i · E_t)` vanishes when
/// `E_t = 0` regardless of `β_mode` and face geometry.
#[test]
fn port_modal_rhs_zero_e_t_returns_zero() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let beta = beta_te10_wr90_at_10ghz();
    let rhs = assemble_port_modal_rhs(face_vertices, outward_normal, beta, Vector3::zeros());

    for i in 0..3 {
        let value = rhs[i];
        assert_eq!(
            value.re, 0.0,
            "rhs[{i}].re must be exactly zero with E_t = 0; got {}",
            value.re
        );
        assert_eq!(
            value.im, 0.0,
            "rhs[{i}].im must be exactly zero with E_t = 0; got {}",
            value.im
        );
    }
}

/// E2 DoD criterion 5: with a uniform incident `E_t = x̂` on the unit
/// right triangle (vertices `(0,0,0), (1,0,0), (0,1,0)`), each RHS
/// entry matches the closed-form `2 j β · (A / 3) · (t_i · E_t)`.
///
/// Edge tangents are `t_0 = (1,0,0)`, `t_1 = (−1,1,0)`, `t_2 = (0,−1,0)`,
/// so the per-edge projections are `t_0 · E_t = 1`, `t_1 · E_t = −1`,
/// `t_2 · E_t = 0`. Face area `A = 0.5`. The prefactor is
/// `2 j β · A / 3 = j β / 3`.
#[test]
fn port_modal_rhs_uniform_e_t_dot_normal() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let beta = beta_te10_wr90_at_10ghz();
    let e_t = Vector3::new(1.0, 0.0, 0.0);

    let rhs = assemble_port_modal_rhs(face_vertices, outward_normal, beta, e_t);

    // Expected closed-form values: 2 j β · (A / 3) · (t_i · E_t).
    let face_area = 0.5;
    let prefactor_im = 2.0 * beta * face_area / 3.0; // imaginary part magnitude

    let expected_im = [
        prefactor_im,  // t_0 · E_t = 1
        -prefactor_im, // t_1 · E_t = −1
        0.0,           // t_2 · E_t = 0
    ];

    let tol = 1e-10 * beta;
    for i in 0..3 {
        let value = rhs[i];
        assert!(
            value.re.abs() < tol,
            "rhs[{i}].re must be ~0; got {}",
            value.re
        );
        assert!(
            (value.im - expected_im[i]).abs() < tol,
            "rhs[{i}].im mismatch: expected {}, got {} (tol = {tol:e})",
            expected_im[i],
            value.im
        );
    }
}

/// E2 DoD criterion 6: the wave-port modal RHS lies in the range of
/// the wave-port face block. The 3×3 Gram matrix
/// `G[i][j] = (n̂ × t_i) · (n̂ × t_j)` has a one-dimensional null space
/// spanned by `(1, 1, 1)`, because the three edge tangents of any
/// triangle sum to zero (`Σ t_i = 0`). For `B · e_mode = b_mode / 2`
/// to be solvable, `b_mode` must be orthogonal to this null space —
/// equivalently `Σ b_mode,i = 0`. The RHS formula
/// `b_i = 2 j β · (A / 3) · (t_i · E_t)` satisfies this identity
/// exactly because `Σ_i (t_i · E_t) = (Σ_i t_i) · E_t = 0`.
///
/// This is the face-level necessary self-consistency condition that
/// makes the per-face system `B · e_mode = b_mode / 2` solvable. The
/// assembly layer (Phase 4.fem.eig.2 step E3) inherits this property
/// when scattering per-face contributions into the global driven
/// system.
#[test]
fn port_modal_rhs_complex_symmetric_with_block() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let beta = beta_te10_wr90_at_10ghz();

    // Try several non-trivial incident E_t directions and verify the
    // RHS lies in the range of the wave-port face block (i.e. its
    // entries sum to zero).
    let e_t_choices = [
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.7, -0.3, 0.0),
        // A non-tangential component is dropped by the t_i · E_t dot
        // product; the in-face component (1, 1, 0) still satisfies the
        // sum-to-zero identity.
        Vector3::new(1.0, 1.0, 0.5),
    ];

    let tol = 1e-10 * beta;
    for e_t in &e_t_choices {
        let rhs = assemble_port_modal_rhs(face_vertices, outward_normal, beta, *e_t);
        let sum_im = rhs[0].im + rhs[1].im + rhs[2].im;
        let sum_re = rhs[0].re + rhs[1].re + rhs[2].re;
        assert!(
            sum_re.abs() < tol,
            "Σ rhs[i].re must be 0 (range-of-B condition); got {sum_re:e} for E_t = {e_t:?}"
        );
        assert!(
            sum_im.abs() < tol,
            "Σ rhs[i].im must be 0 (range-of-B condition); got {sum_im:e} for E_t = {e_t:?}"
        );

        // Cross-check: confirm `(1, 1, 1)` is indeed a null vector of
        // the wave-port face block by computing `B · (1, 1, 1)` and
        // asserting every entry is zero.
        let block = assemble_port_face_block(face_vertices, outward_normal, beta, 1.0);
        for i in 0..3 {
            let row_sum_re = block[(i, 0)].re + block[(i, 1)].re + block[(i, 2)].re;
            let row_sum_im = block[(i, 0)].im + block[(i, 1)].im + block[(i, 2)].im;
            assert!(
                row_sum_re.abs() < tol && row_sum_im.abs() < tol,
                "(B · (1,1,1))[{i}] must be zero; got ({row_sum_re:e}, {row_sum_im:e})"
            );
        }
    }
}
