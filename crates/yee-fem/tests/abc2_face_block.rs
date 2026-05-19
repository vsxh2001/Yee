//! Phase 4.fem.eig.3 step F3 — unit tests for
//! [`yee_fem::element::assemble_abc2_face_block`].
//!
//! Gate test inventory (per the F3 brief):
//!
//! 1. [`abc2_complex_symmetric`] — the 2nd-order ABC face block stays
//!    complex-symmetric (`B == B^T`, NOT Hermitian). Both `R_1` and
//!    `R_2` are real-symmetric Gram matrices with scalar (imaginary
//!    `R_1`, real `R_2`) prefactors, so the composite preserves
//!    `B == B^T` to machine precision.
//! 2. [`abc2_first_order_term_matches_1st_order_helper`] — at the
//!    high-`k₀` limit, the `R_2 / (2 k₀)` correction is suppressed by
//!    `1/k₀` relative to `k₀ · R_1`, so the 2nd-order block tends to
//!    the 1st-order block as `k₀ → ∞`. This pins the sign and
//!    normalisation of the `R_1` contribution against the existing
//!    [`yee_fem::element::assemble_abc_face_block`] helper.
//! 3. [`abc2_curl_term_imaginary_zero`] — the `R_2` curl correction
//!    contributes only to the **real** part of the block, because its
//!    scalar prefactor `−1/(2 k₀)` is real and the `R_2` Gram form is
//!    real. Combined with the imaginary `+ j k₀` prefactor on `R_1`,
//!    the 2nd-order block has Re ≠ 0 *and* Im ≠ 0 — this is the load-
//!    bearing distinguishing feature from the 1st-order block, which is
//!    purely imaginary.
//! 4. [`abc2_curl_correction_scales_inversely_with_k0`] — the
//!    Engquist–Majda 1979 frequency-scaling identity: at high frequency
//!    the curl correction is smaller relative to `R_1` (1/k₀ scaling).
//!    Concretely, `‖Re(B(2 k₀))‖ / ‖Im(B(2 k₀))‖` is approximately
//!    `1/4` × `‖Re(B(k₀))‖ / ‖Im(B(k₀))‖` (the real part scales as
//!    `1/k₀` and the imaginary part scales as `k₀`, so the ratio
//!    scales as `1/k₀²`).
//!
//! References:
//! * Engquist, B. and Majda, A., "Radiation boundary conditions for the
//!   numerical simulation of waves", *Math. Comp.* 31 (1977),
//!   pp. 629–651, and *IEEE Trans. Antennas Propag.* 27(5) (1979)
//!   p. 661, eq. 9 — the 2nd-order ABC derivation.
//! * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
//!   3rd ed., Wiley 2014, §10.4 (reflection-floor tables).
//! * Phase 4.fem.eig.3 spec
//!   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
//!   §4.2.

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::element::{assemble_abc_face_block, assemble_abc2_face_block};

/// Canonical unit right-triangle face in the xy-plane with outward
/// normal `+ẑ`. Vertices in CCW order seen from `+ẑ`.
fn unit_right_triangle() -> [Vector3<f64>; 3] {
    [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ]
}

/// `k₀` at 10 GHz in vacuum (`λ ≈ 0.03 m`).
fn k0_10ghz() -> f64 {
    2.0 * PI / 0.03
}

/// F3 DoD criterion 1: the 2nd-order ABC face block is complex-symmetric
/// (`B == B^T`), NOT Hermitian. Both `R_1` and `R_2` are real-symmetric
/// Gram matrices and the prefactors are scalars.
#[test]
fn abc2_complex_symmetric() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let block = assemble_abc2_face_block(face_vertices, outward_normal, k0_10ghz(), 1.0);

    let asymmetry = (block - block.transpose()).norm();
    assert!(
        asymmetry < 1e-10,
        "2nd-order ABC face block must be complex-symmetric (B == B^T); \
         got ||B - B^T|| = {asymmetry:e}"
    );
}

/// F3 DoD criterion 2: at the high-`k₀` asymptote, the imaginary part
/// of `assemble_abc2_face_block` matches the imaginary part of
/// `assemble_abc_face_block` (the `R_1` contribution is shared by both
/// helpers; the additional `R_2` term contributes only real entries —
/// see [`abc2_curl_term_imaginary_zero`]). We test this directly by
/// asserting `Im(abc2) == Im(abc1)` exactly.
#[test]
fn abc2_first_order_term_matches_1st_order_helper() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let k0 = k0_10ghz();

    let abc1 = assemble_abc_face_block(face_vertices, outward_normal, k0, 1.0);
    let abc2 = assemble_abc2_face_block(face_vertices, outward_normal, k0, 1.0);

    // Imaginary parts must agree bit-for-bit: the F3 implementation
    // re-uses assemble_abc_face_block for the R_1 contribution and the
    // R_2 correction has a real (not imaginary) prefactor.
    for i in 0..3 {
        for j in 0..3 {
            let d_im = abc2[(i, j)].im - abc1[(i, j)].im;
            assert!(
                d_im.abs() < 1e-12,
                "Im(abc2[{i}][{j}]) must match Im(abc1[{i}][{j}]); \
                 got Δ = {d_im:e}"
            );
        }
    }
}

/// F3 DoD criterion 3: the `R_2` curl correction has prefactor
/// `−1/(2 k₀)` (real, no `j`), so the curl-correction part contributes
/// **real** block entries. Combined with the purely-imaginary 1st-order
/// part, the total 2nd-order block has Re ≠ 0 AND Im ≠ 0.
///
/// This is the load-bearing distinguishing feature from the 1st-order
/// helper, which is purely imaginary.
#[test]
fn abc2_curl_term_imaginary_zero() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let k0 = k0_10ghz();

    let abc1 = assemble_abc_face_block(face_vertices, outward_normal, k0, 1.0);
    let abc2 = assemble_abc2_face_block(face_vertices, outward_normal, k0, 1.0);

    // The difference abc2 - abc1 should be PURELY REAL (only the
    // R_2 curl correction differs).
    let mut max_im_diff = 0.0_f64;
    let mut max_re_diff = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let d = abc2[(i, j)] - abc1[(i, j)];
            max_im_diff = max_im_diff.max(d.im.abs());
            max_re_diff = max_re_diff.max(d.re.abs());
        }
    }
    assert!(
        max_im_diff < 1e-12,
        "abc2 - abc1 must be purely real (R_2 prefactor is real); \
         got max |Im(Δ)| = {max_im_diff:e}"
    );
    assert!(
        max_re_diff > 1e-6,
        "abc2 - abc1 must be non-trivially real (R_2 contributes); \
         got max |Re(Δ)| = {max_re_diff:e}"
    );

    // And the full block has non-trivial Re *and* Im.
    let mut max_re = 0.0_f64;
    let mut max_im = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            max_re = max_re.max(abc2[(i, j)].re.abs());
            max_im = max_im.max(abc2[(i, j)].im.abs());
        }
    }
    assert!(
        max_re > 1e-6,
        "2nd-order block must have Re ≠ 0 (R_2 contribution); \
         got max |Re| = {max_re:e}"
    );
    assert!(
        max_im > 1e-6,
        "2nd-order block must have Im ≠ 0 (R_1 contribution); \
         got max |Im| = {max_im:e}"
    );
}

/// F3 DoD criterion 4: at high frequency the curl correction is smaller
/// relative to `R_1` (the Engquist–Majda 1979 1/k₀ scaling identity).
///
/// `Im(B)` scales as `k₀` (1st-order prefactor `+ j k₀ · A / μ_r`);
/// `Re(B)` scales as `1/k₀` (2nd-order prefactor `−A / (2 k₀ μ_r)`).
/// So `‖Re‖ / ‖Im‖` should scale as `1/k₀²`.
///
/// Concretely: doubling `k₀` should quarter the `‖Re‖ / ‖Im‖` ratio to
/// within Gauss-quadrature accuracy.
#[test]
fn abc2_curl_correction_scales_inversely_with_k0() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);

    let k0_a = k0_10ghz();
    let k0_b = 2.0 * k0_a;

    let block_a = assemble_abc2_face_block(face_vertices, outward_normal, k0_a, 1.0);
    let block_b = assemble_abc2_face_block(face_vertices, outward_normal, k0_b, 1.0);

    // Pick the (0, 0) diagonal entry which is non-zero on the unit
    // right triangle. Re scales as 1/k₀, Im scales as k₀.
    let re_a = block_a[(0, 0)].re;
    let re_b = block_b[(0, 0)].re;
    let im_a = block_a[(0, 0)].im;
    let im_b = block_b[(0, 0)].im;

    // Re ratio: should be (k₀_a / k₀_b) = 0.5.
    let re_ratio = re_b / re_a;
    assert!(
        (re_ratio - 0.5).abs() < 1e-10,
        "Re part should scale as 1/k₀: ratio at 2·k₀ vs k₀ expected 0.5, \
         got {re_ratio}"
    );

    // Im ratio: should be (k₀_b / k₀_a) = 2.
    let im_ratio = im_b / im_a;
    assert!(
        (im_ratio - 2.0).abs() < 1e-10,
        "Im part should scale as k₀: ratio at 2·k₀ vs k₀ expected 2.0, \
         got {im_ratio}"
    );

    // Combined: ‖Re‖/‖Im‖ at 2·k₀ should be (1/4) × that at k₀.
    let ratio_a = re_a.abs() / im_a.abs();
    let ratio_b = re_b.abs() / im_b.abs();
    let drop = ratio_b / ratio_a;
    assert!(
        (drop - 0.25).abs() < 1e-10,
        "‖Re‖/‖Im‖ should scale as 1/k₀² (Engquist-Majda 1979); \
         at 2·k₀ expected 0.25 × that at k₀, got drop = {drop}"
    );
}
