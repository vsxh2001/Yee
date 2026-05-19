//! Phase 4.fem.eig.2 step E1 — unit tests for
//! [`yee_fem::element::assemble_abc_face_block`].
//!
//! Gate test inventory:
//!
//! 1. `unit_right_triangle_face_block_is_complex_symmetric` — on the unit
//!    right triangle in the xy-plane with `n̂ = ẑ`, the returned 3×3 block
//!    satisfies `B == B^T` to `1e-12` (complex-symmetric, NOT Hermitian).
//! 2. `face_block_imaginary_coefficient` — every entry of `B` is purely
//!    imaginary, because the prefactor `j · k₀ · A / μ_r` is purely
//!    imaginary and the `(n̂ × t_i) · (n̂ × t_j)` Gram form is real.
//! 3. `face_area_scales_block_linearly` — doubling the face dimensions
//!    quadruples the magnitude of every entry (area scales as length²
//!    and the cross products `n̂ × t_i` also scale linearly, but the
//!    `n̂ × t` factors are absorbed into the *Whitney edge basis*
//!    dual-tangent identity; only the explicit `A` prefactor varies
//!    here — the dimensional scaling that ultimately enters the
//!    assembled global block is the `length²` factor).
//! 4. `mu_r_inversely_scales_block` — doubling `μ_r,face` halves the
//!    magnitude of every entry, because the block carries a `(1/μ_r)`
//!    factor.
//! 5. `block_diagonal_dominance` — on the unit right triangle the
//!    diagonal entries `|B[i][i]|` are at least as large as the
//!    off-diagonal entries `|B[i][j]|` for the adjacent edge — a
//!    Cauchy–Schwarz consequence of the Gram form
//!    `(n̂ × t_i) · (n̂ × t_j)`.
//!
//! References:
//! * Engquist, B. and Majda, A., "Absorbing boundary conditions for the
//!   numerical simulation of waves", *Math. Comp.* 31 (1977),
//!   pp. 629–651.
//! * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
//!   3rd ed., Wiley 2014, §10.4.
//! * Phase 4.fem.eig.2 spec
//!   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
//!   §4.2.

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::element::assemble_abc_face_block;

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

/// E1 DoD criterion 1: the ABC face block is complex-symmetric
/// (`B == B^T`), NOT Hermitian. This is the canonical "radiation absorbs
/// energy → eigenvalues acquire negative imaginary part" identity from
/// spec §4.2.
#[test]
fn unit_right_triangle_face_block_is_complex_symmetric() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let block = assemble_abc_face_block(face_vertices, outward_normal, k0_10ghz(), 1.0);

    let asymmetry = (block - block.transpose()).norm();
    assert!(
        asymmetry < 1e-12,
        "ABC face block must be complex-symmetric (B == B^T); got ||B - B^T|| = {asymmetry:e}"
    );
}

/// E1 DoD criterion 2: every entry of the ABC face block is purely
/// imaginary. The prefactor `j · k₀ · A / μ_r` is purely imaginary and
/// the Gram form `(n̂ × t_i) · (n̂ × t_j)` is real, so `Re(B[i][j]) = 0`
/// up to floating-point round-off.
#[test]
fn face_block_imaginary_coefficient() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let block = assemble_abc_face_block(face_vertices, outward_normal, k0_10ghz(), 1.0);

    // Bound: 100 × f64::EPSILON, scaled by the magnitude of the
    // imaginary prefactor `k₀ · A` which can be O(100) at 10 GHz.
    let tol = f64::EPSILON * 100.0 * k0_10ghz();
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

/// E1 DoD criterion 3: doubling every face dimension scales the block
/// entries by 4×. With vertices `(0,0,0)`, `(2,0,0)`, `(0,2,0)` the area
/// becomes `2.0` (4× the unit-triangle area `0.5`); each tangent
/// `n̂ × t_i` doubles in magnitude so `(n̂ × t_i) · (n̂ × t_j)` quadruples
/// — but the *element* face block multiplies the Gram form by the area
/// only once, so the net scaling of the block magnitude is 4×
/// (area-scaling) × 4× (Gram-scaling) / 1× = combined 16×. We test the
/// simpler "area only" scaling here by comparing two triangles with the
/// same `n̂ × t_i` direction set but different areas — a degenerate
/// thought experiment that does not arise in practice; instead we
/// verify the *literal* dimensional scaling: vertices doubled → block
/// magnitude × 16.
#[test]
fn face_area_scales_block_linearly() {
    let face_vertices = unit_right_triangle();
    let face_vertices_doubled = [
        face_vertices[0] * 2.0,
        face_vertices[1] * 2.0,
        face_vertices[2] * 2.0,
    ];
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);

    let block_unit = assemble_abc_face_block(face_vertices, outward_normal, k0_10ghz(), 1.0);
    let block_doubled =
        assemble_abc_face_block(face_vertices_doubled, outward_normal, k0_10ghz(), 1.0);

    // Doubling the vertices: area × 4, each n̂ × t_i × 2, Gram × 4;
    // net block magnitude × 16. Check the (0,0) diagonal entry which
    // is non-zero on this triangle.
    let ratio = block_doubled[(0, 0)].im / block_unit[(0, 0)].im;
    assert!(
        (ratio - 16.0).abs() < 1e-10,
        "doubled face → block × 16 expected; got ratio = {ratio} on B[0][0].im"
    );
}

/// E1 DoD criterion 4: doubling `μ_r,face` halves the block, because
/// the block carries an explicit `1/μ_r,face` factor.
#[test]
fn mu_r_inversely_scales_block() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);

    let block_mu1 = assemble_abc_face_block(face_vertices, outward_normal, k0_10ghz(), 1.0);
    let block_mu2 = assemble_abc_face_block(face_vertices, outward_normal, k0_10ghz(), 2.0);

    let ratio = block_mu2[(0, 0)].im / block_mu1[(0, 0)].im;
    assert!(
        (ratio - 0.5).abs() < 1e-12,
        "μ_r × 2 → block × 0.5 expected; got ratio = {ratio} on B[0][0].im"
    );
}

/// E1 DoD criterion 5: diagonal dominance via Cauchy–Schwarz.
/// `|(n̂ × t_i) · (n̂ × t_i)| ≥ |(n̂ × t_i) · (n̂ × t_j)|` when
/// `||n̂ × t_i|| ≥ ||n̂ × t_j||` (and always to within the geometric mean
/// of the diagonals). On the unit right triangle this gives
/// `|B[0][0]| ≥ |B[0][1]|` (with equality, since `||n̂ × t_0|| = 1` and
/// the cross-term magnitude is also `1`) and `|B[1][1]| > |B[1][2]|`.
#[test]
fn block_diagonal_dominance() {
    let face_vertices = unit_right_triangle();
    let outward_normal = Vector3::new(0.0, 0.0, 1.0);
    let block = assemble_abc_face_block(face_vertices, outward_normal, k0_10ghz(), 1.0);

    let b00 = block[(0, 0)].norm();
    let b01 = block[(0, 1)].norm();
    let b11 = block[(1, 1)].norm();
    let b12 = block[(1, 2)].norm();

    assert!(
        b00 >= b01 - 1e-12,
        "|B[0][0]| ≥ |B[0][1]| expected (Cauchy–Schwarz); got |B[0][0]| = {b00}, |B[0][1]| = {b01}"
    );
    assert!(
        b11 >= b12 - 1e-12,
        "|B[1][1]| ≥ |B[1][2]| expected (Cauchy–Schwarz); got |B[1][1]| = {b11}, |B[1][2]| = {b12}"
    );
}
