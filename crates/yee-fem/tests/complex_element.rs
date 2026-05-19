//! Phase 4.fem.eig.1 step D1 — complex-coefficient lift of
//! [`yee_fem::element::assemble_tet_element_complex`].
//!
//! Gate test inventory:
//!
//! 1. `complex_matches_real_for_pure_real_eps_mu` — pure-real `ε`, `μ`
//!    inputs produce `Im(K) = Im(M) = 0` and `Re(K)`, `Re(M)` bit-for-bit
//!    identical to the Phase 4.fem.eig.0 real
//!    [`yee_fem::element::assemble_tet_element`] output.
//! 2. `imaginary_eps_produces_imaginary_M_diagonal` — pure-imaginary
//!    `ε = j` keeps `K_local` real (curl-curl is independent of ε) and
//!    rotates the entire `M_local` into the imaginary axis.
//! 3. `lossy_drude_like_eps_keeps_M_complex_K_real` — complex `ε` with
//!    a Drude-like negative-imaginary part populates both `Re(M)` and
//!    `Im(M)` with the expected signs; `K_local` stays real-only.
//! 4. `lossy_mu_propagates_to_K_only` — complex `μ` populates both
//!    `Re(K)` and `Im(K)` (through `1/μ`); `M_local` stays real-only.
//! 5. `complex_local_blocks_are_symmetric` — `K_local` and `M_local`
//!    are complex *symmetric* (transposed, not Hermitian) per ADR-0039
//!    / spec §11 for any (eps, mu).
//!
//! References:
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
//!   §6 (element-layer changes).
//! * `docs/src/decisions/0039-phase-4-fem-eig-1-dispersive-scope.md`.

#![allow(non_snake_case)]
//
// Test names below intentionally carry capital `K` / `M` to match the
// `K_local` / `M_local` field names in the source. The
// `#![allow(non_snake_case)]` keeps clippy's snake_case lint happy on
// the few lines where the capitals appear inside the test names.

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_fem::element::{LOCAL_EDGES, assemble_tet_element, assemble_tet_element_complex};

/// The unit reference tetrahedron used in every closed-form check:
/// `v_0 = (0,0,0)`, `v_1 = (1,0,0)`, `v_2 = (0,1,0)`, `v_3 = (0,0,1)`.
fn unit_reference_tet() -> [Vector3<f64>; 4] {
    [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ]
}

/// A non-axis-aligned tet with positive signed volume, to exercise a
/// non-trivial barycentric-gradient code path.
fn tilted_tet() -> [Vector3<f64>; 4] {
    [
        Vector3::new(0.1, 0.2, 0.3),
        Vector3::new(1.2, 0.3, 0.1),
        Vector3::new(0.4, 1.5, 0.2),
        Vector3::new(0.3, 0.4, 1.7),
    ]
}

/// D1 gate criterion: pure-real `ε`, `μ` inputs to the complex path
/// reproduce the v0 real path bit-for-bit.
#[test]
fn complex_matches_real_for_pure_real_eps_mu() {
    for tet in [unit_reference_tet(), tilted_tet()] {
        for eps_r in [1.0_f64, 2.5, 9.8] {
            for mu_r in [1.0_f64, 1.7] {
                let real_elem = assemble_tet_element(tet, eps_r, mu_r);
                let complex_elem = assemble_tet_element_complex(
                    tet,
                    Complex64::new(eps_r, 0.0),
                    Complex64::new(mu_r, 0.0),
                );

                for alpha in 0..6 {
                    for beta in 0..6 {
                        // Im parts identically zero for real (eps, mu).
                        assert_eq!(
                            complex_elem.k_local[(alpha, beta)].im,
                            0.0,
                            "K_local[{alpha},{beta}].im != 0 for real input \
                             (eps_r = {eps_r}, mu_r = {mu_r})",
                        );
                        assert_eq!(
                            complex_elem.m_local[(alpha, beta)].im,
                            0.0,
                            "M_local[{alpha},{beta}].im != 0 for real input \
                             (eps_r = {eps_r}, mu_r = {mu_r})",
                        );
                        // Re parts bit-for-bit identical to v0.
                        let re_k = complex_elem.k_local[(alpha, beta)].re;
                        let re_m = complex_elem.m_local[(alpha, beta)].re;
                        let v0_k = real_elem.k_local[(alpha, beta)];
                        let v0_m = real_elem.m_local[(alpha, beta)];
                        assert_eq!(
                            re_k, v0_k,
                            "K_local[{alpha},{beta}] mismatch: complex.re = {re_k}, \
                             v0 = {v0_k} (eps_r = {eps_r}, mu_r = {mu_r})",
                        );
                        assert_eq!(
                            re_m, v0_m,
                            "M_local[{alpha},{beta}] mismatch: complex.re = {re_m}, \
                             v0 = {v0_m} (eps_r = {eps_r}, mu_r = {mu_r})",
                        );
                    }
                }
            }
        }
    }
}

/// Pure-imaginary `ε = j` keeps `K_local` real (curl-curl is independent
/// of ε) and rotates `M_local` from the real to the imaginary axis.
#[test]
fn imaginary_eps_produces_imaginary_M_diagonal() {
    let tet = unit_reference_tet();
    let eps = Complex64::new(0.0, 1.0);
    let mu = Complex64::new(1.0, 0.0);
    let real_elem = assemble_tet_element(tet, 1.0, 1.0);
    let complex_elem = assemble_tet_element_complex(tet, eps, mu);

    for alpha in 0..6 {
        for beta in 0..6 {
            // K_local is curl-curl and is independent of ε. With μ=1 it
            // equals the v0 real K_local bit-for-bit.
            assert_eq!(
                complex_elem.k_local[(alpha, beta)],
                Complex64::new(real_elem.k_local[(alpha, beta)], 0.0),
                "K_local[{alpha},{beta}] must equal v0 K_local when μ=1 \
                 and ε is purely imaginary",
            );
            // M_local = ε · M_real. For ε = j the result is purely
            // imaginary with Im = v0 M_local entry.
            assert_eq!(
                complex_elem.m_local[(alpha, beta)].re,
                0.0,
                "M_local[{alpha},{beta}].re must be zero for ε = j",
            );
            assert_eq!(
                complex_elem.m_local[(alpha, beta)].im,
                real_elem.m_local[(alpha, beta)],
                "M_local[{alpha},{beta}].im must equal v0 M_local for ε = j",
            );
        }
    }
}

/// Drude-like complex `ε = 3.78 − 0.5j` (lossy fused-silica
/// approximation from fem-eig-002): both `Re(M)` and `Im(M)` are
/// populated with the expected `ε · M_real` decomposition; `K_local`
/// stays real (`μ` is real).
#[test]
fn lossy_drude_like_eps_keeps_M_complex_K_real() {
    let tet = tilted_tet();
    let eps = Complex64::new(3.78, -0.5);
    let mu = Complex64::new(1.0, 0.0);

    // Reference: the v0 real `M_local` at ε_r = 1, μ_r = 1.
    let unit_elem = assemble_tet_element(tet, 1.0, 1.0);
    let lossy_elem = assemble_tet_element_complex(tet, eps, mu);

    for alpha in 0..6 {
        for beta in 0..6 {
            // K: independent of ε, real because μ is real.
            assert_eq!(
                lossy_elem.k_local[(alpha, beta)].im,
                0.0,
                "K_local[{alpha},{beta}].im must be zero for real μ",
            );
            let expected_k = unit_elem.k_local[(alpha, beta)];
            assert!(
                (lossy_elem.k_local[(alpha, beta)].re - expected_k).abs() < 1e-14,
                "K_local[{alpha},{beta}].re = {} expected {expected_k}",
                lossy_elem.k_local[(alpha, beta)].re,
            );

            // M: ε · M_real, decomposed component-wise.
            let m_ref = unit_elem.m_local[(alpha, beta)];
            let expected = eps * Complex64::new(m_ref, 0.0);
            let got = lossy_elem.m_local[(alpha, beta)];
            assert!(
                (got - expected).norm() < 1e-13,
                "M_local[{alpha},{beta}] = {got} expected {expected}",
            );
        }
    }
}

/// Complex `μ` populates both `Re(K)` and `Im(K)` (through `1/μ`);
/// `M_local` stays real (`ε` is real).
#[test]
fn lossy_mu_propagates_to_K_only() {
    let tet = unit_reference_tet();
    let eps = Complex64::new(1.0, 0.0);
    let mu = Complex64::new(2.0, 0.3);

    let unit_elem = assemble_tet_element(tet, 1.0, 1.0);
    let lossy_elem = assemble_tet_element_complex(tet, eps, mu);

    let inv_mu = Complex64::new(1.0, 0.0) / mu;

    for alpha in 0..6 {
        for beta in 0..6 {
            // K: (1/μ) · K_real(μ_r=1). The complex inv_mu factor must
            // appear cleanly.
            let k_real_unit = unit_elem.k_local[(alpha, beta)];
            let expected = inv_mu * Complex64::new(k_real_unit, 0.0);
            let got = lossy_elem.k_local[(alpha, beta)];
            assert!(
                (got - expected).norm() < 1e-13,
                "K_local[{alpha},{beta}] = {got} expected {expected}",
            );

            // M: real (ε real).
            assert_eq!(
                lossy_elem.m_local[(alpha, beta)].im,
                0.0,
                "M_local[{alpha},{beta}].im must be zero for real ε",
            );
            assert_eq!(
                lossy_elem.m_local[(alpha, beta)].re,
                unit_elem.m_local[(alpha, beta)],
                "M_local[{alpha},{beta}].re must equal v0 M_local",
            );
        }
    }
}

/// Both `K_local` and `M_local` are complex *symmetric* (transposed,
/// not Hermitian) per ADR-0039 / spec §11. This is the load-bearing
/// invariant for the D5 Hellmann–Feynman derivative.
#[test]
fn complex_local_blocks_are_symmetric() {
    // Use a fully-lossy material (both ε and μ complex) to exercise
    // every entry of both blocks.
    let tet = tilted_tet();
    let eps = Complex64::new(3.78, -0.5);
    let mu = Complex64::new(2.0, 0.3);
    let elem = assemble_tet_element_complex(tet, eps, mu);

    for alpha in 0..6 {
        for beta in 0..alpha {
            let k_ab = elem.k_local[(alpha, beta)];
            let k_ba = elem.k_local[(beta, alpha)];
            assert!(
                (k_ab - k_ba).norm() < 1e-13,
                "K_local complex-symmetry violated at ({alpha},{beta}): \
                 K[{alpha},{beta}] = {k_ab}, K[{beta},{alpha}] = {k_ba}",
            );

            let m_ab = elem.m_local[(alpha, beta)];
            let m_ba = elem.m_local[(beta, alpha)];
            assert!(
                (m_ab - m_ba).norm() < 1e-13,
                "M_local complex-symmetry violated at ({alpha},{beta}): \
                 M[{alpha},{beta}] = {m_ab}, M[{beta},{alpha}] = {m_ba}",
            );
        }
    }
    // LOCAL_EDGES sanity (referenced in module docs).
    assert_eq!(LOCAL_EDGES.len(), 6);
}
