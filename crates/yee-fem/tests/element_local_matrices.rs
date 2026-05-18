//! Integration tests for the 6-edge first-order Nedelec tet element
//! local matrices (Phase 4.fem.eig.0 step T3).
//!
//! Reference: Jin, *FEM in Electromagnetics* 3rd ed. §9.4 (eq. 9.43,
//! quadrature table).

use nalgebra::{SMatrix, Vector3};
use yee_fem::element::{LOCAL_EDGES, assemble_tet_element};

/// The unit reference tetrahedron used in every closed-form check:
/// `v_0 = (0,0,0)`, `v_1 = (1,0,0)`, `v_2 = (0,1,0)`, `v_3 = (0,0,1)`.
///
/// Signed volume `V = 1/6`. The barycentric gradients are
/// `∇λ_0 = (−1,−1,−1)`, `∇λ_1 = (1,0,0)`, `∇λ_2 = (0,1,0)`,
/// `∇λ_3 = (0,0,1)`.
fn unit_reference_tet() -> [Vector3<f64>; 4] {
    [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ]
}

#[test]
fn local_edges_canonical_ordering() {
    // Sanity-check the constant; assembly relies on this exact order.
    assert_eq!(
        LOCAL_EDGES,
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    );
}

#[test]
fn unit_reference_tet_volume_is_one_sixth() {
    // Build the element with ε_r = μ_r = 1 and recover V from K's
    // diagonal: for the reference tet ∇×N_{01} = 2 ∇λ_0 × ∇λ_1 =
    // 2(−1,−1,−1)×(1,0,0) = 2(0, −1, 1) → |curl|² = 8. K[0,0] = V · 8.
    let elem = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let v_from_k = elem.k_local[(0, 0)] / 8.0;
    assert!(
        (v_from_k - 1.0 / 6.0).abs() < 1e-12,
        "expected V = 1/6, got {v_from_k}"
    );
}

#[test]
fn partition_of_unity_gradients_sum_to_zero() {
    // Σ_i ∇λ_i = ∇(Σ_i λ_i) = ∇(1) = 0.
    //
    // The element layer doesn't expose `∇λ_i` directly, but we can
    // recover the gradient identity through the curls:
    //
    //   Σ_α (∇×N_α) · u   for any vector u   = 0
    //
    // is NOT true in general — what IS true is that for any choice of
    // i, Σ_{j ≠ i} ∇λ_i × ∇λ_j = ∇λ_i × (Σ_{j ≠ i} ∇λ_j) = −∇λ_i ×
    // ∇λ_i = 0 (using the partition-of-unity identity).
    //
    // So: sum over the three edges incident to vertex 0 of (∇×N_{0j})
    // — picking the local sign so each runs *out of* vertex 0 — should
    // be zero. In LOCAL_EDGES the three edges (0,1), (0,2), (0,3) all
    // have local-low endpoint 0, so the canonical sum is:
    //
    //   curl_{01} + curl_{02} + curl_{03}
    //     = 2 ∇λ_0 × (∇λ_1 + ∇λ_2 + ∇λ_3)
    //     = 2 ∇λ_0 × (−∇λ_0) = 0.
    //
    // We test this through the curl-curl block: pick e = [1,1,1,0,0,0]
    // (the three vertex-0 edges, canonical signs) and check
    // K_local e = 0 to working precision (because (∇×N_α)·sum = 0 for
    // every α).
    let elem = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let mut e_vec = SMatrix::<f64, 6, 1>::zeros();
    e_vec[0] = 1.0;
    e_vec[1] = 1.0;
    e_vec[2] = 1.0;
    let k_e = elem.k_local * e_vec;
    assert!(
        k_e.norm() < 1e-12,
        "partition-of-unity / gradient-kernel check failed; ‖K · e‖ = {}",
        k_e.norm()
    );
}

#[test]
fn k_local_is_symmetric() {
    let elem = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let asymm = (elem.k_local - elem.k_local.transpose()).norm();
    assert!(asymm < 1e-12, "K_local not symmetric: ‖K − Kᵀ‖ = {asymm}");
}

#[test]
fn m_local_is_symmetric() {
    let elem = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let asymm = (elem.m_local - elem.m_local.transpose()).norm();
    assert!(asymm < 1e-12, "M_local not symmetric: ‖M − Mᵀ‖ = {asymm}");
}

#[test]
fn m_local_is_positive_definite() {
    // M_local is a Gram matrix of linearly-independent basis functions
    // weighted by ε_r > 0, so its smallest eigenvalue must be > 0. We
    // use a small dense symmetric eigen on the 6×6 block.
    let elem = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let dense = nalgebra::DMatrix::from_iterator(6, 6, elem.m_local.iter().copied());
    let eig = nalgebra::SymmetricEigen::new(dense);
    let min_eig = eig
        .eigenvalues
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    assert!(
        min_eig > 1e-10,
        "M_local not positive definite: smallest eigenvalue = {min_eig}"
    );
}

#[test]
fn k_local_kernel_is_dimension_three() {
    // The kernel of the curl-curl bilinear form on the six Nedelec
    // edges of one tet is the discrete gradient subspace, which on a
    // tet has dimension `n_verts − 1 = 3` (one scalar potential
    // gradient per non-anchored vertex).
    let elem = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let dense = nalgebra::DMatrix::from_iterator(6, 6, elem.k_local.iter().copied());
    let eig = nalgebra::SymmetricEigen::new(dense);
    let zero_eig_count = eig.eigenvalues.iter().filter(|&&e| e.abs() < 1e-9).count();
    assert_eq!(
        zero_eig_count,
        3,
        "expected 3 zero eigenvalues (gradient kernel), got {zero_eig_count}; \
         eigenvalues = {:?}",
        eig.eigenvalues.as_slice()
    );
}

#[test]
fn eps_mu_scaling() {
    // Free-space sanity check: ε_r = μ_r = 1 vs ε_r = μ_r = 2.
    //   M_local scales linearly with ε_r → 2× larger at ε_r = 2.
    //   K_local scales with 1/μ_r        → 1/2× at μ_r = 2.
    let base = assemble_tet_element(unit_reference_tet(), 1.0, 1.0);
    let doubled = assemble_tet_element(unit_reference_tet(), 2.0, 2.0);

    let m_ratio_err = (doubled.m_local - 2.0 * base.m_local).norm();
    let k_ratio_err = (doubled.k_local - 0.5 * base.k_local).norm();

    assert!(
        m_ratio_err < 1e-12,
        "M_local ε_r-scaling failed: ‖M(ε=2) − 2 M(ε=1)‖ = {m_ratio_err}"
    );
    assert!(
        k_ratio_err < 1e-12,
        "K_local μ_r-scaling failed: ‖K(μ=2) − 0.5 K(μ=1)‖ = {k_ratio_err}"
    );
}

#[test]
fn uniform_dilation_scales_k_by_inverse_alpha_and_m_by_alpha() {
    // Scale all four vertices by α = 2:
    //   ∇λ_i scales by 1/α        (gradients are inverse-length).
    //   curl ∝ ∇λ × ∇λ            → scales by 1/α².
    //   V    ∝ α³.
    //   K^e ∝ V · |curl|²         → α³ · 1/α⁴ = 1/α       → scales by 1/α.
    //   N_{ij} = λ_i ∇λ_j − …      → scales by 1/α.
    //   ∫ N·N dV ∝ V · |N|²       → α³ · 1/α² = α.        → M ∝ α.
    //
    // (Sanity for the reader: the plan's T3 DoD bullet says "K scales
    // by α, M by α³", which is at odds with the entry-level
    // derivation above. The derivation is correct on inspection: the
    // eigenvalue ratio K/M scales as (1/α) / α = 1/α², so a tet twice
    // as big resonates at half the wavenumber — the physically
    // correct result. The plan wording is surfaced as a finding in
    // the agent report.)
    let base_verts = unit_reference_tet();
    let alpha = 2.0;
    let scaled_verts = [
        base_verts[0] * alpha,
        base_verts[1] * alpha,
        base_verts[2] * alpha,
        base_verts[3] * alpha,
    ];

    let base = assemble_tet_element(base_verts, 1.0, 1.0);
    let scaled = assemble_tet_element(scaled_verts, 1.0, 1.0);

    let k_err = (scaled.k_local - (1.0 / alpha) * base.k_local).norm();
    let m_err = (scaled.m_local - alpha * base.m_local).norm();

    assert!(
        k_err < 1e-12,
        "K dilation scaling: ‖K(2v) − (1/2) K(v)‖ = {k_err}"
    );
    assert!(
        m_err < 1e-12,
        "M dilation scaling: ‖M(2v) − 2 M(v)‖ = {m_err}"
    );

    // Cross-check: the eigenvalue ratio K/M should scale as 1/α². We
    // pick the trace ratio as a cheap scalar proxy (full eig is
    // overkill here; the entry-level scaling above is the rigorous
    // statement).
    let ratio_base = base.k_local.trace() / base.m_local.trace();
    let ratio_scaled = scaled.k_local.trace() / scaled.m_local.trace();
    assert!(
        (ratio_scaled - ratio_base / (alpha * alpha)).abs() < 1e-10,
        "tr(K)/tr(M) should scale by 1/α² = 1/{}: base = {ratio_base}, scaled = {ratio_scaled}",
        alpha * alpha
    );
}

#[test]
fn vertex_permutation_preserves_eigenspectrum() {
    // Permuting (v_1, v_2, v_3) (an even permutation that preserves
    // orientation: (1,2,3) → (2,3,1)) relabels the local edges but
    // leaves the spectrum of K and M invariant.
    let verts = unit_reference_tet();
    let permuted = [verts[0], verts[2], verts[3], verts[1]];

    let base = assemble_tet_element(verts, 1.0, 1.0);
    let perm = assemble_tet_element(permuted, 1.0, 1.0);

    let base_k_dense = nalgebra::DMatrix::from_iterator(6, 6, base.k_local.iter().copied());
    let perm_k_dense = nalgebra::DMatrix::from_iterator(6, 6, perm.k_local.iter().copied());
    let base_eig: Vec<f64> = {
        let mut e: Vec<f64> = nalgebra::SymmetricEigen::new(base_k_dense)
            .eigenvalues
            .iter()
            .copied()
            .collect();
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        e
    };
    let perm_eig: Vec<f64> = {
        let mut e: Vec<f64> = nalgebra::SymmetricEigen::new(perm_k_dense)
            .eigenvalues
            .iter()
            .copied()
            .collect();
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        e
    };

    for (a, b) in base_eig.iter().zip(perm_eig.iter()) {
        assert!(
            (a - b).abs() < 1e-10,
            "K eigenvalues differ under vertex permutation: {a} vs {b}\n\
             base = {base_eig:?}\nperm = {perm_eig:?}"
        );
    }
}
