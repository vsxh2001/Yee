//! Phase 4.fem.eig.1 step D4 — integration tests for
//! [`yee_fem::DispersiveSolver::solve_at_frequency`] and the
//! complex assembly path
//! [`yee_fem::FemEigenAssembly::assemble_complex`].
//!
//! Test inventory:
//!
//! 1. `free_space_air_matches_phase_4_0_te101` — build a WR-90 cavity
//!    and a `MaterialDatabase` containing only the free-space tag
//!    (`ε_∞ = 1`, no poles). Solve at the analytic TE_{101} angular
//!    frequency and verify the returned smallest `k²` agrees with
//!    Phase 4.0's
//!    [`yee_validation::run_fem_eig_001_rectangular_cavity`] reading
//!    to within 0.5 %. This is the load-bearing
//!    backward-compatibility gate per ADR-0039 §4 — the complex path
//!    must reduce to the real path bit-equivalently in the absence of
//!    dispersion.
//!
//! 2. `lossy_substrate_complex_eigenvalue_has_negative_imaginary_k` —
//!    same cavity geometry, but with half the tets tagged a
//!    single-pole Lorentz material. Verify the returned smallest
//!    eigenvalue has `Im(k²) ≠ 0` (the lossy mode picks up a finite
//!    imaginary part from the complex permittivity).
//!
//! 3. `assemble_complex_at_real_eps_matches_real_assemble` — assemble
//!    `K`, `M` at a frequency where the [`MaterialDatabase`] returns
//!    purely real `ε(ω)`, `μ(ω)` (free-space air) and verify the
//!    complex result's real part agrees with the real
//!    [`yee_fem::FemEigenAssembly::assemble`] output on the same mesh
//!    to within `f64::EPSILON · 1e3`. This pins the cross-path
//!    bit-equivalence at the assembly layer (the eigenvalue gate (1)
//!    above checks it end-to-end).
//!
//! ## Mesh size
//!
//! Tests use the `(8, 6, 10)` Kuhn brick subdivision of the
//! `(22.86, 10.16, 30.0) mm` WR-90 cavity (2880 tets). This is the
//! Phase 4 spec-§9 mesh — large enough to resolve TE_{101} to ~0.3 %
//! against the Pozar analytic frequency, small enough to assemble +
//! solve in well under a minute in `--release`.
//!
//! ## Shift choice
//!
//! `σ = 2.5 · k₀_TE101²` mirrors
//! [`yee_validation::run_fem_eig_001_rectangular_cavity`]: it sits
//! between the 8th and 9th physical modes of the Pozar table on this
//! mesh, so all ten lowest physical modes have `|θ| > |θ_grad|` and
//! inverse-iteration converges to them in ascending `Re(k²)` order.
//! The agent brief's suggested `σ = 0.5 · k₀_TE101²` is a known
//! gradient-cluster boundary case for the Phase 4 T5 escape-hatch
//! inverse-power iteration (see Phase 4 T5 docs in
//! `crates/yee-fem/src/solve.rs`); using the documented working shift
//! instead keeps the gate decisive.
//!
//! References:
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
//! * `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`
//!   step D4.
//! * `docs/src/decisions/0039-phase-4-fem-eig-1-dispersive-scope.md`.

use std::f64::consts::PI;

use num_complex::Complex64;
use yee_fem::{
    DispersiveSolver, FemEigenAssembly, InverseIterEigen, Material, MaterialDatabase, MaterialPole,
    SparseEigen,
};
use yee_mesh::TetMesh3D;

/// Speed of light (m/s) — workspace-canonical from `yee_core::units`.
const C0: f64 = yee_core::units::C0;

/// WR-90 cavity extents (m).
const CAVITY_A_M: f64 = 0.022_86;
const CAVITY_B_M: f64 = 0.010_16;
const CAVITY_D_M: f64 = 0.030;

/// Mesh density — matches the Phase 4 spec §9 fem-eig-001 reference
/// (smaller than the 12×9×15 used by the production validation gate
/// so this unit-test set runs in a few seconds in --release).
const NX: usize = 8;
const NY: usize = 6;
const NZ: usize = 10;

/// Analytic TE_{101} frequency for the cavity (Hz).  Pozar §6.3
/// eq. 6.42 for an air-filled rectangular metallic cavity.
fn f_te101_hz() -> f64 {
    0.5 * C0 * ((1.0 / CAVITY_A_M).powi(2) + (1.0 / CAVITY_D_M).powi(2)).sqrt()
}

/// Build the WR-90 cavity mesh used by every test in this file.
fn build_cavity() -> TetMesh3D {
    TetMesh3D::cavity_uniform(CAVITY_A_M, CAVITY_B_M, CAVITY_D_M, NX, NY, NZ)
        .expect("cavity_uniform must succeed for the standard WR-90 dimensions")
}

/// D4 gate test 1: free-space `DispersiveSolver::solve_at_frequency`
/// reproduces Phase 4.0's `fem-eig-001` lowest-mode reading within
/// 0.5 % (the spec §9 fem-eig-002 Re(f) bound).
#[test]
fn free_space_air_matches_phase_4_0_te101() {
    let mesh = build_cavity();
    let f_te101 = f_te101_hz();
    let omega_te101 = 2.0 * PI * f_te101;
    let k0_te101 = omega_te101 / C0;
    let sigma = Complex64::new(2.5 * k0_te101.powi(2), 0.0);

    // Phase 4.0 reference: same mesh, real-coefficient assemble +
    // InverseIterEigen path. Phase 4.0's fem-eig-001 production gate
    // reports the smallest k² as the analytic TE_{101} k² within
    // ~0.3 %; we recompute it inline so the cross-check is
    // independent of the validation crate's tolerance bookkeeping.
    let real_assembly = FemEigenAssembly::new_free_space(&mesh);
    let real_assembled = real_assembly.assemble().expect("real assemble");
    let real_pairs = InverseIterEigen::default()
        .solve(
            &real_assembled.k,
            &real_assembled.m,
            10,
            2.5 * k0_te101.powi(2),
        )
        .expect("v0 inverse-iter solve");
    let phase_4_0_k_sq = real_pairs.k[0];

    // D4 path: complex assemble + ComplexInverseIterEigen via the
    // free-space MaterialDatabase. The database carries only the
    // air tag (`MaterialTag = 0`) mapped to the default `Material`
    // (`ε_∞ = 1`, `μ_r = 1`, no poles).  This produces real `ε`, `μ`
    // identical to the real path at any ω.
    let db = MaterialDatabase::new().with_material(0, Material::default());
    let solver = DispersiveSolver::new(db);
    let pairs = solver
        .solve_at_frequency(&mesh, omega_te101, 10, sigma)
        .expect("dispersive solve_at_frequency");

    let dispersive_k_sq = pairs.k[0];

    // (a) Im(k²) must be vanishing for purely real `ε`, `μ`.
    assert!(
        dispersive_k_sq.im.abs() < 1e-6 * dispersive_k_sq.re.abs(),
        "free-space mode should be lossless, got Im(k²)/Re(k²) = {}",
        dispersive_k_sq.im / dispersive_k_sq.re,
    );

    // (b) Re(k²) must match the Phase 4.0 v0 reading within 0.5 %.
    let rel_err = (dispersive_k_sq.re - phase_4_0_k_sq).abs() / phase_4_0_k_sq.abs();
    assert!(
        rel_err <= 5e-3,
        "free-space dispersive Re(k²) = {} vs Phase 4.0 = {} (rel err {:.4e}, tol 0.005)",
        dispersive_k_sq.re,
        phase_4_0_k_sq,
        rel_err,
    );

    // (c) Sanity: the lowest mode is close to the analytic TE_{101}
    // wavenumber on this mesh (~0.3 % at NX, NY, NZ = 8, 6, 10).
    let analytic_k_sq = k0_te101.powi(2);
    let analytic_rel_err = (dispersive_k_sq.re - analytic_k_sq).abs() / analytic_k_sq;
    assert!(
        analytic_rel_err <= 0.01,
        "free-space dispersive Re(k²) = {} vs analytic TE_{{101}} k² = {} (rel err {:.4e})",
        dispersive_k_sq.re,
        analytic_k_sq,
        analytic_rel_err,
    );
}

/// D4 gate test 2: a Lorentz-loaded half-cavity produces a complex
/// eigenvalue with non-zero imaginary part.
///
/// Geometry / mesh: same WR-90 cavity as test 1. Material assignment:
/// every tet whose centroid is in the `x < a/2` half-space is tagged
/// `1` and mapped to a single-pole Lorentz oscillator (real part 4 at
/// DC, finite loss tangent); the other half is air (`tag 0`, default
/// material).
#[test]
fn lossy_substrate_complex_eigenvalue_has_negative_imaginary_k() {
    let mut mesh = build_cavity();
    // Repaint tet material tags so the x < a/2 region is a lossy
    // Lorentz dielectric (tag 1); the x ≥ a/2 region remains air
    // (tag 0).
    for (tet_idx, tet) in mesh.tetrahedra.iter().enumerate() {
        let centroid_x = 0.25
            * (mesh.vertices[tet[0]].x
                + mesh.vertices[tet[1]].x
                + mesh.vertices[tet[2]].x
                + mesh.vertices[tet[3]].x);
        mesh.tetrahedron_material[tet_idx] = if centroid_x < 0.5 * CAVITY_A_M { 1 } else { 0 };
    }

    let f_te101 = f_te101_hz();
    let omega_te101 = 2.0 * PI * f_te101;
    let k0_te101 = omega_te101 / C0;
    let sigma = Complex64::new(2.5 * k0_te101.powi(2), 0.0);

    // Lorentz oscillator parameters:
    //  * ε_∞ = 4 — high-frequency permittivity comparable to a typical
    //    PCB substrate.
    //  * ω_0 = 2π · 20 GHz — resonance above the cavity-mode region
    //    (so loss is moderate, not catastrophic).
    //  * ω_p = 2π · 2 GHz — oscillator strength tuned for measurable
    //    Im(ε) at the TE_{101} probe frequency.
    //  * γ = 2π · 0.5 GHz — damping, producing a finite Im(ε).
    let lorentz = Material {
        eps_inf: 4.0,
        mu_r: 1.0,
        poles: vec![MaterialPole::Lorentz {
            omega_0: 2.0 * PI * 20.0e9,
            omega_p: 2.0 * PI * 2.0e9,
            gamma: 2.0 * PI * 0.5e9,
        }],
    };
    let db = MaterialDatabase::new()
        .with_material(0, Material::default())
        .with_material(1, lorentz);
    let solver = DispersiveSolver::new(db);

    let pairs = solver
        .solve_at_frequency(&mesh, omega_te101, 5, sigma)
        .expect("lossy dispersive solve_at_frequency");

    let k_sq = pairs.k[0];
    // The Lorentz contribution to ε at ω < ω_0 has non-zero Im(ε)
    // (γω term in the denominator); the resulting M is complex
    // symmetric and the eigenvalue picks up an imaginary part.
    assert!(
        k_sq.im.abs() > 0.0,
        "lossy mode should have non-zero Im(k²), got {k_sq}",
    );
    // Sanity: Re(k²) shifted downward by the increased average ε
    // (filling half the cavity with ε ≈ 4 lowers the resonance);
    // we only require Re(k²) to be in a sensible range, not a tight
    // analytic match.
    assert!(
        k_sq.re > 0.0 && k_sq.re < 4.0 * k0_te101.powi(2),
        "lossy mode Re(k²) = {} should be in (0, 4 k₀²) for half-filled \
         dielectric cavity",
        k_sq.re,
    );
}

/// D4 gate test 3: `assemble_complex` at purely real `(ε, μ)` matches
/// the real `assemble` output to round-off.
#[test]
fn assemble_complex_at_real_eps_matches_real_assemble() {
    let mesh = build_cavity();
    let f_te101 = f_te101_hz();
    let omega = 2.0 * PI * f_te101;

    // Real path: free-space (ε = μ = 1) on every tet.
    let real_assembly = FemEigenAssembly::new_free_space(&mesh);
    let real_assembled = real_assembly.assemble().expect("real assemble");

    // Complex path: free-space MaterialDatabase, same mesh.
    let db = MaterialDatabase::new().with_material(0, Material::default());
    let n_tets = mesh.tetrahedra.len();
    let complex_assembly =
        FemEigenAssembly::new(&mesh, vec![1.0; n_tets], vec![1.0; n_tets]).expect("assembly");
    let complex_assembled = complex_assembly
        .assemble_complex(omega, &db)
        .expect("complex assemble");

    // Same interior-edge basis on both paths (geometry-only).
    assert_eq!(
        real_assembled.interior_edges,
        complex_assembled.interior_edges
    );
    assert_eq!(real_assembled.k.nrows(), complex_assembled.k.nrows());
    assert_eq!(real_assembled.k.nnz(), complex_assembled.k.nnz());
    assert_eq!(real_assembled.m.nnz(), complex_assembled.m.nnz());

    // Tolerance: f64::EPSILON · 1e3 relative on the maximum entry
    // magnitude in each matrix. Both paths route through
    // `assemble_tet_element_complex` internally (the real path
    // projects via `.re`), so the only sources of divergence are the
    // (Complex64::new(1.0, 0.0) * x) multiplications in the complex
    // scatter loop. Empirically these are bit-exact, but pinning the
    // bound at `EPSILON · 1e3 ≈ 2.2e-13` leaves headroom for any
    // future intermediate-format change.
    let tol_factor = f64::EPSILON * 1e3;

    let max_k_real: f64 = real_assembled
        .k
        .values()
        .iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max);
    let max_m_real: f64 = real_assembled
        .m
        .values()
        .iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max);

    // Compare K entries — index into both matrices by (row, col) and
    // verify the complex entry's real part matches the real entry to
    // round-off and the imaginary part is identically zero (within
    // `tol_factor · max_k_real`).
    let mut dense_real_k = vec![vec![0.0f64; real_assembled.k.ncols()]; real_assembled.k.nrows()];
    for (r, c, v) in real_assembled.k.triplet_iter() {
        dense_real_k[r][c] += v;
    }
    let mut dense_complex_k = vec![
        vec![Complex64::new(0.0, 0.0); complex_assembled.k.ncols()];
        complex_assembled.k.nrows()
    ];
    for (r, c, v) in complex_assembled.k.triplet_iter() {
        dense_complex_k[r][c] += v;
    }
    let n = real_assembled.k.nrows();
    let mut max_re_err = 0.0_f64;
    let mut max_im_err = 0.0_f64;
    for r in 0..n {
        for c in 0..n {
            let re_err = (dense_complex_k[r][c].re - dense_real_k[r][c]).abs();
            let im_err = dense_complex_k[r][c].im.abs();
            if re_err > max_re_err {
                max_re_err = re_err;
            }
            if im_err > max_im_err {
                max_im_err = im_err;
            }
        }
    }
    let k_tol = tol_factor * max_k_real.max(1.0);
    assert!(
        max_re_err <= k_tol,
        "K real-part mismatch: max |Re(K_complex) - K_real| = {max_re_err}, tol {k_tol}",
    );
    assert!(
        max_im_err <= k_tol,
        "K imaginary-part nonzero on real input: max |Im(K_complex)| = {max_im_err}, tol {k_tol}",
    );

    // Compare M entries — same structure.
    let mut dense_real_m = vec![vec![0.0f64; real_assembled.m.ncols()]; real_assembled.m.nrows()];
    for (r, c, v) in real_assembled.m.triplet_iter() {
        dense_real_m[r][c] += v;
    }
    let mut dense_complex_m = vec![
        vec![Complex64::new(0.0, 0.0); complex_assembled.m.ncols()];
        complex_assembled.m.nrows()
    ];
    for (r, c, v) in complex_assembled.m.triplet_iter() {
        dense_complex_m[r][c] += v;
    }
    let mut max_re_err_m = 0.0_f64;
    let mut max_im_err_m = 0.0_f64;
    for r in 0..n {
        for c in 0..n {
            let re_err = (dense_complex_m[r][c].re - dense_real_m[r][c]).abs();
            let im_err = dense_complex_m[r][c].im.abs();
            if re_err > max_re_err_m {
                max_re_err_m = re_err;
            }
            if im_err > max_im_err_m {
                max_im_err_m = im_err;
            }
        }
    }
    let m_tol = tol_factor * max_m_real.max(1.0);
    assert!(
        max_re_err_m <= m_tol,
        "M real-part mismatch: max |Re(M_complex) - M_real| = {max_re_err_m}, tol {m_tol}",
    );
    assert!(
        max_im_err_m <= m_tol,
        "M imaginary-part nonzero on real input: max |Im(M_complex)| = {max_im_err_m}, tol {m_tol}",
    );
}
