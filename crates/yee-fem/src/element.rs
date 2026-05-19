//! First-order Nedelec (Whitney-1) edge-element local matrices on a
//! single tetrahedron.
//!
//! For a physical tet with vertices `v_0..v_3`, the four scalar
//! barycentric coordinates `λ_i` are the linear nodal Lagrange basis on
//! the tet: `λ_i(v_j) = δ_{ij}` and `Σ_i λ_i ≡ 1`. The first-order
//! Nedelec (Whitney-1) basis on the six edges is (Jin, *FEM in
//! Electromagnetics* 3rd ed. §9.4, eq. 9.43)
//!
//! ```text
//!     N_{ij}(x) = λ_i(x) ∇λ_j  −  λ_j(x) ∇λ_i      for edge (i, j).
//! ```
//!
//! The curl is constant on each tet:
//!
//! ```text
//!     ∇ × N_{ij} = 2 ∇λ_i × ∇λ_j.
//! ```
//!
//! Two element matrices follow:
//!
//! * **Local stiffness** (curl-curl)
//!
//!   ```text
//!       K^e_{αβ} = (1/μ_r) ∫_T (∇×N_α) · (∇×N_β) dV
//!                = (1/μ_r) V (∇×N_α) · (∇×N_β)
//!   ```
//!
//!   exact in closed form because both curls are constant on the tet.
//!
//! * **Local mass** (vector mass)
//!
//!   ```text
//!       M^e_{αβ} = ε_r ∫_T N_α · N_β dV.
//!   ```
//!
//!   The integrand `N_α · N_β` is polynomial of degree 2 in the
//!   barycentric coordinates, so a 4-point Gauss-tet quadrature rule
//!   (exact for degree 2) integrates it exactly. We use the canonical
//!   symmetric rule with barycentric coordinates `(α, β, β, β)` and
//!   permutations, with `α = 0.585410196624969`,
//!   `β = 0.138196601125011`, each weight `w = V/4` (Jin §9.4 quadrature
//!   table; equivalently the four-point rule of Keast 1986).
//!
//! ## Sign / orientation convention
//!
//! The element-level block is emitted in **canonical local-edge
//! orientation** — edges are ordered `[01], [02], [03], [12], [13],
//! [23]` and each `N_{ij}` runs from the lower-indexed endpoint to the
//! higher-indexed endpoint. The element layer makes no attempt to
//! reconcile that with any global edge orientation; that is the
//! assembly layer's job (step T4): when a tet-local edge runs against
//! the global orientation, the corresponding row *and* column of the
//! local block are negated during scatter. Keeping the negation out of
//! the element layer means [`assemble_tet_element`] is a pure function
//! of `(vertices, ε_r, μ_r)` and produces bitwise-identical output for
//! the same tet regardless of how its neighbours' global edges happen
//! to be oriented.
//!
//! ## Phase 4.fem.eig.1 — dispersive `ε(ω)` lift
//!
//! Phase 4.fem.eig.1 (ADR-0039) extends this layer to complex scalar
//! `ε(ω)`, `μ(ω)` coefficients. The barycentric gradients, Nedelec basis
//! curls, and 4-point Gauss-tet quadrature weights are all real and
//! unchanged; only the scalar pre-multiplier on `K_local` (`1/μ`) and
//! on `M_local` (`ε`) becomes complex. This is implemented by
//! [`assemble_tet_element_complex`] which returns a
//! [`NedelecTetElementComplex`] with `SMatrix<Complex64, 6, 6>` blocks.
//! The real entry point [`assemble_tet_element`] is preserved as a thin
//! wrapper that calls the complex path with `Complex64::from` and
//! projects the result to real via `.re`. Real-valued inputs produce
//! real-valued blocks bit-for-bit identical to the Phase 4.fem.eig.0
//! shipped behaviour — see the `element_complex` integration tests.
//!
//! ## References
//!
//! * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
//!   3rd ed., Wiley 2014, §9.4 (Nedelec edge elements on tetrahedra).
//! * Keast, P., "Moderate-degree tetrahedral quadrature formulas",
//!   *Comp. Methods Appl. Mech. Eng.* 55, 1986, pp. 339–348 (degree-2
//!   four-point symmetric rule).

use nalgebra::{Matrix3, SMatrix, SVector, Vector3};
use num_complex::Complex64;

/// The six canonical local edges of a tetrahedron in
/// lower-endpoint-first order.
///
/// Index `α ∈ 0..6` maps to local edge `LOCAL_EDGES[α] = (i, j)` with
/// `i < j`. This is the row / column ordering of [`NedelecTetElement`]
/// matrices.
pub const LOCAL_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// Per-tet element matrices for the first-order Nedelec edge basis
/// (real-coefficient path — Phase 4.fem.eig.0).
///
/// Both fields are `6 × 6` real symmetric blocks indexed by
/// [`LOCAL_EDGES`]. See module docs for the exact bilinear forms.
///
/// For the Phase 4.fem.eig.1 complex-coefficient path see
/// [`NedelecTetElementComplex`].
#[derive(Debug, Clone, Copy)]
pub struct NedelecTetElement {
    /// Local curl-curl stiffness `K^e_{αβ} = (1/μ_r) V (∇×N_α)·(∇×N_β)`.
    pub k_local: SMatrix<f64, 6, 6>,
    /// Local vector mass `M^e_{αβ} = ε_r ∫_T N_α · N_β dV`.
    pub m_local: SMatrix<f64, 6, 6>,
}

/// Per-tet element matrices for the first-order Nedelec edge basis with
/// complex scalar `ε(ω)`, `μ(ω)` coefficients (Phase 4.fem.eig.1).
///
/// Both fields are `6 × 6` complex symmetric (not Hermitian — see
/// ADR-0039 / spec §11) blocks indexed by [`LOCAL_EDGES`]. The blocks
/// reduce bit-for-bit to the real [`NedelecTetElement`] blocks when
/// `eps` and `mu` are pure-real (`Im(eps) = Im(mu) = 0`). See
/// [`assemble_tet_element_complex`] for the construction.
#[derive(Debug, Clone, Copy)]
pub struct NedelecTetElementComplex {
    /// Local curl-curl stiffness
    /// `K^e_{αβ} = (1/μ(ω)) V (∇×N_α)·(∇×N_β)`. Complex symmetric.
    pub k_local: SMatrix<Complex64, 6, 6>,
    /// Local vector mass
    /// `M^e_{αβ} = ε(ω) ∫_T N_α · N_β dV`. Complex symmetric.
    pub m_local: SMatrix<Complex64, 6, 6>,
}

/// 4-point Gauss-tet quadrature, exact for polynomials of degree 2.
///
/// Each row is a barycentric quadrature point `(λ_0, λ_1, λ_2, λ_3)`;
/// the corresponding weight is `V/4` for every point (uniform).
///
/// Coefficients: `α = (5 + 3√5)/20 ≈ 0.585410196624969`,
/// `β = (5 −  √5)/20 ≈ 0.138196601125011`. Keast 1986; Jin §9.4.
const QUAD_ALPHA: f64 = 0.585_410_196_624_969;
const QUAD_BETA: f64 = 0.138_196_601_125_011;

const QUAD_POINTS: [[f64; 4]; 4] = [
    [QUAD_ALPHA, QUAD_BETA, QUAD_BETA, QUAD_BETA],
    [QUAD_BETA, QUAD_ALPHA, QUAD_BETA, QUAD_BETA],
    [QUAD_BETA, QUAD_BETA, QUAD_ALPHA, QUAD_BETA],
    [QUAD_BETA, QUAD_BETA, QUAD_BETA, QUAD_ALPHA],
];

/// Build the 4 barycentric gradients `∇λ_i` and the signed tet volume.
///
/// The barycentric coordinates `λ_i` are the unique linear functions
/// satisfying `λ_i(v_j) = δ_{ij}` and `Σ_i λ_i ≡ 1`. Writing
/// `λ_i(x) = a_i + b_i x + c_i y + d_i z` and stacking the four
/// constraints into the `4 × 4` system
///
/// ```text
///     [1 v_0^x v_0^y v_0^z] [a_0 a_1 a_2 a_3]     [1 0 0 0]
///     [1 v_1^x v_1^y v_1^z] [b_0 b_1 b_2 b_3]  =  [0 1 0 0]
///     [1 v_2^x v_2^y v_2^z] [c_0 c_1 c_2 c_3]     [0 0 1 0]
///     [1 v_3^x v_3^y v_3^z] [d_0 d_1 d_2 d_3]     [0 0 0 1]
/// ```
///
/// gives `∇λ_i = (b_i, c_i, d_i)`. The signed volume is
/// `V = det(M) / 6` where `M` is the left-hand `4 × 4`.
///
/// Returned: `(grads, signed_volume)` where `grads[i]` is `∇λ_i`.
fn barycentric_gradients_and_volume(vertices: &[Vector3<f64>; 4]) -> ([Vector3<f64>; 4], f64) {
    // Signed volume V = (1/6) (v1 − v0) · ((v2 − v0) × (v3 − v0)).
    let e1 = vertices[1] - vertices[0];
    let e2 = vertices[2] - vertices[0];
    let e3 = vertices[3] - vertices[0];
    let signed_volume = e1.dot(&e2.cross(&e3)) / 6.0;

    // Closed-form face-normal formulas (Jin §9.4):
    //
    //     ∇λ_0 = (1/(6V)) ( (v_2 − v_1) × (v_3 − v_1) )  (face 1-2-3)
    //
    // and analogously for the other three vertices, with the sign of
    // the cross product chosen so that ∇λ_i points *toward* v_i. The
    // sign convention below is the one used in Jin eq. 9.42 — the face
    // opposite vertex i is traversed in CCW order seen from v_i, and
    // the result is divided by 6V (not 3V; the factor of 2 comes from
    // ½ |cross| = face area).
    //
    // We use a vectorised formulation: for each vertex i, ∇λ_i =
    // (v_{i+1} − v_{i+3}) × (v_{i+2} − v_{i+3}) / (6V) with indices mod
    // 4 and the sign flipped for odd permutations of (1,2,3,0). To
    // avoid sign-table bookkeeping we instead solve the 3×3 linear
    // system that defines the gradients directly.
    //
    // The 3×3 system uses any three vertices as the row basis:
    // (v_j − v_0) · ∇λ_i = δ_{ij} − δ_{i0}, j ∈ {1, 2, 3}.
    // Equivalently the matrix [e1 e2 e3]^T (rows) has columns equal to
    // ∇λ_1, ∇λ_2, ∇λ_3 when right-inverted; ∇λ_0 = −Σ_{i>0} ∇λ_i by the
    // partition-of-unity identity Σ_i ∇λ_i = ∇(1) = 0.
    let rows = Matrix3::from_rows(&[e1.transpose(), e2.transpose(), e3.transpose()]);
    let rows_inv = rows
        .try_inverse()
        .expect("degenerate tet: edge vectors are linearly dependent");
    let grad_1 = rows_inv.column(0).into_owned();
    let grad_2 = rows_inv.column(1).into_owned();
    let grad_3 = rows_inv.column(2).into_owned();
    let grad_0 = -(grad_1 + grad_2 + grad_3);

    ([grad_0, grad_1, grad_2, grad_3], signed_volume)
}

/// Assemble the `6 × 6` Nedelec local stiffness + mass block for a
/// single tetrahedron with **complex scalar** `ε(ω)`, `μ(ω)`
/// (Phase 4.fem.eig.1).
///
/// `vertices` must be ordered so that the signed volume is positive
/// (`yee-mesh::TetMesh3D::new` enforces this; callers passing
/// hand-rolled vertices should ensure the same). The returned matrices
/// are emitted in canonical local-edge orientation per [`LOCAL_EDGES`]
/// — sign flips against global-edge orientation are the assembly
/// layer's job.
///
/// The local matrices are
///
/// ```text
///     K^e_{αβ} = (1/μ(ω)) · V · (∇×N_α) · (∇×N_β)
///     M^e_{αβ} = ε(ω) · ∫_T N_α · N_β dV
/// ```
///
/// with the same barycentric / Nedelec / 4-point Gauss-tet machinery as
/// the real-coefficient path; only the scalar pre-multiplier becomes
/// complex. For real `eps` and `mu` (imaginary part exactly zero) the
/// returned blocks reduce bit-for-bit to the real
/// [`assemble_tet_element`] output — this is the load-bearing
/// backward-compatibility invariant per ADR-0039 §6.
///
/// # Panics
///
/// Panics if the tet is degenerate (all four vertices coplanar; the
/// barycentric gradient system is singular). Real meshes go through
/// [`yee_mesh::TetMesh3D::new`] which rejects such tets at
/// construction; this entry point trusts its caller.
pub fn assemble_tet_element_complex(
    vertices: [Vector3<f64>; 4],
    eps: Complex64,
    mu: Complex64,
) -> NedelecTetElementComplex {
    let (grads, signed_volume) = barycentric_gradients_and_volume(&vertices);
    let volume = signed_volume.abs();

    // ---- Local stiffness: K^e_{αβ} = (1/μ(ω)) V (∇×N_α)·(∇×N_β) ----
    //
    // ∇×N_{ij} = 2 ∇λ_i × ∇λ_j is constant on the tet and **real**, so
    // the integral is exactly `volume × (curl_α · curl_β)`. Only the
    // outer `(1/μ)` scalar carries the complex frequency dependence.
    let mut curls = [Vector3::<f64>::zeros(); 6];
    for (alpha, &(i, j)) in LOCAL_EDGES.iter().enumerate() {
        curls[alpha] = 2.0 * grads[i].cross(&grads[j]);
    }

    let mut k_local = SMatrix::<Complex64, 6, 6>::zeros();
    let inv_mu = Complex64::new(1.0, 0.0) / mu;
    for alpha in 0..6 {
        for beta in 0..6 {
            let real_entry = volume * curls[alpha].dot(&curls[beta]);
            k_local[(alpha, beta)] = inv_mu * Complex64::new(real_entry, 0.0);
        }
    }

    // ---- Local mass: M^e_{αβ} = ε(ω) ∫_T N_α · N_β dV --------------
    //
    // 4-point Gauss-tet quadrature, exact for polynomials of degree 2.
    // The basis vectors `N_α(λ)` and their dot products are real; only
    // the outer `ε(ω)` scalar carries the complex frequency
    // dependence. Accumulate the real integral first, then multiply.
    let weight = volume / 4.0;
    let mut m_real = SMatrix::<f64, 6, 6>::zeros();

    for qp in &QUAD_POINTS {
        // Evaluate all six basis vectors at this quadrature point.
        let mut basis = [Vector3::<f64>::zeros(); 6];
        for (alpha, &(i, j)) in LOCAL_EDGES.iter().enumerate() {
            basis[alpha] = qp[i] * grads[j] - qp[j] * grads[i];
        }
        for alpha in 0..6 {
            for beta in 0..6 {
                m_real[(alpha, beta)] += weight * basis[alpha].dot(&basis[beta]);
            }
        }
    }

    let mut m_local = SMatrix::<Complex64, 6, 6>::zeros();
    for alpha in 0..6 {
        for beta in 0..6 {
            m_local[(alpha, beta)] = eps * Complex64::new(m_real[(alpha, beta)], 0.0);
        }
    }

    NedelecTetElementComplex { k_local, m_local }
}

/// Assemble the per-face 1st-order Engquist–Majda ABC contribution
/// (Phase 4.fem.eig.2 step E1).
///
/// On an ABC-tagged exterior triangular face with outward normal `n̂`,
/// the Engquist–Majda 1977 radiation condition
///
/// ```text
///     n̂ × ∇×E   =   −j k₀  n̂ × (n̂ × E)
/// ```
///
/// substituted into the curl-curl variational form's surface integral
/// yields the per-face stiffness contribution (Jin, *FEM in
/// Electromagnetics* 3rd ed. §10.4, eq. 10.28)
///
/// ```text
///     K_ABC^{e,face}_{ij}  =  + j k₀ · (1/μ_r,face) · ∫_face
///                                (n̂ × N_i) · (n̂ × N_j)  dS.
/// ```
///
/// For the first-order Nedelec / Whitney-1 face basis the basis
/// vector `N_i` restricted to the triangular face reduces to a constant
/// edge tangent `t_i = v_{(i+1) mod 3} − v_i` (the dual of the edge in
/// the Whitney complex), so `n̂ × N_i` is constant over the face and the
/// surface integral evaluates exactly to
///
/// ```text
///     ∫_face (n̂ × N_i) · (n̂ × N_j) dS = A · (n̂ × t_i) · (n̂ × t_j)
/// ```
///
/// where `A = 0.5 · ||t_0 × t_1||` is the triangle area. The returned
/// `3 × 3` block is therefore
///
/// ```text
///     B[i][j] = j · k₀ · (A / μ_r,face) · (n̂ × t_i) · (n̂ × t_j),
/// ```
///
/// indexed by the three face edges `(0→1, 1→2, 2→0)` in the canonical
/// CCW traversal of the face vertices. The block is **complex-symmetric**
/// (`B == B^T`, NOT Hermitian) because the imaginary prefactor `j k₀` is
/// scalar and the real `(n̂ × t_i) · (n̂ × t_j)` Gram form is symmetric.
/// Adding ABC face contributions promotes the otherwise-real
/// closed-cavity stiffness matrix to complex-symmetric — the same
/// mathematical fact that lets the ABC absorb outgoing waves (the
/// imaginary part carries the radiation resistance).
///
/// ## Sign / orientation convention
///
/// The block is emitted in **canonical local-edge orientation** —
/// each edge `i` runs from `face_vertices[i]` to
/// `face_vertices[(i + 1) % 3]`. The element layer makes no attempt to
/// reconcile that with any global edge orientation; that is the
/// assembly layer's job (Phase 4.fem.eig.2 step E3): when a face-local
/// edge runs against the global orientation, the corresponding row
/// *and* column of this block are negated during scatter. Keeping the
/// negation out of the element layer means
/// [`assemble_abc_face_block`] is a pure function of
/// `(face_vertices, outward_normal, k0, mu_r_face)`.
///
/// ## References
///
/// * Engquist, B. and Majda, A., "Absorbing boundary conditions for the
///   numerical simulation of waves", *Math. Comp.* 31 (1977),
///   pp. 629–651 — canonical 1st-order ABC derivation.
/// * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
///   3rd ed., Wiley 2014, §10.4 (ABC face contributions).
/// * Phase 4.fem.eig.2 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
///   §4.2 — the bilinear form this helper implements.
pub fn assemble_abc_face_block(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    k0: f64,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3> {
    // Three edge tangents of the triangular face, in CCW order:
    //     t_i = v_{(i+1) mod 3} − v_i.
    let t = [
        face_vertices[1] - face_vertices[0],
        face_vertices[2] - face_vertices[1],
        face_vertices[0] - face_vertices[2],
    ];

    // Face area A = 0.5 · ||t_0 × t_1||. The cross product of any two
    // edge tangents of a planar triangle yields a vector of magnitude
    // 2 · A (twice the triangle area), so dividing by 2 recovers A.
    let face_area = 0.5 * t[0].cross(&t[1]).norm();

    // (n̂ × N_i) is constant per face for the Whitney-1 face basis;
    // for the edge-tangent dual basis this is exactly (n̂ × t_i).
    let n_cross_t = [
        outward_normal.cross(&t[0]),
        outward_normal.cross(&t[1]),
        outward_normal.cross(&t[2]),
    ];

    // Outer prefactor: j · k₀ · (A / μ_r,face). Purely imaginary.
    let prefactor = Complex64::new(0.0, 1.0) * Complex64::new(k0 * face_area / mu_r_face, 0.0);

    let mut block = SMatrix::<Complex64, 3, 3>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            let gram_entry = n_cross_t[i].dot(&n_cross_t[j]);
            block[(i, j)] = prefactor * Complex64::new(gram_entry, 0.0);
        }
    }

    block
}

/// Assemble the per-face modal wave-port stiffness contribution
/// (Phase 4.fem.eig.2 step E2).
///
/// On a wave-port-tagged exterior triangular face with outward normal
/// `n̂`, the modal wave-port boundary condition (Jin, *FEM in
/// Electromagnetics* 3rd ed. §10.5; Pozar, *Microwave Engineering* 4th
/// ed. §3.3) contributes a per-face stiffness term
///
/// ```text
///     K_port^{e,face}_{ij}  =  + j β_mode · (1/μ_r,face) · ∫_face
///                                  (n̂ × N_i) · (n̂ × N_j)  dS
/// ```
///
/// to the global complex stiffness matrix. The structure is identical
/// to [`assemble_abc_face_block`] with the wave-port modal propagation
/// constant `β_mode` replacing the free-space wavenumber `k₀`. For the
/// dominant TE_{10} mode of a rectangular waveguide,
/// `β_mode = sqrt(k₀² ε_r μ_r − (π/a)²)`; below cutoff `β_mode` is
/// purely imaginary and the caller may decide to skip the assembly.
/// `β_mode` is computed externally by the caller (typically from a
/// `NumericalCrossSection` eigensolver dispatched per swept frequency).
///
/// As in the ABC case, the Whitney-1 face basis `N_i` restricted to the
/// triangular face reduces to the constant edge tangent
/// `t_i = v_{(i+1) mod 3} − v_i`, so `n̂ × N_i = n̂ × t_i` is constant
/// per face and the surface integral evaluates exactly to
///
/// ```text
///     B[i][j] = j · β_mode · (A / μ_r,face) · (n̂ × t_i) · (n̂ × t_j),
/// ```
///
/// where `A = 0.5 · ||t_0 × t_1||` is the triangle area. The returned
/// block is **complex-symmetric** (`B == B^T`, NOT Hermitian) — the
/// imaginary prefactor `j β_mode` is scalar and the real
/// `(n̂ × t_i) · (n̂ × t_j)` Gram form is symmetric. When
/// `β_mode = 0` (cutoff) every entry is identically zero.
///
/// ## Sign / orientation convention
///
/// As in [`assemble_abc_face_block`], the block is emitted in
/// **canonical local-edge orientation** — each edge `i` runs from
/// `face_vertices[i]` to `face_vertices[(i + 1) % 3]`. Local-to-global
/// orientation flips are the assembly layer's job (Phase 4.fem.eig.2
/// step E3); this element helper is a pure function of
/// `(face_vertices, outward_normal, beta_mode, mu_r_face)`.
///
/// ## References
///
/// * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
///   3rd ed., Wiley 2014, §10.5 (wave-port modal decomposition).
/// * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012,
///   §3.3 (waveguide TE/TM modes, propagation constants).
/// * Phase 4.fem.eig.2 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
///   §4.3 — the bilinear form this helper implements.
pub fn assemble_port_face_block(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_mode: f64,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3> {
    // Three edge tangents of the triangular face, in CCW order:
    //     t_i = v_{(i+1) mod 3} − v_i.
    let t = [
        face_vertices[1] - face_vertices[0],
        face_vertices[2] - face_vertices[1],
        face_vertices[0] - face_vertices[2],
    ];

    // Face area A = 0.5 · ||t_0 × t_1||.
    let face_area = 0.5 * t[0].cross(&t[1]).norm();

    // (n̂ × N_i) is constant per face for the Whitney-1 face basis;
    // for the edge-tangent dual basis this is exactly (n̂ × t_i).
    let n_cross_t = [
        outward_normal.cross(&t[0]),
        outward_normal.cross(&t[1]),
        outward_normal.cross(&t[2]),
    ];

    // Outer prefactor: j · β_mode · (A / μ_r,face). Purely imaginary.
    // When β_mode = 0 (cutoff), the prefactor is zero and the block
    // vanishes identically — no special-case branch needed.
    let prefactor =
        Complex64::new(0.0, 1.0) * Complex64::new(beta_mode * face_area / mu_r_face, 0.0);

    let mut block = SMatrix::<Complex64, 3, 3>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            let gram_entry = n_cross_t[i].dot(&n_cross_t[j]);
            block[(i, j)] = prefactor * Complex64::new(gram_entry, 0.0);
        }
    }

    block
}

/// Assemble the per-face modal wave-port right-hand-side contribution
/// (Phase 4.fem.eig.2 step E2).
///
/// On a wave-port-tagged exterior triangular face, the driven FEM
/// system carries a per-face right-hand-side contribution encoding the
/// incident modal current (Jin, *FEM in Electromagnetics* 3rd ed.
/// §10.5, eq. 10.74; Pozar, *Microwave Engineering* 4th ed. §3.3):
///
/// ```text
///     b_port,i  =  + 2 j β_mode  ·  ∫_face  N_i · E_t_mode  dS.
/// ```
///
/// The leading factor of `2` is the matched-port double-amplitude
/// convention from Pozar §3.3: for a matched port driving the incident
/// mode at amplitude `a_inc = 1`, the total tangential E-field at the
/// port boundary is `E_inc + E_refl = 2 · E_inc` because a perfectly
/// matched modal termination absorbs the outgoing wave but the boundary
/// itself sees twice the incident amplitude. Any modulation by the
/// caller's incident amplitude is folded into the supplied
/// `mode_e_t_at_centroid` (i.e. the caller passes
/// `a_inc · e_mode(x_c, y_c)`).
///
/// For the first-order Whitney-1 face basis the basis vector
/// `N_i` restricted to the triangular face is treated under the same
/// lumped edge-tangent approximation as the ABC and wave-port
/// face-block helpers above — `N_i|_face ≈ t_i / ||t_i||² · ||t_i||`
/// in the dual sense — so the integrand `N_i · E_t_mode` is
/// approximated by its face-centroid sample. The face-centroid
/// quadrature evaluates to
///
/// ```text
///     ∫_face N_i · E_t_mode dS  ≈  (A / 3) · (t_i · E_t_mode),
/// ```
///
/// where the `1/3` factor is the Whitney-1 lumped edge basis weight at
/// the face centroid. Substituting:
///
/// ```text
///     b_i  =  2 j β_mode · (A / 3) · (t_i · E_t_mode).
/// ```
///
/// ## CCCCCCCCC normalisation note
///
/// The lumped `t_i / 3` weighting is **not** the exact Whitney-1
/// basis-at-centroid identity
/// `N_i(centroid) = (1/3) · (∇λ_b − ∇λ_a)`. The lumped form is paired
/// with the dual approximation in
/// [`crate::open_boundary::OpenBoundarySolver::extract_s11`]'s
/// `e_t_at_face_centroid` so the round-trip
/// modal-RHS-then-modal-projection cancellation is preserved at the
/// lumped level. The CCCCCCCCC scaling fix is in `extract_s11`, which
/// divides the inner product by the modal self-inner-product computed
/// via the same lumped quadrature. A future Phase 4.fem.eig.2.0.1
/// refinement (ADR-0040 §C-3) will lift this RHS, the centroid
/// reconstruction, and the per-Gauss-point modal sampling to the exact
/// Whitney basis identity in a single coupled change.
///
/// The returned `SVector<Complex64, 3>` is indexed by the three face
/// edges `(0→1, 1→2, 2→0)` in the canonical CCW traversal of the face
/// vertices; the assembly layer (Phase 4.fem.eig.2 step E3) is
/// responsible for the local-to-global orientation flips per shared
/// edge.
///
/// `mode_e_t_at_centroid` is the tangential E-field of the incident
/// mode at the face centroid (typically `a_inc · e_mode(x_c, y_c)`
/// where `e_mode` is sourced from
/// `yee_mom::eigensolver::NumericalCrossSection::e_tangential_at`).
/// Its component along the face normal is dropped — only the
/// tangential projection contributes via the dot product with `t_i`,
/// which lies in the face plane.
///
/// ## Sign / orientation convention
///
/// As with the face-block helpers above, the RHS is emitted in
/// **canonical local-edge orientation** — each edge `i` runs from
/// `face_vertices[i]` to `face_vertices[(i + 1) % 3]`. The assembly
/// layer applies the local-to-global sign flip during scatter.
///
/// ## References
///
/// * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
///   3rd ed., Wiley 2014, §10.5 (wave-port modal decomposition),
///   eq. 10.74 (incident-wave RHS).
/// * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012,
///   §3.3 — the matched-port `E_inc + E_refl = 2 · E_inc` convention
///   that motivates the factor-of-`2` prefactor.
/// * Phase 4.fem.eig.2 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
///   §4.3 — the wave-port forcing bilinear form this helper implements.
pub fn assemble_port_modal_rhs(
    face_vertices: [Vector3<f64>; 3],
    _outward_normal: Vector3<f64>,
    beta_mode: f64,
    mode_e_t_at_centroid: Vector3<f64>,
) -> SVector<Complex64, 3> {
    // Three edge tangents of the triangular face, in CCW order:
    //     t_i = v_{(i+1) mod 3} − v_i.
    let t = [
        face_vertices[1] - face_vertices[0],
        face_vertices[2] - face_vertices[1],
        face_vertices[0] - face_vertices[2],
    ];

    // Face area A = 0.5 · ||t_0 × t_1||.
    let face_area = 0.5 * t[0].cross(&t[1]).norm();

    // Outer prefactor: 2 j β_mode · (A / 3). Purely imaginary.
    let prefactor =
        Complex64::new(0.0, 1.0) * Complex64::new(2.0 * beta_mode * face_area / 3.0, 0.0);

    let mut rhs = SVector::<Complex64, 3>::zeros();
    for i in 0..3 {
        let dot = t[i].dot(&mode_e_t_at_centroid);
        rhs[i] = prefactor * Complex64::new(dot, 0.0);
    }

    rhs
}

/// Three-point Gauss-quadrature points on the reference triangle in
/// barycentric coordinates. Each point carries weight `A / 3` where `A`
/// is the triangle area. Together they integrate polynomials up to
/// degree 2 exactly on the reference triangle (Strang & Fix 1973;
/// equivalent to the "second-order" rule in Cowper 1973).
///
/// Each row is `(λ_0, λ_1, λ_2)` for one Gauss point. The three points
/// `(2/3, 1/6, 1/6)`, `(1/6, 2/3, 1/6)`, `(1/6, 1/6, 2/3)` are the
/// canonical permutation-symmetric set placed at the edge midpoints'
/// reflections.
const TRI_GAUSS_3PT_BARY: [[f64; 3]; 3] = [
    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
];

/// Compute the three in-plane barycentric gradients `∇λ_a, ∇λ_b, ∇λ_c`
/// for a triangular face in 3-space and the face area `A`.
///
/// For a triangle with vertices `(v_0, v_1, v_2)` in CCW order seen from
/// the outward-normal side, the in-plane barycentric coordinate gradient
/// is (Bossavit 1988; Jin §8.4)
///
/// ```text
///     ∇λ_a = (v_b − v_c) × n̂ / (2 A),
/// ```
///
/// where `(a, b, c)` is a cyclic permutation of `(0, 1, 2)` and
/// `n̂` is the outward unit normal. The gradient lies in the face plane
/// (`∇λ_a · n̂ = 0`) and points toward `v_a` (in the half-plane bounded
/// by edge `b → c`). Each `∇λ_a` is constant across the face because the
/// barycentric coordinates are linear in space.
///
/// Returns `(grads, area)` where `grads[a]` is `∇λ_a` and
/// `area = 0.5 · ||(v_1 − v_0) × (v_2 − v_0)||`.
fn face_barycentric_gradients_and_area(
    face_vertices: &[Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
) -> ([Vector3<f64>; 3], f64) {
    let v0 = face_vertices[0];
    let v1 = face_vertices[1];
    let v2 = face_vertices[2];

    // Face area from the cross product of any two edges.
    let face_area = 0.5 * (v1 - v0).cross(&(v2 - v0)).norm();

    // Normalise the outward normal (caller may supply a non-unit vector).
    let n_norm = outward_normal.norm();
    let n_hat = if n_norm > 0.0 {
        outward_normal / n_norm
    } else {
        outward_normal
    };

    let inv_two_a = if face_area > 0.0 {
        1.0 / (2.0 * face_area)
    } else {
        0.0
    };

    // ∇λ_a = (v_b − v_c) × n̂ / (2 A) with (a, b, c) cyclic.
    let grad_0 = (v1 - v2).cross(&n_hat) * inv_two_a;
    let grad_1 = (v2 - v0).cross(&n_hat) * inv_two_a;
    let grad_2 = (v0 - v1).cross(&n_hat) * inv_two_a;

    ([grad_0, grad_1, grad_2], face_area)
}

/// Assemble the per-face wave-port stiffness block evaluated at
/// 3-point Gauss quadrature with the **exact Whitney-1 basis**
/// (Phase 4.fem.eig.3 F1).
///
/// This is the coupled-Whitney upgrade of
/// [`assemble_port_face_block`], lifting the lumped
/// `N_i ≈ t_i` proxy to the exact Whitney-1 identity
/// `N_i(ξ) = λ_a(ξ) ∇λ_b − λ_b(ξ) ∇λ_a` evaluated at the three Gauss
/// points
///
/// ```text
///     ξ_g ∈ { (2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 1/6, 2/3) }
/// ```
///
/// in barycentric coordinates on the reference triangle (each weighted
/// `w_g = A / 3`). The 3×3 block entries are
///
/// ```text
///     B[i][j]  =  j · β_mode · (1/μ_r,face) · Σ_g  w_g
///                   · (n̂ × N_i(ξ_g)) · (n̂ × N_j(ξ_g)).
/// ```
///
/// Indexed by the three directed edges `i = 0, 1, 2` with endpoints
/// `(a, b) = (i, (i+1) mod 3)` in CCW order — same canonical local-edge
/// orientation as [`assemble_port_face_block`]; the assembly layer
/// applies any local-to-global sign flip during scatter.
///
/// The block is **complex-symmetric** (`B == B^T`, NOT Hermitian)
/// because the imaginary prefactor `j β_mode` is scalar and the real
/// `(n̂ × N_i) · (n̂ × N_j)` Gram form is symmetric. When `β_mode = 0`
/// (modal cutoff) every entry is identically zero.
///
/// `β_mode` is taken as `Complex64` for full generality (e.g. below-cutoff
/// evanescent regime where `β_mode = j α`); the v2 entry point with
/// `f64` propagation constant lifts to `Complex64::new(β, 0)` at the
/// caller boundary.
///
/// ## Why this differs from the lumped centroid path
///
/// At the face centroid `ξ_c = (1/3, 1/3, 1/3)` the exact Whitney-1
/// identity gives `N_i(ξ_c) = (1/3)(∇λ_b − ∇λ_a)`, **not** the lumped
/// proxy `t_i / 3` used by [`assemble_port_face_block`]. The two
/// vectors agree only on an equilateral triangle; on every Kuhn-
/// decomposed face the lumped proxy mis-evaluates `N_i(centroid)` and
/// drives the round-trip modal-RHS-then-projection cancellation away
/// from the Pozar §3.3 matched-port identity. F1's coupled fix
/// (paired with [`assemble_port_face_rhs_gauss_pts`] in the RHS and the
/// `e_t_at_face_gauss_pts` projection helper in `OpenBoundarySolver`)
/// preserves the round-trip identity at the exact-basis level.
///
/// ## References
///
/// * Bossavit, A., "Whitney forms: a class of finite elements for
///   three-dimensional computations in electromagnetism", *IEE Proc.*
///   135-A (1988), pp. 493–500.
/// * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
///   3rd ed., Wiley 2014, §10.5 (wave-port modal decomposition).
/// * Cowper, G. R., "Gaussian quadrature formulas for triangles",
///   *Int. J. Numer. Meth. Eng.* 7 (1973), pp. 405–408 — the 3-point
///   rule used here.
/// * Phase 4.fem.eig.3 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
///   §4.1 — the bilinear form this helper implements.
pub fn assemble_port_face_block_gauss_pts(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_mode: Complex64,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3> {
    let (grads, face_area) = face_barycentric_gradients_and_area(&face_vertices, outward_normal);

    // Normalise outward normal once for the (n̂ × ·) operations below.
    let n_norm = outward_normal.norm();
    let n_hat = if n_norm > 0.0 {
        outward_normal / n_norm
    } else {
        outward_normal
    };

    // Outer prefactor: j · β_mode · (1 / μ_r,face). Complex.
    let prefactor = Complex64::new(0.0, 1.0) * beta_mode * Complex64::new(1.0 / mu_r_face, 0.0);

    // Per-Gauss-point quadrature weight w_g = A / 3.
    let w_g = face_area / 3.0;

    let mut block = SMatrix::<Complex64, 3, 3>::zeros();

    for bary in &TRI_GAUSS_3PT_BARY {
        // Evaluate all three Whitney-1 edge basis functions at this
        // Gauss point. Edge i runs from vertex a = i to vertex
        // b = (i + 1) mod 3, so N_i(ξ) = λ_a ∇λ_b − λ_b ∇λ_a.
        let mut basis = [Vector3::<f64>::zeros(); 3];
        for (i, basis_i) in basis.iter_mut().enumerate() {
            let a = i;
            let b = (i + 1) % 3;
            *basis_i = bary[a] * grads[b] - bary[b] * grads[a];
        }

        // Pre-compute (n̂ × N_i) per edge.
        let n_cross_n = [
            n_hat.cross(&basis[0]),
            n_hat.cross(&basis[1]),
            n_hat.cross(&basis[2]),
        ];

        for i in 0..3 {
            for j in 0..3 {
                let gram_entry = n_cross_n[i].dot(&n_cross_n[j]);
                block[(i, j)] += prefactor * Complex64::new(w_g * gram_entry, 0.0);
            }
        }
    }

    block
}

/// Assemble the per-face wave-port right-hand-side contribution at
/// 3-point Gauss quadrature with the **exact Whitney-1 basis**
/// (Phase 4.fem.eig.3 F1).
///
/// Companion of [`assemble_port_face_block_gauss_pts`]. The caller
/// pre-evaluates the modal tangential E-field at the three Gauss
/// points on the reference triangle (typically by sampling a
/// `NumericalCrossSection::e_tangential_at` or evaluating an analytic
/// modal profile at the corresponding world-space points). The RHS
/// entries are
///
/// ```text
///     b_i  =  2 j β_mode  ·  Σ_g  w_g · N_i(ξ_g) · E_t_mode(ξ_g),
/// ```
///
/// with `w_g = A / 3` and the same Whitney-1 identity
/// `N_i(ξ) = λ_a(ξ) ∇λ_b − λ_b(ξ) ∇λ_a` as the stiffness block. The
/// factor of `2` is the matched-port double-amplitude convention from
/// Pozar §3.3 (the boundary sees `E_inc + E_refl = 2 · E_inc` at a
/// matched termination).
///
/// `modal_e_t_at_gauss_pts[g]` is the **tangential** incident-mode
/// E-field at the world-space point corresponding to barycentric Gauss
/// point `g` (already scaled by the caller's incident amplitude
/// `a_inc`). Any out-of-plane component is dropped by the dot product
/// with `N_i(ξ_g)` (which lies in the face plane by construction).
///
/// The returned `SVector<Complex64, 3>` is indexed by the three face
/// edges `(0 → 1, 1 → 2, 2 → 0)` in canonical CCW traversal — same
/// orientation convention as [`assemble_port_modal_rhs`]. The assembly
/// layer applies any local-to-global sign flip during scatter.
///
/// `β_mode` is taken as `Complex64` for the same reason as the
/// stiffness block helper above.
///
/// ## References
///
/// * Phase 4.fem.eig.3 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
///   §4.1 — the RHS bilinear form this helper implements.
/// * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §3.3.
pub fn assemble_port_face_rhs_gauss_pts(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_mode: Complex64,
    modal_e_t_at_gauss_pts: [Vector3<f64>; 3],
) -> SVector<Complex64, 3> {
    let (grads, face_area) = face_barycentric_gradients_and_area(&face_vertices, outward_normal);

    // Outer prefactor: 2 j β_mode. Complex.
    let prefactor = Complex64::new(0.0, 2.0) * beta_mode;

    // Per-Gauss-point quadrature weight w_g = A / 3.
    let w_g = face_area / 3.0;

    let mut rhs = SVector::<Complex64, 3>::zeros();

    for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
        // Whitney-1 edge basis at this Gauss point — identical
        // identity as in the stiffness block above.
        let mut basis = [Vector3::<f64>::zeros(); 3];
        for (i, basis_i) in basis.iter_mut().enumerate() {
            let a = i;
            let b = (i + 1) % 3;
            *basis_i = bary[a] * grads[b] - bary[b] * grads[a];
        }

        let e_t_g = modal_e_t_at_gauss_pts[g];
        for i in 0..3 {
            let dot = basis[i].dot(&e_t_g);
            rhs[i] += prefactor * Complex64::new(w_g * dot, 0.0);
        }
    }

    rhs
}

/// Assemble the per-face 2nd-order Engquist–Majda ABC contribution
/// (Phase 4.fem.eig.3 step F3).
///
/// The 2nd-order Engquist–Majda radiation condition on a planar surface
/// with outward normal `n̂` is (Engquist & Majda 1979, *IEEE Trans.
/// Antennas Propag.* 27(5) p. 661, eq. 9; equivalent forms in Jin §10.4)
///
/// ```text
///     n̂ × ∇×E   =   −j k₀ · n̂×(n̂×E)   +   (1/(2 j k₀)) · ∇_t × (∇_t × E_t).
/// ```
///
/// Substituting into the curl-curl variational form, the bilinear form
/// picks up **two** boundary terms per face — the 1st-order Mur term
/// inherited from [`assemble_abc_face_block`] and a new tangential-curl
/// correction (Jin §10.4 / Engquist–Majda 1979 eq. 9):
///
/// ```text
///     a_ABC2(E, v)  =  + j k₀ · ∫_face (n̂×N_i)·(n̂×N_j) dS    ← 1st-order Mur
///                      − (1/(2 k₀)) · ∫_face (∇_t×N_i)(∇_t×N_j) dS  ← 2nd-order correction
/// ```
///
/// where `∇_t × N_i = n̂ · (∇ × N_i)` is the **scalar** tangential curl
/// (the out-of-plane component of the 3D curl of the in-plane Whitney
/// basis). For first-order Whitney-1 elements on a triangular face,
/// `∇ × N_i = 2 ∇λ_a × ∇λ_b` is parallel to `n̂` (because the in-plane
/// barycentric gradients `∇λ_a, ∇λ_b` both lie in the face plane and
/// their cross product is normal to the plane), so the scalar
/// tangential curl `∇_t × N_i = 2 n̂ · (∇λ_a × ∇λ_b)` is **constant per
/// face**. The curl-correction surface integral is therefore exact and
/// reduces to a rank-1 real-symmetric outer product. The returned
/// `3 × 3` block is
///
/// ```text
///     B[i][j]  =  + j · k₀ · (A / μ_r,face) · R_1[i][j]
///                 − (1 / (2 k₀)) · (A / μ_r,face) · R_2[i][j],
/// ```
///
/// where
///
/// ```text
///     R_1[i][j]  =  (n̂ × t_i) · (n̂ × t_j),       t_i = v_{(i+1) mod 3} − v_i,
///     R_2[i][j]  =  c_i · c_j,                    c_i = 2 · n̂ · (∇λ_a × ∇λ_b),
///                                                 (a, b) = (i, (i+1) mod 3).
/// ```
///
/// The 1st-order term `R_1` here is identical to the Gram form used by
/// [`assemble_abc_face_block`] — the lumped Whitney-1 edge-tangent dual
/// identity `N_i|_face = t_i` (Bossavit) is preserved on the 1st-order
/// part so that `AbcOrder::First` and the imaginary part of
/// `AbcOrder::Second` agree bit-for-bit when the 2nd-order term is
/// dropped.
///
/// ## Note on the surface-curl reduction
///
/// The spec design document writes the 2nd-order Gram form as
/// `(n̂ × ∇×N_i) · (n̂ × ∇×N_j)`; on a planar face `∇ × N_i` is parallel
/// to `n̂`, so the literal `n̂ × ∇×N_i` term is zero. The correct
/// physical reduction — the scalar surface curl `∇_t × N_i = n̂·(∇×N_i)`
/// — yields the rank-1 outer product `R_2[i][j] = c_i · c_j` documented
/// above. This is the form Engquist–Majda 1979 eq. 9 and Jin §10.4
/// derive for first-order Whitney-1 face elements.
///
/// ## Sign / orientation convention
///
/// As with [`assemble_abc_face_block`] and the wave-port face helpers,
/// the block is emitted in **canonical local-edge orientation** — each
/// edge `i` runs from `face_vertices[i]` to
/// `face_vertices[(i + 1) % 3]`. Local-to-global orientation flips are
/// the assembly layer's job; this element-layer helper is a pure
/// function of `(face_vertices, outward_normal, k0, mu_r_face)`.
///
/// ## Block symmetry
///
/// The block is **complex-symmetric** (`B == B^T`, NOT Hermitian): both
/// `R_1` and `R_2` are real-symmetric Gram matrices, and the
/// prefactors `j k₀ · A / μ_r` (imaginary) and `−A / (2 k₀ μ_r)` (real)
/// are scalars. The composite has Im ≠ 0 (from `R_1`) and Re ≠ 0 (from
/// `R_2`).
///
/// ## Frequency scaling
///
/// At high `k₀` the curl correction `−(1/(2 k₀)) R_2` is suppressed by a
/// factor `1/k₀` relative to `+ j k₀ R_1`; at low `k₀` it diverges as
/// `1/k₀`. The Engquist–Majda derivation is asymptotic for
/// `k₀ ≫ k_grazing`; below cutoff the 1st-order ABC is the numerically
/// stable choice (see spec §10 mitigation for the WR-90 band-edge
/// behaviour).
///
/// ## References
///
/// * Engquist, B. and Majda, A., "Radiation boundary conditions for the
///   numerical simulation of waves", *Math. Comp.* 31 (1977),
///   pp. 629–651; and *IEEE Trans. Antennas Propag.* 27(5) (1979)
///   p. 661, eq. 9 — the 2nd-order ABC derivation this helper
///   implements.
/// * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
///   3rd ed., Wiley 2014, §10.4 (1st- and 2nd-order ABC face
///   contributions and reflection-floor tables).
/// * Phase 4.fem.eig.3 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
///   §4.2 — the bilinear form this helper implements.
pub fn assemble_abc2_face_block(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    k0: f64,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3> {
    // ---- 1st-order R_1 contribution. Re-use the 1st-order helper so
    // the `AbcOrder::First` path and the imaginary part of `R_1` in the
    // composite agree bit-for-bit. -----------------------------------
    let abc1 = assemble_abc_face_block(face_vertices, outward_normal, k0, mu_r_face);

    // ---- 2nd-order R_2 curl-correction contribution. --------------------
    //
    // Whitney-1 in-plane barycentric gradients ∇λ_a per face vertex —
    // constant across the face. (a, b, c) cyclic in {0, 1, 2}.
    let (grads, face_area) = face_barycentric_gradients_and_area(&face_vertices, outward_normal);

    // Normalise the outward normal once for the (n̂ · ·) reduction below.
    let n_norm = outward_normal.norm();
    let n_hat = if n_norm > 0.0 {
        outward_normal / n_norm
    } else {
        outward_normal
    };

    // Per-edge scalar surface-curl c_i = ∇_t × N_i = n̂ · (∇ × N_i)
    //                                  = 2 · n̂ · (∇λ_a × ∇λ_b)
    // with (a, b) = (i, (i + 1) mod 3). Constant per face because
    // ∇λ_a, ∇λ_b are constant; the 3D curl `∇ × N_i = 2 ∇λ_a × ∇λ_b`
    // is parallel to `n̂` on a planar face, so its scalar projection
    // onto `n̂` carries the full information.
    let mut c = [0.0_f64; 3];
    for (i, slot) in c.iter_mut().enumerate() {
        let a = i;
        let b = (i + 1) % 3;
        *slot = 2.0 * n_hat.dot(&grads[a].cross(&grads[b]));
    }

    // 2nd-order scalar prefactor: −(1 / (2 k₀)) · (A / μ_r,face). REAL.
    // (Engquist–Majda 1979 eq. 9 — the curl term has a real prefactor;
    // the composite block is complex-symmetric with Re from R_2 and Im
    // from R_1.)
    let curl_prefactor = -face_area / (2.0 * k0 * mu_r_face);

    let mut block = abc1;
    for i in 0..3 {
        for j in 0..3 {
            let r2_entry = c[i] * c[j];
            block[(i, j)] += Complex64::new(curl_prefactor * r2_entry, 0.0);
        }
    }

    block
}

/// Assemble the `6 × 6` Nedelec local stiffness + mass block for a
/// single tetrahedron with scalar real `ε_r`, `μ_r`.
///
/// `vertices` must be ordered so that the signed volume is positive
/// (`yee-mesh::TetMesh3D::new` enforces this; callers passing
/// hand-rolled vertices should ensure the same). The returned matrices
/// are emitted in canonical local-edge orientation per [`LOCAL_EDGES`]
/// — sign flips against global-edge orientation are the assembly
/// layer's job.
///
/// This is a thin wrapper over [`assemble_tet_element_complex`]:
/// the real `ε_r`, `μ_r` are lifted to `Complex64`, the complex local
/// matrices are computed, and the result is projected to real via
/// `.re`. For real inputs `Im(K_local) ≡ 0` and `Im(M_local) ≡ 0`, so
/// the projection is lossless and the returned matrices are bit-for-bit
/// identical to the Phase 4.fem.eig.0 shipped implementation.
///
/// # Panics
///
/// Panics if the tet is degenerate (all four vertices coplanar; the
/// barycentric gradient system is singular). Real meshes go through
/// [`yee_mesh::TetMesh3D::new`] which rejects such tets at
/// construction; this entry point trusts its caller.
pub fn assemble_tet_element(
    vertices: [Vector3<f64>; 4],
    eps_r: f64,
    mu_r: f64,
) -> NedelecTetElement {
    let elem_complex = assemble_tet_element_complex(
        vertices,
        Complex64::new(eps_r, 0.0),
        Complex64::new(mu_r, 0.0),
    );

    // Project to real. For real `eps_r`, `mu_r` the imaginary parts are
    // identically zero (within FP round-off of the `1/mu_r` reciprocal
    // step, but since `Im(mu_r) = 0` exactly the imaginary component
    // does not appear) so the `.re` projection is bit-for-bit lossless.
    let mut k_local = SMatrix::<f64, 6, 6>::zeros();
    let mut m_local = SMatrix::<f64, 6, 6>::zeros();
    for alpha in 0..6 {
        for beta in 0..6 {
            k_local[(alpha, beta)] = elem_complex.k_local[(alpha, beta)].re;
            m_local[(alpha, beta)] = elem_complex.m_local[(alpha, beta)].re;
        }
    }

    NedelecTetElement { k_local, m_local }
}
