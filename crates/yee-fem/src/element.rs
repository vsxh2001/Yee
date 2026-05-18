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
//! ## References
//!
//! * Jin, J.-M., *The Finite Element Method in Electromagnetics*,
//!   3rd ed., Wiley 2014, §9.4 (Nedelec edge elements on tetrahedra).
//! * Keast, P., "Moderate-degree tetrahedral quadrature formulas",
//!   *Comp. Methods Appl. Mech. Eng.* 55, 1986, pp. 339–348 (degree-2
//!   four-point symmetric rule).

use nalgebra::{Matrix3, SMatrix, Vector3};

/// The six canonical local edges of a tetrahedron in
/// lower-endpoint-first order.
///
/// Index `α ∈ 0..6` maps to local edge `LOCAL_EDGES[α] = (i, j)` with
/// `i < j`. This is the row / column ordering of [`NedelecTetElement`]
/// matrices.
pub const LOCAL_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// Per-tet element matrices for the first-order Nedelec edge basis.
///
/// Both fields are `6 × 6` real symmetric blocks indexed by
/// [`LOCAL_EDGES`]. See module docs for the exact bilinear forms.
#[derive(Debug, Clone, Copy)]
pub struct NedelecTetElement {
    /// Local curl-curl stiffness `K^e_{αβ} = (1/μ_r) V (∇×N_α)·(∇×N_β)`.
    pub k_local: SMatrix<f64, 6, 6>,
    /// Local vector mass `M^e_{αβ} = ε_r ∫_T N_α · N_β dV`.
    pub m_local: SMatrix<f64, 6, 6>,
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
/// single tetrahedron with scalar real `ε_r`, `μ_r`.
///
/// `vertices` must be ordered so that the signed volume is positive
/// (`yee-mesh::TetMesh3D::new` enforces this; callers passing
/// hand-rolled vertices should ensure the same). The returned matrices
/// are emitted in canonical local-edge orientation per [`LOCAL_EDGES`]
/// — sign flips against global-edge orientation are the assembly
/// layer's job.
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
    let (grads, signed_volume) = barycentric_gradients_and_volume(&vertices);
    let volume = signed_volume.abs();

    // ---- Local stiffness: K^e_{αβ} = (1/μ_r) V (∇×N_α)·(∇×N_β) ------
    //
    // ∇×N_{ij} = 2 ∇λ_i × ∇λ_j is constant on the tet, so the integral
    // is exactly volume × (curl_α · curl_β).
    let mut curls = [Vector3::<f64>::zeros(); 6];
    for (alpha, &(i, j)) in LOCAL_EDGES.iter().enumerate() {
        curls[alpha] = 2.0 * grads[i].cross(&grads[j]);
    }

    let mut k_local = SMatrix::<f64, 6, 6>::zeros();
    let inv_mu = 1.0 / mu_r;
    for alpha in 0..6 {
        for beta in 0..6 {
            k_local[(alpha, beta)] = inv_mu * volume * curls[alpha].dot(&curls[beta]);
        }
    }

    // ---- Local mass: M^e_{αβ} = ε_r ∫_T N_α · N_β dV --------------
    //
    // 4-point Gauss-tet quadrature, exact for polynomials of degree 2.
    // At each barycentric quadrature point λ = (λ_0, λ_1, λ_2, λ_3):
    //
    //     N_{ij}(λ) = λ_i ∇λ_j − λ_j ∇λ_i      (Vector3<f64>)
    //
    // and the integrand `N_α · N_β` is degree 2 in the λ_k. The weight
    // is V/4 per point.
    let weight = volume / 4.0;
    let mut m_local = SMatrix::<f64, 6, 6>::zeros();

    for qp in &QUAD_POINTS {
        // Evaluate all six basis vectors at this quadrature point.
        let mut basis = [Vector3::<f64>::zeros(); 6];
        for (alpha, &(i, j)) in LOCAL_EDGES.iter().enumerate() {
            basis[alpha] = qp[i] * grads[j] - qp[j] * grads[i];
        }
        for alpha in 0..6 {
            for beta in 0..6 {
                m_local[(alpha, beta)] += weight * basis[alpha].dot(&basis[beta]);
            }
        }
    }
    m_local *= eps_r;

    NedelecTetElement { k_local, m_local }
}
