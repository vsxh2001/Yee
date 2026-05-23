//! Nedelec edge-element + nodal-Lagrange element-matrix assembly for the
//! 2-D wave-port cross-section eigenproblem.
//!
//! The discretisation of the transverse-vector Helmholtz equation on a
//! guided-wave cross-section (Jin, *FEM in Electromagnetics* 3rd ed.
//! §8.5) gives a generalized eigenproblem
//!
//! ```text
//!     A x = β² B x
//! ```
//!
//! where for a TE-mode-capable formulation on a homogeneously-filled
//! PEC waveguide the **transverse-only** form is exact (the longitudinal
//! `E_z` is identically zero, enforced by the PEC walls and the
//! homogeneous interior):
//!
//! ```text
//!     A[i,j] = k_0² ∫_Ω ε_r N_i · N_j dA  −  ∫_Ω (1/μ_r) (∇×N_i)(∇×N_j) dA
//!     B[i,j] = ∫_Ω ε_r N_i · N_j dA
//! ```
//!
//! Here `N_i` are the first-order curl-conforming Nedelec (Whitney-1)
//! edge basis functions:
//!
//! ```text
//!     N_e = ℓ_e σ_e (λ_a ∇λ_b − λ_b ∇λ_a)
//! ```
//!
//! with `(a, b)` the edge endpoints (local CCW traversal — see
//! [`super::mesh`]), `σ_e ∈ {+1, −1}` the per-triangle orientation sign
//! relative to the canonical global-edge direction, and `ℓ_e` the edge
//! length. The curl is constant on each triangle:
//!
//! ```text
//!     ∇ × N_e = ℓ_e σ_e / A · ẑ
//! ```
//!
//! Element matrices use closed-form integrals over linear triangles
//! (`∫ λ_p λ_q dA = (A/12)(1 + δ_{pq})`, `∇λ_p · ∇λ_q = (b_p b_q + c_p c_q)/(4A²)`).
//!
//! **Mixed formulation (Phase 1.3.1.1 step 5).** [`assemble_mixed`]
//! builds the full mixed `(E_t, E_z)` Lee-Sun-Cendes block pencil from
//! the longitudinal nodal-Lagrange element matrices [`local_a_zz`] /
//! [`local_b_zz`] and the `1/μ_r`-weighted edge-node coupling
//! [`local_b_ze`], consumed by [`super::solve_dense_mixed`] for
//! inhomogeneous (dielectric-loaded / microstrip) cross-sections. The
//! transverse-only [`assemble_transverse`] / [`super::solve_dense`] are
//! retained as the homogeneous reference and its regression tests; on a
//! homogeneous guide the two agree to machine precision (the dominant
//! mode is purely transverse). See [`AssembledMixed`] for the block
//! layout and the load-bearing-coupling discussion.

use nalgebra::{Complex, DMatrix};
use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::{MaterialTag, TriMesh2D};

use super::mesh::EdgeTable;

/// Per-triangle FEM scalars used by every element matrix.
///
/// `b[i] = y_{i+1} - y_{i+2}`, `c[i] = x_{i+2} - x_{i+1}` (indices mod 3,
/// CCW vertex order) so that the gradient of barycentric coordinate `λ_i`
/// is `∇λ_i = (b[i], c[i]) / (2A)`.
#[derive(Debug, Clone, Copy)]
struct TriGeom {
    area: f64,
    /// `b[i] = y_{(i+1)%3} - y_{(i+2)%3}`
    b: [f64; 3],
    /// `c[i] = x_{(i+2)%3} - x_{(i+1)%3}`
    c: [f64; 3],
    /// Edge lengths in local order (edge `e` opposite local vertex `e`).
    edge_len: [f64; 3],
}

impl TriGeom {
    fn from_mesh(mesh: &TriMesh2D, tri_idx: usize) -> Self {
        let tri = mesh.triangles[tri_idx];
        let v: [[f64; 2]; 3] = [
            mesh.vertices[tri[0]],
            mesh.vertices[tri[1]],
            mesh.vertices[tri[2]],
        ];
        let area = mesh.area(tri_idx);
        let mut b = [0.0; 3];
        let mut c = [0.0; 3];
        let mut edge_len = [0.0; 3];
        for i in 0..3 {
            let i1 = (i + 1) % 3;
            let i2 = (i + 2) % 3;
            b[i] = v[i1][1] - v[i2][1];
            c[i] = v[i2][0] - v[i1][0];
            // Edge `i` connects local vertex `i1` to `i2` (per Jin §8.5
            // / `mesh.rs` convention: edge i is opposite vertex i).
            let dx = v[i2][0] - v[i1][0];
            let dy = v[i2][1] - v[i1][1];
            edge_len[i] = (dx * dx + dy * dy).sqrt();
        }
        Self {
            area,
            b,
            c,
            edge_len,
        }
    }

    /// `∇λ_p · ∇λ_q` (constant on the triangle).
    fn grad_dot(&self, p: usize, q: usize) -> f64 {
        (self.b[p] * self.b[q] + self.c[p] * self.c[q]) / (4.0 * self.area * self.area)
    }
}

/// `∫_T λ_p λ_q dA = (A/12)(1 + δ_{pq})` — the standard linear-triangle
/// mass integral.
fn int_lambda_lambda(area: f64, p: usize, q: usize) -> f64 {
    let delta = if p == q { 1.0 } else { 0.0 };
    area * (1.0 + delta) / 12.0
}

/// Local edge endpoints, matching [`super::mesh::EdgeTable::build`].
/// Edge `e` opposite local vertex `e`, traversed CCW.
const LOCAL_EDGES: [[usize; 2]; 3] = [[1, 2], [2, 0], [0, 1]];

/// Local Nedelec curl-curl stiffness on a single triangle:
///
/// ```text
///   S^e[i,j] = ∫_T (1/μ_r) (∇×N_i)(∇×N_j) dA
///            = σ_i σ_j ℓ_i ℓ_j / (μ_r · A)
/// ```
///
/// `signs` carry the per-triangle orientation of each local edge (`+1`
/// or `−1`) against the canonical global-edge direction.
// 3×3 dense local matrix assembly is naturally written with explicit
// (i, j) indices; `enumerate`-then-index zip would read worse.
#[allow(clippy::needless_range_loop)]
fn local_a_ee_curl(geom: &TriGeom, mu_r: Complex64, signs: [f64; 3]) -> [[Complex64; 3]; 3] {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    let inv_mu_a = Complex::new(1.0 / geom.area, 0.0) / mu_r;
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = inv_mu_a
                * Complex::new(
                    signs[i] * signs[j] * geom.edge_len[i] * geom.edge_len[j],
                    0.0,
                );
        }
    }
    out
}

/// Local Nedelec mass matrix on a single triangle:
///
/// ```text
///   T^e[i,j] = ∫_T ε_r N_i · N_j dA
/// ```
///
/// where `N_e = ℓ_e σ_e (λ_a ∇λ_b − λ_b ∇λ_a)` with `(a, b)` the local
/// edge endpoints. Expanding the dot product and integrating
/// term-by-term using `∫ λ_p λ_q dA = (A/12)(1+δ_{pq})` yields a
/// closed-form 3×3 dense block.
#[allow(clippy::needless_range_loop)]
fn local_b_ee_mass(geom: &TriGeom, eps_r: Complex64, signs: [f64; 3]) -> [[Complex64; 3]; 3] {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    let area = geom.area;
    for i in 0..3 {
        let [ai, bi] = LOCAL_EDGES[i];
        for j in 0..3 {
            let [aj, bj] = LOCAL_EDGES[j];
            // (λ_{ai} ∇λ_{bi} − λ_{bi} ∇λ_{ai}) · (λ_{aj} ∇λ_{bj} − λ_{bj} ∇λ_{aj})
            // =   λ_{ai} λ_{aj} (∇λ_{bi}·∇λ_{bj})
            //   − λ_{ai} λ_{bj} (∇λ_{bi}·∇λ_{aj})
            //   − λ_{bi} λ_{aj} (∇λ_{ai}·∇λ_{bj})
            //   + λ_{bi} λ_{bj} (∇λ_{ai}·∇λ_{aj})
            let s = int_lambda_lambda(area, ai, aj) * geom.grad_dot(bi, bj)
                - int_lambda_lambda(area, ai, bj) * geom.grad_dot(bi, aj)
                - int_lambda_lambda(area, bi, aj) * geom.grad_dot(ai, bj)
                + int_lambda_lambda(area, bi, bj) * geom.grad_dot(ai, aj);
            let coeff = geom.edge_len[i] * geom.edge_len[j] * signs[i] * signs[j];
            out[i][j] = eps_r * Complex::new(coeff * s, 0.0);
        }
    }
    out
}

/// Local nodal-Lagrange gradient-gradient stiffness on a single triangle:
///
/// ```text
///   A_zz^e[i,j] = ∫_T (1/μ_r) ∇L_i · ∇L_j dA = (1/μ_r) · A · ∇λ_i · ∇λ_j
/// ```
///
/// Used by the mixed Lee-Sun-Cendes formulation
/// ([`assemble_mixed`]). Unused by the transverse-only WR-90
/// validation path.
#[allow(clippy::needless_range_loop)]
fn local_a_zz(geom: &TriGeom, mu_r: Complex64) -> [[Complex64; 3]; 3] {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    let inv_mu = Complex::new(1.0, 0.0) / mu_r;
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = inv_mu * Complex::new(geom.area * geom.grad_dot(i, j), 0.0);
        }
    }
    out
}

/// Local nodal-Lagrange mass on a single triangle:
///
/// ```text
///   B_zz^e[i,j] = ∫_T ε_r L_i L_j dA = ε_r (A/12)(1 + δ_{ij})
/// ```
///
/// Used by the mixed Lee-Sun-Cendes formulation
/// ([`assemble_mixed`]). Unused by the transverse-only WR-90
/// validation path.
#[allow(clippy::needless_range_loop)]
fn local_b_zz(geom: &TriGeom, eps_r: Complex64) -> [[Complex64; 3]; 3] {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = eps_r * Complex::new(int_lambda_lambda(geom.area, i, j), 0.0);
        }
    }
    out
}

/// Local edge-node coupling on a single triangle:
///
/// ```text
///   B_ze^e[i,j] = ∫_T (1/μ_r) ∇L_i · N_j dA
/// ```
///
/// `N_j = ℓ_j σ_j (λ_{a_j} ∇λ_{b_j} − λ_{b_j} ∇λ_{a_j})` and
/// `∇L_i = ∇λ_i` is constant, so
/// `∫ ∇λ_i · N_j dA = ℓ_j σ_j [(∫ λ_{a_j} dA) (∇λ_i · ∇λ_{b_j})
///                              − (∫ λ_{b_j} dA) (∇λ_i · ∇λ_{a_j})]`
/// with `∫ λ_p dA = A/3`.
///
/// **Weight convention (Phase 1.3.1.1 step-5 review fix — `1/μ_r`, not
/// `ε_r`).** The Lee-Sun-Cendes inter-block coupling arises from the
/// **curl-curl** term of the transverse vector Helmholtz functional
/// (`∫(1/μ_r)(∇_t E_z + jβ E_t)·(…)` — the `jβ ∇_t E_z · N` cross term),
/// so it carries the `1/μ_r` weight, matching [`local_a_zz`]. The
/// originally-staged `ε_r` weight made the coupling a *divergence-penalty*
/// term (`∫ε_r ∇L·E_t`), which the curl-curl eigenvector is **exactly
/// orthogonal** to in the ε_r-mass inner product (Boffi-Brezzi-Demkowicz:
/// the Whitney-1 curl kernel is the Whitney-0 gradient space, and the
/// eigenvector is `T = ∫ε_r N·N`-orthogonal to gradients). With the `ε_r`
/// weight the coupling was therefore **structurally inert** for the
/// dominant mode (`B_zt e_t = 0` to machine precision, `E_z ≡ 0` for any
/// piecewise-constant fill) — so it never participated and could not be
/// validated, which is exactly the step-5-review coverage gap. The
/// `1/μ_r` weight is *not* annihilated (`∫(1/μ_r)∇L·E_t ≠ 0` when ε_r
/// varies, since `E_t` is ε_r-orthogonal, not `1/μ_r`-orthogonal, to
/// gradients), so the dominant mode of an inhomogeneous guide genuinely
/// develops `E_z ≠ 0`. On a homogeneous guide with `μ_r = ε_r = 1` both
/// weights coincide and give `E_z ≡ 0`, preserving the WR-90 canary.
///
/// Used by the mixed Lee-Sun-Cendes formulation ([`assemble_mixed`]).
/// Unused by the transverse-only WR-90 validation path.
#[allow(clippy::needless_range_loop)]
fn local_b_ze(geom: &TriGeom, mu_r: Complex64, signs: [f64; 3]) -> [[Complex64; 3]; 3] {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    let inv_mu = Complex::new(1.0, 0.0) / mu_r;
    let third_area = geom.area / 3.0;
    for i in 0..3 {
        for j in 0..3 {
            let [aj, bj] = LOCAL_EDGES[j];
            let s = third_area * (geom.grad_dot(i, bj) - geom.grad_dot(i, aj));
            out[i][j] = inv_mu * Complex::new(geom.edge_len[j] * signs[j] * s, 0.0);
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────
// Phase 1.3.1.1 step 5.5 — second-order (p=2) elements.
//
// At second order the Nedelec curl is no longer constant on a triangle, so
// the curl-curl / mass / coupling element integrals cannot use the
// first-order closed forms. Every p=2 integral instead goes through a 2-D
// triangle Gauss rule exact for the integrand degree (the p=2 mass / nodal
// mass integrand is quartic in the barycentric coordinates, so a degree-4
// rule is required). The basis functions are evaluated point-wise from the
// barycentric coordinates and their (constant) gradients, exactly as in the
// `local_b_ze` independent-quadrature unit test, and the same point-wise
// definitions feed both the production element matrices and the J1
// independent-quadrature pin.
// ─────────────────────────────────────────────────────────────────────────

/// A 6-point degree-4 symmetric triangle Gauss rule (Dunavant 1985),
/// returned as `(λ₀, λ₁, λ₂, weight)` tuples with weights summing to **1**
/// (the reference-triangle convention — multiply by the physical area to
/// integrate). Exact for bivariate polynomials up to total degree 4, which
/// covers every p=2 element integrand: the quartic Nedelec / nodal mass
/// (`∫ε_r N·N`, `∫ε_r L·L`), the quadratic curl-curl / gradient stiffness,
/// and the cubic edge-node coupling `∫(1/μ_r)∇L·N`.
///
/// Hand-rolled (no `Cargo.toml` dependency, per the step-5.5 lane). The
/// weights/points are pinned exact-on-monomials in
/// [`tests::tri_gauss_deg4_integrates_quartics_exactly`].
pub(crate) fn tri_gauss_deg4() -> [(f64, f64, f64, f64); 6] {
    // Orbit 1: (a, a, 1−2a), weight w_a, three cyclic permutations.
    let a = 0.445_948_490_915_965;
    let wa = 0.223_381_589_678_011;
    // Orbit 2: (b, b, 1−2b), weight w_b, three cyclic permutations.
    let b = 0.091_576_213_509_771;
    let wb = 0.109_951_743_655_322;
    [
        (a, a, 1.0 - 2.0 * a, wa),
        (a, 1.0 - 2.0 * a, a, wa),
        (1.0 - 2.0 * a, a, a, wa),
        (b, b, 1.0 - 2.0 * b, wb),
        (b, 1.0 - 2.0 * b, b, wb),
        (1.0 - 2.0 * b, b, b, wb),
    ]
}

/// The eight second-order Nedelec-first-kind basis vectors of a triangle,
/// evaluated at a point given by its barycentric coordinates `lam` with
/// the (constant, per-triangle) barycentric gradients `g[i] = ∇λ_i`.
///
/// **Layout (the row/column order of every p=2 `E_t` element block):**
/// * `0..6` — two functions per local edge `e` (`LOCAL_EDGES[e] = (a, b)`):
///   - `2e`   : the Whitney / rotational function `W_e = λ_a ∇λ_b − λ_b ∇λ_a`
///     (the first-order edge function — its tangential trace on edge `e`
///     is constant);
///   - `2e+1` : the gradient function `G_e = λ_a ∇λ_b + λ_b ∇λ_a
///     = ∇(λ_a λ_b)` (the added second-order edge DoF).
/// * `6..8` — the two interior ("face") functions, each with **zero
///   tangential trace on all three edges** (so they never couple between
///   elements and carry no orientation sign):
///   - `6` : `F₀ = λ₀ (λ₁ ∇λ₂ − λ₂ ∇λ₁)`
///   - `7` : `F₁ = λ₁ (λ₂ ∇λ₀ − λ₀ ∇λ₂)`
///
/// Together these span the Nedelec-first-kind order-2 space (`dim = 8` on a
/// triangle: 2/edge + 2 interior). Completeness (rank 8) and the
/// tangential-trace property are pinned in
/// [`tests::p2_basis_has_full_rank_eight`] and
/// [`tests::p2_interior_functions_have_zero_tangential_trace`].
///
/// The **edge-function orientation sign** (`±1` per edge, from the global
/// edge-direction reconciliation) multiplies the *whole* edge function
/// (both `W_e` and `G_e`) the same way the first-order Whitney sign does —
/// it is applied by the caller during scatter, not here, so this evaluator
/// is a pure function of geometry. The interior functions take no sign.
#[inline]
fn p2_nedelec_basis(lam: [f64; 3], g: [[f64; 2]; 3]) -> [[f64; 2]; 8] {
    let mut out = [[0.0; 2]; 8];
    // Edge functions: Whitney + gradient, two per edge.
    for (e, &[a, b]) in LOCAL_EDGES.iter().enumerate() {
        let w = [
            lam[a] * g[b][0] - lam[b] * g[a][0],
            lam[a] * g[b][1] - lam[b] * g[a][1],
        ];
        let grad = [
            lam[a] * g[b][0] + lam[b] * g[a][0],
            lam[a] * g[b][1] + lam[b] * g[a][1],
        ];
        out[2 * e] = w;
        out[2 * e + 1] = grad;
    }
    // Interior functions F_k = λ_k (λ_{k+1} ∇λ_{k+2} − λ_{k+2} ∇λ_{k+1}),
    // k = 0, 1 (two independent of the three cyclic candidates).
    for (slot, k) in [6usize, 7].into_iter().zip([0usize, 1]) {
        let k1 = (k + 1) % 3;
        let k2 = (k + 2) % 3;
        out[slot] = [
            lam[k] * (lam[k1] * g[k2][0] - lam[k2] * g[k1][0]),
            lam[k] * (lam[k1] * g[k2][1] - lam[k2] * g[k1][1]),
        ];
    }
    out
}

/// `∇ × N` (the scalar `ẑ`-component) for each of the eight p=2 Nedelec
/// basis vectors at barycentric point `lam`, with constant gradients
/// `g[i] = ∇λ_i`.
///
/// For a planar vector field `N = (N_x, N_y)` the curl is the scalar
/// `∂N_y/∂x − ∂N_x/∂y`. With `λ` linear (`∇λ` constant) the curls follow
/// from `∇×(φ ∇ψ) = ∇φ × ∇ψ` (a scalar in 2-D, `g_φ × g_ψ ≡ g_φ.x g_ψ.y −
/// g_φ.y g_ψ.x`):
/// * `∇×W_e = 2 (∇λ_a × ∇λ_b)` — **constant** (the Whitney curl);
/// * `∇×G_e = ∇×∇(λ_aλ_b) = 0` — a pure gradient is curl-free;
/// * `∇×F_k = ∇×(λ_k(λ_{k1}∇λ_{k2} − λ_{k2}∇λ_{k1}))`
///   `= ∇λ_k×(λ_{k1}∇λ_{k2}−λ_{k2}∇λ_{k1}) + λ_k(∇λ_{k1}×∇λ_{k2} −
///   ∇λ_{k2}×∇λ_{k1})`, which is **linear** in `λ` (non-constant — the
///   reason p=2 needs quadrature).
///
/// `cross(u, v) = u.x v.y − u.y v.x`.
#[inline]
fn p2_nedelec_curl(lam: [f64; 3], g: [[f64; 2]; 3]) -> [f64; 8] {
    let cross = |u: [f64; 2], v: [f64; 2]| u[0] * v[1] - u[1] * v[0];
    let mut out = [0.0; 8];
    for (e, &[a, b]) in LOCAL_EDGES.iter().enumerate() {
        // ∇×W_e = 2 ∇λ_a × ∇λ_b (constant); ∇×G_e = 0.
        out[2 * e] = 2.0 * cross(g[a], g[b]);
        out[2 * e + 1] = 0.0;
    }
    for (slot, k) in [6usize, 7].into_iter().zip([0usize, 1]) {
        let k1 = (k + 1) % 3;
        let k2 = (k + 2) % 3;
        // ∇×(λ_k v) where v = λ_{k1}∇λ_{k2} − λ_{k2}∇λ_{k1}:
        //   = ∇λ_k × v  +  λ_k (∇×v),
        //   ∇×v = ∇λ_{k1}×∇λ_{k2} − ∇λ_{k2}×∇λ_{k1} = 2 ∇λ_{k1}×∇λ_{k2}.
        let v = [
            lam[k1] * g[k2][0] - lam[k2] * g[k1][0],
            lam[k1] * g[k2][1] - lam[k2] * g[k1][1],
        ];
        out[slot] = cross(g[k], v) + lam[k] * 2.0 * cross(g[k1], g[k2]);
    }
    out
}

/// The six quadratic nodal-Lagrange shape functions of a triangle at
/// barycentric point `lam`.
///
/// **Layout (the row/column order of every p=2 `E_z` block):**
/// * `0..3` — vertex nodes: `L_i = λ_i (2 λ_i − 1)`;
/// * `3..6` — edge-midpoint nodes, one per local edge `e`
///   (`LOCAL_EDGES[e] = (a, b)`): `L_{3+e} = 4 λ_a λ_b`.
///
/// (Standard P2 nodal basis: `L_node(node_j) = δ` at the 6 nodes.)
#[inline]
fn p2_nodal_basis(lam: [f64; 3]) -> [f64; 6] {
    let mut out = [0.0; 6];
    for i in 0..3 {
        out[i] = lam[i] * (2.0 * lam[i] - 1.0);
    }
    for (e, &[a, b]) in LOCAL_EDGES.iter().enumerate() {
        out[3 + e] = 4.0 * lam[a] * lam[b];
    }
    out
}

/// Gradients `∇L` of the six quadratic nodal-Lagrange shape functions at
/// barycentric point `lam`, with constant `g[i] = ∇λ_i`. Same layout as
/// [`p2_nodal_basis`]:
/// * vertex: `∇L_i = (4 λ_i − 1) ∇λ_i`;
/// * midpoint of edge `(a, b)`: `∇L = 4 (λ_a ∇λ_b + λ_b ∇λ_a)`.
#[inline]
fn p2_nodal_grad(lam: [f64; 3], g: [[f64; 2]; 3]) -> [[f64; 2]; 6] {
    let mut out = [[0.0; 2]; 6];
    for i in 0..3 {
        let s = 4.0 * lam[i] - 1.0;
        out[i] = [s * g[i][0], s * g[i][1]];
    }
    for (e, &[a, b]) in LOCAL_EDGES.iter().enumerate() {
        out[3 + e] = [
            4.0 * (lam[a] * g[b][0] + lam[b] * g[a][0]),
            4.0 * (lam[a] * g[b][1] + lam[b] * g[a][1]),
        ];
    }
    out
}

/// Per-triangle barycentric gradients `∇λ_i = (b_i, c_i)/(2A)` as `[f64; 2]`
/// vectors, the constant data every p=2 point-wise evaluator needs.
#[inline]
fn bary_grads(geom: &TriGeom) -> [[f64; 2]; 3] {
    let inv = 1.0 / (2.0 * geom.area);
    [
        [geom.b[0] * inv, geom.c[0] * inv],
        [geom.b[1] * inv, geom.c[1] * inv],
        [geom.b[2] * inv, geom.c[2] * inv],
    ]
}

/// Local p=2 Nedelec curl-curl stiffness `A_tt^e[i,j] =
/// ∫_T (1/μ_r)(∇×N_i)(∇×N_j) dA`, the 8×8 second-order analogue of
/// [`local_a_ee_curl`]. The curl is linear (non-constant), so the integral
/// is taken by the degree-4 [`tri_gauss_deg4`] rule (exact for the
/// quadratic integrand). `signs[e] ∈ {±1}` is the per-edge orientation,
/// applied to **both** edge DoFs of edge `e` (slots `2e`, `2e+1`); interior
/// slots `6, 7` take sign `+1`.
fn local_a_ee_curl_p2(geom: &TriGeom, mu_r: Complex64, signs: [f64; 3]) -> [[Complex64; 8]; 8] {
    let g = bary_grads(geom);
    let sgn = p2_edge_signs(signs);
    let inv_mu = Complex::new(1.0, 0.0) / mu_r;
    let mut acc = [[0.0f64; 8]; 8];
    for (l0, l1, l2, w) in tri_gauss_deg4() {
        let curl = p2_nedelec_curl([l0, l1, l2], g);
        let wa = w * geom.area;
        for i in 0..8 {
            for j in 0..8 {
                acc[i][j] += wa * sgn[i] * sgn[j] * curl[i] * curl[j];
            }
        }
    }
    finish_complex8(acc, inv_mu)
}

/// Local p=2 Nedelec mass `B_tt^e[i,j] = ∫_T ε_r N_i·N_j dA`, the 8×8
/// second-order analogue of [`local_b_ee_mass`]. The integrand is quartic,
/// so the degree-4 rule integrates it exactly. Orientation handling matches
/// [`local_a_ee_curl_p2`].
fn local_b_ee_mass_p2(geom: &TriGeom, eps_r: Complex64, signs: [f64; 3]) -> [[Complex64; 8]; 8] {
    let g = bary_grads(geom);
    let sgn = p2_edge_signs(signs);
    let mut acc = [[0.0f64; 8]; 8];
    for (l0, l1, l2, w) in tri_gauss_deg4() {
        let n = p2_nedelec_basis([l0, l1, l2], g);
        let wa = w * geom.area;
        for i in 0..8 {
            for j in 0..8 {
                let dot = n[i][0] * n[j][0] + n[i][1] * n[j][1];
                acc[i][j] += wa * sgn[i] * sgn[j] * dot;
            }
        }
    }
    finish_complex8(acc, eps_r)
}

/// Local p=2 nodal gradient-gradient stiffness `A_zz^e[i,j] =
/// ∫_T (1/μ_r) ∇L_i·∇L_j dA`, the 6×6 second-order analogue of
/// [`local_a_zz`]. `∇L` is linear; the degree-4 rule is exact.
fn local_a_zz_p2(geom: &TriGeom, mu_r: Complex64) -> [[Complex64; 6]; 6] {
    let g = bary_grads(geom);
    let inv_mu = Complex::new(1.0, 0.0) / mu_r;
    let mut acc = [[0.0f64; 6]; 6];
    for (l0, l1, l2, w) in tri_gauss_deg4() {
        let gl = p2_nodal_grad([l0, l1, l2], g);
        let wa = w * geom.area;
        for i in 0..6 {
            for j in 0..6 {
                acc[i][j] += wa * (gl[i][0] * gl[j][0] + gl[i][1] * gl[j][1]);
            }
        }
    }
    finish_complex6(acc, inv_mu)
}

/// Local p=2 nodal mass `B_zz^e[i,j] = ∫_T ε_r L_i L_j dA`, the 6×6
/// second-order analogue of [`local_b_zz`]. The integrand is quartic; the
/// degree-4 rule is exact.
fn local_b_zz_p2(geom: &TriGeom, eps_r: Complex64) -> [[Complex64; 6]; 6] {
    let mut acc = [[0.0f64; 6]; 6];
    for (l0, l1, l2, w) in tri_gauss_deg4() {
        let l = p2_nodal_basis([l0, l1, l2]);
        let wa = w * geom.area;
        for i in 0..6 {
            for j in 0..6 {
                acc[i][j] += wa * l[i] * l[j];
            }
        }
    }
    finish_complex6(acc, eps_r)
}

/// Local p=2 edge-node coupling `B_ze^e[i_node][j_edge] =
/// ∫_T (1/μ_r) ∇L_i·N_j dA` (nodal-row / Nedelec-col), the 6×8
/// second-order analogue of [`local_b_ze`]. Integrand is cubic; the
/// degree-4 rule is exact. The Nedelec orientation sign multiplies the
/// edge-function columns (interior columns take `+1`); the nodal rows take
/// no sign. The `1/μ_r` weight is the curl-curl cross term, matching
/// [`local_a_zz_p2`] and the first-order [`local_b_ze`].
fn local_b_ze_p2(geom: &TriGeom, mu_r: Complex64, signs: [f64; 3]) -> [[Complex64; 8]; 6] {
    let g = bary_grads(geom);
    let sgn = p2_edge_signs(signs);
    let inv_mu = Complex::new(1.0, 0.0) / mu_r;
    let mut acc = [[0.0f64; 8]; 6];
    for (l0, l1, l2, w) in tri_gauss_deg4() {
        let gl = p2_nodal_grad([l0, l1, l2], g);
        let n = p2_nedelec_basis([l0, l1, l2], g);
        let wa = w * geom.area;
        for i in 0..6 {
            for j in 0..8 {
                let dot = gl[i][0] * n[j][0] + gl[i][1] * n[j][1];
                acc[i][j] += wa * sgn[j] * dot;
            }
        }
    }
    let mut out = [[Complex64::new(0.0, 0.0); 8]; 6];
    for i in 0..6 {
        for j in 0..8 {
            out[i][j] = inv_mu * Complex::new(acc[i][j], 0.0);
        }
    }
    out
}

/// Expand a 3-edge orientation-sign triple into the 8-slot p=2 `E_t` sign
/// vector — the local→global basis transformation for each DoF slot.
///
/// **The two edge DoFs transform differently under edge reversal** (the
/// step-5.5 risk (b) subtlety):
/// * the Whitney slot `2e` (`W_e = λ_a∇λ_b − λ_b∇λ_a`) is **odd** under
///   swapping the edge endpoints `a ↔ b`, so it carries the orientation
///   sign `signs[e] ∈ {±1}` — exactly as the first-order Whitney DoF does;
/// * the gradient slot `2e+1` (`G_e = ∇(λ_a λ_b)`) is **even** (`λ_a λ_b`
///   is symmetric in `a, b`), so its tangential trace on the shared edge is
///   single-valued *without* any sign flip → sign **`+1` always**,
///   independent of orientation.
///
/// Mixing these up (signing the gradient DoF) would make the global
/// gradient-edge DoF double-valued and corrupt the assembly; it is caught
/// by the J3 homogeneous-TE10 anchor (a wrong sign there fails to reproduce
/// the analytic β). The two interior slots (`6, 7`) take `+1` (interior
/// functions have zero tangential trace, so no global orientation).
#[inline]
fn p2_edge_signs(signs: [f64; 3]) -> [f64; 8] {
    [signs[0], 1.0, signs[1], 1.0, signs[2], 1.0, 1.0, 1.0]
}

/// Multiply a real 8×8 accumulator by a complex material weight, producing
/// the complex element block (lossless path keeps the imaginary part zero).
#[inline]
#[allow(clippy::needless_range_loop)]
fn finish_complex8(acc: [[f64; 8]; 8], weight: Complex64) -> [[Complex64; 8]; 8] {
    let mut out = [[Complex64::new(0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            out[i][j] = weight * Complex::new(acc[i][j], 0.0);
        }
    }
    out
}

/// Multiply a real 6×6 accumulator by a complex material weight.
#[inline]
#[allow(clippy::needless_range_loop)]
fn finish_complex6(acc: [[f64; 6]; 6], weight: Complex64) -> [[Complex64; 6]; 6] {
    let mut out = [[Complex64::new(0.0, 0.0); 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = weight * Complex::new(acc[i][j], 0.0);
        }
    }
    out
}

/// Result of [`assemble_transverse`]: the generalized eigenproblem
/// `S e_t = k_c² T e_t` reduced to the interior-edge DoFs by elimination
/// of the PEC-boundary tangential `E_t` DoFs.
///
/// The propagation constant follows from the cutoff eigenvalue via
/// `β² = k_0² − k_c²`. The eigenproblem itself is **frequency-independent**:
/// `S` and `T` depend only on geometry and material; `k_0²` enters only
/// through the post-solve mapping `(k_c², ε_r) → β`. This is the
/// formulation that lets the spurious gradient null-space cluster at
/// `k_c² = 0` (easy to filter out) instead of `k_c² = k_0²` (impossible
/// to distinguish from the physical mode in a single-precision solve).
///
/// Retained as the transverse-only reference (and exercised by the
/// homogeneous-guide regression tests) after
/// [`crate::ports::NumericalCrossSection::solve`] switched to the mixed
/// [`AssembledMixed`] / [`super::solve_dense_mixed`] path in Phase
/// 1.3.1.1 step 5.
#[allow(dead_code)]
pub(crate) struct AssembledTransverse {
    /// Curl-curl stiffness matrix `S[i,j] = ∫ (1/μ_r) (∇×N_i)(∇×N_j) dA`.
    pub s: DMatrix<Complex64>,
    /// ε_r-weighted Nedelec mass matrix `T_ε[i,j] = ∫ ε_r N_i · N_j dA`.
    pub t: DMatrix<Complex64>,
    /// **Unweighted** Nedelec mass matrix `T_1[i,j] = ∫ N_i · N_j dA`
    /// (ε_r ≡ 1). This is the RHS metric of the β-direct generalized
    /// eigenproblem `(k_0² T_ε − S) x = β² T_1 x` (Phase 1.3.1.1 step 5.2):
    /// the physical transverse Helmholtz equation
    /// `∇×(1/μ_r ∇×E_t) = (k_0² ε_r − β²) E_t` puts ε_r only on the
    /// `k_0²` side, so the `−β² E_t` term carries the **unweighted** mass.
    /// The earlier `S x = k_c² T_ε x` / `β² = k_0² − k_c²` extraction was
    /// algebraically correct only for `ε_r ≡ 1` (then `T_1 = T_ε`); for any
    /// `ε_r ≠ 1` it under-counted the dielectric.
    pub t1: DMatrix<Complex64>,
    /// Map from interior-edge DoF index to the global edge index, so
    /// post-solve eigenvector components can be located back on the mesh.
    #[allow(dead_code)] // consumed by the Phase 1.3.1.1 step 5 eigenvector recovery
    pub interior_to_global: Vec<usize>,
    /// Largest ε_r magnitude seen during assembly. Retained as a diagnostic
    /// of the fill contrast; the β-extraction itself no longer needs it
    /// (the β-direct form `(k_0² T_ε − S) x = β² T_1 x` carries ε_r through
    /// the `T_ε` operator, not a scalar correction). Lossy / heterogeneous
    /// complex ε_r is Phase 1.3.1.2.
    #[allow(dead_code)]
    pub eps_r_max_re: f64,
}

/// Assemble the transverse-only (`E_t`-block) generalized eigenproblem
/// `S e_t = k_c² T e_t` on the supplied mesh.
///
/// PEC boundary edges (those incident on exactly one triangle) are
/// eliminated by Dirichlet condition `E_t = 0`. The returned matrices
/// are dense `n × n` with `n` = interior edge count.
///
/// Material data (`eps_r`, `mu_r`) is keyed by [`MaterialTag`] and
/// looked up per triangle via [`TriMesh2D::triangle_material`]. Missing
/// keys default to vacuum (ε_r = μ_r = 1). The assembly is
/// **frequency-independent**; the caller maps the resulting cutoff
/// eigenvalue to a propagation constant via `β² = k_0² − k_c²`.
///
/// Retained as the transverse-only reference path after
/// [`crate::ports::NumericalCrossSection::solve`] switched to
/// [`assemble_mixed`] in Phase 1.3.1.1 step 5; exercised by the
/// homogeneous-guide regression tests.
#[allow(dead_code)]
pub(crate) fn assemble_transverse(
    mesh: &TriMesh2D,
    eps_r: &HashMap<MaterialTag, Complex64>,
    mu_r: &HashMap<MaterialTag, Complex64>,
    edge_table: &EdgeTable,
) -> AssembledTransverse {
    // Interior-edge DoF numbering: map global edge -> interior DoF index,
    // and remember the reverse for eigenvector recovery.
    let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; edge_table.n_edges()];
    let mut interior_to_global: Vec<usize> = Vec::new();
    for (gid, &is_bnd) in edge_table.is_boundary.iter().enumerate() {
        if !is_bnd {
            interior_dof_of_edge[gid] = Some(interior_to_global.len());
            interior_to_global.push(gid);
        }
    }
    let n = interior_to_global.len();

    let zero = Complex64::new(0.0, 0.0);
    let mut s = DMatrix::from_element(n, n, zero);
    let mut t = DMatrix::from_element(n, n, zero);
    let mut t1 = DMatrix::from_element(n, n, zero);
    let mut eps_r_max_re: f64 = 1.0;

    let default_one = Complex64::new(1.0, 0.0);

    for (tri_idx, conn) in edge_table.tri_edges.iter().enumerate() {
        let geom = TriGeom::from_mesh(mesh, tri_idx);
        let tag = mesh.triangle_material[tri_idx];
        let eps = *eps_r.get(&tag).unwrap_or(&default_one);
        let mu = *mu_r.get(&tag).unwrap_or(&default_one);
        if eps.re > eps_r_max_re {
            eps_r_max_re = eps.re;
        }

        let s_local = local_a_ee_curl(&geom, mu, conn.sign);
        let t_local = local_b_ee_mass(&geom, eps, conn.sign);
        // Unweighted mass (ε_r ≡ 1) for the β-direct RHS metric T_1.
        let t1_local = local_b_ee_mass(&geom, default_one, conn.sign);

        // Scatter into globals (Dirichlet-eliminated): skip rows/cols
        // whose global edge is on the PEC boundary.
        for li in 0..3 {
            let gi = conn.global_edge[li];
            let Some(ii) = interior_dof_of_edge[gi] else {
                continue;
            };
            for lj in 0..3 {
                let gj = conn.global_edge[lj];
                let Some(jj) = interior_dof_of_edge[gj] else {
                    continue;
                };
                s[(ii, jj)] += s_local[li][lj];
                t[(ii, jj)] += t_local[li][lj];
                t1[(ii, jj)] += t1_local[li][lj];
            }
        }
    }

    AssembledTransverse {
        s,
        t,
        t1,
        interior_to_global,
        eps_r_max_re,
    }
}

/// Result of [`assemble_mixed`]: the **mixed `(E_t, E_z)`**
/// Lee-Sun-Cendes (1991) block generalized eigenproblem with
/// `x = [E_t; E_z]`. The staged blocks are a **stiffness** `A` and **two**
/// mass matrices — the ε_r-weighted [`Self::b`] and the unweighted
/// [`Self::b1`] — which [`super::solve_dense_mixed`] combines into the
/// β-direct pencil `(k_0² B − A) x = β² B_1 x`.
///
/// **Eigenvalue convention (Phase 1.3.1.1 step 5.2 — β-direct).** The
/// solver now extracts `β²` as the **direct eigenvalue** of
/// `(k_0² B − A) x = β² B_1 x`, where `A = diag(A_tt, A_zz)` is the pure
/// curl/gradient stiffness (no `k_0²` term — see [`local_a_zz`] /
/// [`local_b_zz`]), `B` is the ε_r-weighted block-mass + coupling, and
/// `B_1` is the **unweighted** block-mass + coupling. This replaces the
/// step-5 cutoff parameterization `A x = k_c² B x` followed by
/// `β² = k_0² − k_c²` (vacuum `k_0`), which was algebraically correct only
/// for `ε_r ≡ 1` (where `B_1 ≡ B`, so the two pencils coincide and the
/// dominant-mode β is unchanged — the DoD-V1 canary). For `ε_r ≠ 1` the
/// old extraction under-counted the dielectric; the β-direct form puts ε_r
/// on the `k_0² B` side and the `−β²` (unweighted-mass) term on the RHS,
/// matching the transverse Helmholtz equation
/// `∇×(1/μ_r ∇×E_t) = (k_0² ε_r − β²) E_t`. The pencil remains
/// **frequency-independent** in `A`, `B`, `B_1`; `k_0²` enters only at the
/// solve.
///
/// **Block layout** (edge DoFs stacked above vertex DoFs):
///
/// ```text
///     ┌ A_tt   0  ┐ ┌E_t┐         ┌ B_tt   B_tz ┐ ┌E_t┐
///     │           │ │   │ = k_c²  │             │ │   │
///     └  0   A_zz ┘ └E_z┘         └ B_zt   B_zz ┘ └E_z┘
/// ```
///
/// with the blocks from the staged element matrices:
/// * `A_tt` = curl-curl stiffness ([`local_a_ee_curl`], `∫(1/μ_r)(∇×N)(∇×N)`);
/// * `A_zz` = nodal gradient stiffness ([`local_a_zz`], `∫(1/μ_r)∇L·∇L`);
/// * `B_tt` = Nedelec mass ([`local_b_ee_mass`], `∫ε_r N·N`);
/// * `B_zz` = nodal mass ([`local_b_zz`], `∫ε_r L·L`);
/// * `B_tz = B_ztᵀ` = edge-node coupling from [`local_b_ze`]
///   (`B_ze[i_vert][j_edge] = ∫(1/μ_r) ∇L_i·N_j`, so the **edge-row /
///   vertex-col** block `B_tz` is the transpose of the assembled
///   `B_ze`, and `B_zt` is `B_ze` itself). The coupling carries the
///   `1/μ_r` weight (the curl-curl cross term), **not** the `ε_r` weight
///   of the originally-staged matrix — see [`local_b_ze`] for the
///   step-5-review correction and why the `ε_r` weight was inert.
///
/// **Homogeneous decoupling.** On a homogeneous guide with
/// `μ_r = ε_r = 1` the coupling acting on the transverse eigenvector
/// vanishes: `B_zt e_t = ∫(1/μ_r)∇L_j·E_t dA = ∫∇L_j·E_t dA`, and the
/// Nedelec curl-curl eigenvector is exactly orthogonal to nodal
/// gradients in the `T = ∫ε_r N·N` inner product
/// (Boffi-Brezzi-Demkowicz), which for `ε_r = 1` is the same integral.
/// So `B_zt e_t = 0` and `[e_t; 0]` solves the mixed pencil with the
/// *same* `k_c²` as the transverse pencil (the DoD-V1 regression canary).
/// On an **inhomogeneous** guide (`ε_r` varying, `μ_r = 1`) the `1/μ_r`
/// coupling is **no longer** annihilated — `∫∇L·E_t ≠ 0` because `E_t`
/// is `ε_r`-orthogonal, not plain-orthogonal, to gradients — so the
/// dominant mode genuinely develops `E_z ≠ 0` and the coupling block is
/// load-bearing (guarded by the horizontal-slab test). NB: the canary
/// alone guards only `A_tt`/`A_zz`/`B_tt`/`B_zz` placement and the
/// vertex elimination; the **coupling** sign/scale/transpose is guarded
/// by the dedicated `local_b_ze` pin + the load-bearing tests, since on
/// the homogeneous guide the coupling never multiplies a nonzero `E_z`.
pub(crate) struct AssembledMixed {
    /// Block-stiffness matrix `A = diag(A_tt, A_zz)`, size `n × n` with
    /// `n = n_t + n_z`. Edge DoFs occupy `0..n_t`, vertex DoFs `n_t..n`.
    pub a: DMatrix<Complex64>,
    /// ε_r-weighted block-mass matrix `B = [[B_tt, B_tz], [B_zt, B_zz]]`,
    /// size `n × n`. `B_tt = ∫ε_r N·N`, `B_zz = ∫ε_r L·L`; the off-diagonal
    /// coupling `B_tz = B_ztᵀ = ∫(1/μ_r)∇L·N` is ε-independent.
    pub b: DMatrix<Complex64>,
    /// **Unweighted** block-mass matrix `B_1`, identical to [`Self::b`]
    /// except the diagonal mass blocks use `ε_r ≡ 1`
    /// (`B_tt → ∫N·N`, `B_zz → ∫L·L`); the `1/μ_r`-weighted coupling block
    /// is ε-independent and therefore unchanged. This is the RHS metric of
    /// the β-direct mixed pencil `(k_0² B − A) x = β² B_1 x` (Phase 1.3.1.1
    /// step 5.2): the `−β²` term of the transverse vector-Helmholtz system
    /// carries no ε_r, so its mass is unweighted. On a homogeneous guide
    /// (`ε_r ≡ 1`) `B_1 ≡ B` and the pencil reduces to the step-5 form
    /// (`β² = k_0² − k_c²`), preserving the DoD-V1 canary; for `ε_r ≠ 1`
    /// the two differ and the β-direct form removes the ε_r=1-only bias of
    /// the old `β² = k_0² − k_c²` extraction.
    pub b1: DMatrix<Complex64>,
    /// Map from interior-edge DoF index (`0..n_t`) to the global edge
    /// index, for scattering the `E_t` eigenvector components back onto
    /// the mesh. Identical in meaning to
    /// [`AssembledTransverse::interior_to_global`].
    pub interior_to_global_edges: Vec<usize>,
    /// Map from interior-vertex DoF index (`0..n_z`) to the global
    /// vertex index, for scattering the `E_z` eigenvector components
    /// back onto the mesh.
    pub interior_to_global_verts: Vec<usize>,
    /// Number of interior-edge (`E_t`) DoFs. Equals
    /// `interior_to_global_edges.len()`.
    pub n_t: usize,
    /// Number of interior-vertex (`E_z`) DoFs. Equals
    /// `interior_to_global_verts.len()`.
    pub n_z: usize,
}

/// Assemble the mixed `(E_t, E_z)` Lee-Sun-Cendes generalized
/// eigenproblem `A x = k_c² B x` on the supplied mesh.
///
/// The transverse `E_t` block reuses the same Nedelec stiffness / mass
/// as [`assemble_transverse`]; the longitudinal `E_z` block adds the
/// staged nodal-Lagrange matrices ([`local_a_zz`], [`local_b_zz`]) plus
/// the edge-node coupling ([`local_b_ze`]). See [`AssembledMixed`] for
/// the block layout, eigenvalue convention, and the homogeneous-
/// decoupling property.
///
/// PEC walls impose homogeneous Dirichlet on **both** the tangential
/// `E_t` (boundary-edge elimination, as in [`assemble_transverse`])
/// **and** the longitudinal `E_z` (boundary-vertex elimination, via
/// [`super::mesh::EdgeTable::boundary_vertices`]). Material data is
/// keyed by [`MaterialTag`] and looked up per triangle; missing keys
/// default to vacuum.
pub(crate) fn assemble_mixed(
    mesh: &TriMesh2D,
    eps_r: &HashMap<MaterialTag, Complex64>,
    mu_r: &HashMap<MaterialTag, Complex64>,
    edge_table: &EdgeTable,
) -> AssembledMixed {
    // Interior-edge DoF numbering (mirrors `assemble_transverse`).
    let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; edge_table.n_edges()];
    let mut interior_to_global_edges: Vec<usize> = Vec::new();
    for (gid, &is_bnd) in edge_table.is_boundary.iter().enumerate() {
        if !is_bnd {
            interior_dof_of_edge[gid] = Some(interior_to_global_edges.len());
            interior_to_global_edges.push(gid);
        }
    }
    let n_t = interior_to_global_edges.len();

    // Interior-vertex DoF numbering: drop PEC boundary vertices, mirror
    // of the edge elimination above.
    let boundary_vertex = edge_table.boundary_vertices(mesh.n_verts());
    let mut interior_dof_of_vert: Vec<Option<usize>> = vec![None; mesh.n_verts()];
    let mut interior_to_global_verts: Vec<usize> = Vec::new();
    for (vid, &is_bnd) in boundary_vertex.iter().enumerate() {
        if !is_bnd {
            interior_dof_of_vert[vid] = Some(interior_to_global_verts.len());
            interior_to_global_verts.push(vid);
        }
    }
    let n_z = interior_to_global_verts.len();
    let n = n_t + n_z;

    let zero = Complex64::new(0.0, 0.0);
    let mut a = DMatrix::from_element(n, n, zero);
    let mut b = DMatrix::from_element(n, n, zero);
    // Unweighted block-mass B_1: same structure as B but the diagonal mass
    // blocks use ε_r ≡ 1. The β-direct RHS metric (Phase 1.3.1.1 step 5.2).
    let mut b1 = DMatrix::from_element(n, n, zero);

    let default_one = Complex64::new(1.0, 0.0);

    for (tri_idx, conn) in edge_table.tri_edges.iter().enumerate() {
        let geom = TriGeom::from_mesh(mesh, tri_idx);
        let tag = mesh.triangle_material[tri_idx];
        let eps = *eps_r.get(&tag).unwrap_or(&default_one);
        let mu = *mu_r.get(&tag).unwrap_or(&default_one);

        let a_tt = local_a_ee_curl(&geom, mu, conn.sign);
        let b_tt = local_b_ee_mass(&geom, eps, conn.sign);
        let a_zz = local_a_zz(&geom, mu);
        let b_zz = local_b_zz(&geom, eps);
        // Unweighted (ε_r ≡ 1) mass blocks for B_1.
        let b_tt1 = local_b_ee_mass(&geom, default_one, conn.sign);
        let b_zz1 = local_b_zz(&geom, default_one);
        // B_ze[i_vert][j_edge] = ∫ (1/μ_r) ∇L_i · N_j (vertex-row /
        // edge-col). The 1/μ_r weight (matching A_zz) is the load-bearing
        // curl-curl coupling; see `local_b_ze`'s docstring for why the
        // ε_r weight was inert. The coupling is ε-independent, so it is the
        // SAME in B and B_1.
        let b_ze = local_b_ze(&geom, mu, conn.sign);

        let tri = mesh.triangles[tri_idx];

        // --- transverse-transverse block (edge / edge) ---
        for li in 0..3 {
            let Some(ii) = interior_dof_of_edge[conn.global_edge[li]] else {
                continue;
            };
            for lj in 0..3 {
                let Some(jj) = interior_dof_of_edge[conn.global_edge[lj]] else {
                    continue;
                };
                a[(ii, jj)] += a_tt[li][lj];
                b[(ii, jj)] += b_tt[li][lj];
                b1[(ii, jj)] += b_tt1[li][lj];
            }
        }

        // --- longitudinal-longitudinal block (vertex / vertex) ---
        for li in 0..3 {
            let Some(ii) = interior_dof_of_vert[tri[li]] else {
                continue;
            };
            let gi = n_t + ii;
            for lj in 0..3 {
                let Some(jj) = interior_dof_of_vert[tri[lj]] else {
                    continue;
                };
                let gj = n_t + jj;
                a[(gi, gj)] += a_zz[li][lj];
                b[(gi, gj)] += b_zz[li][lj];
                b1[(gi, gj)] += b_zz1[li][lj];
            }
        }

        // --- coupling blocks (edge ↔ vertex), mass side only ---
        // B_ze[lv][le] sits at global (vertex-row n_t+iv, edge-col ie):
        // that is the B_zt block. Its transpose populates B_tz at
        // (edge-row ie, vertex-col n_t+iv). Assembling both halves keeps
        // B symmetric (B_tz = B_ztᵀ), which the solve requires. The
        // coupling is ε-independent, so it populates BOTH B and B_1
        // identically (on a homogeneous guide B_1 ≡ B, preserving the
        // canary; the coupling is the curl-curl cross term, present on the
        // β-direct RHS through B_1's coupling sub-block).
        for lv in 0..3 {
            let Some(iv) = interior_dof_of_vert[tri[lv]] else {
                continue;
            };
            let row_v = n_t + iv;
            for le in 0..3 {
                let Some(ie) = interior_dof_of_edge[conn.global_edge[le]] else {
                    continue;
                };
                let c = b_ze[lv][le];
                b[(row_v, ie)] += c; // B_zt
                b[(ie, row_v)] += c; // B_tz = B_ztᵀ
                b1[(row_v, ie)] += c;
                b1[(ie, row_v)] += c;
            }
        }
    }

    AssembledMixed {
        a,
        b,
        b1,
        interior_to_global_edges,
        interior_to_global_verts,
        n_t,
        n_z,
    }
}

/// Assemble the **second-order (p=2)** mixed `(E_t, E_z)` Lee-Sun-Cendes
/// block generalized eigenproblem on the supplied mesh (Phase 1.3.1.1 step
/// 5.5).
///
/// Produces the **same** [`AssembledMixed`] `(A, B, B₁)` block structure as
/// the first-order [`assemble_mixed`] — just a larger pencil (`n_t` = 2 ×
/// interior-edge + 2 × triangle interior; `n_z` = interior-vertex +
/// interior-edge-midpoint) — so the order-agnostic
/// [`super::solve_dense_mixed`] consumes it unchanged. The transverse block
/// uses the p=2 Nedelec-first-kind order-2 element matrices
/// ([`local_a_ee_curl_p2`] / [`local_b_ee_mass_p2`]); the longitudinal block
/// uses quadratic nodal-Lagrange ([`local_a_zz_p2`] / [`local_b_zz_p2`]);
/// the coupling is the p=2 [`local_b_ze_p2`] (vertex/midpoint-row,
/// edge/interior-col), carrying the `1/μ_r` weight exactly as at first order.
///
/// The DoF bookkeeping (edge 2-DoF orientation, interior face DoFs,
/// midpoint nodes, PEC elimination) is owned by [`super::mesh::P2DofMap`].
/// **First-order stays the default** ([`assemble_mixed`]); this path is
/// selected only via [`super::ElementOrder::Second`] for the high-contrast
/// inhomogeneous case.
///
/// **`interior_to_global_*` caveat.** At p=2 the transverse interior (face)
/// DoFs have no owning global edge and the nodal midpoint DoFs no owning
/// vertex; their entries in [`AssembledMixed::interior_to_global_edges`] /
/// [`AssembledMixed::interior_to_global_verts`] are `usize::MAX` sentinels.
/// The β² eigenvalue is fully meaningful (the only thing the step-5.5 gates
/// consume); the first-order edge-scatter field reconstruction in
/// [`crate::ports::NumericalCrossSection::solve`] is **not** wired for p=2
/// (that path always calls the first-order [`assemble_mixed`]).
#[allow(dead_code)] // selected via ElementOrder::Second; consumed by the J3/J4 lib tests
pub(crate) fn assemble_mixed_p2(
    mesh: &TriMesh2D,
    eps_r: &HashMap<MaterialTag, Complex64>,
    mu_r: &HashMap<MaterialTag, Complex64>,
    edge_table: &EdgeTable,
) -> AssembledMixed {
    let dofs = super::mesh::P2DofMap::build(mesh, edge_table);
    let n_t = dofs.n_t;
    let n_z = dofs.n_z;
    let n = n_t + n_z;

    let zero = Complex64::new(0.0, 0.0);
    let mut a = DMatrix::from_element(n, n, zero);
    let mut b = DMatrix::from_element(n, n, zero);
    let mut b1 = DMatrix::from_element(n, n, zero);

    let default_one = Complex64::new(1.0, 0.0);

    for (tri_idx, conn) in edge_table.tri_edges.iter().enumerate() {
        let geom = TriGeom::from_mesh(mesh, tri_idx);
        let tag = mesh.triangle_material[tri_idx];
        let eps = *eps_r.get(&tag).unwrap_or(&default_one);
        let mu = *mu_r.get(&tag).unwrap_or(&default_one);

        let a_tt = local_a_ee_curl_p2(&geom, mu, conn.sign);
        let b_tt = local_b_ee_mass_p2(&geom, eps, conn.sign);
        let b_tt1 = local_b_ee_mass_p2(&geom, default_one, conn.sign);
        let a_zz = local_a_zz_p2(&geom, mu);
        let b_zz = local_b_zz_p2(&geom, eps);
        let b_zz1 = local_b_zz_p2(&geom, default_one);
        // B_ze_p2[i_node][j_edge] = ∫(1/μ_r) ∇L_i · N_j (ε-independent).
        let b_ze = local_b_ze_p2(&geom, mu, conn.sign);

        let tdof = &dofs.tri_t_dofs[tri_idx];
        let zdof = &dofs.tri_z_dofs[tri_idx];

        // --- transverse-transverse block (8×8) ---
        for li in 0..8 {
            let Some(ii) = tdof[li] else { continue };
            for lj in 0..8 {
                let Some(jj) = tdof[lj] else { continue };
                a[(ii, jj)] += a_tt[li][lj];
                b[(ii, jj)] += b_tt[li][lj];
                b1[(ii, jj)] += b_tt1[li][lj];
            }
        }

        // --- longitudinal-longitudinal block (6×6, vertex/midpoint nodes) ---
        for li in 0..6 {
            let Some(ii) = zdof[li] else { continue };
            let gi = n_t + ii;
            for lj in 0..6 {
                let Some(jj) = zdof[lj] else { continue };
                let gj = n_t + jj;
                a[(gi, gj)] += a_zz[li][lj];
                b[(gi, gj)] += b_zz[li][lj];
                b1[(gi, gj)] += b_zz1[li][lj];
            }
        }

        // --- coupling blocks (node ↔ Nedelec); assemble both halves so B
        // stays symmetric (B_tz = B_ztᵀ), as the solve requires. The coupling
        // is ε-independent → populates B and B_1 identically. ---
        for li in 0..6 {
            let Some(iv) = zdof[li] else { continue };
            let row_v = n_t + iv;
            for lj in 0..8 {
                let Some(je) = tdof[lj] else { continue };
                let c = b_ze[li][lj];
                b[(row_v, je)] += c; // B_zt
                b[(je, row_v)] += c; // B_tz = B_ztᵀ
                b1[(row_v, je)] += c;
                b1[(je, row_v)] += c;
            }
        }
    }

    AssembledMixed {
        a,
        b,
        b1,
        interior_to_global_edges: dofs.t_dof_edge,
        interior_to_global_verts: dofs.z_dof_vert,
        n_t,
        n_z,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigensolver::mesh::EdgeTable;

    fn unit_right_triangle() -> TriMesh2D {
        // Single CCW triangle with legs 1.0 / 1.0 — area = 0.5.
        TriMesh2D::new(
            vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            vec![[0, 1, 2]],
            None,
            None,
        )
        .unwrap()
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn local_curl_matrix_symmetric_on_unit_triangle() {
        let mesh = unit_right_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let s = local_a_ee_curl(&geom, Complex64::new(1.0, 0.0), table.tri_edges[0].sign);
        // Real-symmetric for the lossless case.
        for i in 0..3 {
            for j in 0..3 {
                let diff = (s[i][j] - s[j][i]).norm();
                assert!(
                    diff < 1e-12,
                    "S[{i},{j}] = {:?}, S[{j},{i}] = {:?}",
                    s[i][j],
                    s[j][i]
                );
            }
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn local_mass_matrix_symmetric_on_unit_triangle() {
        let mesh = unit_right_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let t = local_b_ee_mass(&geom, Complex64::new(1.0, 0.0), table.tri_edges[0].sign);
        for i in 0..3 {
            for j in 0..3 {
                let diff = (t[i][j] - t[j][i]).norm();
                assert!(diff < 1e-12, "T not symmetric at ({i},{j}): diff={diff}");
            }
        }
        // Diagonal entries must be strictly positive (it's a Gram matrix).
        for i in 0..3 {
            assert!(
                t[i][i].re > 0.0,
                "T[{i},{i}] = {:?} should have positive real part",
                t[i][i]
            );
        }
    }

    #[test]
    fn assembled_s_t_dimensions_match_interior_edges() {
        // Two-triangle unit square: 5 edges, 4 boundary, 1 interior.
        let mesh = TriMesh2D::new(
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            vec![[0, 1, 2], [0, 2, 3]],
            None,
            None,
        )
        .unwrap();
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = assemble_transverse(&mesh, &eps, &mu, &table);
        assert_eq!(asm.s.nrows(), 1);
        assert_eq!(asm.s.ncols(), 1);
        assert_eq!(asm.t.nrows(), 1);
        assert_eq!(asm.t.ncols(), 1);
        assert_eq!(asm.interior_to_global.len(), 1);
    }

    #[test]
    fn local_a_zz_and_b_zz_symmetric_positive_diagonal() {
        // Architectural smoke test for the longitudinal-block helpers.
        // They are unused by the WR-90 TE10 path but are part of the
        // mixed-formulation contract with the spec.
        let mesh = unit_right_triangle();
        let geom = TriGeom::from_mesh(&mesh, 0);
        let azz = local_a_zz(&geom, Complex64::new(1.0, 0.0));
        let bzz = local_b_zz(&geom, Complex64::new(1.0, 0.0));
        for i in 0..3 {
            for j in 0..3 {
                assert!((azz[i][j] - azz[j][i]).norm() < 1e-12);
                assert!((bzz[i][j] - bzz[j][i]).norm() < 1e-12);
            }
            // bzz is a Gram matrix; azz is positive-semidefinite (only
            // constant-λ kernel is in the nullspace), but its on-diagonal
            // entries are still strictly positive.
            assert!(bzz[i][i].re > 0.0);
            assert!(azz[i][i].re > 0.0);
        }
    }

    /// Structured `nx × ny` quad-grid of CCW triangles spanning
    /// `[0, a] × [0, b]` — the same fixture the WR-90 gate uses. All
    /// triangles share material tag `tag`.
    fn rectangular_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
        let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                vertices.push([a * (i as f64) / (nx as f64), b * (j as f64) / (ny as f64)]);
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut triangles = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let v00 = idx(i, j);
                let v10 = idx(i + 1, j);
                let v11 = idx(i + 1, j + 1);
                let v01 = idx(i, j + 1);
                triangles.push([v00, v10, v11]);
                triangles.push([v00, v11, v01]);
            }
        }
        TriMesh2D::new(vertices, triangles, None, None).unwrap()
    }

    #[test]
    fn assemble_mixed_dimensions_and_symmetry() {
        // A 4×4-quad WR-90-shaped mesh: 25 vertices, perimeter = 16 of
        // them are boundary → 9 interior vertices; interior edges from
        // the transverse path. The mixed pencil is (n_t + n_z) square
        // and both A, B are symmetric (real-symmetric for the lossless
        // case).
        let a = 22.86e-3;
        let b = 10.16e-3;
        let mesh = rectangular_mesh(a, b, 4, 4);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = assemble_mixed(&mesh, &eps, &mu, &table);

        // 5×5 vertex grid → 9 interior vertices.
        assert_eq!(asm.n_z, 9, "interior-vertex count");
        let n = asm.n_t + asm.n_z;
        assert_eq!(asm.a.nrows(), n);
        assert_eq!(asm.a.ncols(), n);
        assert_eq!(asm.b.nrows(), n);
        assert_eq!(asm.b.ncols(), n);
        assert_eq!(asm.interior_to_global_edges.len(), asm.n_t);
        assert_eq!(asm.interior_to_global_verts.len(), asm.n_z);

        // Real-symmetric (Gram structure; the coupling is assembled as
        // B_tz = B_ztᵀ explicitly). Both mass matrices B and the
        // unweighted B_1 must be symmetric.
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (asm.a[(i, j)] - asm.a[(j, i)]).norm() < 1e-12,
                    "A not symmetric at ({i},{j})"
                );
                assert!(
                    (asm.b[(i, j)] - asm.b[(j, i)]).norm() < 1e-12,
                    "B not symmetric at ({i},{j})"
                );
                assert!(
                    (asm.b1[(i, j)] - asm.b1[(j, i)]).norm() < 1e-12,
                    "B_1 not symmetric at ({i},{j})"
                );
            }
        }

        // On this homogeneous (air, ε_r ≡ 1) mesh the unweighted B_1 must
        // equal the ε_r-weighted B bit-for-bit — the structural guarantee
        // that the β-direct pencil reduces to the step-5 form on a
        // homogeneous guide (DoD-V1 canary at the matrix level).
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (asm.b[(i, j)] - asm.b1[(i, j)]).norm() < 1e-12,
                    "homogeneous B and B_1 must coincide at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn b1_unweighted_differs_from_b_under_dielectric_loading() {
        // The β-direct fix hinges on B_1 (unweighted mass) differing from B
        // (ε_r-weighted) when ε_r ≠ 1. On a dielectric-loaded mesh the mass
        // diagonal blocks must differ; the ε-independent coupling must not.
        let a = 22.86e-3;
        let b = 10.16e-3;
        let mesh = horizontal_slab_mesh(a, b, 6, 6);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        eps.insert(1u32, Complex64::new(10.2, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        mu.insert(1u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = assemble_mixed(&mesh, &eps, &mu, &table);
        let n_t = asm.n_t;
        let n = n_t + asm.n_z;

        // Some diagonal mass entry in the dielectric region must differ
        // (B carries ε_r = 10.2, B_1 carries 1).
        let mut any_mass_diff = false;
        for i in 0..n {
            if (asm.b[(i, i)] - asm.b1[(i, i)]).norm() > 1e-9 {
                any_mass_diff = true;
                break;
            }
        }
        assert!(
            any_mass_diff,
            "B_1 must differ from B on the diagonal under dielectric loading"
        );

        // The off-diagonal edge↔vertex coupling block is ε-independent, so
        // B and B_1 must agree there exactly.
        for i in 0..n_t {
            for j in n_t..n {
                assert!(
                    (asm.b[(i, j)] - asm.b1[(i, j)]).norm() < 1e-12,
                    "coupling block must be identical in B and B_1 at ({i},{j})"
                );
            }
        }
    }

    /// Horizontal-slab WR-90 mesh: lower-y half tagged 1, rest tagged 0.
    /// Mirrors the integration-test fixture; used by the B_1 contrast test.
    fn horizontal_slab_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
        let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                vertices.push([a * (i as f64) / (nx as f64), b * (j as f64) / (ny as f64)]);
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut triangles = Vec::with_capacity(2 * nx * ny);
        let mut tags = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let v00 = idx(i, j);
                let v10 = idx(i + 1, j);
                let v11 = idx(i + 1, j + 1);
                let v01 = idx(i, j + 1);
                let yc = b * ((j as f64) + 0.5) / (ny as f64);
                let tag = if yc < b / 2.0 { 1u32 } else { 0u32 };
                triangles.push([v00, v10, v11]);
                tags.push(tag);
                triangles.push([v00, v11, v01]);
                tags.push(tag);
            }
        }
        TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
    }

    #[test]
    fn local_b_ze_column_sums_vanish_partition_of_unity() {
        // The per-element coupling B_ze[i_vert][j_edge] = ∫(1/μ_r) ∇L_i·N_j
        // does NOT vanish entry-wise, but its column sums (Σ_i over the
        // three vertices, for each edge j) vanish because
        // Σ_i ∇L_i = ∇(Σ_i L_i) = ∇(1) = 0. This is the element-local
        // partition-of-unity fingerprint, NOT a statement that the
        // coupling is inert: globally the 1/μ_r coupling is load-bearing
        // on inhomogeneous guides (see the horizontal-slab test). The
        // earlier ε_r-weighted coupling WAS globally inert for the
        // dominant mode (annihilated by the divergence-free transverse
        // eigenvector) — the step-5-review bug this weight fix corrects.
        let mesh = unit_right_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let bze = local_b_ze(&geom, Complex64::new(1.0, 0.0), table.tri_edges[0].sign);
        // Not all-zero entry-wise.
        let any_nonzero = bze.iter().flatten().any(|z| z.norm() > 1e-9);
        assert!(any_nonzero, "B_ze should have nonzero entries");
        // Column sums vanish: Σ_i B_ze[i][j] = 0 for each edge-column j.
        let mut col_sum = [Complex64::new(0.0, 0.0); 3];
        for row in &bze {
            for (j, entry) in row.iter().enumerate() {
                col_sum[j] += *entry;
            }
        }
        for (j, sum) in col_sum.iter().enumerate() {
            assert!(
                sum.norm() < 1e-12,
                "B_ze column {j} sum = {sum:?} should vanish (partition of unity)"
            );
        }
    }

    #[test]
    fn local_b_ze_coupling_smoke_test() {
        // Smoke test: just confirms the helper compiles and returns
        // finite values. The coupling block is unused by the WR-90
        // validation path.
        let mesh = unit_right_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let bze = local_b_ze(&geom, Complex64::new(1.0, 0.0), table.tri_edges[0].sign);
        for row in &bze {
            for entry in row {
                assert!(entry.re.is_finite() && entry.im.is_finite());
            }
        }
    }

    #[test]
    fn local_b_ze_matches_independent_quadrature_sign_and_scale() {
        // Step-5-review P1-1 guard for the highest-risk item: pin every
        // entry of the coupling block `B_ze[i][j] = ∫(1/μ_r)∇L_i·N_j`
        // against an INDEPENDENT 3-point edge-midpoint quadrature (exact
        // for the linear `λ × const` integrand). A sign / scale /
        // edge-endpoint-ordering error in `local_b_ze` changes individual
        // entries and is caught here — the homogeneous β canary cannot
        // catch it because its dominant eigenvector has `E_z = 0`, so the
        // coupling never multiplies a nonzero vector.
        //
        // Uses a deliberately GENERIC triangle (non-right, non-unit,
        // mixed-sign edge orientations) so the test exercises the σ and ℓ
        // factors, not just the trivial unit-triangle symmetry.
        let mesh = TriMesh2D::new(
            vec![[0.2, 0.1], [1.3, 0.0], [0.4, 1.1]],
            vec![[0, 1, 2]],
            None,
            None,
        )
        .unwrap();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let signs = table.tri_edges[0].sign;
        // μ_r = 1 so the weight is unity and the quadrature is a pure
        // geometric reference.
        let bze = local_b_ze(&geom, Complex64::new(1.0, 0.0), signs);

        let v: [[f64; 2]; 3] = [mesh.vertices[0], mesh.vertices[1], mesh.vertices[2]];
        let area = geom.area;
        // ∇λ_i = (b_i, c_i)/(2A) — recompute independently.
        let grad = |i: usize| -> [f64; 2] {
            let i1 = (i + 1) % 3;
            let i2 = (i + 2) % 3;
            let b = v[i1][1] - v[i2][1];
            let c = v[i2][0] - v[i1][0];
            [b / (2.0 * area), c / (2.0 * area)]
        };
        // Barycentric coords of an arbitrary point via sub-triangle areas.
        let lam_at = |p: [f64; 2]| -> [f64; 3] {
            let sub = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
                0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
            };
            [
                sub(p, v[1], v[2]) / area,
                sub(v[0], p, v[2]) / area,
                sub(v[0], v[1], p) / area,
            ]
        };
        // Nedelec basis N_j at point p.
        let n_at = |p: [f64; 2], j: usize| -> [f64; 2] {
            let [a, b] = LOCAL_EDGES[j];
            let lam = lam_at(p);
            let ga = grad(a);
            let gb = grad(b);
            let s = signs[j] * geom.edge_len[j];
            [
                s * (lam[a] * gb[0] - lam[b] * ga[0]),
                s * (lam[a] * gb[1] - lam[b] * ga[1]),
            ]
        };
        // 3-point midpoint rule (weight A/3 at each edge midpoint).
        let mids: [[f64; 2]; 3] = [
            [0.5 * (v[1][0] + v[2][0]), 0.5 * (v[1][1] + v[2][1])],
            [0.5 * (v[2][0] + v[0][0]), 0.5 * (v[2][1] + v[0][1])],
            [0.5 * (v[0][0] + v[1][0]), 0.5 * (v[0][1] + v[1][1])],
        ];
        for (i, row) in bze.iter().enumerate() {
            let gi = grad(i);
            for (j, entry) in row.iter().enumerate() {
                let mut quad = 0.0;
                for &mp in &mids {
                    let nj = n_at(mp, j);
                    quad += (area / 3.0) * (gi[0] * nj[0] + gi[1] * nj[1]);
                }
                assert!(
                    (entry.re - quad).abs() < 1e-12,
                    "B_ze[{i}][{j}] = {} disagrees with independent quadrature {quad}",
                    entry.re
                );
                assert!(entry.im.abs() < 1e-15, "lossless → real");
            }
        }
    }

    // ── Phase 1.3.1.1 step 5.5 — p=2 element matrices (J1) ──────────────
    //
    // These are the J1 correctness anchors: each p=2 element matrix is pinned
    // against an INDEPENDENT, higher-order quadrature evaluation of the same
    // point-wise basis (mirroring the first-order `local_b_ze` pin), plus the
    // basis is checked for completeness (rank 8) and tangential conformity
    // (interior functions vanish tangentially on ∂T). A wrong quadrature,
    // sign, or basis is caught here BEFORE any eigensolve (the J3/J4 ladder).

    /// Generic non-right, non-unit triangle exercising the σ and ℓ factors.
    fn generic_triangle() -> TriMesh2D {
        TriMesh2D::new(
            vec![[0.2, 0.1], [1.3, 0.0], [0.4, 1.1]],
            vec![[0, 1, 2]],
            None,
            None,
        )
        .unwrap()
    }

    /// An INDEPENDENT degree-5, 7-point symmetric triangle quadrature
    /// (Dunavant 1985) — a different rule from the degree-4 6-point
    /// production [`tri_gauss_deg4`] (different points, weights, and count),
    /// exact to total degree 5 (one above the highest p=2 integrand degree
    /// of 4). Agreement of a degree-4 production matrix with this degree-5
    /// reference certifies the production rule rather than re-deriving it.
    /// Returns `(λ₀, λ₁, λ₂, weight)`, weights summing to 1.
    fn independent_tri_quad() -> Vec<(f64, f64, f64, f64)> {
        let w0 = 0.225;
        let a1 = 0.470_142_064_105_115;
        let w1 = 0.132_394_152_788_506;
        let a2 = 0.101_286_507_323_456;
        let w2 = 0.125_939_180_544_827;
        vec![
            (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, w0),
            (a1, a1, 1.0 - 2.0 * a1, w1),
            (a1, 1.0 - 2.0 * a1, a1, w1),
            (1.0 - 2.0 * a1, a1, a1, w1),
            (a2, a2, 1.0 - 2.0 * a2, w2),
            (a2, 1.0 - 2.0 * a2, a2, w2),
            (1.0 - 2.0 * a2, a2, a2, w2),
        ]
    }

    #[test]
    fn tri_gauss_deg4_integrates_quartics_exactly() {
        // Pin the hand-rolled degree-4 rule: every barycentric monomial
        // λ₀^p λ₁^q λ₂^r with p+q+r ≤ 4 integrates to the exact reference-
        // triangle value ∫ λ₀^p λ₁^q λ₂^r dA / A = p! q! r! · 2 / (p+q+r+2)!
        // (weights normalised to sum 1, so the area factor drops out).
        let rule = tri_gauss_deg4();
        let wsum: f64 = rule.iter().map(|&(_, _, _, w)| w).sum();
        assert!((wsum - 1.0).abs() < 1e-14, "weights must sum to 1: {wsum}");
        let fact = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0];
        for p in 0..=4usize {
            for q in 0..=(4 - p) {
                let r = 4 - p - q; // test the top total-degree shell (=4)
                let exact = fact[p] * fact[q] * fact[r] * 2.0 / fact[p + q + r + 2];
                let approx: f64 = rule
                    .iter()
                    .map(|&(l0, l1, l2, w)| {
                        w * l0.powi(p as i32) * l1.powi(q as i32) * l2.powi(r as i32)
                    })
                    .sum();
                assert!(
                    (approx - exact).abs() < 1e-13,
                    "deg-4 rule wrong on λ0^{p}λ1^{q}λ2^{r}: got {approx}, exact {exact}"
                );
            }
        }
    }

    #[test]
    fn p2_basis_has_full_rank_eight() {
        // Completeness (risk (a)/(b) mitigation): the 8 p=2 Nedelec basis
        // vectors must be linearly independent (span N₁(order 2), dim 8).
        // Sample them at many interior points, stack the 2-component values
        // into a tall (2·npts)×8 matrix, and assert its singular values are
        // all bounded away from zero (rank 8). A wrong / degenerate interior
        // pair collapses the rank and is caught here, before any eigensolve.
        use nalgebra::DMatrix;
        let mesh = generic_triangle();
        let geom = TriGeom::from_mesh(&mesh, 0);
        let g = bary_grads(&geom);
        // 10 scattered interior barycentric points.
        let pts = [
            [0.5, 0.3, 0.2],
            [0.2, 0.5, 0.3],
            [0.3, 0.2, 0.5],
            [0.6, 0.1, 0.3],
            [0.1, 0.6, 0.3],
            [0.3, 0.1, 0.6],
            [0.34, 0.33, 0.33],
            [0.7, 0.2, 0.1],
            [0.15, 0.25, 0.6],
            [0.45, 0.45, 0.1],
        ];
        let mut m = DMatrix::<f64>::zeros(2 * pts.len(), 8);
        for (pi, lam) in pts.iter().enumerate() {
            let n = p2_nedelec_basis(*lam, g);
            for (j, nj) in n.iter().enumerate() {
                m[(2 * pi, j)] = nj[0];
                m[(2 * pi + 1, j)] = nj[1];
            }
        }
        let svd = m.singular_values();
        let smin = svd.iter().cloned().fold(f64::INFINITY, f64::min);
        let smax = svd.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            smin > 1e-9 * smax,
            "p=2 Nedelec basis is rank-deficient (σ_min/σ_max = {:.3e}); the 8 functions \
             must be linearly independent (span N₁(2))",
            smin / smax
        );
    }

    #[test]
    fn p2_interior_functions_have_zero_tangential_trace() {
        // Conformity (risk (b) mitigation): the two interior ("face")
        // functions must have ZERO tangential component along all three
        // edges — that is what makes them purely interior (no inter-element
        // coupling, no orientation sign). Sample each edge at several
        // parameter values and assert N·t̂ ≈ 0 for slots 6 and 7. (The edge
        // functions, by contrast, have nonzero tangential trace on their own
        // edge — also checked, as a sanity counterpoint.)
        let mesh = generic_triangle();
        let v: [[f64; 2]; 3] = [mesh.vertices[0], mesh.vertices[1], mesh.vertices[2]];
        let geom = TriGeom::from_mesh(&mesh, 0);
        let g = bary_grads(&geom);
        for (e, &[a, b]) in LOCAL_EDGES.iter().enumerate() {
            let ta = [v[b][0] - v[a][0], v[b][1] - v[a][1]];
            let tlen = (ta[0] * ta[0] + ta[1] * ta[1]).sqrt();
            let that = [ta[0] / tlen, ta[1] / tlen];
            let mut edge_trace_max = 0.0_f64;
            for s in [0.1, 0.25, 0.5, 0.75, 0.9] {
                // Point on edge e: λ_a = 1−s, λ_b = s, λ_other = 0.
                let mut lam = [0.0; 3];
                lam[a] = 1.0 - s;
                lam[b] = s;
                let n = p2_nedelec_basis(lam, g);
                // Interior functions: tangential trace must vanish.
                for slot in [6usize, 7] {
                    let tan = n[slot][0] * that[0] + n[slot][1] * that[1];
                    assert!(
                        tan.abs() < 1e-12,
                        "interior function {slot} has nonzero tangential trace {tan} on edge {e}"
                    );
                }
                // The Whitney slot of this edge should NOT vanish tangentially
                // (sanity: the edge DoFs do carry the tangential field).
                let w_tan = n[2 * e][0] * that[0] + n[2 * e][1] * that[1];
                edge_trace_max = edge_trace_max.max(w_tan.abs());
            }
            assert!(
                edge_trace_max > 1e-6,
                "edge {e} Whitney function has ~zero tangential trace on its own edge — basis bug"
            );
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn p2_mass_matrix_matches_independent_quadrature() {
        // J1 PRIMARY PIN (mirrors `local_b_ze_matches_independent_quadrature`
        // at p=2): every entry of the 8×8 p=2 Nedelec mass `B_tt = ∫ε_r N·N`
        // must agree with an INDEPENDENT high-order quadrature of the same
        // point-wise basis (built differently from the production rule). A
        // quadrature, sign, or scale error changes individual entries and is
        // caught here. Generic triangle so σ/ℓ are exercised; μ irrelevant.
        let mesh = generic_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let signs = table.tri_edges[0].sign;
        let sgn = p2_edge_signs(signs);
        let b = local_b_ee_mass_p2(&geom, Complex64::new(1.0, 0.0), signs);
        let g = bary_grads(&geom);
        let quad = independent_tri_quad();
        for i in 0..8 {
            for j in 0..8 {
                let mut acc = 0.0;
                for &(l0, l1, l2, w) in &quad {
                    let n = p2_nedelec_basis([l0, l1, l2], g);
                    let dot = n[i][0] * n[j][0] + n[i][1] * n[j][1];
                    acc += w * geom.area * sgn[i] * sgn[j] * dot;
                }
                assert!(
                    (b[i][j].re - acc).abs() < 1e-12,
                    "B_tt^p2[{i}][{j}] = {} disagrees with independent quadrature {acc}",
                    b[i][j].re
                );
                assert!(b[i][j].im.abs() < 1e-15, "lossless → real");
                assert!(
                    (b[i][j] - b[j][i]).norm() < 1e-12,
                    "B_tt^p2 must be symmetric"
                );
            }
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn p2_curl_matrix_matches_independent_quadrature() {
        // J1 pin for the curl-curl stiffness `A_tt = ∫(1/μ)(∇×N)(∇×N)` at p=2
        // (curl is NON-CONSTANT — the whole reason p=2 needs quadrature).
        // Pinned against the independent rule + the closed-form constant-curl
        // sub-block (the 3 Whitney slots have constant curl 2∇λ_a×∇λ_b, the 3
        // gradient slots have zero curl).
        let mesh = generic_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let signs = table.tri_edges[0].sign;
        let sgn = p2_edge_signs(signs);
        let a = local_a_ee_curl_p2(&geom, Complex64::new(1.0, 0.0), signs);
        let g = bary_grads(&geom);
        let quad = independent_tri_quad();
        for i in 0..8 {
            for j in 0..8 {
                let mut acc = 0.0;
                for &(l0, l1, l2, w) in &quad {
                    let curl = p2_nedelec_curl([l0, l1, l2], g);
                    acc += w * geom.area * sgn[i] * sgn[j] * curl[i] * curl[j];
                }
                assert!(
                    (a[i][j].re - acc).abs() < 1e-12,
                    "A_tt^p2[{i}][{j}] = {} disagrees with independent quadrature {acc}",
                    a[i][j].re
                );
                assert!(
                    (a[i][j] - a[j][i]).norm() < 1e-12,
                    "A_tt^p2 must be symmetric"
                );
            }
            // The 3 gradient slots (2e+1) are curl-free → zero rows/cols.
            for e in 0..3 {
                assert!(
                    a[2 * e + 1][2 * e + 1].norm() < 1e-12,
                    "gradient edge function {} must have zero curl-curl self-energy",
                    2 * e + 1
                );
            }
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn p2_nodal_matrices_match_independent_quadrature() {
        // J1 pin for the quadratic nodal `A_zz = ∫(1/μ)∇L·∇L` and
        // `B_zz = ∫ε_r L·L` against the independent rule. Also pins the P2
        // nodal interpolation property `L_node(node_k) = δ` at the 6 nodes.
        let mesh = generic_triangle();
        let geom = TriGeom::from_mesh(&mesh, 0);
        let azz = local_a_zz_p2(&geom, Complex64::new(1.0, 0.0));
        let bzz = local_b_zz_p2(&geom, Complex64::new(1.0, 0.0));
        let g = bary_grads(&geom);
        let quad = independent_tri_quad();
        for i in 0..6 {
            for j in 0..6 {
                let (mut a_acc, mut b_acc) = (0.0, 0.0);
                for &(l0, l1, l2, w) in &quad {
                    let gl = p2_nodal_grad([l0, l1, l2], g);
                    let l = p2_nodal_basis([l0, l1, l2]);
                    a_acc += w * geom.area * (gl[i][0] * gl[j][0] + gl[i][1] * gl[j][1]);
                    b_acc += w * geom.area * l[i] * l[j];
                }
                assert!(
                    (azz[i][j].re - a_acc).abs() < 1e-12,
                    "A_zz^p2[{i}][{j}] = {} vs independent {a_acc}",
                    azz[i][j].re
                );
                assert!(
                    (bzz[i][j].re - b_acc).abs() < 1e-12,
                    "B_zz^p2[{i}][{j}] = {} vs independent {b_acc}",
                    bzz[i][j].re
                );
                assert!((azz[i][j] - azz[j][i]).norm() < 1e-12, "A_zz^p2 symmetric");
                assert!((bzz[i][j] - bzz[j][i]).norm() < 1e-12, "B_zz^p2 symmetric");
            }
            assert!(bzz[i][i].re > 0.0, "B_zz^p2 Gram diagonal positive");
        }
        // P2 nodal interpolation: node coordinates are the 3 vertices +
        // 3 edge midpoints. L_node(node_k) = δ_{node,k}.
        let node_bary = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.5, 0.5], // midpoint edge 0 (verts 1,2)
            [0.5, 0.0, 0.5], // midpoint edge 1 (verts 2,0)
            [0.5, 0.5, 0.0], // midpoint edge 2 (verts 0,1)
        ];
        for (k, nb) in node_bary.iter().enumerate() {
            let l = p2_nodal_basis(*nb);
            for (node, lv) in l.iter().enumerate() {
                let expect = if node == k { 1.0 } else { 0.0 };
                assert!(
                    (lv - expect).abs() < 1e-12,
                    "P2 nodal L_{node}(node {k}) = {lv}, expected {expect}"
                );
            }
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn p2_coupling_matrix_matches_independent_quadrature() {
        // J1 pin for the p=2 edge-node coupling `B_ze = ∫(1/μ)∇L_i·N_j`
        // (6 nodal rows × 8 Nedelec cols) against the independent rule —
        // the highest-risk block (sign/scale/transpose), the p=2 analogue of
        // the first-order `local_b_ze` pin.
        let mesh = generic_triangle();
        let table = EdgeTable::build(&mesh);
        let geom = TriGeom::from_mesh(&mesh, 0);
        let signs = table.tri_edges[0].sign;
        let sgn = p2_edge_signs(signs);
        let bze = local_b_ze_p2(&geom, Complex64::new(1.0, 0.0), signs);
        let g = bary_grads(&geom);
        let quad = independent_tri_quad();
        for i in 0..6 {
            for j in 0..8 {
                let mut acc = 0.0;
                for &(l0, l1, l2, w) in &quad {
                    let gl = p2_nodal_grad([l0, l1, l2], g);
                    let n = p2_nedelec_basis([l0, l1, l2], g);
                    acc += w * geom.area * sgn[j] * (gl[i][0] * n[j][0] + gl[i][1] * n[j][1]);
                }
                assert!(
                    (bze[i][j].re - acc).abs() < 1e-12,
                    "B_ze^p2[{i}][{j}] = {} disagrees with independent quadrature {acc}",
                    bze[i][j].re
                );
                assert!(bze[i][j].im.abs() < 1e-15, "lossless → real");
            }
        }
    }

    // ── Phase 1.3.1.1 step 5.5 — p=2 DoF map + assembly (J2) ────────────

    #[test]
    fn p2_dof_map_counts_two_tri_unit_square() {
        // Two-triangle unit square: 5 edges (4 boundary, 1 interior), 4
        // vertices (all on the perimeter → all boundary), 2 triangles. At
        // p=2: n_t = 2·(interior edges) + 2·(triangles) = 2·1 + 2·2 = 6;
        // n_z = (interior verts) + (interior-edge midpoints) = 0 + 1 = 1.
        let mesh = TriMesh2D::new(
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            vec![[0, 1, 2], [0, 2, 3]],
            None,
            None,
        )
        .unwrap();
        let table = EdgeTable::build(&mesh);
        let dofs = crate::eigensolver::mesh::P2DofMap::build(&mesh, &table);
        assert_eq!(dofs.n_t, 6, "p=2 transverse DoF count");
        assert_eq!(dofs.n_z, 1, "p=2 nodal DoF count");
        assert_eq!(dofs.t_dof_edge.len(), dofs.n_t);
        assert_eq!(dofs.z_dof_vert.len(), dofs.n_z);
        // Interior face DoFs (2 per triangle) must always be assigned (never
        // PEC-eliminated) — every triangle's slots 6,7 are Some.
        for td in &dofs.tri_t_dofs {
            assert!(
                td[6].is_some() && td[7].is_some(),
                "interior face DoFs always present"
            );
        }
    }

    #[test]
    fn assemble_mixed_p2_dimensions_and_symmetry() {
        // p=2 pencil on a 4×4-quad WR-90 mesh: (A, B, B_1) are (n_t+n_z)
        // square and symmetric (real-symmetric, lossless); on this
        // homogeneous (air, ε_r≡1) mesh B and B_1 must coincide bit-for-bit
        // (the matrix-level homogeneous canary — the β-direct pencil reduces
        // to the cutoff form), exactly mirroring the first-order
        // `assemble_mixed_dimensions_and_symmetry`.
        let a = 22.86e-3;
        let b = 10.16e-3;
        let mesh = rectangular_mesh(a, b, 4, 4);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = assemble_mixed_p2(&mesh, &eps, &mu, &table);
        let n = asm.n_t + asm.n_z;
        assert_eq!(asm.a.nrows(), n);
        assert_eq!(asm.b.nrows(), n);
        assert_eq!(asm.b1.nrows(), n);
        // p=2 pencil is strictly larger than the first-order one.
        let asm1 = assemble_mixed(&mesh, &eps, &mu, &table);
        assert!(
            asm.n_t > asm1.n_t && asm.n_z > asm1.n_z,
            "p=2 must add DoFs: p2 (n_t={}, n_z={}) vs p1 (n_t={}, n_z={})",
            asm.n_t,
            asm.n_z,
            asm1.n_t,
            asm1.n_z
        );
        // Symmetry tolerances are relative to the matrix scale (the curl
        // block reaches ~1/(μ·A)·ℓ² ≈ 1e5 on WR-90, so an absolute 1e-18 is
        // unphysical). B ≡ B_1 on a homogeneous mesh is an exact (same code
        // path, ε_r ≡ 1) equality, held to a tight scale-relative bound.
        let scale_a = asm
            .a
            .iter()
            .map(|z| z.norm())
            .fold(0.0_f64, f64::max)
            .max(1.0);
        let scale_b = asm
            .b
            .iter()
            .map(|z| z.norm())
            .fold(0.0_f64, f64::max)
            .max(1.0);
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (asm.a[(i, j)] - asm.a[(j, i)]).norm() < 1e-10 * scale_a,
                    "A^p2 symmetric at ({i},{j})"
                );
                assert!(
                    (asm.b[(i, j)] - asm.b[(j, i)]).norm() < 1e-10 * scale_b,
                    "B^p2 symmetric at ({i},{j})"
                );
                assert!(
                    (asm.b1[(i, j)] - asm.b1[(j, i)]).norm() < 1e-10 * scale_b,
                    "B_1^p2 symmetric at ({i},{j})"
                );
                assert!(
                    (asm.b[(i, j)] - asm.b1[(i, j)]).norm() < 1e-12 * scale_b,
                    "homogeneous B and B_1 must coincide at ({i},{j})"
                );
            }
        }
    }

    // ── Phase 1.3.1.1 step 5.5 — J3 correctness anchor (DoD-4) ──────────

    /// Solve the dominant-mode β² for a homogeneous air-filled WR-90 mesh at
    /// the given element order, returning β (rad/m). Shared by the J3 anchor.
    fn solve_homogeneous_beta(nx: usize, ny: usize, order: super::super::ElementOrder) -> f64 {
        let a = 22.86e-3;
        let b = 10.16e-3;
        let freq_hz = 10.0e9;
        let mesh = rectangular_mesh(a, b, nx, ny);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = match order {
            super::super::ElementOrder::First => assemble_mixed(&mesh, &eps, &mu, &table),
            super::super::ElementOrder::Second => assemble_mixed_p2(&mesh, &eps, &mu, &table),
        };
        let sol =
            crate::eigensolver::solve_dense_mixed(&asm, freq_hz).expect("homogeneous mixed solve");
        sol.beta_sq.re.max(0.0).sqrt()
    }

    #[test]
    fn p2_homogeneous_wr90_te10_at_least_as_accurate_as_p1() {
        // J3 CORRECTNESS ANCHOR (DoD-4) — the sharp, no-singularity check
        // that the p=2 element matrices + DoF map + solve are all correct.
        // On the homogeneous (air-filled) WR-90 the dominant mode is the
        // analytic TE10 β = √(k₀² − (π/a)²); p=2 must reproduce it AT LEAST
        // AS ACCURATELY AS p=1 on the same mesh. A wrong p=2 element matrix,
        // a wrong edge-DoF orientation sign, or a wrong global assembly fails
        // here — this anchor MUST pass before the high-contrast J4 case.
        let a = 22.86e-3;
        let k0 = std::f64::consts::TAU * 10.0e9 / yee_core::units::C0;
        let kx = std::f64::consts::PI / a;
        let beta_te10 = (k0 * k0 - kx * kx).sqrt();

        // Same 6×6 mesh the WR-90 gate uses.
        let beta_p1 = solve_homogeneous_beta(6, 6, super::super::ElementOrder::First);
        let beta_p2 = solve_homogeneous_beta(6, 6, super::super::ElementOrder::Second);
        let rel_p1 = (beta_p1 - beta_te10).abs() / beta_te10;
        let rel_p2 = (beta_p2 - beta_te10).abs() / beta_te10;
        eprintln!(
            "J3 homogeneous WR-90 TE10 (6×6): analytic β {beta_te10:.6}, \
             p1 β {beta_p1:.6} (rel {rel_p1:.3e}), p2 β {beta_p2:.6} (rel {rel_p2:.3e})"
        );
        // p=2 must be a valid propagating mode within 1 % (the WR-90 gate's
        // own tolerance) AND no worse than p=1 on the identical mesh (DoD-4:
        // "at least as accurately as p=1"). A small numerical slack guards
        // against a tie being flipped by round-off.
        assert!(
            rel_p2 < 0.01,
            "p=2 homogeneous TE10 β {beta_p2} must match analytic {beta_te10} within 1 % \
             (rel {rel_p2:.3e})"
        );
        assert!(
            rel_p2 <= rel_p1 + 1e-9,
            "DoD-4: p=2 TE10 error {rel_p2:.3e} must be ≤ p=1 error {rel_p1:.3e} on the same mesh \
             — a regression here means the p=2 element matrices / DoF sign / assembly are wrong"
        );
    }
}
