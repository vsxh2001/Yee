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
//! **Mixed-formulation status.** Architecturally the assembly also
//! produces the longitudinal nodal-Lagrange (`E_z`) element matrices
//! [`local_a_zz`] and [`local_b_zz`], plus the edge-node coupling
//! [`local_b_ze`], so a full Bardi-Biró mixed eigensolve can be wired
//! in once non-trivial dielectric stacks need it. The Phase 1.3.1.1
//! step 3 dense fallback ([`super::solve_dense`]) uses only the
//! transverse block, which is exact for the WR-90 TE10 validation
//! case and bit-stable on homogeneously filled PEC cross-sections.

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
/// Used by the mixed Lee-Sun-Cendes formulation. Unused by the
/// transverse-only WR-90 validation path; kept for the architectural
/// contract with the spec and exercised by a smoke test.
#[allow(dead_code, clippy::needless_range_loop)] // wired into the eigensolve in Phase 1.3.1.1 step 5
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
/// Unused by the transverse-only WR-90 validation path; kept for the
/// architectural contract with the spec.
#[allow(dead_code, clippy::needless_range_loop)] // wired into the eigensolve in Phase 1.3.1.1 step 5
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
///   B_ze^e[i,j] = ∫_T ε_r ∇L_i · N_j dA
/// ```
///
/// Linear in `λ` and constant in `∇λ`, so the integral reduces to
/// `(A/3) ε_r ℓ_j σ_j (∇λ_i · ∇λ_{b_j} − ∇λ_i · ∇λ_{a_j})` — wait, more
/// carefully: `N_j = ℓ_j σ_j (λ_{a_j} ∇λ_{b_j} − λ_{b_j} ∇λ_{a_j})` and
/// `∇L_i = ∇λ_i` is constant, so
/// `∫ ∇λ_i · N_j dA = ℓ_j σ_j [(∫ λ_{a_j} dA) (∇λ_i · ∇λ_{b_j})
///                              − (∫ λ_{b_j} dA) (∇λ_i · ∇λ_{a_j})]`
/// and `∫ λ_p dA = A/3`.
///
/// Unused by the transverse-only WR-90 validation path; kept for the
/// architectural contract with the spec.
#[allow(dead_code, clippy::needless_range_loop)] // wired into the eigensolve in Phase 1.3.1.1 step 5
fn local_b_ze(geom: &TriGeom, eps_r: Complex64, signs: [f64; 3]) -> [[Complex64; 3]; 3] {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    let third_area = geom.area / 3.0;
    for i in 0..3 {
        for j in 0..3 {
            let [aj, bj] = LOCAL_EDGES[j];
            let s = third_area * (geom.grad_dot(i, bj) - geom.grad_dot(i, aj));
            out[i][j] = eps_r * Complex::new(geom.edge_len[j] * signs[j] * s, 0.0);
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
pub(crate) struct AssembledTransverse {
    /// Curl-curl stiffness matrix `S[i,j] = ∫ (1/μ_r) (∇×N_i)(∇×N_j) dA`.
    pub s: DMatrix<Complex64>,
    /// ε_r-weighted Nedelec mass matrix `T[i,j] = ∫ ε_r N_i · N_j dA`.
    pub t: DMatrix<Complex64>,
    /// Map from interior-edge DoF index to the global edge index, so
    /// post-solve eigenvector components can be located back on the mesh.
    #[allow(dead_code)] // consumed by the Phase 1.3.1.1 step 5 eigenvector recovery
    pub interior_to_global: Vec<usize>,
    /// Largest ε_r magnitude seen during assembly. Recorded for the
    /// caller-side `β² = k_0² − k_c²` translation (the ε_r weighting in
    /// `T` is folded into `k_c²` since `T` is ε_r-weighted, so the
    /// relation is exactly `β² = k_0² − k_c²` when ε_r is real and
    /// uniform; lossy / heterogeneous ε_r is Phase 1.3.1.2).
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
            }
        }
    }

    AssembledTransverse {
        s,
        t,
        interior_to_global,
        eps_r_max_re,
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
}
