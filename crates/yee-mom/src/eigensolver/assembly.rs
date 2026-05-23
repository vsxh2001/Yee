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

/// Result of [`assemble_mixed`]: the **mixed `(E_t, E_z)`**
/// Lee-Sun-Cendes (1991) block generalized eigenproblem
/// `A x = k_c² B x` with `x = [E_t; E_z]`.
///
/// **Eigenvalue convention.** This pencil is assembled in the
/// **cutoff-wavenumber** parameterization `k_c² = k_0² − β²`, identical
/// to [`AssembledTransverse`] / [`super::solve_dense`]. The propagation
/// constant follows from `β² = k_0² − k_c²` at the post-solve mapping;
/// the pencil itself is **frequency-independent** (the staged
/// longitudinal element matrices carry no `k_0²` term — see
/// [`local_a_zz`] / [`local_b_zz`] — which is the decisive evidence that
/// the `k_c²` parameterization is the one they were built for, not the
/// `β²` parameterization quoted in the design-spec §3 prose).
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
    /// Block-mass matrix `B = [[B_tt, B_tz], [B_zt, B_zz]]`, size `n × n`.
    pub b: DMatrix<Complex64>,
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
        // B_ze[i_vert][j_edge] = ∫ (1/μ_r) ∇L_i · N_j (vertex-row /
        // edge-col). The 1/μ_r weight (matching A_zz) is the load-bearing
        // curl-curl coupling; see `local_b_ze`'s docstring for why the
        // ε_r weight was inert.
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
            }
        }

        // --- coupling blocks (edge ↔ vertex), mass side only ---
        // B_ze[lv][le] sits at global (vertex-row n_t+iv, edge-col ie):
        // that is the B_zt block. Its transpose populates B_tz at
        // (edge-row ie, vertex-col n_t+iv). Assembling both halves keeps
        // B symmetric (B_tz = B_ztᵀ), which the Cholesky-symmetrised
        // solve requires.
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
            }
        }
    }

    AssembledMixed {
        a,
        b,
        interior_to_global_edges,
        interior_to_global_verts,
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
        // B_tz = B_ztᵀ explicitly).
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
            }
        }
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
}
