//! Global FEM assembly — first-order Nedelec curl-curl stiffness `K` and
//! vector mass `M` on a tetrahedral mesh, with PEC tangential-`E`-zero
//! Dirichlet elimination by row/column drop.
//!
//! ## Pipeline
//!
//! 1. Walk [`yee_mesh::TetMesh3D::tetrahedra`]. For each tet build the six
//!    canonical local edges per [`crate::element::LOCAL_EDGES`] and resolve
//!    each to a *global* edge index in a deduplicated edge list keyed by
//!    its canonical (lower-endpoint-first) vertex pair.
//! 2. Record the per-tet **local-vs-global orientation sign** per edge:
//!    `+1` if the local edge runs from the lower to the higher global
//!    vertex index, `-1` otherwise. Mirrors the 2-D analog in
//!    `crates/yee-mom/src/eigensolver/mesh.rs`.
//! 3. Call [`crate::element::assemble_tet_element`] on each tet to obtain
//!    the canonical-orientation 6×6 [`crate::element::NedelecTetElement`]
//!    block. Scatter the block into an [`nalgebra_sparse::coo::CooMatrix`]
//!    entry-summing builder, multiplying each row *and* column by its
//!    orientation sign (`s_α s_β` on `(α, β)`).
//! 4. Classify each global edge as boundary or interior. An edge is on the
//!    PEC boundary iff it lies on at least one tet-face that is shared by
//!    fewer than two tets. (A face is identified by its sorted triple of
//!    vertex indices.) Build `interior_edges: Vec<usize>` mapping
//!    interior-DoF index → full-edge index so eigenvectors can be lifted
//!    back to the full edge basis.
//! 5. Convert the COO builder to [`nalgebra_sparse::csr::CsrMatrix`] and
//!    drop the boundary rows and columns from both `K` and `M`.
//!
//! ## Why row/column drop and not penalty
//!
//! Row/column elimination keeps the resulting spectra clean: every
//! returned eigenvalue is a physical resonance of the interior problem.
//! Penalty methods either bias the spectrum by `O(1/penalty)` or push
//! spurious modes into the cluster near `K e = k² M e`'s zero kernel,
//! both of which the Phase 4 plan explicitly rules out (v0 spec §6,
//! design rationale).
//!
//! ## Sign convention
//!
//! Mirrors the 2-D analog in `crates/yee-mom/src/eigensolver/mesh.rs`:
//! the **canonical global orientation** of an edge is from its
//! lower-indexed endpoint to its higher-indexed endpoint. The element
//! layer ([`crate::element::assemble_tet_element`]) emits its 6×6 block
//! in **canonical local-edge orientation** — local edges are listed in
//! [`crate::element::LOCAL_EDGES`] which itself is lower-endpoint-first.
//! Because the tet vertex tuple `[v0, v1, v2, v3]` is not necessarily
//! sorted (it merely has positive signed volume), a local edge `(i, j)`
//! with `i < j` *in local indices* may still run against the global
//! orientation when `tet[i] > tet[j]` in *global* indices. The scatter
//! step multiplies row `α` and column `β` by their respective signs to
//! reconcile the two orientations.
//!
//! ## Boundary-edge classifier
//!
//! Per the Phase 4 plan T2/T4 escape hatch: build a face-incidence map
//! `HashMap<[usize; 3], usize>` (sorted vertex triplet → tet-count). An
//! edge is boundary iff at least one of its incident **faces** has
//! tet-count `< 2`. This is the same classifier the plan documents at
//! `yee-mesh::TetMesh3D::boundary_edges()`; that method is **not yet
//! implemented** at the worktree base SHA (Phase 4 T2 shipped only the
//! constructor and signed-volume / centroid accessors), so the
//! classifier lives here for now. The same definition will move to
//! `yee-mesh` when T2's edges()/boundary_edges() are finished —
//! flagged as an out-of-lane finding in the Track KKKKKK report.

use std::collections::{HashMap, HashSet};

use nalgebra_sparse::{coo::CooMatrix, csr::CsrMatrix};
use num_complex::Complex64;
use yee_mesh::TetMesh3D;

use crate::element::{LOCAL_EDGES, assemble_tet_element, assemble_tet_element_complex};
use crate::material::MaterialDatabase;

/// Output of [`FemEigenAssembly::assemble`] — the interior-DoF-reduced
/// curl-curl stiffness `K`, vector mass `M`, and the lift map from
/// interior-DoF index to full-edge index for eigenvector recovery.
///
/// `K` and `M` are real-symmetric and have identical sparsity pattern.
/// They are returned in [`nalgebra_sparse::csr::CsrMatrix`] form so the
/// downstream sparse generalized eigensolve (Phase 4 T5) can consume
/// them with no further format conversion.
#[derive(Debug, Clone)]
pub struct AssembledMatrices {
    /// Curl-curl stiffness `K[i,j] = Σ_e (1/μ_r^e) ∫_T_e (∇×N_i)·(∇×N_j) dV`
    /// reduced to interior edges by PEC Dirichlet row/column elimination.
    pub k: CsrMatrix<f64>,
    /// Vector mass `M[i,j] = Σ_e ε_r^e ∫_T_e N_i · N_j dV`, reduced to
    /// the same interior-DoF basis as [`Self::k`].
    pub m: CsrMatrix<f64>,
    /// Lift map: `interior_edges[i]` is the full-edge index of the
    /// interior DoF at row/column `i`. Used to scatter eigenvector
    /// components back to the full edge basis at solve time.
    pub interior_edges: Vec<usize>,
}

/// Output of [`FemEigenAssembly::assemble_complex`] — the
/// interior-DoF-reduced **complex** curl-curl stiffness `K(ω)`, vector
/// mass `M(ω)`, and the lift map from interior-DoF index to full-edge
/// index for eigenvector recovery.
///
/// `K` and `M` are complex symmetric (not Hermitian — see ADR-0039 §6
/// and `solve.rs` module docs) with identical sparsity pattern. They
/// are returned in [`nalgebra_sparse::csr::CsrMatrix`] of [`Complex64`]
/// so the downstream Phase 4.fem.eig.1 complex sparse eigensolve
/// ([`crate::solve::ComplexInverseIterEigen`]) can consume them with
/// no further format conversion.
///
/// For inputs whose per-tet `ε(ω)`, `μ(ω)` are purely real, the
/// imaginary parts of every entry are bit-for-bit zero and `.re` agrees
/// with [`AssembledMatrices`] from the real path on the same mesh; this
/// is the backward-compatibility invariant exercised by the
/// `dispersive_solve` integration tests.
#[derive(Debug, Clone)]
pub struct AssembledMatricesComplex {
    /// Curl-curl stiffness
    /// `K[i,j] = Σ_e (1/μ_e(ω)) ∫_T_e (∇×N_i)·(∇×N_j) dV`
    /// reduced to interior edges by PEC Dirichlet row/column
    /// elimination.
    pub k: CsrMatrix<Complex64>,
    /// Vector mass `M[i,j] = Σ_e ε_e(ω) ∫_T_e N_i · N_j dV`, reduced to
    /// the same interior-DoF basis as [`Self::k`].
    pub m: CsrMatrix<Complex64>,
    /// Lift map: `interior_edges[i]` is the full-edge index of the
    /// interior DoF at row/column `i`. Used to scatter eigenvector
    /// components back to the full edge basis at solve time. Identical
    /// to [`AssembledMatrices::interior_edges`] for the same mesh —
    /// orientation and boundary classification are geometric, not
    /// material-dependent.
    pub interior_edges: Vec<usize>,
}

/// Per-tet material-aware global FEM assembler.
///
/// The assembler is parameterised on a borrowed [`yee_mesh::TetMesh3D`]
/// and two per-tet material vectors (`eps_r`, `mu_r`). For uniform
/// free-space the [`Self::new_free_space`] convenience constructor seeds
/// both with `1.0` for every tet; callers with non-trivial materials
/// supply explicit length-`N_tets` vectors via [`Self::new`].
///
/// The struct is intentionally cheap to build — material vectors are
/// owned `Vec<f64>` so the assembler does not need to keep the caller's
/// `MaterialTag` lookup table alive — and idempotent: repeated calls to
/// [`Self::assemble`] produce bit-identical output for the same inputs.
#[derive(Debug, Clone)]
pub struct FemEigenAssembly<'m> {
    /// Borrowed mesh; the assembler does not mutate it.
    mesh: &'m TetMesh3D,
    /// Length-`N_tets` relative permittivity, one entry per tet.
    eps_r: Vec<f64>,
    /// Length-`N_tets` relative permeability, one entry per tet.
    mu_r: Vec<f64>,
}

impl<'m> FemEigenAssembly<'m> {
    /// Build an assembler with explicit per-tet `ε_r` and `μ_r`.
    ///
    /// Both vectors must have length equal to `mesh.tetrahedra.len()`,
    /// otherwise the assembler returns [`yee_core::Error::Invalid`] on
    /// [`Self::assemble`]. We validate at construction so the error
    /// surfaces as close to the caller as possible.
    ///
    /// # Errors
    ///
    /// Returns [`yee_core::Error::Invalid`] if either material vector
    /// length does not match `mesh.tetrahedra.len()`.
    pub fn new(
        mesh: &'m TetMesh3D,
        eps_r: Vec<f64>,
        mu_r: Vec<f64>,
    ) -> Result<Self, yee_core::Error> {
        let n_tets = mesh.tetrahedra.len();
        if eps_r.len() != n_tets {
            return Err(yee_core::Error::Invalid(format!(
                "eps_r length {} does not match tet count {n_tets}",
                eps_r.len()
            )));
        }
        if mu_r.len() != n_tets {
            return Err(yee_core::Error::Invalid(format!(
                "mu_r length {} does not match tet count {n_tets}",
                mu_r.len()
            )));
        }
        Ok(Self { mesh, eps_r, mu_r })
    }

    /// Convenience constructor for an air-filled (free-space) cavity:
    /// every tet has `ε_r = μ_r = 1.0`. Matches the fem-eig-001
    /// validation gate (Pozar §6.3 air cavity).
    pub fn new_free_space(mesh: &'m TetMesh3D) -> Self {
        let n_tets = mesh.tetrahedra.len();
        Self {
            mesh,
            eps_r: vec![1.0; n_tets],
            mu_r: vec![1.0; n_tets],
        }
    }

    /// Assemble the global sparse `K`, `M` and apply PEC Dirichlet
    /// elimination.
    ///
    /// See module docs for the pipeline. The output is an
    /// [`AssembledMatrices`] holding interior-DoF-reduced `K`, `M` (both
    /// in CSR form, same sparsity, real-symmetric) and the lift map for
    /// eigenvector recovery.
    ///
    /// # Errors
    ///
    /// Returns [`yee_core::Error::Invalid`] if the mesh is empty (no
    /// tetrahedra) or — in principle — if any tet has degenerate
    /// geometry. The latter is rejected at
    /// [`yee_mesh::TetMesh3D::new`], so a well-formed mesh will not
    /// trigger it here.
    pub fn assemble(&self) -> Result<AssembledMatrices, yee_core::Error> {
        if self.mesh.tetrahedra.is_empty() {
            return Err(yee_core::Error::Invalid(
                "cannot assemble FEM matrices on a mesh with zero tetrahedra".to_string(),
            ));
        }

        // ---- Step 1+2: build the global edge table and per-tet
        // local-to-global edge map with orientation signs --------------
        let table = TetEdgeTable::build(self.mesh);
        let n_edges = table.edges.len();

        // ---- Step 4 (pre-scatter): classify each global edge as
        // boundary or interior, and build the lift map ----------------
        let interior_edges: Vec<usize> = (0..n_edges).filter(|&e| !table.is_boundary[e]).collect();
        let n_interior = interior_edges.len();
        // Inverse map: global-edge-index -> interior DoF index (or None
        // if the edge is on the boundary and thus eliminated).
        let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; n_edges];
        for (dof, &gid) in interior_edges.iter().enumerate() {
            interior_dof_of_edge[gid] = Some(dof);
        }

        // ---- Step 3+5: scatter signed 6×6 local blocks into COO,
        // skipping boundary DoFs as we go (cheaper than building the
        // full-edge matrix and then slicing) --------------------------
        let mut k_coo: CooMatrix<f64> = CooMatrix::new(n_interior, n_interior);
        let mut m_coo: CooMatrix<f64> = CooMatrix::new(n_interior, n_interior);

        for (tet_idx, conn) in table.tet_edges.iter().enumerate() {
            let tet = &self.mesh.tetrahedra[tet_idx];
            let vertices = [
                self.mesh.vertices[tet[0]],
                self.mesh.vertices[tet[1]],
                self.mesh.vertices[tet[2]],
                self.mesh.vertices[tet[3]],
            ];
            let eps_r = self.eps_r[tet_idx];
            let mu_r = self.mu_r[tet_idx];
            let elem = assemble_tet_element(vertices, eps_r, mu_r);

            for alpha in 0..6 {
                let gi = conn.global_edge[alpha];
                let Some(ii) = interior_dof_of_edge[gi] else {
                    continue;
                };
                let sa = conn.sign[alpha];
                for beta in 0..6 {
                    let gj = conn.global_edge[beta];
                    let Some(jj) = interior_dof_of_edge[gj] else {
                        continue;
                    };
                    let sb = conn.sign[beta];
                    let signed = sa * sb;
                    // CooMatrix's `push` is entry-summing on conversion to
                    // CSR, so multiple contributions from neighbouring
                    // tets to the same `(ii, jj)` are accumulated.
                    k_coo.push(ii, jj, signed * elem.k_local[(alpha, beta)]);
                    m_coo.push(ii, jj, signed * elem.m_local[(alpha, beta)]);
                }
            }
        }

        let k = CsrMatrix::from(&k_coo);
        let m = CsrMatrix::from(&m_coo);

        Ok(AssembledMatrices {
            k,
            m,
            interior_edges,
        })
    }

    /// Assemble the global **complex** sparse `K(ω)`, `M(ω)` at angular
    /// frequency `omega` (rad/s) and apply PEC Dirichlet elimination.
    ///
    /// Walks the same per-tet pipeline as [`Self::assemble`] but looks
    /// up `(ε(ω), μ(ω))` per-tet from the supplied [`MaterialDatabase`]
    /// keyed by [`yee_mesh::TetMesh3D::tetrahedron_material`], and
    /// calls [`crate::element::assemble_tet_element_complex`] instead
    /// of the real path. The scatter, orientation signs, boundary-edge
    /// classifier, and PEC row/column elimination are geometric (not
    /// material-dependent) and are reused verbatim from the real path.
    ///
    /// For materials whose `ε(ω)`, `μ(ω)` are purely real at the
    /// supplied `omega` (e.g. [`crate::material::Material::default`]
    /// free-space response) the returned blocks have
    /// `Im(K[i,j]) = Im(M[i,j]) = 0` bit-for-bit and `.re` agrees with
    /// the real [`Self::assemble`] output on the same mesh. The
    /// `assemble_complex_at_real_eps_matches_real_assemble` integration
    /// test pins this invariant — it is the load-bearing
    /// backward-compatibility check per ADR-0039 §4.
    ///
    /// # Arguments
    ///
    /// * `omega` — trial angular frequency (rad/s). Real-valued by
    ///   design: the Phase 4.fem.eig.1 Newton tracker (plan step D5)
    ///   sweeps the linearised eigenproblem along the **real** ω axis
    ///   and lets the complex eigenvalue `k²` carry the imaginary part.
    ///   A complex trial frequency would force `assemble_tet_element`
    ///   to evaluate `ε`, `μ` at complex argument, which the
    ///   single-pole [`crate::material::Material::eps_at`] surface does
    ///   not support in v1.
    /// * `db` — [`MaterialDatabase`] keyed by [`yee_mesh::MaterialTag`].
    ///   Unregistered tags fall back to free-space `ε = μ = 1`.
    ///
    /// # Errors
    ///
    /// Returns [`yee_core::Error::Invalid`] if the mesh is empty (no
    /// tetrahedra). The per-tet material lookup itself does not fail —
    /// unregistered tags return free-space `ε = μ = 1` per
    /// [`MaterialDatabase::eps_at`] / [`MaterialDatabase::mu_at`].
    pub fn assemble_complex(
        &self,
        omega: f64,
        db: &MaterialDatabase,
    ) -> Result<AssembledMatricesComplex, yee_core::Error> {
        if self.mesh.tetrahedra.is_empty() {
            return Err(yee_core::Error::Invalid(
                "cannot assemble FEM matrices on a mesh with zero tetrahedra".to_string(),
            ));
        }

        // ---- Step 1+2: build the global edge table and per-tet
        // local-to-global edge map with orientation signs.  This
        // mirrors [`Self::assemble`] line-for-line because the
        // geometry / topology is material-independent.
        let table = TetEdgeTable::build(self.mesh);
        let n_edges = table.edges.len();

        // ---- Step 4 (pre-scatter): classify each global edge as
        // boundary or interior, and build the lift map.
        let interior_edges: Vec<usize> = (0..n_edges).filter(|&e| !table.is_boundary[e]).collect();
        let n_interior = interior_edges.len();
        let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; n_edges];
        for (dof, &gid) in interior_edges.iter().enumerate() {
            interior_dof_of_edge[gid] = Some(dof);
        }

        // ---- Step 3+5: scatter signed 6×6 complex local blocks into
        // COO, skipping boundary DoFs as we go.  The orientation sign
        // `sa * sb ∈ {-1, +1}` is real and applied as a real scalar to
        // the complex block, matching the real path bit-for-bit on the
        // sign convention.
        let mut k_coo: CooMatrix<Complex64> = CooMatrix::new(n_interior, n_interior);
        let mut m_coo: CooMatrix<Complex64> = CooMatrix::new(n_interior, n_interior);

        for (tet_idx, conn) in table.tet_edges.iter().enumerate() {
            let tet = &self.mesh.tetrahedra[tet_idx];
            let vertices = [
                self.mesh.vertices[tet[0]],
                self.mesh.vertices[tet[1]],
                self.mesh.vertices[tet[2]],
                self.mesh.vertices[tet[3]],
            ];
            let tag = self.mesh.tetrahedron_material[tet_idx];
            let eps_omega = db.eps_at(tag, omega);
            let mu_omega = db.mu_at(tag, omega);
            let elem = assemble_tet_element_complex(vertices, eps_omega, mu_omega);

            for alpha in 0..6 {
                let gi = conn.global_edge[alpha];
                let Some(ii) = interior_dof_of_edge[gi] else {
                    continue;
                };
                let sa = conn.sign[alpha];
                for beta in 0..6 {
                    let gj = conn.global_edge[beta];
                    let Some(jj) = interior_dof_of_edge[gj] else {
                        continue;
                    };
                    let sb = conn.sign[beta];
                    let signed = Complex64::new(sa * sb, 0.0);
                    k_coo.push(ii, jj, signed * elem.k_local[(alpha, beta)]);
                    m_coo.push(ii, jj, signed * elem.m_local[(alpha, beta)]);
                }
            }
        }

        let k = CsrMatrix::from(&k_coo);
        let m = CsrMatrix::from(&m_coo);

        Ok(AssembledMatricesComplex {
            k,
            m,
            interior_edges,
        })
    }

    /// Assemble the global **complex** sparse `K(ω)`, `M(ω)` at angular
    /// frequency `omega` using an **explicit PEC-edge set** for the
    /// Dirichlet elimination, rather than the geometric boundary-face
    /// classifier.
    ///
    /// This is the entry point consumed by Phase 4.fem.eig.2
    /// [`crate::OpenBoundarySolver`] when the caller has tagged some
    /// exterior faces ABC or wave-port: those faces' edges must remain
    /// in the interior-DoF basis (where the boundary-term scatter then
    /// adds the corresponding face block / RHS contribution), while
    /// edges on caller-tagged PEC faces are eliminated by row/column
    /// drop.
    ///
    /// The function returns the **edge table's `is_boundary` array**
    /// for the open-boundary caller — i.e. the geometric boundary-edge
    /// set computed from the mesh's face-incidence map — alongside the
    /// PEC-reduced complex `K(ω)`, `M(ω)`, and the interior-DoF lift
    /// map. The caller's [`crate::OpenBoundarySolver`] uses the
    /// (`is_boundary` ∩ ¬`pec_edges`) set to identify ABC + wave-port
    /// edges and scatter the corresponding face-block / RHS
    /// contributions on top.
    ///
    /// When `pec_edges` is the same as the geometric boundary-edge set
    /// (i.e. all-PEC), this function is bit-for-bit equivalent to
    /// [`Self::assemble_complex`].
    ///
    /// # Arguments
    ///
    /// * `omega` — trial angular frequency (rad/s).
    /// * `db` — [`MaterialDatabase`] keyed by [`yee_mesh::MaterialTag`].
    /// * `pec_edges` — set of global edge indices to eliminate by
    ///   row/column drop. Edges not in this set remain as interior
    ///   DoFs even if they geometrically lie on the mesh exterior.
    ///
    /// # Errors
    ///
    /// Returns [`yee_core::Error::Invalid`] on an empty mesh, same as
    /// [`Self::assemble_complex`].
    pub fn assemble_complex_with_pec_edges(
        &self,
        omega: f64,
        db: &MaterialDatabase,
        pec_edges: &HashSet<usize>,
    ) -> Result<AssembledMatricesComplex, yee_core::Error> {
        if self.mesh.tetrahedra.is_empty() {
            return Err(yee_core::Error::Invalid(
                "cannot assemble FEM matrices on a mesh with zero tetrahedra".to_string(),
            ));
        }

        // ---- Step 1+2: build the global edge table identical to
        // assemble_complex (geometry is material-independent).
        let table = TetEdgeTable::build(self.mesh);
        let n_edges = table.edges.len();

        // ---- Step 4 (pre-scatter): classify each global edge as
        // interior or PEC by the caller-supplied pec_edges set. Edges
        // not in pec_edges remain as interior DoFs (including
        // exterior-but-non-PEC faces — ABC and wave-port — whose
        // boundary-term contributions are scattered by the caller).
        let interior_edges: Vec<usize> = (0..n_edges).filter(|e| !pec_edges.contains(e)).collect();
        let n_interior = interior_edges.len();
        let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; n_edges];
        for (dof, &gid) in interior_edges.iter().enumerate() {
            interior_dof_of_edge[gid] = Some(dof);
        }

        // ---- Step 3+5: scatter signed 6×6 complex local blocks into
        // COO, skipping PEC DoFs as we go.
        let mut k_coo: CooMatrix<Complex64> = CooMatrix::new(n_interior, n_interior);
        let mut m_coo: CooMatrix<Complex64> = CooMatrix::new(n_interior, n_interior);

        for (tet_idx, conn) in table.tet_edges.iter().enumerate() {
            let tet = &self.mesh.tetrahedra[tet_idx];
            let vertices = [
                self.mesh.vertices[tet[0]],
                self.mesh.vertices[tet[1]],
                self.mesh.vertices[tet[2]],
                self.mesh.vertices[tet[3]],
            ];
            let tag = self.mesh.tetrahedron_material[tet_idx];
            let eps_omega = db.eps_at(tag, omega);
            let mu_omega = db.mu_at(tag, omega);
            let elem = assemble_tet_element_complex(vertices, eps_omega, mu_omega);

            for alpha in 0..6 {
                let gi = conn.global_edge[alpha];
                let Some(ii) = interior_dof_of_edge[gi] else {
                    continue;
                };
                let sa = conn.sign[alpha];
                for beta in 0..6 {
                    let gj = conn.global_edge[beta];
                    let Some(jj) = interior_dof_of_edge[gj] else {
                        continue;
                    };
                    let sb = conn.sign[beta];
                    let signed = Complex64::new(sa * sb, 0.0);
                    k_coo.push(ii, jj, signed * elem.k_local[(alpha, beta)]);
                    m_coo.push(ii, jj, signed * elem.m_local[(alpha, beta)]);
                }
            }
        }

        let k = CsrMatrix::from(&k_coo);
        let m = CsrMatrix::from(&m_coo);

        Ok(AssembledMatricesComplex {
            k,
            m,
            interior_edges,
        })
    }

    /// Borrow the mesh this assembler was built against. Useful for
    /// downstream solve / post-processing steps that need to look up
    /// vertex coordinates by edge index.
    pub fn mesh(&self) -> &'m TetMesh3D {
        self.mesh
    }

    /// Read-only borrow of the per-tet `ε_r` array. Length matches
    /// `mesh.tetrahedra.len()`.
    pub fn eps_r(&self) -> &[f64] {
        &self.eps_r
    }

    /// Read-only borrow of the per-tet `μ_r` array. Length matches
    /// `mesh.tetrahedra.len()`.
    pub fn mu_r(&self) -> &[f64] {
        &self.mu_r
    }
}

/// Canonical edge key — lower-endpoint-first vertex pair. Mirrors
/// [`yee_mom::eigensolver::mesh::EdgeKey`] one dimension up.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EdgeKey {
    from: usize,
    to: usize,
}

impl EdgeKey {
    fn new(a: usize, b: usize) -> Self {
        if a < b {
            Self { from: a, to: b }
        } else {
            Self { from: b, to: a }
        }
    }
}

/// Per-tet local→global edge map and orientation signs.
#[derive(Debug, Clone, Copy)]
struct TetEdgeConnectivity {
    /// `global_edge[α]` is the global edge index for local edge `α` of
    /// the tet (per [`crate::element::LOCAL_EDGES`]).
    global_edge: [usize; 6],
    /// `sign[α]` is `+1` if the local-edge direction matches the
    /// canonical global orientation, `-1` otherwise.
    sign: [f64; 6],
}

/// Edge table for a [`yee_mesh::TetMesh3D`]. Holds the deduplicated
/// edge list, the per-edge boundary flag, and the per-tet local→global
/// edge connectivity with orientation signs.
#[derive(Debug, Clone)]
struct TetEdgeTable {
    /// All distinct edges in canonical (`from < to`) orientation.
    edges: Vec<EdgeKey>,
    /// `is_boundary[e]` is `true` iff edge `e` lies on at least one
    /// boundary face (a face shared by exactly one tet).
    is_boundary: Vec<bool>,
    /// Per-tet local→global edge map + orientation signs.
    tet_edges: Vec<TetEdgeConnectivity>,
}

impl TetEdgeTable {
    /// Build the edge table in `O(n_tets)` with a `HashMap` for
    /// edge-key dedup and a parallel `HashMap` for face-incidence
    /// counting (boundary-edge classifier).
    fn build(mesh: &TetMesh3D) -> Self {
        // ---- Edge dedup --------------------------------------------
        let mut edge_map: HashMap<EdgeKey, usize> = HashMap::new();
        let mut edges: Vec<EdgeKey> = Vec::new();
        let mut tet_edges: Vec<TetEdgeConnectivity> = Vec::with_capacity(mesh.tetrahedra.len());

        for tet in &mesh.tetrahedra {
            let mut global_edge = [0usize; 6];
            let mut sign = [0.0f64; 6];
            for (alpha, &(li, lj)) in LOCAL_EDGES.iter().enumerate() {
                let a = tet[li];
                let b = tet[lj];
                let key = EdgeKey::new(a, b);
                let idx = *edge_map.entry(key).or_insert_with(|| {
                    edges.push(key);
                    edges.len() - 1
                });
                global_edge[alpha] = idx;
                // Local direction is `a -> b` (local `li -> lj` with
                // `li < lj` since LOCAL_EDGES is lower-first). Canonical
                // global direction is `from -> to` with `from < to`.
                // Sign is `+1` iff `a < b` in global vertex indices.
                sign[alpha] = if a < b { 1.0 } else { -1.0 };
            }
            tet_edges.push(TetEdgeConnectivity { global_edge, sign });
        }

        // ---- Boundary-edge classifier (per the Phase 4 plan T2
        // escape-hatch and T4 brief): build a face-incidence map keyed
        // by sorted vertex triplet. A face is "boundary" iff it appears
        // in exactly one tet. An edge is boundary iff *any* of its
        // incident faces is boundary. -------------------------------
        let mut face_count: HashMap<[usize; 3], usize> = HashMap::new();
        // Each tet has 4 faces: the three vertices opposite to each
        // local vertex `i ∈ 0..4`. Build sorted-triplet keys so faces
        // shared by two tets collide regardless of vertex permutation.
        const TET_FACES: [[usize; 3]; 4] = [
            [1, 2, 3], // face opposite v0
            [0, 2, 3], // face opposite v1
            [0, 1, 3], // face opposite v2
            [0, 1, 2], // face opposite v3
        ];
        for tet in &mesh.tetrahedra {
            for &[a, b, c] in &TET_FACES {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_count.entry(key).or_insert(0) += 1;
            }
        }

        // An edge `(u, v)` is boundary iff some boundary face contains
        // both `u` and `v`. Equivalently, walk every face and, if the
        // face has count `< 2`, mark its three edges (each pair of
        // its three vertices) as boundary.
        let mut is_boundary: Vec<bool> = vec![false; edges.len()];
        for (face, &count) in &face_count {
            if count >= 2 {
                continue;
            }
            // Three undirected edges of this boundary face.
            for &(la, lb) in &[(0usize, 1usize), (0, 2), (1, 2)] {
                let key = EdgeKey::new(face[la], face[lb]);
                if let Some(&idx) = edge_map.get(&key) {
                    is_boundary[idx] = true;
                }
            }
        }

        Self {
            edges,
            is_boundary,
            tet_edges,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// Reference unit tet `[(0,0,0), (1,0,0), (0,1,0), (0,0,1)]`.
    fn single_tet() -> TetMesh3D {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 3]];
        TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
    }

    /// Two tets sharing a triangular face (the face (0, 1, 2) at z=0).
    /// Vertex 3 is above the plane, vertex 4 is below.
    fn two_tets_shared_face() -> TetMesh3D {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.0, 0.0, -1.0),
        ];
        let tetrahedra = vec![[0, 1, 2, 3], [0, 1, 2, 4]];
        TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
    }

    #[test]
    fn edge_table_single_tet_has_six_boundary_edges() {
        let mesh = single_tet();
        let table = TetEdgeTable::build(&mesh);
        assert_eq!(table.edges.len(), 6, "a single tet has 6 distinct edges");
        // Every edge of a free tet is on a boundary face (each of the 4
        // faces is shared by 1 tet only), so every edge is boundary.
        assert!(
            table.is_boundary.iter().all(|&b| b),
            "every edge of a free tet must be a boundary edge"
        );
    }

    #[test]
    fn edge_table_two_tets_shared_face_dedup_to_nine_edges() {
        let mesh = two_tets_shared_face();
        let table = TetEdgeTable::build(&mesh);
        // The two tets share the face (0, 1, 2). That face contributes 3
        // edges to the global list, and the remaining 6 unique edges
        // come from the two "tip" connections: (0,3), (1,3), (2,3) and
        // (0,4), (1,4), (2,4). Total = 3 + 6 = 9.
        assert_eq!(
            table.edges.len(),
            9,
            "two tets sharing one triangular face have 9 distinct edges"
        );
    }
}
