//! Phase 4.fem.eig.2 step E3 — open-boundary FEM driven solver.
//!
//! This module wraps the Phase 4.fem.eig.0/1 closed-cavity assembly path
//! ([`crate::FemEigenAssembly::assemble_complex`]) with **face-kind
//! classification** and **per-face boundary-term scatter** so the FEM
//! stack can solve the driven open-boundary system
//!
//! ```text
//!     ( K(ω) − k₀² M(ω)  +  Σ_ABC  j k₀ B_ABC  +  Σ_port j β B_port ) e
//!         =  Σ_port b_port,
//! ```
//!
//! at a single angular frequency ω, returning the complex edge-DoF
//! solution vector on the interior-edge basis (PEC Dirichlet eliminated).
//!
//! ## Face classification and PEC precedence
//!
//! On construction, [`OpenBoundarySolver::new`] iterates over every
//! exterior face of the supplied [`yee_mesh::TetMesh3D`] and assigns one
//! of three [`FaceKind`] tags:
//!
//! - [`FaceKind::Pec`] — tangential-`E`-zero Dirichlet, eliminated from
//!   the global system by row/column drop (the default for any face the
//!   caller does not tag).
//! - [`FaceKind::Abc`] — 1st-order Engquist–Majda absorbing boundary;
//!   contributes a per-face `+ j k₀ B_ABC` block to the global complex
//!   stiffness matrix (Phase 4.fem.eig.2 step E1).
//! - [`FaceKind::WavePort`] — modal wave-port; contributes a per-face
//!   `+ j β B_port` block to the global stiffness matrix and a per-face
//!   `+ 2 j β · ∫ N_i · e_t dS` contribution to the global RHS vector
//!   (Phase 4.fem.eig.2 step E2).
//!
//! **PEC precedence (spec §10 risk #5):** an edge that lies on the
//! intersection of a PEC face and a wave-port face is forced PEC. The
//! PEC sidewall must enforce tangential-`E` zero on shared edges; the
//! wave-port modal source on those edges would conflict with the PEC
//! tangential-zero condition. The classifier therefore collects the
//! union of edges-on-PEC-faces *first*, then the wave-port face scatter
//! skips any edge in that PEC Dirichlet set during the modal-RHS and
//! face-block scatter steps. Tested by
//! `pec_precedence_over_waveport_at_shared_edges` in
//! `crates/yee-fem/tests/open_boundary_assembly.rs`.
//!
//! ## Sign / orientation convention
//!
//! Face vertices are ordered so the outward normal of the parent tet
//! points *away* from the tet interior. For the exterior faces of a
//! [`yee_mesh::TetMesh3D`] this is enforced by the face classifier
//! `ExteriorFaceTable::build`: for each boundary face, the three
//! face vertices are emitted in CCW order as seen from outside the mesh
//! (i.e. `(v_a, v_b, v_c) × (v_b − v_a) · n̂_out > 0`).
//!
//! Edge tangents on the face follow the canonical CCW traversal:
//! `t_i = v_{(i+1) mod 3} − v_i` for `i ∈ {0, 1, 2}`. Each face-local
//! edge is matched to a global edge index via the same
//! lower-endpoint-first canonical orientation used by
//! [`crate::assembly::FemEigenAssembly::assemble_complex`]; the
//! orientation sign `s_i ∈ {-1, +1}` is applied at scatter time to row
//! AND column (or row alone for the RHS).
//!
//! ## Pipeline at a single ω
//!
//! 1. Construct the per-tet `(ε(ω), μ(ω))` via the stored
//!    [`crate::MaterialDatabase`] and call
//!    [`crate::FemEigenAssembly::assemble_complex`] to obtain the
//!    PEC-reduced complex sparse `K(ω)` and `M(ω)` along with the
//!    interior-edge lift map.
//! 2. Form the closed-cavity driven core
//!    `A(ω) = K(ω) − k₀² M(ω)` with `k₀ = ω / c`.
//! 3. For each ABC face: call
//!    [`crate::element::assemble_abc_face_block`] with the face's outward
//!    normal, `k₀`, and free-space `μ_r = 1`; scatter the 3×3 block into
//!    `A(ω)` at the interior-DoF indices of the three face edges,
//!    applying the per-edge orientation sign.
//! 4. For each wave-port face: compute the modal `β_mode(ω)` via the
//!    caller-supplied [`PortDefinition::beta_mode`] closure, evaluate
//!    the modal `e_t(centroid)` via [`PortDefinition::modal_e_t`], then
//!    scatter both
//!    [`crate::element::assemble_port_face_block`] and
//!    [`crate::element::assemble_port_modal_rhs`] contributions.
//! 5. Solve `A(ω) e = b` once via `faer::sparse::Lu<usize, Complex64>`
//!    (the same surface Phase 4.fem.eig.1 already exercises).
//!
//! ## API placeholder vs spec §6
//!
//! The spec ships a `WavePortFace` type carrying a
//! `NumericalCrossSection` from `yee-mom`'s Phase 1.3.1.1 cross-section
//! eigensolver. The E3 lane is restricted to `crates/yee-fem/**`, so we
//! ship a closure-based [`PortDefinition`] surface instead — the caller
//! supplies `beta_mode(ω)` and `modal_e_t(x)` as Rust closures, which
//! either evaluate analytic profiles (TE_{10} `ŷ sin(π x / a)` for
//! WR-90) or wrap a `NumericalCrossSection` accessor on the consumer
//! side. The spec's `NumericalCrossSection` integration lands in a
//! follow-up cross-lane PR (Phase 4.fem.eig.2.0.2 per spec §13).
//!
//! ## References
//!
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
//!   §4 (theory), §6 (API surface).
//! * `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md`
//!   step E3.
//! * `docs/src/decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md`.

use std::collections::{HashMap, HashSet};

use faer::linalg::solvers::SolveCore;
use faer::sparse::{SparseColMat, Triplet, linalg::solvers::Lu};
use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::Error;
use yee_core::units::C0;
use yee_mesh::TetMesh3D;

use crate::assembly::FemEigenAssembly;
use crate::element::{
    LOCAL_EDGES, assemble_abc_face_block, assemble_port_face_block, assemble_port_modal_rhs,
};
use crate::material::MaterialDatabase;

/// Identifier of a wave-port descriptor in
/// [`OpenBoundarySolver::ports`]. Stable across calls to
/// [`OpenBoundarySolver::solve_at_frequency`].
pub type PortId = usize;

/// Tag classifying an exterior mesh face for the open-boundary FEM
/// pipeline.
///
/// See module-level docs for the assembly contribution each kind emits.
/// The default for any unannotated exterior face is [`FaceKind::Pec`];
/// PEC tangential-`E`-zero is enforced by Dirichlet row/column drop
/// during [`crate::FemEigenAssembly::assemble_complex`], so a fully-PEC
/// caller produces a complex matrix that matches the Phase 4.fem.eig.0
/// closed-cavity assembly bit-for-bit.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FaceKind {
    /// Perfect electric conductor — tangential-`E`-zero Dirichlet on
    /// every edge of the face. Edges on PEC faces are eliminated from
    /// the global system by row/column drop. Takes precedence over
    /// [`FaceKind::WavePort`] on shared edges (spec §10 risk #5).
    Pec,
    /// 1st-order Engquist–Majda absorbing boundary. Contributes a
    /// per-face `+ j k₀ B_ABC` block to the global complex stiffness
    /// matrix via [`crate::element::assemble_abc_face_block`].
    Abc,
    /// Modal wave-port with the descriptor at index `PortId` in
    /// [`OpenBoundarySolver::ports`]. Contributes a per-face
    /// `+ j β B_port` block to the stiffness matrix and a per-face
    /// modal-current contribution to the RHS vector.
    WavePort(PortId),
}

/// Caller-supplied wave-port descriptor.
///
/// One [`PortDefinition`] per physical port. The two closures together
/// specify the port's modal behaviour:
///
/// - [`Self::beta_mode`] returns the modal propagation constant `β(ω)`
///   at angular frequency `ω`. For a TE_{10} mode on rectangular
///   waveguide of broad-wall `a`, this is
///   `β(ω) = sqrt((ω/c)² − (π/a)²)`. The caller is responsible for
///   handling the below-cutoff regime (returning `0` or a small real
///   value); the assembly path passes the returned value through
///   verbatim to [`crate::element::assemble_port_face_block`] and
///   [`crate::element::assemble_port_modal_rhs`].
/// - [`Self::modal_e_t`] returns the tangential incident-mode E-field
///   at a world-space point on the port face, already scaled by the
///   caller's incident amplitude `a_inc`. Phase 4.fem.eig.2 v0 samples
///   this at the face centroid; per-Gauss-point sampling is a
///   Phase 4.fem.eig.2.0.1 refinement.
///
/// The wave-port face is identified by the face-classification list
/// [`OpenBoundarySolver::face_kinds`] — every face tagged
/// `FaceKind::WavePort(p)` for the same `p` contributes to port `p`'s
/// stiffness + RHS scatter.
pub struct PortDefinition {
    /// Modal propagation constant at angular frequency `ω` (rad/s).
    /// Returns `β(ω)` in rad/m. Real-valued; the caller is responsible
    /// for clipping below-cutoff (`β = 0`) if applicable.
    pub beta_mode: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    /// Tangential modal `E_t(x)` already scaled by the incident
    /// amplitude `a_inc`. The argument is a world-space point on the
    /// port face (typically the face centroid); the function should
    /// return the **tangential** incident-mode E-field at that point.
    /// Components normal to the face are dropped by the dot product
    /// inside [`crate::element::assemble_port_modal_rhs`].
    pub modal_e_t: Box<dyn Fn(Vector3<f64>) -> Vector3<f64> + Send + Sync>,
}

impl std::fmt::Debug for PortDefinition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PortDefinition")
            .field("beta_mode", &"<fn(f64) -> f64>")
            .field("modal_e_t", &"<fn(Vector3<f64>) -> Vector3<f64>>")
            .finish()
    }
}

/// Output of [`OpenBoundarySolver::assemble_driven_system`] — the
/// fully-assembled driven linear system at a single frequency, with
/// the lift map needed to scatter the solution back to the full edge
/// basis.
///
/// This is the introspection surface for unit tests of the assembled
/// matrix shape, sparsity, and entry values. Production callers go
/// through [`OpenBoundarySolver::solve_at_frequency`] which factors and
/// back-substitutes in one step.
#[derive(Debug, Clone)]
pub struct DrivenSystem {
    /// Complex sparse driven matrix `A(ω) = K(ω) − k₀² M(ω) +
    /// boundary terms` on the interior-DoF basis (PEC edges eliminated).
    pub matrix: SparseColMat<usize, Complex64>,
    /// Complex RHS vector `b(ω) = Σ_port modal-current contributions`
    /// on the same interior-DoF basis.
    pub rhs: Vec<Complex64>,
    /// Lift map: `interior_edges[i]` is the global edge index of the
    /// interior-DoF at row/column `i`. Identical to
    /// [`crate::AssembledMatricesComplex::interior_edges`].
    pub interior_edges: Vec<usize>,
    /// Inverse lift map: `interior_dof_of_edge[gid]` is `Some(dof)` if
    /// global edge `gid` is an interior DoF, or `None` if it is
    /// PEC-eliminated. Sized to cover every global edge index in the
    /// mesh's edge table.
    pub interior_dof_of_edge: Vec<Option<usize>>,
}

/// Phase 4.fem.eig.2 open-boundary FEM driven solver.
///
/// Wraps a borrowed [`yee_mesh::TetMesh3D`], a [`MaterialDatabase`], a
/// per-exterior-face [`FaceKind`] tagging, and a vector of
/// [`PortDefinition`] descriptors. Construction performs face
/// classification once (no per-frequency cost); each call to
/// [`Self::solve_at_frequency`] assembles the complex driven system
/// `A(ω) e = b(ω)` and returns the interior-edge-indexed complex
/// solution.
///
/// The struct intentionally does not own the `mesh` (lifetime `'m`) so
/// the caller can rebuild it cheaply across multiple solver
/// configurations.
pub struct OpenBoundarySolver<'m> {
    /// Borrowed mesh.
    mesh: &'m TetMesh3D,
    /// Per-tet material database keyed by [`yee_mesh::MaterialTag`].
    material_db: MaterialDatabase,
    /// Per-exterior-face classification. The length and ordering match
    /// the exterior-face list produced by `ExteriorFaceTable::build`.
    /// The caller passes the list in the same order
    /// [`OpenBoundarySolver::new`] resolves the exterior face order, so
    /// `face_kinds[i]` tags `exterior_faces.faces[i]`.
    face_kinds: Vec<FaceKind>,
    /// Wave-port descriptors; index = [`PortId`] referenced by
    /// `FaceKind::WavePort(p)` in [`Self::face_kinds`].
    ports: Vec<PortDefinition>,
    /// Pre-computed exterior-face table — face list with outward
    /// normals, three-vertex tuples, three global-edge indices per face,
    /// and per-edge orientation signs against the global canonical
    /// (lower-endpoint-first) orientation.
    exterior_faces: ExteriorFaceTable,
    /// Global edges classified as PEC (Dirichlet-eliminated) — the
    /// union of every edge on every `FaceKind::Pec`-tagged face.
    /// Computed once at construction; consumed by every
    /// [`Self::solve_at_frequency`] call.
    pec_global_edges: HashSet<usize>,
}

impl<'m> OpenBoundarySolver<'m> {
    /// Build an [`OpenBoundarySolver`] from a mesh, face-kind tagging,
    /// port descriptors, and material database.
    ///
    /// # Arguments
    ///
    /// * `mesh` — tet mesh with stable exterior-face ordering. The
    ///   exterior-face list is computed by walking
    ///   [`yee_mesh::TetMesh3D::tetrahedra`] and emitting every face
    ///   shared by exactly one tet.
    /// * `face_kinds` — per-exterior-face classification. Length must
    ///   match the exterior-face count; index `i` tags exterior face
    ///   `i` (the iteration order is canonical — see
    ///   `ExteriorFaceTable::build`). Unannotated callers may pass
    ///   `vec![FaceKind::Pec; n_exterior_faces]` to reproduce the
    ///   closed-cavity boundary condition.
    /// * `ports` — wave-port descriptors. Indexed by the `PortId` in
    ///   each `FaceKind::WavePort(p)` tag. May be empty if no
    ///   wave-port faces are present.
    /// * `material_db` — per-tet `(ε(ω), μ(ω))` lookup, consumed by the
    ///   inner [`FemEigenAssembly::assemble_complex`] call.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if `face_kinds.len()` does not match
    /// the mesh's exterior-face count, or if any
    /// `FaceKind::WavePort(p)` references a `p >= ports.len()`. The
    /// face classifier itself does not fail on a well-formed mesh.
    pub fn new(
        mesh: &'m TetMesh3D,
        face_kinds: Vec<FaceKind>,
        ports: Vec<PortDefinition>,
        material_db: MaterialDatabase,
    ) -> Result<Self, Error> {
        if mesh.tetrahedra.is_empty() {
            return Err(Error::Invalid(
                "OpenBoundarySolver::new: mesh has zero tetrahedra".to_string(),
            ));
        }

        let exterior_faces = ExteriorFaceTable::build(mesh);

        if face_kinds.len() != exterior_faces.faces.len() {
            return Err(Error::Invalid(format!(
                "OpenBoundarySolver::new: face_kinds length {} does not match \
                 exterior-face count {}",
                face_kinds.len(),
                exterior_faces.faces.len()
            )));
        }

        for (i, kind) in face_kinds.iter().enumerate() {
            if let FaceKind::WavePort(p) = *kind
                && p >= ports.len()
            {
                return Err(Error::Invalid(format!(
                    "OpenBoundarySolver::new: face {i} is WavePort({p}) but \
                     only {} ports are defined",
                    ports.len()
                )));
            }
        }

        // PEC Dirichlet set: union of all global edges of every PEC
        // face. Built once at construction so the per-frequency solve
        // path is a constant-time lookup.
        let mut pec_global_edges: HashSet<usize> = HashSet::new();
        for (i, kind) in face_kinds.iter().enumerate() {
            if matches!(kind, FaceKind::Pec) {
                let face = &exterior_faces.faces[i];
                for &gid in &face.global_edges {
                    pec_global_edges.insert(gid);
                }
            }
        }

        Ok(Self {
            mesh,
            material_db,
            face_kinds,
            ports,
            exterior_faces,
            pec_global_edges,
        })
    }

    /// Number of exterior faces classified by this solver. Equal to
    /// `face_kinds.len()`.
    pub fn n_exterior_faces(&self) -> usize {
        self.face_kinds.len()
    }

    /// Read-only borrow of the per-exterior-face [`FaceKind`] tagging
    /// supplied at construction.
    pub fn face_kinds(&self) -> &[FaceKind] {
        &self.face_kinds
    }

    /// Read-only borrow of the PEC global-edge set — the union of all
    /// global-edge indices that lie on at least one [`FaceKind::Pec`]
    /// face. Used internally to enforce PEC precedence over wave-port
    /// modal source on shared edges.
    pub fn pec_global_edges(&self) -> &HashSet<usize> {
        &self.pec_global_edges
    }

    /// Read-only borrow of the borrowed mesh.
    pub fn mesh(&self) -> &'m TetMesh3D {
        self.mesh
    }

    /// Read-only borrow of the wave-port descriptors.
    pub fn ports(&self) -> &[PortDefinition] {
        &self.ports
    }

    /// Assemble the driven open-boundary system at angular frequency
    /// `omega` without solving it. Returns the complex sparse driven
    /// matrix `A(ω)`, the complex RHS vector `b(ω)`, and the
    /// interior-edge lift map.
    ///
    /// This is the introspection surface for unit tests that need to
    /// verify the assembled matrix shape or sparsity (e.g.
    /// "all-PEC matches the closed-cavity assemble", "ABC face
    /// introduces complex boundary terms", "single wave-port modal RHS
    /// is non-zero"). Production callers should go through
    /// [`Self::solve_at_frequency`] which factors and back-substitutes
    /// in one step.
    pub fn assemble_driven_system(&self, omega: f64) -> Result<DrivenSystem, Error> {
        if omega <= 0.0 {
            return Err(Error::Invalid(format!(
                "OpenBoundarySolver::assemble_driven_system: omega = {omega} must be positive"
            )));
        }

        // ---- Step 1: assemble the PEC-reduced complex K(ω), M(ω).
        // Use the explicit-PEC-edges variant so only edges on PEC-tagged
        // faces are eliminated; edges on ABC / wave-port faces remain
        // as interior DoFs and receive the boundary-term scatter
        // below.
        let n_tets = self.mesh.tetrahedra.len();
        let assembly = FemEigenAssembly::new(self.mesh, vec![1.0; n_tets], vec![1.0; n_tets])?;
        let assembled = assembly.assemble_complex_with_pec_edges(
            omega,
            &self.material_db,
            &self.pec_global_edges,
        )?;

        let n_interior = assembled.interior_edges.len();
        // Inverse map: global-edge-index → interior-DoF index. `None`
        // means the edge is PEC-eliminated. Sized by both the
        // interior-edge lift map and the exterior-face table to ensure
        // any PEC edge with a higher global index than every interior
        // edge still fits.
        let max_interior = assembled.interior_edges.iter().copied().max().unwrap_or(0);
        let max_face = self.total_global_edge_count();
        let n_global_edges = max_interior.max(max_face.saturating_sub(1)) + 1;
        let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; n_global_edges];
        for (dof, &gid) in assembled.interior_edges.iter().enumerate() {
            interior_dof_of_edge[gid] = Some(dof);
        }

        // ---- Step 2: form the closed-cavity driven core
        // A(ω) = K(ω) − k₀² M(ω). Both K and M live on the interior-DoF
        // basis already; the subtraction is a sparse axpy. We collect
        // entries into a triplet list because we will append face-block
        // contributions next and let `try_new_from_triplets` accumulate
        // duplicates.
        let k0 = omega / C0;
        let k0_sq = k0 * k0;
        let mut triplets: Vec<Triplet<usize, usize, Complex64>> =
            Vec::with_capacity(assembled.k.nnz() + assembled.m.nnz());
        for (row, col, &val) in assembled.k.triplet_iter() {
            triplets.push(Triplet::new(row, col, val));
        }
        let k0_sq_c = Complex64::new(k0_sq, 0.0);
        for (row, col, &val) in assembled.m.triplet_iter() {
            triplets.push(Triplet::new(row, col, -k0_sq_c * val));
        }

        // RHS vector on the interior-DoF basis. Starts at zero;
        // wave-port faces accumulate modal-current contributions.
        let mut rhs: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_interior];

        // ---- Step 3+4: scatter per-face boundary-term contributions
        // (ABC and wave-port). PEC faces contribute nothing — they are
        // already eliminated by the row/column drop inside
        // `assemble_complex`. ----------------------------------------
        for (i, kind) in self.face_kinds.iter().enumerate() {
            let face = &self.exterior_faces.faces[i];
            match *kind {
                FaceKind::Pec => {
                    // No assembly contribution. Edges on this face are
                    // tangential-E-zero by the Dirichlet elimination
                    // already applied to K(ω) and M(ω) above.
                }
                FaceKind::Abc => {
                    self.scatter_abc_face(face, k0, &interior_dof_of_edge, &mut triplets);
                }
                FaceKind::WavePort(p) => {
                    let port = &self.ports[p];
                    let beta = (port.beta_mode)(omega);
                    let centroid = face.centroid(self.mesh);
                    let e_t = (port.modal_e_t)(centroid);
                    self.scatter_port_face(
                        face,
                        beta,
                        e_t,
                        &interior_dof_of_edge,
                        &mut triplets,
                        &mut rhs,
                    );
                }
            }
        }

        let matrix = SparseColMat::try_new_from_triplets(n_interior, n_interior, &triplets)
            .map_err(|e| {
                Error::Numerical(format!(
                    "OpenBoundarySolver::assemble_driven_system: failed to build driven matrix: {e:?}"
                ))
            })?;

        Ok(DrivenSystem {
            matrix,
            rhs,
            interior_edges: assembled.interior_edges,
            interior_dof_of_edge,
        })
    }

    /// Solve the driven open-boundary system
    /// `A(ω) e = b(ω)` at a single angular frequency, where
    ///
    /// ```text
    ///     A(ω) = K(ω) − k₀² M(ω)
    ///            + Σ_ABC  j k₀ B_ABC^{face}
    ///            + Σ_port j β  B_port^{face}
    /// ```
    ///
    /// and
    ///
    /// ```text
    ///     b(ω) = Σ_port  + 2 j β · ∫_face N_i · e_t(x) dS.
    /// ```
    ///
    /// Returns the complex edge-DoF solution vector on the **interior**
    /// edge basis (PEC edges are eliminated by the inner
    /// [`FemEigenAssembly::assemble_complex`] row/column drop).
    /// S-parameter extraction (project the FEM solution against the
    /// modal profile on each port face) lands in
    /// Phase 4.fem.eig.2 step E4.
    ///
    /// # Arguments
    ///
    /// * `omega` — real-valued angular frequency (rad/s). Real-valued
    ///   for the same reason
    ///   [`FemEigenAssembly::assemble_complex`] is — the per-tet
    ///   [`MaterialDatabase::eps_at`] lookup evaluates on a real
    ///   frequency.
    ///
    /// # Errors
    ///
    /// Propagates [`Error`] from
    /// [`FemEigenAssembly::assemble_complex`] (mesh / material shape
    /// mismatch, empty mesh) and adds an
    /// [`Error::Numerical`] variant when the sparse LU on the driven
    /// matrix fails.
    pub fn solve_at_frequency(&self, omega: f64) -> Result<Vec<Complex64>, Error> {
        let system = self.assemble_driven_system(omega)?;
        let n_interior = system.rhs.len();

        let lu: Lu<usize, Complex64> = system.matrix.sp_lu().map_err(|e| {
            Error::Numerical(format!(
                "OpenBoundarySolver::solve_at_frequency: sparse LU of driven matrix failed: {e:?}"
            ))
        })?;

        let mut rhs_mat = faer::Mat::<Complex64>::zeros(n_interior, 1);
        for (i, &b_i) in system.rhs.iter().enumerate() {
            rhs_mat[(i, 0)] = b_i;
        }
        lu.solve_in_place_with_conj(faer::Conj::No, rhs_mat.as_mut());

        let mut out = vec![Complex64::new(0.0, 0.0); n_interior];
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = rhs_mat[(i, 0)];
        }
        Ok(out)
    }

    /// Scatter the per-face ABC `+ j k₀ B_ABC` block into the driven
    /// triplet list. PEC edges on the face are silently skipped — they
    /// are eliminated by the global row/column drop applied to `K(ω)`
    /// and `M(ω)` by [`FemEigenAssembly::assemble_complex`].
    fn scatter_abc_face(
        &self,
        face: &ExteriorFace,
        k0: f64,
        interior_dof_of_edge: &[Option<usize>],
        triplets: &mut Vec<Triplet<usize, usize, Complex64>>,
    ) {
        let block = assemble_abc_face_block(face.world_vertices(self.mesh), face.normal, k0, 1.0);
        for i in 0..3 {
            let gi = face.global_edges[i];
            let Some(ii) = interior_dof_of_edge[gi] else {
                continue;
            };
            let si = face.signs[i];
            for j in 0..3 {
                let gj = face.global_edges[j];
                let Some(jj) = interior_dof_of_edge[gj] else {
                    continue;
                };
                let sj = face.signs[j];
                let sign = Complex64::new(si * sj, 0.0);
                triplets.push(Triplet::new(ii, jj, sign * block[(i, j)]));
            }
        }
    }

    /// Scatter the per-face wave-port `+ j β B_port` stiffness block
    /// and `+ 2 j β · ∫ N_i · e_t dS` RHS contribution into the driven
    /// triplet list and RHS vector. PEC edges on the face are silently
    /// skipped (the PEC-precedence rule from spec §10 risk #5).
    fn scatter_port_face(
        &self,
        face: &ExteriorFace,
        beta: f64,
        e_t: Vector3<f64>,
        interior_dof_of_edge: &[Option<usize>],
        triplets: &mut Vec<Triplet<usize, usize, Complex64>>,
        rhs: &mut [Complex64],
    ) {
        let face_vertices = face.world_vertices(self.mesh);

        // Stiffness contribution: + j β B_port^{face}.
        let block = assemble_port_face_block(face_vertices, face.normal, beta, 1.0);
        for i in 0..3 {
            let gi = face.global_edges[i];
            // PEC precedence: skip edges that lie on a PEC face even if
            // they also lie on this wave-port face. Without this guard
            // the modal source would conflict with the PEC tangential-
            // zero condition (spec §10 risk #5).
            if self.pec_global_edges.contains(&gi) {
                continue;
            }
            let Some(ii) = interior_dof_of_edge[gi] else {
                continue;
            };
            let si = face.signs[i];
            for j in 0..3 {
                let gj = face.global_edges[j];
                if self.pec_global_edges.contains(&gj) {
                    continue;
                }
                let Some(jj) = interior_dof_of_edge[gj] else {
                    continue;
                };
                let sj = face.signs[j];
                let sign = Complex64::new(si * sj, 0.0);
                triplets.push(Triplet::new(ii, jj, sign * block[(i, j)]));
            }
        }

        // RHS contribution: + 2 j β · ∫ N_i · e_t dS (face-centroid
        // quadrature; per element-layer docs).
        let rhs_block = assemble_port_modal_rhs(face_vertices, face.normal, beta, e_t);
        for i in 0..3 {
            let gi = face.global_edges[i];
            if self.pec_global_edges.contains(&gi) {
                continue;
            }
            let Some(ii) = interior_dof_of_edge[gi] else {
                continue;
            };
            let si = face.signs[i];
            let sign = Complex64::new(si, 0.0);
            rhs[ii] += sign * rhs_block[i];
        }
    }
}

impl OpenBoundarySolver<'_> {
    /// Total global-edge count for this mesh — the maximum global-edge
    /// index referenced by any tet plus one. Used internally to size
    /// the `global-edge → interior-DoF` inverse lookup table; the
    /// `interior_edges` lift map alone may under-size the table when
    /// PEC-eliminated edges happen to carry the largest global index.
    fn total_global_edge_count(&self) -> usize {
        let mut max_idx = 0usize;
        for face in &self.exterior_faces.faces {
            for &gid in &face.global_edges {
                if gid > max_idx {
                    max_idx = gid;
                }
            }
        }
        max_idx + 1
    }
}

// ---------------------------------------------------------------------
// Exterior-face table — local to this module
// ---------------------------------------------------------------------

/// Canonical edge key — lower-endpoint-first vertex pair. Peer of
/// [`crate::assembly`]'s private `EdgeKey`; reproduced here so the
/// exterior-face classifier can map face edges to the same global edge
/// index space used by [`FemEigenAssembly::assemble_complex`].
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

/// One exterior face of a [`yee_mesh::TetMesh3D`]. Carries the three
/// vertex indices in CCW order seen from outside the mesh, the outward
/// normal, the three global edge indices (matching the
/// [`crate::assembly`] global-edge numbering exactly), and the three
/// per-edge orientation signs.
#[derive(Debug, Clone, Copy)]
struct ExteriorFace {
    /// The three vertex indices of the face, in CCW order seen from
    /// outside the mesh (i.e. the outward normal is given by
    /// `((v1 − v0) × (v2 − v0)).normalize()`).
    vertices: [usize; 3],
    /// Pre-computed outward unit normal.
    normal: Vector3<f64>,
    /// Global edge indices of the three face edges, in the canonical
    /// face-local edge order `(v0 → v1, v1 → v2, v2 → v0)`.
    global_edges: [usize; 3],
    /// Per-edge orientation sign: `+1` if the face-local edge direction
    /// matches the canonical global orientation (lower-vertex-first),
    /// `-1` otherwise. Applied during scatter to reconcile the
    /// canonical-orientation element-layer output with the global
    /// orientation.
    signs: [f64; 3],
}

impl ExteriorFace {
    /// World-space coordinates of the three face vertices in the order
    /// matching [`Self::vertices`].
    fn world_vertices(&self, mesh: &TetMesh3D) -> [Vector3<f64>; 3] {
        [
            mesh.vertices[self.vertices[0]],
            mesh.vertices[self.vertices[1]],
            mesh.vertices[self.vertices[2]],
        ]
    }

    /// Face centroid — arithmetic mean of the three face vertices.
    fn centroid(&self, mesh: &TetMesh3D) -> Vector3<f64> {
        let v = self.world_vertices(mesh);
        (v[0] + v[1] + v[2]) / 3.0
    }
}

/// Exterior-face table for a [`yee_mesh::TetMesh3D`].
///
/// Built once per [`OpenBoundarySolver::new`]; iterates over every tet,
/// identifies faces shared by exactly one tet (the boundary faces),
/// orients them outward, and resolves the three face edges to the same
/// global-edge index space used by
/// [`FemEigenAssembly::assemble_complex`].
#[derive(Debug, Clone)]
struct ExteriorFaceTable {
    /// All exterior faces of the mesh, in canonical enumeration order
    /// — see [`Self::build`] for the exact ordering.
    faces: Vec<ExteriorFace>,
}

impl ExteriorFaceTable {
    /// Build the exterior-face table from a [`yee_mesh::TetMesh3D`].
    ///
    /// Algorithm (mirrors [`crate::assembly::TetEdgeTable::build`]'s
    /// boundary-edge classifier one dimension up):
    ///
    /// 1. Build the global edge table by walking every tet's six
    ///    canonical local edges per [`crate::element::LOCAL_EDGES`] and
    ///    de-duplicating with a `HashMap<EdgeKey, usize>`. This is the
    ///    same global-edge numbering [`FemEigenAssembly`] consumes.
    /// 2. Walk every tet's four faces (the face opposite each local
    ///    vertex). For each face, build a sorted-vertex-triplet key and
    ///    a count of how many tets reference it. Faces with count `1`
    ///    are exterior.
    /// 3. For each exterior face, orient it outward: compute the face
    ///    centroid and the parent tet's centroid; the outward normal
    ///    points away from the tet centroid. If the CCW-ordered face
    ///    normal points the wrong way, swap two vertices.
    /// 4. Resolve the three face edges to global edge indices and
    ///    record the orientation sign per edge.
    ///
    /// The exterior-face enumeration order is the order in which the
    /// algorithm encounters new boundary faces while walking the tet
    /// list — deterministic, stable across runs for a given mesh.
    fn build(mesh: &TetMesh3D) -> Self {
        // ---- Step 1: build the global edge table identical to
        // crate::assembly::TetEdgeTable. ------------------------------
        let mut edge_map: HashMap<EdgeKey, usize> = HashMap::new();
        for tet in &mesh.tetrahedra {
            for &(li, lj) in LOCAL_EDGES.iter() {
                let a = tet[li];
                let b = tet[lj];
                let key = EdgeKey::new(a, b);
                let next_idx = edge_map.len();
                edge_map.entry(key).or_insert(next_idx);
            }
        }

        // ---- Step 2: face-incidence map. ----------------------------
        const TET_FACES: [[usize; 3]; 4] = [
            [1, 2, 3], // face opposite v0
            [0, 2, 3], // face opposite v1
            [0, 1, 3], // face opposite v2
            [0, 1, 2], // face opposite v3
        ];

        // Map from sorted-vertex-triplet → list of (tet_idx,
        // local_face_idx). Boundary faces appear exactly once.
        let mut face_map: HashMap<[usize; 3], Vec<(usize, usize)>> = HashMap::new();
        for (tet_idx, tet) in mesh.tetrahedra.iter().enumerate() {
            for (face_local, &[a, b, c]) in TET_FACES.iter().enumerate() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                face_map.entry(key).or_default().push((tet_idx, face_local));
            }
        }

        // ---- Step 3+4: walk every tet and emit exterior faces in
        // canonical (tet-then-local-face) order. The face_map serves
        // only as the "is this face boundary?" lookup.
        let mut faces: Vec<ExteriorFace> = Vec::new();
        let mut seen: HashSet<[usize; 3]> = HashSet::new();
        for (tet_idx, tet) in mesh.tetrahedra.iter().enumerate() {
            for (face_local, &[a, b, c]) in TET_FACES.iter().enumerate() {
                let raw_verts = [tet[a], tet[b], tet[c]];
                let mut key = raw_verts;
                key.sort_unstable();
                let occurrences = face_map.get(&key).map(|v| v.len()).unwrap_or(0);
                if occurrences != 1 || !seen.insert(key) {
                    continue;
                }

                // Tet centroid for the outward-normal check.
                let v_tet = [
                    mesh.vertices[tet[0]],
                    mesh.vertices[tet[1]],
                    mesh.vertices[tet[2]],
                    mesh.vertices[tet[3]],
                ];
                let tet_centroid = (v_tet[0] + v_tet[1] + v_tet[2] + v_tet[3]) / 4.0;

                // Initial face vertex order from TET_FACES.
                let mut verts = raw_verts;
                let p = [
                    mesh.vertices[verts[0]],
                    mesh.vertices[verts[1]],
                    mesh.vertices[verts[2]],
                ];
                let face_centroid = (p[0] + p[1] + p[2]) / 3.0;
                let face_to_tet = tet_centroid - face_centroid;
                let candidate_normal = (p[1] - p[0]).cross(&(p[2] - p[0]));
                // Outward iff candidate_normal points away from the
                // tet centroid (i.e. opposite to face_to_tet).
                if candidate_normal.dot(&face_to_tet) > 0.0 {
                    // Wrong orientation: swap two vertices to flip the
                    // cross-product sign.
                    verts.swap(1, 2);
                }

                let p = [
                    mesh.vertices[verts[0]],
                    mesh.vertices[verts[1]],
                    mesh.vertices[verts[2]],
                ];
                let mut normal = (p[1] - p[0]).cross(&(p[2] - p[0]));
                let n_norm = normal.norm();
                if n_norm > 0.0 {
                    normal /= n_norm;
                }

                // Resolve the three face edges to global edge indices
                // and compute the per-edge orientation sign.
                let edges_local = [(0usize, 1usize), (1, 2), (2, 0)];
                let mut global_edges = [0usize; 3];
                let mut signs = [0.0f64; 3];
                for (i, &(la, lb)) in edges_local.iter().enumerate() {
                    let a_g = verts[la];
                    let b_g = verts[lb];
                    let key = EdgeKey::new(a_g, b_g);
                    let gid = *edge_map.get(&key).expect(
                        "face edge missing from global edge table — \
                         exterior face references an edge not visited by any tet \
                         (bug in ExteriorFaceTable::build vs assembly::TetEdgeTable::build)",
                    );
                    global_edges[i] = gid;
                    // Local direction is a_g → b_g; canonical global
                    // direction is from-to with from < to. Sign is +1
                    // iff a_g < b_g.
                    signs[i] = if a_g < b_g { 1.0 } else { -1.0 };
                }

                let _ = tet_idx; // tet index is implicit in enumeration order
                let _ = face_local;
                faces.push(ExteriorFace {
                    vertices: verts,
                    normal,
                    global_edges,
                    signs,
                });
            }
        }

        Self { faces }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two tets sharing a triangular face — same fixture as
    /// `crate::assembly::tests::two_tets_shared_face` so any drift in
    /// global edge numbering between this module's
    /// [`ExteriorFaceTable`] and `crate::assembly::TetEdgeTable` is
    /// visible.
    fn two_tets_shared_face_mesh() -> TetMesh3D {
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
    fn exterior_face_table_two_tets_has_six_exterior_faces() {
        let mesh = two_tets_shared_face_mesh();
        let table = ExteriorFaceTable::build(&mesh);
        // Two tets sharing one face — each tet contributes 3 exterior
        // faces (the 4 faces minus the shared one). Total = 6.
        assert_eq!(
            table.faces.len(),
            6,
            "two tets sharing one face should produce 6 exterior faces"
        );
    }

    #[test]
    fn exterior_face_outward_normals_point_away_from_mesh_centroid() {
        let mesh = two_tets_shared_face_mesh();
        let mesh_centroid = mesh
            .vertices
            .iter()
            .copied()
            .fold(Vector3::zeros(), |a, b| a + b)
            / (mesh.vertices.len() as f64);
        let table = ExteriorFaceTable::build(&mesh);
        for face in &table.faces {
            let face_centroid = face.centroid(&mesh);
            let outward_dir = face_centroid - mesh_centroid;
            let dot = face.normal.dot(&outward_dir);
            // Both `outward_dir` and `face.normal` should point away
            // from the mesh interior; their dot product should be
            // non-negative. The strict inequality holds for every face
            // of this fixture because no face passes through the mesh
            // centroid.
            assert!(
                dot >= -1e-12,
                "face normal {:?} should point away from mesh centroid; \
                 got dot = {dot} with outward_dir = {outward_dir:?}",
                face.normal
            );
        }
    }
}
