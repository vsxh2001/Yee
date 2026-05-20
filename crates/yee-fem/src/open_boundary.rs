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
//! ## S-parameter extraction (Phase 4.fem.eig.2 step E4 + CCCCCCCCC fix)
//!
//! [`OpenBoundarySolver::sweep`] runs the per-frequency driven solve at
//! every `ω` in the supplied list and extracts the diagonal scattering
//! matrix entries `S_{p,p}(ω)` via modal projection on each port face.
//! Per Pozar §3.3 / Jin §10.7, with normalised incident amplitude
//! `a_inc_p = 1` and modal self-inner-product
//! `M_pp = ⟨e_mode_p, e_mode_p⟩_port` computed via the same face-
//! centroid quadrature:
//!
//! ```text
//!     b_p(ω)     =   ⟨ E_FEM,t , e_mode_p ⟩_port / M_pp  −  a_inc_p,
//!     S_{p,p}(ω) =   b_p(ω) / a_inc_p.
//! ```
//!
//! The modal normalisation by `M_pp` (CCCCCCCCC scaling fix) replaces
//! the original spec §4.3 formula `b_p = 2 ⟨E_FEM, e_mode⟩ − a_inc`,
//! which had implicitly assumed `M_pp = 1/2`. With the standard
//! Pozar §3.3 orthonormalisation `M_pp ≈ 1` used by the driver
//! (`crates/yee-validation/src/lib.rs::fem_eig_003_modal_e_t_te10`),
//! the un-normalised formula saturated `|S_{11}|` at 1.0 even on a
//! matched-port total field (BBBBBBBBB E5 finding). The corrected
//! formula recovers the Pozar §3.3 matched-port identity
//! `S_{11} ≈ 0` for `E_FEM,t ≈ a_inc · e_mode` and the PEC-reflection
//! identity `|S_{11}| ≈ 1` for `E_FEM,t ≈ 0`.
//!
//! The modal projection is computed by face-centroid quadrature: for
//! every port face `f` carrying the port `p` tag,
//!
//! ```text
//!     ⟨ E_FEM,t , e_mode_p ⟩_port
//!         ≈  Σ_face  A_face · ( E_FEM,t(centroid_f) · e_mode_p(centroid_f) ),
//! ```
//!
//! and `E_FEM,t(centroid_f)` is reconstructed from the per-edge complex
//! DoFs of the three face edges by evaluating the Whitney-1 face basis
//! at the centroid:
//!
//! ```text
//!     E_FEM,t(centroid)  =  Σ_{i ∈ face_edges}  s_i · e_i · (t_i / 3),
//! ```
//!
//! where `t_i = v_{(i+1) mod 3} − v_i` is the canonical face-edge
//! tangent, `s_i ∈ {-1, +1}` is the local-to-global orientation sign,
//! `e_i` is the interior-DoF complex amplitude (or `0` if edge `i` is
//! PEC-eliminated), and the `1/3` is the Whitney-1 edge basis value at
//! the centroid (each edge basis integrates to `A/3` against a constant
//! test function — see [`crate::element::assemble_port_modal_rhs`] for
//! the dual formulation). Cross-port `S_{p,q}` for `p ≠ q` is deferred
//! to Phase 4.fem.eig.2.0.2 per spec §13.
//!
//! ## References
//!
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
//!   §4 (theory), §6 (API surface).
//! * `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md`
//!   step E3, step E4 (S-parameter extraction).
//! * `docs/src/decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md`.
//! * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §3.3
//!   — wave-port modal characterisation and `S_{11}` extraction
//!   convention.
//! * Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd
//!   ed., Wiley 2014, §10.5 — modal decomposition for FEM wave-port
//!   driven analysis.

use std::collections::{HashMap, HashSet};

use faer::linalg::solvers::SolveCore;
use faer::sparse::{SparseColMat, Triplet, linalg::solvers::Lu};
use nalgebra::{DMatrix, Vector3};
use num_complex::Complex64;
use yee_core::Error;
use yee_core::units::C0;
use yee_mesh::TetMesh3D;

use nalgebra::SMatrix;

use crate::assembly::FemEigenAssembly;
use crate::element::{
    LOCAL_EDGES, assemble_abc_face_block, assemble_abc2_face_block, assemble_port_face_block,
    assemble_port_face_block_gauss_pts, assemble_port_face_rhs_gauss_pts, assemble_port_modal_rhs,
    assemble_tet_element_complex, assemble_tet_element_complex_anisotropic,
};
use crate::material::MaterialDatabase;

/// Three-point Gauss-quadrature barycentric coordinates on the
/// reference triangle (mirror of `element::TRI_GAUSS_3PT_BARY` —
/// kept private here so the open-boundary helper can sample modal
/// profiles at the same Gauss-point world-space positions as the
/// element-layer F1 helpers consume).
///
/// Each row is `(λ_0, λ_1, λ_2)` for one Gauss point; the
/// corresponding weight is `A / 3` (uniform).
const TRI_GAUSS_3PT_BARY: [[f64; 3]; 3] = [
    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
];

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

/// Selects the open-boundary truncation kernel on
/// [`FaceKind::Abc`]-tagged exterior faces.
///
/// The default is [`AbcOrder::First`], which reproduces the
/// Phase 4.fem.eig.2 v2 + CCCCCCCCC behaviour bit-for-bit: every ABC
/// face contributes the 1st-order Mur block `+ j k₀ · (A / μ_r) · R_1`
/// via [`crate::element::assemble_abc_face_block`]. The reflection
/// floor for a TE plane wave at normal incidence is `~ −40 dB`
/// (Jin §10.4, Table 10.1).
///
/// Selecting [`AbcOrder::Second`] augments the bilinear form with the
/// tangential-curl correction `−(1 / (2 k₀)) · (A / μ_r) · R_2` from
/// Engquist–Majda 1979 eq. 9, lowering the normal-incidence reflection
/// floor to `~ −60 dB`. The 2nd-order block is computed by
/// [`crate::element::assemble_abc2_face_block`]; the curl correction
/// has a **real** scalar prefactor while the 1st-order part stays
/// purely imaginary, so the composite block is complex-symmetric with
/// non-trivial real *and* imaginary content.
///
/// [`AbcOrder::CfsPml`] (Phase 4.fem.eig.3.5) replaces the
/// surface-integral Engquist–Majda kernel with a volumetric CFS-PML
/// (Roden–Gedney 2000) buffer-layer absorber. The PML is a thin shell
/// of additional tetrahedra outside the original cavity volume, in
/// which the constitutive tensor `ε(ω)` becomes the stretched-coordinate
/// form `ε · Λ(ω)` and absorbs off-normal and evanescent modal content
/// that the local Engquist–Majda operators cannot. See [`PmlConfig`]
/// for the grading-parameter surface.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub enum AbcOrder {
    /// 1st-order Engquist–Majda absorbing boundary. The v0 / v2 default
    /// — reproduces the shipped behaviour bit-for-bit. Reflection floor
    /// `~ −40 dB` at normal incidence (Jin §10.4).
    #[default]
    First,
    /// 2nd-order Engquist–Majda absorbing boundary. Adds the
    /// tangential-curl correction term per Engquist–Majda 1979 eq. 9.
    /// Reflection floor `~ −60 dB` at normal incidence (Jin §10.4).
    Second,
    /// CFS-PML (Roden–Gedney 2000) volumetric buffer-layer absorber
    /// (Phase 4.fem.eig.3.5). Replaces the surface-integral
    /// Engquist–Majda kernel with a thin shell of extra tetrahedra
    /// outside the original cavity in which `ε → ε · Λ(ω)` with
    /// stretched-coordinate factor
    /// `s_α(ω) = κ_α + σ_α / (α_α + j ω ε_0)`. The variant payload
    /// carries the grading parameters; see [`PmlConfig`].
    CfsPml(PmlConfig),
}

/// CFS-PML grading-parameter configuration (Phase 4.fem.eig.3.5;
/// Roden–Gedney 2000, *IEEE MWCL* 10:5).
///
/// One [`PmlConfig`] applies symmetrically to every PML-tagged face.
/// The grading parameters control the polynomial-graded conductivity
/// `σ_α(d) = σ_max · (d/D)^m`, coordinate stretching
/// `κ_α(d) = 1 + (κ_max − 1) · (d/D)^m`, and the CFS frequency-shift
/// `α_α(d) = α_max`. The depth `d ∈ [0, D]` is measured inward from the
/// PML's outer truncation surface; at the inner boundary (`d = 0`) both
/// `σ` and `κ − 1` vanish so the material is continuous with the
/// cavity interior, eliminating surface-reflection spurious modes.
///
/// Default values follow Roden–Gedney 2000 §III + Table I for microwave
/// waveguide benchmarks. `sigma_max` and `alpha_max` use sentinel
/// zeros at construction; the
/// [`OpenBoundarySolver::with_cfs_pml`] builder recomputes them from
/// the band-centre `ω` and mean tet edge length using the recommended
/// `σ_max ≈ (m + 1) / (150 π · h · √ε_r)` and `α_max ≈ ω₀ ε_0`
/// formulae. Callers wanting full control can populate either field
/// explicitly; the [`Self::resolved`] helper returns a fully-populated
/// copy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PmlConfig {
    /// PML shell thickness in tet layers. Default 6 (Roden–Gedney 2000
    /// §III recommends 6 to 10 cells for microwave applications).
    pub thickness_cells: usize,
    /// Maximum conductivity (S/m) at the outer truncation surface.
    /// `0.0` is the sentinel meaning "recompute from frequency and
    /// mesh"; [`Self::resolved`] fills this in.
    pub sigma_max: f64,
    /// CFS frequency-shift parameter `α_max` (rad·s⁻¹ × ε₀, i.e. the
    /// same units as `j ω ε_0`). `0.0` is the sentinel; [`Self::resolved`]
    /// sets it to `ω₀ ε_0` per Roden–Gedney 2000 §IV.
    pub alpha_max: f64,
    /// Coordinate-stretching parameter `κ_max`. Default 5.0 per
    /// Roden–Gedney 2000 Table I for microwave-waveguide benchmarks.
    pub kappa_max: f64,
    /// Polynomial grading order `m` for `σ_α(d) = σ_max · (d/D)^m`.
    /// Default 3. Values 2, 3, 4 are typical; higher orders steepen the
    /// gradient near the outer truncation surface and the inner cavity
    /// boundary.
    pub m: usize,
}

impl Default for PmlConfig {
    fn default() -> Self {
        Self {
            thickness_cells: 6,
            sigma_max: 0.0,
            alpha_max: 0.0,
            kappa_max: 5.0,
            m: 3,
        }
    }
}

impl PmlConfig {
    /// Resolve sentinel `0.0` values for `sigma_max` and `alpha_max`
    /// against a band-centre frequency `freq_hz` and a mean PML tet
    /// edge length `h_cell` (m).
    ///
    /// Per Roden–Gedney 2000 §III + §IV:
    ///
    /// * `σ_max ≈ (m + 1) / (150 π · h_cell · √ε_r)` with `ε_r = 1`
    ///   (PML built outside an air-filled cavity).
    /// * `α_max ≈ 2 π · freq_hz · ε_0` (band-centre rule of thumb).
    ///
    /// Non-zero `sigma_max` / `alpha_max` are passed through verbatim
    /// (callers retain full control). `kappa_max`, `thickness_cells`,
    /// `m` are never touched.
    pub fn resolved(self, freq_hz: f64, h_cell: f64) -> Self {
        let mut out = self;
        if out.sigma_max == 0.0 {
            let m_plus_1 = (self.m as f64) + 1.0;
            // ε_r = 1 (PML against air).
            out.sigma_max = m_plus_1 / (150.0 * std::f64::consts::PI * h_cell.max(1e-12));
        }
        if out.alpha_max == 0.0 {
            // ω₀ ε_0 with ω₀ = 2 π f₀.
            out.alpha_max = 2.0 * std::f64::consts::PI * freq_hz * yee_core::units::EPS0;
        }
        out
    }
}

/// Caller-supplied designation of which exterior faces become
/// PML-fronted (Phase 4.fem.eig.3.5).
///
/// One [`PmlRegion`] is consumed by
/// [`OpenBoundarySolver::with_cfs_pml`]; if absent, the builder
/// defaults to "every [`FaceKind::Abc`]-tagged face is PML-fronted",
/// which is the standard one-shell-per-ABC pattern used by
/// fem-eig-003 and fem-eig-006. The per-axis [`PmlConfig`] travels in
/// the [`Self::config`] payload; if `None`, the solver-level default
/// applies.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PmlRegion {
    /// Faces in the original mesh's [`FaceKind::Abc`] set that should
    /// be replaced with a PML shell. An empty `faces` set is "every
    /// ABC face becomes PML-fronted" (the spec §6 default).
    pub faces: Vec<FaceKind>,
    /// Per-axis PML configuration. If `None`, the solver-level
    /// [`PmlConfig`] default applies.
    pub config: Option<PmlConfig>,
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

/// Frequency-swept diagonal S-parameter matrix
/// (Phase 4.fem.eig.2 step E4 output of [`OpenBoundarySolver::sweep`]).
///
/// Carries only the **diagonal** entries `S_{p,p}(ω)` per port, one
/// sweep vector per port. Cross-port `S_{p,q}` for `p ≠ q` is deferred
/// to Phase 4.fem.eig.2.0.2 per spec §13 (single-incident-mode-per-port
/// driven analysis ships single-port S-parameters only; multi-port
/// scattering requires per-port driven sweeps with cross-projection,
/// which is out of scope for v0).
///
/// # Layout
///
/// `s_pp.len() == n_ports` and `s_pp[p].len() == omegas.len()` for every
/// `p ∈ [0, n_ports)`. `s_pp[p][k]` is the per-port `S_{p,p}` at
/// `omegas[k]`.
#[derive(Debug, Clone)]
pub struct SParameters {
    /// Real-valued angular frequencies (rad/s) at which the sweep was
    /// evaluated; matches the order of the slice passed to
    /// [`OpenBoundarySolver::sweep`].
    pub omegas: Vec<f64>,
    /// Per-port diagonal S-parameter sweep vectors.
    /// `s_pp[p][k] = S_{p,p}(omegas[k])`. Length of the outer vector
    /// equals the number of [`PortDefinition`]s registered with the
    /// solver; length of each inner vector equals `omegas.len()`.
    pub s_pp: Vec<Vec<Complex64>>,
}

/// Frequency-swept full multi-port S-parameter matrix
/// (Phase 4.fem.eig.3 step F5 output of
/// [`OpenBoundarySolver::sweep_matrix`]).
///
/// Carries the **complete** `n_ports × n_ports` scattering matrix per
/// swept frequency. Entry `s[k][(q, p)]` is `S_{q,p}(omegas[k])` — the
/// modal amplitude received at port `q` when port `p` is driven with
/// `a_inc_p = 1` and every other port is matched (`a_inc_q = 0` for
/// `q ≠ p`).
///
/// # Layout
///
/// `s.len() == omegas.len()` and every `s[k]` is an
/// `(n_ports × n_ports)` complex dense matrix. Indexing follows the
/// nalgebra convention `(row, col) = (q, p)`: rows index the receive
/// port, columns index the excited port.
///
/// # Per-frequency cost model
///
/// Per spec §7, the driven matrix `A(ω)` is **independent** of which
/// port is excited (every wave-port face contributes its
/// `+ j β B_port` stiffness block unconditionally — only the RHS
/// carries the `a_inc_p` selection). The implementation therefore
/// factors `A(ω)` once per frequency and back-substitutes once per
/// excited port, giving an asymptotic per-frequency cost of
/// `O(LU(N) + n_ports · BS(N))` rather than the naive
/// `O(n_ports · LU(N))`.
///
/// # References
///
/// * Phase 4.fem.eig.3 spec
///   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
///   §4.3 (multi-port column-extraction convention) and §7 (LU-factor
///   reuse).
/// * Sheen, D. M., Ali, S. M., Abouzahra, M. D., Katehi, P. B. L.,
///   "Application of the three-dimensional finite-difference time-domain
///   method to the analysis of planar microstrip circuits",
///   *IEEE Trans. Microwave Theory Tech.* 38(7) (1990), pp. 849-857
///   — eq. 7 column extraction.
/// * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §4.3
///   — reciprocity `S_{p,q} = S_{q,p}` for lossless multi-ports.
#[derive(Debug, Clone)]
pub struct SParametersMatrix {
    /// Real-valued angular frequencies (rad/s) at which the sweep was
    /// evaluated; matches the order of the slice passed to
    /// [`OpenBoundarySolver::sweep_matrix`].
    pub omegas: Vec<f64>,
    /// Per-frequency `n_ports × n_ports` complex S-parameter matrix.
    /// `s[k][(q, p)]` is `S_{q,p}(omegas[k])` — response at port `q`
    /// driven by port `p`. Length of the outer vector equals
    /// `omegas.len()`; every inner matrix has shape
    /// `(n_ports, n_ports)`.
    pub s: Vec<DMatrix<Complex64>>,
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
    /// If `true`, scatter wave-port faces and project the FEM solution
    /// using the **exact Whitney-1 basis** at 3-point Gauss quadrature
    /// (Phase 4.fem.eig.3 F1+F2). Default `false` reproduces the v2 +
    /// CCCCCCCCC lumped-centroid behaviour bit-for-bit. Toggled by
    /// [`Self::with_coupled_whitney`].
    coupled_whitney: bool,
    /// Selects the Engquist–Majda ABC bilinear form on
    /// [`FaceKind::Abc`]-tagged faces (Phase 4.fem.eig.3 F4). Default
    /// [`AbcOrder::First`] reproduces the v2 1st-order Mur behaviour
    /// bit-for-bit. Toggled by [`Self::with_abc_order`].
    abc_order: AbcOrder,
    /// CFS-PML per-tet classification array (Phase 4.fem.eig.3.5 P4).
    /// `None` for the v3 default surface-integral path; `Some` after
    /// [`Self::with_cfs_pml`] has been invoked. When populated, every
    /// non-[`crate::PmlClass::Interior`] tet contributes its mass /
    /// stiffness block through the anisotropic-ε per-tet helper
    /// ([`crate::assemble_tet_element_complex_anisotropic`]) with the
    /// stretched-coordinate `Λ(ω)` factor; interior tets use the
    /// scalar [`crate::assemble_tet_element_complex`] bit-for-bit
    /// identical to the v3 path.
    pml_classes: Option<Vec<crate::PmlClass>>,
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
            coupled_whitney: false,
            abc_order: AbcOrder::First,
            pml_classes: None,
        })
    }

    /// Toggle the **coupled exact-Whitney-1** wave-port path
    /// (Phase 4.fem.eig.3 F1+F2).
    ///
    /// When `coupled = false` (the default), the wave-port scatter path
    /// uses the v2 lumped face-centroid quadrature
    /// ([`crate::element::assemble_port_face_block`] +
    /// [`crate::element::assemble_port_modal_rhs`]) and
    /// [`Self::extract_s11`]'s `E_FEM`-reconstruction uses the lumped
    /// `t_i / 3` proxy — bit-for-bit identical to the v2 + CCCCCCCCC
    /// shipped behaviour.
    ///
    /// When `coupled = true`, both the modal RHS and the FEM-side
    /// projection are lifted to the **exact Whitney-1 identity**
    /// `N_i(ξ) = λ_a(ξ) ∇λ_b − λ_b(ξ) ∇λ_a` evaluated at the three
    /// Gauss points
    ///
    /// ```text
    ///     ξ_g ∈ { (2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 1/6, 2/3) }
    /// ```
    ///
    /// on the reference triangle (each weighted `A / 3`). The two
    /// paths are changed together so the modal round-trip
    /// cancellation that Pozar §3.3 / Jin §10.5 derives is preserved
    /// at the exact-basis level, not the lumped level.
    ///
    /// The stiffness face block computed via the F1 entry point
    /// [`crate::element::assemble_port_face_block_gauss_pts`] is also
    /// substituted for the v2 lumped block; for a planar face the two
    /// stiffness blocks agree numerically (the Gauss-rule sum is
    /// degree-2 exact and the integrand `(n̂ × N_i) · (n̂ × N_j)` is
    /// linear × linear = degree 2), so the coupled stiffness path
    /// produces the same matrix to round-off — only the RHS and the
    /// `extract_s11` projection actually differ between the two paths.
    pub fn with_coupled_whitney(mut self, coupled: bool) -> Self {
        self.coupled_whitney = coupled;
        self
    }

    /// Read-only borrow of the `coupled_whitney` flag.
    pub fn coupled_whitney(&self) -> bool {
        self.coupled_whitney
    }

    /// Set the Engquist–Majda ABC bilinear-form order on
    /// [`FaceKind::Abc`]-tagged faces (Phase 4.fem.eig.3 F4).
    ///
    /// Default [`AbcOrder::First`] reproduces the v2 1st-order Mur
    /// behaviour bit-for-bit — every ABC face contributes
    /// [`crate::element::assemble_abc_face_block`] unchanged. Selecting
    /// [`AbcOrder::Second`] augments each ABC face's stiffness block
    /// with the tangential-curl correction term from Engquist–Majda
    /// 1979 eq. 9 via [`crate::element::assemble_abc2_face_block`]; the
    /// reflection floor for a TE plane wave at normal incidence drops
    /// from `~ −40 dB` to `~ −60 dB` (Jin §10.4, Table 10.1).
    ///
    /// The 2nd-order correction has a real scalar prefactor while the
    /// 1st-order part stays purely imaginary, so an `AbcOrder::Second`
    /// face contributes both Re and Im entries to the driven matrix on
    /// the face-edge rows/columns. The composite block is still
    /// complex-symmetric (`B == B^T`); `faer::sparse::Lu<usize,
    /// Complex64>` handles it unchanged from v2.
    pub fn with_abc_order(mut self, order: AbcOrder) -> Self {
        self.abc_order = order;
        self
    }

    /// Read-only borrow of the `abc_order` configuration.
    pub fn abc_order(&self) -> AbcOrder {
        self.abc_order
    }

    /// Configure the open-boundary solver for **CFS-PML** volumetric
    /// truncation (Phase 4.fem.eig.3.5 P4; Roden–Gedney 2000).
    ///
    /// The caller has already extended the mesh with PML brick shells
    /// via [`crate::extend_mesh_with_pml`]; this builder accepts the
    /// resulting per-tet [`crate::PmlClass`] classification (one entry
    /// per tet in the extended mesh, length must match
    /// `self.mesh.tetrahedra.len()`) and the [`PmlConfig`] carrying
    /// the grading parameters.
    ///
    /// Effect on assembly:
    ///
    /// * Sets `abc_order` to [`AbcOrder::CfsPml(config)`] so the
    ///   surface-integral Engquist–Majda kernel is **not** applied to
    ///   ABC-tagged faces — the volumetric PML absorbs in the bulk.
    ///   ABC-tagged faces become effectively "transparent" (their
    ///   surface integral is zero in the spec §3.2 limit `Λ(d=0) = I`).
    /// * Every per-tet stiffness + mass block is now computed via the
    ///   anisotropic helper [`crate::assemble_tet_element_complex_anisotropic`]
    ///   with a diagonal `Λ(ω)` factor. Interior tets get
    ///   `Λ = I` and the result matches the v3 scalar path bit-for-bit.
    /// * For PML tets, the diagonal `Λ(ω)` follows the
    ///   stretched-coordinate identity per spec §3.1:
    ///   `Λ = diag(s_y s_z / s_x, s_z s_x / s_y, s_x s_y / s_z)`,
    ///   with `s_α(ω) = κ_α(d_α) + σ_α(d_α) / (α_α + j ω ε_0)`. The
    ///   per-axis depths `d_α` are read from the [`crate::PmlClass`]
    ///   variant payload (always non-negative; v3.5 emits at most one
    ///   non-zero axis per tet).
    ///
    /// `with_cfs_pml` is **mutually exclusive** with
    /// `with_abc_order(First | Second)` — calling both is fine
    /// (later call wins) but only `CfsPml` triggers the volumetric
    /// path; surface-integral ABC kernels are not applied alongside
    /// the PML shell.
    ///
    /// # Arguments
    ///
    /// * `config` — grading parameters. Use [`PmlConfig::default`]
    ///   plus [`PmlConfig::resolved`] to populate the sentinel
    ///   `sigma_max` / `alpha_max` from a band-centre frequency and
    ///   mean tet edge length. Roden–Gedney 2000 §III/IV defaults
    ///   apply otherwise.
    /// * `pml_classes` — per-tet classification (length = number of
    ///   tets in `self.mesh`).
    ///
    /// # Panics
    ///
    /// Panics if `pml_classes.len()` does not match the mesh's tet
    /// count. The panic is preferred over a `Result` here because the
    /// classes array is produced by [`crate::extend_mesh_with_pml`]
    /// against the same extended mesh; a length mismatch is a caller
    /// bug, not a runtime failure mode.
    pub fn with_cfs_pml(mut self, config: PmlConfig, pml_classes: Vec<crate::PmlClass>) -> Self {
        assert_eq!(
            pml_classes.len(),
            self.mesh.tetrahedra.len(),
            "with_cfs_pml: pml_classes length {} does not match mesh tet count {}",
            pml_classes.len(),
            self.mesh.tetrahedra.len(),
        );
        self.abc_order = AbcOrder::CfsPml(config);
        self.pml_classes = Some(pml_classes);
        self
    }

    /// Read-only borrow of the CFS-PML per-tet classification array,
    /// `None` if the v3 surface-integral path is active.
    pub fn pml_classes(&self) -> Option<&[crate::PmlClass]> {
        self.pml_classes.as_deref()
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

    /// Geometric centroids of the exterior faces, in the same canonical
    /// order as [`Self::face_kinds`].
    ///
    /// Callers building face-kind tagging from geometric criteria (e.g.
    /// "tag the face at `z = 0` ABC, the face at `z = d` WavePort") can
    /// use this accessor to identify faces by centroid position **after
    /// constructing the solver** with a placeholder tagging — note that
    /// in practice the convention is to build the centroid list first
    /// (via a temporary all-PEC solver), classify by centroid position,
    /// then rebuild the solver with the proper tagging.
    pub fn exterior_face_centroids(&self) -> Vec<Vector3<f64>> {
        self.exterior_faces
            .faces
            .iter()
            .map(|f| f.centroid(self.mesh))
            .collect()
    }

    /// Outward unit normals of the exterior faces, in the same canonical
    /// order as [`Self::face_kinds`]. Peer of
    /// [`Self::exterior_face_centroids`].
    pub fn exterior_face_normals(&self) -> Vec<Vector3<f64>> {
        self.exterior_faces.faces.iter().map(|f| f.normal).collect()
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

        // Phase 4.fem.eig.3.5 P4 dispatch: if CFS-PML has been wired
        // via `with_cfs_pml`, the assembly path is the volumetric
        // anisotropic-ε path (per-tet `Λ(ω)` factor on the PML shell;
        // interior tets reduce to the scalar path bit-for-bit).
        // Surface-integral ABC scatter is skipped — the PML absorbs in
        // the bulk and the original Abc-tagged faces become smooth
        // material interfaces (continuity is preserved by the
        // polynomial grading `σ(d=0) = 0`).
        if let AbcOrder::CfsPml(config) = self.abc_order
            && self.pml_classes.is_some()
        {
            return self.assemble_driven_system_pml(omega, config);
        }
        // If CfsPml is set but pml_classes is not (caller forgot to
        // call with_cfs_pml), surface the configuration mismatch.
        if matches!(self.abc_order, AbcOrder::CfsPml(_)) {
            return Err(Error::Invalid(
                "OpenBoundarySolver::assemble_driven_system: \
                 abc_order = CfsPml but pml_classes is None — call \
                 `with_cfs_pml(config, classes)` to wire the PML path"
                    .to_string(),
            ));
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
                    if self.coupled_whitney {
                        // Phase 4.fem.eig.3 F1+F2 path: sample the
                        // modal profile at the three Gauss points on
                        // the reference triangle and call the exact-
                        // Whitney face-block + RHS helpers.
                        let face_vertices = face.world_vertices(self.mesh);
                        let mut e_t_gauss = [Vector3::<f64>::zeros(); 3];
                        for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
                            let p_g = bary[0] * face_vertices[0]
                                + bary[1] * face_vertices[1]
                                + bary[2] * face_vertices[2];
                            e_t_gauss[g] = (port.modal_e_t)(p_g);
                        }
                        self.scatter_port_face_gauss(
                            face,
                            beta,
                            e_t_gauss,
                            &interior_dof_of_edge,
                            &mut triplets,
                            &mut rhs,
                        );
                    } else {
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

    /// Frequency-sweep driven solve with diagonal S-parameter extraction
    /// (Phase 4.fem.eig.2 step E4).
    ///
    /// For each `ω` in `omegas`:
    ///
    /// 1. Call [`Self::solve_at_frequency`] to obtain the interior-DoF
    ///    complex solution vector `e_interior(ω)`.
    /// 2. For each wave-port `p`: project `e_interior` onto port `p`'s
    ///    modal profile via face-centroid quadrature (see module-level
    ///    docs for the formula), giving the modal reflection amplitude
    ///    `b_p(ω) = 2 ⟨E_FEM,t, e_mode_p⟩_port − a_inc_p`.
    /// 3. With normalised incident amplitude `a_inc_p = 1`,
    ///    `S_{p,p}(ω) = b_p(ω)`.
    ///
    /// The returned [`SParameters`] carries one diagonal sweep vector
    /// per port: `s_pp[p]` is a `Vec<Complex64>` of length `omegas.len()`
    /// indexing the per-frequency `S_{p,p}` value. Cross-port `S_{p,q}`
    /// for `p ≠ q` lands in Phase 4.fem.eig.2.0.2 per spec §13.
    ///
    /// # Arguments
    ///
    /// * `omegas` — non-empty slice of real-valued angular frequencies
    ///   (rad/s). Every entry must be positive; below-cutoff frequencies
    ///   are passed through to [`Self::solve_at_frequency`] verbatim, and
    ///   the wave-port modal contribution is whatever
    ///   [`PortDefinition::beta_mode`] returns (the caller is
    ///   responsible for the below-cutoff branch).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if `omegas` is empty. Propagates any
    /// error from [`Self::solve_at_frequency`] (mesh / material shape
    /// mismatch, sparse LU failure) verbatim.
    pub fn sweep(&self, omegas: &[f64]) -> Result<SParameters, Error> {
        if omegas.is_empty() {
            return Err(Error::Invalid(
                "OpenBoundarySolver::sweep: omegas slice is empty".to_string(),
            ));
        }

        let n_ports = self.ports.len();
        let mut s_pp: Vec<Vec<Complex64>> = vec![Vec::with_capacity(omegas.len()); n_ports];

        for &omega in omegas {
            // Re-assemble the driven system so we have access to the
            // interior-edge lift map needed for modal projection.
            let system = self.assemble_driven_system(omega)?;
            let n_interior = system.rhs.len();

            // Sparse LU factor + back-substitute, mirroring
            // solve_at_frequency.
            let lu: Lu<usize, Complex64> = system.matrix.sp_lu().map_err(|e| {
                Error::Numerical(format!(
                    "OpenBoundarySolver::sweep: sparse LU of driven matrix at \
                     omega = {omega} failed: {e:?}"
                ))
            })?;

            let mut rhs_mat = faer::Mat::<Complex64>::zeros(n_interior, 1);
            for (i, &b_i) in system.rhs.iter().enumerate() {
                rhs_mat[(i, 0)] = b_i;
            }
            lu.solve_in_place_with_conj(faer::Conj::No, rhs_mat.as_mut());

            let e_interior: Vec<Complex64> = (0..n_interior).map(|i| rhs_mat[(i, 0)]).collect();

            // Extract S_{p,p}(ω) for every port.
            for (p, s_vec) in s_pp.iter_mut().enumerate().take(n_ports) {
                let s_pp_omega = self.extract_s11(p, omega, &e_interior, &system)?;
                s_vec.push(s_pp_omega);
            }
        }

        Ok(SParameters {
            omegas: omegas.to_vec(),
            s_pp,
        })
    }

    /// Frequency-sweep driven solve returning the full multi-port
    /// `S_{p,q}` matrix (Phase 4.fem.eig.3 step F5).
    ///
    /// For each `ω` in `omegas`:
    ///
    /// 1. Assemble the driven system `A(ω) e = b(ω)` once via
    ///    [`Self::assemble_driven_system`] — the matrix is independent
    ///    of which port is excited, so we factor it via
    ///    `faer::sparse::Lu<usize, Complex64>` exactly once per
    ///    frequency.
    /// 2. For each excited port `p ∈ 0..n_ports`: build a port-specific
    ///    RHS in which **only** port `p`'s modal contribution is
    ///    included (`a_inc_p = 1`, `a_inc_q = 0` for `q ≠ p`); the
    ///    other ports' face stiffness blocks remain in the matrix so
    ///    the matched-port condition is enforced naturally by the wave-
    ///    port bilinear form (Pozar §3.3 / Jin §10.5). Back-substitute
    ///    against the cached LU factor to obtain
    ///    `e_interior(ω; driven by p)`.
    /// 3. For each receive port `q ∈ 0..n_ports`: extract
    ///    `S_{q,p}(ω) = ⟨E_FEM, e_mode_q⟩_port / M_qq − a_inc_q`, where
    ///    `a_inc_q = δ_{q,p}` (i.e. `1` on the diagonal, `0`
    ///    off-diagonal). Pack into the entry `s[k][(q, p)]` of the
    ///    output matrix.
    ///
    /// # Arguments
    ///
    /// * `omegas` — non-empty slice of real-valued angular frequencies
    ///   (rad/s). Every entry must be positive; below-cutoff behaviour
    ///   is governed by each port's
    ///   [`PortDefinition::beta_mode`] closure (same convention as
    ///   [`Self::sweep`]).
    ///
    /// # Returns
    ///
    /// [`SParametersMatrix`] carrying one `n_ports × n_ports` complex
    /// dense matrix per swept frequency. Entry `(q, p)` follows the
    /// Sheen–Ali–Abouzahra–Katehi 1990 column-extraction convention —
    /// `S_{q,p}` is the modal amplitude at port `q` when port `p` is
    /// driven.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if `omegas` is empty. Propagates any
    /// error from [`Self::assemble_driven_system`] (mesh / material
    /// shape mismatch) and surfaces an [`Error::Numerical`] variant if
    /// the sparse LU of the driven matrix fails at any frequency or if
    /// any port's modal self-inner-product is numerically zero.
    pub fn sweep_matrix(&self, omegas: &[f64]) -> Result<SParametersMatrix, Error> {
        if omegas.is_empty() {
            return Err(Error::Invalid(
                "OpenBoundarySolver::sweep_matrix: omegas slice is empty".to_string(),
            ));
        }

        let n_ports = self.ports.len();
        let mut s_out: Vec<DMatrix<Complex64>> = Vec::with_capacity(omegas.len());

        for &omega in omegas {
            // Assemble the driven system once at this frequency. The
            // matrix carries every port's `+ j β B_port` stiffness
            // block regardless of which port is excited, so the same
            // factor handles every excited-port RHS in the inner loop.
            // We discard `system.rhs` because it bundles every port's
            // modal-current contribution simultaneously; per spec §4.3
            // the multi-port extraction needs the RHS to isolate one
            // excited port at a time.
            let system = self.assemble_driven_system(omega)?;
            let n_interior = system.rhs.len();

            let lu: Lu<usize, Complex64> = system.matrix.sp_lu().map_err(|e| {
                Error::Numerical(format!(
                    "OpenBoundarySolver::sweep_matrix: sparse LU of driven matrix at \
                     omega = {omega} failed: {e:?}"
                ))
            })?;

            // Allocate the n_ports × n_ports output for this frequency.
            let mut s_k = DMatrix::<Complex64>::zeros(n_ports, n_ports);

            for p in 0..n_ports {
                // Build the per-excited-port RHS. Only port `p`
                // contributes its modal-current term (a_inc_p = 1);
                // all other ports contribute zero RHS (a_inc_q = 0)
                // because the modal-current scatter is linear in
                // `e_t` and `e_t = a_inc · e_mode`.
                let rhs_p = self.build_rhs_for_excited_port(
                    omega,
                    p,
                    &system.interior_dof_of_edge,
                    n_interior,
                )?;

                // Back-substitute against the cached LU factor.
                let mut rhs_mat = faer::Mat::<Complex64>::zeros(n_interior, 1);
                for (i, &b_i) in rhs_p.iter().enumerate() {
                    rhs_mat[(i, 0)] = b_i;
                }
                lu.solve_in_place_with_conj(faer::Conj::No, rhs_mat.as_mut());

                let e_interior: Vec<Complex64> = (0..n_interior).map(|i| rhs_mat[(i, 0)]).collect();

                // Extract S_{q, p} for every receive port q. The
                // `extract_s_qp` helper subtracts `a_inc_q = δ_{q,p}`
                // (1 on the diagonal, 0 off-diagonal) per the
                // Sheen et al. 1990 convention.
                for q in 0..n_ports {
                    let a_inc_q = if q == p { 1.0 } else { 0.0 };
                    let s_qp = self.extract_s_qp(q, a_inc_q, &e_interior, &system)?;
                    s_k[(q, p)] = s_qp;
                }
            }

            s_out.push(s_k);
        }

        Ok(SParametersMatrix {
            omegas: omegas.to_vec(),
            s: s_out,
        })
    }

    /// Build a port-specific RHS vector for the multi-port sweep
    /// (Phase 4.fem.eig.3 F5 helper).
    ///
    /// Returns a fresh RHS in which **only** wave-port faces tagged
    /// `WavePort(excited_port)` contribute their modal-current scatter;
    /// every other port's RHS contribution is zero (`a_inc_q = 0` for
    /// `q ≠ excited_port`). The matrix stays unchanged because the
    /// wave-port bilinear form `+ j β B_port` is intrinsic to the
    /// boundary condition — every wave-port face contributes its
    /// stiffness block regardless of whether it is driven or matched.
    ///
    /// This is a thin wrapper around the same scatter helpers used by
    /// [`Self::assemble_driven_system`] (`scatter_port_face` for the
    /// v2 lumped path; `scatter_port_face_gauss` for the
    /// coupled-Whitney path), with the matrix-side accumulation
    /// discarded. PEC-precedence handling and per-edge orientation
    /// signs are inherited unchanged from the v2 scatter path.
    fn build_rhs_for_excited_port(
        &self,
        omega: f64,
        excited_port: PortId,
        interior_dof_of_edge: &[Option<usize>],
        n_interior: usize,
    ) -> Result<Vec<Complex64>, Error> {
        if excited_port >= self.ports.len() {
            return Err(Error::Invalid(format!(
                "OpenBoundarySolver::build_rhs_for_excited_port: \
                 excited_port = {excited_port} out of range (n_ports = {})",
                self.ports.len()
            )));
        }

        let mut rhs: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_interior];
        // Scratch matrix-side triplets — discarded after the helper
        // returns. The scatter helpers accept a `&mut Vec<Triplet>`
        // unconditionally; we feed them a local sink so the RHS path
        // is exercised without polluting the (already-factored) matrix.
        let mut sink: Vec<Triplet<usize, usize, Complex64>> = Vec::new();

        for (i, kind) in self.face_kinds.iter().enumerate() {
            if let FaceKind::WavePort(p) = *kind
                && p == excited_port
            {
                let face = &self.exterior_faces.faces[i];
                let port = &self.ports[p];
                let beta = (port.beta_mode)(omega);
                if self.coupled_whitney {
                    let face_vertices = face.world_vertices(self.mesh);
                    let mut e_t_gauss = [Vector3::<f64>::zeros(); 3];
                    for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
                        let p_g = bary[0] * face_vertices[0]
                            + bary[1] * face_vertices[1]
                            + bary[2] * face_vertices[2];
                        e_t_gauss[g] = (port.modal_e_t)(p_g);
                    }
                    self.scatter_port_face_gauss(
                        face,
                        beta,
                        e_t_gauss,
                        interior_dof_of_edge,
                        &mut sink,
                        &mut rhs,
                    );
                } else {
                    let centroid = face.centroid(self.mesh);
                    let e_t = (port.modal_e_t)(centroid);
                    self.scatter_port_face(
                        face,
                        beta,
                        e_t,
                        interior_dof_of_edge,
                        &mut sink,
                        &mut rhs,
                    );
                }
            }
        }

        Ok(rhs)
    }

    /// Extract `S_{q,p}(ω)` for a generic receive port `q` with
    /// caller-supplied `a_inc_q` (Phase 4.fem.eig.3 F5 helper).
    ///
    /// Implements the Sheen–Ali–Abouzahra–Katehi 1990 eq. 7 column
    /// extraction
    ///
    /// ```text
    ///     b_q  =  ⟨ E_FEM , e_mode_q ⟩_port / M_qq  −  a_inc_q,
    ///     S_{q,p}  =  b_q / a_inc_p   =   b_q          (a_inc_p = 1),
    /// ```
    ///
    /// where `a_inc_q = δ_{q,p}` is the incident amplitude at the
    /// receive port (`1` on the diagonal entry, `0` off-diagonal).
    /// The modal-projection inner product and `M_qq` self-inner-
    /// product use the same face-centroid (or three-point Gauss-
    /// quadrature, under [`Self::coupled_whitney`]) integration as
    /// [`Self::extract_s11`].
    ///
    /// This is the multi-port generalisation of [`Self::extract_s11`]:
    /// the single-port case `q = p, a_inc_q = 1` reproduces
    /// `extract_s11(p, ...)` bit-for-bit.
    ///
    /// # Arguments
    ///
    /// * `port_id` — receive port `q` (matched against
    ///   `FaceKind::WavePort(q)` tags).
    /// * `a_inc_q` — incident amplitude at port `q`. Conventionally
    ///   `1.0` if `q` is the driven port, `0.0` otherwise.
    /// * `e_interior` — interior-DoF complex solution from the per-
    ///   excited-port back-substitution.
    /// * `system` — driven system returned by
    ///   [`Self::assemble_driven_system`] (provides the lift map).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if `port_id` is out of range.
    /// Returns [`Error::Numerical`] if the modal self-inner-product
    /// `M_qq` is numerically zero.
    fn extract_s_qp(
        &self,
        port_id: PortId,
        a_inc_q: f64,
        e_interior: &[Complex64],
        system: &DrivenSystem,
    ) -> Result<Complex64, Error> {
        if port_id >= self.ports.len() {
            return Err(Error::Invalid(format!(
                "OpenBoundarySolver::extract_s_qp: port_id = {port_id} \
                 out of range (n_ports = {})",
                self.ports.len()
            )));
        }

        let port = &self.ports[port_id];
        let mut inner_product = Complex64::new(0.0, 0.0);
        let mut mode_self_inner = 0.0_f64;

        for (i, kind) in self.face_kinds.iter().enumerate() {
            if let FaceKind::WavePort(p) = *kind
                && p == port_id
            {
                let face = &self.exterior_faces.faces[i];
                let face_vertices = face.world_vertices(self.mesh);

                let t0 = face_vertices[1] - face_vertices[0];
                let t1 = face_vertices[2] - face_vertices[1];
                let face_area = 0.5 * t0.cross(&t1).norm();

                if self.coupled_whitney {
                    let e_fem_g =
                        self.e_t_at_face_gauss_pts(face, e_interior, &system.interior_dof_of_edge);
                    let w_g = face_area / 3.0;
                    for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
                        let p_g = bary[0] * face_vertices[0]
                            + bary[1] * face_vertices[1]
                            + bary[2] * face_vertices[2];
                        let e_mode_g = (port.modal_e_t)(p_g);
                        let e_fem = e_fem_g[g];
                        let dot_g = e_fem.x * Complex64::new(e_mode_g.x, 0.0)
                            + e_fem.y * Complex64::new(e_mode_g.y, 0.0)
                            + e_fem.z * Complex64::new(e_mode_g.z, 0.0);
                        inner_product += Complex64::new(w_g, 0.0) * dot_g;
                        mode_self_inner += w_g * e_mode_g.dot(&e_mode_g);
                    }
                } else {
                    let centroid = face.centroid(self.mesh);
                    let e_mode = (port.modal_e_t)(centroid);
                    let e_fem =
                        self.e_t_at_face_centroid(face, e_interior, &system.interior_dof_of_edge);
                    let face_dot = e_fem.x * Complex64::new(e_mode.x, 0.0)
                        + e_fem.y * Complex64::new(e_mode.y, 0.0)
                        + e_fem.z * Complex64::new(e_mode.z, 0.0);
                    inner_product += Complex64::new(face_area, 0.0) * face_dot;
                    mode_self_inner += face_area * e_mode.dot(&e_mode);
                }
            }
        }

        if mode_self_inner <= f64::EPSILON {
            return Err(Error::Numerical(format!(
                "OpenBoundarySolver::extract_s_qp: modal self-inner-product \
                 ⟨e_mode, e_mode⟩_port = {mode_self_inner} is numerically \
                 zero for port {port_id}; cannot normalise extraction"
            )));
        }

        let m_qq = Complex64::new(mode_self_inner, 0.0);
        let a_inc_c = Complex64::new(a_inc_q, 0.0);
        Ok(inner_product / m_qq - a_inc_c)
    }

    /// Extract `S_{p,p}(ω)` for a single port from an interior-DoF
    /// complex solution vector (Phase 4.fem.eig.2 step E4).
    ///
    /// Implements the modal projection
    ///
    /// ```text
    ///     ⟨ E_FEM,t , e_mode_p ⟩_port
    ///         =  Σ_face  A_face · ( E_FEM,t(centroid) · e_mode_p(centroid) ),
    /// ```
    ///
    /// summed over every exterior face tagged
    /// [`FaceKind::WavePort`]`(p)`. `E_FEM,t(centroid)` is reconstructed
    /// from the per-edge interior DoFs via the Whitney-1 face basis
    /// evaluated at the centroid (see the private
    /// `e_t_at_face_centroid` helper for the closed-form basis-at-
    /// centroid weighting).
    ///
    /// ## Modal-normalisation correction (CCCCCCCCC)
    ///
    /// The original Phase 4.fem.eig.2 E4 formula
    /// `b_p = 2 ⟨E_FEM, e_mode⟩ − a_inc` (spec §4.3) implicitly assumed
    /// `⟨e_mode, e_mode⟩_port = 1/2`, but the [`PortDefinition::modal_e_t`]
    /// contract in this crate (and the WR-90 TE_{10} profile used by
    /// `fem-eig-003` / `crates/yee-validation`) carries the standard
    /// orthonormalisation `⟨e_mode, e_mode⟩_port = 1` (Pozar §3.3).
    /// With the un-normalised formula, even a matched-port total field
    /// `E_FEM = a_inc · e_mode` yielded `|S_{11}| = 1` instead of the
    /// expected `0`, causing the fem-eig-003 sweep to saturate at
    /// `|S_{11}| = 1.0` (BBBBBBBBB E5 finding).
    ///
    /// The corrected extraction normalises the projection by the modal
    /// self-inner-product computed via the same face-centroid quadrature:
    ///
    /// ```text
    ///     M_pp   =   ⟨ e_mode_p , e_mode_p ⟩_port,
    ///     b_p    =   ⟨ E_FEM,t , e_mode_p ⟩_port / M_pp  −  a_inc,
    ///     S_{p,p} =  b_p / a_inc.
    /// ```
    ///
    /// With the standard `a_inc = 1` and a matched-port total field
    /// `E_FEM ≈ a_inc · e_mode`, the corrected formula gives
    /// `b_p ≈ a_inc − a_inc = 0` — the expected matched-port identity.
    /// With a fully-reflective PEC termination behind the port and no
    /// modal amplitude landing at the port (`E_FEM ≈ 0`),
    /// `b_p ≈ 0 − a_inc = −a_inc`, recovering `|S_{11}| ≈ 1` (full
    /// reflection). Both end-cases match Pozar §3.3 and Jin §10.7.
    ///
    /// `M_pp` is computed once per call; for a degenerate fixture with
    /// `M_pp ≈ 0` (no port faces or modal profile zero everywhere) the
    /// function returns [`Error::Numerical`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if `port_id` is out of range.
    /// Returns [`Error::Numerical`] if the modal self-inner-product
    /// `M_pp` is numerically zero (which would otherwise produce a
    /// division-by-zero `Inf` in the extracted `S_{p,p}`).
    pub fn extract_s11(
        &self,
        port_id: PortId,
        _omega: f64,
        e_interior: &[Complex64],
        system: &DrivenSystem,
    ) -> Result<Complex64, Error> {
        if port_id >= self.ports.len() {
            return Err(Error::Invalid(format!(
                "OpenBoundarySolver::extract_s11: port_id = {port_id} \
                 out of range (n_ports = {})",
                self.ports.len()
            )));
        }

        let port = &self.ports[port_id];
        let mut inner_product = Complex64::new(0.0, 0.0);
        let mut mode_self_inner = 0.0_f64;

        for (i, kind) in self.face_kinds.iter().enumerate() {
            if let FaceKind::WavePort(p) = *kind
                && p == port_id
            {
                let face = &self.exterior_faces.faces[i];
                let face_vertices = face.world_vertices(self.mesh);

                // Face area A_face = 0.5 · ||t_0 × t_1||.
                let t0 = face_vertices[1] - face_vertices[0];
                let t1 = face_vertices[2] - face_vertices[1];
                let face_area = 0.5 * t0.cross(&t1).norm();

                if self.coupled_whitney {
                    // Phase 4.fem.eig.3 F1+F2 path: project E_FEM and
                    // the modal profile at the same three Gauss points
                    // used by `assemble_port_face_rhs_gauss_pts`. The
                    // round-trip cancellation now holds at the exact
                    // Whitney-1 basis level (Pozar §3.3 matched-port
                    // identity), not just the CCCCCCCCC M_pp level.
                    let e_fem_g =
                        self.e_t_at_face_gauss_pts(face, e_interior, &system.interior_dof_of_edge);
                    let w_g = face_area / 3.0;
                    for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
                        let p_g = bary[0] * face_vertices[0]
                            + bary[1] * face_vertices[1]
                            + bary[2] * face_vertices[2];
                        let e_mode_g = (port.modal_e_t)(p_g);
                        let e_fem = e_fem_g[g];
                        let dot_g = e_fem.x * Complex64::new(e_mode_g.x, 0.0)
                            + e_fem.y * Complex64::new(e_mode_g.y, 0.0)
                            + e_fem.z * Complex64::new(e_mode_g.z, 0.0);
                        inner_product += Complex64::new(w_g, 0.0) * dot_g;
                        mode_self_inner += w_g * e_mode_g.dot(&e_mode_g);
                    }
                } else {
                    // v2 + CCCCCCCCC lumped-centroid path. Bit-for-bit
                    // unchanged from the Phase 4.fem.eig.2 shipped
                    // behaviour.
                    let centroid = face.centroid(self.mesh);
                    let e_mode = (port.modal_e_t)(centroid);

                    let e_fem =
                        self.e_t_at_face_centroid(face, e_interior, &system.interior_dof_of_edge);

                    let face_dot = e_fem.x * Complex64::new(e_mode.x, 0.0)
                        + e_fem.y * Complex64::new(e_mode.y, 0.0)
                        + e_fem.z * Complex64::new(e_mode.z, 0.0);
                    inner_product += Complex64::new(face_area, 0.0) * face_dot;
                    mode_self_inner += face_area * e_mode.dot(&e_mode);
                }
            }
        }

        if mode_self_inner <= f64::EPSILON {
            return Err(Error::Numerical(format!(
                "OpenBoundarySolver::extract_s11: modal self-inner-product \
                 ⟨e_mode, e_mode⟩_port = {mode_self_inner} is numerically \
                 zero for port {port_id}; cannot normalise S_{{11}} \
                 extraction (modal profile vanishes on every port face?)"
            )));
        }

        // S_{p,p} = b_p / a_inc with a_inc = 1 and
        // b_p = ⟨E_FEM, e_mode⟩_port / M_pp − a_inc (CCCCCCCCC).
        let a_inc = Complex64::new(1.0, 0.0);
        let m_pp = Complex64::new(mode_self_inner, 0.0);
        let b_p = inner_product / m_pp - a_inc;
        Ok(b_p / a_inc)
    }

    /// Reconstruct the tangential `E`-field at a port face's centroid
    /// from the global interior-DoF complex solution vector
    /// (Phase 4.fem.eig.2 step E4 helper).
    ///
    /// For the Whitney-1 face basis the per-edge basis function `N_i`
    /// evaluated at the centroid is treated as the lumped edge-tangent
    /// proxy `t_i / 3`, where `t_i = v_{(i+1) mod 3} − v_i` is the
    /// canonical face-edge tangent. The face-centroid FEM E-field is
    /// therefore
    ///
    /// ```text
    ///     E_FEM,t(centroid)  =  Σ_{i ∈ face_edges}  s_i · e_i · (t_i / 3),
    /// ```
    ///
    /// where `s_i ∈ {-1, +1}` is the local-to-global orientation sign
    /// and `e_i` is the interior-DoF amplitude (or `0` if edge `i` is
    /// PEC-eliminated).
    ///
    /// ## CCCCCCCCC scaling note
    ///
    /// The `t_i / 3` lumped weighting is **not** the exact Whitney-1
    /// basis-at-centroid identity
    /// `N_i(centroid) = (1/3)(∇λ_b − ∇λ_a)`; in general
    /// `(∇λ_b − ∇λ_a) ≠ t_i`. The lumped form is retained here to
    /// match the dual approximation already in
    /// [`crate::element::assemble_port_modal_rhs`], so the round-trip
    /// modal-RHS-then-modal-projection cancellation that the spec
    /// §4.3 derivation relies on is preserved at the lumped level. The
    /// CCCCCCCCC scaling fix lives in [`Self::extract_s11`], which
    /// divides the inner product by the modal self-inner-product `M_pp`
    /// computed via the same lumped quadrature; that ratio is what
    /// retires the `|S_{11}| = 1` saturation. A future Phase
    /// 4.fem.eig.2.0.1 refinement (ADR-0040 §C-3) will lift both the
    /// element-layer RHS and this reconstruction to the exact Whitney
    /// basis identity in a single coupled change — independently
    /// validated against the cross-section eigensolver per-Gauss-point
    /// modal sampling.
    fn e_t_at_face_centroid(
        &self,
        face: &ExteriorFace,
        e_interior: &[Complex64],
        interior_dof_of_edge: &[Option<usize>],
    ) -> nalgebra::Vector3<Complex64> {
        let face_vertices = face.world_vertices(self.mesh);
        let t = [
            face_vertices[1] - face_vertices[0],
            face_vertices[2] - face_vertices[1],
            face_vertices[0] - face_vertices[2],
        ];

        let mut e_t = nalgebra::Vector3::<Complex64>::zeros();
        for (i, t_i) in t.iter().enumerate() {
            let gi = face.global_edges[i];
            let Some(dof) = interior_dof_of_edge[gi] else {
                continue;
            };
            let coeff = e_interior[dof];
            let sign = Complex64::new(face.signs[i], 0.0);
            let weight = sign * coeff * Complex64::new(1.0 / 3.0, 0.0);
            e_t.x += weight * Complex64::new(t_i.x, 0.0);
            e_t.y += weight * Complex64::new(t_i.y, 0.0);
            e_t.z += weight * Complex64::new(t_i.z, 0.0);
        }
        e_t
    }

    /// Scatter the per-face ABC `+ j k₀ B_ABC` block into the driven
    /// triplet list. PEC edges on the face are silently skipped — they
    /// are eliminated by the global row/column drop applied to `K(ω)`
    /// and `M(ω)` by [`FemEigenAssembly::assemble_complex`].
    ///
    /// Branches on [`Self::abc_order`]: [`AbcOrder::First`] (default)
    /// calls [`crate::element::assemble_abc_face_block`] for v2
    /// bit-for-bit behaviour; [`AbcOrder::Second`] calls
    /// [`crate::element::assemble_abc2_face_block`] which adds the
    /// Engquist–Majda 1979 eq. 9 tangential-curl correction term.
    fn scatter_abc_face(
        &self,
        face: &ExteriorFace,
        k0: f64,
        interior_dof_of_edge: &[Option<usize>],
        triplets: &mut Vec<Triplet<usize, usize, Complex64>>,
    ) {
        let block = match self.abc_order {
            AbcOrder::First => {
                assemble_abc_face_block(face.world_vertices(self.mesh), face.normal, k0, 1.0)
            }
            AbcOrder::Second => {
                assemble_abc2_face_block(face.world_vertices(self.mesh), face.normal, k0, 1.0)
            }
            AbcOrder::CfsPml(_) => {
                // Unreachable: `assemble_driven_system` returns
                // `Error::NotEnabled` before reaching this scatter when
                // `abc_order` is `CfsPml`. The volumetric PML path
                // wires in at Phase 4.fem.eig.3.5 P4.
                unreachable!(
                    "scatter_abc_face called with AbcOrder::CfsPml; \
                     guarded out by assemble_driven_system early return"
                );
            }
        };
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

    /// Coupled-Whitney variant of [`Self::scatter_port_face`]
    /// (Phase 4.fem.eig.3 F1+F2). Scatters the wave-port face block +
    /// RHS computed via the exact Whitney-1 basis at three Gauss
    /// points on the reference triangle. PEC-precedence and per-edge
    /// orientation-sign handling are identical to v2.
    fn scatter_port_face_gauss(
        &self,
        face: &ExteriorFace,
        beta: f64,
        e_t_gauss: [Vector3<f64>; 3],
        interior_dof_of_edge: &[Option<usize>],
        triplets: &mut Vec<Triplet<usize, usize, Complex64>>,
        rhs: &mut [Complex64],
    ) {
        let face_vertices = face.world_vertices(self.mesh);
        let beta_c = Complex64::new(beta, 0.0);

        // Stiffness contribution via the exact-Whitney Gauss-pt block.
        let block = assemble_port_face_block_gauss_pts(face_vertices, face.normal, beta_c, 1.0);
        for i in 0..3 {
            let gi = face.global_edges[i];
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

        // RHS contribution via the exact-Whitney Gauss-pt RHS helper.
        let rhs_block =
            assemble_port_face_rhs_gauss_pts(face_vertices, face.normal, beta_c, e_t_gauss);
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

    /// Reconstruct the tangential `E`-field at the three reference-
    /// triangle Gauss points of a port face from the global interior-
    /// DoF complex solution vector (Phase 4.fem.eig.3 F2 helper).
    ///
    /// For each Gauss point `ξ_g ∈ {(2/3, 1/6, 1/6), (1/6, 2/3, 1/6),
    /// (1/6, 1/6, 2/3)}` the FEM-side tangential E-field is the sum
    /// over the three face edges
    ///
    /// ```text
    ///     E_FEM,t(ξ_g)  =  Σ_{i ∈ face_edges}  s_i · e_i · N_i(ξ_g),
    /// ```
    ///
    /// where `s_i ∈ {-1, +1}` is the local-to-global orientation sign,
    /// `e_i` is the interior-DoF amplitude (or `0` if edge `i` is
    /// PEC-eliminated), and `N_i(ξ_g)` is the **exact** Whitney-1 edge
    /// basis at the Gauss point, computed from the in-plane
    /// barycentric gradients `∇λ_a, ∇λ_b, ∇λ_c` and the Whitney
    /// identity `N_i = λ_a ∇λ_b − λ_b ∇λ_a`.
    ///
    /// Pairs with [`Self::scatter_port_face_gauss`] in the coupled-
    /// Whitney path enabled by [`Self::with_coupled_whitney`]; the two
    /// helpers share the same Gauss-point set, the same per-face
    /// gradient construction, and the same Whitney-1 basis identity,
    /// preserving the modal-RHS-then-projection round-trip
    /// cancellation that Pozar §3.3 / Jin §10.5 derives at the
    /// exact-basis level.
    fn e_t_at_face_gauss_pts(
        &self,
        face: &ExteriorFace,
        e_interior: &[Complex64],
        interior_dof_of_edge: &[Option<usize>],
    ) -> [nalgebra::Vector3<Complex64>; 3] {
        let face_vertices = face.world_vertices(self.mesh);

        // In-plane barycentric gradients ∇λ_a (same identity as the
        // element-layer F1 helpers). For a triangle with vertices
        // (v_0, v_1, v_2) in CCW order seen from +n̂:
        //
        //     ∇λ_a = (v_b − v_c) × n̂ / (2 A),
        //
        // with (a, b, c) cyclic.
        let v0 = face_vertices[0];
        let v1 = face_vertices[1];
        let v2 = face_vertices[2];
        let face_area = 0.5 * (v1 - v0).cross(&(v2 - v0)).norm();

        let n_norm = face.normal.norm();
        let n_hat = if n_norm > 0.0 {
            face.normal / n_norm
        } else {
            face.normal
        };

        let inv_two_a = if face_area > 0.0 {
            1.0 / (2.0 * face_area)
        } else {
            0.0
        };
        let grad = [
            (v1 - v2).cross(&n_hat) * inv_two_a,
            (v2 - v0).cross(&n_hat) * inv_two_a,
            (v0 - v1).cross(&n_hat) * inv_two_a,
        ];

        // Per-edge weighted DoF amplitude s_i · e_i with PEC-eliminated
        // edges contributing zero.
        let mut edge_coeff = [Complex64::new(0.0, 0.0); 3];
        for (i, slot) in edge_coeff.iter_mut().enumerate() {
            let gi = face.global_edges[i];
            if let Some(dof) = interior_dof_of_edge[gi] {
                let coeff = e_interior[dof];
                let sign = Complex64::new(face.signs[i], 0.0);
                *slot = sign * coeff;
            }
        }

        let mut out = [nalgebra::Vector3::<Complex64>::zeros(); 3];
        for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
            // Exact Whitney-1 basis N_i(ξ_g) for i = 0, 1, 2 (edge i
            // runs from a = i to b = (i + 1) mod 3).
            let mut basis = [nalgebra::Vector3::<f64>::zeros(); 3];
            for (i, basis_i) in basis.iter_mut().enumerate() {
                let a = i;
                let b = (i + 1) % 3;
                *basis_i = bary[a] * grad[b] - bary[b] * grad[a];
            }

            let mut e_t = nalgebra::Vector3::<Complex64>::zeros();
            for (i, basis_i) in basis.iter().enumerate() {
                let c = edge_coeff[i];
                e_t.x += c * Complex64::new(basis_i.x, 0.0);
                e_t.y += c * Complex64::new(basis_i.y, 0.0);
                e_t.z += c * Complex64::new(basis_i.z, 0.0);
            }
            out[g] = e_t;
        }

        out
    }

    /// Assemble the CFS-PML driven open-boundary system at angular
    /// frequency `omega` (Phase 4.fem.eig.3.5 P4).
    ///
    /// Mirrors [`Self::assemble_driven_system`] but uses per-tet
    /// anisotropic-ε assembly via
    /// [`crate::assemble_tet_element_complex_anisotropic`] on PML
    /// tets (the `Λ(ω)` stretched-coordinate factor is computed
    /// per-axis per-tet) and the scalar
    /// [`crate::assemble_tet_element_complex`] on cavity-interior
    /// tets (bit-for-bit unchanged from v3). The boundary-term scatter
    /// for [`FaceKind::Abc`] is **skipped** — the PML absorbs in the
    /// bulk and the surface integral is identically zero in the
    /// `Λ(d=0) = I` continuity limit (spec §3.2).
    fn assemble_driven_system_pml(
        &self,
        omega: f64,
        config: PmlConfig,
    ) -> Result<DrivenSystem, Error> {
        let pml_classes = self.pml_classes.as_ref().ok_or_else(|| {
            Error::Invalid(
                "assemble_driven_system_pml: pml_classes is None — \
                 call `with_cfs_pml(config, classes)` first"
                    .to_string(),
            )
        })?;

        // Resolve sentinel config against the band-centre frequency
        // and a mean tet edge length read from the mesh.
        let freq_hz = omega / (2.0 * std::f64::consts::PI);
        let h_cell = self.mean_tet_edge_length();
        let cfg = config.resolved(freq_hz, h_cell);

        // ---- 1. Build the per-tet edge connectivity table (mirror of
        // assembly::TetEdgeTable, reproduced here to keep this method
        // self-contained inside the open_boundary lane). The same edge
        // numbering and orientation conventions are used as
        // FemEigenAssembly so the scalar path produces a matching
        // matrix structure on interior tets.
        let edge_table = PmlAssemblyEdgeTable::build(self.mesh);
        let n_edges = edge_table.edges.len();

        // ---- 2. Classify edges as PEC vs interior. PEC edges are
        // eliminated by row/column drop.
        let interior_edges: Vec<usize> = (0..n_edges)
            .filter(|e| !self.pec_global_edges.contains(e))
            .collect();
        let n_interior = interior_edges.len();
        let mut interior_dof_of_edge: Vec<Option<usize>> = vec![None; n_edges];
        for (dof, &gid) in interior_edges.iter().enumerate() {
            interior_dof_of_edge[gid] = Some(dof);
        }

        // ---- 3. Per-tet scatter into the driven matrix triplet list.
        // We accumulate `A = K - k0^2 M` directly (rather than separate
        // K and M COO matrices) because per-frequency assembly is the
        // only consumer.
        let k0 = omega / C0;
        let k0_sq = k0 * k0;
        let k0_sq_c = Complex64::new(k0_sq, 0.0);
        let mut triplets: Vec<Triplet<usize, usize, Complex64>> = Vec::new();

        for (tet_idx, conn) in edge_table.tet_edges.iter().enumerate() {
            let tet = &self.mesh.tetrahedra[tet_idx];
            let vertices = [
                self.mesh.vertices[tet[0]],
                self.mesh.vertices[tet[1]],
                self.mesh.vertices[tet[2]],
                self.mesh.vertices[tet[3]],
            ];
            let tag = self.mesh.tetrahedron_material[tet_idx];
            let eps_omega = self.material_db.eps_at(tag, omega);
            let mu_omega = self.material_db.mu_at(tag, omega);

            let class = pml_classes[tet_idx];
            let elem = if class.is_interior() {
                assemble_tet_element_complex(vertices, eps_omega, mu_omega)
            } else {
                let lam = pml_stretching_lambda(class, &cfg, omega);
                let mut eps_tensor = SMatrix::<Complex64, 3, 3>::zeros();
                let mut mu_inv_tensor = SMatrix::<Complex64, 3, 3>::zeros();
                let mu_inv = Complex64::new(1.0, 0.0) / mu_omega;
                for d in 0..3 {
                    eps_tensor[(d, d)] = eps_omega * lam[d];
                    // For Roden-Gedney 2000 §II the same Λ applies to
                    // both ε and μ; μ_inv therefore picks up `1 / Λ_d`.
                    let lam_d_inv = Complex64::new(1.0, 0.0) / lam[d];
                    mu_inv_tensor[(d, d)] = mu_inv * lam_d_inv;
                }
                assemble_tet_element_complex_anisotropic(vertices, eps_tensor, mu_inv_tensor)?
            };

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
                    // A_{αβ} = K_{αβ} - k0^2 M_{αβ}
                    let a_entry = signed
                        * (elem.k_local[(alpha, beta)] - k0_sq_c * elem.m_local[(alpha, beta)]);
                    triplets.push(Triplet::new(ii, jj, a_entry));
                }
            }
        }

        // ---- 4. Wave-port face scatter (same as the v3 path). ABC
        // faces are SKIPPED — the PML absorbs in the bulk.
        let mut rhs: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_interior];
        for (i, kind) in self.face_kinds.iter().enumerate() {
            let face = &self.exterior_faces.faces[i];
            match *kind {
                FaceKind::Pec | FaceKind::Abc => {
                    // PEC: row/column drop above already handled.
                    // Abc: surface-integral kernel suppressed; PML
                    //      handles absorption volumetrically.
                }
                FaceKind::WavePort(p) => {
                    let port = &self.ports[p];
                    let beta = (port.beta_mode)(omega);
                    if self.coupled_whitney {
                        let face_vertices = face.world_vertices(self.mesh);
                        let mut e_t_gauss = [Vector3::<f64>::zeros(); 3];
                        for (g, bary) in TRI_GAUSS_3PT_BARY.iter().enumerate() {
                            let p_g = bary[0] * face_vertices[0]
                                + bary[1] * face_vertices[1]
                                + bary[2] * face_vertices[2];
                            e_t_gauss[g] = (port.modal_e_t)(p_g);
                        }
                        self.scatter_port_face_gauss(
                            face,
                            beta,
                            e_t_gauss,
                            &interior_dof_of_edge,
                            &mut triplets,
                            &mut rhs,
                        );
                    } else {
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
        }

        let matrix = SparseColMat::try_new_from_triplets(n_interior, n_interior, &triplets)
            .map_err(|e| {
                Error::Numerical(format!(
                    "OpenBoundarySolver::assemble_driven_system_pml: \
                     failed to build driven matrix: {e:?}"
                ))
            })?;

        Ok(DrivenSystem {
            matrix,
            rhs,
            interior_edges,
            interior_dof_of_edge,
        })
    }

    /// Mean edge length of every tet in the mesh (m). Used by the
    /// CFS-PML config resolver to populate sentinel `sigma_max` /
    /// `alpha_max` parameters. Walks every tet's six local edges and
    /// averages their world-space lengths.
    fn mean_tet_edge_length(&self) -> f64 {
        let mut sum = 0.0_f64;
        let mut count = 0_usize;
        for tet in &self.mesh.tetrahedra {
            for &(li, lj) in LOCAL_EDGES.iter() {
                let v_a = self.mesh.vertices[tet[li]];
                let v_b = self.mesh.vertices[tet[lj]];
                sum += (v_b - v_a).norm();
                count += 1;
            }
        }
        if count == 0 {
            1.0e-3
        } else {
            sum / (count as f64)
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
// CFS-PML helpers (Phase 4.fem.eig.3.5 P4)
// ---------------------------------------------------------------------

/// Evaluate the diagonal `Λ(ω)` stretched-coordinate factor at a tet
/// classified by `class`.
///
/// Implements spec §3.1 (Roden–Gedney 2000 §II):
///
/// ```text
///     Λ = diag( s_y · s_z / s_x,
///               s_z · s_x / s_y,
///               s_x · s_y / s_z ),
///     s_α(ω) = κ_α(d_α) + σ_α(d_α) / (α_α + j ω ε_0).
/// ```
///
/// For a v3.5 single-axis PML tet (e.g. `PmlClass::PmlX`), only one
/// `s_α` deviates from `1`; the others collapse to unity. Interior
/// tets are guarded out by the caller — this helper always returns
/// the unit diagonal `[1, 1, 1]` for `PmlClass::Interior` so callers
/// can avoid a branch in their own scatter loop.
fn pml_stretching_lambda(class: crate::PmlClass, cfg: &PmlConfig, omega: f64) -> [Complex64; 3] {
    use crate::PmlClass;
    let m = cfg.m as f64;
    let d_total = (cfg.thickness_cells as f64).max(1.0); // sentinel guard
    // d_max for a tet centroid in the t-th brick layer is roughly
    // (t - 0.5) cells. We pick D = thickness_cells * h_cell at the
    // outer truncation surface and parameterise σ via d / D.
    // Here we pass `d / D` directly from the PmlClass payload, where
    // `d` is already in metres; we infer `D` from the cfg by
    // estimating the cell size as `d_max_observed / (thickness - 1 + 0.5)`.
    // To keep this self-contained, treat `d_total` as a normalised
    // count and substitute the per-cell depth d_cells = d_metres /
    // h_cell. We compute h_cell from the omega-resolved sigma_max
    // backward — but `cfg.resolved` already used h_cell. The cleanest
    // recovery is to scale by thickness_cells: `(d / D) ≈ d_cells /
    // thickness_cells`. We can derive d_cells if we know the lattice
    // spacing. For v3.5 we make the simplifying choice: pass d in
    // metres and assume the PML inner depth runs 0 → D_metres where
    // D_metres ≡ thickness_cells * h_cell_estimate. Since cfg's
    // `sigma_max` was already calibrated against `h_cell`, the
    // normalised depth `(d / D_metres) = (d / (thickness * h_cell))`.
    // We recover h_cell from the sigma_max - back-formula:
    //
    //     sigma_max = (m + 1) / (150 π h_cell)  ⇒
    //     h_cell    = (m + 1) / (150 π sigma_max).
    //
    // (This is robust to a caller passing a custom sigma_max; the
    // estimator just becomes a self-consistent "what h_cell does this
    // σ_max imply" rather than re-querying the mesh.)
    let h_cell = if cfg.sigma_max > 0.0 {
        (m + 1.0) / (150.0 * std::f64::consts::PI * cfg.sigma_max)
    } else {
        1.0e-3 // 1 mm fallback if grading parameters degenerate
    };
    let d_max = d_total * h_cell;

    let s_for = |d_alpha: f64| -> Complex64 {
        if d_alpha <= 0.0 || d_max <= 0.0 {
            return Complex64::new(1.0, 0.0);
        }
        // Polynomial grading: σ(d) = σ_max (d/D)^m, κ(d) = 1 + (κ_max - 1) (d/D)^m.
        let ratio = (d_alpha / d_max).clamp(0.0, 1.0);
        let pow = ratio.powf(m);
        let sigma_d = cfg.sigma_max * pow;
        let kappa_d = 1.0 + (cfg.kappa_max - 1.0) * pow;
        // s(ω) = κ + σ / (α + j ω ε_0)
        let denom = Complex64::new(cfg.alpha_max, omega * yee_core::units::EPS0);
        // Avoid division by zero at DC + zero alpha.
        if denom.norm_sqr() <= f64::MIN_POSITIVE {
            return Complex64::new(kappa_d, 0.0);
        }
        Complex64::new(kappa_d, 0.0) + Complex64::new(sigma_d, 0.0) / denom
    };

    let (sx, sy, sz) = match class {
        PmlClass::Interior => (
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
        ),
        PmlClass::PmlX { d } => (s_for(d), Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0)),
        PmlClass::PmlY { d } => (Complex64::new(1.0, 0.0), s_for(d), Complex64::new(1.0, 0.0)),
        PmlClass::PmlZ { d } => (Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0), s_for(d)),
        PmlClass::PmlXY { d_x, d_y } => (s_for(d_x), s_for(d_y), Complex64::new(1.0, 0.0)),
        PmlClass::PmlYZ { d_y, d_z } => (Complex64::new(1.0, 0.0), s_for(d_y), s_for(d_z)),
        PmlClass::PmlZX { d_z, d_x } => (s_for(d_x), Complex64::new(1.0, 0.0), s_for(d_z)),
        PmlClass::PmlXYZ { d_x, d_y, d_z } => (s_for(d_x), s_for(d_y), s_for(d_z)),
    };

    let one = Complex64::new(1.0, 0.0);
    // Λ = diag(s_y s_z / s_x, s_z s_x / s_y, s_x s_y / s_z)
    [
        if sx.norm_sqr() > 0.0 {
            sy * sz / sx
        } else {
            one
        },
        if sy.norm_sqr() > 0.0 {
            sz * sx / sy
        } else {
            one
        },
        if sz.norm_sqr() > 0.0 {
            sx * sy / sz
        } else {
            one
        },
    ]
}

/// Edge connectivity table mirroring `assembly::TetEdgeTable` so the
/// PML assembly path can scatter signed local blocks into the same
/// global-edge index space without depending on assembly's private
/// types.
#[derive(Debug, Clone)]
struct PmlAssemblyEdgeTable {
    edges: Vec<EdgeKey>,
    /// Per-tet local→global edge map + orientation signs (parallel to
    /// `assembly::TetEdgeConnectivity`).
    tet_edges: Vec<PmlTetEdgeConnectivity>,
}

#[derive(Debug, Clone, Copy)]
struct PmlTetEdgeConnectivity {
    global_edge: [usize; 6],
    sign: [f64; 6],
}

impl PmlAssemblyEdgeTable {
    fn build(mesh: &TetMesh3D) -> Self {
        let mut edge_map: HashMap<EdgeKey, usize> = HashMap::new();
        let mut edges: Vec<EdgeKey> = Vec::new();
        let mut tet_edges: Vec<PmlTetEdgeConnectivity> = Vec::with_capacity(mesh.tetrahedra.len());

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
                sign[alpha] = if a < b { 1.0 } else { -1.0 };
            }
            tet_edges.push(PmlTetEdgeConnectivity { global_edge, sign });
        }
        Self { edges, tet_edges }
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
