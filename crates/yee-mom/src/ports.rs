//! Port abstractions for planar MoM.
//!
//! Phase 1.3 walking skeleton: a [`Port`] trait plus two implementations.
//!
//! * [`DeltaGapPort`] — preserves the Phase 1.0 / mom-001 behaviour
//!   bit-for-bit: a 1 V (or user-supplied) impulse across every RWG edge
//!   tagged with `tag`, Galerkin-tested into `b[k] = V × length_k`.
//! * [`WavePort`] — **partially placeholder.** For Phase 1.3.1.0, the
//!   rectangular-waveguide TE10 closed-form mode is shipped via
//!   [`RectangularWaveguideTe10`] and selected with
//!   [`WavePort::with_rectangular_te10`]. Without a modal spec the
//!   [`WavePort::rhs`] is identical to a [`DeltaGapPort`] at the same
//!   voltage and tag — i.e. a uniform mode on a TEM-like cross-section.
//!   Microstrip / CPW / arbitrary cross-section modes via a numerical
//!   2-D eigensolver are deferred to Phase 1.3.1.1.
//!
//! The trait is `pub(crate)` because Phase 1.3 only exposes ports through
//! the high-level [`crate::PlanarMoM::run`] entry point. A public surface
//! will follow once the wave-port modal solver lands and the API is stable.
//! [`RectangularWaveguideTe10`] is `pub` so callers and integration tests
//! can query its analytic cutoff / β / Z / mode profile directly.

#![allow(dead_code)]
// Phase 1.3.0 wires the `Port` trait in but only the `DeltaGapPort` ::rhs /
// ::port_current / ::port_voltage path is reached by `PlanarMoM::run`. The
// `tag()` accessor and `WavePort` are part of the API surface for Phase
// 1.3.1+ and exercised by unit tests in this module.

use crate::basis::RwgBasis;
use faer::Mat;
use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::{MaterialTag, TriMesh2D};

#[doc(inline)]
pub use crate::eigensolver::ElementOrder;

/// Excitation source contract for the MoM solver.
///
/// A port owns:
/// * the RWG edge set it drives (identified by [`Port::tag`]),
/// * the construction of the Galerkin-tested RHS column vector at a given
///   frequency ([`Port::rhs`]),
/// * the extraction of the scalar port current from the solved current
///   vector ([`Port::port_current`]),
/// * the reference voltage used to compute `Z_in = V / I` ([`Port::port_voltage`]).
///
/// The solver in [`crate::solve`] is parameterised over `&dyn Port` so new
/// port models (lumped R/L/C, wave ports with modal weighting, edge ports,
/// etc.) can be added without touching the assembly/solve pipeline.
pub(crate) trait Port {
    /// Tag identifying which RWG edges belong to this port. Resolved by
    /// [`RwgBasis::port_basis_indices`].
    fn tag(&self) -> u32;
    /// Build the RHS column vector for this port at `freq_hz`.
    ///
    /// Returned vector has length `basis.n_basis()`; non-port entries are
    /// zero by construction.
    fn rhs(&self, basis: &RwgBasis, freq_hz: f64) -> Mat<Complex64>;
    /// Given the post-LU current vector `i`, extract the scalar port
    /// current in amperes — the Galerkin projection of the basis-weighted
    /// current back onto the port mode.
    fn port_current(&self, basis: &RwgBasis, i: &Mat<Complex64>) -> Complex64;
    /// Reference voltage `V` for `Z_in = V / I_port`.
    fn port_voltage(&self) -> Complex64;
}

/// Classical delta-gap source. Identical to the Phase 1.0 `delta_gap_rhs`:
/// `b[k] = V × length_k` on every edge with `tag()`, zero elsewhere.
pub(crate) struct DeltaGapPort {
    /// Port tag — must match the mesh's triangle-tag scheme used by
    /// [`RwgBasis::port_basis_indices`].
    pub tag: u32,
    /// Driving voltage. Phase 1.0 / mom-001 uses `1.0 + 0i`.
    pub voltage: Complex64,
}

impl Port for DeltaGapPort {
    fn tag(&self) -> u32 {
        self.tag
    }
    fn rhs(&self, basis: &RwgBasis, _freq_hz: f64) -> Mat<Complex64> {
        let n = basis.n_basis();
        let mut b = Mat::<Complex64>::zeros(n, 1);
        for k in basis.port_basis_indices(self.tag) {
            b[(k, 0)] = self.voltage * Complex64::new(basis.edges[k].length, 0.0);
        }
        b
    }
    fn port_current(&self, basis: &RwgBasis, i: &Mat<Complex64>) -> Complex64 {
        let mut total = Complex64::new(0.0, 0.0);
        for k in basis.port_basis_indices(self.tag) {
            total += Complex64::new(basis.edges[k].length, 0.0) * i[(k, 0)];
        }
        total
    }
    fn port_voltage(&self) -> Complex64 {
        self.voltage
    }
}

/// TEM-mode-weighted smoothed delta-gap port (Track WWWWWWW / P1 fix).
///
/// The classical [`DeltaGapPort`] drives every RWG basis function at the
/// port with `b_k = V · length_k`. Track TTTTTTT's port-edge diagnostic
/// (`tests/mom_002_port_edge_diagnostic.rs`) showed that on the
/// IIIIIII-reframed mom-002 strip mesh this uniform drive couples to a
/// single longitudinal-edge mode (alternating per-edge currents,
/// `|i|` ratios `~5×` between longitudinal-edge and diagonal-edge
/// basis functions) rather than the dominant quasi-TEM microstrip
/// mode. Track QQQQQQQ separately exonerated the Sommerfeld kernel
/// (`ε_eff_solver = 3.385` vs Hammerstad-Jensen `3.32`, +1.83 % error),
/// so the residual `|Im(Z)| ≈ 674 Ω` capacitive bias at 1 GHz is a
/// **port-excitation modeling** problem, not a kernel problem.
///
/// This port type weights the delta-gap RHS by an analytic Maxwell
/// transverse-mode envelope
///
/// ```text
///   w_TEM(y) = sqrt(2 / (π · (1 − (2 y / w)²)))
/// ```
///
/// peaked at the strip edges (`|2 y / w| → 1`) and minimum at strip
/// centre (`y = 0`). `y_k` is the y-coordinate of the basis function's
/// edge midpoint and `w` is the strip width (the
/// [`Self::strip_width_m`] field). The denominator is regularised
/// away from the singularity by clipping `|2 y / w| ≤ 1 − ε` with
/// `ε = 1e-3` so the weights stay finite on a uniform-y mesh that
/// places edges exactly on the strip edges. The `√(2/π)` normalisation
/// is the analytic `∫ w_TEM² dy / w = 1` factor for a Maxwell
/// `1/√(1−u²)` density on `u ∈ [-1, 1]`, but the result is
/// scale-invariant — `Z_in = V / I` is unchanged by a uniform rescale
/// of the weights — so the constant matters only for diagnostic
/// comparability across ports.
///
/// The same `w_TEM` weighting is applied **symmetrically** to the
/// port-current extraction so that `Z_in = V / I_port` retains the
/// Galerkin inner-product structure: a pure delta-gap recovers
/// bit-for-bit when `w_TEM ≡ 1` (the `strip_width_m → ∞` limit).
///
/// On a single-column wire port (the mom-001 dipole geometry) every
/// edge midpoint sits at `y = 0`, so `w_TEM(0) = √(2/π)` is a constant
/// uniform rescale and `Z_in` is **unchanged** from the [`DeltaGapPort`]
/// answer. This is the key property that lets mom-002 pick up the
/// TEM-smoothed port without disturbing the mom-001 NEC-4 gate.
pub(crate) struct TemSmoothedPort {
    /// Port tag — must match the mesh's triangle-tag scheme used by
    /// [`RwgBasis::port_basis_indices`].
    pub tag: u32,
    /// Driving voltage `V`. mom-002 uses `1.0 + 0i`.
    pub voltage: Complex64,
    /// Strip width `w` (meters). The transverse-mode envelope
    /// `w_TEM(y) = √(2 / (π · (1 − (2 y / w)²)))` evaluates the
    /// edge-singular Maxwell density on `y ∈ [-w/2, w/2]`. On the
    /// mom-002 IIIIIII reframe this is `2.94 mm`. A `0` width
    /// degenerates to the [`DeltaGapPort`] form (every weight is
    /// `√(2/π)`).
    pub strip_width_m: f64,
}

impl TemSmoothedPort {
    /// Evaluate the analytic Maxwell transverse-mode envelope at a
    /// single edge-midpoint y-coordinate. The
    /// `√(2 / (π · (1 − (2 y / w)²)))` form is the integrable
    /// edge-singular density for a thin strip carrying a quasi-TEM
    /// current; see Harrington, *Time-Harmonic Electromagnetic Fields*,
    /// §5.5 (Maxwell envelope) and the `tests/port_tem_smoothed_rhs.rs`
    /// validation gate.
    ///
    /// The denominator is clipped at `|u| ≤ 1 − ε` with `ε = 1e-3` so
    /// edges that sit exactly on the strip boundary (which a
    /// uniform-y mesh produces at `j = 0` and `j = n_width`) don't
    /// divide-by-zero. A zero or non-finite `w` falls back to the
    /// uniform `√(2/π)` weight — i.e. the [`DeltaGapPort`] form modulo
    /// a uniform rescale, which leaves `Z_in` unchanged.
    fn weight(&self, y: f64) -> f64 {
        let w = self.strip_width_m;
        if !w.is_finite() || w <= 0.0 {
            return (2.0 / std::f64::consts::PI).sqrt();
        }
        let u = (2.0 * y / w).abs().min(1.0 - 1e-3);
        (2.0 / (std::f64::consts::PI * (1.0 - u * u))).sqrt()
    }

    /// Sample the y-coordinate of the `k`-th basis function's
    /// shared-edge midpoint. The Maxwell envelope is a function of
    /// the transverse coordinate, so this is the point where the
    /// weight is evaluated.
    fn edge_y_midpoint(basis: &RwgBasis, k: usize) -> f64 {
        let edge = &basis.edges[k];
        let p0 = basis.mesh.vertices[edge.v0 as usize];
        let p1 = basis.mesh.vertices[edge.v1 as usize];
        0.5 * (p0.y + p1.y)
    }
}

impl Port for TemSmoothedPort {
    fn tag(&self) -> u32 {
        self.tag
    }
    fn rhs(&self, basis: &RwgBasis, _freq_hz: f64) -> Mat<Complex64> {
        let n = basis.n_basis();
        let mut b = Mat::<Complex64>::zeros(n, 1);
        for k in basis.port_basis_indices(self.tag) {
            let y = Self::edge_y_midpoint(basis, k);
            let w_tem = self.weight(y);
            b[(k, 0)] = self.voltage * Complex64::new(basis.edges[k].length * w_tem, 0.0);
        }
        b
    }
    fn port_current(&self, basis: &RwgBasis, i: &Mat<Complex64>) -> Complex64 {
        // Symmetric weighting on the Galerkin projection: the same
        // `length_k · w_TEM(y_k)` factor that built `b_k` extracts the
        // port current `I_port = Σ length_k · w_TEM(y_k) · i_k`. This
        // preserves the inner-product structure of `Z_in = V / I`.
        let mut total = Complex64::new(0.0, 0.0);
        for k in basis.port_basis_indices(self.tag) {
            let y = Self::edge_y_midpoint(basis, k);
            let w_tem = self.weight(y);
            total += Complex64::new(basis.edges[k].length * w_tem, 0.0) * i[(k, 0)];
        }
        total
    }
    fn port_voltage(&self) -> Complex64 {
        self.voltage
    }
}

/// Modal distribution attached to a [`WavePort`].
///
/// Phase 1.3.1.0 ships [`ModalDistribution::Te10`] for rectangular
/// waveguides; the default [`ModalDistribution::Uniform`] preserves the
/// Phase 1.3.0 delta-gap-equivalent placeholder behaviour. Phase 1.3.1.1
/// step 0-1 adds [`ModalDistribution::Numerical2D`] as a typed slot for
/// the numerical 2-D cross-section eigensolver; the eigensolve itself
/// lands in Phase 1.3.1.1 step 2-5, until which point the RHS falls
/// back to the uniform / delta-gap-equivalent behaviour of
/// [`ModalDistribution::Uniform`].
pub enum ModalDistribution {
    /// Uniform (TEM-like) dominant-mode amplitude across the port edges.
    /// Galerkin-tested into edge-length-weighted RWG basis functions this
    /// reduces to the delta-gap form, preserving bit-for-bit equivalence
    /// with [`DeltaGapPort`] for the mom-001 gate.
    Uniform,
    /// Closed-form TE10 mode of a rectangular waveguide cross-section.
    /// See [`RectangularWaveguideTe10`].
    Te10(RectangularWaveguideTe10),
    /// Numerical 2-D cross-section mode, computed by FEM eigensolve on
    /// an externally supplied [`TriMesh2D`]. See [`NumericalCrossSection`].
    ///
    /// **Status (Phase 1.3.1.1 step 7, this commit):** the eigensolve
    /// and Nedelec interpolation are wired in. After
    /// [`NumericalCrossSection::solve`] runs, the cached
    /// `mode_profile` is sampled at each port-edge midpoint via
    /// [`NumericalCrossSection::e_tangential_at`] and projected onto
    /// the edge's tangent unit vector — see [`WavePort::rhs`] for the
    /// formula. Before [`NumericalCrossSection::solve`] has run, the
    /// RHS is identically zero (the documented programmer-error
    /// path); callers must call `solve` first.
    ///
    /// Boxed because a [`NumericalCrossSection`] carries a full
    /// [`TriMesh2D`] (potentially hundreds of triangles) — keeping the
    /// payload behind an indirection prevents the enum from being
    /// dominated by a single rarely-built variant. Constructed via
    /// [`WavePort::with_numerical_cross_section`].
    Numerical2D(Box<NumericalCrossSection>),
}

/// Numerical 2-D cross-section mode for a [`WavePort`].
///
/// Phase 1.3.1.1 (full): runs a vector Helmholtz FEM eigensolve over a
/// `TriMesh2D` of the port cross-section to extract the dominant
/// quasi-TEM / quasi-TE mode, then caches `β` (propagation constant)
/// and `Z_w` (wave impedance) for use during RHS assembly.
///
/// Phase 1.3.1.1 step 0-1 (this commit): ships the type and the
/// [`WavePort::with_numerical_cross_section`] builder so the public API
/// freezes before the assembly code is written. [`Self::solve`] is a
/// stub that returns `Error::Unimplemented`. The `β` / `Z_w` cache
/// fields are `None` until a successful solve fills them.
///
/// Material data is keyed by [`MaterialTag`] (matching the tag in
/// [`TriMesh2D::triangle_material`]) — the caller supplies one
/// permittivity / permeability per distinct tag rather than per
/// triangle, which is how dielectric stack-ups are conventionally
/// described.
pub struct NumericalCrossSection {
    /// 2-D triangular mesh of the port cross-section.
    pub mesh: TriMesh2D,
    /// Complex relative permittivity per material tag.
    pub eps_r: HashMap<MaterialTag, Complex64>,
    /// Complex relative permeability per material tag.
    pub mu_r: HashMap<MaterialTag, Complex64>,
    /// Propagation constant `β` cached at the most recent
    /// [`Self::solve`] frequency. `None` before any successful solve.
    pub beta: Option<Complex64>,
    /// Wave impedance `Z_w` cached at the most recent [`Self::solve`]
    /// frequency. `None` before any successful solve.
    pub z_w: Option<Complex64>,
    /// Dominant-mode eigenvector (Nedelec edge-DoF amplitudes,
    /// **global-edge** indexing — already scattered out from the
    /// interior-edge DoF set with Dirichlet boundary edges set to 0).
    /// Cached on a successful [`Self::solve`]. `None` otherwise.
    ///
    /// Real-valued on the Phase 1.3.1.1 step 3 lossless path; stored
    /// as `Complex64` to mirror the future-proofed `beta` / `z_w`
    /// API and the assembly module's complex storage.
    ///
    /// Used by [`Self::e_tangential_at`] (Phase 1.3.1.1 step 7) to
    /// build the wave-port RHS by interpolating the modal `E_t`
    /// field at port-edge midpoints in the MoM-side mesh.
    pub mode_profile: Option<Vec<Complex64>>,
    /// Dominant-mode **longitudinal** field `E_z` (nodal-Lagrange
    /// vertex-DoF amplitudes, **global-vertex** indexing — already
    /// scattered out from the interior-vertex DoF set with Dirichlet
    /// boundary vertices set to 0). Cached on a successful
    /// [`Self::solve`]. `None` otherwise.
    ///
    /// Filled by the Phase 1.3.1.1 step 5 mixed `(E_t, E_z)`
    /// Lee-Sun-Cendes solve. On a homogeneous (air-filled) guide the
    /// dominant mode is purely transverse so this is ~zero; on an
    /// inhomogeneous (dielectric-loaded / microstrip) cross-section it
    /// carries the genuine longitudinal field that couples through the
    /// dielectric interface. Real-valued on the lossless path; stored
    /// as `Complex64` to mirror [`Self::mode_profile`].
    pub mode_profile_ez: Option<Vec<Complex64>>,
    /// Per-triangle local→global edge map (the `EdgeTable::tri_edges`
    /// payload) cached on a successful [`Self::solve`]. Needed by
    /// [`Self::e_tangential_at`] to interpolate the Nedelec edge basis
    /// without rebuilding the edge table on every sample. `None`
    /// otherwise.
    pub(crate) tri_edges_cache: Option<Vec<TriEdgesCacheEntry>>,
    /// Finite-element polynomial order for the cross-section mixed
    /// `(E_t, E_z)` eigensolve. Defaults to [`ElementOrder::First`] (the
    /// validated Whitney-1 Nedelec + linear-nodal path with full field
    /// reconstruction). [`ElementOrder::Second`] (Phase 1.3.1.1 step 5.6,
    /// set via [`Self::with_element_order`]) selects the second-order
    /// family for the high-contrast inhomogeneous case; see [`Self::solve`]
    /// for the p=2 field-reconstruction caveat.
    pub element_order: ElementOrder,
    /// When `true`, [`Self::solve`] dispatches to
    /// [`crate::eigensolver::solve_dense_mixed_quasi_tem`] instead of the
    /// default [`crate::eigensolver::solve_dense_mixed`].
    ///
    /// Set via [`Self::with_quasi_tem`]. Default: `false` (closed-guide
    /// path, bit-identical to previous behaviour).
    ///
    /// **Why a separate flag?** The closed-guide `solve_dense_mixed` selects
    /// the dominant mode by largest-β² among above-cutoff-floor candidates —
    /// correct for a rectangular waveguide above cutoff but it discards the
    /// microstrip quasi-TEM mode (`k_c² ≈ 0`, no cutoff). The quasi-TEM path
    /// uses a TEM-scale β-direct shift-invert ladder instead, validated by
    /// `tests/eigensolver_microstrip_quasi_tem.rs` (Phase 1.3.1.2, ADR-0060).
    /// Keeping the flag explicit preserves all existing closed-guide callers
    /// bit-for-bit; only callers that explicitly opt into quasi-TEM get the
    /// new path.
    use_quasi_tem: bool,
}

/// Compact per-triangle cache entry for the [`NumericalCrossSection`]
/// eigenmode-interpolation path: `(triangle_index, global_edge_indices,
/// orientation_signs)`. Mirrors a `pub(crate)` shape of the
/// `eigensolver::mesh::TriEdgeConnectivity` payload — kept here as a
/// crate-local type alias so the public [`NumericalCrossSection`]
/// struct does not re-export the assembly module's internals.
pub(crate) type TriEdgesCacheEntry = (usize, [usize; 3], [f64; 3]);

impl NumericalCrossSection {
    /// Build a cross-section mode descriptor with empty caches. The
    /// eigensolve is deferred to [`Self::solve`]; until that call (and
    /// until Phase 1.3.1.1 step 2-5 implements the eigensolve), `beta`
    /// and `z_w` remain `None`.
    pub fn new(
        mesh: TriMesh2D,
        eps_r: HashMap<MaterialTag, Complex64>,
        mu_r: HashMap<MaterialTag, Complex64>,
    ) -> Self {
        Self {
            mesh,
            eps_r,
            mu_r,
            beta: None,
            z_w: None,
            mode_profile: None,
            mode_profile_ez: None,
            tri_edges_cache: None,
            element_order: ElementOrder::First,
            use_quasi_tem: false,
        }
    }

    /// Select the finite-element polynomial order for the cross-section
    /// eigensolve (Phase 1.3.1.1 step 5.6). Non-breaking builder:
    /// [`ElementOrder::First`] is the default. [`ElementOrder::Second`]
    /// uses the validated p=2 element family (Nedelec-first-kind order-2 +
    /// quadratic nodal-Lagrange) for the high-contrast inhomogeneous case.
    ///
    /// **p=2 caveat.** At second order [`Self::solve`] caches the
    /// propagation constant `β` (the validated quantity) and a closed-form
    /// TE wave impedance `Z_w = η₀k₀/β`, but does **not** populate the
    /// [`Self::mode_profile`] / [`Self::e_tangential_at`] field
    /// reconstruction (the p=2 interior-face / edge-midpoint DoFs have no
    /// owning global edge / vertex, so the first-order edge-scatter does
    /// not apply). A `Numerical2D` wave-port RHS therefore still requires
    /// the first-order order; this knob exists so the validated p=2 β path
    /// is reachable end-to-end (e.g. for dispersion / β studies).
    ///
    /// **The p=2 `Z_w = η₀k₀/β` closed form is valid only on a homogeneous
    /// (`μ_r ≡ 1`) guide.** On an inhomogeneous (dielectric-loaded)
    /// cross-section it is **systematically incorrect** — the first-order
    /// path computes the energy-averaged `Z_w` (see [`Self::solve`]), but
    /// the p=2 path does not. Use the p=2 knob for `β` / dispersion studies,
    /// not for inhomogeneous `Z_w`.
    pub fn with_element_order(mut self, order: ElementOrder) -> Self {
        self.element_order = order;
        self
    }

    /// Enable the quasi-TEM mode-selection path for this cross-section
    /// eigensolve (Phase 1.3.1.2 / ADR-0060 adoption into `NumericalCrossSection`).
    ///
    /// Non-breaking builder: omitting it leaves the default closed-guide
    /// [`crate::eigensolver::solve_dense_mixed`] path bit-identical to all
    /// prior callers. Calling `.with_quasi_tem()` switches the dispatch to
    /// [`crate::eigensolver::solve_dense_mixed_quasi_tem`] — the validated
    /// TEM-scale β-direct shift-invert path (FR-4 50 Ω microstrip `ε_eff`
    /// within 1.2% of Hammerstad-Jensen, loose ≤10% box-truncation tol).
    ///
    /// **Applicable only to `First`-order (p=1) element meshes.** The p=2
    /// path does not define a quasi-TEM variant (out of scope — the p=2
    /// knob exists for closed-guide β / dispersion studies, not for
    /// microstrip excitation). If both `.with_element_order(Second)` and
    /// `.with_quasi_tem()` are called, `solve` returns an `Unimplemented`
    /// error rather than silently using the wrong solver.
    ///
    /// Use this when the cross-section carries a quasi-TEM mode (microstrip,
    /// CPW, strip-line) rather than a cutoff-bearing closed-waveguide mode.
    pub fn with_quasi_tem(mut self) -> Self {
        self.use_quasi_tem = true;
        self
    }

    /// Run the 2-D eigensolve at `freq_hz`.
    ///
    /// **Phase 1.3.1.1 step 5 (this commit):** assembles the full mixed
    /// `(E_t, E_z)` Lee-Sun-Cendes block generalized eigenproblem
    /// `A x = k_c² B x` on [`Self::mesh`] (via
    /// [`crate::eigensolver::assemble_mixed`], with `eps_r` / `mu_r`
    /// looked up per [`yee_mesh::MaterialTag`]), dense-solves the dominant
    /// quasi-TEM mode (via [`crate::eigensolver::solve_dense_mixed`] — the
    /// `B⁻¹A` non-symmetric path, since the block pencil is symmetric
    /// *indefinite*), caches `β = √(k_0² − k_c²)`, and extracts a
    /// **numerical** wave impedance `Z_w` off the solved eigenvector
    /// (see below) — replacing the earlier TE-mode `η₀k₀/β`
    /// approximation.
    ///
    /// On a **homogeneous** (air-filled) cross-section the longitudinal
    /// `E_z` block decouples and the result reproduces the transverse-
    /// only solve to machine precision (the WR-90 TE10 gate stays green
    /// within its 1 % tolerance — DoD-V1). On an **inhomogeneous**
    /// (dielectric-loaded / microstrip) cross-section the `E_z` coupling
    /// through the dielectric interface shifts `β` and the mode shape —
    /// the capability this solver exists for.
    ///
    /// **Numerical Z_w.** For the dominant mode the energy-averaged wave
    /// impedance is
    ///
    /// ```text
    ///   Z_w = (ω μ₀ / β) · ( ∫_Ω |E_t|² dA ) / ( ∫_Ω (1/μ_r) |E_t|² dA )
    /// ```
    ///
    /// from the modal field relation `H_t = (β / ω μ₀ μ_r)(ẑ × E_t)`
    /// (so the local wave impedance is `ω μ₀ μ_r / β` and the above is
    /// its `|E_t|²`-energy-weighted average). On a homogeneous guide
    /// `μ_r ≡ 1` this reduces **exactly** to `ω μ₀ / β = η₀ k₀ / β`,
    /// the TE-mode value the closed-form
    /// [`RectangularWaveguideTe10::wave_impedance`] gives — the DoD-V3
    /// regression guard. The two transverse-energy integrals are
    /// computed by Nedelec quadrature off the cached eigenvector
    /// (`local_b_ee_mass`-style closed forms).
    ///
    /// The dense path is `O(n³)` in `n = n_t + n_z`; viable for coarse
    /// cross-sections (≤ a few hundred DoF). The WR-90 case lands at
    /// `n ≈ 121`. Sparse shift-and-invert is a later step.
    pub fn solve(&mut self, freq_hz: f64) -> yee_core::Result<()> {
        use crate::eigensolver::{
            ElementOrder,
            assembly::{assemble_mixed, assemble_mixed_p2},
            mesh::EdgeTable,
            solve_dense_mixed, solve_dense_mixed_quasi_tem,
        };

        // Quasi-TEM + Second-order: unsupported combination. The p=2 knob
        // exists for closed-guide β/dispersion studies; the quasi-TEM
        // path is p=1 only (the validated path, Phase 1.3.1.2 / ADR-0060).
        if self.use_quasi_tem && matches!(self.element_order, ElementOrder::Second) {
            return Err(yee_core::Error::Unimplemented(
                "NumericalCrossSection: quasi-TEM path is p=1 only; \
                 Second-order + quasi-TEM is out of scope",
            ));
        }

        let table = EdgeTable::build(&self.mesh);
        let asm = match self.element_order {
            ElementOrder::First => assemble_mixed(&self.mesh, &self.eps_r, &self.mu_r, &table),
            ElementOrder::Second => assemble_mixed_p2(&self.mesh, &self.eps_r, &self.mu_r, &table),
        };
        // Dispatch: quasi-TEM path uses the TEM-scale β-direct shift-invert
        // ladder (Phase 1.3.1.2 / ADR-0060); the default closed-guide path
        // is unchanged for all existing callers.
        let sol = if self.use_quasi_tem {
            solve_dense_mixed_quasi_tem(&asm, freq_hz)?
        } else {
            solve_dense_mixed(&asm, freq_hz)?
        };
        // β = √(β²). Lossless inputs give real β² ≥ 0; take the
        // principal square root (positive real branch).
        let beta_sq = sol.beta_sq;
        let beta = if beta_sq.im.abs() < 1e-9 * beta_sq.re.abs() {
            Complex64::new(beta_sq.re.max(0.0).sqrt(), 0.0)
        } else {
            beta_sq.sqrt()
        };
        self.beta = Some(beta);

        // ── p=2 path (Phase 1.3.1.1 step 5.6): β is the validated quantity;
        // the first-order edge-scatter field reconstruction does NOT apply
        // (interior-face / edge-midpoint DoFs carry `usize::MAX` sentinels in
        // `interior_to_global_*`). Cache β + a closed-form TE wave impedance
        // and leave the field-profile caches `None`; see
        // `with_element_order`'s caveat. ──
        if matches!(self.element_order, ElementOrder::Second) {
            let omega = std::f64::consts::TAU * freq_hz;
            let k0 = omega / yee_core::units::C0;
            self.z_w = Some(Complex64::new(yee_core::units::ETA0 * k0, 0.0) / beta);
            self.mode_profile = None;
            self.mode_profile_ez = None;
            self.tri_edges_cache = None;
            return Ok(());
        }

        // Scatter the interior-edge DoF eigenvector out to global-edge
        // indexing (E_t). PEC boundary edges have E_t = 0 by Dirichlet
        // elimination so they remain zero in the global profile.
        let mut global_mode = vec![Complex64::new(0.0, 0.0); table.n_edges()];
        for (i_dof, &gid) in asm.interior_to_global_edges.iter().enumerate() {
            global_mode[gid] = sol.e_t[i_dof];
        }
        self.mode_profile = Some(global_mode);

        // Scatter the interior-vertex DoF eigenvector out to global-
        // vertex indexing (E_z). PEC boundary vertices are 0 by Dirichlet
        // elimination.
        let mut global_ez = vec![Complex64::new(0.0, 0.0); self.mesh.n_verts()];
        for (i_dof, &vid) in asm.interior_to_global_verts.iter().enumerate() {
            global_ez[vid] = sol.e_z[i_dof];
        }
        self.mode_profile_ez = Some(global_ez);

        // Cache per-triangle (tri_idx, global_edge_indices, orientation
        // signs) so e_tangential_at can interpolate without rebuilding
        // the edge table. Compact tuple form avoids a public re-export
        // of the crate-private `TriEdgeConnectivity` type.
        let tri_edges: Vec<TriEdgesCacheEntry> = table
            .tri_edges
            .iter()
            .enumerate()
            .map(|(t, c)| (t, c.global_edge, c.sign))
            .collect();
        self.tri_edges_cache = Some(tri_edges.clone());

        // Numerical wave impedance Z_w (see method docstring): the
        // |E_t|²-energy-weighted modal wave impedance, which reduces to
        // `η₀ k₀ / β` on the homogeneous guide and shifts with the
        // dielectric fill otherwise.
        let omega = std::f64::consts::TAU * freq_hz;
        let mu0 = yee_core::units::MU0;
        // Accumulate ∫|E_t|² dA and ∫(1/μ_r)|E_t|² dA over the mesh by
        // Nedelec quadrature off the cached eigenvector.
        let (energy_plain, energy_inv_mu) = self.transverse_energy_integrals(&table);
        let z_w = if beta.norm() > 0.0 && energy_inv_mu > 0.0 {
            Complex64::new(omega * mu0 / beta.re * (energy_plain / energy_inv_mu), 0.0)
        } else {
            // Degenerate fallback to the closed-form TE relation so a
            // single-DoF / pathological mesh still caches a finite value.
            let k0 = omega / yee_core::units::C0;
            Complex64::new(yee_core::units::ETA0 * k0, 0.0) / beta
        };
        self.z_w = Some(z_w);

        Ok(())
    }

    /// Compute the two transverse-energy integrals
    /// `(∫_Ω |E_t|² dA, ∫_Ω (1/μ_r) |E_t|² dA)` over the cross-section
    /// by Nedelec quadrature off the cached `mode_profile` eigenvector.
    ///
    /// Used by the numerical `Z_w` extraction in [`Self::solve`]. The
    /// `∫|E_t|² dA` integrand is `N_i · N_j` (the unweighted Nedelec
    /// mass), assembled from the same closed-form linear-triangle moments
    /// the element-matrix code uses; the `1/μ_r`-weighted variant divides
    /// each triangle's contribution by its `μ_r` so the ratio reduces to
    /// `1` on a homogeneous non-magnetic guide.
    fn transverse_energy_integrals(
        &self,
        table: &crate::eigensolver::mesh::EdgeTable,
    ) -> (f64, f64) {
        let Some(profile) = &self.mode_profile else {
            return (0.0, 0.0);
        };
        let default_one = Complex64::new(1.0, 0.0);
        let local_edges: [[usize; 2]; 3] = [[1, 2], [2, 0], [0, 1]];
        let mut energy_plain = 0.0_f64;
        let mut energy_inv_mu = 0.0_f64;
        for (t_idx, conn) in table.tri_edges.iter().enumerate() {
            let tri = self.mesh.triangles[t_idx];
            let v: [[f64; 2]; 3] = [
                self.mesh.vertices[tri[0]],
                self.mesh.vertices[tri[1]],
                self.mesh.vertices[tri[2]],
            ];
            let area = self.mesh.area(t_idx);
            let tag = self.mesh.triangle_material[t_idx];
            let mu = self.mu_r.get(&tag).unwrap_or(&default_one).re;
            // Local FEM scalars (Jin convention).
            let mut bb = [0.0; 3];
            let mut cc = [0.0; 3];
            let mut ell = [0.0; 3];
            for i in 0..3 {
                let i1 = (i + 1) % 3;
                let i2 = (i + 2) % 3;
                bb[i] = v[i1][1] - v[i2][1];
                cc[i] = v[i2][0] - v[i1][0];
            }
            for (e, &[a, b]) in local_edges.iter().enumerate() {
                let dx = v[b][0] - v[a][0];
                let dy = v[b][1] - v[a][1];
                ell[e] = (dx * dx + dy * dy).sqrt();
            }
            let grad_dot =
                |p: usize, q: usize| (bb[p] * bb[q] + cc[p] * cc[q]) / (4.0 * area * area);
            let int_ll = |p: usize, q: usize| area * (1.0 + if p == q { 1.0 } else { 0.0 }) / 12.0;
            // Local |E_t|² = Σ_ij amp_i amp_j ∫ N_i·N_j dA.
            let mut tri_energy = 0.0_f64;
            for li in 0..3 {
                let [ai, bi] = local_edges[li];
                let amp_i = profile[conn.global_edge[li]].re * conn.sign[li] * ell[li];
                for lj in 0..3 {
                    let [aj, bj] = local_edges[lj];
                    let amp_j = profile[conn.global_edge[lj]].re * conn.sign[lj] * ell[lj];
                    // ∫ (λ_ai∇λ_bi − λ_bi∇λ_ai)·(λ_aj∇λ_bj − λ_bj∇λ_aj) dA
                    let m = int_ll(ai, aj) * grad_dot(bi, bj)
                        - int_ll(ai, bj) * grad_dot(bi, aj)
                        - int_ll(bi, aj) * grad_dot(ai, bj)
                        + int_ll(bi, bj) * grad_dot(ai, aj);
                    tri_energy += amp_i * amp_j * m;
                }
            }
            energy_plain += tri_energy;
            energy_inv_mu += tri_energy / mu;
        }
        (energy_plain, energy_inv_mu)
    }

    /// Evaluate the dominant-mode transverse electric field
    /// `E_t = (E_x, E_y)` at cross-section coordinate `(x, y)` by
    /// Nedelec interpolation.
    ///
    /// Locates the mesh triangle containing `(x, y)`, computes
    /// barycentric coordinates, and sums the three local Nedelec edge
    /// basis functions weighted by the cached eigenvector components.
    /// Returns `[0.0, 0.0]` if `(x, y)` lies outside the mesh or
    /// before [`Self::solve`] has been called.
    ///
    /// **Sign / scale convention.** The eigenvector returned by the
    /// dense eigensolve is determined up to a global scalar
    /// (eigenvectors of `M y = λ y` are scale-free). The Phase 1.3.1.1
    /// step 3 path fixes the global sign so the largest-magnitude DoF
    /// is positive but leaves the scale arbitrary — callers comparing
    /// against an analytic reference must normalize. The Phase 1.3.1.1
    /// step 7 wave-port RHS [`WavePort::rhs`] consumes this directly
    /// in its `Numerical2D` arm and inherits the same scale-freedom;
    /// the impedance `Z_in = V / I` is scale-invariant under any
    /// global rescaling of the modal RHS, so this is benign for the
    /// scattering / `Z_in` extraction. The downstream
    /// `tests/wave_port_numerical_te10.rs` validation gate
    /// explicitly renormalizes both `b_num` and `b_analytic` to unit
    /// L2 norm before computing the 1 % agreement.
    pub fn e_tangential_at(&self, x: f64, y: f64) -> [f64; 2] {
        let Some(profile) = &self.mode_profile else {
            return [0.0, 0.0];
        };
        let Some(tri_edges) = &self.tri_edges_cache else {
            return [0.0, 0.0];
        };
        for (t_idx, global_edge, sign) in tri_edges {
            let tri = self.mesh.triangles[*t_idx];
            let v: [[f64; 2]; 3] = [
                self.mesh.vertices[tri[0]],
                self.mesh.vertices[tri[1]],
                self.mesh.vertices[tri[2]],
            ];
            let area = self.mesh.area(*t_idx);
            // Barycentric coordinates of (x, y) wrt CCW triangle v0,v1,v2.
            // λ_i = ((b_i, c_i) · ((x,y) - v_origin) + const) / (2A),
            // computed directly via sub-triangle areas.
            let sub_area = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| -> f64 {
                0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
            };
            let p = [x, y];
            let lam0 = sub_area(p, v[1], v[2]) / area;
            let lam1 = sub_area(v[0], p, v[2]) / area;
            let lam2 = sub_area(v[0], v[1], p) / area;
            // Point lies inside iff all three barycentric coordinates
            // are non-negative (small tolerance for floating-point noise
            // on triangle boundaries).
            let eps = -1e-12;
            if lam0 < eps || lam1 < eps || lam2 < eps {
                continue;
            }
            // ∇λ_i = (b_i, c_i) / (2A), with the Jin convention
            // b[i] = y_{i+1} - y_{i+2}, c[i] = x_{i+2} - x_{i+1}.
            let mut bb = [0.0; 3];
            let mut cc = [0.0; 3];
            for i in 0..3 {
                let i1 = (i + 1) % 3;
                let i2 = (i + 2) % 3;
                bb[i] = v[i1][1] - v[i2][1];
                cc[i] = v[i2][0] - v[i1][0];
            }
            let grad = |i: usize| -> [f64; 2] { [bb[i] / (2.0 * area), cc[i] / (2.0 * area)] };
            let lam = [lam0, lam1, lam2];
            // Local edge endpoints per `eigensolver::mesh::LOCAL_EDGES`:
            // edge `e` opposite local vertex `e`, traversed CCW.
            // Edge 0: v1 → v2, edge 1: v2 → v0, edge 2: v0 → v1.
            let local_edges: [[usize; 2]; 3] = [[1, 2], [2, 0], [0, 1]];
            let mut e_field = [0.0f64; 2];
            for (le, &[a, b]) in local_edges.iter().enumerate() {
                // Edge length matches the canonical `EdgeKey::new` ordering
                // (smaller→larger vertex) since lengths are direction-
                // independent.
                let dx = v[b][0] - v[a][0];
                let dy = v[b][1] - v[a][1];
                let ell = (dx * dx + dy * dy).sqrt();
                let sigma = sign[le];
                let gid = global_edge[le];
                let amp = profile[gid].re;
                let ga = grad(a);
                let gb = grad(b);
                // N_e = ℓ σ (λ_a ∇λ_b − λ_b ∇λ_a)
                let nx = ell * sigma * (lam[a] * gb[0] - lam[b] * ga[0]);
                let ny = ell * sigma * (lam[a] * gb[1] - lam[b] * ga[1]);
                e_field[0] += amp * nx;
                e_field[1] += amp * ny;
            }
            return e_field;
        }
        [0.0, 0.0]
    }
}

/// Closed-form TE10 mode of a rectangular waveguide of inner dimensions
/// `a × b` (with `a > b`, conventional).
///
/// Phase 1.3.1.0: analytic mode only. A numerical 2-D eigensolver for
/// arbitrary cross-sections is Phase 1.3.1.1.
///
/// Reference: Pozar, *Microwave Engineering* 4th ed., §3.3.
///
/// The mode profile in the cross-section, with the convention that `x`
/// runs along the long dimension `a` and `y` along the short dimension
/// `b`, is
///
/// ```text
///   E_y(x, y, z) = E_0 sin(π x / a) exp(-j β z)
///   H_x(x, y, z) = -(E_0 / Z_TE10) sin(π x / a) exp(-j β z)
/// ```
///
/// with cutoff `f_c = c / (2 a √ε_r)`, phase constant
/// `β_10 = sqrt(k² - (π/a)²)`, and wave impedance
/// `Z_TE10 = η / sqrt(1 - (f_c/f)²)` where `η = η_0 / √ε_r` is the
/// intrinsic impedance of the fill medium. Below cutoff the mode is
/// evanescent; Phase 1.3.1.0 rejects that regime by returning `NaN`
/// from [`Self::beta`] and [`Self::wave_impedance`].
pub struct RectangularWaveguideTe10 {
    /// Long inner dimension of the waveguide cross-section (meters).
    pub a: f64,
    /// Short inner dimension (meters).
    pub b: f64,
    /// Relative permittivity of the fill medium. Use `1.0` for air.
    pub eps_r: f64,
}

impl RectangularWaveguideTe10 {
    /// Cutoff frequency `f_c = c / (2 a √ε_r)` for the TE10 mode.
    pub fn cutoff_hz(&self) -> f64 {
        let c = yee_core::units::C0 / self.eps_r.sqrt();
        c / (2.0 * self.a)
    }

    /// Phase constant `β_10` at frequency `freq_hz`. Returns `NaN` at or
    /// below cutoff (the evanescent regime is out of scope for Phase
    /// 1.3.1.0; callers should reject this case before driving a sweep).
    pub fn beta(&self, freq_hz: f64) -> f64 {
        let k = std::f64::consts::TAU * freq_hz * self.eps_r.sqrt() / yee_core::units::C0;
        let kc = std::f64::consts::PI / self.a;
        if k <= kc {
            return f64::NAN;
        }
        (k * k - kc * kc).sqrt()
    }

    /// Wave impedance `Z_TE10` at frequency `freq_hz`. Returns `NaN` at or
    /// below cutoff.
    pub fn wave_impedance(&self, freq_hz: f64) -> f64 {
        let eta = yee_core::units::ETA0 / self.eps_r.sqrt();
        let fc = self.cutoff_hz();
        if freq_hz <= fc {
            return f64::NAN;
        }
        eta / (1.0 - (fc / freq_hz).powi(2)).sqrt()
    }

    /// Sample the transverse `E_y` modal profile at cross-section
    /// coordinates `(x, y)`.
    ///
    /// Returns `sin(π x / a)`; the TE10 mode is uniform in `y`. Domain:
    /// `x ∈ [0, a]`, `y ∈ [0, b]`. Returns `0` outside the cross-section
    /// rectangle, including on the conducting walls.
    pub fn e_y_profile(&self, x: f64, y: f64) -> f64 {
        if x < 0.0 || x > self.a || y < 0.0 || y > self.b {
            return 0.0;
        }
        (std::f64::consts::PI * x / self.a).sin()
    }
}

/// Wave port — 1D modal source on a tagged edge set.
///
/// **Status (Phase 1.3.1.0):** The closed-form rectangular-waveguide TE10
/// mode is shipped via [`Self::with_rectangular_te10`]. A `WavePort`
/// constructed without a modal spec defaults to [`ModalDistribution::Uniform`]
/// and is bit-for-bit equivalent to [`DeltaGapPort`] at the same voltage
/// and tag, preserving the Phase 1.3.0 behaviour and the mom-001 gate.
/// Microstrip / CPW / arbitrary cross-sections still degenerate to the
/// uniform placeholder pending the Phase 1.3.1.1 numerical 2-D eigensolver.
pub struct WavePort {
    /// Port tag — matches the mesh tagging scheme.
    pub tag: u32,
    /// Modal-source reference voltage. With [`ModalDistribution::Uniform`]
    /// this maps directly onto the equivalent delta-gap drive amplitude.
    pub voltage: Complex64,
    /// Phase-velocity factor for the lowest-order mode on the port cross
    /// section, as a fraction of `c₀`. Used by the (currently unmodulated)
    /// `β = ω / (c₀ · v_factor)` term that will dress the RHS once
    /// long-section propagation is wired through. Pre-1.3.1.0 callers
    /// continue to set this directly; the rectangular-waveguide path
    /// instead derives `β` from [`RectangularWaveguideTe10::beta`].
    pub mode_phase_velocity_factor: f64,
    /// Modal distribution. Defaults to [`ModalDistribution::Uniform`];
    /// set via [`Self::with_rectangular_te10`] for a TE10 rectangular
    /// waveguide.
    pub modal_distribution: ModalDistribution,
}

impl WavePort {
    /// Attach a closed-form rectangular-waveguide TE10 mode to this
    /// wave-port. Non-breaking builder: omitting it leaves the default
    /// [`ModalDistribution::Uniform`], which is bit-for-bit equivalent to
    /// a [`DeltaGapPort`] at the same voltage and tag.
    pub fn with_rectangular_te10(mut self, mode: RectangularWaveguideTe10) -> Self {
        self.modal_distribution = ModalDistribution::Te10(mode);
        self
    }

    /// Attach a numerical 2-D cross-section mode to this wave-port.
    ///
    /// Callers must invoke [`NumericalCrossSection::solve`] on `mode`
    /// **before** placing it on the port — `WavePort::rhs` needs the
    /// cached `mode_profile` to sample the Nedelec eigenmode at port-
    /// edge midpoints. A `Numerical2D` port whose mode has not been
    /// solved produces an all-zero RHS (the documented degenerate
    /// path); see [`ModalDistribution::Numerical2D`] for the contract.
    pub fn with_numerical_cross_section(mut self, mode: NumericalCrossSection) -> Self {
        self.modal_distribution = ModalDistribution::Numerical2D(Box::new(mode));
        self
    }
}

impl Port for WavePort {
    fn tag(&self) -> u32 {
        self.tag
    }
    fn rhs(&self, basis: &RwgBasis, freq_hz: f64) -> Mat<Complex64> {
        // Phase 1.3.1.0: the wave-port modal weighting differs from
        // delta-gap by distributing the source across port edges
        // according to the mode field. For a uniform (TEM-like)
        // dominant-mode approximation the distribution is uniform, so
        // we degenerate to delta-gap. For TE10 on a rectangular
        // waveguide we evaluate sin(π x / a) at the port-edge midpoint
        // using the edge's centroid x-coordinate, weight by edge
        // length, and scale by the analytic β / Z_TE10 implicit in the
        // mode shape. The β term computed below for the uniform case
        // remains reserved for the Phase 1.3.1.1 propagation correction.
        let n = basis.n_basis();
        let mut b = Mat::<Complex64>::zeros(n, 1);
        match &self.modal_distribution {
            ModalDistribution::Uniform => {
                let omega = std::f64::consts::TAU * freq_hz;
                let beta = omega / (yee_core::units::C0 * self.mode_phase_velocity_factor);
                let _ = beta; // reserved for Phase 1.3.1.1
                for k in basis.port_basis_indices(self.tag) {
                    b[(k, 0)] = self.voltage * Complex64::new(basis.edges[k].length, 0.0);
                }
            }
            ModalDistribution::Te10(mode) => {
                // Sample the analytic TE10 E_y profile at each port-edge
                // midpoint. The convention is that the long dimension `a`
                // runs along the mesh's local x axis; callers must align
                // the port mesh accordingly. The profile is real and
                // bounded by `sin(π x / a) ∈ [0, 1]`, so the resulting
                // RHS is a real-scaled version of the delta-gap RHS.
                for k in basis.port_basis_indices(self.tag) {
                    let edge = &basis.edges[k];
                    let p0 = basis.mesh.vertices[edge.v0 as usize];
                    let p1 = basis.mesh.vertices[edge.v1 as usize];
                    let mid_x = 0.5 * (p0.x + p1.x);
                    let mid_y = 0.5 * (p0.y + p1.y);
                    let profile = mode.e_y_profile(mid_x, mid_y);
                    b[(k, 0)] = self.voltage * Complex64::new(edge.length * profile, 0.0);
                }
            }
            ModalDistribution::Numerical2D(mode) => {
                // Phase 1.3.1.1 step 7: the numerical-cross-section
                // wave-port RHS samples the Nedelec eigenmode at each
                // port-edge midpoint and projects it onto the edge's
                // tangent unit vector, weighted by edge length and the
                // driving voltage.
                //
                //   b[k] = V · ℓ_k · (E_t(x_mid, y_mid) · t̂_k)
                //
                // The cross-section coordinate convention matches the
                // [`ModalDistribution::Te10`] arm: the MoM-side mesh's
                // (x, y) coordinates are taken as the cross-section
                // (x, y) coordinates, so the same RWG port-edge mesh
                // can be paired with either a closed-form TE10 or a
                // numerical 2-D eigenmode by swapping the modal
                // distribution.
                //
                // Sign convention: the dominant eigenvector returned
                // by [`crate::eigensolver::solve_dense`] has its
                // largest-magnitude DoF pinned positive; this
                // corresponds to the positive-going wave (`β > 0`)
                // selected by the smallest-strictly-positive `k_c²`
                // branch. Callers comparing to analytic references
                // typically renormalize both sides to unit L2 — see
                // `tests/wave_port_numerical_te10.rs`.
                //
                // If the modal profile has not been solved yet, the
                // mode field is zero everywhere and the resulting RHS
                // is zero. This is the documented degenerate path
                // (rather than the legacy uniform fallback) because a
                // post-`new` / pre-`solve` `NumericalCrossSection` is
                // a programmer error — the solver expects `solve`
                // to have run.
                if mode.mode_profile.is_none() {
                    // No cached profile — falls back to all-zero RHS
                    // for the Numerical2D arm. mom-001 and the existing
                    // wave-port tests do not exercise this arm.
                    // Returned `b` is already zero-initialized.
                    return b;
                }
                for k in basis.port_basis_indices(self.tag) {
                    let edge = &basis.edges[k];
                    let p0 = basis.mesh.vertices[edge.v0 as usize];
                    let p1 = basis.mesh.vertices[edge.v1 as usize];
                    let mid_x = 0.5 * (p0.x + p1.x);
                    let mid_y = 0.5 * (p0.y + p1.y);
                    let e_field = mode.e_tangential_at(mid_x, mid_y);
                    // Edge tangent unit vector in the cross-section
                    // (x, y) plane. Use (v0 → v1) direction with the
                    // 2-D projection.
                    let tx = p1.x - p0.x;
                    let ty = p1.y - p0.y;
                    let tn = (tx * tx + ty * ty).sqrt();
                    let (tux, tuy) = if tn > 0.0 {
                        (tx / tn, ty / tn)
                    } else {
                        (0.0, 0.0)
                    };
                    let projection = e_field[0] * tux + e_field[1] * tuy;
                    b[(k, 0)] = self.voltage * Complex64::new(edge.length * projection, 0.0);
                }
            }
        }
        b
    }
    fn port_current(&self, basis: &RwgBasis, i: &Mat<Complex64>) -> Complex64 {
        let mut total = Complex64::new(0.0, 0.0);
        for k in basis.port_basis_indices(self.tag) {
            total += Complex64::new(basis.edges[k].length, 0.0) * i[(k, 0)];
        }
        total
    }
    fn port_voltage(&self) -> Complex64 {
        self.voltage
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::RwgBasis;
    use nalgebra::Vector3;
    use yee_mesh::TriMesh;

    fn two_tri_mesh_with_port() -> TriMesh {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
            Vector3::new(0.1, 0.1, 0.0),
            Vector3::new(0.0, 0.1, 0.0),
        ];
        let triangles = vec![[0u32, 1, 2], [0u32, 2, 3]];
        let tags = vec![1u32, 2u32];
        TriMesh::new(vertices, triangles, tags).unwrap()
    }

    #[test]
    fn delta_gap_port_rhs_matches_legacy() {
        let basis = RwgBasis::from_mesh(two_tri_mesh_with_port()).unwrap();
        let port = DeltaGapPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
        };
        let b = port.rhs(&basis, 1.0e9);
        for k in basis.port_basis_indices(1) {
            let expected = basis.edges[k].length;
            assert!((b[(k, 0)].re - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn wave_port_rhs_matches_delta_gap_for_unit_mode() {
        // Phase 1.3.0: WavePort with mode_phase_velocity_factor = 1.0 degenerates to DeltaGapPort.
        let basis = RwgBasis::from_mesh(two_tri_mesh_with_port()).unwrap();
        let dg = DeltaGapPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
        };
        let wp = WavePort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
            mode_phase_velocity_factor: 1.0,
            modal_distribution: ModalDistribution::Uniform,
        };
        let b1 = dg.rhs(&basis, 1.0e9);
        let b2 = wp.rhs(&basis, 1.0e9);
        let n = basis.n_basis();
        for k in 0..n {
            assert!((b1[(k, 0)] - b2[(k, 0)]).norm() < 1e-15);
        }
    }

    fn unit_square_cross_section() -> NumericalCrossSection {
        // Trivial 2-tri cross-section spanning the unit square; the
        // contents don't matter for the stub-equivalence test (the
        // eigensolve is not run), only that the type constructs.
        let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = vec![[0, 1, 2], [0, 2, 3]];
        let mesh = TriMesh2D::new(vertices, triangles, None, None).unwrap();
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        NumericalCrossSection::new(mesh, eps, mu)
    }

    #[test]
    fn numerical_cross_section_new_initializes_caches_to_none() {
        let m = unit_square_cross_section();
        assert!(m.beta.is_none());
        assert!(m.z_w.is_none());
    }

    #[test]
    fn numerical_cross_section_solve_fills_caches_on_unit_square() {
        // Phase 1.3.1.1 step 3: the eigensolve now runs on the
        // unit-square 2-tri fixture. With only 1 interior edge the
        // result is not physically meaningful, but the call must
        // succeed without panicking and populate the caches with
        // finite values — or return an `Error::Numerical` if the
        // single-DoF problem admits no positive β². Either is
        // acceptable as a smoke gate; what we explicitly forbid is
        // the old `Unimplemented` stub.
        let mut m = unit_square_cross_section();
        match m.solve(10e9) {
            Ok(()) => {
                let beta = m.beta.expect("β should be cached on success");
                let z_w = m.z_w.expect("Z_w should be cached on success");
                assert!(beta.re.is_finite() && beta.im.is_finite());
                assert!(z_w.re.is_finite() && z_w.im.is_finite());
            }
            Err(yee_core::Error::Numerical(_)) => {
                // Degenerate single-DoF case — no physically valid β².
                // Caches must remain unfilled on a failed solve.
                assert!(m.beta.is_none());
                assert!(m.z_w.is_none());
            }
            other => panic!("unexpected solve outcome: {other:?}"),
        }
    }

    #[test]
    fn wave_port_numerical_rhs_is_zero_before_solve() {
        // Phase 1.3.1.1 step 7: the Numerical2D arm now requires the
        // eigensolve to have run before it can sample the mode at port-
        // edge midpoints. Without a cached `mode_profile` the RHS
        // degenerates to all-zeros — programmer-error path (the test
        // covers the documented zero-fallback). The previous "uniform
        // fallback" behavior was a Phase 1.3.1.1 step 0-1 stub that
        // existed only because the eigensolve was unimplemented.
        // mom-001 / Phase 1.3.0 paths are untouched: they use
        // `ModalDistribution::Uniform` or `DeltaGapPort` directly.
        let basis = RwgBasis::from_mesh(two_tri_mesh_with_port()).unwrap();
        let numerical = WavePort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
            mode_phase_velocity_factor: 1.0,
            modal_distribution: ModalDistribution::Numerical2D(Box::new(
                unit_square_cross_section(),
            )),
        };
        let b_numerical = numerical.rhs(&basis, 1.0e9);
        let n = basis.n_basis();
        for k in 0..n {
            assert!(
                b_numerical[(k, 0)].norm() < 1e-15,
                "pre-solve Numerical2D RHS must be zero at k={k}, got {:?}",
                b_numerical[(k, 0)]
            );
        }
    }

    #[test]
    fn tem_smoothed_port_degenerate_width_equals_uniform_rescale() {
        // With zero strip_width, the Maxwell envelope collapses to a
        // uniform `√(2/π)` and the port acts as a delta-gap port at a
        // rescaled voltage. The Z_in = V / I_port computation is
        // invariant under a uniform RHS rescale (the same factor
        // multiplies b and is removed by the symmetric port_current
        // extraction), so the answer must equal the delta-gap result
        // bit-for-bit on a single-column wire mesh.
        let basis = RwgBasis::from_mesh(two_tri_mesh_with_port()).unwrap();
        let dg = DeltaGapPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
        };
        let tem = TemSmoothedPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
            strip_width_m: 0.0,
        };
        let b_dg = dg.rhs(&basis, 1.0e9);
        let b_tem = tem.rhs(&basis, 1.0e9);
        let scale = (2.0 / std::f64::consts::PI).sqrt();
        let n = basis.n_basis();
        for k in 0..n {
            let expected = b_dg[(k, 0)] * Complex64::new(scale, 0.0);
            assert!(
                (b_tem[(k, 0)] - expected).norm() < 1e-12,
                "k={k}: tem b = {:?}, expected scale·b_dg = {:?}",
                b_tem[(k, 0)],
                expected
            );
        }
    }

    #[test]
    fn tem_smoothed_port_weight_peaks_at_strip_edges() {
        // The Maxwell envelope `1/√(1−(2y/w)²)` is minimum at y = 0
        // (the strip centre) and diverges at y → ±w/2. Sanity-check
        // the monotone behaviour and the clipped finite limit.
        let port = TemSmoothedPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
            strip_width_m: 1.0,
        };
        let w0 = port.weight(0.0);
        let w_mid = port.weight(0.25);
        let w_edge = port.weight(0.5);
        assert!(w0 > 0.0);
        assert!(w_mid > w0);
        assert!(w_edge > w_mid);
        // ε = 1e-3 clip means w_edge is finite — explicit bound is
        // sqrt(2/π · 1/(2ε − ε²)) ≈ sqrt(1/π · 1/ε) ≈ 17.8
        assert!(w_edge.is_finite());
        assert!(w_edge < 1e2, "w_edge = {w_edge} should be clipped");
    }

    #[test]
    fn wave_port_with_numerical_cross_section_builder_sets_variant() {
        let wp = WavePort {
            tag: 7,
            voltage: Complex64::new(1.0, 0.0),
            mode_phase_velocity_factor: 1.0,
            modal_distribution: ModalDistribution::Uniform,
        }
        .with_numerical_cross_section(unit_square_cross_section());
        assert!(matches!(
            wp.modal_distribution,
            ModalDistribution::Numerical2D(_)
        ));
        assert_eq!(wp.tag, 7);
    }

    /// Build a coarse **shielded-microstrip** cross-section in the
    /// strip-as-hole PEC style, matching the construction used in
    /// `tests/eigensolver_microstrip_quasi_tem.rs`.
    ///
    /// The outer PEC box spans `[0, wb] × [0, hb]` with an `nx × ny`
    /// structured grid. Cells whose centroid falls inside the signal-strip
    /// rectangle `(x ∈ [xc−w/2, xc+w/2], y ∈ [h_sub, h_sub+t])` are
    /// **omitted** so that their boundary edges become mesh-boundary = PEC
    /// Dirichlet — forming the inner PEC signal-strip conductor. FR-4
    /// substrate fills `y < h_sub` (tag 1, ε_r = 4.4); air above (tag 0).
    ///
    /// Grid kept to 8 × 8 (≈ 100–120 DoF after hole) so the `O(n³)` dense
    /// eigensolve is fast in a unit test context.
    fn shielded_microstrip_ports_test(nx: usize, ny: usize) -> NumericalCrossSection {
        // 50 Ω FR-4 geometry: w/h ≈ 1.9, h_sub = 1.6 mm.
        let eps_r_sub = 4.4_f64;
        let h_sub = 1.6e-3_f64;
        let w_strip = 1.9 * h_sub;
        let wb = 8.0 * w_strip;
        let hb = 7.0 * h_sub;
        let xs: Vec<f64> = (0..=nx).map(|i| wb * (i as f64) / (nx as f64)).collect();
        let ys: Vec<f64> = (0..=ny).map(|j| hb * (j as f64) / (ny as f64)).collect();
        let xc = wb / 2.0;
        let (sx0, sx1) = (xc - w_strip / 2.0, xc + w_strip / 2.0);
        // Strip thickness ≈ 1 cell so the strip boundary has at least one
        // cell of interior edges (the hole border becomes PEC automatically).
        let t_strip = hb / (ny as f64);
        let (sy0, sy1) = (h_sub, h_sub + t_strip);
        let in_strip = |cx: f64, cy: f64| {
            cx > sx0 - 1e-12 && cx < sx1 + 1e-12 && cy > sy0 - 1e-12 && cy < sy1 + 1e-12
        };
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
        for &y in &ys {
            for &x in &xs {
                vertices.push([x, y]);
            }
        }
        let mut triangles = Vec::new();
        let mut tags: Vec<u32> = Vec::new();
        for j in 0..ny {
            let yc = 0.5 * (ys[j] + ys[j + 1]);
            for i in 0..nx {
                let xcell = 0.5 * (xs[i] + xs[i + 1]);
                if in_strip(xcell, yc) {
                    continue; // omit strip cells → their edges become PEC
                }
                let v00 = idx(i, j);
                let v10 = idx(i + 1, j);
                let v11 = idx(i + 1, j + 1);
                let v01 = idx(i, j + 1);
                let tag = if yc < h_sub { 1u32 } else { 0u32 };
                triangles.push([v00, v10, v11]);
                tags.push(tag);
                triangles.push([v00, v11, v01]);
                tags.push(tag);
            }
        }
        let mesh = TriMesh2D::new(vertices, triangles, None, Some(tags))
            .expect("shielded microstrip mesh invariants");
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        eps.insert(1u32, Complex64::new(eps_r_sub, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        mu.insert(1u32, Complex64::new(1.0, 0.0));
        NumericalCrossSection::new(mesh, eps, mu)
    }

    /// Phase 1.3.1.2 + B (ADR-0060): `with_quasi_tem()` builder sets the flag
    /// and `solve()` on a shielded-microstrip cross-section succeeds,
    /// populates `mode_profile`, and `e_tangential_at` at an interior point
    /// returns a finite non-NaN `[f64; 2]`.
    ///
    /// Grid: 12 × 6 (≈ 130 DoF) — coarse enough for a fast unit test but
    /// sufficient for the shift-invert ladder to resolve the quasi-TEM mode
    /// (the strip-as-hole spans 3–4 cells in x at 12 columns across 8 strip
    /// widths). The validated gate uses 20 × 10; this is the minimal scale.
    #[test]
    fn with_quasi_tem_builder_sets_flag_and_solve_succeeds() {
        // Confirm the flag-defaulting (closed-guide default = false).
        let m = shielded_microstrip_ports_test(12, 6);
        assert!(!m.use_quasi_tem, "use_quasi_tem must default to false");
        // Confirm the builder flips it.
        let m2 = shielded_microstrip_ports_test(12, 6).with_quasi_tem();
        assert!(m2.use_quasi_tem, "with_quasi_tem() must set the flag");

        // Smoke: quasi-TEM solve on the shielded-microstrip cross-section
        // must succeed and populate mode_profile. At 2 GHz the box modes are
        // all below cutoff, so the quasi-TEM path finds the dominant
        // zero-cutoff mode.
        let mut m3 = shielded_microstrip_ports_test(12, 6).with_quasi_tem();
        let res = m3.solve(2.0e9);
        assert!(
            res.is_ok(),
            "with_quasi_tem() solve must succeed on shielded microstrip, got: {res:?}"
        );
        assert!(
            m3.mode_profile.is_some(),
            "mode_profile must be populated after successful quasi-TEM solve"
        );
        // e_tangential_at at a representative interior point (substrate midpoint)
        // must return a finite non-NaN pair — the Nedelec interpolation ran.
        let h_sub = 1.6e-3_f64;
        let w_strip = 1.9 * h_sub;
        let wb = 8.0 * w_strip;
        let xc = wb / 2.0;
        let et = m3.e_tangential_at(xc, h_sub / 2.0);
        assert!(
            et[0].is_finite() && et[1].is_finite(),
            "e_tangential_at must return finite values, got {et:?}"
        );
        // At least one component must be non-zero (the mode is not trivially
        // zero inside the domain for a physical quasi-TEM).
        let norm = (et[0] * et[0] + et[1] * et[1]).sqrt();
        assert!(
            norm > 0.0,
            "e_tangential_at must be non-zero at interior point"
        );
    }

    /// Phase 1.3.1.2 + B: `with_quasi_tem()` + `with_element_order(Second)`
    /// returns `Unsupported` — the p=2 quasi-TEM combination is explicitly
    /// out of scope (p=1 only for the validated excitation path).
    #[test]
    fn quasi_tem_plus_second_order_returns_unsupported() {
        let mut m = shielded_microstrip_ports_test(4, 4)
            .with_quasi_tem()
            .with_element_order(ElementOrder::Second);
        match m.solve(2.0e9) {
            Err(yee_core::Error::Unimplemented(_)) => { /* expected */ }
            other => panic!("quasi-TEM + Second-order must return Unimplemented, got: {other:?}"),
        }
    }
}
