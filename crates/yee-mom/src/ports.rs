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

/// Modal distribution attached to a [`WavePort`].
///
/// Phase 1.3.1.0 ships [`ModalDistribution::Te10`] for rectangular
/// waveguides; the default [`ModalDistribution::Uniform`] preserves the
/// Phase 1.3.0 delta-gap-equivalent placeholder behaviour. A future
/// `Numerical2D` variant will hold the result of a 2-D cross-section
/// eigensolve once Phase 1.3.1.1 lands.
pub enum ModalDistribution {
    /// Uniform (TEM-like) dominant-mode amplitude across the port edges.
    /// Galerkin-tested into edge-length-weighted RWG basis functions this
    /// reduces to the delta-gap form, preserving bit-for-bit equivalence
    /// with [`DeltaGapPort`] for the mom-001 gate.
    Uniform,
    /// Closed-form TE10 mode of a rectangular waveguide cross-section.
    /// See [`RectangularWaveguideTe10`].
    Te10(RectangularWaveguideTe10),
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
}
