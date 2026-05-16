//! Port abstractions for planar MoM.
//!
//! Phase 1.3 walking skeleton: a [`Port`] trait plus two implementations.
//!
//! * [`DeltaGapPort`] — preserves the Phase 1.0 / mom-001 behaviour
//!   bit-for-bit: a 1 V (or user-supplied) impulse across every RWG edge
//!   tagged with `tag`, Galerkin-tested into `b[k] = V × length_k`.
//! * [`WavePort`] — **API placeholder only** for Phase 1.3.0. The
//!   modal-distribution logic that distinguishes a real wave port from a
//!   delta-gap (a non-uniform mode amplitude across the port cross-section,
//!   produced by a 1D eigenmode solve on the port edges) is deferred to
//!   Phase 1.3.1. At Phase 1.3.0 the [`WavePort::rhs`] is identical to a
//!   [`DeltaGapPort`] at the same voltage and tag — i.e. a uniform mode on
//!   a TEM-like cross-section with `mode_phase_velocity_factor = 1.0`. The
//!   `mode_phase_velocity_factor` field and the frequency-dependent
//!   propagation constant `β = ω / (c₀ · v_factor)` are wired up but not yet
//!   applied to the RHS distribution; they exist so the call sites and the
//!   trait surface settle in 1.3.0 and only the internal eigenmode plumbing
//!   needs to change in 1.3.1.
//!
//! The trait is `pub(crate)` because Phase 1.3 only exposes ports through
//! the high-level [`crate::PlanarMoM::run`] entry point. A public surface
//! will follow once the wave-port modal solver lands and the API is stable.

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

/// Wave port — 1D modal source on a tagged edge set.
///
/// **Phase 1.3.0 status: API placeholder.** The RHS today is identical to a
/// [`DeltaGapPort`] at the same voltage and tag because the Phase 1.3.0
/// approximation is a uniform (TEM-like) dominant mode whose Galerkin
/// projection onto edge-length-weighted RWG basis functions reduces to the
/// delta-gap form. The full eigenmode-solve and non-uniform modal weighting
/// will land in Phase 1.3.1 (microstrip / CPW / waveguide cross-sections).
pub(crate) struct WavePort {
    /// Port tag — matches the mesh tagging scheme.
    pub tag: u32,
    /// Modal-source reference voltage. With a uniform mode at
    /// `mode_phase_velocity_factor = 1.0` this maps directly onto the
    /// equivalent delta-gap drive amplitude.
    pub voltage: Complex64,
    /// Phase-velocity factor for the lowest-order mode on the port cross
    /// section, as a fraction of `c₀`. Phase 1.3.0 uses `1.0` (free-space /
    /// TEM); Phase 1.3.1 will compute this from a 1D eigenmode solve and
    /// use it to build the modal field distribution applied to the RHS.
    pub mode_phase_velocity_factor: f64,
}

impl Port for WavePort {
    fn tag(&self) -> u32 {
        self.tag
    }
    fn rhs(&self, basis: &RwgBasis, freq_hz: f64) -> Mat<Complex64> {
        // Phase 1.3.0: the wave-port modal weighting differs from delta-gap
        // by distributing the source across port edges according to the
        // mode field. For a uniform mode (TEM dominant-mode approximation),
        // the distribution is uniform, so we degenerate to delta-gap with
        // a frequency-dependent mode-amplitude scaling derived from the
        // TEM phase velocity. The β term below is computed so the call
        // site is in place for Phase 1.3.1 but does not yet modulate the
        // returned RHS — preserving bit-for-bit equivalence with the
        // delta-gap path for the mom-001 gate.
        let n = basis.n_basis();
        let mut b = Mat::<Complex64>::zeros(n, 1);
        let omega = std::f64::consts::TAU * freq_hz;
        let beta = omega / (yee_core::units::C0 * self.mode_phase_velocity_factor);
        let _ = beta; // reserved for Phase 1.3.1 propagation-correction term
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
        };
        let b1 = dg.rhs(&basis, 1.0e9);
        let b2 = wp.rhs(&basis, 1.0e9);
        let n = basis.n_basis();
        for k in 0..n {
            assert!((b1[(k, 0)] - b2[(k, 0)]).norm() < 1e-15);
        }
    }
}
