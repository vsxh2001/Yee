//! Free-space scalar Green's function, singularity-subtracted form, and
//! the [`Greens`] trait abstraction over Green's-function flavours.
//!
//! The [`Greens`] trait splits the mixed-potential integral equation (MPIE)
//! kernel into two scalar Green's functions:
//!
//! * `G^A(r, r')` — used in the magnetic-vector-potential term
//!   `j ω μ₀ ⟨f_m, A_n⟩`. In free space this is the canonical
//!   `exp(-j k₀ R) / (4 π R)`. In multilayer media it splits into the TE
//!   (parallel) component carrying the tangential-current contribution.
//! * `G^Φ(r, r')` — used in the scalar-potential term
//!   `(1 / (j ω ε₀)) ⟨∇·f_m, φ_n⟩`. In free space this is identical to
//!   `G^A`; in multilayer media it carries the TM (perpendicular) component
//!   tied to surface charge.
//!
//! `FreeSpaceGreen` collapses both onto the same scalar kernel; the
//! distinction only matters once `MultilayerGreens` (sibling module
//! `multilayer`) lands.
//!
//! Dead-code allowances: this module landed in Phase 1.0 Task 5 ahead of
//! the impedance fill that consumes it (Task 8); Phase 1.1 extends it with
//! the `Greens` trait. Clippy with `-D warnings` flags every yet-unused
//! symbol as an error, so the struct and its associated items are
//! explicitly tagged `#[allow(dead_code)]` at the module boundary.
#![allow(dead_code)]

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::units::C0;

/// Abstraction over the MPIE Green's-function kernel.
///
/// Implementors expose the wave number `k0` and characteristic impedance
/// `eta0` of the background medium, plus four scalar Green's evaluations:
/// the full vector/scalar pair and their singularity-subtracted (`_smooth`)
/// variants suitable for use inside Duffy regularisation at coincident or
/// nearly-coincident integration points.
///
/// For free space the vector and scalar variants are identical; for the
/// multilayer kernels they differ because the TE / TM partitioning of the
/// substrate Sommerfeld integrals (or their DCIM approximation) projects
/// onto each potential differently.
pub trait Greens {
    /// Background wave number (Hz·s/m → rad/m once multiplied by ω).
    fn k0(&self) -> Complex64;
    /// Free-space wave impedance √(μ₀/ε₀), in ohms.
    fn eta0(&self) -> f64;
    /// Scalar Green's function `G^A(r, r')` for the vector-potential term
    /// of the MPIE. For free space this equals `exp(-j k₀ R) / (4 π R)`.
    fn scalar_vector(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64;
    /// Scalar Green's function `G^Φ(r, r')` for the scalar-potential term.
    /// For free space this also equals `exp(-j k₀ R) / (4 π R)`; the two
    /// only differ in multilayer media where TE and TM modes split.
    fn scalar_scalar(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64;
    /// Singularity-subtracted form of [`Greens::scalar_vector`], finite at
    /// `R = 0`. Provided for use as a graceful fallback inside Duffy
    /// quadrature when an inner Gauss point coincides bit-exactly with the
    /// outer anchor (the Jacobian vanishes there, so the value is weighted
    /// by zero — only its finiteness matters).
    fn scalar_vector_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64;
    /// Singularity-subtracted form of [`Greens::scalar_scalar`], finite at
    /// `R = 0`. See [`Greens::scalar_vector_smooth`] for usage notes.
    fn scalar_scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64;
}

/// Free-space scalar Green's function kernel. The MPIE vector- and
/// scalar-potential terms collapse onto a single `exp(-j k₀ R) / (4 π R)`
/// here; see the [`Greens`] trait for the multilayer split.
pub struct FreeSpaceGreen {
    /// Background wave number `k₀ = ω / c` (real-valued for lossless free
    /// space, but stored as `Complex64` so the same arithmetic kernels can
    /// drive lossy / multilayer extensions).
    pub k0: Complex64,
    /// Free-space wave impedance √(μ₀/ε₀) in ohms.
    pub eta0: f64,
}

impl FreeSpaceGreen {
    /// Construct the free-space kernel at `freq_hz`, computing
    /// `k₀ = 2 π f / c₀` and pulling the canonical `η₀` from
    /// [`yee_core::units`].
    pub fn new(freq_hz: f64) -> Self {
        let omega = std::f64::consts::TAU * freq_hz;
        let k0 = Complex64::new(omega / C0, 0.0);
        let eta0 = yee_core::units::ETA0;
        Self { k0, eta0 }
    }

    /// Scalar Green's function G(R) = exp(-j k0 R) / (4 π R).
    ///
    /// # Panics
    ///
    /// Panics if `r1 == r2` (singular). Use `scalar_smooth` for the
    /// singularity-subtracted form valid at R = 0.
    pub fn scalar(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        let r = (r1 - r2).norm();
        assert!(
            r > 0.0,
            "scalar Green's function singular at r1 == r2; use scalar_smooth"
        );
        let k0 = self.k0.re;
        Complex64::from_polar(1.0 / (4.0 * std::f64::consts::PI * r), -k0 * r)
    }

    /// Singularity-subtracted scalar Green's function: G - 1/(4 π R).
    /// Finite at R = 0 — returns `-j k0 / (4 π)` in that limit.
    pub fn scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        let r = (r1 - r2).norm();
        if r == 0.0 {
            return Complex64::new(0.0, -self.k0.re / (4.0 * std::f64::consts::PI));
        }
        let g = self.scalar(r1, r2);
        g - Complex64::new(1.0 / (4.0 * std::f64::consts::PI * r), 0.0)
    }
}

impl Greens for FreeSpaceGreen {
    fn k0(&self) -> Complex64 {
        self.k0
    }
    fn eta0(&self) -> f64 {
        self.eta0
    }
    fn scalar_vector(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.scalar(r1, r2)
    }
    fn scalar_scalar(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.scalar(r1, r2)
    }
    fn scalar_vector_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.scalar_smooth(r1, r2)
    }
    fn scalar_scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.scalar_smooth(r1, r2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Use a frequency derived from the exact C0 constant so that λ = 2 m and
    // k0 = π hold exactly — this keeps the assertion tolerances tight while
    // remaining faithful to the documented λ/k0 relationship.
    const F_LAMBDA_2M: f64 = C0 / 2.0;

    #[test]
    fn wave_number_from_frequency() {
        let g = FreeSpaceGreen::new(F_LAMBDA_2M);
        // λ = c / f = 2 m ; k0 = 2π/λ = π
        let expected = std::f64::consts::PI;
        assert!((g.k0.re - expected).abs() < 1e-9 * expected);
    }

    #[test]
    fn scalar_amplitude_at_quarter_wavelength() {
        let g = FreeSpaceGreen::new(F_LAMBDA_2M); // λ = 2 m
        let r1 = Vector3::new(0.0, 0.0, 0.0);
        let r2 = Vector3::new(0.5, 0.0, 0.0); // R = λ/4
        let v = g.scalar(r1, r2);
        let expected_mag = 1.0 / (4.0 * std::f64::consts::PI * 0.5);
        assert!((v.norm() - expected_mag).abs() < 1e-12 * expected_mag);
        // phase = -k0 R = -π/2
        let expected_phase = -std::f64::consts::FRAC_PI_2;
        assert!((v.arg() - expected_phase).abs() < 1e-12);
    }

    #[test]
    fn scalar_smooth_limit_at_zero() {
        let g = FreeSpaceGreen::new(F_LAMBDA_2M);
        let r = Vector3::new(0.0, 0.0, 0.0);
        let v = g.scalar_smooth(r, r);
        let expected = Complex64::new(0.0, -g.k0.re / (4.0 * std::f64::consts::PI));
        assert!((v - expected).norm() < 1e-12);
    }

    #[test]
    fn free_space_trait_matches_inherent() {
        // The Greens trait impl for FreeSpaceGreen must dispatch to the
        // same scalar G(R) for both vector and scalar potentials.
        let g = FreeSpaceGreen::new(F_LAMBDA_2M);
        let r1 = Vector3::new(0.0, 0.0, 0.0);
        let r2 = Vector3::new(0.3, 0.4, 0.0);
        let direct = g.scalar(r1, r2);
        let via_v = <FreeSpaceGreen as Greens>::scalar_vector(&g, r1, r2);
        let via_s = <FreeSpaceGreen as Greens>::scalar_scalar(&g, r1, r2);
        assert!((direct - via_v).norm() < 1e-14);
        assert!((direct - via_s).norm() < 1e-14);
        assert_eq!(<FreeSpaceGreen as Greens>::k0(&g), g.k0);
        assert_eq!(<FreeSpaceGreen as Greens>::eta0(&g), g.eta0);
    }
}
