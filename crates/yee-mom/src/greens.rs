//! Free-space scalar Green's function and singularity-subtracted form.
//!
//! Dead-code allowances: this module lands in Phase 1.0 Task 5 ahead of the
//! impedance fill that consumes it (Task 8). Clippy with `-D warnings` flags
//! every yet-unused symbol as an error, so the struct and its associated
//! items are explicitly tagged `#[allow(dead_code)]` at the module boundary.
//! The allow will be removed implicitly once `fill.rs` references these
//! items.
#![allow(dead_code)]

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::units::C0;

pub(crate) struct FreeSpaceGreen {
    pub k0: Complex64,
    pub eta0: f64,
}

impl FreeSpaceGreen {
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
}
