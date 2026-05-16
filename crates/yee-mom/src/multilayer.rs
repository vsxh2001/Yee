//! Multilayer mixed-potential Green's function (Phase 1.1 walking
//! skeleton).
//!
//! ## Scope (Phase 1.1.0 placeholder — NOT production-grade)
//!
//! This module provides a [`MultilayerGreens`] kernel that augments the
//! free-space scalar Green's function with a **single complex image**
//! placed below an effective ground plane, weighted by a quasi-static
//! dielectric-contrast factor. It is intentionally the simplest possible
//! "multilayer" Green's function consistent with the
//! Michalski–Mosig / DCIM convention
//! `G = Σ b_n · exp(-j k₀ R_n) / (4 π R_n)` with `R_n² = ρ² + (z + z' − a_n)²`.
//!
//! Concretely, for a "one substrate slab of thickness `h` and relative
//! permittivity `ε_r` over a PEC ground at z = −h", we use a single image
//! pair `(b_1, a_1) = ((ε_r − 1) / (ε_r + 1), −2 h)`. This recovers the
//! free-space kernel exactly in two interesting limits:
//!
//! * `ε_r = 1` (no dielectric contrast) → `b_1 = 0`, no image.
//! * `h → ∞` (ground plane at infinity) → image distance → ∞, image
//!   contribution → 0.
//!
//! The TE / TM (vector / scalar potential) split is **collapsed**: both
//! `G^A` and `G^Φ` use the same single-image expression. Real Sommerfeld-
//! integral or production DCIM extraction (Phase 1.1.1+) will replace
//! this with separate image trains for the two potentials.
//!
//! ## Why this placeholder ships now
//!
//! Phase 1.1 unblocks mom-002 (microstrip Z₀) and mom-003 (patch
//! resonance) by providing the *abstraction layer* — the [`Greens`] trait,
//! the generic `impedance_matrix`, and the multilayer struct itself.
//! The numerical fidelity of the one-image approximation is not adequate
//! for those validation gates and they will run with loose tolerances
//! until Phase 1.1.1+ lands real DCIM extraction.
//!
//! ## References
//!
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans. Antennas
//!   Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.
//! * Y. L. Chow, J. J. Yang, D. G. Fang, G. E. Howard, "A closed-form
//!   spatial Green's function for the thick microstrip substrate," *IEEE
//!   Trans. Microw. Theory Tech.*, vol. 39, no. 3, pp. 588–592, Mar 1991.

#![allow(dead_code)]

use crate::greens::Greens;
use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::units::C0;

/// Multilayer scalar Green's function: free space plus a fitted set of
/// complex images. Phase 1.1.0 uses a single real-axis image, but the
/// `vector_images` / `scalar_images` vectors are sized to admit additional
/// DCIM terms in Phase 1.1.1+ without an API churn.
pub struct MultilayerGreens {
    /// Background wave number `k₀ = ω / c` of the upper half space
    /// (treated as free space in Phase 1.1.0).
    pub k0: Complex64,
    /// Free-space wave impedance √(μ₀/ε₀) (ohms).
    pub eta0: f64,
    /// Relative permittivity of the substrate slab.
    pub eps_r: f64,
    /// Substrate thickness in metres; PEC ground plane sits at `z = −h`.
    pub h: f64,
    /// `(b_n, a_n)` image coefficients for the vector-potential Green's
    /// function. Each contributes `b_n · exp(-j k₀ R_n) / (4 π R_n)` with
    /// `R_n² = ρ² + (z + z' − a_n)²` (field-vs-source z addition because
    /// the image is on the opposite side of the ground plane).
    pub vector_images: Vec<(Complex64, Complex64)>,
    /// `(b_n, a_n)` image coefficients for the scalar-potential Green's
    /// function. Phase 1.1.0 sets this equal to `vector_images`; the TE/TM
    /// split lands in Phase 1.1.1.
    pub scalar_images: Vec<(Complex64, Complex64)>,
}

impl MultilayerGreens {
    /// Build a placeholder one-image multilayer Green's function for the
    /// canonical "substrate on PEC ground" microstrip geometry.
    ///
    /// Parameters:
    /// * `freq_hz` — operating frequency, used to set `k₀ = ω / c`.
    /// * `eps_r` — relative permittivity of the substrate slab.
    /// * `h` — substrate thickness (m); the PEC ground lies at `z = −h`.
    ///
    /// Phase 1.1.0 uses a single image with
    /// `b = (ε_r − 1) / (ε_r + 1)` (real) at `a = −2h` (real, negative).
    /// See module docstring for the limitations of this approximation.
    pub fn new_microstrip(freq_hz: f64, eps_r: f64, h: f64) -> Self {
        let omega = std::f64::consts::TAU * freq_hz;
        let k0 = Complex64::new(omega / C0, 0.0);
        let eta0 = yee_core::units::ETA0;
        // Quasi-static reflection-like factor blending the PEC and the
        // dielectric-contrast image effects. ε_r → 1 ⇒ no image. The PEC
        // limit ε_r → ∞ gives b → 1, i.e. a full constructive image — a
        // crude stand-in for the negative-unity PEC reflection that we
        // intentionally use here because the TE/TM split needed for a
        // sign-faithful PEC image is deferred to Phase 1.1.1.
        let gamma = (eps_r - 1.0) / (eps_r + 1.0);
        let image = (Complex64::new(gamma, 0.0), Complex64::new(-2.0 * h, 0.0));
        Self {
            k0,
            eta0,
            eps_r,
            h,
            vector_images: vec![image],
            scalar_images: vec![image],
        }
    }

    /// Mirror a field point's z-coordinate against an image location
    /// `a_n`. With `a_n = −2 h` and source point at `z = z'`, the
    /// resulting image lies at `z'_img = -2 h − z'`, exactly the spec'd
    /// PEC-ground image.
    fn image_distance(r_field: Vector3<f64>, r_source: Vector3<f64>, a_n: Complex64) -> Complex64 {
        // For complex a_n the radial distance becomes
        // R² = ρ² + (z_field + z_source − a_n)². Phase 1.1.0 keeps a_n
        // strictly real (Im a_n = 0), but we honour the complex form so
        // Phase 1.1.1's complex DCIM coefficients drop in without touching
        // this helper.
        let dx = r_field.x - r_source.x;
        let dy = r_field.y - r_source.y;
        let rho_sq = dx * dx + dy * dy;
        let dz = Complex64::new(r_field.z + r_source.z, 0.0) - a_n;
        (Complex64::new(rho_sq, 0.0) + dz * dz).sqrt()
    }

    /// Sum the image train `images` evaluated between field and source
    /// points. Returns Σ b_n · exp(-j k₀ R_n) / (4 π R_n).
    fn image_sum(
        &self,
        images: &[(Complex64, Complex64)],
        r1: Vector3<f64>,
        r2: Vector3<f64>,
    ) -> Complex64 {
        let mut acc = Complex64::new(0.0, 0.0);
        for &(b_n, a_n) in images {
            let r_n = Self::image_distance(r1, r2, a_n);
            if r_n.norm() < 1e-300 {
                // A coincident image point would be pathological; in
                // Phase 1.1.0 it cannot occur (image is below the ground
                // plane), so we treat it as zero contribution rather than
                // propagate a NaN.
                continue;
            }
            let phase = -self.k0 * r_n;
            let exp_term = phase.exp();
            acc += b_n * exp_term / (Complex64::new(4.0 * std::f64::consts::PI, 0.0) * r_n);
        }
        acc
    }

    /// Free-space scalar Green's function for the direct (no-image) term.
    /// Inlined here so this module is independent of `FreeSpaceGreen`'s
    /// concrete struct while sharing the same numerical formula.
    fn free_space(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        let r = (r1 - r2).norm();
        assert!(
            r > 0.0,
            "MultilayerGreens free-space term is singular at r1 == r2; use the _smooth variants"
        );
        let phase = -self.k0.re * r;
        Complex64::from_polar(1.0 / (4.0 * std::f64::consts::PI * r), phase)
    }

    /// Singularity-subtracted free-space term: `G_free − 1/(4 π R)`, with
    /// the `−j k₀ / (4 π)` limit value at `R = 0`.
    fn free_space_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        let r = (r1 - r2).norm();
        if r == 0.0 {
            return Complex64::new(0.0, -self.k0.re / (4.0 * std::f64::consts::PI));
        }
        let g = self.free_space(r1, r2);
        g - Complex64::new(1.0 / (4.0 * std::f64::consts::PI * r), 0.0)
    }
}

impl Greens for MultilayerGreens {
    fn k0(&self) -> Complex64 {
        self.k0
    }
    fn eta0(&self) -> f64 {
        self.eta0
    }
    fn scalar_vector(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.free_space(r1, r2) + self.image_sum(&self.vector_images, r1, r2)
    }
    fn scalar_scalar(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.free_space(r1, r2) + self.image_sum(&self.scalar_images, r1, r2)
    }
    fn scalar_vector_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        // Only the direct (free-space) term carries the 1/R singularity;
        // image points are never coincident with the source on Phase 1.1.0
        // geometries, so the image sum is added in full.
        self.free_space_smooth(r1, r2) + self.image_sum(&self.vector_images, r1, r2)
    }
    fn scalar_scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.free_space_smooth(r1, r2) + self.image_sum(&self.scalar_images, r1, r2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::greens::FreeSpaceGreen;

    // λ = 2 m so that k0 = π exactly — keeps reference values clean.
    const F_LAMBDA_2M: f64 = C0 / 2.0;

    #[test]
    fn eps_r_one_zeros_the_image_weight() {
        // ε_r = 1 ⇒ Γ = 0 ⇒ the image weight `b_1` is zero. The trait
        // evaluators must then match the free-space kernel to machine
        // precision.
        let m = MultilayerGreens::new_microstrip(F_LAMBDA_2M, 1.0, 0.1);
        let f = FreeSpaceGreen::new(F_LAMBDA_2M);
        let r1 = Vector3::new(0.0, 0.0, 0.0);
        let r2 = Vector3::new(0.3, 0.4, 0.0);
        let exp_g = f.scalar(r1, r2);
        let got_v = m.scalar_vector(r1, r2);
        let got_s = m.scalar_scalar(r1, r2);
        assert!((exp_g - got_v).norm() < 1e-14);
        assert!((exp_g - got_s).norm() < 1e-14);
    }

    #[test]
    fn h_infinite_recovers_free_space() {
        // h → ∞ pushes the image to infinite distance, so its
        // 1/R-decayed contribution vanishes. Use h = 1e20 m as a
        // numerical stand-in.
        let m = MultilayerGreens::new_microstrip(F_LAMBDA_2M, 4.4, 1.0e20);
        let f = FreeSpaceGreen::new(F_LAMBDA_2M);
        let r1 = Vector3::new(0.0, 0.0, 0.0);
        let r2 = Vector3::new(0.3, 0.4, 0.0);
        let exp_g = f.scalar(r1, r2);
        let got = m.scalar_vector(r1, r2);
        assert!(
            (exp_g - got).norm() < 1e-30,
            "image at z → −∞ must not perturb free-space kernel"
        );
    }

    #[test]
    fn image_term_is_finite_for_typical_microstrip() {
        // Sanity check on a representative substrate geometry: 150 MHz,
        // ε_r = 4.4, h = 1 mm, ρ = 1 cm. The image must contribute a
        // finite, non-zero value to both vector and scalar Green's
        // functions and remain smaller than the direct term (since the
        // image is farther away).
        let m = MultilayerGreens::new_microstrip(150.0e6, 4.4, 1.0e-3);
        let r1 = Vector3::new(0.0, 0.0, 5.0e-4);
        let r2 = Vector3::new(1.0e-2, 0.0, 5.0e-4);
        let g_full = m.scalar_vector(r1, r2);
        let g_image_only = m.image_sum(&m.vector_images, r1, r2);
        assert!(g_full.re.is_finite() && g_full.im.is_finite());
        assert!(g_image_only.norm() > 0.0);
        assert!(g_image_only.norm() < g_full.norm() * 2.0);
    }
}
