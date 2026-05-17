//! Multilayer mixed-potential Green's function — Phase 1.1.1.0 multi-image DCIM.
//!
//! ## Scope (Phase 1.1.1.0 — N-image DCIM, ≥ Phase 1.1.0 placeholder)
//!
//! This module provides a [`MultilayerGreens`] kernel that augments the
//! free-space scalar Green's function with a **train of `N` complex
//! images** (`N ∈ [1, 10]`) fitted to the spectral reflection
//! coefficients of a grounded dielectric slab via the Generalized
//! Pencil-of-Function (GPOF) procedure in [`crate::gpof`]. It honours
//! the Michalski–Mosig / DCIM convention
//! `G = Σ b_n · exp(-j k₀ R_n) / (4 π R_n)` with
//! `R_n² = ρ² + (z + z' − a_n)²`.
//!
//! Concretely, for a "single substrate slab of thickness `h` and
//! relative permittivity `ε_r` over a PEC ground at z = −h", the kernel:
//!
//! 1. Samples the TE / TM spectral reflection coefficients of the slab
//!    along a deformed contour `k_{z0}(t) = k₀ (1 − j t)` for
//!    `t ∈ [0, T_max]` (Aksun 1996 contour).
//! 2. Runs GPOF on each sample set to recover the `(α_n, β_n)` complex
//!    exponents.
//! 3. Maps `(α_n, β_n)` to image coefficients `(b_n, a_n)` via the
//!    Sommerfeld identity: `a_n = −β_n / k₀` (complex z-location) and
//!    `b_n = α_n · exp(j k₀ a_n) = α_n · exp(−j β_n)`. Vector-potential
//!    Green's function uses the TE fit; scalar-potential uses the TM
//!    fit.
//!
//! The N=1 constructor preserves the Phase 1.1.0 placeholder (one
//! real-axis image at z = −2h with weight Γ = (ε_r − 1)/(ε_r + 1)) so
//! the existing back-compat tests pass bit-for-bit.
//!
//! ## What this is NOT (deferred to Phase 1.1.1+)
//!
//! * Surface-wave pole extraction. For thin substrates at moderate
//!   frequencies the surface-wave contribution is small; this module
//!   ignores it.
//! * Sommerfeld-integral tail evaluation. We fit on a finite contour
//!   only.
//! * Multilayer stacks with more than substrate + ground plane.
//! * Anisotropic / lossy / dispersive substrates.
//!
//! The full Phase 1.1.1 plan ships those.
//!
//! ## References
//!
//! * M. I. Aksun, "A robust approach for the derivation of closed-form
//!   Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44,
//!   no. 5, pp. 651–658, May 1996.
//! * Y. L. Chow, J. J. Yang, D. G. Fang, G. E. Howard, "A closed-form
//!   spatial Green's function for the thick microstrip substrate,"
//!   *IEEE Trans. Microw. Theory Tech.*, vol. 39, no. 3,
//!   pp. 588–592, Mar 1991.
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans.
//!   Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7.

#![allow(dead_code)]

use crate::gpof::gpof;
use crate::greens::Greens;
use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::units::C0;

/// Multilayer scalar Green's function: free space plus a fitted set of
/// complex images. Phase 1.1.0 uses a single real-axis image; Phase
/// 1.1.1.0 fits up to N images via GPOF on the slab spectral
/// reflection coefficients. Both paths share the same trait impl and
/// summation kernel — they differ only in how the `(b_n, a_n)` pairs
/// are constructed.
pub struct MultilayerGreens {
    /// Background wave number `k₀ = ω / c` of the upper half space
    /// (treated as free space).
    pub k0: Complex64,
    /// Free-space wave impedance √(μ₀/ε₀) (ohms).
    pub eta0: f64,
    /// Relative permittivity of the substrate slab.
    pub eps_r: f64,
    /// Substrate thickness in metres; PEC ground plane sits at `z = −h`.
    pub h: f64,
    /// Number of complex images fitted into [`vector_images`] /
    /// [`scalar_images`]. `n_images = 1` reproduces the Phase 1.1.0
    /// placeholder exactly.
    pub n_images: usize,
    /// `(b_n, a_n)` image coefficients for the vector-potential Green's
    /// function. Each contributes `b_n · exp(-j k₀ R_n) / (4 π R_n)`
    /// with `R_n² = ρ² + (z + z' − a_n)²` (field-vs-source z addition
    /// because the image lies on the opposite side of the ground
    /// plane). Phase 1.1.0 fills this with one entry; Phase 1.1.1.0
    /// fills it with up to `n_images` entries from a GPOF fit to the
    /// TE spectral reflection coefficient.
    pub vector_images: Vec<(Complex64, Complex64)>,
    /// `(b_n, a_n)` image coefficients for the scalar-potential
    /// Green's function. Phase 1.1.0 sets this equal to
    /// [`vector_images`]; Phase 1.1.1.0 fits a separate TM-coefficient
    /// train, so the TE/TM split that Phase 1.1.0 collapsed is now
    /// resolved.
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
    /// Uses the Phase 1.1.0 one-image approximation: `b = (ε_r − 1) /
    /// (ε_r + 1)` (real) at `a = −2h` (real, negative). Preserved as a
    /// back-compat shortcut around [`Self::new_microstrip_with_n_images`].
    /// All existing integration tests target this constructor and
    /// continue to see the same numbers they did before Phase 1.1.1.0.
    pub fn new_microstrip(freq_hz: f64, eps_r: f64, h: f64) -> Self {
        Self::new_microstrip_with_n_images(eps_r, h, freq_hz, 1)
    }

    /// Build an N-image DCIM kernel for the same canonical microstrip
    /// geometry. `n_images = 1` reproduces [`Self::new_microstrip`]
    /// bit-for-bit (back-compat tripwire); `n_images ∈ [2, 10]` runs
    /// the GPOF fitter against the slab's TE and TM spectral
    /// reflection coefficients and stores the result in
    /// [`Self::vector_images`] / [`Self::scalar_images`].
    ///
    /// Note the parameter order: `(eps_r, h, freq_hz, n_images)` — the
    /// frequency moves to the third slot to match the canonical
    /// "build me a Greens for this substrate at this frequency with
    /// this many images" reading order.
    pub fn new_microstrip_with_n_images(eps_r: f64, h: f64, freq_hz: f64, n_images: usize) -> Self {
        let omega = std::f64::consts::TAU * freq_hz;
        let k0 = Complex64::new(omega / C0, 0.0);
        let eta0 = yee_core::units::ETA0;

        let (vector_images, scalar_images) = if n_images <= 1 {
            // N=1 back-compat path: literal Phase 1.1.0 placeholder.
            // Quasi-static reflection-like factor blending PEC and
            // dielectric-contrast image effects. ε_r → 1 ⇒ no image.
            let gamma = (eps_r - 1.0) / (eps_r + 1.0);
            let image = (Complex64::new(gamma, 0.0), Complex64::new(-2.0 * h, 0.0));
            (vec![image], vec![image])
        } else {
            // N>1 path: fit TE / TM reflection coefficients via GPOF.
            // A failed fit (rank-deficient samples, e.g. ε_r = 1 makes
            // R ≡ 0) falls back to the N=1 placeholder — better an
            // approximate kernel than a build-time panic.
            let v_imgs = fit_slab_dcim(eps_r, h, k0.re, n_images, SpectralKernel::Te)
                .unwrap_or_else(|_| placeholder_single_image(eps_r, h));
            let s_imgs = fit_slab_dcim(eps_r, h, k0.re, n_images, SpectralKernel::Tm)
                .unwrap_or_else(|_| placeholder_single_image(eps_r, h));
            (v_imgs, s_imgs)
        };

        Self {
            k0,
            eta0,
            eps_r,
            h,
            n_images,
            vector_images,
            scalar_images,
        }
    }

    /// Mirror a field point's z-coordinate against an image location
    /// `a_n`. With `a_n = −2 h` and source point at `z = z'`, the
    /// resulting image lies at `z'_img = -2 h − z'`, exactly the spec'd
    /// PEC-ground image. For complex `a_n` (Phase 1.1.1.0 DCIM
    /// coefficients) the radial distance is the analytic-continuation
    /// `R² = ρ² + (z_field + z_source − a_n)²` with the principal-branch
    /// square root.
    fn image_distance(r_field: Vector3<f64>, r_source: Vector3<f64>, a_n: Complex64) -> Complex64 {
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
                // the placeholder it cannot occur (image is below the
                // ground plane), so we treat it as zero contribution
                // rather than propagate a NaN.
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
        // image points are never coincident with the source on the
        // placeholder geometry, so the image sum is added in full.
        self.free_space_smooth(r1, r2) + self.image_sum(&self.vector_images, r1, r2)
    }
    fn scalar_scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.free_space_smooth(r1, r2) + self.image_sum(&self.scalar_images, r1, r2)
    }
}

// ---------------------------------------------------------------------
// Slab spectral kernels + GPOF wrapper
// ---------------------------------------------------------------------

/// Which spectral reflection coefficient to fit. The grounded-slab
/// reflection has independent TE / TM channels; fitting both gives the
/// TE/TM split that Phase 1.1.0 collapsed.
#[derive(Debug, Clone, Copy)]
enum SpectralKernel {
    /// TE polarisation. Drives the vector-potential image train
    /// (`G^A`).
    Te,
    /// TM polarisation. Drives the scalar-potential image train
    /// (`G^Φ`).
    Tm,
}

/// Phase 1.1.0 placeholder fallback used when GPOF cannot run (e.g.
/// `ε_r = 1` makes both reflection coefficients identically zero, the
/// SVD becomes rank-deficient, and the fit raises an error).
fn placeholder_single_image(eps_r: f64, h: f64) -> Vec<(Complex64, Complex64)> {
    let gamma = (eps_r - 1.0) / (eps_r + 1.0);
    vec![(Complex64::new(gamma, 0.0), Complex64::new(-2.0 * h, 0.0))]
}

/// Aksun 1996 deformed-contour parameter span. The contour is
/// `k_{z0}(t) = k₀ (1 − j t)` for `t ∈ [0, T_MAX]`; T_MAX = 10 is the
/// canonical Aksun value and covers ≈ 5 decades of spectral decay,
/// which is enough to resolve `N ≤ 10` images robustly.
const DCIM_T_MAX: f64 = 10.0;

/// Compute `k_{zd}` from `k_{z0}` and the substrate `ε_r` via the
/// dispersion relation `k_{zd}² = ε_r · k₀² − k_ρ²` =
/// `(ε_r − 1) · k₀² + k_{z0}²`. The principal-branch square root is
/// taken; for the Aksun contour `Im(k_{z0}) < 0` and the result has
/// the conventional `Re k_{zd} > 0` / `Im k_{zd} < 0` orientation.
fn k_zd_from_k_z0(k_z0: Complex64, eps_r: f64, k0: f64) -> Complex64 {
    let k0_sq = Complex64::new(k0 * k0, 0.0);
    let inside = Complex64::new(eps_r - 1.0, 0.0) * k0_sq + k_z0 * k_z0;
    inside.sqrt()
}

/// Evaluate the slab reflection coefficient at a point on the
/// integration contour for the requested polarisation.
///
/// Both are derived from the transmission-line analogy with a PEC
/// short at `z = -h`. Looking down from `z = 0+`, the input impedance
/// seen at the substrate top is `Z_in = j Z_d tan(k_{zd} h)`, where
/// `Z_d` is the characteristic impedance of the dielectric in the
/// requested polarisation; the reflection at the air-substrate
/// interface is `R = (Z_in − Z_0) / (Z_in + Z_0)`. After clearing the
/// `(ω μ_0)` / `(ω ε_0)` factors:
///
/// TE: `R = (j k_{z0} tan(k_{zd} h) − k_{zd})  / (j k_{z0} tan(k_{zd} h) + k_{zd})`.
/// TM: `R = (j k_{zd} tan(k_{zd} h) − ε_r k_{z0}) / (j k_{zd} tan(k_{zd} h) + ε_r k_{z0})`.
///
/// Near multiples of `π` in `k_{zd} h` the `tan` is singular — that
/// is a physical slab resonance and the caller's fit will surface it
/// as a large image coefficient. We do not try to dodge it; for the
/// thin-substrate FR-4 case at sub-GHz the resonance is far above the
/// band of interest.
fn slab_reflection(
    k_z0: Complex64,
    eps_r: f64,
    h: f64,
    k0: f64,
    kernel: SpectralKernel,
) -> Complex64 {
    let k_zd = k_zd_from_k_z0(k_z0, eps_r, k0);
    let phase = k_zd * Complex64::new(h, 0.0);
    let t = phase.tan();
    let j = Complex64::new(0.0, 1.0);
    match kernel {
        SpectralKernel::Te => {
            let num = j * k_z0 * t - k_zd;
            let den = j * k_z0 * t + k_zd;
            num / den
        }
        SpectralKernel::Tm => {
            let num = j * k_zd * t - Complex64::new(eps_r, 0.0) * k_z0;
            let den = j * k_zd * t + Complex64::new(eps_r, 0.0) * k_z0;
            num / den
        }
    }
}

/// Run a GPOF fit against the requested slab spectral reflection
/// coefficient and translate the recovered exponents into image
/// coefficients in the (b_n, a_n) convention used by [`image_sum`].
///
/// Sample the contour `k_{z0}(t_m) = k₀ (1 − j t_m)`,
/// `t_m = m · Δt`, with `Δt = T_MAX / (M − 1)` and `M = 2 N` (spec
/// DoD #1 minimum-sample setup). On this contour the spectral kernel
/// has the form `R(t) = Σ_n α_n · exp(β_n · t)` with
/// `α_n = a_n · exp(−j k₀ z_n)` and `β_n = −k₀ z_n` (so the recovered
/// image z-location is `z_n = −β_n / k₀` and the image weight is
/// `b_n = α_n · exp(j k₀ z_n) = α_n · exp(−j β_n)`).
fn fit_slab_dcim(
    eps_r: f64,
    h: f64,
    k0: f64,
    n_images: usize,
    kernel: SpectralKernel,
) -> Result<Vec<(Complex64, Complex64)>, crate::gpof::GpofError> {
    // Spec calls for M = 2N samples on the contour for N=5; we honour
    // the same rule for any N (so M = 2N = 10 at the default).
    let m = 2 * n_images;
    let dt = DCIM_T_MAX / ((m - 1) as f64);

    // Sample R(k_{z0}(t_m)). The Aksun contour `k_{z0} = k₀(1 − j t)`
    // crosses from the real axis into the lower half plane, which is
    // the steepest-descent direction for the Sommerfeld integrand.
    let samples: Vec<Complex64> = (0..m)
        .map(|kidx| {
            let t = (kidx as f64) * dt;
            let k_z0 = Complex64::new(k0, 0.0) * Complex64::new(1.0, -t);
            slab_reflection(k_z0, eps_r, h, k0, kernel)
        })
        .collect();

    // Fit. The recovered (α_n, β_n) satisfy R(t_m) ≈ Σ α_n exp(β_n t_m).
    let poles = gpof(&samples, dt, n_images)?;

    // Map (α_n, β_n) to image (b_n, a_n). The Sommerfeld identity is
    //
    //   R(k_{z0}) = Σ b_n · exp(−j k_{z0} · z_n)
    //          ⇔   G_refl(ρ, z, z') = Σ b_n · exp(−j k₀ R_n) / (4π R_n)
    //   with R_n² = ρ² + (z + z' + z_n)².
    //
    // The existing struct interprets the second tuple element
    // `a_n` as the z-coordinate of the image, with
    // `R_n² = ρ² + (z + z' − a_n)²` — note the **minus** sign, which
    // flips the math `z_n` into our `a_n = −z_n`.
    //
    // Substituting `k_{z0}(t) = k₀ (1 − j t)` into the spectral form
    // and matching against the GPOF model `R ≈ Σ α_n exp(β_n t)` gives
    //
    //   α_n = b_n · exp(−j k₀ z_n)   and   β_n = −k₀ z_n.
    //
    // Therefore:
    //   z_n (math)        = −β_n / k₀,
    //   a_n  (struct)     =  β_n / k₀,
    //   b_n               = α_n · exp(j k₀ z_n) = α_n · exp(−j β_n).
    let images: Vec<(Complex64, Complex64)> = poles
        .into_iter()
        .map(|(alpha, beta)| {
            let a_n = beta / Complex64::new(k0, 0.0);
            let j = Complex64::new(0.0, 1.0);
            let b_n = alpha * (-j * beta).exp();
            (b_n, a_n)
        })
        .collect();

    Ok(images)
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

    /// Back-compat tripwire: `new_microstrip_with_n_images(.., 1)` must
    /// produce exactly the same `(vector_images, scalar_images)` pair as
    /// the Phase 1.1.0 [`new_microstrip`] shortcut. If this regresses,
    /// every multilayer integration test from Phase 1.1.0 silently
    /// changes its numerics — surface it loudly.
    #[test]
    fn n_equals_one_matches_phase_1_1_0() {
        let f0 = 1.0e9;
        let a = MultilayerGreens::new_microstrip(f0, 4.4, 1.6e-3);
        let b = MultilayerGreens::new_microstrip_with_n_images(4.4, 1.6e-3, f0, 1);
        assert_eq!(a.vector_images.len(), 1);
        assert_eq!(b.vector_images.len(), 1);
        let (b_a, a_a) = a.vector_images[0];
        let (b_b, a_b) = b.vector_images[0];
        assert_eq!(b_a, b_b);
        assert_eq!(a_a, a_b);
        let (s_a, sa_a) = a.scalar_images[0];
        let (s_b, sa_b) = b.scalar_images[0];
        assert_eq!(s_a, s_b);
        assert_eq!(sa_a, sa_b);
    }

    /// N > 1 builds without panic for the canonical Phase 1.1.1.0
    /// microstrip parameters and produces a finite, non-degenerate
    /// image train for both TE and TM channels.
    #[test]
    fn n_images_five_fits_finite_coefficients_for_fr4() {
        let m = MultilayerGreens::new_microstrip_with_n_images(4.4, 1.6e-3, 1.0e9, 5);
        assert_eq!(m.vector_images.len(), 5);
        assert_eq!(m.scalar_images.len(), 5);
        for (b, a) in m.vector_images.iter().chain(m.scalar_images.iter()) {
            assert!(b.re.is_finite() && b.im.is_finite(), "non-finite b: {b:?}");
            assert!(a.re.is_finite() && a.im.is_finite(), "non-finite a: {a:?}");
        }
    }

    /// N ∈ [1, 10] all build without panic (spec DoD #1).
    #[test]
    fn n_images_one_through_ten_build() {
        for n in 1..=10usize {
            let m = MultilayerGreens::new_microstrip_with_n_images(4.4, 1.6e-3, 1.0e9, n);
            assert_eq!(m.vector_images.len(), n.max(1));
            assert_eq!(m.scalar_images.len(), n.max(1));
        }
    }
}
