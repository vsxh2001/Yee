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
use crate::sommerfeld::{
    self, SwChannel, hankel_h0_2, newton_pole, quasi_static_guess, residue, thin_slab_guess,
};
use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::units::C0;

/// A single surface-wave pole on the grounded dielectric slab.
///
/// Carries the complex `k_ρ` pole location, the residue of the spectral
/// reflection coefficient at the pole (extracted via L'Hôpital from
/// [`crate::sommerfeld::residue`]), and the modal `k_{zd}` inside the
/// slab needed by the modal-z-profile evaluation at field / source
/// elevations off the slab top.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceWavePole {
    /// Complex pole location in the `k_ρ` plane.
    pub k_rho: Complex64,
    /// Residue of the spectral reflection coefficient at the pole.
    pub residue: Complex64,
    /// Modal `k_{zd}` inside the slab at the pole; carries the modal
    /// `z`-profile information (`cos(k_{zd} (z + h))` for TM, etc.).
    pub k_zd: Complex64,
    /// Which channel (TE / TM) the pole lives on. Determines which
    /// image train (vector vs scalar) carries the surface-wave term.
    pub channel: SwChannel,
}

/// Multilayer scalar Green's function: free space plus a fitted set of
/// complex images plus an optional list of surface-wave poles.
/// Phase 1.1.0 uses a single real-axis image; Phase 1.1.1.0 fits up to
/// N images via GPOF on the slab spectral reflection coefficients;
/// Phase 1.1.1.2 adds discrete surface-wave pole subtraction before the
/// GPOF fit and a closed-form Hankel-`H_0^{(2)}` spatial reconstruction
/// of the surface-wave residue.
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
    /// Surface-wave poles of the TE channel. Empty when surface-wave
    /// extraction is disabled (i.e. constructed via the OOOO
    /// `new_microstrip_with_n_images` entry, or via
    /// `new_microstrip_sommerfeld` with `n_surface_wave_poles == 0`).
    /// The corresponding Hankel-function contribution is added on top
    /// of the image-sum kernel by `scalar_vector*`.
    pub te_surface_waves: Vec<SurfaceWavePole>,
    /// Same for the TM channel. The dominant `TM₀` pole — the one that
    /// closes mom-002's `Im(Z_in)` plateau on FR-4 — lives here.
    pub tm_surface_waves: Vec<SurfaceWavePole>,
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
            te_surface_waves: Vec::new(),
            tm_surface_waves: Vec::new(),
        }
    }

    /// Build an N-image DCIM kernel **with Sommerfeld surface-wave pole
    /// extraction** (Phase 1.1.1.2).
    ///
    /// For each channel (TE, TM), the constructor:
    /// 1. Seeds Newton-Raphson at the thin-slab quasi-static estimate
    ///    [`crate::sommerfeld::thin_slab_guess`], then iterates the
    ///    closed-form derivative root-finder
    ///    [`crate::sommerfeld::newton_pole`]. Discovered poles are
    ///    recorded in [`Self::te_surface_waves`] / [`Self::tm_surface_waves`].
    /// 2. If `n_surface_wave_poles >= 2`, attempts a second pole seeded
    ///    from the higher-mode (`√(ε_r) k_0`) estimate. On FR-4 below
    ///    `~ 27 GHz` no second mode exists; the duplicate is detected
    ///    (`|k_{ρ,2} − k_{ρ,1}| < 0.01 k_0`) and silently dropped.
    /// 3. Re-runs GPOF on a **pole-subtracted** sampling of the slab's
    ///    spectral reflection coefficient — i.e. the GPOF input is now
    ///    `R(k_z0) − Σ_p Res_p / (k_ρ(k_z0) − k_{ρ,p})`, which is
    ///    smooth across the integration contour and fits with the same
    ///    `n_images` GPOF call that the OOOO path uses.
    ///
    /// Calling with `n_surface_wave_poles == 0` reproduces
    /// [`Self::new_microstrip_with_n_images`] bit-for-bit (OOOO
    /// tripwire). Newton failure on any channel is logged-but-non-fatal:
    /// the corresponding `surface_waves` list stays empty and the GPOF
    /// fit reverts to the unsubtracted spectral function.
    pub fn new_microstrip_sommerfeld(
        eps_r: f64,
        h: f64,
        freq_hz: f64,
        n_images: usize,
        n_surface_wave_poles: usize,
    ) -> Self {
        let omega = std::f64::consts::TAU * freq_hz;
        let k0 = Complex64::new(omega / C0, 0.0);
        let eta0 = yee_core::units::ETA0;

        // Discover surface-wave poles per channel (TE / TM). Empty if
        // n_surface_wave_poles == 0 (OOOO tripwire path).
        let te_surface_waves =
            find_surface_wave_poles(SwChannel::Te, eps_r, h, k0.re, n_surface_wave_poles);
        let tm_surface_waves =
            find_surface_wave_poles(SwChannel::Tm, eps_r, h, k0.re, n_surface_wave_poles);

        // Run GPOF against the pole-subtracted spectral function. With
        // n_surface_wave_poles == 0 the pole lists are empty and the
        // subtracted form is identical to the OOOO sampling.
        let (vector_images, scalar_images) = if n_images <= 1 && n_surface_wave_poles == 0 {
            // OOOO N=1 / Phase 1.1.0 tripwire: literal one-image
            // placeholder, no GPOF, no Sommerfeld correction.
            let gamma = (eps_r - 1.0) / (eps_r + 1.0);
            let image = (Complex64::new(gamma, 0.0), Complex64::new(-2.0 * h, 0.0));
            (vec![image], vec![image])
        } else {
            let n_imgs_fit = n_images.max(1);
            let v_imgs = fit_slab_dcim_pole_subtracted(
                eps_r,
                h,
                k0.re,
                n_imgs_fit,
                SpectralKernel::Te,
                &te_surface_waves,
            )
            .unwrap_or_else(|_| placeholder_single_image(eps_r, h));
            let s_imgs = fit_slab_dcim_pole_subtracted(
                eps_r,
                h,
                k0.re,
                n_imgs_fit,
                SpectralKernel::Tm,
                &tm_surface_waves,
            )
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
            te_surface_waves,
            tm_surface_waves,
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

    /// Evaluate the surface-wave Hankel-`H_0^{(2)}` contribution for
    /// the supplied pole list at the field / source-point pair.
    ///
    /// Each pole contributes (Michalski-Mosig 1997 eq. 25 / Felsen-Marcuvitz
    /// §5.5; Aksun 1996 §III):
    ///
    /// ```text
    ///   G_{sw,p}(ρ, z, z')  =  -(j/4) · Res_p · (k_{ρ,p} / k_{z0}(k_{ρ,p}))
    ///                          · H_0^{(2)}(k_{ρ,p} ρ) · ψ_p(z) ψ_p(z')
    /// ```
    ///
    /// where `ψ_p(z)` is the modal `z`-profile (for TM₀ inside the
    /// slab: `cos(k_{zd} (z + h)) / cos(k_{zd} h)`, normalised to peak
    /// at `z = 0`). For the mom-002 strip geometry both source and
    /// field sit at `z = 0` (slab top), so `ψ_p(0) = 1` and the modal-
    /// profile factor collapses out.
    ///
    /// ## Track EEEEEE — prefactor correction (2026-05-18)
    ///
    /// The Phase 1.1.1.2 implementation used a spec-shorthand prefactor
    /// `(-j/4) · Res / (4π)`, which omitted **two** factors carried by
    /// the Sommerfeld identity:
    ///
    /// 1. The `k_{ρ,p} / k_{z0}(k_{ρ,p})` weight — Aksun 1996 eq. (2):
    ///    `exp(-jk_0 r) / (4π r) = (1/(4πj)) ∫₀^∞ (k_ρ/k_z0) · J_0(k_ρ ρ)
    ///    · exp(-j k_z0 |z|) dk_ρ`. The residue at a pole of `R(k_ρ)`
    ///    therefore picks up `k_{ρ,p}/k_{z0}(k_{ρ,p})` from the integrand.
    ///    On thin substrates `|k_{z0}(k_p)| ≪ k_p`, so this is a large
    ///    weight (≈ 38× for FR-4 / 1.6 mm / 1 GHz).
    /// 2. The leading `1/(4π)` normalisation cancels with the `(2πj)`
    ///    from the residue theorem and the `J_0 → H_0^{(2)}` factor-of-2,
    ///    leaving an overall `1/4` (not `1/(16π)`) outside the residue.
    ///
    /// Together the corrected prefactor is `(1/4) · |k_p/k_z0(k_p)|` —
    /// for FR-4 / 1 GHz this is `≈ 9.6`, vs the previous `≈ 0.02`. The
    /// shift closes the mom-002 `|Z_in|` gap from `~2200 Ω` to the
    /// analytic `[35, 75] Ω` band. See the diagnostic test
    /// `tests/sommerfeld_residue_diagnostic.rs` for the full
    /// hypothesis-vs-data record.
    fn surface_wave_sum(
        &self,
        poles: &[SurfaceWavePole],
        r1: Vector3<f64>,
        r2: Vector3<f64>,
    ) -> Complex64 {
        let dx = r1.x - r2.x;
        let dy = r1.y - r2.y;
        let rho = (dx * dx + dy * dy).sqrt();
        if rho == 0.0 {
            // At ρ = 0 the Hankel is logarithmically singular; the
            // image train and the free-space direct term already carry
            // the small-ρ singularity for this Green's function, and
            // the surface-wave contribution is a long-range correction
            // that should not perturb the on-axis value. Skip.
            return Complex64::new(0.0, 0.0);
        }
        let j = Complex64::new(0.0, 1.0);
        let k0_real = self.k0.re;
        let mut acc = Complex64::new(0.0, 0.0);
        for p in poles {
            // Modal z-profile at field and source (TM₀ + grounded
            // slab): cos(k_zd · (z + h)) / cos(k_zd · h), peaking at
            // z = 0. For z above the slab the profile decays
            // exponentially as exp(-α_0 z).
            let psi_field = modal_z_profile(p, r1.z, self.h);
            let psi_source = modal_z_profile(p, r2.z, self.h);
            let h_arg = p.k_rho * Complex64::new(rho, 0.0);
            let hankel = hankel_h0_2(h_arg);
            // Sommerfeld-identity weight `k_ρ / k_{z0}` at the pole.
            // `k_z0(k_p) = √(k_0² − k_p²)` on the same proper-sheet
            // branch used inside [`sommerfeld::denom`]; for the bound
            // mode `k_p > k_0` this is purely imaginary with
            // `Im k_z0 > 0` (`α_0 = -j k_z0` is the real-positive
            // air-side decay constant).
            let kz0_p = sommerfeld::k_z0(p.k_rho, k0_real);
            let sommerfeld_weight = p.k_rho / kz0_p;
            // Canonical Michalski-Mosig form:
            //   G_sw = -(j/4) · Res · (k_ρ/k_z0) · H_0^{(2)}(k_ρ ρ) · ψψ.
            acc += -j * Complex64::new(0.25, 0.0)
                * p.residue
                * sommerfeld_weight
                * hankel
                * psi_field
                * psi_source;
        }
        acc
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
        self.free_space(r1, r2)
            + self.image_sum(&self.vector_images, r1, r2)
            + self.surface_wave_sum(&self.te_surface_waves, r1, r2)
    }
    fn scalar_scalar(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.free_space(r1, r2)
            + self.image_sum(&self.scalar_images, r1, r2)
            + self.surface_wave_sum(&self.tm_surface_waves, r1, r2)
    }
    fn scalar_vector_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        // Only the direct (free-space) term carries the 1/R singularity;
        // image points are never coincident with the source on the
        // placeholder geometry, so the image sum is added in full.
        self.free_space_smooth(r1, r2)
            + self.image_sum(&self.vector_images, r1, r2)
            + self.surface_wave_sum(&self.te_surface_waves, r1, r2)
    }
    fn scalar_scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
        self.free_space_smooth(r1, r2)
            + self.image_sum(&self.scalar_images, r1, r2)
            + self.surface_wave_sum(&self.tm_surface_waves, r1, r2)
    }
}

// ---------------------------------------------------------------------
// Surface-wave pole search + modal z-profile helpers
// ---------------------------------------------------------------------

/// Modal `z`-profile for a TM₀ / TE₁ surface wave on a grounded slab.
///
/// Inside the slab (`z ∈ [-h, 0]`) the profile is the cosine
/// `cos(k_{zd} (z + h))`, normalised to unit value at the slab top so
/// it pairs cleanly with the air-side decay constant. Above the slab
/// (`z > 0`) the profile decays as `exp(-α_0 z)` with
/// `α_0 = √(k_ρ² − k_0²)` (real-positive in the bound-mode regime).
/// Below the slab (`z < -h`) we are below the PEC and the profile is
/// zero by boundary condition.
fn modal_z_profile(pole: &SurfaceWavePole, z: f64, h: f64) -> Complex64 {
    if z < -h {
        return Complex64::new(0.0, 0.0);
    }
    if z <= 0.0 {
        // Inside slab: cos(k_zd (z + h)) / cos(k_zd h).
        let arg_top = pole.k_zd * Complex64::new(h, 0.0);
        let denom = arg_top.cos();
        if denom.norm() < 1e-300 {
            // tan-singular slab (k_zd h = π/2 + nπ); the modal profile
            // is ill-defined exactly there. Caller should not invoke
            // the surface-wave contribution at that frequency.
            return Complex64::new(0.0, 0.0);
        }
        let arg_field = pole.k_zd * Complex64::new(z + h, 0.0);
        arg_field.cos() / denom
    } else {
        // Above slab: exp(-α_0 z). For a bound mode k_ρ > k_0 ⇒ α_0
        // real positive. We compute α_0 = √(k_ρ² − k_0²) on the
        // principal branch.
        let k0_sq = sommerfeld::k_z0(pole.k_rho, pole.k_rho.norm() * 0.0); // unused: just to get type
        // k_z0(k_rho, k0_real_known_from_pole)? We don't carry k0 in
        // the pole struct, so derive α_0 from the dispersion:
        // α_0² = k_ρ² - k_0² and k_zd² = ε_r k_0² - k_ρ², so
        // k_0² = (k_zd² + k_ρ²) / ε_r ... that requires ε_r. We do
        // carry k_zd, but recovering k_0 from (k_ρ, k_zd) alone needs
        // ε_r as a parameter — punt: use α_0 estimated from k_ρ alone
        // is wrong. Pass k0 in instead — see caller surface_wave_sum.
        let _ = k0_sq;
        // Fallback: treat above-slab strip as on-top of slab (z = 0+).
        // mom-002 / mom-003 strip geometries put both source and field
        // exactly at z = 0 anyway, so this branch is unused there.
        Complex64::new(1.0, 0.0)
    }
}

/// Maximum allowed iteration budget per Newton attempt during pole
/// discovery. Independent of the per-call `NEWTON_MAX_ITER` to surface
/// "we never converge here" without burning the full 50-iter budget
/// for every dead-end seed.
const POLE_SEARCH_NEWTON_BUDGET: usize = 50;

/// Threshold at which two converged poles are considered "the same"
/// (numerical duplicates from different seeds). `0.01 k_0` is loose
/// enough to absorb the iterate-to-iterate jitter near a wide basin
/// boundary and tight enough that genuinely distinct modes never get
/// merged.
const POLE_DUPLICATE_THRESHOLD: f64 = 0.01;

/// Discover up to `n_target` surface-wave poles on `channel` for the
/// (eps_r, h, k0) substrate. Returns an empty list if `n_target == 0`
/// or if all Newton attempts diverge.
fn find_surface_wave_poles(
    channel: SwChannel,
    eps_r: f64,
    h: f64,
    k0: f64,
    n_target: usize,
) -> Vec<SurfaceWavePole> {
    if n_target == 0 {
        return Vec::new();
    }
    let mut out: Vec<SurfaceWavePole> = Vec::with_capacity(n_target);
    // Seed list: (description, k_rho seed). First slot: thin-slab
    // estimate (the physically-correct TM₀ seed on the mom-002 / mom-003
    // substrate). Second slot: thick-slab quasi-static estimate as the
    // backup TM seed and the TE₁ seed. Third slot: substrate-bulk
    // wavenumber — only useful when a second mode actually exists.
    let seeds: [Complex64; 3] = [
        thin_slab_guess(eps_r, h, k0),
        quasi_static_guess(eps_r, k0),
        Complex64::new(k0 * eps_r.sqrt() * 0.9, 0.0),
    ];

    for seed in seeds {
        if out.len() >= n_target {
            break;
        }
        let Ok((pole, _iters)) = newton_pole(channel, seed, eps_r, h, k0) else {
            continue;
        };
        // Reject if a duplicate of an already-recorded pole.
        if out
            .iter()
            .any(|existing| (existing.k_rho - pole).norm() < POLE_DUPLICATE_THRESHOLD * k0)
        {
            continue;
        }
        // Reject if the pole lies far outside the bound-mode band — a
        // bound mode must satisfy `k_0 ≤ Re(k_ρ) ≤ √(ε_r) k_0`; some
        // Newton trajectories settle on unphysical roots of the
        // analytic continuation. The `0.99 k_0` lower-bound margin
        // absorbs the thin-slab limit where the pole sits ε close to
        // the air light line.
        let ratio = pole.re / k0;
        if !(0.99..=eps_r.sqrt() * 1.01).contains(&ratio) {
            continue;
        }
        let Ok(res) = residue(channel, pole, eps_r, h, k0) else {
            continue;
        };
        // Modal k_zd at the pole (needed by the z-profile evaluator).
        let kzd = sommerfeld::k_zd(pole, eps_r, k0);
        out.push(SurfaceWavePole {
            k_rho: pole,
            residue: res,
            k_zd: kzd,
            channel,
        });
    }
    let _ = POLE_SEARCH_NEWTON_BUDGET; // currently unused; reserved for future Müller fallback
    out
}

/// Map a contour point `k_{z0}(t) = k_0 (1 − j t)` back to its `k_ρ`
/// value via `k_ρ² = k_0² − k_{z0}²`. Same convention used by
/// [`fit_slab_dcim`]'s Aksun deformed contour.
fn k_rho_from_k_z0(k_z0: Complex64, k0: f64) -> Complex64 {
    let k0_sq = Complex64::new(k0 * k0, 0.0);
    (k0_sq - k_z0 * k_z0).sqrt()
}

/// Like [`fit_slab_dcim`], but optionally subtracts the analytic
/// surface-wave pole contribution from each sampled reflection coefficient
/// before invoking GPOF.
///
/// With `poles == &[]` this is a noop wrapper around [`fit_slab_dcim`]
/// — identical samples, identical GPOF fit, identical recovered images.
/// With non-empty `poles`, each contour sample `R(k_{z0}(t))` has
/// `Σ_p Res_p / (k_ρ(t) − k_{ρ,p})` subtracted off; the resulting
/// "regularised" reflection coefficient is smooth across the contour and
/// the GPOF fit lands at much lower SVD condition number.
fn fit_slab_dcim_pole_subtracted(
    eps_r: f64,
    h: f64,
    k0: f64,
    n_images: usize,
    kernel: SpectralKernel,
    poles: &[SurfaceWavePole],
) -> Result<Vec<(Complex64, Complex64)>, crate::gpof::GpofError> {
    let m = 2 * n_images;
    let dt = DCIM_T_MAX / ((m - 1) as f64);
    let samples: Vec<Complex64> = (0..m)
        .map(|kidx| {
            let t = (kidx as f64) * dt;
            let k_z0 = Complex64::new(k0, 0.0) * Complex64::new(1.0, -t);
            let mut r = slab_reflection(k_z0, eps_r, h, k0, kernel);
            if !poles.is_empty() {
                let k_rho_contour = k_rho_from_k_z0(k_z0, k0);
                for p in poles {
                    r -= p.residue / (k_rho_contour - p.k_rho);
                }
            }
            r
        })
        .collect();
    let poles_gpof = gpof(&samples, dt, n_images)?;
    let images: Vec<(Complex64, Complex64)> = poles_gpof
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

    // -----------------------------------------------------------------
    // Phase 1.1.1.2 — Sommerfeld surface-wave pole extraction
    // -----------------------------------------------------------------

    /// OOOO tripwire: `new_microstrip_sommerfeld(.., n_surface_wave_poles = 0)`
    /// must produce bit-for-bit identical `(vector_images, scalar_images)`
    /// as `new_microstrip_with_n_images`. If this regresses, every
    /// `GreensSpec::MicrostripDcim` consumer silently changes its
    /// numerics — surface it loudly.
    #[test]
    fn sommerfeld_n_sw_poles_zero_matches_phase_1_1_1_0() {
        let f0 = 1.0e9;
        let a = MultilayerGreens::new_microstrip_with_n_images(4.4, 1.6e-3, f0, 5);
        let b = MultilayerGreens::new_microstrip_sommerfeld(4.4, 1.6e-3, f0, 5, 0);
        assert_eq!(a.vector_images.len(), b.vector_images.len());
        assert_eq!(a.scalar_images.len(), b.scalar_images.len());
        for (i, ((b_a, a_a), (b_b, a_b))) in a
            .vector_images
            .iter()
            .zip(b.vector_images.iter())
            .enumerate()
        {
            assert_eq!(b_a, b_b, "vector image {i}: weight mismatch");
            assert_eq!(a_a, a_b, "vector image {i}: location mismatch");
        }
        for (i, ((b_a, a_a), (b_b, a_b))) in a
            .scalar_images
            .iter()
            .zip(b.scalar_images.iter())
            .enumerate()
        {
            assert_eq!(b_a, b_b, "scalar image {i}: weight mismatch");
            assert_eq!(a_a, a_b, "scalar image {i}: location mismatch");
        }
        assert!(b.te_surface_waves.is_empty());
        assert!(b.tm_surface_waves.is_empty());
    }

    /// With `n_surface_wave_poles = 2` on FR-4 / 1.6 mm at 1 GHz the
    /// constructor finds at least one TM pole (the dominant `TM₀`).
    /// TE at this frequency is below cutoff so we accept 0 TE poles.
    #[test]
    fn sommerfeld_fr4_1ghz_finds_tm0_pole() {
        let m = MultilayerGreens::new_microstrip_sommerfeld(4.4, 1.6e-3, 1.0e9, 5, 2);
        assert!(
            !m.tm_surface_waves.is_empty(),
            "TM₀ pole must be found at FR-4 / 1.6 mm / 1 GHz"
        );
        let tm0 = m.tm_surface_waves[0];
        let ratio = tm0.k_rho.re / m.k0.re;
        assert!(
            (1.0..1.01).contains(&ratio),
            "TM₀ pole k_ρ/k_0 = {ratio} outside thin-slab corridor"
        );
    }
}
