//! mom-002 reflection convention sanity diagnostic — Track SSSSSS.
//!
//! ## Why this file exists
//!
//! Three prior diagnostics have ruled out the leading candidates for the
//! ~30× `|Z_in|` gap that mom-002 still posts on FR-4 / `h = 1.6 mm` /
//! `w = 2.94 mm` / 1 GHz:
//!
//! * **Track EEEEEE** corrected the surface-wave reconstruction prefactor
//!   in [`yee_mom::multilayer::MultilayerGreens::surface_wave_sum`] to the
//!   canonical Michalski-Mosig 1997 form
//!   `G_sw = -(j/4) · Res · (k_ρ/k_z0) · H_0^{(2)}(k_ρ ρ) · ψ²` — closes
//!   the bulk of the magnitude gap but leaves a ~30× residual.
//!
//! * **Track JJJJJJ** swept the strip length `L ∈ [15, 80] mm` and found
//!   only a ~5% monotonic trend — not enough to close 30×.
//!
//! * **Track PPPPPP** ran a GPOF-residual / kernel-decomposition check
//!   and found the pole-subtraction is a no-op on the Aksun contour
//!   (`‖R̃‖/‖R‖ ≈ 1`) but the four-term hand-sum agrees bit-exact with
//!   `scalar_scalar`. The residue itself is the suspect: PPPPPP recorded
//!   `|Res| ≈ 2.82e-2` versus `|R| ~ 1` on the contour, suggesting the
//!   `residue()` output may be the wrong quantity by a sign, a sheet, or
//!   a `2π` factor.
//!
//! ## Hypotheses tested here
//!
//! * **(S1) Sign / sheet bug.** `residue()` and [`slab_reflection`] may
//!   carry inconsistent sign conventions (one uses `α_0 = -j·k_z0`, the
//!   other expresses the same reflection coefficient with the opposite
//!   second-term sign). If true, the L'Hôpital ratio `N(k_p)/D'(k_p)`
//!   measures the residue of the *wrong* meromorphic function.
//!
//! * **(S2) Improper Riemann sheet.** The branch
//!   `k_{z0} = √(k_0² − k_ρ²)` has two sheets. The proper sheet for a
//!   bound surface wave puts `Im(k_z0) > 0` in the `k_ρ > k_0` regime —
//!   i.e. `exp(j·k_z0·z)` decays into the air half-space. If the Newton
//!   solver lands on the improper sheet (`Im(k_z0) < 0` → radiation /
//!   leaky-wave), the residue sign flips relative to the convention used
//!   by [`slab_reflection`].
//!
//! * **(S3) Missing `2π` factor.** The Sommerfeld identity carries a
//!   `1/(2π)` from the `J_0` Fourier transform; the residue theorem
//!   carries a `2π` from `∮ dk_ρ = 2π·j·Res`. If the cancellation in the
//!   spec is off by exactly `2π`, `|Res|` is `~6.28×` too small —
//!   approaching but not closing the observed ~30× gap.
//!
//! ## Diagnostic method
//!
//! The cleanest, model-free way to check `residue()`'s output is to
//! compute the residue **numerically** via a contour integral around the
//! pole. For a simple pole of `R(k_ρ)` at `k_p`, the residue theorem
//! gives
//!
//! ```text
//!   Res = (1 / (2π·i)) · ∮ R(k_ρ) dk_ρ.
//! ```
//!
//! Parameterising the contour as a circle of radius `δ` centred at `k_p`
//! (`k_ρ(θ) = k_p + δ·e^{i·θ}`, `dk_ρ = i·δ·e^{i·θ} dθ`) and applying
//! the trapezoidal rule with `N` equally-spaced samples yields the
//! discrete approximation used below:
//!
//! ```text
//!   Res ≈ (δ / N) · Σ_n R(k_p + δ·e^{i·2π·n/N}) · e^{i·2π·n/N}.
//! ```
//!
//! `δ = 1e-4 · k_0` is small enough that higher-order Laurent terms are
//! negligible (their contribution scales as `δ²`) but large enough to
//! stay above the double-precision floor in the reflection-coefficient
//! evaluation. `N = 64` over-resolves the trapezoidal rule — at this
//! sample count the discrete approximation is exact to machine epsilon
//! for any analytic function on the disc, including the pole part.
//!
//! Comparing `Res` (from [`yee_mom::__internal::sommerfeld::residue`])
//! with `Res_numerical` (from the contour integral) directly tests S1
//! and S3: any sign / sheet / `2π` mismatch shows up as a non-unity
//! ratio in either phase or magnitude. The Riemann-sheet check on
//! `k_z0(k_p)` directly tests S2.
//!
//! ## References
//!
//! * Track EEEEEE prefactor-correction record:
//!   `crates/yee-mom/tests/sommerfeld_residue_diagnostic.rs`.
//! * Track JJJJJJ Hankel-tail / extent sweep:
//!   `crates/yee-mom/tests/mom_002_extent_sensitivity.rs`.
//! * Track PPPPPP GPOF residual / kernel decomposition:
//!   `crates/yee-mom/tests/mom_002_h2_gpof_diagnostic.rs`.
//! * M. I. Aksun, "A robust approach for the derivation of closed-form
//!   Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44,
//!   no. 5, pp. 651–658, May 1996.
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans.
//!   Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.

use num_complex::Complex64;
use yee_mom::__internal::sommerfeld::{
    SwChannel, d_tm, k_z0, k_zd, newton_pole, residue, thin_slab_guess,
};

const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const F_HZ: f64 = 1.0e9;

/// Number of points on the small circle around the pole. 64 trapezoidal
/// samples are ~machine-epsilon for any analytic integrand on a disc;
/// the residue is the `e^{i·θ}`-mode of `R(k_p + δ·e^{i·θ})` and the
/// trapezoidal rule is spectrally accurate for periodic integrands.
const N_CONTOUR_SAMPLES: usize = 64;

/// Radius of the contour around the pole, in units of `k_0`. Small
/// enough that the `O(δ²)` Laurent terms are negligible; large enough
/// that `R` evaluation stays well above the double-precision floor.
const CONTOUR_RADIUS_OVER_K0: f64 = 1.0e-4;

fn k0_at(freq_hz: f64) -> f64 {
    std::f64::consts::TAU * freq_hz / yee_core::units::C0
}

/// Slab reflection coefficient on the TM channel, accepting `k_z0`
/// directly (the same form [`yee_mom::multilayer::slab_reflection`] uses
/// internally). Duplicated here because `slab_reflection` is
/// `pub(crate)` and not exposed through `__internal`.
///
/// ```text
///   R(k_z0) = (j · k_zd · tan(k_zd · h) − ε_r · k_z0)
///           / (j · k_zd · tan(k_zd · h) + ε_r · k_z0)
/// ```
///
/// with `k_zd² = (ε_r − 1) · k_0² + k_z0²` (principal branch).
fn tm_reflection_from_k_z0(k_z0_pt: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let k0_sq = Complex64::new(k0 * k0, 0.0);
    let inside = Complex64::new(eps_r - 1.0, 0.0) * k0_sq + k_z0_pt * k_z0_pt;
    let k_zd_pt = inside.sqrt();
    let phase = k_zd_pt * Complex64::new(h, 0.0);
    let t = phase.tan();
    let j = Complex64::new(0.0, 1.0);
    let num = j * k_zd_pt * t - Complex64::new(eps_r, 0.0) * k_z0_pt;
    let den = j * k_zd_pt * t + Complex64::new(eps_r, 0.0) * k_z0_pt;
    num / den
}

/// Slab reflection coefficient evaluated as a function of `k_ρ`. Picks
/// the principal-branch `k_z0 = √(k_0² − k_ρ²)`, the same convention
/// the [`yee_mom::__internal::sommerfeld::k_z0`] helper uses. This is
/// the meromorphic function whose residue we want to extract; it has a
/// simple pole at the TM₀ surface-wave wavenumber `k_p`.
fn tm_reflection_from_k_rho(k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let k0_sq = Complex64::new(k0 * k0, 0.0);
    let k_z0_pt = (k0_sq - k_rho * k_rho).sqrt();
    tm_reflection_from_k_z0(k_z0_pt, eps_r, h, k0)
}

/// Numerical residue of `f` at `pole` via the trapezoidal rule on a
/// circle of radius `delta_k`:
///
/// ```text
///   Res ≈ (δ / N) · Σ_n f(pole + δ·e^{i·2π·n/N}) · e^{i·2π·n/N}.
/// ```
///
/// Derived from `Res = (1/(2π·i)) · ∮ f dk_ρ` with the parameterisation
/// `k_ρ = pole + δ·e^{i·θ}`, `dk_ρ = i·δ·e^{i·θ} dθ`. The factors of
/// `i` cancel, leaving the real-valued `δ/N` prefactor and the
/// `e^{i·θ}` weight on each sample.
fn numerical_residue<F>(f: F, pole: Complex64, delta_k: f64, n_samples: usize) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let mut acc = Complex64::new(0.0, 0.0);
    let two_pi_over_n = std::f64::consts::TAU / (n_samples as f64);
    for n in 0..n_samples {
        let theta = (n as f64) * two_pi_over_n;
        let e_i_theta = Complex64::from_polar(1.0, theta);
        let z = pole + Complex64::new(delta_k, 0.0) * e_i_theta;
        acc += f(z) * e_i_theta;
    }
    acc * Complex64::new(delta_k / (n_samples as f64), 0.0)
}

/// Print the reflection convention diagnostic table and the verdict.
///
/// Marked `#[ignore]` so the suite never runs it by default; invoke
/// explicitly via
///
/// ```text
/// cargo test -p yee-mom --release --test mom_002_reflection_convention \
///     -- --ignored --nocapture
/// ```
///
/// to dump the table and the S1/S2/S3 verdict.
#[test]
#[ignore = "diagnostic: prints residue convention / Riemann-sheet / 2π-factor numerics for mom-002"]
fn mom_002_reflection_convention_diagnostic() {
    let k0 = k0_at(F_HZ);
    let seed = thin_slab_guess(EPS_R, H_SUBSTRATE_M, k0);
    let (pole, iters) =
        newton_pole(SwChannel::Tm, seed, EPS_R, H_SUBSTRATE_M, k0).expect("Newton converges");
    let resid_d = d_tm(pole, EPS_R, H_SUBSTRATE_M, k0).norm();
    let res = residue(SwChannel::Tm, pole, EPS_R, H_SUBSTRATE_M, k0).expect("residue finite");

    // Sample R(k_ρ) magnitude away from the pole for an "ambient" scale.
    // Used in the verdict to put |2π · Res| in context.
    let r_ambient =
        tm_reflection_from_k_rho(Complex64::new(0.5 * k0, 0.0), EPS_R, H_SUBSTRATE_M, k0);

    // Numerical residue via the trapezoidal rule on a circle of radius
    // δ = 1e-4 · k_0 around the converged pole.
    let delta_k = CONTOUR_RADIUS_OVER_K0 * k0;
    let res_num = numerical_residue(
        |z| tm_reflection_from_k_rho(z, EPS_R, H_SUBSTRATE_M, k0),
        pole,
        delta_k,
        N_CONTOUR_SAMPLES,
    );

    // Riemann-sheet check on k_z0 at the pole.
    let kz0_at_pole = k_z0(pole, k0);
    let kzd_at_pole = k_zd(pole, EPS_R, k0);
    let on_proper_sheet = kz0_at_pole.im > 0.0;

    // Ratio Res / Res_num — both magnitude and phase reveal sign / sheet
    // / 2π mismatches between the two extractions.
    let ratio = if res_num.norm() > 0.0 {
        res / res_num
    } else {
        Complex64::new(f64::NAN, f64::NAN)
    };

    eprintln!("--- mom-002 reflection convention diagnostic ---");
    eprintln!(
        "FR-4 / ε_r = {EPS_R}, h = {} mm, f = {} GHz, TM channel",
        H_SUBSTRATE_M * 1e3,
        F_HZ * 1e-9,
    );
    eprintln!("k_0      = {:.6e} rad/m", k0);
    eprintln!("Newton   = {iters} iters, |D|={:.2e}", resid_d);
    eprintln!(
        "k_p      = {:.6e} + j·{:.6e} rad/m  (k_p/k_0 = {:.6} + j·{:.6})",
        pole.re,
        pole.im,
        pole.re / k0,
        pole.im / k0,
    );
    eprintln!(
        "contour  = circle of radius δ = {:.2e} · k_0 = {:.4e} rad/m, N = {} samples",
        CONTOUR_RADIUS_OVER_K0, delta_k, N_CONTOUR_SAMPLES,
    );

    eprintln!();
    eprintln!(
        "Res from residue():     re = {:>14.6e}, im = {:>14.6e},  |Res|     = {:.6e}",
        res.re,
        res.im,
        res.norm()
    );
    eprintln!(
        "Res from contour:       re = {:>14.6e}, im = {:>14.6e},  |Res_num| = {:.6e}",
        res_num.re,
        res_num.im,
        res_num.norm()
    );
    eprintln!(
        "   ratio Res / Res_num: re = {:>14.6e}, im = {:>14.6e},  |.|       = {:.6e}   arg = {:.6e} rad",
        ratio.re,
        ratio.im,
        ratio.norm(),
        ratio.arg(),
    );

    eprintln!();
    eprintln!(
        "k_z0(k_p):              re = {:>14.6e}, im = {:>14.6e}  → on PROPER sheet (Im>0)? {}",
        kz0_at_pole.re,
        kz0_at_pole.im,
        if on_proper_sheet { "yes" } else { "no" },
    );
    eprintln!(
        "k_zd(k_p):              re = {:>14.6e}, im = {:>14.6e}",
        kzd_at_pole.re, kzd_at_pole.im,
    );

    eprintln!();
    let two_pi = std::f64::consts::TAU;
    eprintln!(
        "|2π · Res|              = {:.6e}    (compare to |R| sample magnitude away from pole ~ {:.4e})",
        two_pi * res.norm(),
        r_ambient.norm(),
    );
    eprintln!("|2π · Res_num|          = {:.6e}", two_pi * res_num.norm());

    // -------------------------------------------------------------------
    // Verdict.
    // -------------------------------------------------------------------
    // Tolerances:
    //   * |ratio|  ≈ 1.0 ± 1% → magnitudes agree.
    //   * arg(ratio) ≈ 0      → phases agree (sign / sheet matches).
    //   * |ratio| ≈ 2π or 1/(2π) → S3 (missing 2π factor).
    //   * arg(ratio) ≈ ±π     → S1 (sign flip).
    //   * arg(ratio) ≈ ±π/2   → S1 with an imaginary-unit mismatch (j vs −j).
    //
    // S2 is checked independently of the ratio: it is true iff
    // `Im(k_z0(k_p)) < 0`.
    let ratio_mag = ratio.norm();
    let ratio_arg = ratio.arg();
    let mag_agrees = (ratio_mag - 1.0).abs() < 0.05;
    let phase_agrees = ratio_arg.abs() < 0.05;

    let s1_detected = !mag_agrees || !phase_agrees;
    let s2_detected = !on_proper_sheet;
    let s3_detected = (ratio_mag - two_pi).abs() / two_pi < 0.05
        || (ratio_mag - 1.0 / two_pi).abs() * two_pi < 0.05;

    eprintln!();
    eprintln!("Verdict:");
    eprintln!(
        "  S1 (sign/sheet bug):     {}",
        if s1_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!(
        "      reason: |ratio| = {:.4} (want ≈ 1.0 ± 5%), arg(ratio) = {:.4} rad (want ≈ 0 ± 0.05 rad)",
        ratio_mag, ratio_arg,
    );
    eprintln!(
        "  S2 (improper sheet):     {}",
        if s2_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!(
        "      reason: Im(k_z0(k_p)) = {:.4e} (want > 0 for proper sheet)",
        kz0_at_pole.im,
    );
    eprintln!(
        "  S3 (missing 2π factor):  {}",
        if s3_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!(
        "      reason: |ratio| = {:.4} (want ≈ 2π ≈ 6.283 or ≈ 1/(2π) ≈ 0.159 for S3)",
        ratio_mag,
    );
}
