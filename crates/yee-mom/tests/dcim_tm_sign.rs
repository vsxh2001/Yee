//! DCIM TM-channel sign gate — Track DDDDDDD (M3 root-cause fix).
//!
//! ## Why this file exists
//!
//! Track YYYYYY (commit `9b140d4`) measured the leading DCIM image of the
//! `MultilayerGreens::scalar_images` train on FR-4 / `h = 1.6 mm` / 1 GHz
//! as
//!
//! ```text
//!   b_leading_TM ≈ -0.92 + j 0.11    at    a_leading_TM ≈ -1.08 mm.
//! ```
//!
//! The analytic Chow-1991 / Aksun-1996 leading image for a horizontal
//! current source above a grounded dielectric slab is
//!
//! ```text
//!   b ≈ (ε_r − 1) / (ε_r + 1) ≈ 0.629    at    a ≈ -2 h = -3.2 mm.
//! ```
//!
//! — **positive real**. The TE channel (`MultilayerGreens::vector_images`)
//! had the correct sign already: its leading image is the PEC ground at
//! `b = -1` at `a = -2 h = -3.2 mm`.
//!
//! Track SSSSSS (`mom_002_reflection_convention.rs`) traced the underlying
//! cause to a `−2×` discrepancy in `sommerfeld::residue()`; Track TTTTTT
//! (commit `43c01aa`) applied the Michalski-Mosig 1997 eq. 19 sign + factor-of-2
//! correction to `residue()`, which fixed the residue function in
//! isolation. The residue check passes (ratio `≈ 1.00`).
//!
//! But the **sampling form** that GPOF actually fits — `slab_reflection`'s
//! TM branch — was still Pozar's TM Fresnel form
//!
//! ```text
//!   R_TM_Pozar = (j · k_zd · tan(k_zd h) − ε_r · k_z0)
//!              / (j · k_zd · tan(k_zd h) + ε_r · k_z0),
//! ```
//!
//! which is the **negative** of the Michalski-Mosig 1991 / 1997 scalar-
//! potential spectral kernel
//!
//! ```text
//!   F_q = -R_TM_Pozar
//!       = (ε_r · k_z0 − j · k_zd · tan(k_zd h))
//!       / (ε_r · k_z0 + j · k_zd · tan(k_zd h)).
//! ```
//!
//! The DCIM expansion of the **scalar-potential** Green's function fits
//! `F_q`, not `R_TM_Pozar`. With the Pozar form fed to GPOF, the leading
//! image carried the wrong sign and ended up tracking the **PEC short**
//! at `R(k_z0 = k_0) → -1` instead of the dielectric-correction image at
//! `b = +(ε_r − 1) / (ε_r + 1)`.
//!
//! Track DDDDDDD's fix flips the TM branch of `slab_reflection` to the
//! `F_q` convention. The TM `SurfaceWavePole.residue` is correspondingly
//! negated at the storage site (in `find_surface_wave_poles`) so the
//! surface-wave Hankel reconstruction in `surface_wave_sum` continues to
//! see the residue of the *new* spectral kernel.
//!
//! ## What this gate checks
//!
//! 1. The leading TM-channel DCIM image on FR-4 has positive real
//!    coefficient (sign of `(ε_r − 1) / (ε_r + 1)`).
//! 2. The magnitude is within `[0.1, 5] ×` the analytic `0.629` target.
//! 3. The depth is within `[0.2, 5] ×` the analytic `-2 h = -3.2 mm`
//!    target.
//!
//! Identical sanity criteria to Track YYYYYY's M3 probe, but expressed
//! as an assert-on-fail gate test rather than an `#[ignore]`d
//! diagnostic.
//!
//! ## References
//!
//! * Track YYYYYY — `crates/yee-mom/tests/mom_002_mpie_audit.rs` (M3 probe).
//! * Track TTTTTT — `crates/yee-mom/src/sommerfeld.rs::residue` (Michalski-
//!   Mosig 1997 eq. 19 residue formula).
//! * Y. L. Chow et al., "A closed-form spatial Green's function for the
//!   thick microstrip substrate," *IEEE Trans. Microw. Theory Tech.*,
//!   vol. 39, no. 3, pp. 588–592, Mar 1991.
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans.
//!   Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.

use yee_mom::__internal::MultilayerGreens;

const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const F_HZ: f64 = 1.0e9;
const N_DCIM_IMAGES: usize = 5;

/// Gate: the leading TM-channel DCIM image is positive-real on FR-4 /
/// 1 GHz / 1.6 mm at depth on the order of `-2 h`.
#[test]
fn dcim_tm_leading_image_is_positive_real_on_fr4() {
    let greens =
        MultilayerGreens::new_microstrip_with_n_images(EPS_R, H_SUBSTRATE_M, F_HZ, N_DCIM_IMAGES);

    // Pick the scalar-image entry with the largest |b| — that is the
    // "leading image" in the Chow / Aksun decomposition.
    let leading = greens
        .scalar_images
        .iter()
        .copied()
        .max_by(|(b1, _), (b2, _)| {
            b1.norm()
                .partial_cmp(&b2.norm())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("scalar_images must be non-empty");

    let (b_lead, a_lead) = leading;
    eprintln!(
        "dcim_tm_sign: leading TM image  b = {:.6e} + j·{:.6e}  a = {:.6e} + j·{:.6e}",
        b_lead.re, b_lead.im, a_lead.re, a_lead.im,
    );

    let b_expected = (EPS_R - 1.0) / (EPS_R + 1.0);
    let a_expected = -2.0 * H_SUBSTRATE_M;

    // Sign: the dielectric-correction image carries `(ε_r − 1) /
    // (ε_r + 1) > 0`. The pre-fix code returned `b ≈ -0.92` here.
    assert!(
        b_lead.re > 0.0,
        "TM leading-image Re(b) must be positive (sign of (ε_r-1)/(ε_r+1)), got {:.4e}",
        b_lead.re,
    );

    // Magnitude — keep YYYYYY's [0.1×, 5×] sanity band so a future
    // formulation tweak (e.g. dropping the PEC image into the DCIM
    // fit) does not falsely regress.
    let mag_ratio = b_lead.re.abs() / b_expected.abs();
    assert!(
        (0.1..=5.0).contains(&mag_ratio),
        "TM leading-image |Re(b)| / 0.629 = {:.3e} outside [0.1, 5.0] band",
        mag_ratio,
    );

    // Depth — within YYYYYY's [0.2×, 5×] sanity band of the analytic
    // `-2 h`.
    let depth_ratio = a_lead.re.abs() / a_expected.abs();
    assert!(
        (0.2..=5.0).contains(&depth_ratio),
        "TM leading-image |Re(a)| / 2 h = {:.3e} outside [0.2, 5.0] band",
        depth_ratio,
    );
}

/// Cross-check: the TE-channel leading image continues to track the
/// PEC ground image at `b = -1` / `a = -2 h`. If this regresses, the
/// TM sign-flip accidentally bled into the TE branch.
#[test]
fn dcim_te_leading_image_still_at_pec_minus_one() {
    let greens =
        MultilayerGreens::new_microstrip_with_n_images(EPS_R, H_SUBSTRATE_M, F_HZ, N_DCIM_IMAGES);

    let leading = greens
        .vector_images
        .iter()
        .copied()
        .max_by(|(b1, _), (b2, _)| {
            b1.norm()
                .partial_cmp(&b2.norm())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("vector_images must be non-empty");

    let (b_lead, a_lead) = leading;
    eprintln!(
        "dcim_te_check: leading TE image  b = {:.6e} + j·{:.6e}  a = {:.6e} + j·{:.6e}",
        b_lead.re, b_lead.im, a_lead.re, a_lead.im,
    );

    let a_expected = -2.0 * H_SUBSTRATE_M;

    // TE PEC ground image: b = -1, a = -2 h. Loose magnitude band
    // (0.5–1.5) plus correct sign.
    assert!(
        b_lead.re < 0.0,
        "TE leading-image Re(b) must be negative (PEC reflection coefficient -1), got {:.4e}",
        b_lead.re,
    );
    assert!(
        (0.5..=1.5).contains(&b_lead.re.abs()),
        "TE leading-image |Re(b)| = {:.4e} outside [0.5, 1.5] PEC band",
        b_lead.re.abs(),
    );

    let depth_ratio = a_lead.re.abs() / a_expected.abs();
    assert!(
        (0.5..=1.5).contains(&depth_ratio),
        "TE leading-image |Re(a)| / 2 h = {:.3e} outside [0.5, 1.5] band",
        depth_ratio,
    );
}
