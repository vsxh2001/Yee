//! mom-002 Hypothesis 2 — GPOF residual fit on the pole-subtracted Aksun
//! contour. Track PPPPPP root-cause diagnostic.
//!
//! ## Why this file exists
//!
//! Two prior diagnostics have eliminated the leading candidates for the
//! ~30× `|Z_in|` gap that mom-002 still posts (`|Z_in| ≈ 2232 Ω` vs
//! analytic `Z_0 ≈ 51 Ω` on FR-4 / `h = 1.6 mm` / `w = 2.94 mm` / 1 GHz):
//!
//! * **Track EEEEEE** corrected the surface-wave reconstruction prefactor
//!   in [`crate::multilayer::MultilayerGreens::surface_wave_sum`] to the
//!   canonical Michalski-Mosig 1997 form
//!   `G_sw = -(j/4) · Res · (k_ρ/k_z0) · H_0^{(2)}(k_ρ ρ) · ψ²` — but the
//!   residue contribution still moves `Z_in` by far less than the spec
//!   predicted.
//!
//! * **Track JJJJJJ** swept the strip length `L ∈ [15, 80] mm` to test
//!   whether the slow `ρ^{-1/2}` Hankel tail is being truncated by the
//!   mesh. The trend is monotonic but only ≈ 5% over that range, far too
//!   small to close the remaining factor-of-30.
//!
//! The remaining hypothesis the EEEEEE module docstring flagged is
//! **Hypothesis 2 — GPOF residual fit**: if the surface-wave pole
//! contribution `R_pole(k_ρ) = Res_p / (k_ρ − k_p)` is too small to
//! materially change the sampled `R(k_ρ)` on the Aksun contour, then the
//! GPOF re-fit of the residual `R̃ = R − R_pole` is essentially the same
//! fit the OOOO (no-pole-subtraction) path runs — and the DCIM image
//! train silently fails to absorb the surface-wave's spatial contribution.
//!
//! Independently, this file also records a falsifiability check on the
//! kernel itself: if the full Sommerfeld kernel
//! (`new_microstrip_sommerfeld(.., n_images = 5, n_surface_wave_poles = 1)`)
//! disagrees with a hand-summed four-term decomposition
//! (free-space + DCIM image train + surface-wave Hankel), then there is
//! a deeper sign / convention / sheet bug between
//! [`MultilayerGreens::scalar_scalar`] and the underlying
//! [`crate::sommerfeld`] residue / Hankel routines that no amount of
//! GPOF re-fitting will fix.
//!
//! ## What this diagnostic prints
//!
//! 1. The TM₀ pole `k_p`, `k_z0(k_p)`, `k_zd(k_p)`, and the residue
//!    `Res_p` at FR-4 / `h = 1.6 mm` / 1 GHz.
//!
//! 2. A 32-row table of the spectral reflection coefficient `R(k_ρ)`
//!    sampled along the Aksun deformed contour
//!    `k_z0(t) = k_0 (1 − j t)`, `t ∈ [0, 10]`, together with the
//!    analytic pole contribution `R_pole(k_ρ) = Res_p / (k_ρ − k_p)`
//!    and the residual `R̃ = R − R_pole`. The norms
//!    `‖R‖₂`, `‖R_pole‖₂`, `‖R̃‖₂` over the 32 contour samples
//!    quantify whether pole subtraction is materially shifting the
//!    GPOF target.
//!
//! 3. Hand-summed four-term decomposition of `G^Φ(r1, r2)` at a probe
//!    pair `r1 = (0, 0, 0)`, `r2 = (15 mm, 0, 0)` (both on the slab top):
//!    free-space + DCIM image train + surface-wave Hankel residue.
//!    Compared against the kernel's own [`Greens::scalar_scalar`]
//!    evaluation of `new_microstrip_sommerfeld(.., 5, 1)`, which sums
//!    the same three pieces internally. The two **must** agree to GPOF
//!    round-off; any larger gap is a kernel bug, not a fit quality
//!    issue.
//!
//! ## Verdict logic
//!
//! * If `‖R̃‖ / ‖R‖ ≈ 1` (pole subtraction barely changes the GPOF
//!   target) → Hypothesis 2 is **WEAK**: the pole is too small to
//!   matter, but that itself is consistent with the residue being
//!   the right size and the gap living elsewhere (e.g. modal-profile
//!   convention or the slab reflection-coefficient normalisation that
//!   feeds the residue).
//!
//! * If `‖R̃‖ / ‖R‖ ≪ 1` (pole captures most of `R` on the contour) →
//!   Hypothesis 2 is **STRONG**: GPOF fits a small smooth residual and
//!   the DCIM image train absorbs only the smooth part of the kernel,
//!   leaving the long-range Hankel decay to the surface-wave term. In
//!   that regime the kernel decomposition should land correctly, and
//!   the `|Z_in|` gap is some other physical mismatch (e.g. the strip
//!   width interacts with the lateral integration of `H_0^{(2)}` in a
//!   way the Galerkin quadrature undercounts).
//!
//! * If `|G^Φ_sommerfeld − G^Φ_hand|` is large (rel > 1e-3) → there is
//!   a **kernel inconsistency** between the surface_wave_sum and
//!   image_sum paths inside [`MultilayerGreens`]. This would invalidate
//!   the EEEEEE prefactor fix and require revisiting the
//!   `(k_p / k_z0)` weight, the `(−j/4)` prefactor, the `ψ²` modal
//!   profile at `z = 0`, or the choice of branch for `k_z0(k_p)`.
//!
//! * If `|G^Φ_sommerfeld − G^Φ_hand|` is at the round-off floor
//!   (rel < 1e-6) → the kernel is **self-consistent**; the gap is not
//!   in the residue extraction or the spatial reconstruction, and the
//!   next diagnostic must look at the Galerkin / port-current side.
//!
//! ## Escape hatch (per Track PPPPPP brief)
//!
//! The brief originally specified an explicit GPOF re-fit of the
//! residual `R̃`. The crate's `gpof::gpof` function is `pub(crate)`
//! and not exposed through `__internal`, so a literal residual re-fit
//! cannot be driven from an integration test without out-of-lane
//! changes to `src/lib.rs` (the brief restricts edits to a single
//! tests file). Per the escape hatch, the qualitative sample-norm
//! comparison `‖R‖ / ‖R̃‖` carries the same information for the
//! hypothesis: if the residual is "much smaller" the GPOF re-fit
//! lands on smooth data; if not, it does not. The kernel
//! decomposition diagnostic (step 4–5 in the brief) is the more
//! important deliverable and is fully reachable through
//! [`Greens::scalar_scalar`] alone.
//!
//! ## References
//!
//! * Track EEEEEE prefactor-correction record:
//!   `crates/yee-mom/tests/sommerfeld_residue_diagnostic.rs`.
//! * Track JJJJJJ Hankel-tail sweep:
//!   `crates/yee-mom/tests/mom_002_extent_sensitivity.rs`.
//! * M. I. Aksun, "A robust approach for the derivation of closed-form
//!   Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44,
//!   no. 5, pp. 651–658, May 1996.
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans.
//!   Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_mom::__internal::sommerfeld::{
    SwChannel, d_tm, hankel_h0_2, k_z0, k_zd, newton_pole, residue, thin_slab_guess,
};
use yee_mom::__internal::{Greens, MultilayerGreens};

const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const F_HZ: f64 = 1.0e9;
const N_CONTOUR_SAMPLES: usize = 32;
const AKSUN_T_MAX: f64 = 10.0;

fn k0_at(freq_hz: f64) -> f64 {
    std::f64::consts::TAU * freq_hz / yee_core::units::C0
}

/// Reflection coefficient on the Aksun deformed contour for the TM
/// channel. Duplicates [`yee_mom::multilayer::slab_reflection`] so the
/// diagnostic does not need to reach into private module internals.
fn tm_reflection(k_z0_pt: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
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

/// Map a contour point `k_z0(t) = k_0 (1 − j t)` back to `k_ρ` via the
/// upper-half-plane branch of `k_ρ = √(k_0² − k_z0²)`.
fn k_rho_on_contour(t: f64, k0: f64) -> Complex64 {
    let k_z0_pt = Complex64::new(k0, 0.0) * Complex64::new(1.0, -t);
    let k0_sq = Complex64::new(k0 * k0, 0.0);
    (k0_sq - k_z0_pt * k_z0_pt).sqrt()
}

/// Hand-sum the free-space + DCIM-image + surface-wave-Hankel terms for
/// the scalar-potential Green's function at `(r1, r2)`, using the public
/// fields of [`MultilayerGreens`].
///
/// This mirrors what `scalar_scalar` does internally; if both sums agree
/// the kernel is internally consistent at the four-term level.
fn hand_summed_scalar_scalar(
    g: &MultilayerGreens,
    r1: Vector3<f64>,
    r2: Vector3<f64>,
) -> (Complex64, Complex64, Complex64, Complex64) {
    // 1. Free space.
    let r = (r1 - r2).norm();
    let phase = -g.k0.re * r;
    let free_space = Complex64::from_polar(1.0 / (4.0 * std::f64::consts::PI * r), phase);

    // 2. DCIM image train (scalar / TM channel).
    let mut image_sum = Complex64::new(0.0, 0.0);
    for &(b_n, a_n) in &g.scalar_images {
        let dx = r1.x - r2.x;
        let dy = r1.y - r2.y;
        let rho_sq = dx * dx + dy * dy;
        let dz = Complex64::new(r1.z + r2.z, 0.0) - a_n;
        let r_n = (Complex64::new(rho_sq, 0.0) + dz * dz).sqrt();
        if r_n.norm() < 1e-300 {
            continue;
        }
        let ph = -g.k0 * r_n;
        let exp_term = ph.exp();
        image_sum += b_n * exp_term / (Complex64::new(4.0 * std::f64::consts::PI, 0.0) * r_n);
    }

    // 3. Surface-wave Hankel residue contribution (TM₀ pole only).
    //    G_sw = -(j/4) · Res · (k_ρ / k_z0(k_ρ)) · H_0^{(2)}(k_ρ ρ) · ψ_field · ψ_source
    //    with ψ(0) = 1 (both points sit at z = 0 on the slab top), so
    //    the modal-z-profile factor collapses out. This is the canonical
    //    Michalski-Mosig 1997 form that Track EEEEEE installed.
    let dx = r1.x - r2.x;
    let dy = r1.y - r2.y;
    let rho = (dx * dx + dy * dy).sqrt();
    let mut sw_sum = Complex64::new(0.0, 0.0);
    if rho > 0.0 {
        let j = Complex64::new(0.0, 1.0);
        for p in &g.tm_surface_waves {
            let kz0_p = k_z0(p.k_rho, g.k0.re);
            let weight = p.k_rho / kz0_p;
            let h_arg = p.k_rho * Complex64::new(rho, 0.0);
            let hk = hankel_h0_2(h_arg);
            // ψ(0) = 1 by the modal-z-profile normalisation; the slab
            // top is the strip plane.
            let psi = Complex64::new(1.0, 0.0);
            sw_sum += -j * Complex64::new(0.25, 0.0) * p.residue * weight * hk * psi * psi;
        }
    }

    let total = free_space + image_sum + sw_sum;
    (free_space, image_sum, sw_sum, total)
}

/// Diagnostic entry point. Marked `#[ignore]` so the suite never runs
/// it by default; invoke explicitly via
///
/// ```text
/// cargo test -p yee-mom --release --test mom_002_h2_gpof_diagnostic \
///     -- --ignored --nocapture
/// ```
///
/// to dump the contour-sample table, the kernel-decomposition
/// comparison, and the hypothesis-2 verdict.
#[test]
#[ignore = "diagnostic: prints GPOF-residual / kernel-decomposition numerics for mom-002 H2"]
fn mom_002_h2_gpof_residual_diagnostic() {
    let k0 = k0_at(F_HZ);
    let seed = thin_slab_guess(EPS_R, H_SUBSTRATE_M, k0);
    let (pole, iters) =
        newton_pole(SwChannel::Tm, seed, EPS_R, H_SUBSTRATE_M, k0).expect("Newton converges");
    let resid_d = d_tm(pole, EPS_R, H_SUBSTRATE_M, k0).norm();
    let res = residue(SwChannel::Tm, pole, EPS_R, H_SUBSTRATE_M, k0).expect("residue finite");
    let kz0_at_pole = k_z0(pole, k0);
    let kzd_at_pole = k_zd(pole, EPS_R, k0);

    eprintln!("--- mom-002 Hypothesis 2 GPOF residual diagnostic ---");
    eprintln!(
        "FR-4 / ε_r = {EPS_R}, h = {} mm, f = {} GHz, TM channel",
        H_SUBSTRATE_M * 1e3,
        F_HZ * 1e-9,
    );
    eprintln!("k_0      = {:.4} rad/m", k0);
    eprintln!("Newton   = {iters} iters, |D|={:.2e}", resid_d);
    eprintln!(
        "k_p      = {:.6}+j{:.6} rad/m  (k_p/k_0 = {:.6}+j{:.6})",
        pole.re,
        pole.im,
        pole.re / k0,
        pole.im / k0,
    );
    eprintln!(
        "k_z0(k_p)= {:.4e}+j{:.4e}  |k_z0|={:.4e}",
        kz0_at_pole.re,
        kz0_at_pole.im,
        kz0_at_pole.norm(),
    );
    eprintln!(
        "k_zd(k_p)= {:.4e}+j{:.4e}  |k_zd|={:.4e}",
        kzd_at_pole.re,
        kzd_at_pole.im,
        kzd_at_pole.norm(),
    );
    eprintln!(
        "Res_p    = {:.4e}+j{:.4e}  |Res|={:.4e}",
        res.re,
        res.im,
        res.norm(),
    );

    // 32-point sweep along the Aksun deformed contour: t ∈ [0, 10].
    eprintln!();
    eprintln!("--- Aksun contour samples (32 points, t ∈ [0, {AKSUN_T_MAX}]) ---");
    eprintln!(
        "{:>5} | {:>22} | {:>12} | {:>12} | {:>14} | {:>10}",
        "t", "k_ρ(t)", "|R|", "|R_pole|", "|R - R_pole|", "ratio",
    );
    eprintln!(
        "{:->6}+{:->24}+{:->14}+{:->14}+{:->16}+{:->12}",
        "", "", "", "", "", ""
    );
    let dt = AKSUN_T_MAX / ((N_CONTOUR_SAMPLES - 1) as f64);
    let mut sum_sq_r = 0.0_f64;
    let mut sum_sq_pole = 0.0_f64;
    let mut sum_sq_resid = 0.0_f64;
    for k in 0..N_CONTOUR_SAMPLES {
        let t = (k as f64) * dt;
        let k_z0_pt = Complex64::new(k0, 0.0) * Complex64::new(1.0, -t);
        let r_full = tm_reflection(k_z0_pt, EPS_R, H_SUBSTRATE_M, k0);
        let k_rho_pt = k_rho_on_contour(t, k0);
        let r_pole = res / (k_rho_pt - pole);
        let r_resid = r_full - r_pole;
        sum_sq_r += r_full.norm_sqr();
        sum_sq_pole += r_pole.norm_sqr();
        sum_sq_resid += r_resid.norm_sqr();
        let ratio = if r_full.norm() > 0.0 {
            r_resid.norm() / r_full.norm()
        } else {
            0.0
        };
        eprintln!(
            "{:>5.2} | {:>9.3e}+j{:>9.3e} | {:>12.4e} | {:>12.4e} | {:>14.4e} | {:>10.3e}",
            t,
            k_rho_pt.re,
            k_rho_pt.im,
            r_full.norm(),
            r_pole.norm(),
            r_resid.norm(),
            ratio,
        );
    }
    let norm_r = sum_sq_r.sqrt();
    let norm_pole = sum_sq_pole.sqrt();
    let norm_resid = sum_sq_resid.sqrt();
    eprintln!();
    eprintln!(
        "Sample-norm summary (ℓ² over 32 contour points):  \
         ‖R‖ = {norm_r:.4e},  ‖R_pole‖ = {norm_pole:.4e},  \
         ‖R − R_pole‖ = {norm_resid:.4e}"
    );
    let resid_fraction = if norm_r > 0.0 {
        norm_resid / norm_r
    } else {
        0.0
    };
    eprintln!("Residual fraction ‖R̃‖ / ‖R‖ = {resid_fraction:.4e}");
    eprintln!("  (= 1.0 → pole subtraction is a no-op; ≪ 1.0 → pole captures most of R)");

    // Note on the GPOF re-fit step: the `gpof::gpof` function is
    // `pub(crate)` inside yee-mom and not surfaced through the
    // `__internal` test helper module. Per the Track PPPPPP escape
    // hatch, we report the residual-fraction summary above instead of
    // attempting a literal GPOF re-fit from this integration-test
    // crate — the qualitative signal is the same: a residual fraction
    // ≪ 1 means the GPOF target collapsed, ≈ 1 means it did not.
    eprintln!();
    eprintln!(
        "(GPOF residual re-fit is `pub(crate)`-gated; the sample-norm \
         residual fraction above proxies the same quantity.)"
    );

    // -------------------------------------------------------------------
    // Kernel decomposition: hand-summed four-term vs full-Sommerfeld vs DCIM-only.
    // -------------------------------------------------------------------
    eprintln!();
    eprintln!("--- G^Φ at probe pair r1 = (0, 0, 0), r2 = (15 mm, 0, 0) ---");
    let r1 = Vector3::new(0.0, 0.0, 0.0);
    let r2 = Vector3::new(15.0e-3, 0.0, 0.0);

    let dcim_only = MultilayerGreens::new_microstrip_sommerfeld(EPS_R, H_SUBSTRATE_M, F_HZ, 5, 0);
    let sommerfeld = MultilayerGreens::new_microstrip_sommerfeld(EPS_R, H_SUBSTRATE_M, F_HZ, 5, 1);

    let g_dcim_only = dcim_only.scalar_scalar(r1, r2);
    let g_sommerfeld = sommerfeld.scalar_scalar(r1, r2);
    let (g_free, g_img, g_sw, g_hand) = hand_summed_scalar_scalar(&sommerfeld, r1, r2);

    eprintln!(
        "free-space term      = {:>12.5e}+j{:>12.5e}   |G| = {:.4e}",
        g_free.re,
        g_free.im,
        g_free.norm(),
    );
    eprintln!(
        "DCIM image train     = {:>12.5e}+j{:>12.5e}   |G| = {:.4e}",
        g_img.re,
        g_img.im,
        g_img.norm(),
    );
    eprintln!(
        "surface-wave Hankel  = {:>12.5e}+j{:>12.5e}   |G| = {:.4e}",
        g_sw.re,
        g_sw.im,
        g_sw.norm(),
    );
    eprintln!();
    eprintln!(
        "DCIM-only (n_sw=0)             = {:>12.5e}+j{:>12.5e}   |G| = {:.4e}",
        g_dcim_only.re,
        g_dcim_only.im,
        g_dcim_only.norm(),
    );
    eprintln!(
        "Sommerfeld (n_sw=1, full)      = {:>12.5e}+j{:>12.5e}   |G| = {:.4e}",
        g_sommerfeld.re,
        g_sommerfeld.im,
        g_sommerfeld.norm(),
    );
    eprintln!(
        "Hand-summed 4-term             = {:>12.5e}+j{:>12.5e}   |G| = {:.4e}",
        g_hand.re,
        g_hand.im,
        g_hand.norm(),
    );
    let delta_kernel = g_sommerfeld - g_hand;
    let rel_kernel = if g_sommerfeld.norm() > 0.0 {
        delta_kernel.norm() / g_sommerfeld.norm()
    } else {
        delta_kernel.norm()
    };
    eprintln!(
        "Δ (Sommerfeld − hand)          = {:>12.5e}+j{:>12.5e}   |Δ| = {:.4e}   rel = {:.4e}",
        delta_kernel.re,
        delta_kernel.im,
        delta_kernel.norm(),
        rel_kernel,
    );

    // -------------------------------------------------------------------
    // Verdict.
    // -------------------------------------------------------------------
    eprintln!();
    eprintln!("--- VERDICT ---");
    eprintln!("HYPOTHESIS 2 (GPOF residual fit):");
    eprintln!("   ‖R̃‖ / ‖R‖ on the 32-point Aksun contour = {resid_fraction:.4e}");
    let h2_verdict: &str = if resid_fraction < 0.1 {
        "STRONG — pole captures > 90% of R on the contour; GPOF re-fit \
         operates on a small smooth residual (good news for the fit)."
    } else if resid_fraction < 0.9 {
        "MIXED — pole captures some but not most of R; GPOF re-fit \
         retains most of the spectral content of the original samples."
    } else {
        "WEAK — pole subtraction barely changes the GPOF target. \
         Either the residue is too small relative to R on the contour, \
         or the pole sits far enough off-contour that the kernel \
         `Res / (k_ρ − k_p)` has little overlap with R(k_ρ) at the \
         sample points."
    };
    eprintln!("   verdict: {h2_verdict}");
    eprintln!();
    eprintln!("KERNEL SELF-CONSISTENCY (Sommerfeld vs hand-summed 4-term):");
    eprintln!(
        "   |Δ G^Φ| = {:.4e}   rel = {rel_kernel:.4e}",
        delta_kernel.norm()
    );
    let kernel_verdict: &str = if rel_kernel < 1.0e-6 {
        "CONSISTENT — the gap is not in the residue / Hankel / image \
         decomposition. Next diagnostic must look elsewhere (e.g. \
         Galerkin quadrature of the long-range Hankel term over the \
         strip, port-current normalisation, or the slab-reflection \
         convention that the residue is extracted from)."
    } else if rel_kernel < 1.0e-3 {
        "BORDERLINE — kernel agrees to GPOF round-off but not to \
         double-precision floor; suspicion is benign numerical \
         differences in how DCIM images are accumulated."
    } else {
        "INCONSISTENT — `scalar_scalar` disagrees materially with the \
         four-term hand-sum. Either the kernel adds a fifth term not \
         accounted for here, or the EEEEEE prefactor / branch / modal \
         profile is being applied differently in `surface_wave_sum` \
         than the spec records. This is the next thing to fix."
    };
    eprintln!("   verdict: {kernel_verdict}");
}
