//! mom-002 ψ_p + α_0 + port-current normalization audit — Track XXXXXX.
//!
//! ## Why this file exists
//!
//! Five prior diagnostics have shrunk-but-not-closed the mom-002 `|Z_in|`
//! gap on FR-4 / `h = 1.6 mm` / `w = 2.94 mm` / 1 GHz:
//!
//! * **Track EEEEEE** (commit `ca0e7bb`) — fixed the Sommerfeld
//!   surface-wave prefactor in
//!   [`yee_mom::multilayer::MultilayerGreens::surface_wave_sum`] to the
//!   canonical Michalski-Mosig 1997 form
//!   `G_sw = -(j/4) · Res · (k_ρ/k_z0) · H_0^{(2)}(k_ρ ρ) · ψ²`. Closes
//!   the bulk of the magnitude gap but leaves a ~30× residual.
//!
//! * **Track JJJJJJ** (commit `4dbeece`) — swept the strip length
//!   `L ∈ [15, 80] mm` and found only a ~5% monotonic trend. Hankel
//!   tail truncation is not the leading residual.
//!
//! * **Track PPPPPP** (commit `d89d0b9`) — GPOF-residual / kernel hand-
//!   sum decomposition. Pole subtraction is a no-op on the Aksun
//!   contour; the four-term hand-sum agrees bit-exact with
//!   `scalar_scalar`. Suspect was localised to `residue()`.
//!
//! * **Track SSSSSS** (commit `0e571b7`) — numerical contour integral
//!   of [`yee_mom::multilayer::slab_reflection`] around the converged
//!   pole detected a `-2×` discrepancy between the closed-form
//!   `residue()` and the contour residue.
//!
//! * **Track TTTTTT** (commit `a4f98a4`) — applied the
//!   Michalski-Mosig 1997 eq. (19) correction
//!   `Res = -N₁(k_p) / (2 · D'(k_p))` to [`yee_mom::sommerfeld::residue`].
//!   mom-002 `|Z_in|` dropped from `~2233 Ω` to `~2215 Ω`. Still ~43×
//!   above the analytic Hammerstad-Jensen `Z_0 ≈ 51 Ω`.
//!
//! The TTTTTT commit body lists three remaining candidates:
//!
//! * **(Y1) ψ_p(z) modal-profile normalization** inside
//!   [`yee_mom::multilayer::MultilayerGreens::surface_wave_sum`]. The
//!   helper `modal_z_profile` uses the textbook form
//!   `cos(k_zd (z + h)) / cos(k_zd h)`, peak-at-z=0 normalised. The
//!   open question is whether the dispersion-relation derivative-
//!   continuity identity at `z = 0` is honoured by this normalisation
//!   at the converged pole, and whether the profile is continuous
//!   across `z = 0`.
//!
//! * **(Y2) α_0 sign convention**. The pole-search uses
//!   `α_0 = -j · k_z0` with `Im k_z0 > 0` in the bound-mode regime, so
//!   `α_0` should be **real positive** (proper outward decay above the
//!   slab). If the residue / dispersion code paths land on the
//!   improper sheet (`Im k_z0 < 0`) the sign of α_0 flips and the
//!   surface-wave reconstruction picks up a global negative.
//!
//! * **(Y3) Galerkin port-current normalization at delta-gap
//!   excitation**. The port current `I_port = Σ_k b_k · i_k` is
//!   evaluated identically across kernels in
//!   [`yee_mom::__internal::z_in_with_greens`], but the value of `i`
//!   itself is kernel-dependent because `Z · i = b` is solved with a
//!   different `Z`. The question is whether the resulting `Z_in` shows
//!   any sign of a kernel-dependent **misnormalisation** — i.e. an
//!   extra factor on `I_port` that should not be there — or whether it
//!   is just the expected physics-dependent shift.
//!
//! ## Diagnostic method
//!
//! Three probes, in independent file regions:
//!
//! 1. **Probe 1 — ψ_p(z) audit.** Evaluate the textbook modal profile
//!    `cos(k_zd (z + h)) / cos(k_zd h)` (inside slab) and `exp(-α_0 z)`
//!    (above slab) at a handful of representative `z` values for the
//!    converged TM₀ pole at FR-4 / 1 GHz. Verify continuity at `z = 0`.
//!    Verify the dispersion-relation derivative-continuity identity
//!    `ε_r α_0 = k_zd · tan(k_zd h)` (TM bound mode) holds at the pole
//!    to machine precision.
//!
//! 2. **Probe 2 — α_0 sign sanity.** Compute `k_z0(k_p)` using the same
//!    `k_z0` helper the pole search uses. In the proper-sheet
//!    convention with `Im k_z0 > 0`, the air-side decay constant
//!    `α_0 = -j · k_z0` must be real positive. Print and check the
//!    sign.
//!
//! 3. **Probe 3 — Port-current normalization.** Build a single-strip
//!    mom-002 problem at `L = 30 mm` and call
//!    [`yee_mom::__internal::z_in_with_greens`] twice — once with the
//!    DCIM-only kernel (`new_microstrip_sommerfeld(.., 5, 0)`) and
//!    once with the Sommerfeld surface-wave kernel
//!    (`new_microstrip_sommerfeld(.., 5, 1)`). Print `Z_in` and the
//!    ratio. The Y3-flag criterion is: a Z-ratio that maps to a
//!    physically implausible `|I_port|` ratio (e.g. `1/Z_dcim`
//!    contradicting the analytic ~51 Ω target).
//!
//!    Per the brief's escape hatch, `|I_port|` itself is not directly
//!    surfaced by `z_in_with_greens` and the lower-level
//!    `impedance_matrix` / `delta_gap_rhs` are `pub(crate)`. We report
//!    `Z_in` only; the Z-ratio is sufficient to flag Y3 (`Z = V/I`, so
//!    a kernel-dependent normalisation on `I` shows as a kernel-
//!    dependent shift on `Z`).
//!
//! ## References
//!
//! * Track EEEEEE — `crates/yee-mom/tests/sommerfeld_residue_diagnostic.rs`.
//! * Track JJJJJJ — `crates/yee-mom/tests/mom_002_extent_sensitivity.rs`.
//! * Track PPPPPP — `crates/yee-mom/tests/mom_002_h2_gpof_diagnostic.rs`.
//! * Track SSSSSS — `crates/yee-mom/tests/mom_002_reflection_convention.rs`.
//! * Track TTTTTT — commit `a4f98a4`
//!   (`Merge Track TTTTTT: yee-mom — fix sommerfeld::residue() sign + factor-of-2`).
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7,
//!   eq. 3.196–3.199.
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans.
//!   Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_mom::__internal::{MultilayerGreens, z_in_with_greens};
use yee_mom::sommerfeld::{SwChannel, k_z0, k_zd, newton_pole, thin_slab_guess};

const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const F_HZ: f64 = 1.0e9;
const STRIP_W_M: f64 = 2.94e-3;
const STRIP_L_M: f64 = 30.0e-3;
const N_LENGTH: usize = 30;
const N_WIDTH: usize = 16;
const N_DCIM_IMAGES: usize = 5;

fn k0_at(freq_hz: f64) -> f64 {
    std::f64::consts::TAU * freq_hz / yee_core::units::C0
}

/// Textbook TM₀ modal `z`-profile inside the grounded dielectric slab,
/// normalised to unit value at the slab top (`z = 0`).
///
/// ```text
///   ψ(z) = cos(k_zd · (z + h)) / cos(k_zd · h),   z ∈ [-h, 0]
///   ψ(z) = exp(-α_0 · z),                         z > 0
/// ```
///
/// At `z = 0` both forms collapse to `1`, so the profile is continuous
/// (provided the dispersion-relation derivative-continuity identity
/// `ε_r α_0 = k_zd · tan(k_zd · h)` holds at the pole; if not, the
/// `∂ψ/∂z` jump at `z = 0` is non-zero — a finding the diagnostic
/// surfaces).
///
/// Replicated inline because the internal `modal_z_profile` helper in
/// [`yee_mom::multilayer`] is not surfaced through `__internal`.
fn psi_textbook(z: f64, k_zd_pole: Complex64, alpha0: Complex64, h: f64) -> Complex64 {
    if z <= 0.0 {
        let arg_top = k_zd_pole * Complex64::new(h, 0.0);
        let denom = arg_top.cos();
        let arg_field = k_zd_pole * Complex64::new(z + h, 0.0);
        arg_field.cos() / denom
    } else {
        (-alpha0 * Complex64::new(z, 0.0)).exp()
    }
}

/// Edge-clustered (Chebyshev-y) strip mesh — bit-for-bit equivalent to
/// `yee_validation::mom_002_strip_mesh_with_spacing` with
/// `StripSpacing::EdgeClustered`. Inlined here so the diagnostic has
/// no cross-lane dependency on `yee-validation` internals (the
/// pattern follows `tests/mom_002_extent_sensitivity.rs`).
fn build_strip_mesh_edge_clustered(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
) -> yee_mesh::TriMesh {
    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(nx * ny);
    let dx = length_m / (n_length as f64);
    let y_nodes: Vec<f64> = (0..=n_width)
        .map(|j| {
            let theta = std::f64::consts::PI * (j as f64) / (n_width as f64);
            -(width_m / 2.0) * theta.cos()
        })
        .collect();
    for i in 0..nx {
        let x = (i as f64) * dx;
        for &y in &y_nodes {
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }
    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_length * n_width);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_length * n_width);
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == 0 {
                1
            } else if i == 1 {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    yee_mesh::TriMesh::new(vertices, triangles, tags).expect("strip mesh invariants")
}

/// Print the three probe tables and the verdict. Marked `#[ignore]` so
/// the suite never runs it by default; invoke explicitly via
///
/// ```text
/// cargo test -p yee-mom --release --test mom_002_psi_port_audit \
///     -- --ignored --nocapture
/// ```
///
/// to dump the table and the Y1 / Y2 / Y3 verdict.
#[test]
#[ignore = "diagnostic: prints ψ_p + α_0 + port-current normalisation audit for mom-002"]
fn mom_002_psi_port_audit_diagnostic() {
    let k0 = k0_at(F_HZ);

    // Converge the TM₀ pole — same path the production
    // `find_surface_wave_poles` helper takes.
    let seed = thin_slab_guess(EPS_R, H_SUBSTRATE_M, k0);
    let (pole, iters) =
        newton_pole(SwChannel::Tm, seed, EPS_R, H_SUBSTRATE_M, k0).expect("Newton converges");

    // Modal `k_zd` and air-side `α_0` at the pole. The dispersion-
    // relation derivative-continuity identity (TM bound mode) is
    //
    //     ε_r α_0  =  k_zd · tan(k_zd · h).
    //
    // Equivalently, this is `d_tm(k_p) = 0` (the Newton residue should
    // be `< NEWTON_TOL`).
    let kzd_pole = k_zd(pole, EPS_R, k0);
    let kz0_pole = k_z0(pole, k0);
    let j = Complex64::new(0.0, 1.0);
    let alpha0 = -j * kz0_pole;

    eprintln!("--- Track XXXXXX: ψ_p audit + α_0 sign + port-current normalization ---");
    eprintln!(
        "FR-4 / ε_r = {EPS_R}, h = {} mm, f = {} GHz, TM channel",
        H_SUBSTRATE_M * 1e3,
        F_HZ * 1e-9,
    );
    eprintln!(
        "Newton: {iters} iters, k_p/k_0 = {:.6} + j·{:.6e}",
        pole.re / k0,
        pole.im / k0,
    );
    eprintln!(
        "k_zd(k_p) = {:.6e} + j·{:.6e},   α_0 = {:.6e} + j·{:.6e}",
        kzd_pole.re, kzd_pole.im, alpha0.re, alpha0.im,
    );

    // -----------------------------------------------------------------
    // Probe 1 — ψ_p(z) profile values
    // -----------------------------------------------------------------
    eprintln!();
    eprintln!("Probe 1 — ψ_p(z) profile values:");
    eprintln!(
        "  {:>14}  {:>26}  {:>26}",
        "z (m)", "ψ_p(z) re", "ψ_p(z) im",
    );

    let zs = [
        ("slab top z = 0", 0.0_f64),
        ("z = -h/4", -H_SUBSTRATE_M / 4.0),
        ("z = -h/2", -H_SUBSTRATE_M / 2.0),
        ("z = -h (ground)", -H_SUBSTRATE_M),
        ("z = +1e-6 (above)", 1.0e-6),
    ];
    let mut psi_values = std::collections::HashMap::new();
    for (label, z) in zs.iter() {
        let psi = psi_textbook(*z, kzd_pole, alpha0, H_SUBSTRATE_M);
        psi_values.insert(*label, psi);
        eprintln!(
            "  {label:<20}  z = {z:>14.6e}  ψ = {:>14.6e} + j·{:>14.6e}",
            psi.re, psi.im
        );
    }

    // Continuity at z=0: the slab-side form ψ(0) = cos(k_zd h) /
    // cos(k_zd h) = 1, and the air-side form lim_{z→0+} exp(-α_0 z) =
    // 1. They must match.
    let psi_zero_slab_side = psi_textbook(-1.0e-12, kzd_pole, alpha0, H_SUBSTRATE_M);
    let psi_zero_air_side = psi_textbook(1.0e-12, kzd_pole, alpha0, H_SUBSTRATE_M);
    let continuity_diff = (psi_zero_slab_side - psi_zero_air_side).norm();
    let continuity_ok = continuity_diff < 1e-6;
    eprintln!(
        "  Continuity at z = 0: ψ(0-) = {:.6e} + j·{:.6e}, ψ(0+) = {:.6e} + j·{:.6e}",
        psi_zero_slab_side.re, psi_zero_slab_side.im, psi_zero_air_side.re, psi_zero_air_side.im,
    );
    eprintln!(
        "    |ψ(0+) - ψ(0-)| = {continuity_diff:.6e}    {{{}}}",
        if continuity_ok { "OK" } else { "MISMATCH" },
    );

    // Dispersion identity at k_p (TM bound mode):
    //   LHS  =  ε_r · α_0
    //   RHS  =  k_zd · tan(k_zd · h)
    // Must agree to machine precision when Newton has converged.
    let lhs = Complex64::new(EPS_R, 0.0) * alpha0;
    let arg = kzd_pole * Complex64::new(H_SUBSTRATE_M, 0.0);
    let rhs = kzd_pole * arg.tan();
    let dispersion_diff = (lhs - rhs).norm();
    let dispersion_ok = dispersion_diff / lhs.norm().max(1e-300) < 1e-6;
    eprintln!("  Dispersion identity ε_r α_0 = k_zd · tan(k_zd · h) at k_p:",);
    eprintln!("    LHS = {:.6e} + j·{:.6e}", lhs.re, lhs.im,);
    eprintln!("    RHS = {:.6e} + j·{:.6e}", rhs.re, rhs.im,);
    eprintln!(
        "    |LHS - RHS| / |LHS| = {:.6e}    {{{}}}",
        dispersion_diff / lhs.norm().max(1e-300),
        if dispersion_ok { "OK" } else { "MISMATCH" },
    );

    // -----------------------------------------------------------------
    // Probe 2 — α_0 sign sanity
    // -----------------------------------------------------------------
    eprintln!();
    eprintln!("Probe 2 — α_0 sign:");
    eprintln!(
        "  k_z0(k_p) = {:.6e} + j·{:.6e}    → α_0 = -j·k_z0 = {:.6e} + j·{:.6e}",
        kz0_pole.re, kz0_pole.im, alpha0.re, alpha0.im,
    );
    // For a bound surface mode: Im k_z0 > 0 (proper sheet),
    // α_0 = -j k_z0 has Re α_0 > 0 (real-positive air-side decay).
    // We tolerate a tiny imaginary part from numerical sqrt rounding.
    let alpha0_real_positive = alpha0.re > 0.0 && alpha0.im.abs() < 1e-6 * alpha0.re.abs().max(1.0);
    let on_proper_sheet = kz0_pole.im > 0.0;
    eprintln!(
        "  on proper sheet (Im k_z0 > 0)? {}    α_0 real-positive? {}",
        if on_proper_sheet { "yes" } else { "no" },
        if alpha0_real_positive { "yes" } else { "no" },
    );
    eprintln!(
        "  {{{}}}",
        if alpha0_real_positive && on_proper_sheet {
            "real-positive air decay (OK)"
        } else {
            "wrong sign / wrong sheet"
        },
    );

    // Also sample slab_reflection-like denominators at a couple of
    // off-pole `k_z0` values to confirm `k_z0` returns the principal
    // branch and the principal branch behaves monotonically over a
    // sane real-axis range.
    eprintln!("  k_z0 sanity samples (off-pole):");
    for &kr_over_k0 in &[0.5_f64, 0.9, 1.0, 1.001, 1.5, 2.0] {
        let kr = Complex64::new(kr_over_k0 * k0, 0.0);
        let kz0 = k_z0(kr, k0);
        eprintln!(
            "    k_ρ/k_0 = {kr_over_k0:>5.3}   k_z0/k_0 = {:.6e} + j·{:.6e}",
            kz0.re / k0,
            kz0.im / k0,
        );
    }

    // -----------------------------------------------------------------
    // Probe 3 — port-current normalisation
    // -----------------------------------------------------------------
    eprintln!();
    eprintln!(
        "Probe 3 — port-current normalization (mom-002, L = {} mm):",
        STRIP_L_M * 1e3
    );

    let mesh = build_strip_mesh_edge_clustered(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH);
    let port_tag = 1u32;

    // Kernel A — DCIM-only, no surface wave.
    let greens_dcim = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ,
        N_DCIM_IMAGES,
        0, // no surface-wave poles
    );
    let z_dcim = z_in_with_greens(&mesh, port_tag, &greens_dcim).expect("DCIM solve");

    // Kernel B — DCIM + Sommerfeld TM₀ surface wave.
    let greens_som = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ,
        N_DCIM_IMAGES,
        1, // one surface-wave pole (TM₀)
    );
    let z_som = z_in_with_greens(&mesh, port_tag, &greens_som).expect("Sommerfeld solve");

    eprintln!(
        "  DCIM-only:    Z_in = {:>12.4e} + j·{:>12.4e}   |Z_in| = {:.4e}",
        z_dcim.re,
        z_dcim.im,
        z_dcim.norm(),
    );
    eprintln!(
        "  Sommerfeld:   Z_in = {:>12.4e} + j·{:>12.4e}   |Z_in| = {:.4e}",
        z_som.re,
        z_som.im,
        z_som.norm(),
    );

    // V_port = 1.0 in both cases (delta-gap), so the ratio of `Z_in`s
    // is the inverse of the `I_port` ratio:
    //   Z_som / Z_dcim  =  I_dcim / I_som.
    //
    // A kernel-independent port-current normalisation should produce
    // a ratio whose magnitude is roughly the same order as the
    // physical `Z_in` shift. A ratio either ~1 (kernel has no effect,
    // which would contradict the EEEEEE prefactor fix) or
    // pathologically large/small (e.g. orders of magnitude beyond the
    // analytic `|Z_dcim| / 51 Ω` ratio) flags Y3.
    let z_ratio = z_som / z_dcim;
    let i_ratio_inferred = Complex64::new(1.0, 0.0) / z_ratio; // I_som / I_dcim, since V=1
    eprintln!(
        "  Z_som / Z_dcim:   |.| = {:.6e}    arg = {:.6e} rad    (= I_dcim / I_som)",
        z_ratio.norm(),
        z_ratio.arg(),
    );
    eprintln!(
        "  inferred I_som/I_dcim:   |.| = {:.6e}    arg = {:.6e} rad",
        i_ratio_inferred.norm(),
        i_ratio_inferred.arg(),
    );

    // Y3 flag: a kernel-dependent port-current normalisation would
    // show up either as a perfectly 1.0 ratio (no surface-wave
    // physics making it through to I — implausible given EEEEEE's
    // fix) or as an extreme ratio that does not correspond to a
    // sensible shift toward the analytic 51 Ω target. The expected
    // behaviour: |Z_dcim| ≈ 2200 Ω and |Z_som| ≈ 2200 Ω (almost
    // unchanged — at L = 30 mm the surface-wave correction is
    // small), so a ratio close to but not exactly 1 is expected.
    // Any ratio above 100 or below 0.01 would be a Y3 smoking gun.
    let y3_extreme = z_ratio.norm() > 100.0 || z_ratio.norm() < 0.01;

    // -----------------------------------------------------------------
    // Verdict
    // -----------------------------------------------------------------
    eprintln!();
    eprintln!("Verdict:");

    // Y1 — ψ normalisation. Detected if either continuity or
    // dispersion-identity check failed.
    let y1_detected = !(continuity_ok && dispersion_ok);
    eprintln!(
        "  Y1 (ψ normalization):     {}",
        if y1_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!("      reason: continuity_ok = {continuity_ok}, dispersion_ok = {dispersion_ok}",);

    // Y2 — α_0 sign. Detected if α_0 is not real-positive or k_z0 is
    // on the improper sheet.
    let y2_detected = !(alpha0_real_positive && on_proper_sheet);
    eprintln!(
        "  Y2 (α_0 sign):            {}",
        if y2_detected {
            "detected"
        } else {
            "not detected"
        },
    );
    eprintln!(
        "      reason: α_0_real_positive = {alpha0_real_positive}, on_proper_sheet = {on_proper_sheet}",
    );

    // Y3 — port-current normalisation. Inconclusive without direct
    // access to `|I_port|` (the `__internal` surface exposes only
    // `Z_in`, per the brief's escape hatch). We flag "detected" only
    // for an extreme `Z_som / Z_dcim` magnitude. Anything in the
    // physically plausible band (0.01 < |.| < 100) is reported as
    // "inconclusive" so future tracks know to plumb `|I_port|`
    // through the test surface before declaring Y3 dead.
    let y3_verdict = if y3_extreme {
        "detected"
    } else {
        "inconclusive"
    };
    eprintln!("  Y3 (port normalization):  {y3_verdict}",);
    eprintln!(
        "      reason: |Z_som / Z_dcim| = {:.6e} (extreme outside [0.01, 100] flags Y3;",
        z_ratio.norm(),
    );
    eprintln!(
        "      direct |I_port| access blocked by `pub(crate)` on `delta_gap_rhs`/`impedance_matrix`)",
    );
}
