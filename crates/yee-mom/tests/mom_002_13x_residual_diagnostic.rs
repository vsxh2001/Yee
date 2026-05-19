//! mom-002 ~13× capacitive-reactance residual diagnostic — Track MMMMMMM.
//!
//! ## Why this file exists
//!
//! Track IIIIIII (commit `65502aa`) landed the ADR-0036 reframe of
//! mom-002: `L = 82 mm` ≈ `λ_eff / 2` half-wave resonator on FR-4 at 1
//! GHz, centered port, uniform y-spacing, 82 × 16 cells. The empirical
//! landing is
//!
//! ```text
//!   Z_in  =  +1.819  +  j (-674.105)   Ω,
//!   |Z_in|  =  674.108  Ω  ≈  13.2 × Z_0.
//! ```
//!
//! Improvement vs the pre-reframe 30-mm strip: `|Z_in|` dropped from
//! `~2569 Ω` (43× Z_0) to `~674 Ω` (13× Z_0). `Re(Z)` is now physically
//! clean and positive. **Dominant residual is the capacitive reactance
//! `Im(Z) = -674 Ω`**.
//!
//! A perfectly-tuned `λ_eff / 2` microstrip line at 1 GHz on FR-4 with
//! `Z_0 ≈ 51 Ω` should land near `Z_in ≈ +Z_0 ≈ 51 Ω` (resistive at exact
//! resonance) for a short-circuited termination, or near a large real
//! `|Z_in| ≈ Z_0² / R_load` for an open-circuited termination. Negative
//! `Im(Z)` of magnitude `~13 × Z_0` indicates the line behaves as
//! **electrically shorter than `λ_eff / 2`** at 1 GHz, i.e. the solver
//! perceives a larger `ε_eff` than the Hammerstad-Jensen analytic
//! `ε_eff ≈ 3.32`, or the strip is just enough off-resonance that the
//! open-circuit-like reactance still dominates.
//!
//! ## Three candidates this file tests
//!
//! * **(R1) `ε_eff` biasing**. The DCIM image train + Sommerfeld
//!   surface-wave kernel could be producing an effective dielectric
//!   constant different from H-J's static `ε_eff`. If
//!   `ε_eff_solver > 3.32`, the strip's `λ_eff / 2 < 82 mm`, the line
//!   looks capacitive at 1 GHz.
//!
//! * **(R2) Resonance-frequency offset**. Sweep `f` near 1 GHz. The
//!   resonance is where `Im(Z)` crosses zero; the offset from 1 GHz
//!   pins the `ε_eff` bias from R1 quantitatively (`f_res / 1 GHz =
//!   √(3.32 / ε_eff_solver)`).
//!
//! * **(R3) Finite-width edge effects**. The 2.94 mm strip is a
//!   non-trivial fraction of substrate thickness 1.6 mm (`w/h ≈ 1.84`).
//!   Edge fringing usually pushes the **other** direction (raises
//!   capacitance, lowers `ε_eff`). Vary `w ∈ {2, 2.94, 4} mm` at the
//!   reference frequency and report the trend on `|Z_in|`.
//!
//! ## Escape-hatch usage
//!
//! The brief allows dropping the frequency sweep to 3 points and the
//! width sweep to 3 widths "if `z_in_with_greens` is slow". On this
//! hardware the production 82×16 mesh measured by
//! [`yee_validation::tests::mom_002_measure_z_in_for_seeding`] takes
//! `~270 s` per solve, so a 5-frequency × 5-width grid is in the hours,
//! not minutes. This diagnostic therefore (a) uses a coarser
//! `41 × 8` cell mesh and (b) walks 3 frequencies × 3 widths. Both
//! reductions surrender precision on the absolute landing — `41 × 8`
//! drops basis count `~3.7×` from the production mesh — but they
//! preserve the **trends** which are the signal R1 / R2 / R3 want.
//! Each verdict block calls this out explicitly.
//!
//! ## References
//!
//! * ADR-0036 — `docs/src/decisions/0036-mom-002-validation-strategy.md`
//!   (reframe to half-wave resonator).
//! * Track IIIIIII commit `65502aa` — `MOM_002_STRIP_LENGTH_M = 82e-3`,
//!   `MOM_002_N_LENGTH = 82`, `MOM_002_N_WIDTH = 16`, centered port,
//!   uniform y-spacing, `MOM_002_Z_IN_MEASURED_OHM = 674.108 Ω`.
//! * Sibling diagnostics (XXXXXX/SSSSSS/TTTTTT/PPPPPP/JJJJJJ/EEEEEE)
//!   under `crates/yee-mom/tests/mom_002_*.rs` for the pattern.
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7 (microstrip
//!   `Z_0` and `ε_eff`), §2.5 (transmission-line input impedance).
//! * E. Hammerstad and Ø. Jensen, "Accurate Models for Microstrip
//!   Computer-Aided Design," *MTT-S Digest*, 1980.

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_mom::__internal::{Greens, MultilayerGreens, z_in_with_greens};

// FR-4 / 1 GHz canonical microstrip parameters (same as
// `yee-validation::MOM_002_*` constants).
const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const STRIP_W_M: f64 = 2.94e-3;
const STRIP_L_M: f64 = 82.0e-3;
const F_HZ_NOMINAL: f64 = 1.0e9;
const Z0_REF: f64 = 50.0;

// Sommerfeld-kernel parameters (match the production headline gate).
const N_DCIM_IMAGES: usize = 5;
const N_SW_POLES: usize = 1;

// Coarse-mesh dimensions for this diagnostic — production is 82×16
// (~270s per solve), so we reduce to 41×8 (~16× cheaper LU) to keep
// the multi-frequency / multi-width sweeps tractable. Both numbers are
// even so the centered-port column straddle still works.
const N_LENGTH: usize = 40;
const N_WIDTH: usize = 8;

// H-J analytic ε_eff for FR-4 / 1.6 mm / w = 2.94 mm. Source: Pozar 4e
// §3.7 eq. 3.195 evaluated at u = w/h = 1.8375. This is the value the
// reframe targeted; if R1 detects a `solver` ε_eff materially above
// this, the kernel is overshooting on the dielectric loading.
const EPS_EFF_HJ_ANALYTIC: f64 = 3.32;

/// Speed-of-light helper — duplicated from sibling diagnostics so the
/// constant lookup stays inline.
fn k0_at(freq_hz: f64) -> f64 {
    std::f64::consts::TAU * freq_hz / yee_core::units::C0
}

/// Build the ADR-0036 centered-uniform strip mesh at the supplied
/// dimensions. Bit-equivalent to
/// `yee_validation::mom_002_strip_mesh_with_spacing(.., StripSpacing::Uniform)`
/// but inlined here so the diagnostic has no cross-lane dependency on
/// `yee-validation` internals (pattern follows
/// `tests/mom_002_psi_port_audit.rs`).
fn build_strip_mesh_centered_uniform(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
) -> yee_mesh::TriMesh {
    assert!(
        n_length >= 4 && n_length.is_multiple_of(2),
        "n_length must be even and >= 4 to host a centered port column"
    );
    assert!(n_width >= 1, "n_width must be >= 1");

    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(nx * ny);
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;

    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..=n_width {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_length * n_width);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_length * n_width);
    let port_left = n_length / 2 - 1;
    let port_right = n_length / 2;
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == port_left {
                1
            } else if i == port_right {
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

/// Solve mom-002 at the supplied geometry / frequency through the
/// production Sommerfeld kernel. Returns `Z_in = V_port / I_port`.
fn z_in_at(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
    freq_hz: f64,
) -> Complex64 {
    let mesh = build_strip_mesh_centered_uniform(length_m, width_m, n_length, n_width);
    let port_tag = 1u32;
    let greens = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        freq_hz,
        N_DCIM_IMAGES,
        N_SW_POLES,
    );
    z_in_with_greens(&mesh, port_tag, &greens).expect("z_in_with_greens")
}

/// Probe R1 — extract `ε_eff_solver` from the Green's-function
/// scalar-potential spatial decay along the strip axis.
///
/// Sample `G_φ(ρ) = scalar_scalar((0, 0, 0), (ρ, 0, 0))` at two well-
/// separated `ρ` values on the slab top. In the surface-wave-dominated
/// far field (`k_ρ ρ ≫ 1`) the kernel reduces to
///
/// ```text
///   G_sw  ~  H_0^{(2)}(k_p ρ)  ~  exp(-j k_p ρ) / √(k_p ρ),
/// ```
///
/// so the **phase** of `G_φ` advances at rate `-k_p` along `ρ` and the
/// extracted phase velocity yields `(k_p / k_0)² = ε_eff_solver`. We
/// fit `k_p` from a linear regression of `arg(G_φ(ρ_k))` against `ρ_k`,
/// which automatically robustifies against the Hankel branch noise the
/// `H_0^{(2)}` evaluator carries near its asymptotic regime.
fn probe_r1_eps_eff_from_greens() -> (f64, bool) {
    eprintln!();
    eprintln!("Probe R1 — solver effective dielectric (from G_φ spatial decay):");

    let k0 = k0_at(F_HZ_NOMINAL);
    let greens = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ_NOMINAL,
        N_DCIM_IMAGES,
        N_SW_POLES,
    );

    // Sample ρ from ~λ_eff / 4 to ~3 λ_eff / 4 along the strip axis.
    // For ε_eff = 3.32, λ_eff at 1 GHz is ~165 mm so ρ ∈ [40, 120] mm.
    // The substrate is h = 1.6 mm and the strip is 82 mm, so 120 mm is
    // already off-strip — we sample G_φ at field points NOT colocated
    // with the strip; this is the kernel function itself, not the
    // basis-function projection. ρ ≥ 6h ensures we're in the
    // surface-wave-dominated far-field regime where the free-space
    // and DCIM-image terms have decayed by a couple of orders of
    // magnitude relative to the TM₀ Hankel.
    let rho_values: Vec<f64> = (0..8).map(|k| 40.0e-3 + (k as f64) * 12.0e-3).collect();
    let r_source = Vector3::new(0.0, 0.0, 0.0);

    eprintln!(
        "  {:>10}  {:>14}  {:>14}  {:>14}",
        "ρ (mm)", "|G_φ|", "arg(G_φ) (rad)", "k_eff*ρ (rad)",
    );
    let mut samples: Vec<(f64, f64)> = Vec::with_capacity(rho_values.len());
    let mut last_arg_unwrapped: Option<f64> = None;
    let mut last_arg_raw: Option<f64> = None;
    for &rho in &rho_values {
        let r_field = Vector3::new(rho, 0.0, 0.0);
        let g = greens.scalar_scalar(r_field, r_source);
        // Phase unwrapping — the Hankel evaluator returns arg in
        // (-π, π], but the underlying e^{-jk_p ρ} drifts monotonically
        // (negatively) with ρ. We walk samples in increasing-ρ order
        // and add 2π whenever the raw arg jumps by more than π upward
        // (the unwrapped sequence should be monotonically decreasing).
        let raw_arg = g.arg();
        let unwrapped = match last_arg_raw {
            None => raw_arg,
            Some(prev_raw) => {
                let prev_un = last_arg_unwrapped.unwrap();
                let mut step = raw_arg - prev_raw;
                while step > std::f64::consts::PI {
                    step -= std::f64::consts::TAU;
                }
                while step < -std::f64::consts::PI {
                    step += std::f64::consts::TAU;
                }
                prev_un + step
            }
        };
        last_arg_raw = Some(raw_arg);
        last_arg_unwrapped = Some(unwrapped);
        eprintln!(
            "  {:>10.3}  {:>14.4e}  {:>14.6}  {:>14.6}",
            rho * 1e3,
            g.norm(),
            unwrapped,
            -unwrapped,
        );
        samples.push((rho, unwrapped));
    }

    // Linear-least-squares fit: arg(ρ) ≈ -k_p · ρ + φ_0. Slope = -k_p.
    let n = samples.len() as f64;
    let sum_x: f64 = samples.iter().map(|(x, _)| *x).sum();
    let sum_y: f64 = samples.iter().map(|(_, y)| *y).sum();
    let sum_xx: f64 = samples.iter().map(|(x, _)| x * x).sum();
    let sum_xy: f64 = samples.iter().map(|(x, y)| x * y).sum();
    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    let k_eff = -slope;
    let eps_eff_solver = (k_eff / k0).powi(2);
    let eps_eff_ratio = eps_eff_solver / EPS_EFF_HJ_ANALYTIC;
    eprintln!(
        "  Linear-fit slope: dφ/dρ = {:.6e} rad/m  →  k_eff = {:.6e} rad/m",
        slope, k_eff,
    );
    eprintln!(
        "  k_eff / k_0 = {:.4},   ε_eff_solver = (k_eff/k_0)² = {:.4}",
        k_eff / k0,
        eps_eff_solver,
    );
    eprintln!("  ε_eff_analytic (Hammerstad-Jensen) = {EPS_EFF_HJ_ANALYTIC}");
    eprintln!(
        "  ratio ε_eff_solver / ε_eff_analytic = {:.4} ({:+.2} %)",
        eps_eff_ratio,
        (eps_eff_ratio - 1.0) * 100.0,
    );

    // Verdict: > 10 % deviation in either direction flags the bias.
    //
    // Direction map:
    //   ε_eff_solver > H-J  ⇒  λ_eff_solver < λ_eff_HJ
    //     ⇒  fixed L = 82 mm is electrically LONGER (more wavelengths fit)
    //     ⇒  resonance (L = λ_eff / 2) moves DOWN in frequency.
    //   ε_eff_solver < H-J  ⇒  λ_eff_solver > λ_eff_HJ
    //     ⇒  fixed L = 82 mm is electrically SHORTER (fewer wavelengths fit)
    //     ⇒  resonance moves UP in frequency, and at 1 GHz the line is
    //         pre-resonance / open-circuit-like / capacitive (Im(Z) < 0) —
    //         exactly the IIIIIII landing.
    let biased = (eps_eff_ratio - 1.0).abs() > 0.10;
    if biased {
        eprintln!(
            "  Verdict: biased — ε_eff_solver disagrees with H-J by > 10 % \
             (sign indicates {})",
            if eps_eff_solver > EPS_EFF_HJ_ANALYTIC {
                "solver overestimates dielectric loading (line electrically LONGER; resonance moves DOWN)"
            } else {
                "solver underestimates dielectric loading (line electrically SHORTER; resonance moves UP; 1 GHz is pre-resonance / capacitive)"
            },
        );
    } else {
        eprintln!("  Verdict: matches — solver ε_eff is within 10 % of H-J");
    }
    (eps_eff_solver, biased)
}

/// Probe R2 — frequency sweep across [0.7, 1.3] GHz at the production
/// 82 mm strip / coarse mesh. Reports `|Z_in|`, `Re(Z)`, `Im(Z)` per
/// frequency and tries to detect the `Im(Z) = 0` zero-crossing (the
/// resonance) by linear interpolation. The escape-hatch 3-point set
/// `{0.7, 1.0, 1.3} GHz` is what we use here — the trend is the signal,
/// not the exact resonance.
fn probe_r2_frequency_sweep() -> (Option<f64>, Vec<(f64, Complex64)>) {
    eprintln!();
    eprintln!(
        "Probe R2 — Frequency sweep ({}x{} coarse mesh, L = {} mm, w = {} mm):",
        N_LENGTH,
        N_WIDTH,
        STRIP_L_M * 1e3,
        STRIP_W_M * 1e3,
    );
    eprintln!(
        "  {:>10}  {:>14}  {:>14}  {:>14}",
        "f (GHz)", "Re(Z) (Ω)", "Im(Z) (Ω)", "|Z| (Ω)",
    );

    let freqs_hz: Vec<f64> = vec![0.7e9, 1.0e9, 1.3e9];
    let mut points: Vec<(f64, Complex64)> = Vec::with_capacity(freqs_hz.len());
    for &f in &freqs_hz {
        let z = z_in_at(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH, f);
        eprintln!(
            "  {:>10.3}  {:>14.4}  {:>14.4}  {:>14.4}",
            f * 1e-9,
            z.re,
            z.im,
            z.norm(),
        );
        points.push((f, z));
    }

    // Detect Im(Z) zero crossing. Walk consecutive (f_i, Im(Z_i))
    // pairs; if Im changes sign with positive slope (capacitive → 0 →
    // inductive going up in f), linearly interpolate to find f_res.
    // For a short open-circuited line, Im(Z) is large-negative below
    // λ/4 and large-positive between λ/4 and λ/2 — the zero-crossing
    // at λ/4 is the parallel resonance; the next (positive-slope)
    // crossing at λ/2 is the series resonance. We want the first
    // positive-slope crossing.
    let mut f_res: Option<f64> = None;
    for win in points.windows(2) {
        let (f0, z0) = win[0];
        let (f1, z1) = win[1];
        if z0.im < 0.0 && z1.im > 0.0 {
            let f_interp = f0 + (f1 - f0) * (-z0.im) / (z1.im - z0.im);
            f_res = Some(f_interp);
            break;
        }
    }

    match f_res {
        Some(f) => {
            let delta_ghz = (f - F_HZ_NOMINAL) * 1e-9;
            let f_ratio = f / F_HZ_NOMINAL;
            // Implied half-wave: λ_eff(f_res) = c / (f_res · √ε_eff).
            // If the line is a half-wave resonator at f_res, then
            // L = λ_eff(f_res) / 2, hence ε_eff_at_fres = (c / (2 L
            // f_res))². We just print this as an alternative bias
            // estimator that doesn't depend on R1's Green's-function
            // sampling.
            let lambda_eff_2_at_fres = yee_core::units::C0 / (2.0 * f * EPS_EFF_HJ_ANALYTIC.sqrt());
            eprintln!(
                "  Im(Z) zero-crossing detected at f_res ≈ {:.4} GHz (Δf = {:+.4} GHz)",
                f * 1e-9,
                delta_ghz,
            );
            eprintln!(
                "  f_res / 1 GHz = {:.4}  →  implied ε_eff = 3.32 · (1 / {:.4})² = {:.4}",
                f_ratio,
                f_ratio,
                EPS_EFF_HJ_ANALYTIC / (f_ratio * f_ratio),
            );
            eprintln!(
                "  At H-J ε_eff = 3.32, λ_eff/2 at f_res would be {:.3} mm; strip length is {:.3} mm.",
                lambda_eff_2_at_fres * 1e3,
                STRIP_L_M * 1e3,
            );
            eprintln!("  Verdict: resonance offset detected");
        }
        None => {
            eprintln!(
                "  Im(Z) does NOT change sign across [0.7, 1.3] GHz — \
                 resonance is outside this sweep range."
            );
            // Report whether Im(Z) is uniformly negative (line still
            // "too short" / capacitive throughout) or uniformly
            // positive (line "too long" / inductive throughout). This
            // tells the direction even without a zero crossing.
            //
            // f_res = c / (2 L √ε_eff), so:
            //   resonance ABOVE 1.3 GHz  ⇒  ε_eff_solver < H-J (smaller
            //     dielectric loading; line electrically shorter at 1 GHz;
            //     Im(Z) < 0 capacitive).
            //   resonance BELOW 0.7 GHz  ⇒  ε_eff_solver > H-J (larger
            //     dielectric loading; line electrically longer at 1 GHz;
            //     Im(Z) > 0 inductive).
            let all_capacitive = points.iter().all(|(_, z)| z.im < 0.0);
            let all_inductive = points.iter().all(|(_, z)| z.im > 0.0);
            if all_capacitive {
                eprintln!(
                    "  All Im(Z) < 0 across the sweep → resonance is ABOVE 1.3 GHz \
                     (line electrically too short at 1 GHz; ε_eff_solver < H-J)"
                );
            } else if all_inductive {
                eprintln!(
                    "  All Im(Z) > 0 across the sweep → resonance is BELOW 0.7 GHz \
                     (line electrically too long at 1 GHz; ε_eff_solver > H-J)"
                );
            }
            eprintln!("  Verdict: resonance offset NOT bracketed by 3-point sweep");
        }
    }
    (f_res, points)
}

/// Probe R3 — width sweep at 1 GHz (the original mom-002 reference
/// frequency, used here regardless of whether R2 located a resonance —
/// the brief notes "the trend is the signal"). Reports `|Z_in|` for
/// each width and flags `width-sensitive` if the spread exceeds a
/// factor of 2.
fn probe_r3_width_sweep() {
    eprintln!();
    eprintln!(
        "Probe R3 — Width sweep at f = {} GHz ({}x{} coarse mesh, L = {} mm):",
        F_HZ_NOMINAL * 1e-9,
        N_LENGTH,
        N_WIDTH,
        STRIP_L_M * 1e3,
    );
    eprintln!(
        "  {:>10}  {:>14}  {:>14}  {:>14}",
        "w (mm)", "Re(Z) (Ω)", "Im(Z) (Ω)", "|Z| (Ω)",
    );

    let widths_m: Vec<f64> = vec![2.0e-3, 2.94e-3, 4.0e-3];
    let mut z_mags: Vec<f64> = Vec::with_capacity(widths_m.len());
    for &w in &widths_m {
        let z = z_in_at(STRIP_L_M, w, N_LENGTH, N_WIDTH, F_HZ_NOMINAL);
        eprintln!(
            "  {:>10.3}  {:>14.4}  {:>14.4}  {:>14.4}",
            w * 1e3,
            z.re,
            z.im,
            z.norm(),
        );
        z_mags.push(z.norm());
    }

    let z_min = z_mags.iter().cloned().fold(f64::INFINITY, f64::min);
    let z_max = z_mags.iter().cloned().fold(0.0_f64, f64::max);
    let spread_ratio = z_max / z_min.max(1e-30);
    eprintln!(
        "  |Z_in| spread: min = {:.2} Ω, max = {:.2} Ω, ratio = {:.3}",
        z_min, z_max, spread_ratio,
    );
    // Width-sensitivity threshold: 2× spread in |Z_in| across a 2x
    // width range is significant; anything less is in the
    // "trend-noise" band.
    if spread_ratio > 2.0 {
        eprintln!("  Verdict: width-sensitive — |Z_in| varies > 2× across w ∈ [2, 4] mm");
    } else {
        eprintln!("  Verdict: width-insensitive — |Z_in| ratio < 2× across w ∈ [2, 4] mm");
    }
}

/// Run the three probes and print the consolidated verdict. Marked
/// `#[ignore]` so the suite never runs it by default; invoke
/// explicitly via:
///
/// ```text
/// cargo test -p yee-mom --release --test mom_002_13x_residual_diagnostic \
///     -- --ignored --nocapture
/// ```
///
/// Wall-time budget on the in-tree hardware: ~60-120 s total (3-freq
/// sweep + 3-width sweep on the coarse 40×8 mesh).
#[test]
#[ignore = "diagnostic: probes ε_eff bias / resonance offset / width sensitivity behind the 13× capacitive residual"]
fn mom_002_13x_residual_diagnostic() {
    eprintln!("--- Track MMMMMMM: mom-002 ~13x capacitive residual diagnostic ---");
    eprintln!(
        "Geometry (post-IIIIIII reframe): L = {} mm, w_nominal = {} mm, \
         centered port, uniform y-spacing. Coarse mesh ({}x{}) used for \
         the sweep loops to keep wall-time under ~2 minutes; the \
         production headline is 82x16.",
        STRIP_L_M * 1e3,
        STRIP_W_M * 1e3,
        N_LENGTH,
        N_WIDTH,
    );
    eprintln!(
        "Substrate: ε_r = {}, h = {} mm, Z_0 (H-J analytic) ≈ 51 Ω, \
         |Z_in| (production at 1 GHz) ≈ 674 Ω = 13.2 × Z_0",
        EPS_R,
        H_SUBSTRATE_M * 1e3,
    );

    let (eps_eff_solver, r1_biased) = probe_r1_eps_eff_from_greens();
    let (f_res_opt, _r2_points) = probe_r2_frequency_sweep();
    probe_r3_width_sweep();

    eprintln!();
    eprintln!("Overall verdict:");
    eprintln!(
        "  R1 (ε_eff biasing):       {}  (ε_eff_solver = {:.4} vs H-J {})",
        if r1_biased { "biased" } else { "matches" },
        eps_eff_solver,
        EPS_EFF_HJ_ANALYTIC,
    );
    match f_res_opt {
        Some(f) => eprintln!(
            "  R2 (resonance offset):    detected (f_res ≈ {:.4} GHz)",
            f * 1e-9
        ),
        None => eprintln!(
            "  R2 (resonance offset):    not bracketed (resonance outside [0.7, 1.3] GHz)"
        ),
    }
    eprintln!("  Dominant cause of capacitive Im(Z) = -674 Ω at 1 GHz:");
    // Decision table:
    //   R1 biased  + R2 offset detected  → R1 / R2 (both consistent
    //     with the same ε_eff drift; the offset is the integrated
    //     consequence of the bias)
    //   R1 not biased + R2 offset        → R2 alone (geometry-only)
    //   R1 biased + R2 not bracketed     → R1 (the bias is large
    //     enough that resonance is outside the sweep window)
    //   neither                          → leave as inconclusive /
    //     consider R3 finite-width or higher-order kernel issues
    let cause = match (r1_biased, f_res_opt.is_some()) {
        (true, true) => {
            "    R1 ε_eff biasing  +  R2 resonance offset (consistent: the offset is \
             the integrated consequence of the bias)"
        }
        (true, false) => {
            "    R1 ε_eff biasing dominates (resonance pushed outside the [0.7, 1.3] GHz \
             sweep window — bias is large)"
        }
        (false, true) => {
            "    R2 resonance offset alone (geometry-only; ε_eff is on-target but L is \
             slightly off from λ_eff / 2 at 1 GHz)"
        }
        (false, false) => {
            "    inconclusive — neither R1 nor R2 detected a clear bias; consider R3 \
             finite-width or higher-order kernel issues"
        }
    };
    eprintln!("{cause}");
    eprintln!(
        "  (Z_0 reference for context: {} Ω; production |Z_in| ≈ 674 Ω ≈ 13.2 × Z_0)",
        Z0_REF,
    );
}
