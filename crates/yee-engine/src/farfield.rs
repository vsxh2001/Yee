//! Absolute far-field products (FS.2b, ADR-0207): gain in dBi from the
//! NTFF pattern amplitude plus the FS.2a port records.
//!
//! ## The normalization chain, audited end to end
//!
//! - `yee_fdtd::NtffState::far_field_at` returns the pattern amplitude
//!   `F(ω, θ, φ)` — Taflove eq. 8.35 with the `e^{−jkr}/r` envelope
//!   dropped — accumulated as a **continuous-transform** DFT (kernel
//!   × dt). The physical field at distance r is `E = F·e^{−jkr}/r`, so
//!   the radiated **energy** per solid angle per rad/s of a pulsed run
//!   is `(1/π)·|F(ω)|²/η₀` (Parseval on the real time signal).
//! - The FS.2a port records give the circuit-side accepted energy. Its
//!   per-frequency density is `(1/π)·(Re[−V_src·I*] − R·|I|²)` with
//!   `V_src, I` the same continuous-transform DFTs — the frequency-domain
//!   twin of the gate-validated `E_emf − E_R` identity (engine-power-001
//!   closure 0.9917), and immune to the β-term pitfall the aperture-side
//!   `v_term·i` accounting measured (ADR-0207).
//! - In the gain ratio the `1/π` and every DFT scale factor cancel:
//!
//!   ```text
//!   G(ω, θ, φ) = 4π · |F|² / (η₀ · (Re[−V_src·I*] − R·|I|²))
//!   ```
//!
//! [`crate::sparams::single_bin_dft`] is a plain sum (no dt), so this
//! module multiplies the port DFTs by dt to match the NTFF convention.

/// Free-space wave impedance η₀ (Ω).
const ETA0: f64 = 376.730_313_668;

/// Per-frequency **accepted power density** from a port's FS.2a record
/// DFTs: `Re[−V_src·I*] − R·|I|²` with `V_src = v_src_dft·dt` and
/// `I = i_dft·dt` (inputs are raw [`crate::sparams::single_bin_dft`]
/// bins). Positive inside the drive band for a port that feeds the
/// field; a non-positive value means the frequency is outside the band
/// or the port is a passive load — a caller error for gain purposes.
pub fn accepted_power_density(
    v_src_dft: (f64, f64),
    i_dft: (f64, f64),
    r_ohm: f64,
    dt_s: f64,
) -> f64 {
    let (vr, vi) = (v_src_dft.0 * dt_s, v_src_dft.1 * dt_s);
    let (ir, ii) = (i_dft.0 * dt_s, i_dft.1 * dt_s);
    // Re[−V·I*] = −(vr·ir + vi·ii); |I|² = ir² + ii².
    -(vr * ir + vi * ii) - r_ohm * (ir * ir + ii * ii)
}

/// Gain in dBi at one direction: `10·log₁₀(4π·|F|²/(η₀·p_acc))` with
/// `e_far_mag = |F|` straight from [`crate::JobResult::far_field`] and
/// `p_acc` from [`accepted_power_density`] at the same frequency.
///
/// # Panics
///
/// Panics if `p_acc ≤ 0` (frequency outside the drive band, or the
/// record does not belong to the driven port).
pub fn gain_dbi(e_far_mag: f64, p_acc: f64) -> f64 {
    assert!(
        p_acc > 0.0 && p_acc.is_finite(),
        "gain_dbi: accepted power density must be positive (got {p_acc:.3e}) — \
         is the frequency inside the drive band and the record the driven port's?"
    );
    let g = 4.0 * std::f64::consts::PI * e_far_mag * e_far_mag / (ETA0 * p_acc);
    10.0 * g.log10()
}

/// A midpoint-rule (θ, φ) sphere raster for [`radiation_efficiency`]:
/// `θ_i = (i+0.5)·π/n_theta`, `φ_j = j·2π/n_phi`, row-major θ-outer —
/// pass straight to `NtffSpec::directions`.
pub fn sphere_grid(n_theta: usize, n_phi: usize) -> Vec<(f64, f64)> {
    let mut dirs = Vec::with_capacity(n_theta * n_phi);
    for i in 0..n_theta {
        let theta = (i as f64 + 0.5) * std::f64::consts::PI / n_theta as f64;
        for j in 0..n_phi {
            let phi = j as f64 * std::f64::consts::TAU / n_phi as f64;
            dirs.push((theta, phi));
        }
    }
    dirs
}

/// Radiation efficiency `η = ∮|F|²dΩ / (η₀·p_acc)` (FS.2c): the gain
/// theorem `∮G dΩ = 4π·η` evaluated by midpoint quadrature over a
/// [`sphere_grid`]-ordered |F| raster. `1` for a lossless antenna (up to
/// the NTFF scale, quadrature, and absorber leakage); drops with
/// substrate/conductor loss.
///
/// # Panics
///
/// Panics if `e_far` is not `n_theta·n_phi` long or `p_acc ≤ 0`.
pub fn radiation_efficiency(e_far: &[f64], n_theta: usize, n_phi: usize, p_acc: f64) -> f64 {
    assert_eq!(e_far.len(), n_theta * n_phi, "raster shape mismatch");
    assert!(p_acc > 0.0 && p_acc.is_finite(), "p_acc must be positive");
    let d_theta = std::f64::consts::PI / n_theta as f64;
    let d_phi = std::f64::consts::TAU / n_phi as f64;
    let mut integral = 0.0;
    for i in 0..n_theta {
        let theta = (i as f64 + 0.5) * d_theta;
        let w = theta.sin() * d_theta * d_phi;
        for j in 0..n_phi {
            let f = e_far[i * n_phi + j];
            integral += f * f * w;
        }
    }
    integral / (ETA0 * p_acc)
}

/// Render a full-sphere pattern as CSV (`theta_deg,phi_deg,e_far,gain_dbi`
/// header + one row per direction) — the FS.2c export artifact. Stable
/// formatting (`{:.6e}`) so exports are byte-checkable.
pub fn pattern_csv(directions: &[(f64, f64)], e_far: &[f64], p_acc: f64) -> String {
    let mut out = String::from("theta_deg,phi_deg,e_far,gain_dbi\n");
    for (&(theta, phi), &f) in directions.iter().zip(e_far) {
        out.push_str(&format!(
            "{:.3},{:.3},{:.6e},{:.3}\n",
            theta.to_degrees(),
            phi.to_degrees(),
            f,
            gain_dbi(f, p_acc),
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_dbi_reference_case_is_exact() {
        // Hand-built isotropic reference: choose V_src, I, R, then set
        // |F| so G = 1 exactly; dt = 1 keeps the DFT scaling trivial.
        let v_src = (-2.0, 0.0);
        let i = (0.01, 0.0);
        let r = 50.0;
        let p = accepted_power_density(v_src, i, r, 1.0);
        // Re[−V·I*] = 0.02; R|I|² = 0.005 → p = 0.015.
        assert!((p - 0.015).abs() < 1e-15);
        let f_mag = (ETA0 * p / (4.0 * std::f64::consts::PI)).sqrt();
        let g = gain_dbi(f_mag, p);
        assert!(g.abs() < 1e-9, "expected 0 dBi, got {g}");
    }

    #[test]
    fn dt_scaling_cancels_in_the_ratio_but_not_in_p() {
        // p scales as dt² (both DFTs), so raw bins from different-dt runs
        // are not comparable — the docs say multiply by dt; verify.
        let p1 = accepted_power_density((-2.0, 0.0), (0.01, 0.0), 50.0, 1.0);
        let p2 = accepted_power_density((-2.0, 0.0), (0.01, 0.0), 50.0, 2.0);
        assert!((p2 / p1 - 4.0).abs() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "accepted power density must be positive")]
    fn passive_record_is_rejected() {
        gain_dbi(1.0, 0.0);
    }

    #[test]
    fn isotropic_radiator_has_unit_efficiency() {
        // G = 1 everywhere → η = 1; the midpoint rule integrates sin θ
        // exactly enough that 24×24 is well under 0.1 %.
        let (nt, np) = (24, 24);
        let p = 0.015;
        let f_mag = (ETA0 * p / (4.0 * std::f64::consts::PI)).sqrt();
        let e: Vec<f64> = vec![f_mag; nt * np];
        let eta = radiation_efficiency(&e, nt, np, p);
        assert!((eta - 1.0).abs() < 1e-3, "η = {eta}");
    }

    #[test]
    fn sphere_grid_shape_and_csv_are_stable() {
        let dirs = sphere_grid(3, 4);
        assert_eq!(dirs.len(), 12);
        assert!((dirs[0].0.to_degrees() - 30.0).abs() < 1e-9);
        let e = vec![1.0e-11; 12];
        let csv = pattern_csv(&dirs, &e, 1.0e-23);
        assert!(csv.starts_with("theta_deg,phi_deg,e_far,gain_dbi\n"));
        assert_eq!(csv.lines().count(), 13);
        // Byte-stability: same inputs, same bytes.
        assert_eq!(csv, pattern_csv(&dirs, &e, 1.0e-23));
    }
}
