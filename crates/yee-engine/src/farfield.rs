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
}
