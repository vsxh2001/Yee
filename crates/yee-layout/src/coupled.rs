//! Static even/odd electrical model for symmetric edge-coupled microstrip.
//!
//! This module adds a closed-form even/odd-mode model for a *symmetric* pair of
//! edge-coupled microstrip lines (two equal-width strips, width `W`, separated
//! by a gap `S`, on a substrate of height `h` and permittivity `ε_r`). It
//! produces the even- and odd-mode characteristic impedances and effective
//! permittivities, plus a coupler-style coupling coefficient. It is the
//! validatable `k` reference the coupled-resonator FDTD driver (Filter Phase
//! F1.1b.1) needs and the initial-dimensioning model (F1.2) uses to pick a
//! gap/width from a target coupling before EM-in-the-loop refinement.
//!
//! Pure `f64`, WASM-safe, no FDTD, no new dependency (ADR-0089 / ADR-0094).
//!
//! ## Model
//!
//! **Kirschning-Jansen quasi-static (zero-frequency) coupled-microstrip
//! model.** Kirschning & Jansen, "Accurate Wide-Range Design Equations for the
//! Frequency-Dependent Characteristic of Parallel Coupled Microstrip Lines,"
//! *IEEE Trans. MTT*, vol. 32, no. 1, pp. 83–90, Jan. 1984 (corr. Nov. 1985).
//! The model builds on the single-line Hammerstad-Jensen `Z₀(u)` and effective
//! permittivity (Hammerstad & Jensen 1980), evaluated at a coupling-modified
//! width for the even mode and corrected by the `Q₁…Q₁₀` functions for the
//! impedances; the static quasi-TEM (DC) limit is taken — full
//! frequency-dispersion is out of scope for this static gate (spec §"Out of
//! scope"). Published accuracy is **better than ≈ 1.4 %** over `0.1 ≤ W/h ≤ 10`,
//! `0.1 ≤ S/h ≤ 10` versus rigorous numerical reference data.
//!
//! The equation set follows the canonical implementation in the QUCS circuit
//! simulator (`qucs-core/src/components/microstrip/mscoupled.cpp`,
//! `analysQuasiStatic`, "Kirschning" branch), which transcribes the
//! Kirschning-Jansen 1984 formulae. The single-line helpers
//! ([`crate::microstrip_width`]'s sibling forms) are reused: the even/odd
//! effective permittivities reuse the Hammerstad-Jensen single-line `εeff(u)`
//! form ([`hj_eps_eff`]) and the impedances reuse the Hammerstad-Jensen
//! single-line `Z₀(u)` form ([`hj_z0_air`]).
//!
//! ## Validation
//!
//! `coupled_microstrip` reproduces the worked example in Steer, *Microwave and
//! RF Design II: Transmission Lines* (3rd ed.), §5.6, Example 5.6.1 (alumina
//! `ε_r = 10`, `W = h = 500 µm` so `W/h = 1`, `S = 250 µm` so `S/h = 0.5`):
//! published `Z₀ₑ = 59 Ω`, `Z₀ₒ = 37 Ω`, `εeff,e = 7.28`, `εeff,o = 5.82`. The
//! model gives `Z₀ₑ ≈ 59.07 Ω`, `Z₀ₒ ≈ 36.96 Ω` (< 0.2 % error), gated by
//! `tests/coupled_001_vs_published.rs`.

use std::f64::consts::PI;

/// Free-space wave impedance `Z₀ = η₀ = 376.730…` Ω, the value the
/// Kirschning-Jansen / Hammerstad-Jensen forms are normalised against (the same
/// constant QUCS uses).
const Z_FREE_SPACE: f64 = 376.730_313_461_770_66;

/// Even/odd-mode electrical parameters of a symmetric edge-coupled microstrip
/// pair, from the static [Kirschning-Jansen] model.
///
/// The *even* mode is the two strips driven in phase (a magnetic wall in the
/// symmetry plane); the *odd* mode is the two strips driven out of phase (an
/// electric wall). Because the odd mode concentrates more field in the
/// substrate-adjacent gap region its effective permittivity is lower and its
/// characteristic impedance is lower than the even mode's, with the split
/// widening as the gap closes.
///
/// [Kirschning-Jansen]: crate::coupled
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupledMicrostrip {
    /// Even-mode characteristic impedance `Z₀ₑ`, ohms.
    pub z0e_ohm: f64,
    /// Odd-mode characteristic impedance `Z₀ₒ`, ohms.
    pub z0o_ohm: f64,
    /// Even-mode effective relative permittivity `εeff,e` (dimensionless).
    pub eps_eff_e: f64,
    /// Odd-mode effective relative permittivity `εeff,o` (dimensionless).
    pub eps_eff_o: f64,
}

/// Hammerstad-Jensen single-line `a(u)·b(εr)` exponent factors used by the
/// effective-permittivity form ([`hj_eps_eff`]).
///
/// ```text
/// a(u)  = 1 + (1/49)·ln[(u⁴ + (u/52)²)/(u⁴ + 0.432)] + (1/18.7)·ln(1 + (u/18.1)³)
/// b(εr) = 0.564·((εr − 0.9)/(εr + 3))^0.053
/// ```
fn hj_ab(u: f64, eps_r: f64) -> (f64, f64) {
    let u2 = u * u;
    let u4 = u2 * u2;
    let a = 1.0 + (u4 + (u / 52.0).powi(2)).ln() / 49.0 - (u4 + 0.432).ln() / 49.0
        + (1.0 + (u / 18.1).powi(3)).ln() / 18.7;
    let b = 0.564 * ((eps_r - 0.9) / (eps_r + 3.0)).powf(0.053);
    (a, b)
}

/// Hammerstad-Jensen single-line effective relative permittivity for a strip of
/// normalised width `u = W/h`:
///
/// ```text
/// εeff(u) = (εr + 1)/2 + (εr − 1)/2 · (1 + 10/u)^(−a(u)·b(εr))
/// ```
///
/// This is the high-accuracy Hammerstad-Jensen form (< 0.2 % over
/// `0.01 ≤ u ≤ 100`); the crate's [`crate::eps_eff`] is the simpler Pozar §3.8
/// form. The Kirschning-Jansen coupled model is built on this HJ form, so it is
/// used here for consistency with the published coefficients.
pub fn hj_eps_eff(u: f64, eps_r: f64) -> f64 {
    let (a, b) = hj_ab(u, eps_r);
    (eps_r + 1.0) / 2.0 + (eps_r - 1.0) / 2.0 * (1.0 + 10.0 / u).powf(-a * b)
}

/// Hammerstad-Jensen single-line characteristic impedance of a strip of
/// normalised width `u = W/h` **in air** (i.e. before dividing by `√εeff`):
///
/// ```text
/// f(u)  = 6 + (2π − 6)·exp(−(30.666/u)^0.7528)
/// Z₀ᵃⁱʳ = (Z_free/2π)·ln(f(u)/u + √(1 + (2/u)²))
/// ```
///
/// Accuracy < 0.01 % for `u ≥ 0.06` (Hammerstad-Jensen 1980). The line
/// impedance on the dielectric is `Z₀ᵃⁱʳ(u)/√εeff(u)`.
pub fn hj_z0_air(u: f64) -> f64 {
    let f_u = 6.0 + (2.0 * PI - 6.0) * (-(30.666 / u).powf(0.7528)).exp();
    Z_FREE_SPACE / (2.0 * PI) * (f_u / u + (1.0 + (2.0 / u).powi(2)).sqrt()).ln()
}

/// Static (zero-frequency) even/odd electrical parameters of a symmetric
/// edge-coupled microstrip pair, via the Kirschning-Jansen quasi-static model.
///
/// See the [module docs](crate::coupled) for the model, its citation, and its
/// published accuracy (≈ 1.4 % over `0.1 ≤ W/h ≤ 10`, `0.1 ≤ S/h ≤ 10`). The
/// returned [`CoupledMicrostrip`] always satisfies `z0e_ohm > z0o_ohm > 0` and
/// `eps_eff_e, eps_eff_o > 0` for physical inputs, with both impedances
/// approaching the single-line value as the gap grows.
///
/// # Arguments
///
/// - `w_m` — strip width `W` (each of the two equal strips), metres.
/// - `s_m` — edge-to-edge coupling gap `S` between the strips, metres.
/// - `h_m` — substrate height `h`, metres.
/// - `eps_r` — substrate relative permittivity `εr`.
///
/// # Panics
///
/// Does not panic for finite positive inputs; degenerate inputs (`w_m ≤ 0`,
/// `s_m ≤ 0`, `h_m ≤ 0`) yield non-finite or unphysical results and are the
/// caller's responsibility to avoid (this is a low-level model helper).
pub fn coupled_microstrip(w_m: f64, s_m: f64, h_m: f64, eps_r: f64) -> CoupledMicrostrip {
    // Normalised width and gap.
    let u = w_m / h_m;
    let g = s_m / h_m;
    let g2 = g * g;
    let g3 = g2 * g;
    let exp_neg_g = (-g).exp();

    // --- Even-mode effective permittivity ---------------------------------
    // Evaluate the single-line HJ εeff at a coupling-modified width `v`.
    let v = u * (20.0 + g2) / (10.0 + g2) + g * exp_neg_g;
    let eps_eff_e = hj_eps_eff(v, eps_r);

    // --- Odd-mode effective permittivity ----------------------------------
    // Built on the single-line HJ εeff at the unmodified width `u`.
    let eps_eff_single = hj_eps_eff(u, eps_r);
    let d = 0.593 + 0.694 * (-0.562 * u).exp();
    let bo = 0.747 * eps_r / (0.15 + eps_r);
    let co = bo - (bo - 0.207) * (-0.414 * u).exp();
    let ao = 0.7287 * (eps_eff_single - (eps_r + 1.0) / 2.0) * (1.0 - (-0.179 * u).exp());
    let eps_eff_o =
        ((eps_r + 1.0) / 2.0 + ao - eps_eff_single) * (-co * g.powf(d)).exp() + eps_eff_single;

    // --- Single-line characteristic impedance on the dielectric -----------
    let zl1 = hj_z0_air(u) / eps_eff_single.sqrt();

    // --- Even-mode characteristic impedance (Q1..Q4) ----------------------
    let q1 = 0.8695 * u.powf(0.194);
    let q2 = 1.0 + 0.7519 * g + 0.189 * g.powf(2.31);
    let q3 = 0.1975
        + (16.6 + (8.4 / g).powi(6)).powf(-0.387)
        + (g.powi(10) / (1.0 + (g / 3.4).powi(10))).ln() / 241.0;
    let q4 = q1 / q2 * 2.0 / (exp_neg_g * u.powf(q3) + (2.0 - exp_neg_g) * u.powf(-q3));
    let z0e_ohm = (eps_eff_single / eps_eff_e).sqrt() * zl1
        / (1.0 - zl1 * eps_eff_single.sqrt() * q4 / Z_FREE_SPACE);

    // --- Odd-mode characteristic impedance (Q5..Q10) ----------------------
    let q5 = 1.794 + 1.14 * (1.0 + 0.638 / (g + 0.517 * g.powf(2.43))).ln();
    let q6 = 0.2305
        + (g.powi(10) / (1.0 + (g / 5.8).powi(10))).ln() / 281.3
        + (1.0 + 0.598 * g.powf(1.154)).ln() / 5.1;
    let q7 = (10.0 + 190.0 * g2) / (1.0 + 82.3 * g3);
    let q8 = (-6.5 - 0.95 * g.ln() - (g / 0.15).powi(5)).exp();
    let q9 = q7.ln() * (q8 + 1.0 / 16.5);
    let q10 = (q2 * q4 - q5 * (u.ln() * q6 * u.powf(-q9)).exp()) / q2;
    let z0o_ohm = (eps_eff_single / eps_eff_o).sqrt() * zl1
        / (1.0 - zl1 * eps_eff_single.sqrt() * q10 / Z_FREE_SPACE);

    CoupledMicrostrip {
        z0e_ohm,
        z0o_ohm,
        eps_eff_e,
        eps_eff_o,
    }
}

/// Coupler-style coupling coefficient `k = (Z₀ₑ − Z₀ₒ)/(Z₀ₑ + Z₀ₒ)`.
///
/// This is the voltage-coupling factor of an edge-coupled-line section
/// (Pozar §7.6): `k → 0` for widely-spaced (weakly coupled) strips and grows
/// toward `1` as the gap closes. For physical inputs from
/// [`coupled_microstrip`] it lies in `(0, 1)` and **decreases monotonically as
/// the gap `S` increases**.
pub fn coupling_coefficient(m: &CoupledMicrostrip) -> f64 {
    (m.z0e_ohm - m.z0o_ohm) / (m.z0e_ohm + m.z0o_ohm)
}
