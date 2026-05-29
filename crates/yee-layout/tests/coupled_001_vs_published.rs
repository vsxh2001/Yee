//! `coupled-001` — published coupled-microstrip even/odd-mode gate.
//!
//! Ground-truth published worked example (do NOT loosen / do NOT invent):
//!
//! > M. Steer, *Microwave and RF Design II: Transmission Lines*, 3rd ed.,
//! > NC State University, 2019, §5.6, **Example 5.6.1** ("Even- and Odd-Mode
//! > Parameters").
//! > <https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electronics/Microwave_and_RF_Design_II_-_Transmission_Lines_(Steer)/05:_Coupled_Lines_and_Applications/5.06:_Formulas_for_Impedance_of_Coupled_Microstrip_Lines>
//! >
//! > Alumina substrate: `ε_r = 10`, `h = 500 µm`; strips `W = 500 µm`
//! > (`W/h = 1.0`), gap `s = 250 µm` (`s/h = 0.5`).
//! >
//! > Published results:
//! >   - `Z₀ₑ = 59 Ω`,   `Z₀ₒ = 37 Ω`
//! >   - `εeff,e = 7.28`, `εeff,o = 5.82`
//!
//! The implemented model is the Kirschning-Jansen quasi-static coupled-
//! microstrip model (Kirschning & Jansen, *IEEE Trans. MTT* 32(1):83–90, 1984;
//! published accuracy ≈ 1.4 %). It reproduces the reference to < 0.5 %, so a
//! 5 % gate is comfortable headroom; this is the model's stated tolerance, not
//! a fudge factor. If a check misses by > 5 % the Q-function transcription or
//! the single-line HJ helper is wrong — fix the math, do NOT loosen the gate.

use yee_layout::{coupled_microstrip, coupling_coefficient};

/// Gate tolerance — well inside the Kirschning-Jansen ≈ 1.4 % model accuracy.
const REL_TOL: f64 = 0.05;

/// Steer Example 5.6.1 geometry.
const EPS_R: f64 = 10.0;
const H_M: f64 = 500.0e-6;
const W_M: f64 = 500.0e-6; // W/h = 1.0
const S_M: f64 = 250.0e-6; // s/h = 0.5

/// Published reference values (Steer Example 5.6.1).
const Z0E_REF: f64 = 59.0;
const Z0O_REF: f64 = 37.0;
const EPS_EFF_E_REF: f64 = 7.28;
const EPS_EFF_O_REF: f64 = 5.82;

fn assert_rel(got: f64, want: f64, label: &str) {
    let rel = (got - want).abs() / want;
    assert!(
        rel < REL_TOL,
        "{label}: got {got:.4}, published {want:.4} (rel err {:.2}% > {:.0}%)",
        rel * 100.0,
        REL_TOL * 100.0
    );
}

#[test]
fn z0e_matches_steer_example_5_6_1() {
    let m = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    assert_rel(m.z0e_ohm, Z0E_REF, "Z0e");
}

#[test]
fn z0o_matches_steer_example_5_6_1() {
    let m = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    assert_rel(m.z0o_ohm, Z0O_REF, "Z0o");
}

#[test]
fn eps_eff_e_matches_steer_example_5_6_1() {
    let m = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    assert_rel(m.eps_eff_e, EPS_EFF_E_REF, "eps_eff_e");
}

#[test]
fn eps_eff_o_matches_steer_example_5_6_1() {
    let m = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    assert_rel(m.eps_eff_o, EPS_EFF_O_REF, "eps_eff_o");
}

/// Sanity on the derived coupling factor: with the published `Z₀ₑ = 59 Ω`,
/// `Z₀ₒ = 37 Ω` the coupler `k = (59 − 37)/(59 + 37) ≈ 0.229`. The model's
/// `coupling_coefficient` must land within the same 5 % of that published-derived
/// value.
#[test]
fn coupling_coefficient_matches_published_derived() {
    let m = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    let k_ref = (Z0E_REF - Z0O_REF) / (Z0E_REF + Z0O_REF);
    assert_rel(coupling_coefficient(&m), k_ref, "k");
}
