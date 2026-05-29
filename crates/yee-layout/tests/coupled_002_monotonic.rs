//! `coupled-002` — coupling-vs-gap monotonicity sanity (pure physics, no data).
//!
//! Textbook coupling-vs-gap law (Pozar §7.6, Hong & Lancaster ch. 5): for a
//! symmetric edge-coupled microstrip pair at fixed width / height / `ε_r`, the
//! even and odd modes always satisfy `Z₀ₑ > Z₀ₒ > 0` and have positive effective
//! permittivities, and the coupler coupling coefficient `k = (Z₀ₑ−Z₀ₒ)/(Z₀ₑ+Z₀ₒ)`
//! is strictly positive and **decreases as the gap widens** (the strips couple
//! more weakly the farther apart they are). This is a property of the physics,
//! not a published data point, so no external reference is cited.

use yee_layout::{coupled_microstrip, coupling_coefficient};

/// FR-4-like substrate; fixed strip width and height, swept gap.
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;

/// Strictly increasing gaps (metres). `k` must strictly decrease across these.
const GAPS_M: [f64; 5] = [0.2e-3, 0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3];

#[test]
fn each_gap_is_physical() {
    for &s in &GAPS_M {
        let m = coupled_microstrip(W_M, s, H_M, EPS_R);
        assert!(
            m.z0o_ohm > 0.0,
            "s = {s:.4e} m: Z0o = {:.3} Ω must be > 0",
            m.z0o_ohm
        );
        assert!(
            m.z0e_ohm > m.z0o_ohm,
            "s = {s:.4e} m: Z0e ({:.3}) must exceed Z0o ({:.3})",
            m.z0e_ohm,
            m.z0o_ohm
        );
        assert!(
            m.eps_eff_e > 0.0 && m.eps_eff_o > 0.0,
            "s = {s:.4e} m: eps_eff_e ({:.3}) and eps_eff_o ({:.3}) must be > 0",
            m.eps_eff_e,
            m.eps_eff_o
        );
    }
}

#[test]
fn coupling_strictly_decreases_with_gap() {
    let ks: Vec<f64> = GAPS_M
        .iter()
        .map(|&s| coupling_coefficient(&coupled_microstrip(W_M, s, H_M, EPS_R)))
        .collect();

    for (i, &k) in ks.iter().enumerate() {
        assert!(
            k > 0.0 && k < 1.0,
            "gap #{i} (s = {:.4e} m): k = {k:.4} must be in (0, 1)",
            GAPS_M[i]
        );
    }
    for w in ks.windows(2) {
        assert!(
            w[0] > w[1],
            "coupling must strictly decrease as gap grows, but k went {:.4} -> {:.4}",
            w[0],
            w[1]
        );
    }
}

/// As the gap grows large the two strips decouple: both even- and odd-mode
/// impedances approach the single (uncoupled) line value and `k → 0`. A
/// directional check that the widest gap is the most weakly coupled.
#[test]
fn widest_gap_is_weakest_coupling() {
    let k_tight = coupling_coefficient(&coupled_microstrip(W_M, GAPS_M[0], H_M, EPS_R));
    let k_wide = coupling_coefficient(&coupled_microstrip(
        W_M,
        *GAPS_M.last().unwrap(),
        H_M,
        EPS_R,
    ));
    assert!(
        k_wide < k_tight,
        "widest gap k ({k_wide:.4}) should be weaker than tightest gap k ({k_tight:.4})"
    );
}
