//! Sommerfeld surface-wave pole search — multi-frequency sanity table.
//!
//! Drives `yee_mom`'s Newton-Raphson pole search across the canonical
//! FR-4 / 1.6 mm substrate at 1, 2.4, and 5 GHz. The DoD #1 corridor in
//! the spec was authored assuming the *thick-slab / strip-wave* limit
//! `k_ρ ≈ k_0 √((ε_r+1)/2) ≈ 1.64 k_0`. Re-derivation against Pozar §3.7
//! shows that for a **thin** substrate (`k_0 h √(ε_r-1) ≪ 1`) the bare
//! grounded-slab TM₀ pole sits *much closer* to `k_0` than that quasi-
//! static value:
//!
//! ```text
//!   α_0  ≈  k_0² · h · (ε_r − 1) / ε_r,
//!   k_ρ  ≈  k_0 · √(1 + α_0² / k_0²).
//! ```
//!
//! For FR-4 / 1.6 mm: `k_ρ / k_0` is `≈ 1.0003` at 1 GHz, `≈ 1.0018` at
//! 2.4 GHz, `≈ 1.0079` at 5 GHz — climbing steadily as the slab becomes
//! electrically thicker. These are the bands this integration test
//! asserts; the higher band quoted in the original spec is preserved
//! only as a comment for future readers comparing notes.

use num_complex::Complex64;
use yee_mom::__internal::sommerfeld::{SwChannel, d_tm, newton_pole, residue, thin_slab_guess};

const EPS_R_FR4: f64 = 4.4;
const H_FR4: f64 = 1.6e-3;

fn k0_at(freq_hz: f64) -> f64 {
    std::f64::consts::TAU * freq_hz / yee_core::units::C0
}

/// Sanity table: Newton seeded at [`thin_slab_guess`] converges in ≤ 15
/// iterations across 1, 2.4, and 5 GHz; converged `|D| < 1e-9`;
/// `k_ρ / k_0` lies in the thin-slab band stated in the module
/// docstring. The DoD #1 escape hatch ("blocked > 25 min → surface and
/// stop") was not triggered — Newton finds the correct (physically-
/// real, not the spec's thick-slab estimate) TM₀ pole in single-digit
/// iterations for every frequency tested.
#[test]
fn pole_search_fr4_at_three_frequencies() {
    let cases: &[(f64, (f64, f64))] = &[
        (1.0e9, (1.0, 1.001)),
        (2.4e9, (1.0, 1.005)),
        (5.0e9, (1.0, 1.020)),
    ];

    for &(f, (lo, hi)) in cases {
        let k0 = k0_at(f);
        let seed = thin_slab_guess(EPS_R_FR4, H_FR4, k0);
        let (pole, iters) = newton_pole(SwChannel::Tm, seed, EPS_R_FR4, H_FR4, k0)
            .unwrap_or_else(|e| panic!("Newton failed at {f:.3e} Hz: {e:?}"));
        let resid = d_tm(pole, EPS_R_FR4, H_FR4, k0).norm();
        let ratio = pole.re / k0;
        assert!(
            resid < 1e-9,
            "f={f:.3e} Hz: |D| at pole = {resid:e} (tol 1e-9)"
        );
        assert!(
            (lo..hi).contains(&ratio),
            "f={f:.3e} Hz: k_ρ/k_0 = {ratio} outside thin-slab band [{lo}, {hi}]"
        );
        assert!(
            iters <= 15,
            "f={f:.3e} Hz: Newton took {iters} iters (budget 15)"
        );
        let _: Complex64 = pole;
    }
}

/// Diagnostic: print the FR-4 TM₀ pole locations at {1, 2.4, 5} GHz —
/// used during Phase 1.1.1.2 development to record exact pole positions
/// for the agent report. Marked `#[ignore]` so it never runs by default
/// (the table above is the validated gate).
#[test]
#[ignore = "diagnostic only: prints exact FR-4 TM₀ pole locations"]
fn pole_diagnostic_table() {
    for &f in &[1.0e9, 2.4e9, 5.0e9] {
        let k0 = k0_at(f);
        let seed = thin_slab_guess(EPS_R_FR4, H_FR4, k0);
        let (pole, iters) = newton_pole(SwChannel::Tm, seed, EPS_R_FR4, H_FR4, k0).expect("conv");
        let resid = d_tm(pole, EPS_R_FR4, H_FR4, k0).norm();
        let r = residue(SwChannel::Tm, pole, EPS_R_FR4, H_FR4, k0).expect("res");
        eprintln!(
            "f={:>4.2} GHz: k_rho/k_0={:.6} k_rho={:.4}+j{:.4} rad/m |D|={:.2e} iters={} residue={:.4e}+j{:.4e}",
            f * 1e-9,
            pole.re / k0,
            pole.re,
            pole.im,
            resid,
            iters,
            r.re,
            r.im,
        );
    }
}
