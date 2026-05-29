//! fdtd-206: Lumped series-LC resonant frequency gate.
//!
//! Validates that the FDTD series-RLC ODE integration reproduces the
//! natural resonant frequency of a lumped LC circuit:
//!
//! ```text
//! f₀ = 1 / (2π √(LC))
//! ```
//!
//! # Geometry
//!
//! 5×5×40 PEC box at dx=1 mm. Series-LC port at cell (2,2,20), z-axis.
//! L = 1 nH, C ≈ 25.33 pF → f₀ = 1 GHz analytic. R = 1 Ω → Q ≈ 6.28.
//!
//! # Protocol
//!
//! Drive phase (6500 steps): inject a GaussianPulse at f₀=1 GHz via the port
//! voltage source. The pulse has 500 MHz bandwidth and peaks at step 200, so
//! by step ~1500 the source has decayed to negligible amplitude.
//! DFT scan on I_L over the late ring-down window [1500, 6500):
//! 1000 bins from 0.5 GHz to 1.5 GHz; find peak frequency.
//!
//! # Gate
//!
//! |f_peak − f₀| / f₀ < 2 %
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-fdtd --test lumped_lc_resonance -- --nocapture
//! ```

use std::f64::consts::PI;

use yee_fdtd::{
    WalkingSkeletonSolver, YeeGrid, boundary,
    lumped::{LumpedRlcPort, SourceWaveform},
};

// Grid
const NX: usize = 5;
const NY: usize = 5;
const NZ: usize = 40;
const DX: f64 = 1.0e-3;

// Circuit
const L_H: f64 = 1.0e-9;
const F0_HZ: f64 = 1.0e9;
const C_F: f64 = 1.0 / (4.0 * PI * PI * F0_HZ * F0_HZ * L_H);
const R_OHM: f64 = 1.0;

// Protocol — GaussianPulse source
const T0_STEPS: usize = 200;
const BW_HZ: f64 = 500.0e6;
const V0: f64 = 1.0e-3;
const N_TOTAL: usize = 6_500;
const DFT_START: usize = 1_500;
const DFT_N_BINS: usize = 1_000;
const DFT_F_LO: f64 = 0.5e9;
const DFT_F_HI: f64 = 1.5e9;
const TOL_REL: f64 = 0.02;

#[test]
fn lumped_lc_resonance_f0_within_two_percent() {
    let grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    let dt = grid.dt;

    let port_cell = (NX / 2, NY / 2, NZ / 2);

    // GaussianPulse excites the LC at resonance; source decays by step ~1500.
    let source = SourceWaveform::GaussianPulse {
        v0: V0,
        f0: F0_HZ,
        bw: BW_HZ,
        t0_steps: T0_STEPS,
    };

    let mut port = LumpedRlcPort::series_rlc(port_cell, R_OHM, L_H, C_F, source);
    let mut solver = WalkingSkeletonSolver::new(grid);

    let mut il_all = Vec::with_capacity(N_TOTAL);
    for n in 0..N_TOTAL {
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        solver.update_e_only();
        port.correct_e(solver.grid_mut(), n, dt);
        il_all.push(port.inductor_current());
        solver.advance_clock();
    }

    eprintln!(
        "Peak I_L: {:.4e} A at step ~{}",
        il_all.iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs())),
        il_all
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0)
    );

    // DFT scan over the late ring-down window (source off after ~step 1500)
    let late = &il_all[DFT_START..];
    let df = (DFT_F_HI - DFT_F_LO) / (DFT_N_BINS as f64 - 1.0);
    let mut peak_amp = 0.0_f64;
    let mut f_peak = DFT_F_LO;
    for k in 0..DFT_N_BINS {
        let f = DFT_F_LO + (k as f64) * df;
        let omega = 2.0 * PI * f;
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for (n, &il) in late.iter().enumerate() {
            let phase = omega * (n as f64) * dt;
            re += il * phase.cos();
            im += il * phase.sin();
        }
        let amp = (re * re + im * im).sqrt();
        if amp > peak_amp {
            peak_amp = amp;
            f_peak = f;
        }
    }

    let rel_err = (f_peak - F0_HZ).abs() / F0_HZ;

    eprintln!(
        "\nfdtd-206: series-LC resonant frequency gate
  f\u{2080}_analytic = {:.4} GHz   (1/(2\u{03c0}\u{221a}LC))
  f_measured  = {:.4} GHz
  rel_err     = {:.3} %     (gate < {} %)
  peak_amp    = {peak_amp:.4e}
",
        F0_HZ * 1e-9,
        f_peak * 1e-9,
        rel_err * 100.0,
        TOL_REL * 100.0,
    );

    assert!(
        peak_amp > 0.0,
        "DFT peak amplitude is zero — LC oscillation did not excite; \
         check GaussianPulse source parameters"
    );
    assert!(
        rel_err < TOL_REL,
        "fdtd-206 FAILED: f_measured={f_peak:.4e} Hz (rel_err={:.2}% > gate {}%)",
        rel_err * 100.0,
        TOL_REL * 100.0,
    );
}
