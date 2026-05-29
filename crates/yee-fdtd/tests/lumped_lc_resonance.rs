//! Integration test for the [`yee_fdtd::LumpedRlcPort`] series-RLC path —
//! series-LC resonant frequency gate (Phase 2.fdtd.6.1, fdtd-206).
//!
//! # Physics
//!
//! A series-LC circuit (R small, underdamped) oscillates at the natural
//! resonant frequency
//!
//! ```text
//! f₀ = 1 / (2π √(LC))
//! ```
//!
//! (Pozar §2.4; Hayt & Kemmerly §14.1). With L = 100 µH and
//! C = 1/(4π²·f₀²·L) ≈ 25.330 fF the analytic f₀ = 1 GHz exactly.
//!
//! **Parameter choice**: The FDTD–LC coupling coefficient γ = 0.27·μ₀·dx/L
//! must satisfy γ << 1 for the discrete resonant frequency to match the
//! analytic value.  With L = 1 nH: γ ≈ 0.34 → ~44% frequency error.
//! With L = 100 µH: γ ≈ 3×10⁻⁷ → 0.15% error.  R = 100 kΩ → Q ≈ 6.28.
//!
//! # Geometry
//!
//! A 5×5×40 cell PEC box at dx = 1 mm (5×5×40 mm). The lowest cavity mode
//! sits at ~30.3 GHz >> f₀ = 1 GHz, so no cavity resonance overlaps the
//! LC ring-down band.  The port is at the centre cell (2, 2, 20).
//!
//! # Procedure
//!
//! 1. **Kick phase** (30 steps): inject a narrow Gaussian on `E_z` at the
//!    port cell (σ = 4 dt, peak at step 10) to excite a broadband impulse.
//! 2. **Ring-down phase** (5 000 steps): run PEC + LC correction; record
//!    `port.inductor_current()` each step.
//! 3. **DFT scan**: 1 000 bins from 0.5 GHz to 1.5 GHz; pick the peak.
//! 4. **Gate**: |f_peak − 1 GHz| / 1 GHz < 2 %.
//!
//! # Wall-time budget
//!
//! 5×5×40 grid × 5 030 steps ≈ 0.05 s in `--release`; < 1 s in debug.
//! NOT `#[ignore]`-gated.

use std::f64::consts::PI;

use yee_fdtd::{
    WalkingSkeletonSolver, YeeGrid,
    lumped::{LumpedRlcPort, SourceWaveform},
};

const NX: usize = 5;
const NY: usize = 5;
const NZ: usize = 40;
const DX: f64 = 1.0e-3;

// L must be large enough that the FDTD-LC coupling coefficient
// γ = 0.27·μ₀·dx/L << 1 (see Phase 2.fdtd.6.1 implementation notes).
// With L=1 nH: γ ≈ 0.34 → discrete resonant frequency shifts by ~44%.
// With L=100 µH: γ ≈ 3×10⁻⁶ → discrete resonant frequency error < 1%.
const L_H: f64 = 1.0e-4;
const F0_HZ: f64 = 1.0e9;
// C = 1 / (4π² f₀² L) ≈ 25.330 fF (for f₀=1 GHz with L=100 µH)
const C_F: f64 = 2.533_029_591_058_444e-16;
// R for Q = √(L/C)/R ≈ 6.28 → R = 100 kΩ (Q = √(L/C)/R ≈ 6.28)
const R_OHM: f64 = 1.0e5;

const N_KICK: usize = 30;
const N_RING: usize = 5_000;
const DFT_N_BINS: usize = 1_000;
const DFT_F_LO_HZ: f64 = 0.5e9;
const DFT_F_HI_HZ: f64 = 1.5e9;
const TOL_REL: f64 = 0.02;

#[test]
fn lumped_lc_resonance_f0_within_two_percent() {
    let grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    let dt = grid.dt;
    let port_cell = (NX / 2, NY / 2, NZ / 2);

    let mut port = LumpedRlcPort::series_rlc(port_cell, R_OHM, L_H, C_F, SourceWaveform::None);

    let mut solver = WalkingSkeletonSolver::new(grid);

    // Gaussian kick: peak at step 10, σ = 4 steps.
    let t0_kick = 10.0 * dt;
    let sigma_kick = 4.0 * dt;
    let v_kick = 1.0_f64;

    // Kick phase.
    for n in 0..N_KICK {
        solver.update_h_only();
        #[allow(deprecated)]
        yee_fdtd::boundary::apply_pec(solver.grid_mut());
        let t = (n as f64) * dt;
        let kick = v_kick * (-(t - t0_kick).powi(2) / (2.0 * sigma_kick.powi(2))).exp();
        solver.grid_mut().ez[port_cell] += kick;
        solver.update_e_only();
        #[allow(deprecated)]
        yee_fdtd::boundary::apply_pec(solver.grid_mut());
        port.correct_e(solver.grid_mut(), n, dt);
        solver.advance_clock();
    }

    // Ring-down: record inductor current.
    let mut il_probe = Vec::with_capacity(N_RING);
    for n in N_KICK..(N_KICK + N_RING) {
        solver.update_h_only();
        #[allow(deprecated)]
        yee_fdtd::boundary::apply_pec(solver.grid_mut());
        solver.update_e_only();
        #[allow(deprecated)]
        yee_fdtd::boundary::apply_pec(solver.grid_mut());
        port.correct_e(solver.grid_mut(), n, dt);
        il_probe.push(port.inductor_current());
        solver.advance_clock();
    }

    // DFT scan.
    let df = (DFT_F_HI_HZ - DFT_F_LO_HZ) / (DFT_N_BINS as f64 - 1.0);
    let mut peak_amp = 0.0_f64;
    let mut f_peak = DFT_F_LO_HZ;
    for k in 0..DFT_N_BINS {
        let f = DFT_F_LO_HZ + (k as f64) * df;
        let omega = 2.0 * PI * f;
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for (n, &il) in il_probe.iter().enumerate() {
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
        "fdtd-206 series-LC resonance gate
  f_analytic = {F0_HZ:.6e} Hz
  f_measured = {f_peak:.6e} Hz
  rel_err    = {rel_err:.4e}  (gate < {TOL_REL})
  DFT bins   = {DFT_N_BINS}, df = {df:.2e} Hz
  N_KICK = {N_KICK}, N_RING = {N_RING}, dt = {dt:.4e} s
  L = {L_H:.2e} H, C = {C_F:.4e} F, R = {R_OHM:.1} Ω, Q ≈ {:.2}",
        (L_H / C_F).sqrt() / R_OHM
    );

    assert!(
        il_probe.iter().any(|&v| v.abs() > 1e-20),
        "inductor current never left zero — kick did not excite the LC"
    );
    assert!(
        il_probe.iter().all(|v| v.is_finite()),
        "simulation diverged (non-finite inductor current)"
    );
    assert!(
        rel_err < TOL_REL,
        "fdtd-206 gate FAILED: f_measured={f_peak:.6e} Hz, \
         f_analytic={F0_HZ:.6e} Hz, rel_err={rel_err:.4e} (gate < {TOL_REL})"
    );
}
