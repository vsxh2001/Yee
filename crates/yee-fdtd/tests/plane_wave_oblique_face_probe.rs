//! Per-face leakage probe for the oblique TF/SF source (Phase 2.fdtd.5.3.2).
//!
//! Phase 2.fdtd.5.3 / 5.3.1 plateaued at 15.6× contrast for the 30°/45°
//! E_phi case — well below the 1000× DoD. GGGGG's finding: dispersion
//! mismatch is NOT the dominant SF leak; there is a sign / cross-section
//! bug in the 12-face oblique stencils that was never exercised at
//! normal incidence.
//!
//! This probe measures leakage in each of the six SF half-spaces
//! independently (lo-x, hi-x, lo-y, hi-y, lo-z, hi-z) so individual
//! faces' stencil bugs can be isolated. Each face's leakage is
//! reported as a fraction of the inside-TF amplitude (so 0.01 = 1%).
//!
//! Pattern mirrors `plane_wave_oblique.rs`.
//!
//! The probe is `#[ignore]`d to keep CI fast; run with
//! `cargo test -p yee-fdtd --release oblique_face_leakage_probe -- --include-ignored --nocapture`
//! to inspect the per-face values directly.

use std::ops::Range;

use yee_fdtd::{CpmlParams, PlaneWaveSource, WalkingSkeletonSolver, YeeGrid};

fn max_abs_field_in_region(
    grid: &YeeGrid,
    is: Range<usize>,
    js: Range<usize>,
    ks: Range<usize>,
) -> f64 {
    let mut peak: f64 = 0.0;
    for i in is.clone() {
        for j in js.clone() {
            for k in ks.clone() {
                if i < grid.ex.shape()[0] && j < grid.ex.shape()[1] && k < grid.ex.shape()[2] {
                    let v = grid.ex[(i, j, k)].abs();
                    if v > peak {
                        peak = v;
                    }
                }
                if i < grid.ey.shape()[0] && j < grid.ey.shape()[1] && k < grid.ey.shape()[2] {
                    let v = grid.ey[(i, j, k)].abs();
                    if v > peak {
                        peak = v;
                    }
                }
                if i < grid.ez.shape()[0] && j < grid.ez.shape()[1] && k < grid.ez.shape()[2] {
                    let v = grid.ez[(i, j, k)].abs();
                    if v > peak {
                        peak = v;
                    }
                }
            }
        }
    }
    peak
}

/// Run a 60³ TF box with oblique 30°/45° E_phi polarization for 400 steps
/// and return `(inside_amp, leak_lo_x, leak_hi_x, leak_lo_y, leak_hi_y,
/// leak_lo_z, leak_hi_z)`. Each leak* value is the peak |E| in the
/// corresponding SF half-space.
fn run_oblique_face_probe(dispersion_match: bool) -> (f64, [f64; 6]) {
    run_oblique_face_probe_steps(dispersion_match, 400)
}

fn run_oblique_face_probe_steps(dispersion_match: bool, n_steps: usize) -> (f64, [f64; 6]) {
    use std::f64::consts::PI;

    const N: usize = 60;
    const DX: f64 = 5.0e-3;
    const FREQ_HZ: f64 = 3.0e9;
    const RAMP: usize = 40;
    const PAD: usize = 8;
    const NPML: usize = 8;

    const I0: usize = 20;
    const I1: usize = 40;
    const J0: usize = 20;
    const J1: usize = 40;
    const K0: usize = 20;
    const K1: usize = 40;

    let grid = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid.dt;
    let cpml_params = CpmlParams::for_grid(&grid, NPML);
    let mut solver = WalkingSkeletonSolver::with_cpml(grid, cpml_params);

    let theta = 30.0_f64.to_radians();
    let phi = 45.0_f64.to_radians();
    let psi = PI / 2.0; // E along e_phi
    let mut pw = PlaneWaveSource::with_oblique_incidence_match(
        I0,
        I1,
        J0,
        J1,
        K0,
        K1,
        theta,
        phi,
        psi,
        FREQ_HZ,
        RAMP,
        DX,
        dt,
        PAD,
        dispersion_match,
    );

    for _ in 0..n_steps {
        solver.step_with_plane_wave(&mut pw);
    }

    let inside_amp = max_abs_field_in_region(
        solver.grid(),
        (I0 + 4)..(I1 - 3),
        (J0 + 4)..(J1 - 3),
        (K0 + 4)..(K1 - 3),
    );

    const CPML_INTERIOR: usize = NPML + 2;
    let lo_x = max_abs_field_in_region(
        solver.grid(),
        CPML_INTERIOR..(I0 - 1),
        J0..(J1 + 1),
        K0..(K1 + 1),
    );
    let hi_x = max_abs_field_in_region(
        solver.grid(),
        (I1 + 2)..(N - CPML_INTERIOR),
        J0..(J1 + 1),
        K0..(K1 + 1),
    );
    let lo_y = max_abs_field_in_region(
        solver.grid(),
        I0..(I1 + 1),
        CPML_INTERIOR..(J0 - 1),
        K0..(K1 + 1),
    );
    let hi_y = max_abs_field_in_region(
        solver.grid(),
        I0..(I1 + 1),
        (J1 + 2)..(N - CPML_INTERIOR),
        K0..(K1 + 1),
    );
    let lo_z = max_abs_field_in_region(
        solver.grid(),
        I0..(I1 + 1),
        J0..(J1 + 1),
        CPML_INTERIOR..(K0 - 1),
    );
    let hi_z = max_abs_field_in_region(
        solver.grid(),
        I0..(I1 + 1),
        J0..(J1 + 1),
        (K1 + 2)..(N - CPML_INTERIOR),
    );

    (inside_amp, [lo_x, hi_x, lo_y, hi_y, lo_z, hi_z])
}

#[test]
#[ignore = "diagnostic: dispersion-match enabled, oblique 30°/45° per-face leakage probe"]
fn oblique_face_leakage_probe_match() {
    let (inside_amp, leaks) = run_oblique_face_probe(true);
    let names = ["lo-x", "hi-x", "lo-y", "hi-y", "lo-z", "hi-z"];
    eprintln!(
        "oblique-30°/45° (dispersion match=true) inside amp = {inside_amp:.6e}"
    );
    for (name, leak) in names.iter().zip(leaks.iter()) {
        let frac = leak / inside_amp.max(1e-30);
        eprintln!(
            "  face {name}: peak |E| = {leak:.6e}  ({:.3}% of TF amplitude)",
            frac * 100.0
        );
    }
    let worst = leaks.iter().cloned().fold(0.0_f64, f64::max);
    let contrast = inside_amp / worst.max(1e-30);
    eprintln!(
        "  worst-face contrast = {contrast:.6e} ({:.2} dB)",
        20.0 * contrast.log10().max(-1000.0)
    );
    // No hard assertion here — diagnostic only. The 30/45 ephi gate
    // test in `plane_wave_oblique.rs` carries the DoD.
    assert!(inside_amp > 0.1);
}

#[test]
#[ignore = "diagnostic: dispersion-match enabled, sweep N_STEPS"]
fn oblique_face_leakage_probe_match_sweep_steps() {
    for &n in &[300usize, 400, 500, 600, 800] {
        let (inside_amp, leaks) = run_oblique_face_probe_steps(true, n);
        let worst = leaks.iter().cloned().fold(0.0_f64, f64::max);
        let lo_only = leaks[0].max(leaks[2]).max(leaks[4]);
        let contrast_worst = inside_amp / worst.max(1e-30);
        let contrast_lo = inside_amp / lo_only.max(1e-30);
        eprintln!(
            "n_steps={n:>4} inside={inside_amp:.6e} worst-face contrast={contrast_worst:.2} lo-only contrast={contrast_lo:.2}"
        );
    }
}

#[test]
#[ignore = "diagnostic: dispersion-match disabled (5.3 baseline) per-face leakage probe"]
fn oblique_face_leakage_probe_no_match() {
    let (inside_amp, leaks) = run_oblique_face_probe(false);
    let names = ["lo-x", "hi-x", "lo-y", "hi-y", "lo-z", "hi-z"];
    eprintln!(
        "oblique-30°/45° (dispersion match=false) inside amp = {inside_amp:.6e}"
    );
    for (name, leak) in names.iter().zip(leaks.iter()) {
        let frac = leak / inside_amp.max(1e-30);
        eprintln!(
            "  face {name}: peak |E| = {leak:.6e}  ({:.3}% of TF amplitude)",
            frac * 100.0
        );
    }
    let worst = leaks.iter().cloned().fold(0.0_f64, f64::max);
    let contrast = inside_amp / worst.max(1e-30);
    eprintln!(
        "  worst-face contrast = {contrast:.6e} ({:.2} dB)",
        20.0 * contrast.log10().max(-1000.0)
    );
    assert!(inside_amp > 0.1);
}
