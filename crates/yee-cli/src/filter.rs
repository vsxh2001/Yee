//! `yee filter synth` handler (Filter Phase F0).
//!
//! Parses a [`yee_filter::FilterSpec`] TOML, synthesizes the lowpass prototype
//! and all-pole coupling matrix, sweeps the closed-form ideal response, writes
//! the S-parameters as a Touchstone `.s2p` via `yee-io`, and prints the
//! spec-mask verdict. Exit 0 on PASS, 1 on FAIL.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use num_complex::Complex64;
use std::process::ExitCode;

use yee_filter::{FilterSpec, check_mask, ideal_response, synthesize};
use yee_io::touchstone::{File, Format, FreqUnit};

/// Number of points in the response sweep written to the Touchstone file.
const SWEEP_POINTS: usize = 401;
/// Sweep span as a multiple of the fractional bandwidth on each side of `f0`.
/// `f0·(1 ± SPAN·FBW/2)` keeps a wide skirt around the passband.
const SPAN_MULT: f64 = 6.0;

/// Run `yee filter synth <spec> [--output <out.s2p>]`.
pub fn run_synth(spec_path: &Path, output: Option<&Path>) -> Result<ExitCode> {
    let text = std::fs::read_to_string(spec_path)
        .with_context(|| format!("failed to read filter spec {}", spec_path.display()))?;
    let spec: FilterSpec = toml::from_str(&text)
        .with_context(|| format!("failed to parse filter spec {}", spec_path.display()))?;

    let proj = synthesize(&spec);

    // ---- prototype g-values ----------------------------------------------
    let g = &proj.prototype.g;
    let n = proj.prototype.order();
    println!("Filter synthesis ({:?}, order N={n})", spec.approximation);
    println!(
        "  f0 = {:.6e} Hz   FBW = {:.4}   Z0 = {} Ohm",
        spec.f0_hz, spec.fbw, spec.z0_ohm
    );
    print!("  prototype g-values: g0={:.4}", g[0]);
    for (i, gi) in g.iter().enumerate().skip(1) {
        print!("  g{i}={gi:.4}");
    }
    println!();

    // ---- coupling matrix + external Q ------------------------------------
    println!(
        "  external Q: Qe_in={:.4}  Qe_out={:.4}",
        proj.coupling.qe_in, proj.coupling.qe_out
    );
    println!("  coupling matrix M (normalized, {n}x{n}):");
    for row in &proj.coupling.m {
        let cells: Vec<String> = row.iter().map(|v| format!("{v:+.4}")).collect();
        println!("    [ {} ]", cells.join("  "));
    }

    // ---- ideal-response sweep --------------------------------------------
    let freqs = sweep_freqs(spec.f0_hz, spec.fbw);
    let s21 = ideal_response(&proj, &freqs);

    // ---- spec-mask verdict -----------------------------------------------
    // Grade over the same sweep so the in-band ripple/RL is well-sampled.
    let report = check_mask(&proj, &freqs);
    println!(
        "  mask: passband ripple {:.3} dB (spec {:.3}), in-band RL {:.3} dB (spec {:.3})",
        report.worst_passband_ripple_db,
        spec.mask.passband_ripple_db,
        report.worst_return_loss_db,
        spec.mask.return_loss_db,
    );
    for (f_hz, achieved, required, met) in &report.stopband {
        println!(
            "  stopband {f_hz:.4e} Hz: rejection {achieved:.2} dB (required {required:.2} dB) {}",
            if *met { "OK" } else { "UNDER" }
        );
    }
    for fail in &report.failures {
        println!("  FAILURE: {fail}");
    }

    // ---- write Touchstone .s2p -------------------------------------------
    let out_path = output
        .map(PathBuf::from)
        .unwrap_or_else(|| spec_path.with_extension("s2p"));
    write_s2p(&out_path, spec.z0_ohm, &freqs, &s21)
        .with_context(|| format!("failed to write Touchstone {}", out_path.display()))?;
    println!("  wrote Touchstone: {}", out_path.display());

    if report.pass {
        println!("VERDICT: PASS");
        Ok(ExitCode::SUCCESS)
    } else {
        println!("VERDICT: FAIL");
        Ok(ExitCode::FAILURE)
    }
}

/// Linear sweep of `SWEEP_POINTS` frequencies centred on `f0`, spanning
/// `f0·(1 ± SPAN_MULT·fbw/2)` (clamped to be strictly positive).
fn sweep_freqs(f0: f64, fbw: f64) -> Vec<f64> {
    let half = SPAN_MULT * fbw / 2.0;
    let lo = (f0 * (1.0 - half)).max(f0 * 1e-3);
    let hi = f0 * (1.0 + half);
    (0..SWEEP_POINTS)
        .map(|i| lo + (hi - lo) * (i as f64) / ((SWEEP_POINTS - 1) as f64))
        .collect()
}

/// Build and write a 2-port Touchstone file for a reciprocal, symmetric
/// lossless filter: `S21 = S12` from the ideal response, `S11 = S22` from
/// `|S11|² = 1 − |S21|²` (zero phase, magnitude model).
fn write_s2p(path: &Path, z0: f64, freqs: &[f64], s21: &[Complex64]) -> yee_io::Result<()> {
    let mut data = Vec::with_capacity(freqs.len());
    for s21f in s21 {
        let s21_mag = s21f.norm().min(1.0);
        let s11_mag = (1.0 - s21_mag * s21_mag).max(0.0).sqrt();
        let s11 = Complex64::new(s11_mag, 0.0);
        let s21c = Complex64::new(s21_mag, 0.0);
        // Row-major n×n: [S11, S12, S21, S22].
        data.push(vec![s11, s21c, s21c, s11]);
    }
    let file = File {
        n_ports: 2,
        z0,
        freq_unit: FreqUnit::Hz,
        format: Format::RealImag,
        freq_hz: freqs.to_vec(),
        data,
        comments: vec![" yee filter synth — ideal closed-form response".to_string()],
    };
    yee_io::touchstone::write(path, &file)
}
