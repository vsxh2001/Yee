//! Gate `cli-s2p-phase` for the **default** `yee filter synth` `.s2p`
//! (no `--q-unloaded`) — ADR-0174 T11 (the ADR-0172 T9 follow-on).
//!
//! T11 swaps the default `.s2p` response from the old
//! `ideal_response().map(lossless_s_pair)` (`|S21|` magnitude with `S11` placed
//! in lossless quadrature → **flat S21 phase**, `arg(S21) ≡ 0`) to the complex
//! coupling-matrix response `coupling_matrix_s_params` (the same complex S the
//! studio's distributed `.s2p` emits → CLI↔studio parity). This gate runs the
//! CLI against the standard 3-pole 0.5 dB Chebyshev bandpass (f0 = 2 GHz,
//! FBW = 0.10, Z0 = 50) **with no `--q-unloaded`**, **reads the written `.s2p`
//! back** via `yee_io::touchstone::read`, and asserts the emitted file's
//! properties (non-circular — the phase comes from the coupling-matrix solve;
//! the test only inspects the written Touchstone):
//!
//!   1. **Round-trip:** the `.s2p` parses cleanly via `yee_io::touchstone`
//!      (a 2-port file with a frequency grid).
//!   2. **Non-flat S21 phase:** `arg(S21)` varies meaningfully across the band
//!      (max − min > 0.5 rad). The OLD flat-phase default had `arg(S21) ≡ 0`;
//!      this is what distinguishes T11's complex response from it.
//!   3. **Passive (lossless):** `|S11|² + |S21|² ≤ 1 + ε` at every frequency,
//!      and `≈ 1` at midband (the coupling-matrix response is lossless by
//!      construction; `yee_io` would reject a super-unitary file on read-back).
//!
//! Pure closed-form (no EM/FDTD), so NOT `#[ignore]`'d.

use std::path::{Path, PathBuf};
use std::process::Command;

use num_complex::Complex64;

/// Standard 3-pole 0.5 dB Chebyshev bandpass (matches `cli_finite_q_s2p.rs`).
const F0_HZ: f64 = 2.0e9;
const FBW: f64 = 0.10;
const Z0_OHM: f64 = 50.0;
const ORDER: usize = 3;
const RIPPLE_DB: f64 = 0.5;

/// Write the standard 3-pole spec TOML into the cargo target tmp dir.
fn write_spec() -> PathBuf {
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_s2p_phase_spec.toml");
    let toml = format!(
        "response = \"Bandpass\"\n\
         f0_hz = {F0_HZ:e}\n\
         fbw = {FBW}\n\
         order = {ORDER}\n\
         z0_ohm = {Z0_OHM}\n\
         \n\
         [approximation.Chebyshev]\n\
         ripple_db = {RIPPLE_DB}\n\
         \n\
         [mask]\n\
         passband_ripple_db = {RIPPLE_DB}\n\
         return_loss_db = 9.0\n\
         stopband = [[2.4e9, 20.0]]\n"
    );
    std::fs::write(&path, toml).expect("write 3-pole spec TOML");
    path
}

/// Run `yee filter synth <spec> --output <out>` (DEFAULT path, no
/// `--q-unloaded`) and return `(success, stdout)`.
fn run_synth_default(out: &Path) -> (bool, String) {
    let spec = write_spec();
    let _ = std::fs::remove_file(out);
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "synth"])
        .arg(&spec)
        .arg("--output")
        .arg(out)
        .output()
        .expect("invoke yee");
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
    )
}

#[test]
fn cli_s2p_phase_default_carries_complex_phase() {
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_s2p_phase_out.s2p");
    let (ok, stdout) = run_synth_default(&out);
    assert!(ok, "default synth exited non-zero; stdout:\n{stdout}");
    assert!(
        out.exists(),
        "default Touchstone {} was not written",
        out.display()
    );

    // ---- (1) round-trip through yee_io::touchstone ----------------------
    let file = yee_io::touchstone::read(&out).expect("read written default .s2p");
    assert_eq!(file.n_ports, 2, ".s2p must be a 2-port file");
    assert!(!file.freq_hz.is_empty(), ".s2p has no frequency points");
    assert_eq!(
        file.data.len(),
        file.freq_hz.len(),
        "data rows must match frequency points"
    );

    // S11 is row-major index [0], S21 is index [1*2 + 0] = 2 of each 2×2 row.
    let s11: Vec<Complex64> = file.data.iter().map(|row| row[0]).collect();
    let s21: Vec<Complex64> = file.data.iter().map(|row| row[2]).collect();

    // ---- (2) non-flat S21 phase -----------------------------------------
    // The OLD default emitted S21 real (`arg ≡ 0`); the coupling-matrix
    // response carries physical phase. Use the in-band points (|S21| not
    // vanishingly small) so the argument is well-defined, then require the
    // phase to swing materially across the band.
    let args_inband: Vec<f64> = s21
        .iter()
        .filter(|z| z.norm() > 1e-3)
        .map(|z| z.arg())
        .collect();
    assert!(
        !args_inband.is_empty(),
        "no in-band S21 points with non-negligible magnitude"
    );
    let arg_min = args_inband.iter().cloned().fold(f64::INFINITY, f64::min);
    let arg_max = args_inband
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let arg_span = arg_max - arg_min;

    // The largest |arg(S21)| across the band — the old flat-phase output had
    // arg(S21) ≡ 0, so this is also strictly positive only for the complex S.
    let max_abs_arg = s21
        .iter()
        .filter(|z| z.norm() > 1e-3)
        .map(|z| z.arg().abs())
        .fold(0.0_f64, f64::max);

    println!(
        "cli-s2p-phase: n_freq = {} | arg(S21) span = {arg_span:.4} rad \
         (min {arg_min:.4}, max {arg_max:.4}) | max|arg(S21)| = {max_abs_arg:.4} rad",
        file.freq_hz.len()
    );

    assert!(
        arg_span > 0.5,
        "default .s2p S21 phase is too flat: arg(S21) spans only {arg_span:.4} rad \
         (≤ 0.5) — the OLD flat-phase output had arg(S21) ≡ 0; T11 should carry \
         physical phase"
    );

    // ---- (3) passive everywhere + lossless (≈1) at midband --------------
    let mut worst_power = 0.0_f64;
    for (i, (a, b)) in s11.iter().zip(s21.iter()).enumerate() {
        let p = a.norm_sqr() + b.norm_sqr();
        worst_power = worst_power.max(p);
        assert!(
            p <= 1.0 + 1e-6,
            "default .s2p is not passive at f = {:.4e} Hz (idx {i}): \
             |S11|²+|S21|² = {p:.6} > 1",
            file.freq_hz[i]
        );
    }

    // At midband the lossless coupling-matrix response conserves power (≈1) —
    // the complement of the finite-Q path's absorption (`cli_finite_q_s2p`).
    let mid_idx = file
        .freq_hz
        .iter()
        .enumerate()
        .min_by(|(_, x), (_, y)| (*x - F0_HZ).abs().partial_cmp(&(*y - F0_HZ).abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let mid_power = s11[mid_idx].norm_sqr() + s21[mid_idx].norm_sqr();
    println!(
        "cli-s2p-phase: worst |S11|²+|S21|² = {worst_power:.6} (≤ 1) | \
         midband |S11|²+|S21|² = {mid_power:.6} (≈ 1, lossless)"
    );
    assert!(
        (mid_power - 1.0).abs() < 1e-3,
        "midband |S11|²+|S21|² = {mid_power:.6} is not ≈ 1 (lossless coupling-matrix \
         response expected)"
    );

    let _ = std::fs::remove_file(&out);
}
