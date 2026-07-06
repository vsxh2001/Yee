//! S-parameter post-processing over job results (S.6, ADR-0183).
//!
//! Pure functions on the per-step probe series a [`crate::JobResult`]
//! carries, so every client of the job protocol — studio, Python, WS,
//! tests — extracts spectra the same way. The workhorse is the two-run
//! transmission ratio (Sheen et al. 1990, adapted to lumped ports):
//! run the bare feed line (reference) and the device (DUT) as two
//! otherwise-identical jobs, then
//! `|S21|(f) = |DFT(dut)(f)| / |DFT(reference)(f)|` — feed-line loss,
//! launch discontinuity, and probe coupling divide out.

/// Complex single-bin DFT of a uniformly sampled series:
/// `X(f) = Σₙ x[n]·e^{−j·2πf·n·dt}`, returned as `(re, im)`.
///
/// The same single-bin correlation the ε_eff gates use, exposed as a
/// reusable function. No windowing and no normalization — ratios of bins
/// taken with the same `series.len()` cancel both.
pub fn single_bin_dft(series: &[f64], dt_s: f64, f_hz: f64) -> (f64, f64) {
    let omega = std::f64::consts::TAU * f_hz;
    let mut re = 0.0;
    let mut im = 0.0;
    for (n, x) in series.iter().enumerate() {
        let phase = omega * n as f64 * dt_s;
        re += x * phase.cos();
        im -= x * phase.sin();
    }
    (re, im)
}

/// Transmission magnitude in dB at each requested frequency:
/// `20·log₁₀(|DFT(dut)(f)| / |DFT(reference)(f)|)`.
///
/// Both series must be sampled with the same `dt_s` (use the
/// [`crate::JobResult::dt_s`] the runs report) and should have the same
/// length so the un-normalized bins cancel. Frequencies must lie inside
/// the drive's spectral band: a reference bin with no drive energy
/// produces `+∞`/NaN, which is a caller error, not a signal.
pub fn transmission_db(dut: &[f64], reference: &[f64], dt_s: f64, freqs_hz: &[f64]) -> Vec<f64> {
    freqs_hz
        .iter()
        .map(|&f| {
            let (dr, di) = single_bin_dft(dut, dt_s, f);
            let (rr, ri) = single_bin_dft(reference, dt_s, f);
            let mag_dut = dr.hypot(di);
            let mag_ref = rr.hypot(ri);
            20.0 * (mag_dut / mag_ref).log10()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn single_bin_dft_recovers_a_known_sinusoid() {
        // A·sin(2πf t) over exactly 8 periods: |X| = A·N/2, phase −90°
        // (sin = (e^{jθ} − e^{−jθ})/2j picks up −j at the +f bin).
        let f = 2.0e9;
        let dt = 1.0 / (f * 64.0); // 64 samples per period
        let n = 64 * 8;
        let a = 3.5;
        let series: Vec<f64> = (0..n)
            .map(|i| a * (TAU * f * i as f64 * dt).sin())
            .collect();
        let (re, im) = single_bin_dft(&series, dt, f);
        let mag = re.hypot(im);
        assert!((mag - a * n as f64 / 2.0).abs() / (a * n as f64 / 2.0) < 1e-9);
        assert!(re.abs() < 1e-6 * mag, "real part should vanish: {re}");
        assert!(im < 0.0, "positive-frequency bin of sin carries −j");
    }

    #[test]
    fn transmission_of_a_half_scaled_copy_is_minus_6_db() {
        let f0 = 2.0e9;
        let dt = 1.0 / (f0 * 64.0);
        let reference: Vec<f64> = (0..1024)
            .map(|i| {
                let t = i as f64 * dt;
                (-((t - 512.0 * dt) / (128.0 * dt)).powi(2)).exp() * (TAU * f0 * t).sin()
            })
            .collect();
        let dut: Vec<f64> = reference.iter().map(|x| 0.5 * x).collect();
        let freqs = [0.8 * f0, f0, 1.2 * f0];
        for db in transmission_db(&dut, &reference, dt, &freqs) {
            assert!((db - 20.0 * 0.5_f64.log10()).abs() < 1e-9, "got {db} dB");
        }
    }
}
