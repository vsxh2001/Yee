//! S-parameter post-processing over job results (S.6, ADR-0183).
//!
//! Pure functions on the per-step probe series a [`crate::JobResult`]
//! carries, so every client of the job protocol — studio, Python, WS,
//! tests — extracts spectra the same way. The workhorse is the two-run
//! transmission ratio (Sheen et al. 1990, adapted to lumped ports):
//! run the bare feed line (reference) and the device (DUT) as two
//! otherwise-identical jobs, then
//! `|S21|(f) = |DFT(dut)(f)| / |DFT(reference)(f)|` — feed-line loss,
//! launch discontinuity, and probe coupling divide out. The same two
//! runs also yield |S11| ([`reflection_db`], S.7): the reference run's
//! port-1 probe is the incident wave, so subtracting it from the DUT
//! run's isolates the device-caused reflection.

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

/// Reflection magnitude in dB at each requested frequency, via
/// incident/reflected separation at the port-1 reference plane (S.7):
/// `reflected(t) = dut_p1(t) − ref_p1(t)`, then
/// `20·log₁₀(|DFT(reflected)(f)| / |DFT(ref_p1)(f)|)`.
///
/// `ref_p1` is the same probe recorded on the reference (bare-line) run —
/// the two runs share launch, line, and grid, so that series **is** the
/// incident wave at the plane, and the sample-wise difference isolates
/// the device-caused reflection. Second-order caveat: the reflected wave
/// re-reflects off the imperfectly matched drive port and passes the
/// plane again; accepted at walking-skeleton tolerance.
///
/// # Panics
///
/// Panics if the series lengths differ — they must come from the same
/// pair of jobs.
pub fn reflection_db(dut_p1: &[f64], ref_p1: &[f64], dt_s: f64, freqs_hz: &[f64]) -> Vec<f64> {
    assert_eq!(
        dut_p1.len(),
        ref_p1.len(),
        "reflection_db: P1 series lengths differ — not the same job pair"
    );
    let reflected: Vec<f64> = dut_p1.iter().zip(ref_p1).map(|(d, r)| d - r).collect();
    transmission_db(&reflected, ref_p1, dt_s, freqs_hz)
}

/// Forward/backward wave split from a three-probe standing-wave fit
/// (S.12, ADR-0189) — a verbatim port of `yee-voxel`'s `fit_standing_wave`
/// (F2.3-h, ADR-0129/0131/0133), the repo's established cure for the
/// standing-wave / over-unity artifact: with three equally spaced phasors
/// `V₀,V₁,V₂` on a line carrying `V(x) = a·e^{−jβx} + b·e^{+jβx}`,
/// `cos(βd) = (V₀+V₂)/(2V₁)` recovers β, then a linear solve splits the
/// waves. All complex numbers are `(re, im)` tuples.
#[derive(Debug, Clone, Copy)]
pub struct WaveSplit {
    /// Forward (+x) wave phasor at the first probe's plane.
    pub fwd: (f64, f64),
    /// Backward (−x) wave phasor at the first probe's plane.
    pub bwd: (f64, f64),
    /// Fitted propagation constant β (rad/m).
    pub beta_rad_m: f64,
    /// |Im cos(βd)| — a consistency residual (→ 0 for a clean fit).
    pub residual: f64,
}

fn cadd(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 + b.0, a.1 + b.1)
}
fn csub(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 - b.0, a.1 - b.1)
}
fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn cdiv(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let n = b.0 * b.0 + b.1 * b.1;
    ((a.0 * b.0 + a.1 * b.1) / n, (a.1 * b.0 - a.0 * b.1) / n)
}
fn cabs(a: (f64, f64)) -> f64 {
    a.0.hypot(a.1)
}
fn cexpj(theta: f64) -> (f64, f64) {
    (theta.cos(), theta.sin())
}

/// Fit `V(x) = a·e^{−jβx} + b·e^{+jβx}` to three equally spaced phasors
/// (`spacing_m` apart, ordered along +x). Returns the split at `v0`'s
/// plane. Degenerate spacings (βd ≈ 0 or π) fall back to all-forward,
/// flagged by the residual.
pub fn fit_standing_wave(
    v0: (f64, f64),
    v1: (f64, f64),
    v2: (f64, f64),
    spacing_m: f64,
) -> WaveSplit {
    let cos_bd = cdiv(cadd(v0, v2), (2.0 * v1.0, 2.0 * v1.1));
    let residual = cos_bd.1.abs();
    let c = cos_bd.0.clamp(-1.0, 1.0);
    let beta_d = c.acos();
    let beta_rad_m = beta_d / spacing_m;

    let p = cexpj(-2.0 * beta_d);
    let q = cexpj(2.0 * beta_d);
    let denom = csub(q, p);
    let (fwd, bwd) = if cabs(denom) > 1e-12 {
        let a = cdiv(csub(cmul(v0, q), v2), denom);
        let b = csub(v0, a);
        (a, b)
    } else {
        (v0, (0.0, 0.0))
    };
    WaveSplit {
        fwd,
        bwd,
        beta_rad_m,
        residual,
    }
}

/// **Directional** transmission magnitude in dB: at each frequency, DFT
/// three equally spaced probes per run, split forward/backward waves via
/// [`fit_standing_wave`], and ratio the **forward** amplitudes —
/// `20·log₁₀(|fwd_dut| / |fwd_ref|)`. Immune to the reflected wave that
/// makes the plain single-probe [`transmission_db`] ripple when ports
/// are imperfectly matched (the ADR-0188 finding). Probe triples must be
/// ordered along +x with the same spacing in both runs.
pub fn directional_transmission_db(
    dut: [&[f64]; 3],
    reference: [&[f64]; 3],
    dt_s: f64,
    spacing_m: f64,
    freqs_hz: &[f64],
) -> Vec<f64> {
    freqs_hz
        .iter()
        .map(|&f| {
            let d: Vec<(f64, f64)> = dut.iter().map(|s| single_bin_dft(s, dt_s, f)).collect();
            let r: Vec<(f64, f64)> = reference
                .iter()
                .map(|s| single_bin_dft(s, dt_s, f))
                .collect();
            let d_split = fit_standing_wave(d[0], d[1], d[2], spacing_m);
            let r_split = fit_standing_wave(r[0], r[1], r[2], spacing_m);
            20.0 * (cabs(d_split.fwd) / cabs(r_split.fwd)).log10()
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

    #[test]
    fn standing_wave_fit_recovers_known_forward_and_backward_waves() {
        // V(x) = a·e^{−jβx} + b·e^{+jβx} with known a, b, β.
        let beta = 120.0; // rad/m
        let d = 5.0e-3; // 5 mm spacing → βd = 0.6 rad
        let a = (0.8, -0.3);
        let b = (0.25, 0.15);
        let v_at = |x: f64| cadd(cmul(a, cexpj(-beta * x)), cmul(b, cexpj(beta * x)));
        let split = fit_standing_wave(v_at(0.0), v_at(d), v_at(2.0 * d), d);
        assert!(
            (split.beta_rad_m - beta).abs() < 1e-9,
            "{}",
            split.beta_rad_m
        );
        assert!(cabs(csub(split.fwd, a)) < 1e-9);
        assert!(cabs(csub(split.bwd, b)) < 1e-9);
        assert!(split.residual < 1e-12);
    }

    #[test]
    fn directional_ratio_ignores_a_backward_wave() {
        // DUT carries 0.5× the reference's forward wave PLUS a strong
        // backward wave; the directional ratio must still read −6.02 dB.
        let f0 = 2.0e9;
        let dt = 1.0 / (f0 * 64.0);
        let beta = 130.0;
        let d = 5.0e-3;
        let n = 1024;
        let series = |a_amp: f64, b_amp: f64, x: f64| -> Vec<f64> {
            (0..n)
                .map(|i| {
                    let t = i as f64 * dt;
                    let env = (-((t - 512.0 * dt) / (170.0 * dt)).powi(2)).exp();
                    env * (a_amp * (TAU * f0 * t - beta * x).sin()
                        + b_amp * (TAU * f0 * t + beta * x).sin())
                })
                .collect()
        };
        let reference: Vec<Vec<f64>> = (0..3).map(|m| series(1.0, 0.0, m as f64 * d)).collect();
        let dut: Vec<Vec<f64>> = (0..3).map(|m| series(0.5, 0.4, m as f64 * d)).collect();
        let freqs = [f0];
        let db = directional_transmission_db(
            [&dut[0], &dut[1], &dut[2]],
            [&reference[0], &reference[1], &reference[2]],
            dt,
            d,
            &freqs,
        );
        assert!(
            (db[0] - 20.0 * 0.5_f64.log10()).abs() < 0.1,
            "directional ratio {} dB, want −6.02",
            db[0]
        );
    }

    #[test]
    fn reflection_of_a_synthetic_quarter_echo_is_minus_12_db() {
        // dut_p1 = incident + 0.25·incident → reflected/incident = 0.25
        // exactly, i.e. −12.04 dB at every in-band frequency.
        let f0 = 2.0e9;
        let dt = 1.0 / (f0 * 64.0);
        let incident: Vec<f64> = (0..1024)
            .map(|i| {
                let t = i as f64 * dt;
                (-((t - 512.0 * dt) / (128.0 * dt)).powi(2)).exp() * (TAU * f0 * t).sin()
            })
            .collect();
        let dut_p1: Vec<f64> = incident.iter().map(|x| 1.25 * x).collect();
        let freqs = [0.8 * f0, f0, 1.2 * f0];
        for db in reflection_db(&dut_p1, &incident, dt, &freqs) {
            assert!((db - 20.0 * 0.25_f64.log10()).abs() < 1e-9, "got {db} dB");
        }
    }
}
