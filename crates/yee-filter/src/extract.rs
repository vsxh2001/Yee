//! Coupling/Qe extraction (Filter Phase F1.1b.0).
//!
//! Pure-DSP inversion that turns an EM/measured response into the two filter
//! design parameters the coupled-resonator synthesis trades in:
//!
//! - [`extract_coupling`] reads a coupling coefficient `k` from the two split
//!   resonance peaks of a synchronously-tuned coupled-resonator pair, and
//! - [`extract_q_ringdown`] reads the loaded/external quality factor `Q` from a
//!   resonator ring-down time series.
//!
//! Both are the inverse of synthesis (measured response → coupling) and are the
//! exact API the F1.1b.1 FDTD driver will call. The module is pure `f64`, takes
//! no new dependency, and is WASM-safe (ADR-0093). It is validated against
//! *analytic* signals only — no FDTD (gates `extract-001` / `extract-002`).
//!
//! # Method
//!
//! Coupling (Pozar §8.8, Hong & Lancaster ch. 8): for two synchronously-tuned
//! resonators at centre `f0` with coupling `k`, the response splits into two
//! peaks at `f_lo = f0/√(1+k)` and `f_hi = f0/√(1−k)`. The coupling inverts
//! from the split peaks as
//!
//! ```text
//! k = (f_hi² − f_lo²) / (f_hi² + f_lo²).
//! ```
//!
//! Q (Pozar §6.1; mirrors `yee-fdtd`'s `cavity_q.rs` decay fit): the decaying
//! upper envelope of the ring-down is fit log-linearly,
//! `ln|env| = a − t/τ`, giving the time constant `τ = −1/slope` and
//! `Q = π · f0 · τ`.

/// The coupling-coefficient extraction: the two split resonance peaks and the
/// coupling `k` inverted from them.
#[derive(Debug, Clone, PartialEq)]
pub struct CouplingExtraction {
    /// Lower split-peak frequency, Hz (`f_lo < f_hi`).
    pub f_lo_hz: f64,
    /// Upper split-peak frequency, Hz (`f_hi > f_lo`).
    pub f_hi_hz: f64,
    /// Coupling coefficient `k = (f_hi² − f_lo²) / (f_hi² + f_lo²)`.
    pub k: f64,
}

/// Extract the coupling coefficient `k` from the two split resonance peaks of a
/// synchronously-tuned coupled-resonator pair (Pozar §8.8).
///
/// Finds the interior local maxima of `mag` (index `i` in `1..n-1` with
/// `mag[i] > mag[i-1]` and `mag[i] > mag[i+1]`), takes the two with the largest
/// magnitude, orders their frequencies `f_lo < f_hi`, and returns
/// `k = (f_hi² − f_lo²) / (f_hi² + f_lo²)`.
///
/// Returns `None` if `freqs_hz.len() != mag.len()`, the length is `< 5`, or
/// fewer than two distinct interior local maxima exist (e.g. a single peak).
pub fn extract_coupling(freqs_hz: &[f64], mag: &[f64]) -> Option<CouplingExtraction> {
    if freqs_hz.len() != mag.len() || mag.len() < 5 {
        return None;
    }

    // Interior local maxima: strictly greater than both neighbours.
    let mut peaks: Vec<usize> = (1..mag.len() - 1)
        .filter(|&i| mag[i] > mag[i - 1] && mag[i] > mag[i + 1])
        .collect();
    if peaks.len() < 2 {
        return None;
    }

    // Keep the two largest-magnitude peaks. Sort descending by magnitude; ties
    // resolve by index so the selection is deterministic.
    peaks.sort_by(|&a, &b| {
        mag[b]
            .partial_cmp(&mag[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });
    let (i0, i1) = (peaks[0], peaks[1]);

    // Order by frequency: f_lo < f_hi.
    let (f_a, f_b) = (freqs_hz[i0], freqs_hz[i1]);
    let (f_lo_hz, f_hi_hz) = if f_a <= f_b { (f_a, f_b) } else { (f_b, f_a) };

    let k = (f_hi_hz * f_hi_hz - f_lo_hz * f_lo_hz) / (f_hi_hz * f_hi_hz + f_lo_hz * f_lo_hz);

    Some(CouplingExtraction {
        f_lo_hz,
        f_hi_hz,
        k,
    })
}

/// Extract the loaded/external quality factor `Q` from a resonator ring-down
/// time series (Pozar §6.1; mirrors `yee-fdtd`'s `cavity_q.rs` decay fit).
///
/// Builds the decaying upper envelope as the interior local maxima of
/// `|samples|` (index `i` in `1..n-1` that is a local maximum of magnitude),
/// skipping any in the first 1/3 of the record (the initial transient, matching
/// `cavity_q.rs`). Fits `ln|env_k| = a − t_k/τ` by ordinary least squares with
/// `t_k = i_k · dt_s`; a negative slope `m < 0` gives `τ = −1/m` and
/// `Q = π · f0_hz · τ`.
///
/// Returns `None` if fewer than 3 strictly-positive envelope points survive,
/// the fit is degenerate (no spread in `t` — all-equal magnitudes), or the
/// fitted slope `m ≥ 0` (no decay).
pub fn extract_q_ringdown(samples: &[f64], dt_s: f64, f0_hz: f64) -> Option<f64> {
    let n = samples.len();
    if n < 3 {
        return None;
    }

    // Skip the first 1/3 (initial transient) before envelope detection.
    let skip = n / 3;

    // Upper envelope: interior local maxima of |samples| past the transient,
    // keeping only strictly-positive magnitudes (ln-fit needs |env| > 0).
    let mut ts: Vec<f64> = Vec::new();
    let mut ln_env: Vec<f64> = Vec::new();
    for i in 1..n - 1 {
        if i < skip {
            continue;
        }
        let m = samples[i].abs();
        let is_local_max = m >= samples[i - 1].abs() && m >= samples[i + 1].abs();
        if is_local_max && m > 0.0 {
            ts.push(i as f64 * dt_s);
            ln_env.push(m.ln());
        }
    }

    if ts.len() < 3 {
        return None;
    }

    // Ordinary least squares: ln|env| = a + m·t (slope m = −1/τ).
    let count = ts.len() as f64;
    let t_mean = ts.iter().sum::<f64>() / count;
    let y_mean = ln_env.iter().sum::<f64>() / count;

    let num: f64 = ts
        .iter()
        .zip(ln_env.iter())
        .map(|(&t, &y)| (t - t_mean) * (y - y_mean))
        .sum();
    let den: f64 = ts.iter().map(|&t| (t - t_mean).powi(2)).sum();

    // Degenerate: no spread in t (e.g. all-equal magnitudes collapse onto one
    // bin, or every envelope sample shares a time).
    if den <= 0.0 {
        return None;
    }

    let slope = num / den;
    // No decay (flat or growing) → not a ring-down.
    if slope >= 0.0 {
        return None;
    }

    let tau = -1.0 / slope;
    Some(std::f64::consts::PI * f0_hz * tau)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Single Lorentzian peak `1 / (1 + (2·Q·(f−fc)/fc)²)`.
    fn lorentzian(f: f64, fc: f64, q: f64) -> f64 {
        let x = 2.0 * q * (f - fc) / fc;
        1.0 / (1.0 + x * x)
    }

    /// `extract-001` — two-Lorentzian split from a known `k_true` round-trips
    /// back to `k_true` within the peak-bin tolerance; a single peak → `None`.
    #[test]
    fn extract_001_coupling() {
        let f0 = 2.0e9_f64;
        let k_true = 0.04_f64;
        // Hong-Lancaster split frequencies.
        let f_lo = f0 / (1.0 + k_true).sqrt();
        let f_hi = f0 / (1.0 - k_true).sqrt();
        let q = 100.0_f64; // peaks clearly separated at k = 0.04.

        // 0.7–1.3·f0 sweep, 401 points.
        let n = 401usize;
        let lo = 0.7 * f0;
        let hi = 1.3 * f0;
        let freqs: Vec<f64> = (0..n)
            .map(|i| lo + (hi - lo) * (i as f64) / ((n - 1) as f64))
            .collect();
        let mag: Vec<f64> = freqs
            .iter()
            .map(|&f| lorentzian(f, f_lo, q) + lorentzian(f, f_hi, q))
            .collect();

        let got = extract_coupling(&freqs, &mag).expect("two peaks → Some");
        assert!(
            got.f_lo_hz < got.f_hi_hz,
            "f_lo {} should be < f_hi {}",
            got.f_lo_hz,
            got.f_hi_hz
        );
        let err = (got.k - k_true).abs();
        eprintln!(
            "extract-001: k_recovered = {:.6}, k_true = {:.6}, |err| = {:.3e} \
             (f_lo = {:.6e}, f_hi = {:.6e})",
            got.k, k_true, err, got.f_lo_hz, got.f_hi_hz
        );
        assert!(
            err <= 1e-2,
            "recovered k = {:.6} vs k_true = {:.6}, |err| = {:.3e} exceeds 1e-2",
            got.k,
            k_true,
            err
        );

        // Negative control: a single Lorentzian has one interior maximum → None.
        let single: Vec<f64> = freqs.iter().map(|&f| lorentzian(f, f0, q)).collect();
        assert!(
            extract_coupling(&freqs, &single).is_none(),
            "single Lorentzian should yield None"
        );
    }

    /// `extract-002` — `Q` from `exp(−t/τ)·sin(2π f0 t)` ring-down within 3 % of
    /// `π f0 τ`; a constant (no decay) → `None`.
    #[test]
    fn extract_002_q_ringdown() {
        let f0 = 2.0e9_f64;
        let tau = 5.0e-9_f64;
        let dt = 1.0 / (40.0 * f0); // ~40 samples per RF cycle.
        // ~8τ of record.
        let n = (8.0 * tau / dt).ceil() as usize;
        let samples: Vec<f64> = (0..n)
            .map(|i| {
                let t = i as f64 * dt;
                (-t / tau).exp() * (2.0 * PI * f0 * t).sin()
            })
            .collect();

        let q = extract_q_ringdown(&samples, dt, f0).expect("decaying ring-down → Some");
        let q_expected = PI * f0 * tau;
        let rel_err = (q - q_expected).abs() / q_expected;
        eprintln!(
            "extract-002: Q_recovered = {:.6}, Q_expected = π·f0·τ = {:.6}, \
             rel_err = {:.4} %",
            q,
            q_expected,
            rel_err * 100.0
        );
        assert!(
            rel_err <= 0.03,
            "recovered Q = {:.6} vs expected {:.6}, rel_err = {:.4} % exceeds 3 %",
            q,
            q_expected,
            rel_err * 100.0
        );

        // Negative control: a constant (no decay) → None.
        let constant = vec![1.0_f64; n];
        assert!(
            extract_q_ringdown(&constant, dt, f0).is_none(),
            "constant samples (no decay) should yield None"
        );
    }
}
