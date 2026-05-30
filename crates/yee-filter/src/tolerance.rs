//! Monte-Carlo **tolerance / yield** analysis (Filter Phase F2.4).
//!
//! Real lumped components have a manufacturing tolerance (E24 ±5 %, E96 ±1 %),
//! so the as-built filter's spec compliance is a *distribution*, not a point.
//! This module quantifies it with a seeded Monte-Carlo **yield**: it snaps each
//! resonator's `L`/`C` to the chosen E-series value (F2.1), then for `M` random
//! trials perturbs every value uniformly within `±tolerance`, rebuilds the LC
//! ladder, evaluates the realized response ([`crate::ladder_s21`]) against the
//! spec mask ([`crate::mask_verdict`]), and reports the **fraction that passes**
//! plus worst-case in-band return loss / stopband rejection across all trials.
//!
//! Pure `f64` + serde, WASM-safe, **NO FDTD and NO `rand` dependency**: the
//! randomness is a tiny in-module **seeded PRNG** (SplitMix64), so results are
//! reproducible (same `seed` → identical yield) and the crate stays dep-free.
//! It builds only on the shipped F2.0 ([`LumpedLadder`] + `ladder_s21`) and F2.1
//! ([`ESeries::nearest`] + per-series tolerance).
//!
//! # Method
//!
//! 1. **Realize:** snap each [`crate::LcResonator`]'s `l_henry`/`c_farad` to the
//!    nearest [`ESeries`] value ([`ESeries::nearest`]).
//! 2. **Sample (`M` seeded trials):** perturb each realized value by a uniform
//!    random factor `1 + tol·(2·u − 1)` with `u ∈ [0,1)` and
//!    `tol = ESeries::tolerance_pct/100`.
//! 3. **Evaluate:** rebuild the perturbed [`LumpedLadder`], sweep `ladder_s21`
//!    over the band, apply [`crate::mask_verdict`] (strict — no realization
//!    slack).
//! 4. **Aggregate:** `yield_fraction = passes / M`, plus the worst in-band
//!    return loss and worst stopband rejection seen across *all* trials.
//!
//! Out of scope (documented follow-ons): correlated / √-N statistics,
//! per-component sensitivity ranking, σ/Cpk, non-uniform (Gaussian) part
//! distributions, parasitics (F2.1b), and FDTD-based yield (would use F2.3).

use serde::{Deserialize, Serialize};

use crate::{ESeries, LumpedLadder, SpecMask, mask_verdict};

/// Number of frequency samples in the per-trial passband/stopband sweep.
///
/// Matches the `lumped_001` realization gate's resolution so the in-band ripple
/// is well sampled; the sweep spans `f0·(1 ± 3·fbw)`.
const SWEEP_POINTS: usize = 801;

/// Half-span of the per-trial sweep, as a multiple of the fractional bandwidth.
/// The sweep runs over `f0·(1 ± SWEEP_SPAN_FBW·fbw)`.
const SWEEP_SPAN_FBW: f64 = 3.0;

/// A tiny seeded **SplitMix64** PRNG.
///
/// Pure-`f64` / `u64`, WASM-safe, dependency-free, and fully reproducible from
/// its seed — the F2.4 yield analysis uses it instead of the `rand` crate so the
/// crate stays dep-free and the gate is deterministic. SplitMix64 is the
/// standard fast splittable generator (the seeding routine recommended for
/// `xoshiro`); it is *not* cryptographic, which is exactly right for a
/// reproducible tolerance Monte-Carlo.
struct Rng {
    /// The 64-bit generator state, advanced by the golden-ratio increment.
    state: u64,
}

impl Rng {
    /// Seed the generator. Any seed is valid (including `0`).
    fn new(seed: u64) -> Self {
        Rng { state: seed }
    }

    /// Next 64-bit output (the SplitMix64 mixing function).
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Next uniform `f64` in `[0, 1)`, using the top 53 bits (the f64 mantissa
    /// width) so every representable value in `[0,1)` is reachable.
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// The result of a Monte-Carlo tolerance / yield analysis.
///
/// Produced by [`monte_carlo_yield`]. `yield_fraction` is the share of trials
/// whose realized response met the spec mask; the worst-case fields are the
/// extremes observed across **all** trials (pass or fail), useful as a margin
/// readout.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct YieldResult {
    /// Fraction of trials that passed the spec mask, in `[0, 1]`
    /// (`passes / n_samples`).
    pub yield_fraction: f64,
    /// Number of Monte-Carlo trials run.
    pub n_samples: usize,
    /// Worst-case (smallest) in-band return loss observed across all trials, dB.
    pub worst_inband_rl_db: f64,
    /// Worst-case (smallest) stopband rejection observed across all trials, dB.
    pub worst_stopband_rej_db: f64,
}

/// Run a seeded Monte-Carlo tolerance / yield analysis on a [`LumpedLadder`].
///
/// **Realizes** the ladder by snapping each resonator's `l_henry`/`c_farad` to
/// the nearest `series` value ([`ESeries::nearest`]), then runs `n_samples`
/// reproducible trials seeded by `seed`: each trial perturbs every realized
/// `L`/`C` by an independent uniform factor `1 + tol·(2·u − 1)` (with
/// `tol = series.tolerance_pct()/100` and `u ∈ [0,1)` from the in-module
/// SplitMix64 PRNG), rebuilds the ladder, and grades its realized `ladder_s21`
/// response against `mask` via [`mask_verdict`] (strict — no realization slack).
///
/// Returns the [`YieldResult`]: `yield_fraction = passes / n_samples` plus the
/// worst in-band return loss and worst stopband rejection seen across all
/// trials. The same `(ladder, series, mask, n_samples, seed)` always yields an
/// identical result (the PRNG is fully determined by `seed`).
///
/// With `n_samples == 0` the `yield_fraction` is `0.0` and the worst-case fields
/// are their neutral extremes (`+∞` rejection, `+∞` return loss reported as the
/// `f64::INFINITY` sentinels), since no trial was evaluated.
pub fn monte_carlo_yield(
    ladder: &LumpedLadder,
    series: ESeries,
    mask: &SpecMask,
    n_samples: usize,
    seed: u64,
) -> YieldResult {
    // --- realize: snap every L/C to the nearest E-series value -------------
    let mut realized = ladder.clone();
    for r in &mut realized.resonators {
        r.l_henry = series.nearest(r.l_henry);
        r.c_farad = series.nearest(r.c_farad);
    }

    let tol = series.tolerance_pct() / 100.0;
    let f0 = realized.f0_hz;
    let fbw = realized.fbw;
    let freqs: Vec<f64> = {
        let lo = (f0 * (1.0 - SWEEP_SPAN_FBW * fbw)).max(f0 * 1e-3);
        let hi = f0 * (1.0 + SWEEP_SPAN_FBW * fbw);
        (0..SWEEP_POINTS)
            .map(|i| lo + (hi - lo) * (i as f64) / ((SWEEP_POINTS - 1) as f64))
            .collect()
    };

    let mut rng = Rng::new(seed);
    let mut passes = 0usize;
    let mut worst_rl = f64::INFINITY;
    let mut worst_rej = f64::INFINITY;

    // Reuse one ladder buffer across trials, overwriting each element value.
    let mut sample = realized.clone();
    for _ in 0..n_samples {
        for (s, base) in sample.resonators.iter_mut().zip(&realized.resonators) {
            s.l_henry = base.l_henry * (1.0 + tol * (2.0 * rng.next_f64() - 1.0));
            s.c_farad = base.c_farad * (1.0 + tol * (2.0 * rng.next_f64() - 1.0));
        }
        // Strict verdict (ripple slack 0.0): real parts are graded at full
        // strength, no realization margin.
        let v = mask_verdict(&sample, mask, f0, fbw, &freqs, 0.0);
        if v.pass {
            passes += 1;
        }
        worst_rl = worst_rl.min(v.worst_return_loss_db);
        worst_rej = worst_rej.min(v.worst_stopband_rej_db);
    }

    let yield_fraction = if n_samples == 0 {
        0.0
    } else {
        passes as f64 / n_samples as f64
    };

    YieldResult {
        yield_fraction,
        n_samples,
        worst_inband_rl_db: worst_rl,
        worst_stopband_rej_db: worst_rej,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_is_deterministic_and_in_range() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..1000 {
            let x = a.next_f64();
            assert_eq!(x, b.next_f64(), "same seed must reproduce the stream");
            assert!((0.0..1.0).contains(&x), "next_f64 out of [0,1): {x}");
        }
        // A different seed gives a different stream.
        let mut c = Rng::new(43);
        assert_ne!(Rng::new(42).next_u64(), c.next_u64());
    }
}
