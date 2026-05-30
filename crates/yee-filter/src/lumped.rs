//! Closed-form lumped-element **LC ladder synthesis** (Filter Phase F2.0).
//!
//! Turns an abstract synthesized lowpass prototype (the `yee-synth` g-values
//! carried on [`FilterProject::prototype`](crate::FilterProject)) into **ideal
//! `L`/`C` element values** for a lumped band-pass filter, by applying the
//! textbook low-pass-prototype → band-pass ladder transform. Pure `f64`,
//! WASM-safe, NO FDTD, NO parts/PCB — this is the *initial* ideal-element
//! realization that seeds the later component-selection (F2.1), lumped-EM
//! (F2.3), and tolerance (F2.4) phases, mirroring the
//! [`crate::dimension`] (`synthesize` → physical-dimensions) pattern on the
//! distributed track.
//!
//! # Method (Pozar §8.3 / Hong & Lancaster ch. 3)
//!
//! From the low-pass prototype `g0, g1, …, gN, g_{N+1}` (already produced by
//! [`crate::synthesize`]), the centre `ω0 = 2π·f0`, the fractional bandwidth
//! `Δ = FBW`, and the system `Z0`, the standard low-pass → band-pass ladder
//! transform maps **each** reactive prototype element `g_k` (`k = 1..N`) to a
//! **series** or **shunt** LC resonator, alternating along the ladder:
//!
//! - **Series-branch resonator** (a series L–C): `L_k = g_k·Z0/(ω0·Δ)`,
//!   `C_k = Δ/(ω0·Z0·g_k)`.
//! - **Shunt-branch resonator** (a parallel L–C): `L_k = Z0·Δ/(ω0·g_k)`,
//!   `C_k = g_k/(ω0·Z0·Δ)`.
//!
//! Either way `L_k·C_k = 1/ω0²`, so **every resonator is tuned to `ω0`** — the
//! defining property the [`lumped_001`](../../tests/lumped_001.rs) gate checks.
//!
//! ## Branch-ordering convention (shunt-first)
//!
//! The first reactive element `g1` is realized as a **shunt** resonator and the
//! branches alternate from there — resonator `k` (1-based) is **shunt** for odd
//! `k` and **series** for even `k`. This is the conventional choice when the
//! prototype is drawn with `g1` as the first shunt element (Pozar fig. 8.25 /
//! Hong & Lancaster fig. 3.x): for the symmetric, equally-terminated all-pole
//! prototypes `yee-synth` produces (`g0 = g_{N+1} = 1` for odd `N`), the dual
//! "series-first" ladder is electrically equivalent, so shunt-first is a
//! well-defined and lossless representative. The choice is documented here
//! rather than made configurable because F2.0 is the band-pass-only walking
//! skeleton; a `series_first` flag is a trivial later addition if a topology
//! ever needs the dual.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use yee_synth::lowpass_to_bandpass;

use crate::{FilterProject, Response, SpecMask};

/// Which ladder branch an [`LcResonator`] sits on.
///
/// A `Series` resonator is a series L–C in the through arm; a `Shunt`
/// resonator is a parallel L–C to ground. The branches alternate along the
/// ladder starting **shunt-first** (see the [module docs](self)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LcBranch {
    /// Series L–C in the through (series) arm of the ladder.
    Series,
    /// Parallel L–C shunting the ladder node to ground.
    Shunt,
}

/// One L–C resonator of a lumped band-pass ladder, tuned to the centre `ω0`.
///
/// `l_henry` and `c_farad` satisfy `L·C = 1/ω0²` by construction (the resonator
/// is tuned to `f0`); see [`synthesize_lumped`] for the transform. The
/// [`branch`](LcResonator::branch) selects whether it is a series or shunt
/// resonator.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LcResonator {
    /// Which ladder branch this resonator sits on.
    pub branch: LcBranch,
    /// Resonator inductance, henries.
    pub l_henry: f64,
    /// Resonator capacitance, farads.
    pub c_farad: f64,
}

/// A synthesized lumped-element LC band-pass ladder.
///
/// Mirrors [`crate::EdgeCoupledDimensions`] on the distributed track: the
/// design centre / bandwidth / impedance plus the ordered list of
/// [`LcResonator`]s (one per reactive prototype element `g1..gN`, alternating
/// branch). Produced by [`synthesize_lumped`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LumpedLadder {
    /// Design centre frequency, Hz (`f0`).
    pub f0_hz: f64,
    /// Fractional bandwidth `Δ = (f2 − f1)/f0`.
    pub fbw: f64,
    /// System reference impedance, Ω (`Z0`).
    pub z0_ohm: f64,
    /// The `N` LC resonators in ladder order (first is **shunt**, then
    /// alternating — see the [module docs](self)).
    pub resonators: Vec<LcResonator>,
}

/// Errors from [`synthesize_lumped`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LumpedError {
    /// The project's [`Response`] is not [`Response::Bandpass`]; F2.0 only
    /// realizes the band-pass low-pass→band-pass transform.
    UnsupportedResponse,
    /// The filter order `N < 1`: there is no reactive element to realize.
    OrderTooSmall,
}

impl std::fmt::Display for LumpedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LumpedError::UnsupportedResponse => write!(
                f,
                "lumped LC ladder synthesis supports only Response::Bandpass"
            ),
            LumpedError::OrderTooSmall => {
                write!(f, "filter order N must be >= 1 to realize an LC resonator")
            }
        }
    }
}

impl std::error::Error for LumpedError {}

/// Synthesize an ideal lumped-element [`LumpedLadder`] from a synthesized
/// [`FilterProject`].
///
/// Applies the low-pass-prototype → band-pass ladder transform (Pozar §8.3) to
/// each reactive prototype element `g_k` (`k = 1..N`), producing one
/// [`LcResonator`] per element, alternating **shunt-first** (resonator `k` is
/// shunt for odd `k`, series for even `k` — see the [module docs](self)). Every
/// resonator is tuned to `ω0 = 2π·f0` (`L_k·C_k = 1/ω0²`).
///
/// Closed-form throughout: no optimizer, no FDTD. The `g0`/`g_{N+1}`
/// terminations of the prototype are the source/load `Z0` and are not
/// themselves resonators, so the ladder has exactly `N` resonators.
///
/// # Errors
///
/// - [`LumpedError::UnsupportedResponse`] if the project's response is not
///   [`Response::Bandpass`] (band-pass only for the F2.0 skeleton).
/// - [`LumpedError::OrderTooSmall`] if the order `N < 1`.
pub fn synthesize_lumped(project: &FilterProject) -> Result<LumpedLadder, LumpedError> {
    if project.spec.response != Response::Bandpass {
        return Err(LumpedError::UnsupportedResponse);
    }

    let n = project.prototype.order();
    if n < 1 {
        return Err(LumpedError::OrderTooSmall);
    }

    let f0 = project.spec.f0_hz;
    let fbw = project.spec.fbw;
    let z0 = project.spec.z0_ohm;
    let omega0 = std::f64::consts::TAU * f0;

    // Prototype g-vector is [g0, g1, …, gN, g_{N+1}]; the reactive elements are
    // g[1..=N]. Realize each as an alternating shunt/series LC resonator.
    let mut resonators = Vec::with_capacity(n);
    for k in 1..=n {
        let g_k = project.prototype.g[k];
        // Shunt-first: k odd → shunt, k even → series (see module docs).
        let branch = if k % 2 == 1 {
            LcBranch::Shunt
        } else {
            LcBranch::Series
        };
        let (l_henry, c_farad) = match branch {
            // Series L–C: L = g·Z0/(ω0·Δ), C = Δ/(ω0·Z0·g).
            LcBranch::Series => (g_k * z0 / (omega0 * fbw), fbw / (omega0 * z0 * g_k)),
            // Parallel (shunt) L–C: L = Z0·Δ/(ω0·g), C = g/(ω0·Z0·Δ).
            LcBranch::Shunt => (z0 * fbw / (omega0 * g_k), g_k / (omega0 * z0 * fbw)),
        };
        resonators.push(LcResonator {
            branch,
            l_henry,
            c_farad,
        });
    }

    Ok(LumpedLadder {
        f0_hz: f0,
        fbw,
        z0_ohm: z0,
        resonators,
    })
}

/// Forward transmission `S21` of the [`LumpedLadder`] at `f_hz`, by cascading
/// each resonator's ABCD matrix between a `Z0` source and a `Z0` load.
///
/// Each shunt resonator contributes a shunt admittance ABCD `[[1, 0], [Y, 1]]`
/// with `Y = jωC + 1/(jωL)`; each series resonator contributes a series
/// impedance ABCD `[[1, Z], [0, 1]]` with `Z = jωL + 1/(jωC)`. The cascade is
/// the ordered matrix product, and with equal `Z0` terminations
/// `S21 = 2 / (A + B/Z0 + C·Z0 + D)` (Pozar eq. 4.74).
///
/// This is an internal realized-response helper, **not** part of the documented
/// public API — it is `#[doc(hidden)] pub` solely so the
/// [`lumped_001`](../../tests/lumped_001.rs) integration gate (a separate crate)
/// can verify the LC realization reproduces the synthesized response.
#[doc(hidden)]
pub fn ladder_s21(ladder: &LumpedLadder, f_hz: f64) -> Complex64 {
    let z0 = Complex64::new(ladder.z0_ohm, 0.0);
    let omega = std::f64::consts::TAU * f_hz;
    let j = Complex64::new(0.0, 1.0);
    let jw = j * omega;

    // Start from the identity ABCD and right-multiply each resonator's matrix.
    let mut a = Complex64::new(1.0, 0.0);
    let mut b = Complex64::new(0.0, 0.0);
    let mut c = Complex64::new(0.0, 0.0);
    let mut d = Complex64::new(1.0, 0.0);

    for res in &ladder.resonators {
        let l = Complex64::new(res.l_henry, 0.0);
        let cap = Complex64::new(res.c_farad, 0.0);
        // Element ABCD [[ea, eb], [ec, ed]].
        let (ea, eb, ec, ed) = match res.branch {
            LcBranch::Series => {
                // Series impedance Z = jωL + 1/(jωC).
                let z = jw * l + Complex64::new(1.0, 0.0) / (jw * cap);
                (
                    Complex64::new(1.0, 0.0),
                    z,
                    Complex64::new(0.0, 0.0),
                    Complex64::new(1.0, 0.0),
                )
            }
            LcBranch::Shunt => {
                // Shunt admittance Y = jωC + 1/(jωL).
                let y = jw * cap + Complex64::new(1.0, 0.0) / (jw * l);
                (
                    Complex64::new(1.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    y,
                    Complex64::new(1.0, 0.0),
                )
            }
        };
        // [a b; c d] := [a b; c d] · [ea eb; ec ed].
        let na = a * ea + b * ec;
        let nb = a * eb + b * ed;
        let nc = c * ea + d * ec;
        let nd = c * eb + d * ed;
        a = na;
        b = nb;
        c = nc;
        d = nd;
    }

    Complex64::new(2.0, 0.0) / (a + b / z0 + c * z0 + d)
}

/// The realized-response spec-mask verdict for a [`LumpedLadder`].
///
/// Carries the strict pass/fail plus the worst-case graded quantities. Produced
/// by [`mask_verdict`], which evaluates the **realized** ladder response
/// ([`ladder_s21`]) — *not* the closed-form ideal response that
/// [`crate::check_mask`] grades. It is the shared verdict logic behind both the
/// `lumped_001` realization gate and the F2.4 [`crate::monte_carlo_yield`]
/// Monte-Carlo tolerance/yield analysis.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MaskVerdict {
    /// Overall verdict: `true` iff in-band ripple ≤ the mask's ripple bound,
    /// in-band return loss ≥ the mask's bound, and every stopband point is met.
    /// The ripple comparison uses `ripple_tol_db` slack (see [`mask_verdict`]).
    pub pass: bool,
    /// Worst-case in-band insertion-loss ripple observed, dB (`max_IL − min_IL`).
    pub worst_passband_ripple_db: f64,
    /// Worst-case (smallest) in-band return loss observed, dB.
    pub worst_return_loss_db: f64,
    /// Worst-case (smallest) stopband rejection across the mask's stopband
    /// points, dB; `+∞` when the mask has no stopband points.
    pub worst_stopband_rej_db: f64,
    /// Whether any swept frequency fell inside the passband (`|Ω| ≤ 1`). A
    /// verdict with no in-band sample cannot `pass`.
    pub saw_passband: bool,
}

/// Grade a realized [`LumpedLadder`]'s `S21` response against a [`SpecMask`].
///
/// This is the **shared** mask verdict that both the `lumped_001` realization
/// gate and the F2.4 [`crate::monte_carlo_yield`] yield analysis run on a
/// concrete (possibly tolerance-perturbed) ladder. It mirrors the closed-form
/// [`crate::check_mask`] comparison logic, but evaluates the *realized* response
/// [`ladder_s21`] rather than the ideal closed-form transfer function:
///
/// - **In-band** is `|Ω| ≤ 1` under the band-pass map (between the band edges).
///   Insertion loss is `−20·log10(|S21|)`; return loss is `−10·log10(|S11|²)`
///   with the lossless `|S11|² = 1 − |S21|²`. Ripple is `max_IL − min_IL` over
///   the in-band sweep; the worst (smallest) in-band return loss is tracked.
/// - **Stopband** points are graded at their own frequency: rejection
///   `−20·log10(|S21|)` must meet the required minimum.
///
/// `freqs_hz` should sample the passband densely enough to resolve the ripple
/// (the gates use ~800 points over `f0·(1 ± 3·fbw)`). `ripple_tol_db` is the
/// realization slack added to the ripple bound: pass `0.0` for a strict verdict
/// (the F2.4 yield path) or the small `lumped_001` realization margin that
/// absorbs the arithmetic-vs-geometric band-edge mismatch of the narrow-band LC
/// transform. The return-loss and stopband checks carry no slack.
pub fn mask_verdict(
    ladder: &LumpedLadder,
    mask: &SpecMask,
    f0_hz: f64,
    fbw: f64,
    freqs_hz: &[f64],
    ripple_tol_db: f64,
) -> MaskVerdict {
    let mut min_il = f64::INFINITY; // best (smallest) in-band insertion loss
    let mut max_il = f64::NEG_INFINITY; // worst (largest) in-band insertion loss
    let mut worst_rl = f64::INFINITY; // smallest in-band return loss
    let mut saw_passband = false;

    for &f in freqs_hz {
        if f <= 0.0 {
            continue;
        }
        let omega = lowpass_to_bandpass(f, f0_hz, fbw);
        if omega.abs() > 1.0 {
            continue; // out of band; graded by the stopband points instead
        }
        saw_passband = true;
        let s21_mag = ladder_s21(ladder, f).norm();
        let s11_sq = (1.0 - s21_mag * s21_mag).max(0.0);
        let il_db = -20.0 * s21_mag.max(1e-300).log10();
        let rl_db = if s11_sq <= 0.0 {
            f64::INFINITY
        } else {
            -10.0 * s11_sq.log10()
        };
        min_il = min_il.min(il_db);
        max_il = max_il.max(il_db);
        worst_rl = worst_rl.min(rl_db);
    }

    let ripple = if saw_passband {
        (max_il - min_il).max(0.0)
    } else {
        0.0
    };

    let mut pass = saw_passband;
    if saw_passband {
        if ripple > mask.passband_ripple_db + ripple_tol_db + 1e-9 {
            pass = false;
        }
        if worst_rl + 1e-9 < mask.return_loss_db {
            pass = false;
        }
    }

    let mut worst_rej = f64::INFINITY;
    for &(f_hz, required_db) in &mask.stopband {
        let s21_sq_mag = if f_hz <= 0.0 {
            0.0
        } else {
            ladder_s21(ladder, f_hz).norm()
        };
        let rejection_db = -20.0 * s21_sq_mag.max(1e-300).log10();
        worst_rej = worst_rej.min(rejection_db);
        if rejection_db + 1e-9 < required_db {
            pass = false;
        }
    }

    MaskVerdict {
        pass,
        worst_passband_ripple_db: ripple,
        worst_return_loss_db: if worst_rl.is_finite() { worst_rl } else { 0.0 },
        worst_stopband_rej_db: worst_rej,
        saw_passband,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn series_and_shunt_resonators_are_tuned() {
        // A trivial hand-built ladder: every resonator must satisfy L·C·ω0² = 1.
        let f0 = 2.0e9;
        let omega0 = std::f64::consts::TAU * f0;
        let ladder = LumpedLadder {
            f0_hz: f0,
            fbw: 0.1,
            z0_ohm: 50.0,
            resonators: vec![
                LcResonator {
                    branch: LcBranch::Shunt,
                    l_henry: 1.0 / (omega0 * omega0 * 1e-12),
                    c_farad: 1e-12,
                },
                LcResonator {
                    branch: LcBranch::Series,
                    l_henry: 1.0 / (omega0 * omega0 * 2e-12),
                    c_farad: 2e-12,
                },
            ],
        };
        for r in &ladder.resonators {
            let prod = r.l_henry * r.c_farad * omega0 * omega0;
            assert!((prod - 1.0).abs() < 1e-9, "resonator not tuned: {prod}");
        }
        // S21 at f0 is finite and well-formed.
        let s21 = ladder_s21(&ladder, f0);
        assert!(s21.norm().is_finite());
    }
}
