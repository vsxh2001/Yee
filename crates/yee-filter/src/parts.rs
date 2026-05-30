//! E-series component selection + **bill of materials** (Filter Phase F2.1).
//!
//! Maps the ideal `L`/`C` element values of a synthesized
//! [`LumpedLadder`](crate::LumpedLadder) (the just-shipped F2.0 output) onto
//! **real, purchasable standard components** — the nearest IEC 60063 preferred
//! ("E-series") value — and emits a grouped [`Bom`]. Pure `f64`/data + serde,
//! WASM-safe, NO FDTD, NO PCB footprints: this is the *component-choosing* brick
//! that turns ideal reactances into orderable parts and carries the per-part
//! tolerance that the later yield analysis (F2.4) consumes.
//!
//! It mirrors the [`crate::lumped`] module's shape (module-doc + serde structs +
//! `lib.rs` re-export) and consumes its [`LumpedLadder`](crate::LumpedLadder) /
//! [`LcResonator`](crate::LcResonator) types unchanged.
//!
//! # Method (IEC 60063)
//!
//! - **E-series:** the standard preferred values are a geometric series of
//!   `2 sig-fig` mantissas tiled across decades. [`ESeries::E24`] (24 values per
//!   decade, ±5 %) and [`ESeries::E96`] (96 values per decade, ±1 %) are
//!   supported.
//! - **Nearest-value selection** ([`ESeries::nearest`]): for an ideal value `x`,
//!   pick the E-series value minimizing `|log10(chosen) − log10(x)|`. Log-nearest
//!   is the correct metric for a geometric series — the decision boundary between
//!   two adjacent members sits at their *geometric* midpoint `√(mₖ·mₖ₊₁)`, not the
//!   arithmetic one. The decade is read from `x` and the member is searched over
//!   that decade and its two neighbours (so values near a decade boundary resolve
//!   correctly).
//! - **Deviation:** the recorded `deviation_pct = (chosen − ideal)/ideal·100` is
//!   bounded by half the series ratio step (E24 ⇒ ≲ 5 %, E96 ⇒ ≲ 1 %) by
//!   construction.
//! - **BOM** ([`select_components`]): one line per ladder element (an inductor
//!   for each resonator's `l_henry`, a capacitor for each `c_farad`),
//!   nearest-value-selected, then **identical `(kind, chosen_value)` lines are
//!   grouped** into a single line with a summed `qty`. The symmetric prototypes
//!   `yee-synth` produces (`g_k = g_{N+1−k}`) therefore collapse duplicate
//!   resonators into shared part numbers.
//!
//! # Parasitics
//!
//! [`BomLine`] carries OPTIONAL `esr_ohm`/`srf_hz` fields, defaulted `None` for
//! this skeleton; selection here is **value-only**. A real vendor-parts library
//! with measured ESR/SRF is a documented follow-on (F2.1b).

use serde::{Deserialize, Serialize};

use crate::{LcResonator, LumpedLadder};

/// IEC 60063 preferred-value ("E") series.
///
/// The mantissas are the standard 2-sig-fig preferred values for a decade; the
/// [`nearest`](ESeries::nearest) selector tiles them across decades. The two
/// variants carry the conventional component tolerance (see
/// [`tolerance_pct`](ESeries::tolerance_pct)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ESeries {
    /// 24 values per decade, ±5 % tolerance.
    E24,
    /// 96 values per decade, ±1 % tolerance.
    E96,
}

/// E24 mantissas (24 standard preferred values per decade, IEC 60063, ±5 %).
const E24_VALUES: [f64; 24] = [
    1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6,
    6.2, 6.8, 7.5, 8.2, 9.1,
];

/// E96 mantissas (96 standard preferred values per decade, IEC 60063, ±1 %).
///
/// These are the canonical published 2-sig-fig values, identical to
/// `round(10^(k/96)·100)/100` for `k = 0..96`.
const E96_VALUES: [f64; 96] = [
    1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30, 1.33, 1.37, 1.40, 1.43,
    1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74, 1.78, 1.82, 1.87, 1.91, 1.96, 2.00, 2.05, 2.10,
    2.15, 2.21, 2.26, 2.32, 2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09,
    3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12, 4.22, 4.32, 4.42, 4.53,
    4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49, 5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65,
    6.81, 6.98, 7.15, 7.32, 7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76,
];

impl ESeries {
    /// The standard preferred-value mantissas for one decade (`1.0 ≤ m < 10.0`).
    pub fn values_decade(&self) -> &'static [f64] {
        match self {
            ESeries::E24 => &E24_VALUES,
            ESeries::E96 => &E96_VALUES,
        }
    }

    /// Conventional component tolerance for the series, in percent.
    ///
    /// E24 → `5.0` (±5 %), E96 → `1.0` (±1 %). This is the per-part tolerance
    /// carried into the F2.4 yield analysis.
    pub fn tolerance_pct(&self) -> f64 {
        match self {
            ESeries::E24 => 5.0,
            ESeries::E96 => 1.0,
        }
    }

    /// The E-series value (`mantissa × 10^decade`) nearest to `x` in `log10`.
    ///
    /// Minimizes `|log10(chosen) − log10(x)|` — the correct distance for a
    /// geometric series, where adjacent members are separated by a constant
    /// ratio and the decision boundary is their *geometric* midpoint. The decade
    /// is read from `x` and every mantissa is tested in that decade and its two
    /// neighbours (`decade ± 1`), so a value just below a decade boundary can
    /// still snap up to `1.0 × 10^(decade+1)` when that is log-closer.
    ///
    /// `x` must be finite and strictly positive (component values are positive
    /// reactances); for `x ≤ 0` or non-finite `x` the result is `x` unchanged.
    pub fn nearest(&self, x: f64) -> f64 {
        if !x.is_finite() || x <= 0.0 {
            return x;
        }
        let log_x = x.log10();
        let decade = log_x.floor() as i32;
        let mantissas = self.values_decade();

        let mut best = x;
        let mut best_dist = f64::INFINITY;
        for d in [decade - 1, decade, decade + 1] {
            let scale = 10f64.powi(d);
            for &m in mantissas {
                let candidate = m * scale;
                let dist = (candidate.log10() - log_x).abs();
                if dist < best_dist {
                    best_dist = dist;
                    best = candidate;
                }
            }
        }
        best
    }
}

/// Which kind of lumped component a [`BomLine`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CompKind {
    /// An inductor (realizing a resonator's `l_henry`), value in henries.
    Inductor,
    /// A capacitor (realizing a resonator's `c_farad`), value in farads.
    Capacitor,
}

/// One bill-of-materials line: a standard-value part plus its quantity.
///
/// Produced by [`select_components`] from each ladder element. The `ideal_value`
/// is the synthesized reactance; `chosen_value` is the nearest E-series member;
/// `deviation_pct = (chosen − ideal)/ideal·100`. `tolerance_pct` is the series
/// tolerance (carried into the F2.4 yield analysis). `esr_ohm`/`srf_hz` are
/// optional parasitics, `None` in this value-only skeleton (a real vendor-parts
/// library is the F2.1b follow-on). `qty` is the number of identical
/// `(kind, chosen_value)` parts grouped into this line.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BomLine {
    /// Inductor or capacitor.
    pub kind: CompKind,
    /// The ideal synthesized value (H for an inductor, F for a capacitor).
    pub ideal_value: f64,
    /// The nearest E-series value chosen to realize [`ideal_value`](Self::ideal_value).
    pub chosen_value: f64,
    /// Signed deviation of the chosen value from ideal, percent:
    /// `(chosen − ideal)/ideal·100`.
    pub deviation_pct: f64,
    /// The E-series the value was chosen from.
    pub series: ESeries,
    /// The series tolerance, percent (e.g. `5.0` for E24). Feeds F2.4.
    pub tolerance_pct: f64,
    /// Quantity: number of identical `(kind, chosen_value)` parts on this line.
    pub qty: usize,
    /// Optional equivalent series resistance, Ω (`None` in this skeleton).
    pub esr_ohm: Option<f64>,
    /// Optional self-resonant frequency, Hz (`None` in this skeleton).
    pub srf_hz: Option<f64>,
}

/// A grouped bill of materials: standard-value parts + quantities.
///
/// Produced by [`select_components`]. Identical `(kind, chosen_value)` parts are
/// merged into a single [`BomLine`] with a summed [`qty`](BomLine::qty); see
/// [`total_parts`](Bom::total_parts) for the physical part count.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bom {
    /// The grouped BOM lines (one per distinct `(kind, chosen_value)` part).
    pub lines: Vec<BomLine>,
}

impl Bom {
    /// Total number of physical parts: the sum of every line's
    /// [`qty`](BomLine::qty).
    pub fn total_parts(&self) -> usize {
        self.lines.iter().map(|l| l.qty).sum()
    }
}

/// Select standard E-series components for a [`LumpedLadder`] and emit a [`Bom`].
///
/// For each [`LcResonator`] in the ladder, emit an inductor line (for
/// `l_henry`) and a capacitor line (for `c_farad`), each snapped to the nearest
/// `series` value (log-nearest, [`ESeries::nearest`]) with the recorded signed
/// `deviation_pct` and the series `tolerance_pct`; `esr_ohm`/`srf_hz` are left
/// `None`. Identical `(kind, chosen_value)` lines are then **grouped** into a
/// single line with a summed `qty`, so a symmetric ladder (whose duplicate
/// resonators select identical parts) collapses to shared part numbers. Lines
/// are kept in first-encountered order (inductor-then-capacitor, resonator by
/// resonator).
pub fn select_components(ladder: &LumpedLadder, series: ESeries) -> Bom {
    let tolerance_pct = series.tolerance_pct();

    // Build the ungrouped line for one (kind, ideal-value) pair.
    let make_line = |kind: CompKind, ideal_value: f64| -> BomLine {
        let chosen_value = series.nearest(ideal_value);
        let deviation_pct = if ideal_value != 0.0 {
            (chosen_value - ideal_value) / ideal_value * 100.0
        } else {
            0.0
        };
        BomLine {
            kind,
            ideal_value,
            chosen_value,
            deviation_pct,
            series,
            tolerance_pct,
            qty: 1,
            esr_ohm: None,
            srf_hz: None,
        }
    };

    // Emit L then C for each resonator, grouping identical (kind, chosen) parts.
    let mut lines: Vec<BomLine> = Vec::new();
    for &LcResonator {
        l_henry, c_farad, ..
    } in &ladder.resonators
    {
        for line in [
            make_line(CompKind::Inductor, l_henry),
            make_line(CompKind::Capacitor, c_farad),
        ] {
            match lines
                .iter_mut()
                .find(|l| l.kind == line.kind && l.chosen_value == line.chosen_value)
            {
                Some(existing) => existing.qty += 1,
                None => lines.push(line),
            }
        }
    }

    Bom { lines }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearest_is_always_a_series_member() {
        for series in [ESeries::E24, ESeries::E96] {
            for &x in &[1.0, 4.5e-9, 1.04e3, 6.6e-12, 9.5, 3.3e-6, 7.0] {
                let chosen = series.nearest(x);
                let mantissa = chosen / 10f64.powi(chosen.log10().floor() as i32);
                let member = series
                    .values_decade()
                    .iter()
                    .any(|&m| (m - mantissa).abs() < 1e-6);
                assert!(member, "{chosen} (mantissa {mantissa}) not in {series:?}");
            }
        }
    }

    #[test]
    fn tolerances() {
        assert_eq!(ESeries::E24.tolerance_pct(), 5.0);
        assert_eq!(ESeries::E96.tolerance_pct(), 1.0);
    }
}
