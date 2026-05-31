//! Guided technique-recommender (App.2.0, ADR-0136).
//!
//! A **deterministic decision tree** that turns a [`FilterSpec`] into a
//! recommended physical [`RealizationTechnique`], a plain-language **rationale**
//! naming the deciding factor, and a ranked list of **alternatives**. This is
//! the validatable, pure-domain core of the studio's guided "recommend-a-
//! technique" entry (the Nuhertz FilterQuick pattern) — the UI in
//! `yee-studio-web` is a thin consumer of [`recommend_technique`].
//!
//! The tree is keyed on the response class, the centre / cutoff frequency
//! (`f0_hz`), and the fractional bandwidth (`fbw`). The thresholds — the
//! ≈500 MHz distributed-feasibility floor and the 5 % / 20 % fractional-
//! bandwidth bands — are documented engineering judgment (Pozar, *Microwave
//! Engineering* Ch. 8; Hong & Lancaster; Matthaei, Young & Jones) **pinned by a
//! gate** (`tech_001`) so they cannot silently drift.
//!
//! The engine stays pure-domain: it has no notion of which techniques the
//! studio can currently *build* (the "live" vs "Soon" distinction is a UI
//! concern, kept in `yee-studio-web`).
//!
//! ## Example
//!
//! ```
//! use yee_filter::{recommend_technique, RealizationTechnique};
//! # use yee_filter::{Approximation, FilterSpec, Response, SpecMask};
//! # let spec = FilterSpec {
//! #     response: Response::Bandpass,
//! #     approximation: Approximation::Chebyshev { ripple_db: 0.5 },
//! #     f0_hz: 2.4e9,
//! #     fbw: 0.05,
//! #     order: Some(5),
//! #     z0_ohm: 50.0,
//! #     mask: SpecMask { passband_ripple_db: 0.5, return_loss_db: 10.0, stopband: vec![] },
//! # };
//! let rec = recommend_technique(&spec);
//! assert_eq!(rec.primary, RealizationTechnique::EdgeCoupled);
//! assert!(!rec.rationale.is_empty());
//! ```

use crate::{FilterSpec, Response};

/// The distributed-feasibility frequency floor, Hz.
///
/// Below this centre / cutoff frequency a quarter-wave distributed resonator is
/// physically large (λ/4 ≈ 15 cm at 500 MHz on FR-4), so a lumped-element
/// realization is preferred; at or above it distributed microstrip techniques
/// become practical (Pozar Ch. 8).
const DISTRIBUTED_FLOOR_HZ: f64 = 500e6;

/// The wide-bandwidth fractional-bandwidth threshold.
///
/// At or above `fbw = 0.20` the band is wide enough that parallel-coupled
/// (edge-coupled) resonators are the natural distributed choice (Hong &
/// Lancaster).
const WIDE_FBW: f64 = 0.20;

/// The narrow-bandwidth fractional-bandwidth threshold.
///
/// Below `fbw = 0.05` edge-coupled gaps become impractically tight, so a
/// compact high-Q technique (interdigital / combline) is preferred over
/// parallel-coupled lines (Matthaei, Young & Jones).
const NARROW_FBW: f64 = 0.05;

/// A physical realization technique the studio can target.
///
/// Each variant maps 1:1 to a Technique-stage gallery topology in
/// `yee-studio-web`; whether a given technique is currently *buildable* (live)
/// or roadmapped (Soon) is a UI concern, not part of this pure-domain engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RealizationTechnique {
    /// Discrete lumped-element L/C ladder (SMD parts + BOM + tolerance).
    LumpedLc,
    /// Distributed edge-coupled (parallel-coupled) half-wave microstrip lines.
    EdgeCoupled,
    /// U-folded half-wave resonators — a more compact edge-coupled variant.
    Hairpin,
    /// Grounded quarter-wave resonators with end-loading capacitors.
    Combline,
    /// Interleaved grounded quarter-wave fingers — compact and high-Q.
    Interdigital,
    /// Distributed lowpass: alternating high-/low-impedance line sections.
    SteppedImpedance,
}

impl RealizationTechnique {
    /// A short human-readable name for display.
    pub fn name(self) -> &'static str {
        match self {
            RealizationTechnique::LumpedLc => "Lumped LC",
            RealizationTechnique::EdgeCoupled => "Edge-coupled",
            RealizationTechnique::Hairpin => "Hairpin",
            RealizationTechnique::Combline => "Combline",
            RealizationTechnique::Interdigital => "Interdigital",
            RealizationTechnique::SteppedImpedance => "Stepped-impedance",
        }
    }
}

/// A technique recommendation: a primary choice + rationale + ranked
/// alternatives.
///
/// Produced by [`recommend_technique`]. The `primary` is the engine's first
/// choice; `alternatives` is a ranked list of `(technique, one-line tradeoff)`
/// pairs (best first). The `primary` never appears in `alternatives`.
#[derive(Debug, Clone, PartialEq)]
pub struct TechniqueRecommendation {
    /// The recommended realization technique.
    pub primary: RealizationTechnique,
    /// Plain-language rationale that **names the deciding factor** (the response
    /// class, the frequency vs the distributed floor, or the fractional
    /// bandwidth band).
    pub rationale: String,
    /// Ranked alternative techniques, each with a one-line tradeoff note (best
    /// alternative first). Does **not** contain [`primary`](Self::primary).
    pub alternatives: Vec<(RealizationTechnique, String)>,
}

/// Recommend a realization technique from a filter [`FilterSpec`].
///
/// A pure, deterministic decision tree keyed on `spec.response`, `spec.f0_hz`
/// (centre frequency for band filters, cutoff for low/high-pass), and
/// `spec.fbw`. Returns a [`TechniqueRecommendation`] whose `rationale` names the
/// deciding factor and whose `alternatives` are ranked best-first (and never
/// include the primary). Never panics on a valid spec.
///
/// The tree (thresholds pinned by the `tech_001` gate):
///
/// - **Lowpass:** cutoff ≥ 500 MHz → [`SteppedImpedance`](RealizationTechnique::SteppedImpedance)
///   (distributed hi/lo-Z sections); else [`LumpedLc`](RealizationTechnique::LumpedLc)
///   (distributed sections impractically long below ~500 MHz).
/// - **Highpass:** [`LumpedLc`](RealizationTechnique::LumpedLc) (distributed
///   microstrip high-pass is impractical; honest — Yee's distributed techniques
///   are LP/BP-oriented).
/// - **Bandpass / Bandstop:**
///   - `f0` < 500 MHz → [`LumpedLc`](RealizationTechnique::LumpedLc)
///     (distributed resonators too large).
///   - `f0` ≥ 500 MHz: `fbw` ≥ 0.20 → [`EdgeCoupled`](RealizationTechnique::EdgeCoupled);
///     0.05 ≤ `fbw` < 0.20 → [`EdgeCoupled`](RealizationTechnique::EdgeCoupled)
///     primary with a [`Hairpin`](RealizationTechnique::Hairpin) alternative;
///     `fbw` < 0.05 → [`Interdigital`](RealizationTechnique::Interdigital)
///     primary with a [`Combline`](RealizationTechnique::Combline) alternative.
pub fn recommend_technique(spec: &FilterSpec) -> TechniqueRecommendation {
    use RealizationTechnique::*;

    let f0_ghz = spec.f0_hz / 1e9;

    match spec.response {
        Response::Lowpass => {
            if spec.f0_hz >= DISTRIBUTED_FLOOR_HZ {
                TechniqueRecommendation {
                    primary: SteppedImpedance,
                    rationale: format!(
                        "Lowpass at a {f0_ghz:.3} GHz cutoff (≥ 500 MHz): a distributed \
                         stepped-impedance microstrip lowpass (alternating high-/low-Z \
                         sections) is compact and practical here."
                    ),
                    alternatives: vec![(
                        LumpedLc,
                        "Discrete L/C ladder — simpler to tune, but larger parts and \
                         parasitics climb with frequency."
                            .to_string(),
                    )],
                }
            } else {
                TechniqueRecommendation {
                    primary: LumpedLc,
                    rationale: format!(
                        "Lowpass at a {:.1} MHz cutoff (< 500 MHz): distributed line \
                         sections become impractically long, so a discrete lumped L/C \
                         ladder is the natural choice.",
                        spec.f0_hz / 1e6
                    ),
                    alternatives: vec![(
                        SteppedImpedance,
                        "Distributed hi/lo-Z lowpass — viable only once the cutoff is in \
                         the GHz range."
                            .to_string(),
                    )],
                }
            }
        }
        Response::Highpass => TechniqueRecommendation {
            primary: LumpedLc,
            rationale: "Highpass response: a distributed microstrip high-pass is \
                 impractical (Yee's distributed techniques are lowpass/bandpass-oriented), \
                 so a discrete lumped L/C realization is recommended."
                .to_string(),
            alternatives: Vec::new(),
        },
        Response::Bandpass | Response::Bandstop => {
            let stop_note = if spec.response == Response::Bandstop {
                " A distributed band-stop (open-stub) technique is a future option."
            } else {
                ""
            };
            if spec.f0_hz < DISTRIBUTED_FLOOR_HZ {
                TechniqueRecommendation {
                    primary: LumpedLc,
                    rationale: format!(
                        "Band filter centred at {:.1} MHz (< 500 MHz): a quarter-wave \
                         distributed resonator (λ/4 ≈ 15 cm at 500 MHz) would be too large, \
                         so a discrete lumped L/C ladder is recommended.{stop_note}",
                        spec.f0_hz / 1e6
                    ),
                    alternatives: vec![(
                        EdgeCoupled,
                        "Distributed parallel-coupled lines — possible, but the board \
                         becomes large at this frequency."
                            .to_string(),
                    )],
                }
            } else if spec.fbw >= WIDE_FBW {
                TechniqueRecommendation {
                    primary: EdgeCoupled,
                    rationale: format!(
                        "Band filter at {f0_ghz:.3} GHz with a wide {:.0}% fractional \
                         bandwidth (≥ 20%): parallel-coupled (edge-coupled) microstrip \
                         lines handle wide bandwidths comfortably.{stop_note}",
                        spec.fbw * 100.0
                    ),
                    alternatives: vec![(
                        LumpedLc,
                        "Discrete L/C ladder — compact, but high-frequency parts get lossy \
                         and parasitic-limited."
                            .to_string(),
                    )],
                }
            } else if spec.fbw >= NARROW_FBW {
                TechniqueRecommendation {
                    primary: EdgeCoupled,
                    rationale: format!(
                        "Band filter at {f0_ghz:.3} GHz with a moderate {:.1}% fractional \
                         bandwidth (5–20%): parallel-coupled (edge-coupled) microstrip \
                         lines are the standard distributed choice.{stop_note}",
                        spec.fbw * 100.0
                    ),
                    alternatives: vec![(
                        Hairpin,
                        "U-folds the same half-wave resonators for a smaller board, at the \
                         cost of slightly more coupling-design care."
                            .to_string(),
                    )],
                }
            } else {
                TechniqueRecommendation {
                    primary: Interdigital,
                    rationale: format!(
                        "Band filter at {f0_ghz:.3} GHz with a narrow {:.1}% fractional \
                         bandwidth (< 5%): edge-coupled gaps would become impractically \
                         tight, so a compact, high-Q interdigital filter (quarter-wave \
                         coupled lines) is recommended.{stop_note}",
                        spec.fbw * 100.0
                    ),
                    alternatives: vec![(
                        Combline,
                        "Capacitively end-loaded quarter-wave resonators — even more \
                         compact and tunable."
                            .to_string(),
                    )],
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Approximation, SpecMask};

    /// Build a minimal spec for the given response / centre-or-cutoff frequency
    /// / fractional bandwidth (the only fields the recommender reads).
    fn spec(response: Response, f0_hz: f64, fbw: f64) -> FilterSpec {
        FilterSpec {
            response,
            approximation: Approximation::Chebyshev { ripple_db: 0.5 },
            f0_hz,
            fbw,
            order: Some(5),
            z0_ohm: 50.0,
            mask: SpecMask {
                passband_ripple_db: 0.5,
                return_loss_db: 10.0,
                stopband: vec![],
            },
        }
    }

    #[test]
    fn primary_is_never_in_alternatives() {
        // Cover every leaf of the tree.
        let cases = [
            spec(Response::Lowpass, 1e9, 0.0),
            spec(Response::Lowpass, 50e6, 0.0),
            spec(Response::Highpass, 1e9, 0.0),
            spec(Response::Bandpass, 100e6, 0.05),
            spec(Response::Bandpass, 2.4e9, 0.05),
            spec(Response::Bandpass, 2.4e9, 0.25),
            spec(Response::Bandpass, 5e9, 0.02),
            spec(Response::Bandstop, 2.4e9, 0.02),
        ];
        for s in &cases {
            let rec = recommend_technique(s);
            assert!(
                !rec.rationale.trim().is_empty(),
                "rationale must be non-empty for {:?}",
                s.response
            );
            assert!(
                !rec.alternatives.iter().any(|(t, _)| *t == rec.primary),
                "primary {:?} must not appear in alternatives for {:?}",
                rec.primary,
                s.response
            );
            for (_, note) in &rec.alternatives {
                assert!(
                    !note.trim().is_empty(),
                    "alternative note must be non-empty"
                );
            }
        }
    }

    #[test]
    fn distributed_floor_is_inclusive() {
        // Exactly 500 MHz lowpass is distributed (≥ floor).
        assert_eq!(
            recommend_technique(&spec(Response::Lowpass, 500e6, 0.0)).primary,
            RealizationTechnique::SteppedImpedance
        );
        // Just below stays lumped.
        assert_eq!(
            recommend_technique(&spec(Response::Lowpass, 499e6, 0.0)).primary,
            RealizationTechnique::LumpedLc
        );
    }

    #[test]
    fn bandpass_fbw_bands() {
        // 20% boundary is inclusive → EdgeCoupled (wide).
        let wide = recommend_technique(&spec(Response::Bandpass, 2.4e9, 0.20));
        assert_eq!(wide.primary, RealizationTechnique::EdgeCoupled);
        // 5% boundary is inclusive → EdgeCoupled (moderate), Hairpin alt.
        let mid = recommend_technique(&spec(Response::Bandpass, 2.4e9, 0.05));
        assert_eq!(mid.primary, RealizationTechnique::EdgeCoupled);
        assert_eq!(mid.alternatives[0].0, RealizationTechnique::Hairpin);
        // Below 5% → Interdigital, Combline alt.
        let narrow = recommend_technique(&spec(Response::Bandpass, 2.4e9, 0.04));
        assert_eq!(narrow.primary, RealizationTechnique::Interdigital);
        assert_eq!(narrow.alternatives[0].0, RealizationTechnique::Combline);
    }
}
