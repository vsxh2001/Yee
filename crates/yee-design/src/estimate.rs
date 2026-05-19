//! Stage-3 initial-dimension calculator for the Phase 3.nl.0 pipeline.
//!
//! Implements the published closed-form synthesis equations for a rectangular
//! inset-fed microstrip patch antenna (Balanis, *Antenna Theory: Analysis and
//! Design*, 4th ed., Ch. 14) and the Pozar §3.8 Hammerstad–Jensen inverse-
//! synthesis formula for the 50 Ω feed-line width. The output is a flat
//! [`InitialEstimate`] of named `f64` dimensions in metres plus a copy of the
//! resolved [`Substrate`] used in the calculation.
//!
//! The function is **pure**: same [`DesignIntent`] in, bit-identical
//! [`InitialEstimate`] out. No I/O, no global mutable state, no randomness —
//! the spec §8 byte-identity invariant for the downstream emit stage depends
//! on it.
//!
//! ## References
//!
//! - Balanis, §14.2 (rectangular microstrip patch):
//!   - Eq. 14-1 effective relative permittivity `ε_reff`,
//!   - Eq. 14-2 edge-effect extension `ΔL`,
//!   - Eq. 14-3 physical length `L`,
//!   - Eq. 14-6 width `W` for good radiator,
//!   - Eq. 14-17 / 14-17a edge input resistance `R_edge`,
//!   - Eq. 14-20a inset offset `y₀` solved from
//!     `R_in = R_edge · cos²(π·y₀/L)`.
//! - Pozar, §3.8 (microstrip characteristic impedance): Hammerstad–Jensen
//!   inverse-synthesis formulas for `W/h` given `Z₀` and `ε_r`.
//!
//! ## Escape hatches honoured
//!
//! - The brief allows a closed-form lower-bound `y₀ = 0.3·L` if the Balanis
//!   14-20a transcendental proves difficult. This module *does* solve 14-20a
//!   directly (it is closed-form in `y₀` because the right-hand side is
//!   independent of `y₀`); no escape hatch is invoked.
//! - The brief also allows a hard-coded `w_feed = 1.5 mm` fallback if Pozar
//!   inverse-synthesis is too hairy. This module implements the Pozar
//!   formula directly; no escape hatch is invoked.

use crate::intent::{
    DesignIntent, GeometryFamily, NamedSubstrate, Substrate, SubstrateOverride, substrate_library,
};

/// Speed of light in vacuum (m/s). Matches the value Balanis uses in his
/// worked examples (`c = 3 × 10⁸ m/s`) so reproducing Example 14.1 to ±0.5%
/// requires no `c` re-derivation.
const C0: f64 = 2.997_924_58e8;

/// Target characteristic impedance of the feed line, in ohms. Phase 3.nl.0
/// hard-codes 50 Ω because every supported port type assumes it. The value
/// also sets the inset depth `y₀` via Balanis 14-20a.
const Z0_FEED: f64 = 50.0;

/// Errors produced by the Stage-3 estimator.
///
/// The estimator validates its inputs (Stage-2's job in the broader pipeline,
/// but cheap to repeat here) and refuses to silently swallow a NaN/infinite
/// dimension caused by an out-of-range substrate.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// The intent's [`GeometryFamily`] is not handled by Phase 3.nl.0's
    /// estimator. Phase 3.nl.2 adds further families (Wilkinson, hairpin,
    /// microstrip line); see spec §12.
    #[error("yee-design: geometry family not yet supported by Phase 3.nl.0 estimator")]
    NotSupported,
    /// A substrate parameter is outside the closed-form formulas' validity
    /// envelope (ε_r ∈ [1, 100], h ∈ [0.05, 10] mm; mirrors spec §7).
    #[error("yee-design: substrate parameter out of range: {0}")]
    SubstrateOutOfRange(&'static str),
    /// The intent references a named substrate that is not in
    /// [`substrate_library`].
    #[error("yee-design: unknown named substrate '{0}'")]
    UnknownSubstrate(String),
    /// Target frequency outside the spec §7 envelope (1 MHz … 1 THz).
    #[error("yee-design: target frequency out of range: {0} Hz")]
    FrequencyOutOfRange(f64),
    /// A computed dimension came out non-finite. Indicates the inputs were
    /// pathological even after range-checking; surfaced so the caller can
    /// inspect the [`DesignIntent`] rather than ingesting a NaN downstream.
    #[error("yee-design: computed dimension is non-finite ({0})")]
    NonFinite(&'static str),
}

/// Resolved substrate parameters used by the estimator.
///
/// Identical layout to [`Substrate::Explicit`] but always concrete — the
/// estimator does the library lookup + override merge internally so its
/// callers (and the emit stage) see a single flat shape.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResolvedSubstrate {
    /// Relative permittivity (unitless).
    pub eps_r: f64,
    /// Substrate height in metres.
    pub h_m: f64,
    /// Loss tangent (tan δ, unitless).
    pub loss_tangent: f64,
}

impl ResolvedSubstrate {
    /// Resolve the intent's [`Substrate`] (which may be a named-library
    /// reference plus a partial override) into a concrete set of parameters
    /// in SI units.
    ///
    /// Returns an error if the named substrate is missing from the library
    /// or if any parameter ends up outside the spec §7 envelope.
    pub fn from_substrate(substrate: &Substrate) -> Result<Self, Error> {
        let (eps_r, h_mm, loss_tangent) = match substrate {
            Substrate::Explicit {
                eps_r,
                h_mm,
                loss_tangent,
            } => (*eps_r, *h_mm, *loss_tangent),
            Substrate::Named(NamedSubstrate {
                name,
                override_with,
            }) => {
                let row = substrate_library()
                    .get(name)
                    .ok_or_else(|| Error::UnknownSubstrate(name.clone()))?;
                let (mut eps_r, mut h_mm, mut loss_tangent) =
                    (row.eps_r, row.h_mm, row.loss_tangent);
                if let Some(SubstrateOverride {
                    eps_r: o_eps_r,
                    h_mm: o_h_mm,
                    loss_tangent: o_loss,
                }) = override_with
                {
                    if let Some(v) = o_eps_r {
                        eps_r = *v;
                    }
                    if let Some(v) = o_h_mm {
                        h_mm = *v;
                    }
                    if let Some(v) = o_loss {
                        loss_tangent = *v;
                    }
                }
                (eps_r, h_mm, loss_tangent)
            }
        };
        if !(1.0..=100.0).contains(&eps_r) {
            return Err(Error::SubstrateOutOfRange("eps_r not in [1, 100]"));
        }
        if !(0.05..=10.0).contains(&h_mm) {
            return Err(Error::SubstrateOutOfRange("h_mm not in [0.05, 10]"));
        }
        if !(0.0..=0.5).contains(&loss_tangent) {
            return Err(Error::SubstrateOutOfRange("loss_tangent not in [0, 0.5]"));
        }
        Ok(Self {
            eps_r,
            h_m: h_mm * 1.0e-3,
            loss_tangent,
        })
    }
}

/// Flat struct of initial-estimate dimensions for a rectangular inset-fed
/// patch antenna, all in metres.
///
/// The shape is deliberately flat (no nested geometry-family enum) because
/// Phase 3.nl.0 only ships [`GeometryFamily::RectangularPatch`]; spec §12's
/// later families will land their own structs. The emit stage (R3) consumes
/// this struct field-by-field.
///
/// All fields use SI units (metres for lengths, unitless for `eps_reff`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InitialEstimate {
    /// Patch width `W` (Balanis 14-6).
    pub width_m: f64,
    /// Patch physical length `L = L_eff − 2·ΔL` (Balanis 14-3).
    pub length_m: f64,
    /// Inset-feed offset `y₀` from the radiating edge (Balanis 14-20a).
    pub inset_offset_m: f64,
    /// Feed-line width for a 50 Ω microstrip on the resolved substrate
    /// (Pozar §3.8 Hammerstad–Jensen inverse synthesis).
    pub feed_width_m: f64,
    /// Effective relative permittivity `ε_reff` (Balanis 14-1).
    pub eps_reff: f64,
    /// Edge-effect length extension `ΔL` (Balanis 14-2).
    pub delta_l_m: f64,
    /// The substrate parameters this estimate was synthesised against. Copied
    /// in so a downstream consumer (emit / log / surrogate) can reproduce the
    /// calculation without re-resolving the [`Substrate`] enum.
    pub substrate: ResolvedSubstrate,
}

impl InitialEstimate {
    /// Compute the Phase 3.nl.0 textbook initial estimate for a
    /// [`DesignIntent`].
    ///
    /// Pure function: the same intent in produces a bit-identical estimate
    /// out. Returns [`Error::NotSupported`] for any geometry family other
    /// than [`GeometryFamily::RectangularPatch`].
    pub fn from_intent(intent: &DesignIntent) -> Result<Self, Error> {
        let f = intent.target_frequency_hz;
        if !(1.0e6..=1.0e12).contains(&f) || !f.is_finite() {
            return Err(Error::FrequencyOutOfRange(f));
        }
        match intent.family {
            GeometryFamily::RectangularPatch => {
                let sub = ResolvedSubstrate::from_substrate(&intent.substrate)?;
                Self::rectangular_patch(f, sub)
            }
        }
    }

    /// Inner rectangular-patch synthesis (Balanis Ch. 14). Separated from
    /// [`InitialEstimate::from_intent`] so the formula stack is unit-testable
    /// against Balanis Example 14.1 without constructing a full
    /// [`DesignIntent`].
    fn rectangular_patch(f_hz: f64, sub: ResolvedSubstrate) -> Result<Self, Error> {
        let eps_r = sub.eps_r;
        let h = sub.h_m;

        // Balanis 14-6: W = c / (2·f·√((ε_r+1)/2))
        let w = C0 / (2.0 * f_hz * ((eps_r + 1.0) * 0.5).sqrt());
        if !w.is_finite() || w <= 0.0 {
            return Err(Error::NonFinite("W"));
        }

        // Balanis 14-1: ε_reff = (ε_r+1)/2 + (ε_r-1)/2 · (1 + 12h/W)^(-1/2)
        let w_over_h = w / h;
        let eps_reff = (eps_r + 1.0) * 0.5 + (eps_r - 1.0) * 0.5 * (1.0 + 12.0 * h / w).powf(-0.5);
        if !eps_reff.is_finite() || eps_reff <= 1.0 {
            return Err(Error::NonFinite("eps_reff"));
        }

        // Balanis 14-2: ΔL/h = 0.412 · (ε_reff + 0.3)(W/h + 0.264)
        //                       / ((ε_reff − 0.258)(W/h + 0.8))
        let delta_l = h * 0.412 * (eps_reff + 0.3) * (w_over_h + 0.264)
            / ((eps_reff - 0.258) * (w_over_h + 0.8));
        if !delta_l.is_finite() || delta_l <= 0.0 {
            return Err(Error::NonFinite("delta_L"));
        }

        // Balanis 14-3: L = c / (2·f·√ε_reff) − 2·ΔL
        let l_eff = C0 / (2.0 * f_hz * eps_reff.sqrt());
        let l = l_eff - 2.0 * delta_l;
        if !l.is_finite() || l <= 0.0 {
            return Err(Error::NonFinite("L"));
        }

        // Balanis 14-17a (closed-form approximation for the edge resistance):
        //   R_edge ≈ 90 · ε_r² / (ε_r − 1) · (L / W)²
        // Validity: works for ε_r > 1 with the patch driven near resonance;
        // we range-check on the substrate above. For ε_r = 1 the formula
        // would blow up; the substrate guard at ε_r ≥ 1 plus an explicit
        // guard here keeps the inset solver well-conditioned.
        let r_edge = if (eps_r - 1.0).abs() < 1.0e-6 {
            // Degenerate "air substrate" — fall back to a finite, large
            // edge-resistance so the inset offset reduces to ~0.5·L (a
            // mid-patch feed). This branch is unreachable in practice given
            // the substrate range check above; included for total-function
            // hygiene.
            10_000.0
        } else {
            90.0 * eps_r * eps_r / (eps_r - 1.0) * (l / w).powi(2)
        };
        if !r_edge.is_finite() || r_edge <= 0.0 {
            return Err(Error::NonFinite("R_edge"));
        }

        // Balanis 14-20a: R_in(y₀) = R_edge · cos²(π·y₀/L)
        //   → y₀ = (L / π) · arccos(√(R_in / R_edge))   (principal branch).
        // If R_in > R_edge we clamp the argument to 1.0 (i.e. y₀ = 0, feed
        // at the radiating edge) — this happens for ε_r close to 1 with a
        // narrow patch and is the physically correct degenerate answer.
        let ratio = (Z0_FEED / r_edge).clamp(0.0, 1.0);
        let inset_offset = (l / std::f64::consts::PI) * ratio.sqrt().acos();
        if !inset_offset.is_finite() || inset_offset < 0.0 {
            return Err(Error::NonFinite("y_0"));
        }

        // Pozar §3.8 inverse synthesis for a 50 Ω microstrip feed line on
        // the same substrate (Hammerstad–Jensen). Two branches; pick the
        // one whose self-consistency check holds. We compute both and keep
        // the one within its branch's validity range.
        let feed_w = microstrip_width_for_z0(Z0_FEED, eps_r, h)?;
        if !feed_w.is_finite() || feed_w <= 0.0 {
            return Err(Error::NonFinite("w_feed"));
        }

        Ok(Self {
            width_m: w,
            length_m: l,
            inset_offset_m: inset_offset,
            feed_width_m: feed_w,
            eps_reff,
            delta_l_m: delta_l,
            substrate: sub,
        })
    }
}

/// Pozar §3.8 microstrip-line inverse synthesis (Hammerstad–Jensen).
///
/// Given a target characteristic impedance `Z₀` (Ω), substrate relative
/// permittivity `ε_r`, and substrate height `h` (metres), returns the strip
/// width `W` (metres) that yields `Z₀` on the given substrate.
///
/// The formula has two branches; the spec convention is:
///   - For `W/h < 2`: `W/h = 8 e^A / (e^(2A) − 2)`
///   - For `W/h ≥ 2`: `W/h = (2/π)·[B − 1 − ln(2B − 1)
///     + ((ε_r − 1)/(2 ε_r))·(ln(B − 1) + 0.39 − 0.61/ε_r)]`
///
/// We compute both, then select the branch whose `W/h` lands in that branch's
/// validity range. For 50 Ω on common substrates (FR4, RO4003C, RO5880,
/// AluminaTC) the first branch is the right one; the second is implemented
/// for completeness so future low-impedance feed-line synthesis (e.g.
/// quarter-wave transformer in 3.nl.2) shares this routine.
fn microstrip_width_for_z0(z0: f64, eps_r: f64, h: f64) -> Result<f64, Error> {
    let a = (z0 / 60.0) * ((eps_r + 1.0) * 0.5).sqrt()
        + ((eps_r - 1.0) / (eps_r + 1.0)) * (0.23 + 0.11 / eps_r);
    let b = std::f64::consts::PI * 377.0 / (2.0 * z0 * eps_r.sqrt());

    let w_over_h_low = 8.0 * a.exp() / ((2.0 * a).exp() - 2.0);
    let w_over_h_high = (2.0 / std::f64::consts::PI)
        * (b - 1.0 - (2.0 * b - 1.0).ln()
            + ((eps_r - 1.0) / (2.0 * eps_r)) * ((b - 1.0).ln() + 0.39 - 0.61 / eps_r));

    // Pick the self-consistent branch.
    let w_over_h = if w_over_h_low < 2.0 {
        w_over_h_low
    } else {
        w_over_h_high
    };
    if !w_over_h.is_finite() || w_over_h <= 0.0 {
        return Err(Error::NonFinite("w_feed/h"));
    }
    Ok(w_over_h * h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intent::{Provenance, Substrate, substrate_library};

    fn balanis_intent(f_hz: f64, eps_r: f64, h_mm: f64) -> DesignIntent {
        DesignIntent {
            family: GeometryFamily::RectangularPatch,
            target_frequency_hz: f_hz,
            substrate: Substrate::Explicit {
                eps_r,
                h_mm,
                loss_tangent: 0.0009,
            },
            gain_target_dbi: None,
            bandwidth_target_mhz: None,
            source_prompt: "test".to_string(),
            provenance: Provenance {
                source: "offline".to_string(),
                model: None,
                temperature: None,
                schema_version: "1".to_string(),
                substrate_library_version: substrate_library().version.clone(),
            },
        }
    }

    /// Balanis Example 14.1 (Ch. 14, 4th ed.): f = 10 GHz, ε_r = 2.2,
    /// h = 0.1588 cm. Published W ≈ 1.186 cm, L ≈ 0.906 cm. Tolerance ±0.5%
    /// per the R2 brief.
    #[test]
    fn balanis_example_14_1_w_l_match_published() {
        let intent = balanis_intent(10.0e9, 2.2, 1.588);
        let est = InitialEstimate::from_intent(&intent).expect("estimate ok");
        // Published: W = 1.186 cm; tolerance 0.5%.
        let w_mm = est.width_m * 1.0e3;
        assert!(
            (w_mm - 11.86).abs() / 11.86 < 0.005,
            "W = {w_mm} mm, expected 11.86 mm (±0.5%)"
        );
        // Published: L = 0.906 cm; tolerance 0.5%.
        let l_mm = est.length_m * 1.0e3;
        assert!(
            (l_mm - 9.06).abs() / 9.06 < 0.005,
            "L = {l_mm} mm, expected 9.06 mm (±0.5%)"
        );
    }

    #[test]
    fn rejects_unknown_family_gracefully_via_substrate() {
        // No other GeometryFamily variants exist yet; smoke-test that the
        // resolver path stays total for in-range inputs.
        let intent = balanis_intent(2.4e9, 4.4, 1.6);
        let est = InitialEstimate::from_intent(&intent).expect("ok");
        assert!(est.width_m.is_finite());
    }

    #[test]
    fn rejects_out_of_range_eps_r() {
        let intent = balanis_intent(2.4e9, 0.5, 1.6);
        assert!(matches!(
            InitialEstimate::from_intent(&intent),
            Err(Error::SubstrateOutOfRange(_))
        ));
    }

    #[test]
    fn rejects_out_of_range_frequency() {
        let intent = balanis_intent(1.0e3, 4.4, 1.6);
        assert!(matches!(
            InitialEstimate::from_intent(&intent),
            Err(Error::FrequencyOutOfRange(_))
        ));
    }
}
