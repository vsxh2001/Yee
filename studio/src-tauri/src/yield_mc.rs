//! Yield-analysis flow for the studio (FS.5c, ADR-0222): the ADR-0211
//! closed-form patch-resonance testcase behind a form-shaped command.
//!
//! The request is design-centric — the user supplies the design frequency
//! f₀ and fabrication tolerances; the nominal patch length is derived from
//! the same closed form the `surrogate-yield-001` gate certified,
//! `f = c / (2 L √ε_eff)` with the zeroth-order `ε_eff = (ε_r + 1) / 2`.
//! Sampling and CI come straight from `yee_surrogate::yield_estimate`
//! (deterministic splitmix64 + Box-Muller, Wilson 95 % interval), so the
//! studio shows the exact numbers a gate would pin for the same seed.

use serde::{Deserialize, Serialize};
use yee_surrogate::ToleranceSpec;

/// Speed of light in vacuum, m/s.
const C0: f64 = 299_792_458.0;

/// Synchronous-command sample cap: the closed form at 10⁷ samples is
/// still sub-second, and the cap keeps a typo from freezing the shell.
const MAX_SAMPLES: usize = 10_000_000;

/// A yield-analysis request from the spec form.
#[derive(Debug, Clone, Deserialize)]
pub struct YieldRequest {
    /// Design (nominal resonance) frequency, Hz.
    pub f0_hz: f64,
    /// Nominal substrate relative permittivity (default FR-4, the
    /// ADR-0211 testcase value).
    #[serde(default = "default_eps_r")]
    pub eps_r: f64,
    /// Fabrication tolerance (1σ) on the patch length, metres.
    pub sigma_l_m: f64,
    /// Fabrication tolerance (1σ) on ε_r (batch spread).
    pub sigma_eps_r: f64,
    /// Spec half-width, Hz: a sample passes iff its resonance lands
    /// within ±this of `f0_hz`.
    pub spec_halfwidth_hz: f64,
    /// Number of Monte-Carlo samples.
    pub n_samples: usize,
    /// RNG seed — the estimate is bit-identical for identical inputs.
    pub seed: u64,
}

fn default_eps_r() -> f64 {
    4.4
}

/// The yield response: point estimate + explicit Wilson 95 % bounds.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct YieldResponse {
    /// Estimated yield, `n_pass / n_samples`.
    pub yield_frac: f64,
    /// Lower Wilson 95 % bound, clamped to [0, 1].
    pub ci95_lo: f64,
    /// Upper Wilson 95 % bound, clamped to [0, 1].
    pub ci95_hi: f64,
    /// Number of samples that met spec.
    pub n_pass: usize,
    /// Total number of samples drawn.
    pub n_samples: usize,
    /// Derived nominal patch length, metres — the dimension `sigma_l_m`
    /// perturbs, surfaced so the user sees what their tolerance applies to.
    pub length_nominal_m: f64,
}

/// Patch resonance closed form, Hz (ADR-0211 testcase).
fn patch_f_hz(l_m: f64, eps_r: f64) -> f64 {
    C0 / (2.0 * l_m * ((eps_r + 1.0) / 2.0).sqrt())
}

/// Pure yield flow: derive the nominal length, run the deterministic MC,
/// report yield + Wilson bounds. Instant at the default n = 10⁴.
pub fn yield_estimate_impl(req: &YieldRequest) -> Result<YieldResponse, String> {
    if !(req.f0_hz > 0.0) {
        return Err("f0 must be positive".into());
    }
    if !(req.eps_r > 1.0) {
        return Err("eps_r must be > 1".into());
    }
    if !(req.sigma_l_m >= 0.0 && req.sigma_eps_r >= 0.0) {
        return Err("tolerances (sigma) must be non-negative".into());
    }
    if !(req.spec_halfwidth_hz > 0.0) {
        return Err("spec half-width must be positive".into());
    }
    if req.n_samples == 0 || req.n_samples > MAX_SAMPLES {
        return Err(format!("n_samples must be in 1..={MAX_SAMPLES}"));
    }

    let eps_eff = (req.eps_r + 1.0) / 2.0;
    let l0 = C0 / (2.0 * req.f0_hz * eps_eff.sqrt());
    let spec = ToleranceSpec {
        nominal: vec![l0, req.eps_r],
        sigma: vec![req.sigma_l_m, req.sigma_eps_r],
    };
    let (f0, half) = (req.f0_hz, req.spec_halfwidth_hz);
    let est = yee_surrogate::yield_estimate(
        |p| {
            // A non-physical draw (possible only at absurd σ) fails spec
            // instead of feeding the closed form a division by ≤ 0.
            if p[0] <= 0.0 || p[1] <= -1.0 {
                return false;
            }
            (patch_f_hz(p[0], p[1]) - f0).abs() <= half
        },
        &spec,
        req.n_samples,
        req.seed,
    );
    Ok(YieldResponse {
        yield_frac: est.yield_frac,
        ci95_lo: (est.yield_frac - est.ci95_half_width).max(0.0),
        ci95_hi: (est.yield_frac + est.ci95_half_width).min(1.0),
        n_pass: est.n_pass,
        n_samples: est.n_samples,
        length_nominal_m: l0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ADR-0211 gate regime, phrased design-first: f₀ at the closed
    /// form's value for L = 29 mm / FR-4, gate seed, gate tolerances.
    fn adr_0211_request() -> YieldRequest {
        YieldRequest {
            f0_hz: patch_f_hz(29.0e-3, 4.4),
            eps_r: 4.4,
            sigma_l_m: 0.1e-3,
            sigma_eps_r: 0.05,
            spec_halfwidth_hz: 40.0e6,
            n_samples: 20_000,
            seed: 20260711,
        }
    }

    #[test]
    fn adr_0211_regime_yield_reproduced() {
        // surrogate-yield-001 measured brute-force 0.9721 at this seed.
        // The derived L₀ round-trips through f₀ (may differ from the
        // pinned 29 mm by an ULP), so assert the regime, not the bit.
        let resp = yield_estimate_impl(&adr_0211_request()).expect("valid request");
        assert!(
            resp.yield_frac > 0.95 && resp.yield_frac < 0.99,
            "yield {} outside the ADR-0211 regime",
            resp.yield_frac
        );
        assert!(resp.ci95_lo < resp.yield_frac && resp.yield_frac < resp.ci95_hi);
        assert!((resp.length_nominal_m - 29.0e-3).abs() < 1e-9);
        assert_eq!(resp.n_samples, 20_000);
    }

    #[test]
    fn deterministic_in_the_request() {
        let a = yield_estimate_impl(&adr_0211_request()).unwrap();
        let b = yield_estimate_impl(&adr_0211_request()).unwrap();
        assert_eq!(a, b, "same request must reproduce bit-identically");
    }

    #[test]
    fn rejects_invalid_inputs() {
        let ok = adr_0211_request();
        for (label, bad) in [
            ("f0", YieldRequest { f0_hz: 0.0, ..ok.clone() }),
            ("eps_r", YieldRequest { eps_r: 1.0, ..ok.clone() }),
            ("sigma_l", YieldRequest { sigma_l_m: -1e-6, ..ok.clone() }),
            ("halfwidth", YieldRequest { spec_halfwidth_hz: 0.0, ..ok.clone() }),
            ("n=0", YieldRequest { n_samples: 0, ..ok.clone() }),
            ("n>cap", YieldRequest { n_samples: MAX_SAMPLES + 1, ..ok.clone() }),
        ] {
            assert!(yield_estimate_impl(&bad).is_err(), "{label} must be rejected");
        }
    }
}
