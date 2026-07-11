//! Monte-Carlo yield analysis over Gaussian dimension tolerances
//! (FS.5a, ADR-0211).
//!
//! "Yield" is the fraction of manufactured devices meeting spec when each
//! design parameter is perturbed by its fabrication tolerance. The
//! estimator is model-agnostic: the pass/fail closure may evaluate a
//! closed form, a trained [`crate::GaussianProcess`], or a full solver —
//! the estimator only draws parameter samples and counts passes.
//!
//! Randomness is a deterministic in-crate `splitmix64` stream (no `rand`
//! dependency): the same seed reproduces the same estimate bit-for-bit on
//! every platform, which makes yield numbers pinnable in gates and
//! regression-stable in CI (gate `yield-mc-002`).

/// Independent Gaussian tolerance on each design parameter.
///
/// `params[i] ~ N(nominal[i], sigma[i]²)`. This is the universal PCB-fab
/// datum (etch width, substrate thickness, ε_r spread are quoted as ±3σ
/// or ±tolerance bands); correlated or non-Gaussian tolerances are a
/// follow-on behind the same API.
#[derive(Debug, Clone)]
pub struct ToleranceSpec {
    /// Nominal (as-designed) parameter values.
    pub nominal: Vec<f64>,
    /// Standard deviation of each parameter's fabrication tolerance.
    /// Must have the same length as `nominal`; a `sigma` of 0 pins that
    /// parameter to its nominal.
    pub sigma: Vec<f64>,
}

/// Result of a Monte-Carlo yield run.
#[derive(Debug, Clone, PartialEq)]
pub struct YieldEstimate {
    /// Estimated yield: `n_pass / n_samples`.
    pub yield_frac: f64,
    /// Half-width of the Wilson 95 % score interval around `yield_frac`.
    /// (Wilson, not Wald: the Wald interval collapses to zero width at
    /// yield → 0 or 1, exactly where yield analysis operates.)
    pub ci95_half_width: f64,
    /// Number of samples that met spec.
    pub n_pass: usize,
    /// Total number of Monte-Carlo samples drawn.
    pub n_samples: usize,
}

/// Estimate yield by Monte-Carlo: draw `n_samples` parameter vectors from
/// the tolerance distribution and count how many satisfy `pass`.
///
/// Deterministic in `seed` — identical inputs give a bit-identical
/// [`YieldEstimate`]. Panics if `spec.nominal` and `spec.sigma` lengths
/// differ or `n_samples` is 0 (both are caller bugs, not data).
pub fn yield_estimate(
    mut pass: impl FnMut(&[f64]) -> bool,
    spec: &ToleranceSpec,
    n_samples: usize,
    seed: u64,
) -> YieldEstimate {
    assert_eq!(
        spec.nominal.len(),
        spec.sigma.len(),
        "nominal and sigma must have the same length"
    );
    assert!(n_samples > 0, "n_samples must be positive");
    let mut rng = SplitMix64::new(seed);
    let mut normals = BoxMuller::new();
    let mut params = vec![0.0_f64; spec.nominal.len()];
    let mut n_pass = 0usize;
    for _ in 0..n_samples {
        for (p, (&mu, &s)) in params.iter_mut().zip(spec.nominal.iter().zip(&spec.sigma)) {
            *p = mu + s * normals.next(&mut rng);
        }
        if pass(&params) {
            n_pass += 1;
        }
    }
    let (yield_frac, ci95_half_width) = wilson_95(n_pass, n_samples);
    YieldEstimate {
        yield_frac,
        ci95_half_width,
        n_pass,
        n_samples,
    }
}

/// Wilson 95 % score interval: returns (point estimate p̂, half-width of
/// the interval **around p̂**, i.e. max distance from p̂ to either Wilson
/// bound) so `p̂ ± half_width` always covers the Wilson interval.
fn wilson_95(n_pass: usize, n: usize) -> (f64, f64) {
    const Z: f64 = 1.959_963_984_540_054; // Φ⁻¹(0.975)
    let nf = n as f64;
    let p = n_pass as f64 / nf;
    let z2 = Z * Z;
    let denom = 1.0 + z2 / nf;
    let centre = (p + z2 / (2.0 * nf)) / denom;
    let half = (Z / denom) * (p * (1.0 - p) / nf + z2 / (4.0 * nf * nf)).sqrt();
    let lo = (centre - half).max(0.0);
    let hi = (centre + half).min(1.0);
    (p, (p - lo).max(hi - p))
}

/// `splitmix64` — the public-domain 64-bit mixer (Steele/Lea/Flood 2014;
/// Vigna's reference implementation). Tiny, passes BigCrush as a stream
/// generator, and — crucially for gates — has no platform-dependent
/// state: pure wrapping integer arithmetic.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in the open interval (0, 1): 53-bit mantissa, then the
    /// zero (probability 2⁻⁵³) is nudged to the smallest step so
    /// `ln(u)` in Box-Muller never sees 0.
    fn next_open01(&mut self) -> f64 {
        let u = ((self.next_u64() >> 11) as f64) * (1.0 / (1u64 << 53) as f64);
        if u == 0.0 {
            0.5 / (1u64 << 53) as f64
        } else {
            u
        }
    }
}

/// Box-Muller transform: pairs of uniforms → pairs of independent
/// standard normals, buffering the second of each pair.
struct BoxMuller {
    spare: Option<f64>,
}

impl BoxMuller {
    fn new() -> Self {
        Self { spare: None }
    }

    fn next(&mut self, rng: &mut SplitMix64) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        let u1 = rng.next_open01();
        let u2 = rng.next_open01();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = std::f64::consts::TAU * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_matches_reference_vector() {
        // Vigna's reference splitmix64 with seed 1234567:
        // first three outputs (checked against the C reference).
        let mut rng = SplitMix64::new(1234567);
        assert_eq!(rng.next_u64(), 6457827717110365317);
        assert_eq!(rng.next_u64(), 3203168211198807973);
        assert_eq!(rng.next_u64(), 9817491932198370423);
    }

    #[test]
    fn normal_stream_has_sane_moments() {
        let mut rng = SplitMix64::new(42);
        let mut bm = BoxMuller::new();
        let n = 200_000;
        let (mut sum, mut sum2) = (0.0, 0.0);
        for _ in 0..n {
            let z = bm.next(&mut rng);
            sum += z;
            sum2 += z * z;
        }
        let mean = sum / n as f64;
        let var = sum2 / n as f64 - mean * mean;
        // 3σ bounds for n = 2e5: mean ±0.0067, var ±0.0095.
        assert!(mean.abs() < 0.01, "mean {mean}");
        assert!((var - 1.0).abs() < 0.015, "var {var}");
    }

    #[test]
    fn wilson_interval_stays_in_unit_range_at_endpoints() {
        for &(n_pass, n) in &[(0usize, 1000usize), (1000, 1000), (1, 1000), (999, 1000)] {
            let (p, hw) = wilson_95(n_pass, n);
            assert!((0.0..=1.0).contains(&(p - hw).max(0.0)));
            assert!((0.0..=1.0).contains(&(p + hw).min(1.0)));
            assert!(hw > 0.0, "Wilson half-width must not collapse at p={p}");
        }
    }

    #[test]
    fn zero_sigma_pins_parameter_to_nominal() {
        let spec = ToleranceSpec {
            nominal: vec![3.0, 5.0],
            sigma: vec![0.0, 1.0],
        };
        let est = yield_estimate(|p| p[0] == 3.0, &spec, 100, 7);
        assert_eq!(est.n_pass, 100);
    }
}
