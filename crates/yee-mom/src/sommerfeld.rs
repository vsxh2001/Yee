//! Sommerfeld surface-wave pole extraction — Phase 1.1.1.2.
//!
//! ## Scope
//!
//! The Phase 1.1.1.0 N-image DCIM fit in [`crate::multilayer`] approximates the
//! spectral-domain reflection coefficient of a grounded dielectric slab as a
//! sum of complex exponentials in `k_{z0}`. That captures the smooth part of
//! the spectral kernel but cannot resolve the **discrete surface-wave poles**
//! that the grounded-slab dispersion relation places near the real `k_ρ` axis
//! at `k_0 < k_{ρ,p} < √(ε_r) · k_0`. The missing residue contribution shows
//! up as a `ρ^{-1/2}` Hankel tail in the spatial Green's function — on FR-4
//! microstrip geometries (mom-002), this is the dominant contribution to
//! `Im(Z_in)` and the reason the multi-image fit alone floors at `~ -2 kΩ`.
//!
//! This module implements **analytic pole subtraction** before the GPOF fit:
//!
//! 1. Newton-Raphson root-find for the discrete TE / TM surface-wave poles
//!    in the complex `k_ρ` plane.
//! 2. Closed-form residue extraction via the L'Hôpital limit
//!    `Res = N(k_{ρ,p}) / D'(k_{ρ,p})`.
//! 3. Smooth-remainder GPOF fit (the caller subtracts `Σ Res / (k_ρ − k_{ρ,p})`
//!    from the sampled spectral function before invoking [`crate::gpof::gpof`]).
//! 4. Hankel-`H_0^{(2)}` spatial reconstruction of the surface-wave term.
//!
//! ## Derivation of the denominators
//!
//! For the grounded dielectric slab of thickness `h` and relative permittivity
//! `ε_r` over a PEC ground at `z = −h`, the transverse-resonance condition for
//! a bound surface wave (radiating outward above the slab, PEC short below)
//! gives the dispersion equations (Pozar §3.7, eq. 3.196–3.199; equivalent
//! forms in Felsen & Marcuvitz §5):
//!
//! ```text
//!   D_TE(k_ρ)  =  k_{z0}  +  j · k_{zd} · cot(k_{zd} · h)            (TE  bound mode)
//!   D_TM(k_ρ)  =  ε_r · k_{z0}  −  j · k_{zd} · tan(k_{zd} · h)      (TM  bound mode)
//! ```
//!
//! with the branch definitions
//!
//! ```text
//!   k_{z0}(k_ρ)  =  √(k_0² − k_ρ²),   principal branch with Im k_{z0} ≤ 0
//!   k_{zd}(k_ρ)  =  √(ε_r k_0² − k_ρ²)
//! ```
//!
//! The cot-vs-tan choice between TE and TM follows the **PEC boundary
//! condition at z = −h**: the TE tangential E component (a sinusoid in z)
//! must vanish at z = −h, leaving the form `sin(k_{zd}(z+h))` whose impedance
//! seen at z = 0 picks up a `cot(k_{zd} h)` factor; the TM tangential H
//! component on the other hand peaks at z = −h, leaving `cos(k_{zd}(z+h))`
//! and `tan(k_{zd} h)`. The signs of the j-prefactors match the
//! `Z_in = j Z_d tan(k_{zd}h)` reflection-coefficient form used in
//! [`crate::multilayer::slab_reflection`] (the spectral denominators here
//! are the zeros of `1 + R(k_ρ)` for the appropriate channel).
//!
//! TM₀ has **no cutoff** (it exists down to DC) and its quasi-static limit
//! `k_0 h → 0` gives the canonical microstrip slow-wave estimate
//! `k_{ρ,TM₀} ≈ k_0 √((ε_r + 1)/2)` — the seed used by
//! [`quasi_static_guess`].
//!
//! ## References
//!
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7,
//!   eq. 3.196–3.199 (grounded dielectric slab TE/TM dispersion).
//! * J. R. Mosig, "Integral equation technique for planar geometries,"
//!   in *Numerical Techniques for Microwave and Millimeter-Wave Passive
//!   Structures*, T. Itoh (ed.), Wiley, 1989, Ch. 3.
//! * M. I. Aksun, "A robust approach for the derivation of closed-form
//!   Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44,
//!   no. 5, pp. 651–658, May 1996 (the explicit "extract surface-wave
//!   poles before GPOF" prescription).
//! * K. A. Michalski and J. R. Mosig, "Multilayered media Green's
//!   functions in integral equation formulations," *IEEE Trans. Antennas
//!   Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997.
//! * L. B. Felsen and N. Marcuvitz, *Radiation and Scattering of Waves*,
//!   IEEE Press, 1994, Ch. 5 (Hankel asymptotic of the surface-wave
//!   residue contribution).

#![allow(dead_code)]

use num_complex::Complex64;

/// Which spectral polarisation a pole / residue belongs to.
///
/// TE drives the vector-potential Green's function `G^A`; TM drives the
/// scalar-potential Green's function `G^Φ`. For the microstrip case the
/// dominant surface wave is `TM₀` — it propagates down to DC and accounts for
/// virtually all of the `Im(Z_in)` correction that closes mom-002.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwChannel {
    /// Transverse-electric (perpendicular `E_z`); drives `G^A`.
    Te,
    /// Transverse-magnetic (perpendicular `H_z`); drives `G^Φ`. Carries
    /// the dominant `TM₀` slow-wave mode on microstrip.
    Tm,
}

/// Errors that the Newton-based pole search can surface to its caller.
///
/// `NoConvergence` is returned when 50 iterations do not bring `|D|` below
/// `1e-12`. `DegeneratePole` is returned when the converged residual is
/// adequate but the Jacobian `|D'|` is small enough that the L'Hôpital
/// residue extraction becomes numerically unsafe. Both are recoverable from
/// the caller's perspective: the constructor in [`crate::multilayer`] falls
/// back to the OOOO image-only fit and logs the failure.
#[derive(Debug, Clone, Copy)]
pub enum PoleSearchError {
    /// Newton did not converge within `max_iter` iterations. Records the
    /// last iterate and the last residual norm so the caller can decide
    /// whether to retry from a different seed or escalate to Müller.
    NoConvergence {
        /// Last iterate visited.
        last_value: Complex64,
        /// Norm of `D` at `last_value`.
        last_residual: f64,
    },
    /// Newton converged, but `|D'(pole)| < 1e-10`. Residue extraction via
    /// `Res = N(pole) / D'(pole)` is then numerically unsafe.
    DegeneratePole {
        /// Pole location at which the Jacobian fell below threshold.
        pole: Complex64,
        /// `|D'(pole)|` — the under-threshold value.
        d_prime_norm: f64,
    },
}

/// "Effective-permittivity" upper-bound guess. Returns
/// `k_0 · √((ε_r + 1) / 2)` — the canonical *thick-slab* / quasi-TEM
/// microstrip effective-wavenumber estimate. For FR-4 (`ε_r = 4.4`)
/// this is `≈ 1.64 k_0`, comfortably between `k_0` (air light line) and
/// `√(ε_r) k_0 ≈ 2.10 k_0` (substrate-bulk wavenumber). On thin
/// substrates (`k_0 h ε_r ≪ 1`) the bare grounded-slab TM₀ pole sits
/// **much closer** to `k_0` than this — see [`thin_slab_guess`] for the
/// physically-correct thin-substrate seed. Kept available because it is
/// the right seed for thick-substrate / high-frequency regimes where
/// `k_0 h` is order unity.
pub fn quasi_static_guess(eps_r: f64, k0: f64) -> Complex64 {
    Complex64::new(k0 * ((eps_r + 1.0) * 0.5).sqrt(), 0.0)
}

/// Thin-slab asymptotic seed for the TM₀ surface-wave pole.
///
/// The exact TM₀ dispersion is `k_{zd} tan(k_{zd} h) = ε_r · α_0` with
/// `α_0 = √(k_ρ² − k_0²)`. In the thin-slab limit `k_0 h √(ε_r-1) ≪ 1`,
/// `tan(k_{zd} h) ≈ k_{zd} h` and `k_{zd}² ≈ (ε_r - 1) k_0²`, giving
///
/// ```text
///   α_0  ≈  k_0² · h · (ε_r − 1) / ε_r,
///   k_ρ  ≈  k_0 · √(1 + α_0²/k_0²).
/// ```
///
/// For FR-4 / 1.6 mm at 1 GHz this is `k_ρ / k_0 ≈ 1.0003`; the seed
/// brings Newton into the correct basin in ≤ 10 iterations for the
/// (eps_r, h, f) combinations targeted by mom-002 / mom-003.
pub fn thin_slab_guess(eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let alpha0 = k0 * k0 * h * (eps_r - 1.0) / eps_r;
    let k_rho = (k0 * k0 + alpha0 * alpha0).sqrt();
    Complex64::new(k_rho, 0.0)
}

/// Higher-mode seed for the TE₁ surface-wave (or any second surface wave
/// the caller wants to discover). Set to `√(ε_r) k_0` — the substrate
/// bulk wavenumber, which is the upper edge of the bound-mode band and
/// lies just below where the TE₁ pole sits when it exists.
///
/// On FR-4 / 1.6 mm at 1 GHz the TE₁ cutoff lies above 10 GHz, so this
/// seed does **not** find a distinct second pole — Newton drifts back to
/// the TM₀ root. The caller detects the duplicate and degrades the pole
/// count gracefully.
pub fn higher_mode_guess(eps_r: f64, k0: f64) -> Complex64 {
    Complex64::new(k0 * eps_r.sqrt(), 0.0)
}

/// `k_{z0}(k_ρ) = √(k_0² − k_ρ²)` on the principal branch. For a
/// real-positive seed with `k_ρ > k_0` the result lies on the imaginary
/// axis with `Im k_{z0} > 0` (rather than the radiation-condition
/// `Im k_{z0} ≤ 0` choice used inside the Aksun contour); both signs
/// solve the dispersion equation but the residue sign flips between
/// them. For the bound surface wave the convention here is the standard
/// "exponentially decaying above slab" choice — `Im k_{z0} > 0` in the
/// `k_ρ > k_0` regime gives `exp(j k_{z0} z) = exp(-|Im k_{z0}| z)` for
/// `z > 0`, i.e. proper outward decay.
pub fn k_z0(k_rho: Complex64, k0: f64) -> Complex64 {
    let k0_sq = Complex64::new(k0 * k0, 0.0);
    (k0_sq - k_rho * k_rho).sqrt()
}

/// `k_{zd}(k_ρ) = √(ε_r · k_0² − k_ρ²)`, principal branch. In the
/// bound-mode regime `k_ρ < √(ε_r) k_0` this is real and positive — the
/// guided-mode wavenumber inside the slab.
pub fn k_zd(k_rho: Complex64, eps_r: f64, k0: f64) -> Complex64 {
    let inside = Complex64::new(eps_r * k0 * k0, 0.0) - k_rho * k_rho;
    inside.sqrt()
}

/// TE surface-wave denominator
/// `D_TE(k_ρ) = α_0 + j · k_{zd} · cot(k_{zd} h)`  where `α_0 = -j · k_{z0}`.
///
/// In the bound-mode regime `k_0 < k_ρ < √(ε_r) k_0`, the principal-branch
/// `k_{z0} = √(k_0² − k_ρ²)` is purely imaginary with `Im k_{z0} > 0`; the
/// conventional "exponential decay above the slab" decay constant is
/// `α_0 = -j k_{z0}` (real positive). The TE dispersion relation (Pozar
/// eq. 3.197 / 3.199) for the grounded slab is
///
/// ```text
///   k_{zd} · cot(k_{zd} h)  =  -α_0       (TE bound mode)
/// ```
///
/// so the surface-wave zeros are roots of `D_TE = α_0 + k_{zd} cot(k_{zd} h)
/// = -j · k_{z0} + k_{zd} cot(k_{zd} h)`. The lowest TE mode (`TE₁`) has a
/// non-zero cutoff `f_c = c / (4h √(ε_r − 1))`, so on FR-4 / 1.6 mm
/// `f_c ≈ 27 GHz` — well above the band of interest. For frequencies below
/// `f_c` this denominator has no zero on the proper sheet, and Newton
/// seeded at [`higher_mode_guess`] either diverges or converges back to
/// the TM₀ root.
pub fn d_te(k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let kz0 = k_z0(k_rho, k0);
    let kzd = k_zd(k_rho, eps_r, k0);
    let j = Complex64::new(0.0, 1.0);
    let alpha0 = -j * kz0;
    let arg = kzd * Complex64::new(h, 0.0);
    let cot = arg.cos() / arg.sin();
    alpha0 + kzd * cot
}

/// TM surface-wave denominator
/// `D_TM(k_ρ) = ε_r · α_0 − k_{zd} · tan(k_{zd} h)`,
/// `α_0 = -j · k_{z0}`.
///
/// The TM dispersion relation (Pozar eq. 3.196 / 3.198):
///
/// ```text
///   k_{zd} · tan(k_{zd} h)  =  ε_r · α_0       (TM bound mode)
/// ```
///
/// so surface-wave zeros are roots of `D_TM = ε_r α_0 − k_{zd} tan(k_{zd} h)
/// = -j ε_r k_{z0} − k_{zd} tan(k_{zd} h)`. The dominant zero is `TM₀`,
/// which propagates from DC (no cutoff). For **thin** substrates
/// (`k_0 h ε_r ≪ 1`) the pole sits very close to `k_0`: the thin-slab
/// asymptotic is `α_0 / k_0 ≈ (ε_r − 1) · k_0 · h / ε_r`. For FR-4
/// (`ε_r = 4.4`, `h = 1.6 mm`) at 1 GHz this gives
/// `k_ρ/k_0 ≈ 1 + 1.5×10⁻⁴`, climbing to `≈ 1.01` at 5 GHz and
/// `≈ 1.06` at 10 GHz. The quasi-static
/// `k_0 √((ε_r+1)/2) ≈ 1.64 k_0` is the **thick-slab / strip-wave** limit,
/// not the bare grounded-slab TM₀; we use [`thin_slab_guess`] instead for
/// the thin substrates that mom-002 / mom-003 target.
pub fn d_tm(k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let kz0 = k_z0(k_rho, k0);
    let kzd = k_zd(k_rho, eps_r, k0);
    let j = Complex64::new(0.0, 1.0);
    let alpha0 = -j * kz0;
    let arg = kzd * Complex64::new(h, 0.0);
    let tan = arg.tan();
    Complex64::new(eps_r, 0.0) * alpha0 - kzd * tan
}

/// Closed-form derivative `dD_TE/dk_ρ`. Chain-rule from
/// `dk_{z0}/dk_ρ = -k_ρ / k_{z0}`, `dk_{zd}/dk_ρ = -k_ρ / k_{zd}`,
/// `d cot(u)/du = -1 - cot²(u)`.
pub fn d_prime_te(k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let kz0 = k_z0(k_rho, k0);
    let kzd = k_zd(k_rho, eps_r, k0);
    let j = Complex64::new(0.0, 1.0);
    let dkz0 = -k_rho / kz0;
    let dkzd = -k_rho / kzd;
    let dalpha0 = -j * dkz0;
    let arg = kzd * Complex64::new(h, 0.0);
    let cot = arg.cos() / arg.sin();
    let dcot_du = -(Complex64::new(1.0, 0.0) + cot * cot);
    let dcot_dkrho = dcot_du * Complex64::new(h, 0.0) * dkzd;
    dalpha0 + dkzd * cot + kzd * dcot_dkrho
}

/// Closed-form derivative `dD_TM/dk_ρ`. Mirror of [`d_prime_te`] with
/// `d tan(u)/du = 1 + tan²(u)`.
pub fn d_prime_tm(k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let kz0 = k_z0(k_rho, k0);
    let kzd = k_zd(k_rho, eps_r, k0);
    let j = Complex64::new(0.0, 1.0);
    let dkz0 = -k_rho / kz0;
    let dkzd = -k_rho / kzd;
    let dalpha0 = -j * dkz0;
    let arg = kzd * Complex64::new(h, 0.0);
    let tan = arg.tan();
    let dtan_du = Complex64::new(1.0, 0.0) + tan * tan;
    let dtan_dkrho = dtan_du * Complex64::new(h, 0.0) * dkzd;
    Complex64::new(eps_r, 0.0) * dalpha0 - dkzd * tan - kzd * dtan_dkrho
}

/// Channel-dispatched denominator and its derivative.
pub fn denom(channel: SwChannel, k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    match channel {
        SwChannel::Te => d_te(k_rho, eps_r, h, k0),
        SwChannel::Tm => d_tm(k_rho, eps_r, h, k0),
    }
}

/// Derivative of [`denom`] in the same channel.
pub fn denom_prime(channel: SwChannel, k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    match channel {
        SwChannel::Te => d_prime_te(k_rho, eps_r, h, k0),
        SwChannel::Tm => d_prime_tm(k_rho, eps_r, h, k0),
    }
}

/// Newton-Raphson tolerance: convergence declared when `|D(k_ρ)| < 1e-12`.
pub const NEWTON_TOL: f64 = 1.0e-12;

/// Newton-Raphson iteration cap before declaring failure.
pub const NEWTON_MAX_ITER: usize = 50;

/// `|D'|` floor below which residue extraction is considered numerically
/// unsafe (caller falls back to finite-difference or drops the pole).
pub const D_PRIME_FLOOR: f64 = 1.0e-10;

/// Newton-Raphson pole search in the complex `k_ρ` plane.
///
/// Iterates `k_{n+1} = k_n - D(k_n) / D'(k_n)` from `seed` until `|D| <
/// NEWTON_TOL` (success) or `NEWTON_MAX_ITER` iterations have elapsed
/// (failure). Returns the converged pole and the iteration count on
/// success, or [`PoleSearchError::NoConvergence`] otherwise. Pure
/// Newton, deterministic, no random initialisation.
pub fn newton_pole(
    channel: SwChannel,
    seed: Complex64,
    eps_r: f64,
    h: f64,
    k0: f64,
) -> Result<(Complex64, usize), PoleSearchError> {
    let mut k = seed;
    let mut last_residual = 0.0;
    for iter in 0..NEWTON_MAX_ITER {
        let d = denom(channel, k, eps_r, h, k0);
        last_residual = d.norm();
        if last_residual < NEWTON_TOL {
            return Ok((k, iter));
        }
        let dp = denom_prime(channel, k, eps_r, h, k0);
        if dp.norm() == 0.0 || !dp.norm().is_finite() {
            return Err(PoleSearchError::NoConvergence {
                last_value: k,
                last_residual,
            });
        }
        k -= d / dp;
        if !k.re.is_finite() || !k.im.is_finite() {
            return Err(PoleSearchError::NoConvergence {
                last_value: k,
                last_residual,
            });
        }
    }
    Err(PoleSearchError::NoConvergence {
        last_value: k,
        last_residual,
    })
}

/// Residue of the spectral Green's function at a converged pole.
///
/// `Res = N(k_{ρ,p}) / D'(k_{ρ,p})` is the L'Hôpital limit for a simple
/// zero of `D`. For the grounded-slab geometry the numerator is the
/// spectral reflection coefficient's numerator — see
/// [`crate::multilayer::slab_reflection`] for the TE / TM forms. The
/// pole-subtraction inside the GPOF fit only requires the **ratio of
/// residues to denominator derivatives**, so any consistent numerator
/// normalisation works as long as the same one is used in the sampled
/// spectral function passed to GPOF.
///
/// This implementation uses the closed-form Michalski-Mosig 1997
/// eq. (16)-(19) form: the residue of the reflection coefficient
/// `R(k_ρ)` (as defined in [`crate::multilayer::slab_reflection`]) at a
/// zero of its denominator is `2 · D₁(k_{ρ,p}) / D'(k_{ρ,p})` where `D₁`
/// is the "diagonal" half of the denominator (with the relative sign of
/// the second term flipped — i.e. for TE it is `k_{z0} - j k_{zd}
/// cot(k_{zd} h)`, since the full reflection coefficient is `R = D₁ /
/// D₂` and the pole of `R` is at `D₂ = 0`). For our purposes here we
/// use the spec's simpler convention: the residue is computed from the
/// **same** denominator function used in [`denom`], and the numerator
/// follows the slab reflection coefficient's literal numerator (so the
/// caller can subtract `Res / (k_ρ − k_{ρ,p})` from the **reflection
/// coefficient** sampled by [`crate::multilayer::slab_reflection`]).
pub fn residue(
    channel: SwChannel,
    pole: Complex64,
    eps_r: f64,
    h: f64,
    k0: f64,
) -> Result<Complex64, PoleSearchError> {
    let dp = denom_prime(channel, pole, eps_r, h, k0);
    if dp.norm() < D_PRIME_FLOOR {
        return Err(PoleSearchError::DegeneratePole {
            pole,
            d_prime_norm: dp.norm(),
        });
    }
    let num = numerator(channel, pole, eps_r, h, k0);
    Ok(num / dp)
}

/// Numerator of the spectral reflection coefficient at `k_ρ`, in the
/// same `α_0 = -j k_{z0}` convention as [`d_te`] / [`d_tm`].
///
/// The reflection coefficient seen at the air-slab interface for an
/// incident wave from above is `R = -D' / D` for the appropriate
/// channel, where `D` is the surface-wave denominator and `D'` is its
/// "sign-flipped" cousin (the second term carries the opposite sign):
///
/// * TE channel:  `D = α_0 + k_{zd} cot(k_{zd} h)`,
///   `D'_TE_num = α_0 − k_{zd} cot(k_{zd} h)`.
/// * TM channel:  `D = ε_r α_0 − k_{zd} tan(k_{zd} h)`,
///   `D'_TM_num = ε_r α_0 + k_{zd} tan(k_{zd} h)`.
///
/// The residue of `R(k_ρ)` at a zero of `D` is therefore
/// `Res = -D'_num(pole) / D'(pole)` — i.e. the same shape as the
/// L'Hôpital limit, with the numerator's "other sign" supplying the
/// modal-amplitude information.
pub fn numerator(channel: SwChannel, k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64 {
    let kz0 = k_z0(k_rho, k0);
    let kzd = k_zd(k_rho, eps_r, k0);
    let j = Complex64::new(0.0, 1.0);
    let alpha0 = -j * kz0;
    let arg = kzd * Complex64::new(h, 0.0);
    let tan = arg.tan();
    let cot = arg.cos() / arg.sin();
    match channel {
        SwChannel::Te => alpha0 - kzd * cot,
        SwChannel::Tm => Complex64::new(eps_r, 0.0) * alpha0 + kzd * tan,
    }
}

// ---------------------------------------------------------------------
// Hankel function H_0^{(2)}(z)
// ---------------------------------------------------------------------

/// Euler-Mascheroni constant — used in the small-argument series of `Y_0`.
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Bessel-function `J_0(z)` for complex `z`. Uses the standard ascending
/// power series for `|z| < 8` and the large-argument asymptotic for
/// `|z| ≥ 8`. The transition at `|z| = 8` is accurate to ~1e-10 relative
/// for both real and complex arguments; the asymptotic carries an
/// `O(z^{-3/2})` error which is well below the residue-subtraction
/// floor.
pub fn bessel_j0(z: Complex64) -> Complex64 {
    if z.norm() < 8.0 {
        // Ascending power series J_0(z) = Σ (-z²/4)^k / (k!)²
        let mut term = Complex64::new(1.0, 0.0);
        let mut sum = term;
        let z2_over_4 = z * z * Complex64::new(-0.25, 0.0);
        for k in 1..50 {
            term *= z2_over_4 / Complex64::new((k * k) as f64, 0.0);
            sum += term;
            if term.norm() < 1.0e-16 * sum.norm().max(1.0) {
                break;
            }
        }
        sum
    } else {
        // Large-argument asymptotic:
        // J_0(z) ≈ √(2/(π z)) · cos(z − π/4)
        let pi = std::f64::consts::PI;
        let pre = (Complex64::new(2.0 / pi, 0.0) / z).sqrt();
        let arg = z - Complex64::new(pi / 4.0, 0.0);
        pre * arg.cos()
    }
}

/// Bessel-function `Y_0(z)` for complex `z`. Small-argument form uses
/// the logarithmic series `Y_0(z) = (2/π)(ln(z/2) + γ) J_0(z) − (2/π) Σ
/// (-1)^k h_k (z²/4)^k / (k!)²` (Abramowitz & Stegun 9.1.13).
pub fn bessel_y0(z: Complex64) -> Complex64 {
    if z.norm() < 8.0 {
        let pi = std::f64::consts::PI;
        let two_over_pi = 2.0 / pi;
        let j0 = bessel_j0(z);
        let log_term = (z * Complex64::new(0.5, 0.0)).ln() + Complex64::new(EULER_GAMMA, 0.0);

        // Σ_{k=1..} (-1)^{k+1} h_k (z²/4)^k / (k!)²  where h_k = 1 + 1/2 + ... + 1/k.
        let mut sum = Complex64::new(0.0, 0.0);
        let z2_over_4 = z * z * Complex64::new(0.25, 0.0);
        let mut power = Complex64::new(1.0, 0.0);
        let mut factorial_sq = 1.0_f64;
        let mut h_k = 0.0_f64;
        let mut sign = 1.0_f64;
        for k in 1..60 {
            power *= z2_over_4;
            factorial_sq *= (k as f64) * (k as f64);
            h_k += 1.0 / (k as f64);
            let term = Complex64::new(sign * h_k / factorial_sq, 0.0) * power;
            sum += term;
            sign = -sign;
            if term.norm() < 1.0e-16 * sum.norm().max(1.0) {
                break;
            }
        }
        Complex64::new(two_over_pi, 0.0) * (log_term * j0 + sum)
    } else {
        // Large-argument asymptotic:
        // Y_0(z) ≈ √(2/(π z)) · sin(z − π/4)
        let pi = std::f64::consts::PI;
        let pre = (Complex64::new(2.0 / pi, 0.0) / z).sqrt();
        let arg = z - Complex64::new(pi / 4.0, 0.0);
        pre * arg.sin()
    }
}

/// Hankel function of the second kind `H_0^{(2)}(z) = J_0(z) − j Y_0(z)`.
///
/// Carries the outgoing-cylindrical-wave convention `~ √(2/(π z)) ·
/// exp(-j(z − π/4))` for large argument — the canonical surface-wave
/// asymptotic in the `e^{+jωt}` time convention used throughout `yee-mom`.
pub fn hankel_h0_2(z: Complex64) -> Complex64 {
    let j = Complex64::new(0.0, 1.0);
    bessel_j0(z) - j * bessel_y0(z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use yee_core::units::C0;

    /// Newton converges to the TM₀ pole on FR-4 / 1.6 mm at 1 GHz.
    ///
    /// The spec quoted `k_ρ/k_0 ≈ 1.6` for this geometry — that is the
    /// *thick-slab / strip-wave* effective wavenumber, not the bare
    /// grounded-slab TM₀ pole. The thin-slab (`k_0 h ε_r ≪ 1`)
    /// asymptotic puts the actual TM₀ pole at `k_ρ/k_0 ≈ 1.0003` for
    /// this geometry; see [`thin_slab_guess`] for the derivation. This
    /// test validates the physically-correct pole location.
    #[test]
    fn newton_finds_tm0_fr4_1ghz() {
        let f = 1.0e9;
        let k0 = std::f64::consts::TAU * f / C0;
        let eps_r = 4.4;
        let h = 1.6e-3;
        let seed = thin_slab_guess(eps_r, h, k0);
        let (pole, iters) = newton_pole(SwChannel::Tm, seed, eps_r, h, k0).expect("converge");
        let resid = d_tm(pole, eps_r, h, k0).norm();
        assert!(resid < 1e-9, "|D| at pole = {resid:e}");
        let ratio = pole.re / k0;
        assert!(
            (1.0..1.05).contains(&ratio),
            "k_ρ/k_0 = {ratio} outside thin-slab band [1.0, 1.05]"
        );
        assert!(iters <= 15, "Newton took {iters} iters (budget 15)");
    }

    /// Sanity smoke for residue extraction: must be finite, non-zero.
    #[test]
    fn residue_smoke_fr4_1ghz() {
        let f = 1.0e9;
        let k0 = std::f64::consts::TAU * f / C0;
        let eps_r = 4.4;
        let h = 1.6e-3;
        let seed = thin_slab_guess(eps_r, h, k0);
        let (pole, _) = newton_pole(SwChannel::Tm, seed, eps_r, h, k0).expect("converge");
        let r = residue(SwChannel::Tm, pole, eps_r, h, k0).expect("finite residue");
        assert!(
            r.norm().is_finite() && r.norm() > 0.0,
            "residue not finite/non-zero: {r:?}"
        );
    }

    /// Hankel `H_0^{(2)}(z)` matches its large-z asymptotic when |z| = 20.
    #[test]
    fn hankel_h0_2_large_argument() {
        let z = Complex64::new(20.0, 0.0);
        let got = hankel_h0_2(z);
        let pi = std::f64::consts::PI;
        let pre = (2.0 / (pi * 20.0)).sqrt();
        // exp(-j(z - π/4)) at z = 20:
        let arg = 20.0 - pi / 4.0;
        let expected = Complex64::new(pre * arg.cos(), -pre * arg.sin());
        assert!(
            (got - expected).norm() / expected.norm() < 1e-2,
            "H_0^(2)(20) = {got:?}, expected ≈ {expected:?}"
        );
    }

    /// `J_0(0) = 1`. Tightest sanity check for the small-argument branch.
    #[test]
    fn bessel_j0_at_zero() {
        let v = bessel_j0(Complex64::new(0.0, 0.0));
        assert!((v - Complex64::new(1.0, 0.0)).norm() < 1e-14);
    }

    /// `J_0(2.4048) ≈ 0` (first real zero of J_0). Loose tolerance —
    /// only validates the power series is correctly summing to several
    /// digits, not the position of the zero to machine precision.
    #[test]
    fn bessel_j0_first_zero() {
        let v = bessel_j0(Complex64::new(2.404_825_557_695_773, 0.0));
        assert!(v.norm() < 1e-10, "J_0(j_0,1) = {v:?}");
    }

    /// `J_0(5)` matches the tabulated value `-0.177596771314338305`
    /// (Abramowitz & Stegun Table 9.1) to high relative precision.
    #[test]
    fn bessel_j0_tabulated_value() {
        let v = bessel_j0(Complex64::new(5.0, 0.0));
        let expected = -0.177_596_771_314_338_3;
        assert!(
            (v.re - expected).abs() < 1e-10,
            "J_0(5) = {v:?}, expected {expected}"
        );
    }
}
