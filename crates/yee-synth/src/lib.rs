//! # yee-synth
//!
//! Classical filter-synthesis math for the Yee electromagnetic-simulation
//! studio (Filter Phase F0). Pure math: no EM, no I/O.
//!
//! This crate produces the **lowpass-prototype g-values** for Butterworth
//! (maximally flat) and Chebyshev (equi-ripple) responses, the
//! lowpass → bandpass frequency transform, and the all-pole **coupling
//! coefficients + external Q** synthesis from those g-values. Everything
//! downstream (the [`yee-filter`](../yee_filter/index.html) data model, ideal
//! response, and the `yee filter synth` CLI) plugs into the output of this
//! crate.
//!
//! ## References
//!
//! - Pozar, *Microwave Engineering* 4e, §8.3–8.4 (prototype g-values,
//!   Tables 8.3/8.4; the Chebyshev recursion is eq. 8.53).
//! - Matthaei, Young & Jones, *Microwave Filters, Impedance-Matching
//!   Networks…*, Table 4.05-2 (Chebyshev g-values).
//! - Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*,
//!   ch. 8 (coupling coefficients, external Q).
//!
//! ## Example
//!
//! ```
//! use yee_synth::{Approximation, prototype, coupling_design};
//!
//! // 0.5 dB-ripple Chebyshev, order 3.
//! let proto = prototype(Approximation::Chebyshev { ripple_db: 0.5 }, 3);
//! assert!((proto.g[1] - 1.5963).abs() < 1e-3);
//!
//! // Coupling design at 10% fractional bandwidth.
//! let design = coupling_design(&proto, 0.10);
//! assert!((design.k[0] - design.k[1]).abs() < 1e-12); // synchronous symmetry
//! ```

use serde::{Deserialize, Serialize};

/// Filter approximation (response shape) of a lowpass prototype.
///
/// Re-exported by `yee-filter` so a [`crate::FilterSpec`](../yee_filter/struct.FilterSpec.html)
/// can name the approximation directly.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Approximation {
    /// Maximally flat (Butterworth) response.
    Butterworth,
    /// Equi-ripple (Chebyshev) response with the given passband ripple in dB.
    Chebyshev {
        /// Passband ripple `L_Ar` in dB (e.g. `0.5` for a 0.5 dB-ripple filter).
        ripple_db: f64,
    },
}

/// A lowpass-prototype element-value vector.
///
/// `g[0]` is the source termination `g0`, `g[1..=N]` are the reactive element
/// values `g1..gN`, and `g[N+1]` is the load termination `g_{N+1}`. The vector
/// therefore has length `N + 2` for an order-`N` prototype.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Prototype {
    /// Element values `[g0, g1, …, gN, g_{N+1}]`, length `N + 2`.
    pub g: Vec<f64>,
}

impl Prototype {
    /// The filter order `N` (number of reactive elements `g1..gN`).
    pub fn order(&self) -> usize {
        // `g` is `[g0, g1, …, gN, g_{N+1}]`, length `N + 2`.
        self.g.len().saturating_sub(2)
    }
}

/// Compute the lowpass-prototype g-values for the given approximation and
/// order `N`.
///
/// The returned [`Prototype`] holds `[g0, g1, …, gN, g_{N+1}]` (length
/// `N + 2`). See [`Approximation`] for the supported response shapes and the
/// crate-level docs for the references.
///
/// # Panics
///
/// Panics if `order == 0` (a zeroth-order prototype is not defined).
pub fn prototype(approx: Approximation, order: usize) -> Prototype {
    assert!(order >= 1, "filter order must be >= 1, got {order}");
    match approx {
        Approximation::Butterworth => butterworth(order),
        Approximation::Chebyshev { ripple_db } => chebyshev(ripple_db, order),
    }
}

/// Butterworth (maximally flat) prototype, spec §2.1:
/// `g0 = 1`; `g_k = 2·sin((2k−1)·π / (2N))`, `k = 1..N`; `g_{N+1} = 1`.
fn butterworth(n: usize) -> Prototype {
    let nf = n as f64;
    let mut g = Vec::with_capacity(n + 2);
    g.push(1.0); // g0
    for k in 1..=n {
        let kf = k as f64;
        g.push(2.0 * ((2.0 * kf - 1.0) * std::f64::consts::PI / (2.0 * nf)).sin());
    }
    g.push(1.0); // g_{N+1}
    Prototype { g }
}

/// Chebyshev (equi-ripple) prototype, spec §2.2 (Pozar eq. 8.53):
/// ```text
/// β   = ln( coth( L_Ar / 17.37 ) )
/// γ   = sinh( β / (2N) )
/// a_k = sin( (2k−1)·π / (2N) ),  k = 1..N
/// b_k = γ² + sin²( k·π / N ),    k = 1..N
/// g1  = 2·a_1 / γ
/// g_k = (4·a_{k−1}·a_k) / (b_{k−1}·g_{k−1}),  k = 2..N
/// g_{N+1} = 1                   (N odd)
/// g_{N+1} = coth²( β / 4 )      (N even)
/// ```
/// The `17.37 = 40/ln(10)` constant is the standard Pozar form.
fn chebyshev(ripple_db: f64, n: usize) -> Prototype {
    let nf = n as f64;
    let beta = (ripple_db / 17.37).cosh() / (ripple_db / 17.37).sinh(); // coth(L_Ar/17.37)
    let beta = beta.ln();
    let gamma = (beta / (2.0 * nf)).sinh();

    // a_k and b_k for k = 1..N (1-based), stored 0-based in `a` / `b`.
    let a: Vec<f64> = (1..=n)
        .map(|k| {
            let kf = k as f64;
            ((2.0 * kf - 1.0) * std::f64::consts::PI / (2.0 * nf)).sin()
        })
        .collect();
    let b: Vec<f64> = (1..=n)
        .map(|k| {
            let kf = k as f64;
            let s = (kf * std::f64::consts::PI / nf).sin();
            gamma * gamma + s * s
        })
        .collect();

    let mut g = Vec::with_capacity(n + 2);
    g.push(1.0); // g0
    // g1 = 2·a_1 / γ
    g.push(2.0 * a[0] / gamma);
    // g_k = 4·a_{k−1}·a_k / (b_{k−1}·g_{k−1}), k = 2..N
    for k in 2..=n {
        let prev = g[k - 1];
        let gk = 4.0 * a[k - 2] * a[k - 1] / (b[k - 2] * prev);
        g.push(gk);
    }
    // g_{N+1}: 1 for odd N, coth²(β/4) for even N.
    if n % 2 == 1 {
        g.push(1.0);
    } else {
        let coth_quarter = (beta / 4.0).cosh() / (beta / 4.0).sinh();
        g.push(coth_quarter * coth_quarter);
    }
    Prototype { g }
}

/// Estimate the minimum filter order `N` to meet a required stopband rejection.
///
/// Given the approximation (whose ripple/return-loss sets the passband-edge
/// reference), the required rejection `rejection_db` at the stopband ratio
/// `omega_s = ω_s / ω_c`, returns the smallest integer order meeting the spec
/// (spec §2.3):
///
/// - Butterworth:
///   `N ≥ log10((10^{A_s/10} − 1) / (10^{L_Ar/10} − 1)) / (2·log10 Ω_s)`
/// - Chebyshev:
///   `N ≥ acosh(√((10^{A_s/10} − 1) / (10^{L_Ar/10} − 1))) / acosh(Ω_s)`
///
/// For Butterworth there is no passband ripple, so `L_Ar` is taken as the
/// 3 dB band-edge reference (`10^{L_Ar/10} − 1 = 1`).
///
/// # Panics
///
/// Panics if `omega_s <= 1.0` (the stopband must be outside the passband edge).
pub fn min_order(approx: Approximation, rejection_db: f64, omega_s: f64) -> usize {
    assert!(
        omega_s > 1.0,
        "stopband ratio omega_s must be > 1, got {omega_s}"
    );
    let num = 10f64.powf(rejection_db / 10.0) - 1.0;
    let n_real = match approx {
        Approximation::Butterworth => {
            // 3 dB band edge → denominator (10^{L_Ar/10} − 1) = 1.
            num.log10() / (2.0 * omega_s.log10())
        }
        Approximation::Chebyshev { ripple_db } => {
            let den = 10f64.powf(ripple_db / 10.0) - 1.0;
            (num / den).sqrt().acosh() / omega_s.acosh()
        }
    };
    n_real.ceil().max(1.0) as usize
}

/// All-pole, synchronous coupling-design outputs (spec §2.5).
///
/// Holds the inter-resonator coupling coefficients `k`, the input/output
/// external quality factors, and the normalized `N × N` coupling matrix `m`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CouplingDesign {
    /// Inter-resonator coupling coefficients `k_{i,i+1}`, `i = 1..N−1`
    /// (length `N − 1`).
    pub k: Vec<f64>,
    /// Input external quality factor `Qe_in = g0·g1 / FBW`.
    pub qe_in: f64,
    /// Output external quality factor `Qe_out = g_N·g_{N+1} / FBW`.
    pub qe_out: f64,
    /// Normalized `N × N` coupling matrix (synchronous → zero diagonal);
    /// `m[i][i+1] = m[i+1][i] = 1/√(g_i g_{i+1})`, all other entries 0.
    pub m: Vec<Vec<f64>>,
}

/// Synthesize the all-pole coupling coefficients, external Q, and normalized
/// coupling matrix from a lowpass prototype at fractional bandwidth `fbw`
/// (spec §2.5):
///
/// ```text
/// k_{i,i+1} = FBW / √( g_i · g_{i+1} ),  i = 1..N−1
/// Qe_in  = g0·g1 / FBW
/// Qe_out = g_N·g_{N+1} / FBW
/// M[i][i+1] = M[i+1][i] = 1/√(g_i g_{i+1})
/// ```
///
/// # Panics
///
/// Panics if `proto.order() < 1` or `fbw <= 0.0`.
pub fn coupling_design(proto: &Prototype, fbw: f64) -> CouplingDesign {
    let n = proto.order();
    assert!(n >= 1, "prototype order must be >= 1, got {n}");
    assert!(fbw > 0.0, "fractional bandwidth must be > 0, got {fbw}");
    let g = &proto.g; // g[0]=g0 .. g[N+1]

    // k_{i,i+1} = FBW / √(g_i g_{i+1}), i = 1..N−1.
    let mut k = Vec::with_capacity(n.saturating_sub(1));
    for i in 1..n {
        k.push(fbw / (g[i] * g[i + 1]).sqrt());
    }

    let qe_in = g[0] * g[1] / fbw;
    let qe_out = g[n] * g[n + 1] / fbw;

    // Normalized N×N coupling matrix: M[i][i+1] = M[i+1][i] = 1/√(g_i g_{i+1}).
    let mut m = vec![vec![0.0_f64; n]; n];
    for i in 1..n {
        // resonators are 1-based in the formula; matrix is 0-based.
        let val = 1.0 / (g[i] * g[i + 1]).sqrt();
        m[i - 1][i] = val;
        m[i][i - 1] = val;
    }

    CouplingDesign {
        k,
        qe_in,
        qe_out,
        m,
    }
}

/// Lowpass → bandpass frequency map (spec §2.4).
///
/// Given a centre frequency `omega0`, fractional bandwidth `fbw`, and an
/// evaluation angular frequency `omega`, returns the prototype lowpass variable
/// `Ω = (1/FBW)·(ω/ω0 − ω0/ω)`.
///
/// `omega0` and `omega` may be expressed in any consistent units (Hz or rad/s);
/// only the ratio matters.
///
/// # Panics
///
/// Panics if `omega0 <= 0.0`, `omega <= 0.0`, or `fbw <= 0.0`.
pub fn lowpass_to_bandpass(omega: f64, omega0: f64, fbw: f64) -> f64 {
    assert!(omega0 > 0.0, "omega0 must be > 0, got {omega0}");
    assert!(omega > 0.0, "omega must be > 0, got {omega}");
    assert!(fbw > 0.0, "fbw must be > 0, got {fbw}");
    (1.0 / fbw) * (omega / omega0 - omega0 / omega)
}
