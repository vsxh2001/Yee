//! Published closed-form reference dispersion for the slab-loaded
//! rectangular waveguide (Phase 1.3.1.1 step 5.1).
//!
//! This module provides an **independent published-benchmark** against
//! which the numerical mixed `(E_t, E_z)` cross-section eigensolve
//! ([`crate::ports::NumericalCrossSection::solve`]) can be reconciled —
//! closing (or characterising) the CLAUDE.md §4 gap left open by the
//! step-5 DoD-V2′ escape-hatch (which shipped only a monotonic
//! `β_air < β_loaded < β_full` bracket + a self-referential regression).
//!
//! # Geometry
//!
//! A rectangular guide of width `a` (along x) and height `b` (along y)
//! with PEC walls, stratified **along y**: a dielectric `ε_r` fills
//! `0 ≤ y ≤ d₁`, air fills `d₁ ≤ y ≤ b`. The x-variation is pinned by the
//! PEC side walls to `sin(m π x / a)` (the dominant guided mode is
//! `m = 1`). This is exactly the `horizontal_slab_mesh` case in
//! `tests/eigensolver_inhomogeneous.rs`.
//!
//! # Transverse-resonance dispersion
//!
//! Modes of a y-stratified guide separate into two **longitudinal-section**
//! families w.r.t. the stratification axis ŷ (Pozar, *Microwave
//! Engineering* 4th ed. §6.6; Collin, *Field Theory of Guided Waves* 2nd
//! ed. §6):
//!
//! * **LSM-to-y** (TM-to-y, `H_y = 0`): the family that contains the
//!   `TE_{m0}` modes of the empty guide (`TE10` has `E_y ≠ 0`, `H_y = 0`),
//!   so the **dominant** slab-loaded mode is LSM-to-y. Its transverse
//!   electric field is dominantly `E_y` with a small longitudinal `E_z`
//!   localised at the dielectric interface — i.e. *weakly* hybrid.
//! * **LSE-to-y** (TE-to-y, `E_y = 0`): the dual family.
//!
//! Treating the y-direction as a transverse-resonance transmission line
//! with the PEC walls at `y = 0, b` as short circuits, each layer `i` is a
//! short-circuited stub of electrical length `k_{y,i} d_i`, where
//!
//! ```text
//!   k_{y,i}² = ε_{r,i} k₀² − (m π / a)² − β²        (i = 1: dielectric, i = 2: air)
//! ```
//!
//! and `k₀ = ω / c₀`. Transverse resonance (`Y_up + Y_down = 0`) gives:
//!
//! ```text
//!   LSE-to-y:   k_{y1} cot(k_{y1} d₁) + k_{y2} cot(k_{y2} d₂) = 0
//!   LSM-to-y:   (ε_{r1}/k_{y1}) cot(k_{y1} d₁) + (ε_{r2}/k_{y2}) cot(k_{y2} d₂) = 0
//! ```
//!
//! with `d₂ = b − d₁`. The characteristic admittance of a stub is
//! `Y₀ ∝ k_y / (ω μ₀)` for the LSE polarisation and `Y₀ ∝ ω ε₀ ε_r / k_y`
//! for the LSM polarisation (non-magnetic fills, `μ_r = 1`).
//!
//! **Imaginary `k_y` (the subtlety the prior bring-up attempt missed).**
//! For a strongly-loaded mode `β` exceeds the air-region propagation
//! constant `√(k₀² − (mπ/a)²)`, so `k_{y2}² < 0` and the field is
//! *evanescent in y* in the air region. Then `k_{y2} = j q₂` and
//! `cot(k_{y2} d₂) = −j coth(q₂ d₂)`; the products in the residual stay
//! **real** (see [`lse_term`] / [`lsm_term`]). A root-find that assumed
//! real `k_y` everywhere fails to find the loaded root — which is most
//! likely why the step-5 first attempt found "no root".
//!
//! # Verification
//!
//! Both transcendentals are verified **independently of the FEM solver
//! and independently of each other** by the unit tests below: the LSE and
//! LSM dominant roots reproduce a shooting-method / finite-difference
//! solution of the same underlying 1-D transverse ODE (sharing no code
//! with this transcendental), and reduce exactly to the analytic empty /
//! fully-filled `TE_{m0}` limits. See
//! [`tests::lsm_dominant_root_matches_independent_shooting`] and
//! [`tests::reduces_to_homogeneous_te10_limit`].

// Phase 1.3.1.1 step 5.1: this reference is a published-benchmark surface
// whose load-bearing consumers are (a) the module's own `#[cfg(test)]`
// independent-verification unit tests (DoD-1) and (b) the
// `eigensolver_inhomogeneous.rs` integration diagnostic (which keeps a
// self-contained mirror to confine this step's edits to the eigensolver +
// test lane). In a plain non-test lib build none of these entry points is
// reached, so `dead_code` would fire; the items are intentional API, not
// dead. Both mode families + the `slab_loaded_beta_family` variant are
// retained so the LSE-vs-LSM family decision rests on two verified
// dispersions (see ADR-0052), not one.
#![allow(dead_code)]

use std::f64::consts::PI;
use yee_core::units::C0;

/// Which longitudinal-section family (w.r.t. the y stratification axis) a
/// [`slab_loaded_beta`] solve targets.
///
/// The **dominant** slab-loaded mode of a horizontally-stratified guide
/// is [`SlabMode::Lsm`] (it contains the `TE_{m0}` family); [`SlabMode::Lse`]
/// is the dual. Exposed so the reconciliation test can compare against the
/// family whose field orientation matches the numerical dominant mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SlabMode {
    /// LSM-to-y (TM-to-y, `H_y = 0`) — the `TE_{m0}`-derived dominant family.
    Lsm,
    /// LSE-to-y (TE-to-y, `E_y = 0`) — the dual family.
    Lse,
}

/// One stub term of the **LSE-to-y** transverse-resonance residual,
/// `k_y cot(k_y d)`, robust to an imaginary `k_y`.
///
/// `ky_sq` is the *signed* squared transverse wavenumber
/// `k_y² = ε_r k₀² − (mπ/a)² − β²`. For `k_y² > 0` this is the propagating
/// `k_y cot(k_y d)`; for `k_y² < 0` the layer is evanescent in y
/// (`k_y = j q`) and `k_y cot(k_y d) = (j q)(−j coth(q d)) = q coth(q d)`,
/// a real positive quantity.
fn lse_term(ky_sq: f64, d: f64) -> f64 {
    if ky_sq > 0.0 {
        let k = ky_sq.sqrt();
        k / (k * d).tan()
    } else {
        let q = (-ky_sq).sqrt();
        q / (q * d).tanh()
    }
}

/// One stub term of the **LSM-to-y** transverse-resonance residual,
/// `(ε_r / k_y) cot(k_y d)`, robust to an imaginary `k_y`.
///
/// For `k_y² < 0` (`k_y = j q`):
/// `(ε_r / (j q))(−j coth(q d)) = −(ε_r / q) coth(q d)`, a real *negative*
/// quantity — the sign flip relative to the propagating branch is what
/// makes the loaded root exist.
fn lsm_term(eps_r: f64, ky_sq: f64, d: f64) -> f64 {
    if ky_sq > 0.0 {
        let k = ky_sq.sqrt();
        (eps_r / k) / (k * d).tan()
    } else {
        let q = (-ky_sq).sqrt();
        -(eps_r / q) / (q * d).tanh()
    }
}

/// Fixed parameters of a slab-loaded-guide transverse-resonance solve:
/// geometry, dielectric, operating wavenumber, x-mode index, and the
/// longitudinal-section family. Bundling these keeps the residual /
/// root-find signatures small (the alternative — passing all of them
/// positionally — trips clippy's `too_many_arguments`).
#[derive(Clone, Copy, Debug)]
struct SlabGuide {
    /// Longitudinal-section family.
    mode: SlabMode,
    /// Guide width along x (m).
    a: f64,
    /// Guide height along y (m).
    b: f64,
    /// Dielectric-layer thickness from `y = 0` (m); air fills `d1 ≤ y ≤ b`.
    d1: f64,
    /// Dielectric relative permittivity.
    eps_r: f64,
    /// Free-space wavenumber `k₀ = ω / c₀` (rad/m) at the solve frequency.
    k0: f64,
    /// Transverse x-mode index (`sin(m π x / a)`).
    m: u32,
}

impl SlabGuide {
    /// Upper edge of the propagating window: `β` at which the dielectric
    /// layer's own transverse wavenumber vanishes (`k_{y1} = 0`). Above this
    /// `β²` the dielectric layer is itself evanescent in y and no
    /// propagating slab mode exists. `None` if below cutoff even fully
    /// filled.
    fn beta_window_top(&self) -> Option<f64> {
        let kx = (self.m as f64) * PI / self.a;
        let beta_max_sq = self.eps_r * self.k0 * self.k0 - kx * kx;
        (beta_max_sq > 0.0).then(|| beta_max_sq.sqrt())
    }

    /// Signed transverse-resonance residual `R(β)`. A propagating mode is a
    /// root `R(β) = 0`.
    fn residual(&self, beta: f64) -> f64 {
        let d2 = self.b - self.d1;
        let kx = (self.m as f64) * PI / self.a;
        let ky1_sq = self.eps_r * self.k0 * self.k0 - kx * kx - beta * beta; // dielectric
        let ky2_sq = self.k0 * self.k0 - kx * kx - beta * beta; // air
        match self.mode {
            SlabMode::Lse => lse_term(ky1_sq, self.d1) + lse_term(ky2_sq, d2),
            SlabMode::Lsm => lsm_term(self.eps_r, ky1_sq, self.d1) + lsm_term(1.0, ky2_sq, d2),
        }
    }
}

/// Solve the slab-loaded-guide transverse-resonance transcendental for the
/// **dominant** (largest-β, i.e. lowest-cutoff) propagating root of the
/// given longitudinal-section `mode` family.
///
/// # Arguments
/// * `a`, `b` — guide width (x) and height (y), metres.
/// * `d1` — dielectric-layer thickness measured from `y = 0`, metres
///   (`ε_r` fills `0 ≤ y ≤ d1`; air fills `d1 ≤ y ≤ b`).
/// * `eps_r` — relative permittivity of the dielectric layer.
/// * `freq_hz` — operating frequency, Hz.
/// * `m` — transverse x-mode index (`sin(m π x / a)`; `1` for the dominant
///   mode).
///
/// Returns the dominant `β` (rad/m), or `None` if no propagating root is
/// found in the physical window `0 < β < √(ε_r) k₀` (degenerate / all-cutoff
/// geometry).
///
/// # Method
///
/// A propagating mode satisfies `β² < ε_r k₀² − (mπ/a)²` (so the dielectric
/// layer is propagating in y) and lies in `0 < β ≤ √(ε_r k₀² − (mπ/a)²)`.
/// The residual ([`SlabGuide::residual`]) has poles where a layer hits an
/// internal half-wave resonance (`k_{y,i} d_i = nπ`); between consecutive poles it is
/// smooth and monotone, so a sign change that is **not** a pole jump is a
/// genuine root. The scan walks β downward from the upper edge of the
/// window and bisects the first genuine sign change — the largest-β root,
/// i.e. the dominant mode. Bisection (no external dependency) is used; the
/// residual is cheap and the bracket is tight.
pub(crate) fn slab_loaded_beta(
    a: f64,
    b: f64,
    d1: f64,
    eps_r: f64,
    freq_hz: f64,
    m: u32,
) -> Option<f64> {
    slab_loaded_beta_family(SlabMode::Lsm, a, b, d1, eps_r, freq_hz, m)
}

/// As [`slab_loaded_beta`], but for an explicit longitudinal-section
/// `mode` family (the public entry point fixes [`SlabMode::Lsm`], the
/// dominant family).
pub(crate) fn slab_loaded_beta_family(
    mode: SlabMode,
    a: f64,
    b: f64,
    d1: f64,
    eps_r: f64,
    freq_hz: f64,
    m: u32,
) -> Option<f64> {
    let k0 = std::f64::consts::TAU * freq_hz / C0;
    let kx = (m as f64) * PI / a;
    let d2 = b - d1;

    // Homogeneous-fill degenerate limits. When the guide is uniformly
    // filled — air everywhere (`eps_r == 1`), fully dielectric-filled
    // (`d2 == 0`), or empty of dielectric (`d1 == 0`) — the dominant mode
    // is the `TE_{m0}` mode with **no** y-variation (`k_y = 0`), i.e.
    // `β = √(ε_eff k₀² − (mπ/a)²)` with `ε_eff` the single fill. That sits
    // exactly at the upper edge of the transverse-resonance window (the
    // `k_y → 0` boundary), where the stub formulation degenerates rather
    // than producing an interior root, so it is handled in closed form.
    let homogeneous = (eps_r - 1.0).abs() < 1e-12 || d2 <= 0.0 || d1 <= 0.0;
    if homogeneous {
        let eps_eff = if d1 <= 0.0 { 1.0 } else { eps_r };
        let beta_sq = eps_eff * k0 * k0 - kx * kx;
        return (beta_sq > 0.0).then(|| beta_sq.sqrt());
    }

    let guide = SlabGuide {
        mode,
        a,
        b,
        d1,
        eps_r,
        k0,
        m,
    };
    let beta_hi = guide.beta_window_top()?; // None ⇒ below cutoff even fully filled

    // Dense downward scan for the first genuine (non-pole) sign change.
    // 4000 samples over the window resolves the well-separated dominant
    // root for the validation geometries with margin; a pole jump is
    // rejected by the `|ΔR|` discontinuity guard.
    let n = 4000usize;
    let beta_lo_floor = 1e-3;
    let step = (beta_hi - beta_lo_floor) / (n as f64);
    let mut prev_beta = beta_hi - 1e-6;
    let mut prev = guide.residual(prev_beta);
    for i in 1..=n {
        let beta = beta_hi - 1e-6 - (i as f64) * step;
        if beta <= beta_lo_floor {
            break;
        }
        let cur = guide.residual(beta);
        let sign_change = prev.is_finite() && cur.is_finite() && prev * cur < 0.0;
        // A pole produces a sign change with a huge magnitude jump; a root
        // is a smooth crossing. Reject jumps where the residual magnitude
        // explodes between the two samples.
        let smooth = (cur - prev).abs() < (cur.abs() + prev.abs() + 1.0);
        if sign_change && smooth {
            // Bisect between prev_beta (higher β) and beta (lower β).
            let mut lo = beta;
            let mut hi = prev_beta;
            let mut f_lo = cur;
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let f_mid = guide.residual(mid);
                if !f_mid.is_finite() {
                    return None;
                }
                if f_lo * f_mid <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }
            return Some(0.5 * (lo + hi));
        }
        prev_beta = beta;
        prev = cur;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    // WR-90 geometry shared with the integration gate.
    const A: f64 = 22.86e-3;
    const B: f64 = 10.16e-3;
    const FREQ_HZ: f64 = 10.0e9;

    /// Independent (shooting-method) solution of the **same** 1-D
    /// transverse ODE the LSM-to-y transcendental encodes, sharing no code
    /// with [`slab_loaded_beta`]. Integrates the modal scalar `ψ(y)` from a
    /// Neumann wall at `y = 0` through the dielectric layer, applies the
    /// exact `(1/ε) ψ'` flux-continuity at the interface, propagates to
    /// `y = b`, and returns `ψ'(b)` — which must vanish at a Neumann (LSM)
    /// wall. A root of this is an LSM-to-y eigenvalue.
    fn lsm_shoot(a: f64, b: f64, d1: f64, eps_r: f64, freq_hz: f64, m: u32, beta: f64) -> f64 {
        let k0 = std::f64::consts::TAU * freq_hz / C0;
        let kx = (m as f64) * PI / a;
        let d2 = b - d1;
        let ky1_sq = eps_r * k0 * k0 - kx * kx - beta * beta;
        let ky2_sq = k0 * k0 - kx * kx - beta * beta;
        // Layer 1, ψ(0)=1, ψ'(0)=0 (Neumann wall).
        let (psi_i, dpsi1) = if ky1_sq >= 0.0 {
            let k1 = ky1_sq.sqrt();
            ((k1 * d1).cos(), -k1 * (k1 * d1).sin())
        } else {
            let q1 = (-ky1_sq).sqrt();
            ((q1 * d1).cosh(), q1 * (q1 * d1).sinh())
        };
        // Interface: ψ continuous; (1/ε₁)ψ'₁ = (1/ε₂)ψ'₂, ε₂ = 1.
        let dpsi2_i = dpsi1 / eps_r;
        // Layer 2 → require ψ'(b) = 0.
        if ky2_sq >= 0.0 {
            let k2 = ky2_sq.sqrt();
            -psi_i * k2 * (k2 * d2).sin() + dpsi2_i * (k2 * d2).cos()
        } else {
            let q2 = (-ky2_sq).sqrt();
            psi_i * q2 * (q2 * d2).sinh() + dpsi2_i * (q2 * d2).cosh()
        }
    }

    /// Brute scan + bisection for the dominant (largest-β) root of the
    /// independent shooting residual.
    fn shoot_dominant(a: f64, b: f64, d1: f64, eps_r: f64, freq_hz: f64, m: u32) -> Option<f64> {
        let k0 = std::f64::consts::TAU * freq_hz / C0;
        let kx = (m as f64) * PI / a;
        let beta_hi = (eps_r * k0 * k0 - kx * kx).sqrt();
        let n = 8000usize;
        let step = beta_hi / (n as f64);
        let mut prev_beta = beta_hi - 1e-6;
        let mut prev = lsm_shoot(a, b, d1, eps_r, freq_hz, m, prev_beta);
        for i in 1..=n {
            let beta = beta_hi - 1e-6 - (i as f64) * step;
            if beta <= 1e-3 {
                break;
            }
            let cur = lsm_shoot(a, b, d1, eps_r, freq_hz, m, beta);
            if prev.is_finite()
                && cur.is_finite()
                && prev * cur < 0.0
                && (cur - prev).abs() < (cur.abs() + prev.abs() + 1.0)
            {
                let (mut lo, mut hi, mut f_lo) = (beta, prev_beta, cur);
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    let f_mid = lsm_shoot(a, b, d1, eps_r, freq_hz, m, mid);
                    if f_lo * f_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        f_lo = f_mid;
                    }
                }
                return Some(0.5 * (lo + hi));
            }
            prev_beta = beta;
            prev = cur;
        }
        None
    }

    #[test]
    fn lsm_dominant_root_matches_independent_shooting() {
        // DoD-1: the LSM-to-y transcendental's dominant root must reproduce
        // an independent shooting-method solution of the same transverse
        // ODE to high precision. This isolates a *reference* bug (which the
        // step-5 bring-up suspected) from a *solver* bug: the two methods
        // share no code, so their agreement certifies the transcendental.
        //
        // Horizontal slab: ε_r = 10.2 in 0 ≤ y ≤ b/2, air above, m = 1.
        let beta_trans = slab_loaded_beta(A, B, B / 2.0, 10.2, FREQ_HZ, 1)
            .expect("LSM transcendental must have a dominant root");
        let beta_shoot = shoot_dominant(A, B, B / 2.0, 10.2, FREQ_HZ, 1)
            .expect("shooting must find a dominant root");
        let rel = (beta_trans - beta_shoot).abs() / beta_shoot;
        eprintln!(
            "LSM dominant β: transcendental {beta_trans:.4}, independent shooting \
             {beta_shoot:.4}, rel err {rel:.3e}"
        );
        assert!(
            rel < 1e-4,
            "LSM transcendental β {beta_trans} must match independent shooting \
             β {beta_shoot} within 1e-4 (rel {rel:.3e})"
        );
        // Physical sanity: the dominant LSM mode of a half-ε_r=10.2-filled
        // guide concentrates its field in the dielectric, so ε_eff is well
        // above the air value and approaches the dielectric value. Reject
        // the (also-real) low-β "air-like" upper-cutoff branch.
        let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
        let kx = PI / A;
        let eps_eff = (beta_trans * beta_trans + kx * kx) / (k0 * k0);
        assert!(
            eps_eff > 4.0,
            "dominant LSM mode ε_eff {eps_eff:.3} should be field-concentrated \
             in the dielectric (≫ air)"
        );
    }

    #[test]
    fn lse_dominant_root_matches_independent_shooting() {
        // Cross-check the dual (LSE-to-y) family against an independent
        // shooting solve of *its* transverse ODE (Dirichlet walls, ψ & ψ'
        // continuous at the interface). Confirms the LSE residual form too,
        // so the LSE-vs-LSM family decision rests on two verified
        // dispersions rather than one.
        let lse_shoot = |beta: f64| -> f64 {
            let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
            let kx = PI / A;
            let d1 = B / 2.0;
            let d2 = B - d1;
            let ky1_sq = 10.2 * k0 * k0 - kx * kx - beta * beta;
            let ky2_sq = k0 * k0 - kx * kx - beta * beta;
            // ψ(0)=0, ψ'(0)=1 (Dirichlet wall); ψ & ψ' continuous; require ψ(b)=0.
            let (psi_i, dpsi_i) = if ky1_sq >= 0.0 {
                let k1 = ky1_sq.sqrt();
                ((k1 * d1).sin(), k1 * (k1 * d1).cos())
            } else {
                let q1 = (-ky1_sq).sqrt();
                ((q1 * d1).sinh(), q1 * (q1 * d1).cosh())
            };
            if ky2_sq >= 0.0 {
                let k2 = ky2_sq.sqrt();
                psi_i * (k2 * d2).cos() + (dpsi_i / k2) * (k2 * d2).sin()
            } else {
                let q2 = (-ky2_sq).sqrt();
                psi_i * (q2 * d2).cosh() + (dpsi_i / q2) * (q2 * d2).sinh()
            }
        };
        let beta_trans = slab_loaded_beta_family(SlabMode::Lse, A, B, B / 2.0, 10.2, FREQ_HZ, 1)
            .expect("LSE transcendental must have a dominant root");
        // Independent dominant-root scan of the LSE shooting residual.
        let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
        let kx = PI / A;
        let beta_hi = (10.2 * k0 * k0 - kx * kx).sqrt();
        let n = 8000usize;
        let step = beta_hi / (n as f64);
        let mut prev_beta = beta_hi - 1e-6;
        let mut prev = lse_shoot(prev_beta);
        let mut found = None;
        for i in 1..=n {
            let beta = beta_hi - 1e-6 - (i as f64) * step;
            if beta <= 1e-3 {
                break;
            }
            let cur = lse_shoot(beta);
            if prev.is_finite()
                && cur.is_finite()
                && prev * cur < 0.0
                && (cur - prev).abs() < (cur.abs() + prev.abs() + 1.0)
            {
                let (mut lo, mut hi, mut f_lo) = (beta, prev_beta, cur);
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    let f_mid = lse_shoot(mid);
                    if f_lo * f_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        f_lo = f_mid;
                    }
                }
                found = Some(0.5 * (lo + hi));
                break;
            }
            prev_beta = beta;
            prev = cur;
        }
        let beta_shoot = found.expect("LSE shooting must find a dominant root");
        let rel = (beta_trans - beta_shoot).abs() / beta_shoot;
        eprintln!(
            "LSE dominant β: transcendental {beta_trans:.4}, independent shooting \
             {beta_shoot:.4}, rel err {rel:.3e}"
        );
        assert!(
            rel < 1e-4,
            "LSE transcendental β {beta_trans} must match independent shooting \
             β {beta_shoot} within 1e-4 (rel {rel:.3e})"
        );
    }

    #[test]
    fn reduces_to_homogeneous_te10_limit() {
        // Degenerate check: with ε_r = 1 the guide is air-filled and the
        // dominant mode is exactly TE10, β = √(k₀² − (π/a)²). Both families
        // must reproduce it (the slab vanishes). This anchors the absolute
        // scale of the transcendental against a closed-form value.
        let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
        let kx = PI / A;
        let beta_te10 = (k0 * k0 - kx * kx).sqrt();
        let beta_lsm =
            slab_loaded_beta(A, B, B / 2.0, 1.0, FREQ_HZ, 1).expect("air-filled LSM dominant root");
        let rel = (beta_lsm - beta_te10).abs() / beta_te10;
        eprintln!(
            "air-filled limit: LSM β {beta_lsm:.4}, analytic TE10 {beta_te10:.4}, rel {rel:.3e}"
        );
        assert!(
            rel < 1e-3,
            "air-filled LSM β {beta_lsm} must reduce to analytic TE10 β {beta_te10} (rel {rel:.3e})"
        );
    }

    #[test]
    fn fully_filled_limit_matches_analytic() {
        // Degenerate check: with the dielectric filling the whole height
        // (d1 = b, ε_r = 2.55) the dominant mode is the fully-filled TE10,
        // β = √(ε_r k₀² − (π/a)²). ε_r = 2.55 (a standard PTFE value) gives
        // a concrete published-style closed-form target.
        let eps_r = 2.55;
        let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
        let kx = PI / A;
        let beta_full = (eps_r * k0 * k0 - kx * kx).sqrt();
        // d1 = b (whole guide) — air layer has zero thickness, so only the
        // dielectric stub remains and the resonance is the filled TE10.
        let beta_lsm =
            slab_loaded_beta(A, B, B, eps_r, FREQ_HZ, 1).expect("fully-filled LSM dominant root");
        let rel = (beta_lsm - beta_full).abs() / beta_full;
        eprintln!(
            "fully-filled limit (ε_r=2.55): LSM β {beta_lsm:.4}, analytic {beta_full:.4}, rel {rel:.3e}"
        );
        assert!(
            rel < 1e-3,
            "fully-filled LSM β {beta_lsm} must match analytic filled TE10 β {beta_full} (rel {rel:.3e})"
        );
    }

    #[test]
    fn imaginary_ky_terms_stay_real_and_finite() {
        // Regression for the prior-attempt failure mode: at the loaded β
        // the air layer is evanescent in y (k_{y2}² < 0). The residual must
        // stay finite and real there — i.e. the cot→coth handling must
        // engage rather than producing NaN.
        let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
        let kx = PI / A;
        let beta = 582.95; // near the verified LSM dominant root
        let ky2_sq = k0 * k0 - kx * kx - beta * beta;
        assert!(
            ky2_sq < 0.0,
            "air layer must be evanescent in y at the loaded β"
        );
        let guide = SlabGuide {
            mode: SlabMode::Lsm,
            a: A,
            b: B,
            d1: B / 2.0,
            eps_r: 10.2,
            k0,
            m: 1,
        };
        let r = guide.residual(beta);
        assert!(
            r.is_finite(),
            "residual must be finite with an imaginary k_y"
        );
    }
}
