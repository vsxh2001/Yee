//! Phase 4.fem.eig.1 step D3 — `MaterialDatabase` ω → ε(ω), μ(ω) lookup.
//!
//! Gate test inventory:
//!
//! 1. `free_space_material_returns_one` — the [`Material::default`]
//!    free-space tag returns `ε = 1 + 0j` at any ω.
//! 2. `single_drude_pole_recovers_static_limit` — for a Drude material in
//!    the `ω → ∞` limit, `ε → ε_∞` exactly.
//! 3. `lorentz_resonance_peak` — at `ω = ω_0` the Lorentz contribution is
//!    purely imaginary; |Im(ε)| dominates |Re(ε − ε_∞)|.
//! 4. `debye_low_frequency_static_dielectric` — water-like Debye fit
//!    recovers `ε(0) = ε_∞ + Δε`.
//! 5. `database_lookup_by_tag` — two materials registered under distinct
//!    tags return distinct `ε(ω)` at the same `ω`.
//! 6. `database_missing_tag_returns_air` — unregistered tag → `ε = 1`.
//!
//! References:
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
//!   §6 + §7.
//! * `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`
//!   step D3.
//! * `docs/src/decisions/0039-phase-4-fem-eig-1-dispersive-scope.md`.

use std::f64::consts::PI;

use num_complex::Complex64;
use yee_fem::material::{Material, MaterialDatabase, MaterialPole};

#[test]
fn free_space_material_returns_one() {
    let m = Material::default();
    let eps = m.eps_at(2.0 * PI * 1.0e9);
    assert!(
        (eps.re - 1.0).abs() < 1e-15,
        "free-space Re(ε) should be 1, got {}",
        eps.re,
    );
    assert!(
        eps.im.abs() < 1e-15,
        "free-space Im(ε) should be 0, got {}",
        eps.im,
    );

    // μ defaults to 1 + 0j as well.
    let mu = m.mu_at(2.0 * PI * 1.0e9);
    assert!((mu.re - 1.0).abs() < 1e-15);
    assert!(mu.im.abs() < 1e-15);
}

#[test]
fn single_drude_pole_recovers_static_limit() {
    // Drude: ε(ω) = ε_∞ − ω_p² / (ω² + jγω). For ω → ∞, both real and
    // imaginary parts of the contribution → 0, so ε → ε_∞.
    let m = Material {
        eps_inf: 3.5,
        mu_r: 1.0,
        poles: vec![MaterialPole::Drude {
            omega_p: 2.0 * PI * 1.0e10,
            gamma: 2.0 * PI * 1.0e8,
        }],
    };
    let omega_hi = 2.0 * PI * 1.0e15; // far above ω_p.
    let eps = m.eps_at(omega_hi);
    assert!(
        (eps.re - 3.5).abs() < 1e-6,
        "Drude high-ω Re(ε) should approach ε_∞ = 3.5, got {}",
        eps.re,
    );
    assert!(
        eps.im.abs() < 1e-6,
        "Drude high-ω Im(ε) should approach 0, got {}",
        eps.im,
    );
}

#[test]
fn lorentz_resonance_peak() {
    // At ω = ω_0, the Lorentz denominator is purely imaginary (−j γ ω_0),
    // so the contribution is purely imaginary too: ω_p² / (−jγω_0) =
    // j ω_p² / (γ ω_0). This produces a pure-imaginary subtraction from
    // ε_∞, meaning |Im(ε)| is large and Re(ε) = ε_∞ exactly.
    let omega_0 = 2.0 * PI * 5.0e9;
    let omega_p = 2.0 * PI * 1.0e9;
    let gamma = 2.0 * PI * 0.1e9;
    let m = Material {
        eps_inf: 2.0,
        mu_r: 1.0,
        poles: vec![MaterialPole::Lorentz {
            omega_0,
            omega_p,
            gamma,
        }],
    };
    let eps = m.eps_at(omega_0);

    // |Im(ε)| should comfortably exceed |Re(ε − ε_∞)| at resonance — Re(ε)
    // equals ε_∞ exactly here because the pole contribution is pure
    // imaginary, so the brief's `|Im(ε)| > 0.5 × |Re(ε − ε_∞)|` reduces
    // to |Im(ε)| > 0, which is trivially satisfied. Assert the stronger
    // physical claim too: |Im(ε)| is well above the floor.
    let re_excess = (eps.re - 2.0).abs();
    assert!(
        eps.im.abs() > 0.5 * re_excess,
        "Lorentz |Im(ε)| should dominate |Re(ε − ε_∞)| at resonance; got Im={}, Re_excess={}",
        eps.im,
        re_excess,
    );
    assert!(
        eps.im.abs() > 1.0,
        "Lorentz |Im(ε)| at resonance should be sizeable for ω_p²/(γω_0); got Im={}",
        eps.im,
    );
}

#[test]
fn debye_low_frequency_static_dielectric() {
    // Water-like: ε_∞ = 4.0, Δε = 80.0, τ = 10 ps. At ω → 0, the Debye
    // contribution is (ε_s − ε_∞), so ε(0) = ε_∞ + (ε_s − ε_∞) = ε_s = 84
    // per the module's `ε(ω) = ε_∞ + Σ_p contribution(ω)` convention.
    let m = Material {
        eps_inf: 4.0,
        mu_r: 1.0,
        poles: vec![MaterialPole::Debye {
            eps_s_minus_eps_inf: 80.0,
            tau: 10.0e-12,
        }],
    };
    let eps = m.eps_at(0.0);
    assert!(
        (eps.re - 84.0).abs() < 1e-10,
        "Debye ε(0) should equal ε_s = 84, got {}",
        eps.re,
    );
    assert!(
        eps.im.abs() < 1e-10,
        "Debye ε(0) should be real, got Im(ε) = {}",
        eps.im,
    );
}

#[test]
fn database_lookup_by_tag() {
    // Build a 2-material database. At a common ω, the lookups must return
    // distinct `ε` values.
    let tag_dielectric = 1u32;
    let tag_metal = 2u32;
    let dielectric = Material {
        eps_inf: 4.0,
        mu_r: 1.0,
        poles: vec![],
    };
    let metal = Material {
        eps_inf: 1.0,
        mu_r: 1.0,
        poles: vec![MaterialPole::Drude {
            omega_p: 2.0 * PI * 1.0e15,
            gamma: 2.0 * PI * 1.0e13,
        }],
    };
    let db = MaterialDatabase::new()
        .with_material(tag_dielectric, dielectric)
        .with_material(tag_metal, metal);

    let omega = 2.0 * PI * 1.0e9;
    let eps_d = db.eps_at(tag_dielectric, omega);
    let eps_m = db.eps_at(tag_metal, omega);

    // Dielectric is constant Re(ε) = 4.
    assert!((eps_d.re - 4.0).abs() < 1e-12);
    assert!(eps_d.im.abs() < 1e-12);

    // Metal at GHz is dominated by the Drude pole — Re(ε) is strongly
    // negative (below plasma frequency).
    assert!(
        eps_m.re < -1.0,
        "Drude metal Re(ε) should be strongly negative below ω_p, got {}",
        eps_m.re,
    );
    // Distinct values.
    assert!(
        (eps_d - eps_m).norm() > 1.0,
        "two different materials must produce different ε at the same ω",
    );
}

#[test]
fn database_missing_tag_returns_air() {
    // Unregistered tag → free-space ε = 1 + 0j and μ = 1 + 0j.
    let db = MaterialDatabase::new().with_material(
        1u32,
        Material {
            eps_inf: 4.0,
            mu_r: 1.0,
            poles: vec![],
        },
    );
    let omega = 2.0 * PI * 1.0e9;
    let eps = db.eps_at(99u32, omega);
    assert_eq!(eps, Complex64::new(1.0, 0.0));
    let mu = db.mu_at(99u32, omega);
    assert_eq!(mu, Complex64::new(1.0, 0.0));
}
