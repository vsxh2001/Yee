//! FEM-EM brick B3 (ADR-0153) — closure-sanity gate for the quasi-TEM
//! microstrip wave-port closures in [`yee_fem::microstrip_port`].
//!
//! This gate is **closure sanity only** (sub-second, no LU solve). It
//! checks that the analytic `(β(ω), e_t(x))` closure pair is internally
//! consistent and physically sane:
//!
//! 1. [`beta_matches_hammerstad_jensen`] — `β(2π·2 GHz)` equals
//!    `(ω/c)·sqrt(ε_eff)` computed independently from first principles
//!    (a real cross-check, not the closure asserted against itself).
//! 2. [`modal_e_t_is_ez_dominant_in_gap_and_decays_in_air`] — the modal
//!    field is `E_z`-dominant (`|E_z| > |E_x|`, `|E_z| > |E_y|`) at a
//!    sample point inside the trace↔ground gap, and its `|E_z|` decays
//!    at a point well above the trace in air.
//! 3. [`modal_self_inner_product_finite_nonzero`] — the L²
//!    self-inner-product `Σ A_face (e_t·e_t)` over a representative
//!    port-face point set is finite and non-zero (modal normalisation
//!    is sane), mirroring the WR-90
//!    `modal_self_inner_product_matches_orthonormalisation` diagnostic
//!    in `open_boundary_sweep.rs`.
//!
//! The **fidelity** of the analytic modal shape (does driving a line
//! with it recover the Hammerstad-Jensen ε_eff end-to-end?) is **B4**'s
//! gate and is deliberately NOT asserted here — no `sweep` /
//! `sweep_matrix` / LU solve runs in this file.
//!
//! References:
//! * `yee_layout::eps_eff` (validated by `crates/yee-layout`
//!   `tests/geo_002_hammerstad.rs`).
//! * `crates/yee-fem/tests/open_boundary_sweep_matrix.rs`
//!   `beta_te10` / `modal_e_t_te10` — the WR-90 closure template.

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_core::units::C0;
use yee_fem::microstrip_port::{beta_microstrip, microstrip_port, modal_e_t_microstrip};

// Representative FR-4 microstrip geometry (a ~50 Ω line on 1 mm FR-4).
// These match the kind of inputs `geo_002_hammerstad` exercises.
const W: f64 = 1.9e-3;
const H: f64 = 1.0e-3;
const EPS_R: f64 = 4.4;

// ---------------------------------------------------------------------
// Criterion 1 — β(ω) matches first-principles Hammerstad-Jensen.
// ---------------------------------------------------------------------

/// The port's `β(ω)` closure must equal `(ω/c)·sqrt(ε_eff(w,h,εr))`
/// computed from first principles. The expected value is built from
/// the *literal* `ε_eff` algebra here (not by re-calling
/// `yee_layout::eps_eff`), so this is a genuine cross-check of the
/// closure rather than a tautology.
#[test]
fn beta_matches_hammerstad_jensen() {
    let omega = 2.0 * PI * 2.0e9;

    // First-principles ε_eff (Schneider / Hammerstad-Jensen form, the
    // same algebra `yee_layout::eps_eff` implements):
    //   ε_eff = (εr+1)/2 + (εr−1)/2 · (1 + 12 h/W)^(−1/2)
    let eps_eff_expected =
        (EPS_R + 1.0) / 2.0 + (EPS_R - 1.0) / 2.0 * (1.0 + 12.0 * H / W).powf(-0.5);
    let beta_expected = (omega / C0) * eps_eff_expected.sqrt();

    // Via the free function …
    let beta_fn = beta_microstrip(W, H, EPS_R, omega);
    // … and via the PortDefinition's stored closure (modes[0]).
    let port = microstrip_port(W, H, EPS_R);
    let beta_closure = (port.modes[0].beta_mode)(omega);

    assert!(
        (beta_fn - beta_expected).abs() < 1e-9,
        "beta_microstrip = {beta_fn}, first-principles expected = {beta_expected}, \
         |diff| = {:e} exceeds 1e-9",
        (beta_fn - beta_expected).abs()
    );
    assert!(
        (beta_closure - beta_expected).abs() < 1e-9,
        "port closure β = {beta_closure}, first-principles expected = {beta_expected}, \
         |diff| = {:e} exceeds 1e-9",
        (beta_closure - beta_expected).abs()
    );

    // Cross-check the independently-computed ε_eff against the library
    // function as well, so a regression in `yee_layout::eps_eff` that
    // happened to match the closure would still be caught here.
    let eps_eff_lib = yee_layout::eps_eff(W, H, EPS_R);
    assert!(
        (eps_eff_lib - eps_eff_expected).abs() < 1e-12,
        "yee_layout::eps_eff = {eps_eff_lib} disagrees with first-principles \
         {eps_eff_expected}"
    );

    // Sanity: quasi-TEM ε_eff lies between 1 (all-air) and εr (all-
    // dielectric), so β is real, positive, and above the free-space
    // wavenumber.
    assert!(
        beta_expected > omega / C0,
        "quasi-TEM β = {beta_expected} should exceed the free-space wavenumber \
         k0 = {} (ε_eff > 1)",
        omega / C0
    );
}

// ---------------------------------------------------------------------
// Criterion 2 — modal_e_t is E_z-dominant in the gap and decays in air.
// ---------------------------------------------------------------------

/// On a `(x, z)` port face the dominant quasi-TEM component is the
/// substrate-normal `E_z` in the trace↔ground gap. Sample the modal
/// shape:
///
/// * in the gap (between ground `z = 0` and trace top `z = sub_h`,
///   within the trace `x`-window) → `|E_z| > |E_x|` and `|E_z| > |E_y|`;
/// * well above the trace in air (`z ≫ sub_h`) → smaller `|E_z|` than
///   in the gap (the fringing field decays).
#[test]
fn modal_e_t_is_ez_dominant_in_gap_and_decays_in_air() {
    // A point in the middle of the trace↔ground gap, at a representative
    // trace-centred x (the v1 shape is x-uniform, so any x is in-window).
    let x_trace = 2.0e-3;
    let p_gap = Vector3::new(x_trace, 0.0, H / 2.0);
    let e_gap = modal_e_t_microstrip(H, p_gap);

    assert!(
        e_gap.z.abs() > e_gap.x.abs(),
        "in-gap field should be E_z-dominant: |E_z| = {} not > |E_x| = {}",
        e_gap.z.abs(),
        e_gap.x.abs()
    );
    assert!(
        e_gap.z.abs() > e_gap.y.abs(),
        "in-gap field should be E_z-dominant: |E_z| = {} not > |E_y| = {}",
        e_gap.z.abs(),
        e_gap.y.abs()
    );
    assert!(
        e_gap.z.abs() > 0.0,
        "in-gap |E_z| must be non-zero, got {}",
        e_gap.z.abs()
    );

    // A point well above the trace in air: z = sub_h + 3·sub_h = 4·sub_h
    // (three air-decay lengths above the trace top).
    let p_air = Vector3::new(x_trace, 0.0, 4.0 * H);
    let e_air = modal_e_t_microstrip(H, p_air);

    assert!(
        e_air.z.abs() < e_gap.z.abs(),
        "field above the trace in air (|E_z| = {}) should decay below the \
         in-gap value (|E_z| = {})",
        e_air.z.abs(),
        e_gap.z.abs()
    );
    // The decay is exp(−(z−sub_h)/sub_h); three lengths up → e^{−3}.
    let expected_air = (-3.0_f64).exp();
    assert!(
        (e_air.z - expected_air).abs() < 1e-9,
        "air-tail |E_z| = {} should equal exp(−3) = {expected_air} \
         (one-substrate-height decay)",
        e_air.z
    );
}

// ---------------------------------------------------------------------
// Criterion 3 — modal self-inner-product is finite and non-zero.
// ---------------------------------------------------------------------

/// The L² modal self-inner-product `Σ_face A_face · (e_t · e_t)` over a
/// representative port-face point set must be finite and non-zero —
/// otherwise the solver's modal normalisation `S = ⟨E,e_t⟩/⟨e_t,e_t⟩`
/// would divide by zero or NaN. Mirrors the WR-90
/// `modal_self_inner_product_matches_orthonormalisation` diagnostic in
/// `open_boundary_sweep.rs`, but using a representative `(x, z)`
/// port-face sample grid (no mesh / solver constructed → no solve).
#[test]
fn modal_self_inner_product_finite_nonzero() {
    // Representative port-face geometry: a 4 mm (x) × 4 mm (z) box face
    // sampled on an 8 × 8 cell grid (cell-centre Riemann quadrature),
    // matching the kind of mesh a layered_microstrip_mesh(4e-3, 4e-3,
    // …) would present at a y-end cap.
    let box_w = 4.0e-3;
    let box_h = 4.0e-3;
    let nx = 8usize;
    let nz = 8usize;
    let dx = box_w / nx as f64;
    let dz = box_h / nz as f64;
    let area_per_cell = dx * dz;

    let mut inner = 0.0_f64;
    for ix in 0..nx {
        for iz in 0..nz {
            let x = (ix as f64 + 0.5) * dx;
            let z = (iz as f64 + 0.5) * dz;
            let e_t = modal_e_t_microstrip(H, Vector3::new(x, 0.0, z));
            inner += area_per_cell * e_t.dot(&e_t);
        }
    }

    assert!(
        inner.is_finite(),
        "modal self-inner-product must be finite, got {inner}"
    );
    assert!(
        inner > 0.0,
        "modal self-inner-product must be non-zero (modal normalisation \
         would otherwise divide by zero), got {inner}"
    );

    // Lower-bound sanity: the gap region (z ∈ [0, sub_h], |E_z| = 1)
    // alone contributes at least box_w · sub_h to ∫|e_t|² dS, so the
    // total must comfortably exceed a small floor. (Not an upper bound
    // — fidelity is a B4 concern.)
    let gap_floor = box_w * H * 0.5;
    assert!(
        inner > gap_floor,
        "modal self-inner-product {inner} should exceed the gap-region \
         floor {gap_floor}"
    );
}
