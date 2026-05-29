//! `geo-002` — Hammerstad-Jensen microstrip synthesis gate.
//!
//! The published FR-4 50 Ω reference (`ε_r = 4.4`, `h = 1.6 mm`) is
//! `W ≈ 3.0 mm`, `ε_eff ≈ 3.3` (Pozar §3.8 / Hammerstad-Jensen 1980). Both
//! must land within ±5%. These are ground-truth published numbers — if the
//! width misses, the A/B (thin/wide-line) branch selection is wrong; do NOT
//! loosen the tolerance.

use yee_layout::{eps_eff, microstrip_width};

/// Relative tolerance for the published HJ reference.
const REL_TOL: f64 = 0.05;

#[test]
fn fr4_50ohm_width() {
    let z0 = 50.0;
    let eps_r = 4.4;
    let h = 1.6e-3;
    let w = microstrip_width(z0, eps_r, h);
    let want = 3.0e-3;
    let rel = (w - want).abs() / want;
    assert!(
        rel < REL_TOL,
        "microstrip_width(50, 4.4, 1.6mm) = {w:.6e} m, expected ≈ {want:.6e} m (rel err {:.2}% > 5%)",
        rel * 100.0
    );
}

#[test]
fn fr4_50ohm_eps_eff() {
    let z0 = 50.0;
    let eps_r = 4.4;
    let h = 1.6e-3;
    let w = microstrip_width(z0, eps_r, h);
    let ee = eps_eff(w, h, eps_r);
    let want = 3.3;
    let rel = (ee - want).abs() / want;
    assert!(
        rel < REL_TOL,
        "eps_eff(W, 1.6mm, 4.4) = {ee:.4}, expected ≈ {want:.4} (rel err {:.2}% > 5%)",
        rel * 100.0
    );
}

/// The thin-line branch is the one that must fire for FR-4 50 Ω (W/h ≈ 1.9 < 2).
/// A regression that flips to the wide-line branch would push the width well
/// past the ±5% gate; this records the expected ratio explicitly.
#[test]
fn fr4_50ohm_is_thin_line_branch() {
    let w = microstrip_width(50.0, 4.4, 1.6e-3);
    let w_over_h = w / 1.6e-3;
    assert!(
        w_over_h < 2.0,
        "FR-4 50 Ω should be the thin-line (W/h < 2) branch, got W/h = {w_over_h:.4}"
    );
}
