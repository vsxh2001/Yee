//! Source helpers for the FDTD walking skeleton.
//!
//! Phase 2.0 ships a single point-source primitive: a Gaussian-in-time pulse
//! added (soft source) to a chosen cell of `E_z`. Hard sources, modal sources,
//! and lumped ports are Phase 2.1+ work.

use crate::grid::YeeGrid;

/// Add a Gaussian-time pulse to `E_z(i, j, k)`.
///
/// The injected value is `exp(-((t - t0) / sigma)²)` (a unit-amplitude soft
/// source). The caller controls the time stepping; this function simply
/// *adds* the source contribution to the existing field value.
///
/// # Panics
///
/// Panics if `(i, j, k)` is outside the bounds of `E_z`
/// (shape `[nx+1, ny+1, nz]`).
pub fn gaussian_pulse_ez(
    grid: &mut YeeGrid,
    i: usize,
    j: usize,
    k: usize,
    t: f64,
    t0: f64,
    sigma: f64,
) {
    assert!(
        sigma > 0.0 && sigma.is_finite(),
        "gaussian sigma must be positive and finite"
    );
    let arg = (t - t0) / sigma;
    let amplitude = (-arg * arg).exp();
    grid.ez[(i, j, k)] += amplitude;
}
