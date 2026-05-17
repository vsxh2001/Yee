//! Phase 1.3.1.1 step 3 validation gate — WR-90 TE10 numerical eigensolve.
//!
//! Builds a coarse structured triangular mesh of the WR-90 cross-section
//! (a × b = 22.86 mm × 10.16 mm, air-filled) and drives
//! [`NumericalCrossSection::solve`] at 10 GHz. The numerical β must
//! match the analytic [`RectangularWaveguideTe10::beta`] within 1 %,
//! and `Z_w` within 5 %. The 0.1 % / 1 % gates from the design spec
//! are deferred to a refined mesh in Phase 1.3.1.1 step 5.

use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::TriMesh2D;
use yee_mom::ports::{NumericalCrossSection, RectangularWaveguideTe10};

/// Structured `nx × ny` quad-grid WR-90 mesh. Each quad is split along
/// the `(low-x, low-y) → (high-x, high-y)` diagonal into two CCW
/// triangles. All triangles share material tag 0 (air fill).
fn rectangular_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([a * (i as f64) / (nx as f64), b * (j as f64) / (ny as f64)]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;
    let mut triangles = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v11 = idx(i + 1, j + 1);
            let v01 = idx(i, j + 1);
            triangles.push([v00, v10, v11]);
            triangles.push([v00, v11, v01]);
        }
    }
    TriMesh2D::new(vertices, triangles, None, None).unwrap()
}

#[test]
fn eigensolver_wr90_te10_beta_within_1_percent() {
    // WR-90: a = 22.86 mm × b = 10.16 mm, air. 6×6 quads → 72 triangles.
    let a = 22.86e-3;
    let b = 10.16e-3;
    let freq_hz = 10.0e9;

    let mesh = rectangular_mesh(a, b, 6, 6);
    let mut eps_r = HashMap::new();
    eps_r.insert(0u32, Complex64::new(1.0, 0.0));
    let mut mu_r = HashMap::new();
    mu_r.insert(0u32, Complex64::new(1.0, 0.0));

    let mut mode = NumericalCrossSection::new(mesh, eps_r, mu_r);
    mode.solve(freq_hz)
        .expect("WR-90 eigensolve should succeed");

    let beta_num = mode.beta.expect("β cached after solve").re;
    let analytic = RectangularWaveguideTe10 { a, b, eps_r: 1.0 };
    let beta_analytic = analytic.beta(freq_hz);
    let rel_err = (beta_num - beta_analytic).abs() / beta_analytic;
    // Diagnostic always-emit: the values land in `--nocapture` output so
    // the validation log records the achieved error.
    eprintln!(
        "WR-90 TE10 β @ 10 GHz: numerical {beta_num:.6} rad/m, \
         analytic {beta_analytic:.6} rad/m, rel err {rel_err:.6}"
    );
    assert!(
        rel_err < 0.01,
        "WR-90 TE10 β: numerical {beta_num} vs analytic {beta_analytic} \
         (rel err {rel_err:.4}); want < 1 % on a 6×6-quad mesh"
    );
}

#[test]
fn eigensolver_wr90_te10_zw_within_5_percent() {
    let a = 22.86e-3;
    let b = 10.16e-3;
    let freq_hz = 10.0e9;

    let mesh = rectangular_mesh(a, b, 6, 6);
    let mut eps_r = HashMap::new();
    eps_r.insert(0u32, Complex64::new(1.0, 0.0));
    let mut mu_r = HashMap::new();
    mu_r.insert(0u32, Complex64::new(1.0, 0.0));

    let mut mode = NumericalCrossSection::new(mesh, eps_r, mu_r);
    mode.solve(freq_hz)
        .expect("WR-90 eigensolve should succeed");

    let zw_num = mode.z_w.expect("Z_w cached after solve").norm();
    let analytic = RectangularWaveguideTe10 { a, b, eps_r: 1.0 };
    let zw_analytic = analytic.wave_impedance(freq_hz);
    let rel_err = (zw_num - zw_analytic).abs() / zw_analytic;
    assert!(
        rel_err < 0.05,
        "WR-90 TE10 Z_w: numerical {zw_num} Ω vs analytic {zw_analytic} Ω \
         (rel err {rel_err:.4}); want < 5 % on a 6×6-quad mesh"
    );
}
