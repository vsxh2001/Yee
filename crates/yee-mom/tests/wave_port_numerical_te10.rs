//! Phase 1.3.1.1 step 7 validation gate — numerical 2-D wave-port RHS
//! matches the closed-form TE10 RHS on a WR-90 cross-section.
//!
//! The wired Numerical2D arm of [`WavePort::rhs`] samples the cached
//! Nedelec eigenmode at each port-edge midpoint (via
//! [`NumericalCrossSection::e_tangential_at`]) and projects the modal
//! `E_t` field onto the edge's tangent unit vector. On the WR-90 air-
//! filled rectangular waveguide cross-section the dominant numerical
//! mode is TE10 (`E_y = sin(π x / a)`, `E_x = 0`); the gate below
//! checks that the resulting RHS agrees with the closed-form TE10 RHS
//! (the [`ModalDistribution::Te10`] arm) within 1 % L2 once both are
//! renormalized to unit L2.
//!
//! **Mesh setup.** The RWG port face is a 3-D triangulated patch in the
//! z = 0 plane spanning the WR-90 cross-section (a × b = 22.86 mm ×
//! 10.16 mm). Triangle tags alternate by column index so the vertical
//! edges between columns become port edges (port_tag = 1) by the
//! standard RWG tag-mismatch rule. With `nx` columns there are `nx - 1`
//! vertical interfaces, each contributing `ny` y-aligned port edges at
//! distinct x positions — enough to sample the modal profile across
//! the cross-section without degenerating to a single-x line.
//!
//! The 2-D cross-section mesh used by [`NumericalCrossSection`] uses
//! the SAME (x, y) coordinate convention so the port-edge midpoints
//! land inside the eigensolve's domain.

use nalgebra::Vector3;
use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::{TriMesh, TriMesh2D};
use yee_mom::__internal::{RwgBasis, build_basis, wave_port_rhs_for_test};
use yee_mom::ports::{ModalDistribution, NumericalCrossSection, RectangularWaveguideTe10};

/// Build a 3-D RWG mesh in the z = 0 plane spanning `[0, a] × [0, b]`,
/// with alternating-column triangle tags so the vertical inter-column
/// edges become port edges (port_tag = 1).
fn rwg_mesh_alternating_columns(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push(Vector3::new(
                a * (i as f64) / (nx as f64),
                b * (j as f64) / (ny as f64),
                0.0,
            ));
        }
    }
    let idx = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    let mut triangles = Vec::with_capacity(2 * nx * ny);
    let mut tags = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v11 = idx(i + 1, j + 1);
            let v01 = idx(i, j + 1);
            // CCW: v00 -> v10 -> v11, v00 -> v11 -> v01
            triangles.push([v00, v10, v11]);
            triangles.push([v00, v11, v01]);
            // Tag by column parity: even columns -> 1, odd columns -> 2.
            // Tag-mismatch rule makes vertical inter-column edges port edges
            // with port_tag = min(1, 2) = 1.
            let column_tag = if i.is_multiple_of(2) { 1u32 } else { 2u32 };
            tags.push(column_tag);
            tags.push(column_tag);
        }
    }
    TriMesh::new(vertices, triangles, tags).unwrap()
}

/// Build the 2-D cross-section mesh in the SAME (x, y) frame as the
/// RWG port mesh above so the eigenmode interpolant lines up with the
/// port-edge midpoints.
fn cross_section_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
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

fn l2_norm(b: &faer::Mat<Complex64>) -> f64 {
    let mut s = 0.0f64;
    for i in 0..b.nrows() {
        s += b[(i, 0)].norm_sqr();
    }
    s.sqrt()
}

fn l2_diff(b1: &faer::Mat<Complex64>, b2: &faer::Mat<Complex64>) -> f64 {
    let mut s = 0.0f64;
    assert_eq!(b1.nrows(), b2.nrows());
    for i in 0..b1.nrows() {
        let d = b1[(i, 0)] - b2[(i, 0)];
        s += d.norm_sqr();
    }
    s.sqrt()
}

fn scale(b: &faer::Mat<Complex64>, s: Complex64) -> faer::Mat<Complex64> {
    let n = b.nrows();
    let mut out = faer::Mat::<Complex64>::zeros(n, 1);
    for i in 0..n {
        out[(i, 0)] = b[(i, 0)] * s;
    }
    out
}

#[test]
fn wave_port_numerical_matches_te10_within_1_percent() {
    // WR-90 air-filled at 10 GHz. 6×6 mesh: 36 quads → 72 triangles
    // → 84 interior edges in the eigensolve (matches the existing
    // `eigensolver_wr90` gate's mesh density). The RWG port face uses
    // the same 6×6 quad grid so the port-edge midpoints line up with
    // the cross-section mesh's natural vertex spacing.
    let a = 22.86e-3;
    let b = 10.16e-3;
    let freq_hz = 10.0e9;
    let nx = 6;
    let ny = 6;

    // 1) Build the numerical cross-section and solve the eigenproblem.
    let xs_mesh = cross_section_mesh(a, b, nx, ny);
    let mut eps_r = HashMap::new();
    eps_r.insert(0u32, Complex64::new(1.0, 0.0));
    let mut mu_r = HashMap::new();
    mu_r.insert(0u32, Complex64::new(1.0, 0.0));
    let mut mode = NumericalCrossSection::new(xs_mesh, eps_r, mu_r);
    mode.solve(freq_hz)
        .expect("WR-90 numerical eigensolve should succeed");

    // 2) Sample at cross-section centroid to fix the global sign so
    // `E_y > 0` matches the analytic TE10 convention. The dense
    // eigensolve fixes the largest-magnitude DoF positive but the
    // physical sign on `E_y` is not pinned by that convention alone.
    let e_center = mode.e_tangential_at(0.5 * a, 0.5 * b);
    let sign = if e_center[1] >= 0.0 { 1.0 } else { -1.0 };
    if sign < 0.0 {
        // Flip the cached profile in-place to align the convention.
        let profile_flipped: Vec<Complex64> = mode
            .mode_profile
            .as_ref()
            .unwrap()
            .iter()
            .map(|z| -*z)
            .collect();
        mode.mode_profile = Some(profile_flipped);
    }

    // 3) Build the RWG port-face mesh and basis. Alternating-column
    // tags produce y-aligned port edges at multiple x positions, so
    // the RHS samples the modal profile across the cross-section
    // rather than at a single point.
    let port_mesh = rwg_mesh_alternating_columns(a, b, nx, ny);
    let basis: RwgBasis = build_basis(&port_mesh).unwrap();
    let port_indices: Vec<usize> = basis.port_basis_indices(1).collect();
    assert!(
        !port_indices.is_empty(),
        "RWG port mesh must produce at least one port edge"
    );
    eprintln!(
        "WR-90 RHS gate: {} port edges at {} y-positions × {} x-interfaces",
        port_indices.len(),
        ny,
        nx - 1
    );

    // 4) Build both wave-ports and compute their RHS at 10 GHz via
    // the `__internal::wave_port_rhs_for_test` helper, which wraps the
    // crate-private [`yee_mom::ports::Port::rhs`] for integration-test
    // consumers (the trait itself is `pub(crate)`).
    let b_num_raw = wave_port_rhs_for_test(
        &basis,
        1,
        Complex64::new(1.0, 0.0),
        1.0,
        ModalDistribution::Numerical2D(Box::new(mode)),
        freq_hz,
    );
    let b_te10 = wave_port_rhs_for_test(
        &basis,
        1,
        Complex64::new(1.0, 0.0),
        1.0,
        ModalDistribution::Te10(RectangularWaveguideTe10 { a, b, eps_r: 1.0 }),
        freq_hz,
    );

    let norm_te10 = l2_norm(&b_te10);
    let norm_num = l2_norm(&b_num_raw);
    assert!(
        norm_te10 > 0.0,
        "Analytic TE10 RHS L2 norm must be non-zero"
    );
    assert!(
        norm_num > 0.0,
        "Numerical RHS L2 norm must be non-zero (eigensolve produced a vanishing mode?)"
    );

    // 5) Normalize both to unit L2 to remove the eigenvector global-
    // scale freedom: the numerical eigenvector is determined only up
    // to scale, while the analytic profile is the canonical sin(πx/a)
    // ∈ [0, 1]. The comparison metric is the L2 difference of the
    // unit-normalized RHS vectors.
    let b_num_unit = scale(&b_num_raw, Complex64::new(1.0 / norm_num, 0.0));
    let b_te10_unit = scale(&b_te10, Complex64::new(1.0 / norm_te10, 0.0));
    let diff = l2_diff(&b_num_unit, &b_te10_unit);
    let rel = diff; // unit-normalized denominator is 1.

    eprintln!(
        "WR-90 wave-port RHS L2 agreement: numerical-vs-analytic = {:.6} \
         (want < 0.01); n_port_edges = {}",
        rel,
        port_indices.len()
    );

    assert!(
        rel < 0.01,
        "WR-90 numerical-vs-analytic wave-port RHS L2 disagreement: \
         {:.4} (>= 1 % budget) on a {}×{} mesh",
        rel,
        nx,
        ny
    );
}
