//! Track WWWWWWW — TEM-mode smoothed RHS verification gate.
//!
//! ## Why this file exists
//!
//! Track TTTTTTT (`tests/mom_002_port_edge_diagnostic.rs`) ran a P1
//! port-excitation probe on the IIIIIII-reframed mom-002 mesh
//! (`L = 82 mm`, `w = 2.94 mm`, centered uniform `82 × 16`,
//! Sommerfeld kernel at 1 GHz) and found:
//!
//! * the production [`yee_mom::__internal::delta_gap_rhs_for_test`]
//!   excitation drives a piecewise alternating per-edge current
//!   pattern at the port — `|i|` on longitudinal-edge basis functions
//!   ~`1.4` vs `~0.06-0.26` on diagonal-edge basis functions —
//!   instead of a smoothed quasi-TEM transverse mode;
//! * the per-y-bucket transverse-profile sum deviates from a Maxwell
//!   `1/√(1 − (2y/w)²)` envelope by **+580 %** at the most-deviant
//!   bucket.
//!
//! Track QQQQQQQ (`tests/mom_002_beta_eigenmode_probe.rs`) earlier
//! exonerated the kernel — `ε_eff_solver = 3.385` vs Hammerstad-Jensen
//! `3.32` (+1.83 % error). So the residual `|Im(Z)| ≈ 674 Ω`
//! capacitive bias is a **port-excitation modeling** issue, not a
//! kernel issue.
//!
//! Track WWWWWWW's fix replaces the uniform delta-gap RHS with a
//! Maxwell-envelope-weighted smoothed RHS via
//! [`yee_mom::__internal::tem_smoothed_rhs_for_test`] — every port
//! basis function gets `b_k = V · length_k · w_TEM(y_k)` with
//! `w_TEM(y) = sqrt(2 / (π · (1 − (2 y / w)²)))`, and the same
//! weighting is applied symmetrically to the port-current extraction.
//!
//! ## What this gate verifies
//!
//! At the production `82 × 16` strip mesh (same FR-4 substrate, same
//! 1 GHz probe, same centered port placement as the yee-validation
//! headline gate) we
//!
//! 1. Solve the production `Z` matrix once and apply both RHS variants
//!    to the **same** LU factorisation.
//! 2. Walk the port basis functions, bucket by transverse y, sum
//!    `length_k · |i_k|` per bucket, normalise against the centre
//!    bucket, and compute the Max deviation from the analytic Maxwell
//!    envelope (the metric TTTTTTT P1 reported).
//! 3. Report the per-RHS `Z_in = V_port / I_port` so the headline
//!    gate's tripwire shift can be audited against this measurement.
//! 4. Assert the TEM-smoothed RHS reduces the Max deviation by ≥ 5×
//!    versus the delta-gap RHS (the brief's "by at least 5×"
//!    escape-hatch bound). Measured `8.32×` on the `82 × 16` mesh
//!    (`579.82 %` → `69.67 %`).
//!
//! Wall-time budget: ~4 min on a modern desktop, dominated by the
//! dense `3838 × 3838` Sommerfeld-fill + LU. The headline gate
//! (`yee_validation::tests::mom_002_headline_gate_passes`) re-runs
//! the same kernel at the same mesh density so the marginal cost is
//! comparable.
//!
//! ## References
//!
//! * `tests/mom_002_port_edge_diagnostic.rs` (Track TTTTTTT) — P1
//!   diagnosis + Maxwell-envelope deviation metric.
//! * `tests/mom_002_beta_eigenmode_probe.rs` (Track QQQQQQQ) — kernel
//!   exoneration via β-from-Z extraction.
//! * R. F. Harrington, *Time-Harmonic Electromagnetic Fields*,
//!   McGraw-Hill, 1961, §5.5 (Maxwell `1/√(1 − (2y/w)²)`
//!   edge-singularity envelope on a thin strip).
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7 (microstrip
//!   `Z_0` / `ε_eff` reference).

use faer::linalg::solvers::{PartialPivLu, Solve};
use nalgebra::Vector3;
use num_complex::Complex64;
use std::collections::BTreeMap;
use yee_mesh::TriMesh;
use yee_mom::__internal::{
    MultilayerGreens, build_basis, delta_gap_rhs_for_test, impedance_matrix_for_test,
    tem_smoothed_rhs_for_test,
};

// FR-4 / 1 GHz canonical microstrip parameters — match the IIIIIII
// reframed `yee-validation::MOM_002_*` constants exactly so this gate
// and the headline `mom_002_headline_gate_passes` test run on the
// same geometry.
const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const STRIP_W_M: f64 = 2.94e-3;
const STRIP_L_M: f64 = 82.0e-3;
const F_HZ: f64 = 1.0e9;

// Production headline-gate mesh density. Tried `40 × 8` first for
// fast iteration; on that coarser mesh the delta-gap RHS produces
// only `~54 %` Maxwell-envelope deviation and the TEM weighting is
// roughly neutral (the alternating-edge pattern isn't yet sharp
// enough to show the fix). The production `82 × 16` mesh reproduces
// TTTTTTT's `+580 %` baseline cleanly and shows an `8.32 ×` reduction
// under the TEM-smoothed RHS, which is the verdict-bearing
// configuration.
const N_LENGTH: usize = 82;
const N_WIDTH: usize = 16;

// Sommerfeld kernel parameters — match the production headline gate.
const N_DCIM_IMAGES: usize = 5;
const N_SW_POLES: usize = 1;

/// Build a centered-port mesh identical in structure to the
/// `yee_validation::mom_002_strip_mesh_with_spacing` builder under
/// `StripSpacing::Uniform`, parameterised by `(n_length, n_width)`.
/// Lifted from `tests/mom_002_port_edge_diagnostic.rs` so this gate
/// stays self-contained.
fn build_centered_strip_mesh(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
) -> TriMesh {
    assert!(
        n_length >= 4 && n_length.is_multiple_of(2),
        "n_length must be even and ≥ 4"
    );
    assert!(n_width >= 1, "n_width must be ≥ 1");

    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(nx * ny);
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;

    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..=n_width {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }

    let port_left = n_length / 2 - 1;
    let port_right = n_length / 2;
    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_length * n_width);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_length * n_width);
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == port_left {
                1
            } else if i == port_right {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }

    TriMesh::new(vertices, triangles, tags).expect("strip mesh invariants")
}

/// Delta-gap port-current projection: `I_port = Σ length_k · i_k`.
/// Mirrors `DeltaGapPort::port_current` for the test forensics; we
/// can't call the production method directly because `Port` is
/// crate-private.
fn port_current_delta_gap(
    basis: &yee_mom::__internal::RwgBasis,
    i_vec: &faer::Mat<Complex64>,
) -> Complex64 {
    let mut total = Complex64::new(0.0, 0.0);
    for k in basis.port_basis_indices(1) {
        total += Complex64::new(basis.edges[k].length, 0.0) * i_vec[(k, 0)];
    }
    total
}

/// TEM-smoothed port-current projection: `I_port = Σ length_k ·
/// w_TEM(y_k) · i_k` with the same Maxwell envelope used on the way
/// in. Mirrors `TemSmoothedPort::port_current` for the test
/// forensics.
fn port_current_tem_smoothed(
    basis: &yee_mom::__internal::RwgBasis,
    mesh: &TriMesh,
    i_vec: &faer::Mat<Complex64>,
    strip_width_m: f64,
) -> Complex64 {
    let mut total = Complex64::new(0.0, 0.0);
    for k in basis.port_basis_indices(1) {
        let edge = basis.edges[k];
        let p0 = mesh.vertices[edge.v0 as usize];
        let p1 = mesh.vertices[edge.v1 as usize];
        let y = 0.5 * (p0.y + p1.y);
        let u = (2.0 * y / strip_width_m).abs().min(1.0 - 1e-3);
        let w_tem = (2.0 / (std::f64::consts::PI * (1.0 - u * u))).sqrt();
        total += Complex64::new(edge.length * w_tem, 0.0) * i_vec[(k, 0)];
    }
    total
}

/// Compute the per-y-bucket transverse profile and the Max fractional
/// deviation from the analytic Maxwell `1/√(1 − (2y/w)²)` envelope.
///
/// Mirrors the metric reported by Track TTTTTTT's P1 probe so this
/// gate's "before/after" numbers compare apples-to-apples with the
/// `+580 %` baseline diagnosed there. Returns `(max_dev_fraction,
/// per_bucket_sums)`.
fn maxwell_envelope_max_deviation(
    basis: &yee_mom::__internal::RwgBasis,
    mesh: &TriMesh,
    i_vec: &faer::Mat<Complex64>,
    width_m: f64,
    n_width: usize,
) -> (f64, Vec<(f64, f64)>) {
    let port_indices: Vec<usize> = basis.port_basis_indices(1).collect();

    let dy = width_m / (n_width as f64);
    let bucket_y = |y: f64| -> i64 { ((y + width_m / 2.0) / (0.5 * dy)).round() as i64 };

    let mut per_y: BTreeMap<i64, f64> = BTreeMap::new();
    for &k in &port_indices {
        let edge = basis.edges[k];
        let v0 = mesh.vertices[edge.v0 as usize];
        let v1 = mesh.vertices[edge.v1 as usize];
        let ymid = 0.5 * (v0.y + v1.y);
        let mag = i_vec[(k, 0)].norm();
        let key = bucket_y(ymid);
        *per_y.entry(key).or_insert(0.0) += edge.length * mag;
    }

    // Maxwell envelope, regularised away from the singularity.
    let envelope = |y: f64| -> f64 {
        let u = (2.0 * y / width_m).abs().min(1.0 - 1e-3);
        1.0 / (1.0 - u * u).sqrt()
    };

    // Find the bucket closest to y = 0 for normalisation.
    let mut centre_key = 0_i64;
    let mut centre_dist = f64::INFINITY;
    for &k in per_y.keys() {
        let y = (k as f64) * 0.5 * dy - width_m / 2.0;
        if y.abs() < centre_dist {
            centre_dist = y.abs();
            centre_key = k;
        }
    }
    let centre_sum = per_y.get(&centre_key).copied().unwrap_or(1.0).max(1e-30);
    let centre_y = (centre_key as f64) * 0.5 * dy - width_m / 2.0;
    let centre_env = envelope(centre_y).max(1e-30);

    let mut max_dev: f64 = 0.0;
    let mut per_bucket: Vec<(f64, f64)> = Vec::with_capacity(per_y.len());
    for (&key, &sum_li) in &per_y {
        let y = (key as f64) * 0.5 * dy - width_m / 2.0;
        let solver_norm = sum_li / centre_sum;
        let ref_norm = envelope(y) / centre_env;
        per_bucket.push((y, sum_li));
        if ref_norm > 0.0 {
            let ratio = solver_norm / ref_norm;
            let dev = (ratio - 1.0).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
    }

    (max_dev, per_bucket)
}

/// Track WWWWWWW P1-fix gate: assemble the production `Z` once at the
/// `40 × 8` reduced-cost mom-002 mesh, apply the delta-gap and
/// TEM-smoothed RHS variants in turn, and require the TEM-smoothed
/// path to drop the Max Maxwell-envelope deviation by ≥ 5×.
///
/// Per the brief: "If the TEM weighting doesn't reduce the
/// port-current oscillation deviation (TTTTTTT measured +580 %) by at
/// least 5×, surface specifically and stop."
#[test]
fn tem_smoothed_rhs_reduces_port_oscillation() {
    eprintln!("--- Track WWWWWWW: TEM-smoothed RHS verification gate ---");
    eprintln!(
        "  Mesh: L = {} mm, w = {} mm, {} × {}, centered port, uniform y",
        STRIP_L_M * 1e3,
        STRIP_W_M * 1e3,
        N_LENGTH,
        N_WIDTH,
    );

    let mesh = build_centered_strip_mesh(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH);
    let basis = build_basis(&mesh).expect("RwgBasis build");
    eprintln!("  RWG basis count: {}", basis.n_basis());

    // Build the production Sommerfeld kernel once — same numerics as
    // the headline gate at this frequency.
    let greens = MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ,
        N_DCIM_IMAGES,
        N_SW_POLES,
    );

    let t0 = std::time::Instant::now();
    let z = impedance_matrix_for_test(&basis, &greens);
    let lu = PartialPivLu::new(z.as_ref());
    eprintln!(
        "  Z fill + LU factorisation: {:.2} s",
        t0.elapsed().as_secs_f64()
    );

    // Baseline: delta-gap RHS, identical to the Phase 1.0 path.
    let b_dg = delta_gap_rhs_for_test(&basis, 1);
    let i_dg = lu.solve(b_dg.as_ref());
    let (dev_dg, profile_dg) =
        maxwell_envelope_max_deviation(&basis, &mesh, &i_dg, STRIP_W_M, N_WIDTH);

    // Port current (delta-gap projection): Σ length_k · i_k. Z_in =
    // V_port / I_port. The baseline production headline gate reports
    // ≈ 674 Ω for the IIIIIII reframed strip mesh (Track IIIIIII
    // measurement) — sanity check we're hitting the same neighbourhood.
    let i_port_dg = port_current_delta_gap(&basis, &i_dg);
    let z_in_dg = num_complex::Complex64::new(1.0, 0.0) / i_port_dg;

    // Fix: TEM-mode-weighted RHS via Track WWWWWWW's TemSmoothedPort.
    let b_tem = tem_smoothed_rhs_for_test(&basis, 1, STRIP_W_M);
    let i_tem = lu.solve(b_tem.as_ref());
    let (dev_tem, profile_tem) =
        maxwell_envelope_max_deviation(&basis, &mesh, &i_tem, STRIP_W_M, N_WIDTH);

    // Port current (TEM-smoothed Galerkin projection): same
    // `w_TEM(y_k) · length_k` weighting on the way out as on the way
    // in. Z_in = V_port / I_port retains the inner-product structure.
    let i_port_tem = port_current_tem_smoothed(&basis, &mesh, &i_tem, STRIP_W_M);
    let z_in_tem = num_complex::Complex64::new(1.0, 0.0) / i_port_tem;

    eprintln!();
    eprintln!(
        "  Z_in @ delta-gap (baseline):    {:+10.4} + j{:+10.4} Ω, |Z| = {:.4} Ω",
        z_in_dg.re,
        z_in_dg.im,
        z_in_dg.norm()
    );
    eprintln!(
        "  Z_in @ TEM-smoothed (P1 fix):   {:+10.4} + j{:+10.4} Ω, |Z| = {:.4} Ω",
        z_in_tem.re,
        z_in_tem.im,
        z_in_tem.norm()
    );

    // Keep profiles around so they show up in `cargo test -- --nocapture`
    // diagnostic runs. Drop the explicit prints to keep the default
    // test output tight; the headline numbers above are sufficient
    // for the gate verdict.
    let _ = (profile_dg, profile_tem);

    eprintln!();
    eprintln!(
        "  Max Maxwell-envelope deviation @ delta-gap (baseline): {:.2} %",
        dev_dg * 100.0
    );
    eprintln!(
        "  Max Maxwell-envelope deviation @ TEM-smoothed (P1 fix): {:.2} %",
        dev_tem * 100.0
    );
    let reduction = if dev_tem > 0.0 {
        dev_dg / dev_tem
    } else {
        f64::INFINITY
    };
    eprintln!("  Reduction factor: {reduction:.2}×");

    // Sanity: the baseline must actually be oscillation-heavy on this
    // mesh — if the test mesh were so coarse that the delta-gap RHS
    // already looks smooth, the comparison would be vacuous.
    assert!(
        dev_dg > 0.5,
        "baseline delta-gap RHS dev = {dev_dg:.3} is suspiciously small; \
         the {N_LENGTH} × {N_WIDTH} mesh may not reproduce the TTTTTTT \
         P1 pattern — re-tune the mesh density before trusting this gate"
    );

    // Headline gate per the brief: TEM weighting reduces the Maxwell
    // envelope deviation by ≥ 5×.
    let target_reduction = 5.0_f64;
    assert!(
        reduction >= target_reduction,
        "TEM-smoothed RHS only reduced Maxwell-envelope deviation by \
         {reduction:.2}× (baseline {dev_dg:.3} → TEM {dev_tem:.3}); \
         brief escape hatch requires ≥ {target_reduction:.0}×. Surface \
         this as a finding — the deeper fix is a true wave-port \
         (Phase 1.3.1.1 step 5 longitudinal block) out of this lane."
    );
}
