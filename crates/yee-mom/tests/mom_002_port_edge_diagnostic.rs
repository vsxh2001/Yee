//! mom-002 port-excitation + edge-singularity diagnostic — Track TTTTTTT.
//!
//! ## Why this file exists
//!
//! Track QQQQQQQ (the immediate predecessor; see
//! `tests/mom_002_beta_eigenmode_probe.rs`) extracted the strip
//! propagation constant `β` directly from the port-driven current on
//! the IIIIIII reframed mom-002 mesh (`L = 82 mm` centered uniform on
//! FR-4 at 1 GHz) and found
//!
//! ```text
//!   β / k_0          = 1.840
//!   ε_eff_solver     = 3.385
//!   ε_eff_HJ         = 3.32
//!   relative error   = +1.83 %  (within the ±5 % verdict band)
//! ```
//!
//! That **exonerated the kernel**: the strip-eigenmode physics matches
//! Hammerstad-Jensen to under 2 %. The remaining `|Im(Z)| = 674 Ω`
//! capacitive residual at 1 GHz must therefore live in one of two
//! discretisation effects upstream of the kernel itself:
//!
//! * **(P1) Port-excitation modeling.** A centered delta-gap on a
//!   half-wave open-ended strip excites a standing wave with the same
//!   `β` the kernel got right, but the **port admittance**
//!   `Y_port = I_port / V_port` will be off if the delta-gap is being
//!   driven over a single longitudinal edge instead of the strip's
//!   transverse extent. The canonical microstrip-mode launch puts a 1 V
//!   step across the full transverse strip width at `x = L/2`, with a
//!   current distribution that follows the dominant quasi-TEM
//!   eigenmode in `y` (peaked at strip centre with `1/√(1 − (2y/w)²)`
//!   edge-singularity tails). Deviation from this `y`-profile flags an
//!   ill-formed port.
//!
//! * **(P2) Edge-singularity under-resolution.** Strip-edge currents
//!   carry a `1/√d` singularity in the transverse direction `d → 0`
//!   from the strip edge. The IIIIIII production mesh uses
//!   `dy = 2.94 mm / 16 ≈ 184 μm`. If that under-resolves the
//!   singularity, the effective input impedance picks up a spurious
//!   reactance (often capacitive — the under-resolved edge-charge
//!   capacitance lowers the running `Z_0`). Refining `n_width` from 8
//!   to 32+ should converge `Im(Z)` monotonically if this is the
//!   dominant residual.
//!
//! * **(P3) Port placement sensitivity.** If `L = 82 mm` were exactly
//!   `λ_eff / 2` at 1 GHz, the centered port would drive the
//!   half-wave-resonator's voltage antinode and a quarter-wave-shifted
//!   port would see a substantially different `Z_in` (resistive
//!   `~Z_0`). If `Im(Z)` stays large-and-capacitive regardless of port
//!   placement along the strip, the geometry is just off-resonance and
//!   port placement is not the dominant residual.
//!
//! ## Probe construction
//!
//! Three probes against the same `82 × 16` reframed mom-002 mesh
//! (with `n_width` and `port_left` sweeps for P2 / P3):
//!
//! * **P1** — fill the production impedance matrix `Z` via
//!   `__internal::impedance_matrix_for_test`, drive with the
//!   centered delta-gap RHS via `__internal::delta_gap_rhs_for_test`,
//!   solve `i = Z^-1 b`, then group basis functions by their
//!   shared-edge midpoint and report `|i_k|` and `arg(i_k)` for each
//!   port basis function. The transverse profile `Σ_{k ∈ port_at_y}
//!   length_k · |i_k|` is compared against a Maxwell `1/√(1 − (2y/w)²)`
//!   reference envelope. Deviation flags P1.
//!
//! * **P2** — sweep `n_width ∈ {8, 16, 32}` (escape hatch: dropped 48
//!   per the brief's 10-min-per-solve budget) at the same `L = 82 mm`
//!   / centered port / uniform y-spacing. Each point calls
//!   `__internal::z_in_with_greens` once and reports `Re(Z)`, `Im(Z)`,
//!   `|Z_in|`. Monotonic convergence of `Im(Z)` toward a finite limit
//!   flags P2.
//!
//! * **P3** — at the production `82 × 16` mesh, vary `port_left ∈
//!   {30, 40, 50}` (40-41 is the centered port; 30-31 puts the port
//!   ~10 cells off-centre, a quarter-wave shift on the 82 mm strip at
//!   `ε_eff ≈ 3.32`). Each point calls `__internal::z_in_with_greens`
//!   once. Large `Z_in` swing across the three placements flags P3
//!   (port placement matters → geometry is genuinely on-resonance);
//!   flat `Z_in` rules P3 out (geometry is off-resonance).
//!
//! ## References
//!
//! * ADR-0036 — `docs/src/decisions/0036-mom-002-validation-strategy.md`
//!   (IIIIIII reframe to `L = 82 mm` half-wave).
//! * ADR-0037 — `docs/src/decisions/0037-mom-002-r1-metric-retracted.md`
//!   (R1 metric retracted; β-from-Z recommended).
//! * Sibling diagnostic `tests/mom_002_beta_eigenmode_probe.rs` (QQQQQQQ —
//!   exonerated the kernel via the β-from-Z extraction).
//! * Sibling diagnostic `tests/mom_002_13x_residual_diagnostic.rs` (MMMMMMM —
//!   inlined-mesh and `z_in_with_greens` pattern reused here).
//! * R. F. Harrington, *Time-Harmonic Electromagnetic Fields*, McGraw-Hill,
//!   1961, §5.5 (Maxwell `1/√(1−(2y/w)²)` edge-singularity envelope on a
//!   thin strip).
//! * D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7 (microstrip
//!   `Z_0`, `ε_eff`, half-wave resonator open-circuit input impedance).

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_mesh::TriMesh;
use yee_mom::__internal::{
    MultilayerGreens, build_basis, delta_gap_rhs_for_test, impedance_matrix_for_test,
    z_in_with_greens,
};

use faer::linalg::solvers::{PartialPivLu, Solve};

// FR-4 / 1 GHz canonical microstrip parameters — match
// `yee-validation::MOM_002_*` constants.
const EPS_R: f64 = 4.4;
const H_SUBSTRATE_M: f64 = 1.6e-3;
const STRIP_W_M: f64 = 2.94e-3;
const STRIP_L_M: f64 = 82.0e-3;
const F_HZ: f64 = 1.0e9;

// Production headline mesh dimensions per the IIIIIII reframe.
const N_LENGTH: usize = 82;
const N_WIDTH: usize = 16;

// Sommerfeld kernel parameters — match the production headline gate.
const N_DCIM_IMAGES: usize = 5;
const N_SW_POLES: usize = 1;

/// Build the ADR-0036 / IIIIIII centered-uniform strip mesh with a
/// **configurable** port column. Bit-equivalent to
/// `yee_validation::mom_002_strip_mesh_with_spacing(..,
/// StripSpacing::Uniform)` when `port_left = n_length / 2 - 1` (the
/// production centered port); the `port_left` parameter is exposed so
/// Probe P3 can sweep port placement along the strip. Pattern lifted
/// from `tests/mom_002_beta_eigenmode_probe.rs` and
/// `tests/mom_002_13x_residual_diagnostic.rs`.
fn build_strip_mesh_with_port(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
    port_left: usize,
) -> TriMesh {
    assert!(n_length >= 4, "n_length must be >= 4");
    assert!(n_width >= 1, "n_width must be >= 1");
    assert!(
        port_left + 1 < n_length,
        "port_left + 1 must fit inside n_length"
    );

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

    let port_right = port_left + 1;
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

/// Centered-port convenience wrapper (matches the IIIIIII production
/// mesh exactly).
fn build_strip_mesh_centered(
    length_m: f64,
    width_m: f64,
    n_length: usize,
    n_width: usize,
) -> TriMesh {
    assert!(
        n_length.is_multiple_of(2),
        "n_length must be even for a centered port"
    );
    let port_left = n_length / 2 - 1;
    build_strip_mesh_with_port(length_m, width_m, n_length, n_width, port_left)
}

/// Build the production Sommerfeld kernel at `f = F_HZ` — DCIM `N = 5`
/// + 1 TM₀ surface-wave pole, matching the IIIIIII headline gate.
fn build_greens() -> MultilayerGreens {
    MultilayerGreens::new_microstrip_sommerfeld(
        EPS_R,
        H_SUBSTRATE_M,
        F_HZ,
        N_DCIM_IMAGES,
        N_SW_POLES,
    )
}

/// Solve mom-002 at the supplied geometry through the production
/// Sommerfeld kernel and return `Z_in = V_port / I_port`. Thin wrapper
/// around `__internal::z_in_with_greens`.
fn z_in_at_mesh(mesh: &TriMesh) -> Complex64 {
    let greens = build_greens();
    z_in_with_greens(mesh, 1u32, &greens).expect("z_in_with_greens")
}

// ---------------------------------------------------------------------
// Probe P1 — port-current spatial profile
// ---------------------------------------------------------------------

/// Run Probe P1: fill `Z` and `b` through the production helpers, solve
/// `i = Z^-1 b`, walk the port basis functions, and report per-basis
/// `(y, |i_k|, arg(i_k))` plus the transverse profile `Σ |i_k|` per
/// `y`-bucket. Returns a verdict tag plus the most-deviant `y`-bucket
/// fractional error against the Maxwell envelope.
fn run_probe_p1() -> (&'static str, f64) {
    eprintln!();
    eprintln!("Probe P1 — port-current spatial profile at y_columns");
    eprintln!(
        "  Mesh:      L = {} mm, w = {} mm, {} × {}, centered port, uniform y",
        STRIP_L_M * 1e3,
        STRIP_W_M * 1e3,
        N_LENGTH,
        N_WIDTH,
    );

    let mesh = build_strip_mesh_centered(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH);
    let basis = build_basis(&mesh).expect("RwgBasis build");
    eprintln!("  RWG basis count: {}", basis.n_basis());

    let greens = build_greens();
    eprintln!("  Assembling production Z (Duffy-regularised fill)...");
    let t0 = std::time::Instant::now();
    let z = impedance_matrix_for_test(&basis, &greens);
    eprintln!("    Z fill: {:.1} s", t0.elapsed().as_secs_f64());

    let b = delta_gap_rhs_for_test(&basis, 1u32);
    let t1 = std::time::Instant::now();
    let lu = PartialPivLu::new(z.as_ref());
    let i_vec = lu.solve(b.as_ref());
    eprintln!("    LU+solve: {:.1} s", t1.elapsed().as_secs_f64());

    // Walk all basis functions with port_tag == 1. Per
    // `basis.rs::from_mesh`, the port_tag of a port-straddling edge is
    // `min(tag_plus, tag_minus) = min(1, 2) = 1` — so port_tag == 2 is
    // never assigned by the centered-port mesh. We report this fact and
    // walk port_tag == 1 only.
    let port_basis: Vec<usize> = basis.port_basis_indices(1).collect();
    let port_basis_tag2: Vec<usize> = basis.port_basis_indices(2).collect();
    eprintln!(
        "  Port basis: {} edges at port_tag == 1, {} at port_tag == 2",
        port_basis.len(),
        port_basis_tag2.len(),
    );
    eprintln!("    (Note: convention is port_tag = min(tag_plus, tag_minus); the centered-port");
    eprintln!("     mesh produces port_tag == 1 only — every port-straddling edge is enumerated");
    eprintln!("     in the port_tag == 1 sweep.)");

    // Per-edge: y-midpoint, |i|, arg(i), length.
    eprintln!();
    eprintln!(
        "  {:>4}  {:>10}  {:>10}  {:>14}  {:>10}",
        "k", "y (mm)", "len (mm)", "|i_k| (A·m^-1)", "arg (rad)",
    );
    let mut by_y: Vec<(f64, f64, f64, f64, f64)> = Vec::with_capacity(port_basis.len());
    for &k in &port_basis {
        let edge = basis.edges[k];
        let v0 = mesh.vertices[edge.v0 as usize];
        let v1 = mesh.vertices[edge.v1 as usize];
        let ymid = 0.5 * (v0.y + v1.y);
        let ik = i_vec[(k, 0)];
        let mag = ik.norm();
        let arg = ik.arg();
        let len = edge.length;
        by_y.push((ymid, len, mag, arg, k as f64));
        eprintln!(
            "  {:>4}  {:>10.4}  {:>10.4}  {:>14.4e}  {:>10.4}",
            k,
            ymid * 1e3,
            len * 1e3,
            mag,
            arg,
        );
    }

    // Group port basis functions by their shared-edge midpoint y. The
    // structured mesh produces two kinds of port-straddling edges at
    // each y position: the longitudinal-y edge (parallel to x) whose
    // midpoint is at a vertex y-coordinate, and the diagonal edge whose
    // midpoint is at a cell-centre y-coordinate. We bucket both by
    // discretised y at half-`dy` resolution.
    let dy = STRIP_W_M / (N_WIDTH as f64);
    let bucket_y = |y: f64| -> i64 { ((y + STRIP_W_M / 2.0) / (0.5 * dy)).round() as i64 };
    use std::collections::BTreeMap;
    let mut per_y: BTreeMap<i64, (f64, f64)> = BTreeMap::new();
    for &(ymid, len, mag, _arg, _k) in &by_y {
        let key = bucket_y(ymid);
        let e = per_y.entry(key).or_insert((0.0, 0.0));
        // `e.0` = sum of length-weighted |i|. `e.1` = sum of weights.
        e.0 += len * mag;
        e.1 += len;
    }

    eprintln!();
    eprintln!("  Per-y-column |i| sum across all port RWGs (length-weighted):");
    eprintln!(
        "  {:>10}  {:>14}  {:>14}  {:>10}",
        "y (mm)", "Σ len·|i| (A)", "envelope ref", "ratio",
    );

    // Maxwell `1/√(1 − (2y/w)²)` reference envelope, normalised against
    // the sum at strip centre (y = 0). The denominator is regularised
    // away from the singularity by clipping `|2y/w| ≤ 1 − ε`.
    let envelope_at = |y: f64| -> f64 {
        let u = (2.0 * y / STRIP_W_M).abs().min(1.0 - 1e-3);
        1.0 / (1.0 - u * u).sqrt()
    };

    // Find the y closest to centre for normalisation.
    let mut centre_key = 0_i64;
    let mut centre_dist = f64::INFINITY;
    for &k in per_y.keys() {
        let y = (k as f64) * 0.5 * dy - STRIP_W_M / 2.0;
        if y.abs() < centre_dist {
            centre_dist = y.abs();
            centre_key = k;
        }
    }
    let centre_sum = per_y
        .get(&centre_key)
        .map(|(s, _)| *s)
        .filter(|s| *s > 0.0)
        .unwrap_or(1.0);
    let centre_y = (centre_key as f64) * 0.5 * dy - STRIP_W_M / 2.0;
    let centre_env = envelope_at(centre_y).max(1e-30);

    let mut max_dev: f64 = 0.0;
    let mut edge_peaked = false;
    let mut centre_peaked = false;
    // Capture sums for the y-distribution shape verdict.
    let mut profile: Vec<(f64, f64, f64)> = Vec::with_capacity(per_y.len());
    for (&key, &(sum_li, _)) in &per_y {
        let y = (key as f64) * 0.5 * dy - STRIP_W_M / 2.0;
        let env_y = envelope_at(y);
        let ref_norm = env_y / centre_env;
        let solver_norm = sum_li / centre_sum;
        let ratio = if ref_norm > 0.0 {
            solver_norm / ref_norm
        } else {
            f64::NAN
        };
        eprintln!(
            "  {:>10.4}  {:>14.4e}  {:>14.4e}  {:>10.4}",
            y * 1e3,
            sum_li,
            ref_norm,
            ratio,
        );
        profile.push((y, sum_li, ref_norm));
        if ratio.is_finite() {
            let dev = (ratio - 1.0).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
    }

    // Shape verdict: compare strip-edge sums vs centre sums. Skip if
    // we have fewer than 3 buckets (degenerate).
    if profile.len() >= 3 {
        // Sort by |y| descending so [0..2] are the two strongest edge
        // buckets, [last..] are the central ones.
        let mut by_abs = profile.clone();
        by_abs.sort_by(|a, b| b.0.abs().partial_cmp(&a.0.abs()).unwrap());
        let edge_mean = (by_abs[0].1 + by_abs[1].1) / 2.0;
        let centre_mean = profile.iter().map(|(_, s, _)| *s).sum::<f64>() / (profile.len() as f64);
        // Edge-peaked means edge sum noticeably exceeds the mean
        // (Maxwell singularity wins): edge_mean / centre_mean > 1.2.
        // Centre-peaked means centre sum exceeds edge sum (TM₀ Gaussian
        // tail): edge_mean / centre_mean < 0.8. In between flags
        // approximately uniform.
        let edge_ratio = edge_mean / centre_mean.max(1e-30);
        eprintln!();
        eprintln!(
            "  Profile shape: edge/centre mean ratio = {:.3} (TM₀-mode: <0.8, uniform: 0.8-1.2, edge-peaked: >1.2)",
            edge_ratio,
        );
        if edge_ratio > 1.2 {
            edge_peaked = true;
        } else if edge_ratio < 0.8 {
            centre_peaked = true;
        }
    }

    eprintln!();
    eprintln!(
        "  Max deviation from Maxwell envelope: {:.2} %",
        max_dev * 100.0,
    );

    // Asymmetry detection: split per_y into left (y < 0) and right
    // (y > 0) and compare aggregate sums. > 10 % asymmetry on a
    // centred-port mesh is unphysical.
    let mut left_sum = 0.0;
    let mut right_sum = 0.0;
    for &(y, len_mag, _) in &profile {
        if y < 0.0 {
            left_sum += len_mag;
        } else if y > 0.0 {
            right_sum += len_mag;
        }
    }
    let asym = (left_sum - right_sum).abs() / (left_sum + right_sum).max(1e-30);
    eprintln!(
        "  Left/right asymmetry: {:.2} %  (left = {:.3e}, right = {:.3e})",
        asym * 100.0,
        left_sum,
        right_sum,
    );

    let shape = if asym > 0.10 {
        "asymmetric"
    } else if edge_peaked {
        "edge-peaked"
    } else if centre_peaked {
        "TM0-mode"
    } else {
        "uniform"
    };
    eprintln!();
    eprintln!("  y-distribution shape: {shape}");
    // P1 detected if the profile is materially NOT a TM₀-mode shape.
    // The Maxwell `1/√(1−u²)` envelope is the analytic limit for an
    // infinitely-thin strip with quasi-TEM current; significant
    // deviation (> 30 % at any bucket) flags port-excitation modeling.
    // A truly TM₀-mode shape has edge_peaked = true (Maxwell envelope
    // amplifies edges) and low asymmetry.
    let p1_detected = matches!(shape, "asymmetric" | "uniform") || max_dev > 0.30;
    let verdict = if p1_detected {
        "P1 detected: port-excitation modeling off"
    } else {
        "P1 not detected: looks like proper TM0-mode"
    };
    eprintln!("  Verdict: {verdict}");
    (verdict, max_dev)
}

// ---------------------------------------------------------------------
// Probe P2 — n_width refinement sweep
// ---------------------------------------------------------------------

fn run_probe_p2() -> (&'static str, Vec<(usize, Complex64)>) {
    eprintln!();
    eprintln!("Probe P2 — n_width refinement sweep:");
    eprintln!(
        "  Geometry fixed: L = {} mm, centered port, uniform y-spacing, f = {} GHz",
        STRIP_L_M * 1e3,
        F_HZ * 1e-9,
    );

    // n_width sweep. Per the brief's escape hatch: the production fill
    // at n_width = 32 is already ~2× the IIIIIII (n=16) cost, and
    // n_width = 48 grows the basis count and triangle count by a
    // further ~1.5×, pushing a single solve into the 10+ minute
    // range. We keep {8, 16, 32}; the brief allows dropping the 48
    // point when the wall budget tightens. Document the choice so the
    // verdict reflects the actual sweep.
    let n_widths = [8usize, 16, 32];
    let mut samples: Vec<(usize, Complex64)> = Vec::with_capacity(n_widths.len());
    for &nw in &n_widths {
        let mesh = build_strip_mesh_centered(STRIP_L_M, STRIP_W_M, N_LENGTH, nw);
        let t0 = std::time::Instant::now();
        let z_in = z_in_at_mesh(&mesh);
        let dt = t0.elapsed().as_secs_f64();
        eprintln!(
            "  |Z_in| @ n_width = {:>3}:  Z = {:+10.4} + j{:+10.4} Ω,  |Z| = {:>10.4} Ω   ({:.1} s)",
            nw,
            z_in.re,
            z_in.im,
            z_in.norm(),
            dt,
        );
        samples.push((nw, z_in));
    }

    // Trend classification of Im(Z) across the sweep. With only 3
    // points, "monotonic-converging" means Im(Z) is monotone AND the
    // step shrinks (|d2 - d1| < |d1|, where d1 = Im[1]-Im[0], d2 =
    // Im[2]-Im[1]).
    let ims: Vec<f64> = samples.iter().map(|(_, z)| z.im).collect();
    let d1 = ims[1] - ims[0];
    let d2 = ims[2] - ims[1];
    let monotonic = (d1 > 0.0 && d2 > 0.0) || (d1 < 0.0 && d2 < 0.0);
    let converging = monotonic && d2.abs() < d1.abs() * 0.8;
    let trend = if converging {
        "monotonic-converging"
    } else if monotonic {
        "monotonic but not converging"
    } else if d1.signum() != d2.signum() && (d1.abs() + d2.abs()) > 1.0 {
        "non-monotonic"
    } else {
        "flat"
    };
    eprintln!();
    eprintln!(
        "  Im(Z) deltas:  d1 (n=8→16) = {:+8.3} Ω,  d2 (n=16→32) = {:+8.3} Ω",
        d1, d2,
    );
    eprintln!("  Im(Z) trend:      {trend}");
    // P2 verdict: detected only if refinement shows monotonic
    // convergence (the under-resolved edge is actually being resolved
    // and the answer is moving). Flat / non-monotonic Im(Z) rules
    // edge-singularity out as the dominant residual.
    let verdict = if matches!(
        trend,
        "monotonic-converging" | "monotonic but not converging"
    ) {
        "P2 detected: edge-singularity under-resolution"
    } else {
        "P2 not detected: Im(Z) refinement-insensitive"
    };
    eprintln!("  Verdict: {verdict}");
    (verdict, samples)
}

// ---------------------------------------------------------------------
// Probe P3 — port placement sensitivity
// ---------------------------------------------------------------------

fn run_probe_p3() -> (&'static str, Vec<(usize, Complex64)>) {
    eprintln!();
    eprintln!("Probe P3 — port placement sensitivity:");
    eprintln!(
        "  Mesh fixed: L = {} mm, {} × {}, uniform y. port_left ∈ {{30, 40, 50}}.",
        STRIP_L_M * 1e3,
        N_LENGTH,
        N_WIDTH,
    );
    eprintln!(
        "  port_left = 40 reproduces the centered-port IIIIIII reframe; 30 and 50 displace by"
    );
    eprintln!("  ~10 cells (~10 mm) — roughly a quarter-wave on the 82 mm half-wave strip.");

    let port_lefts = [30usize, 40, 50];
    let mut samples: Vec<(usize, Complex64)> = Vec::with_capacity(port_lefts.len());
    for &pl in &port_lefts {
        let mesh = build_strip_mesh_with_port(STRIP_L_M, STRIP_W_M, N_LENGTH, N_WIDTH, pl);
        let t0 = std::time::Instant::now();
        let z_in = z_in_at_mesh(&mesh);
        let dt = t0.elapsed().as_secs_f64();
        let label = if pl == 40 { "  (centered)" } else { "" };
        eprintln!(
            "  Z_in @ port_left = {:>2}:  {:+10.4} + j{:+10.4} Ω{}   ({:.1} s)",
            pl, z_in.re, z_in.im, label, dt,
        );
        samples.push((pl, z_in));
    }

    // Spread: max - min on Im(Z) across the three port placements.
    // > 30 % spread relative to the mean magnitude flags port placement
    // as influencing the result; otherwise the geometry is just enough
    // off-resonance that the port-displacement effect is masked.
    let ims: Vec<f64> = samples.iter().map(|(_, z)| z.im).collect();
    let im_max = ims.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let im_min = ims.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let im_mag_mean = ims.iter().map(|x| x.abs()).sum::<f64>() / (ims.len() as f64);
    let spread = (im_max - im_min) / im_mag_mean.max(1e-30);
    eprintln!();
    eprintln!(
        "  Im(Z) spread: max - min = {:+8.3} Ω,  mean |Im(Z)| = {:.3} Ω,  spread/mean = {:.2}",
        im_max - im_min,
        im_mag_mean,
        spread,
    );
    let detected = spread > 0.30;
    let verdict = if detected {
        "P3 detected: port-placement-sensitive"
    } else {
        "P3 not detected: Im(Z) port-placement-insensitive"
    };
    eprintln!("  Verdict: {verdict}");
    (verdict, samples)
}

// ---------------------------------------------------------------------
// Test driver
// ---------------------------------------------------------------------

/// Run the three probes and print the consolidated verdict block.
///
/// Wall-time budget: dominated by the production Sommerfeld fill +
/// dense LU. Each `82 × 16` Z + solve is in the same order as the
/// production headline gate. P1 = 1 solve, P2 = 3 solves (n_width = 8,
/// 16, 32), P3 = 3 solves — total 7 solves. Per the brief's escape
/// hatch the n_width = 48 point is dropped (≥ 10 min wall on its own).
#[test]
#[ignore = "diagnostic: probes port-excitation + edge-singularity residual on the IIIIIII reframed mom-002"]
fn mom_002_port_edge_diagnostic() {
    eprintln!("--- Track TTTTTTT: mom-002 port-excitation + edge-singularity diagnostic ---");
    eprintln!();
    eprintln!("Predecessor (QQQQQQQ) verdict: kernel exonerated.");
    eprintln!("  β / k_0 = 1.840, ε_eff_solver = 3.385, ε_eff_HJ = 3.32, ");
    eprintln!("  relative error = +1.83 % (within ±5 % verdict band).");
    eprintln!("  Remaining |Im(Z)| = 674 Ω capacitive residual at 1 GHz must");
    eprintln!("  live in port-excitation modeling (P1), edge-singularity");
    eprintln!("  under-resolution (P2), or port-placement effects (P3).");

    let (p1_verdict, _max_dev) = run_probe_p1();
    let (p2_verdict, p2_samples) = run_probe_p2();
    let (p3_verdict, p3_samples) = run_probe_p3();

    eprintln!();
    eprintln!("--- Consolidated verdict ---");
    eprintln!();
    eprintln!("Probe P1 — port-current spatial profile at y_columns:");
    eprintln!("  Verdict: {p1_verdict}");
    eprintln!();
    eprintln!("Probe P2 — n_width refinement sweep:");
    for (nw, z) in &p2_samples {
        eprintln!(
            "  |Z_in| @ n={:>3}:  {:>10.4} Ω   (Re = {:+10.4}, Im = {:+10.4})",
            nw,
            z.norm(),
            z.re,
            z.im,
        );
    }
    eprintln!("  Verdict: {p2_verdict}");
    eprintln!();
    eprintln!("Probe P3 — port placement:");
    for (pl, z) in &p3_samples {
        let tag = if *pl == 40 { "  (centered)" } else { "" };
        eprintln!(
            "  Z_in @ port_left={:>3}:  {:+10.4} + j{:+10.4} Ω{}",
            pl, z.re, z.im, tag,
        );
    }
    eprintln!("  Verdict: {p3_verdict}");
    eprintln!();

    // Dominant-cause classification. Multiple detections combine into
    // "combination"; no detections drops to "residual not isolated".
    let p1_hit = p1_verdict.starts_with("P1 detected");
    let p2_hit = p2_verdict.starts_with("P2 detected");
    let p3_hit = p3_verdict.starts_with("P3 detected");
    let n_hits = [p1_hit, p2_hit, p3_hit].iter().filter(|b| **b).count();
    let dominant = match (n_hits, p1_hit, p2_hit, p3_hit) {
        (0, _, _, _) => "residual not isolated by these probes",
        (1, true, _, _) => "P1",
        (1, _, true, _) => "P2",
        (1, _, _, true) => "P3",
        _ => "combination",
    };
    eprintln!("Dominant cause: {dominant}");

    // Side note for reviewers: the test asserts nothing (it is a
    // forensic probe), so leaving the diagnostics in eprintln! is
    // deliberate. The `#[ignore]` gate keeps this off the default
    // `cargo test` matrix.
}
