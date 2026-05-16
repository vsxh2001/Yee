//! Phase 1.0 dipole validation against NEC-4 finite-radius reference.

#[path = "fixtures/mod.rs"]
mod fixtures;

use num_complex::Complex64;
use yee_core::{FreqRange, Solver};
use yee_mom::{PlanarMoM, SParameters};

const Z0_REF: f64 = 50.0;
// NEC-4 reference at a=5mm, L=1m, half-wave, lateral surface, delta-gap.
const REFERENCE_RE: f64 = 87.0;
const REFERENCE_IM: f64 = 41.0;
const TOL_RE_REL: f64 = 0.05;
const TOL_IM_REL: f64 = 0.10;
const N_AXIAL: usize = 24;
const N_AROUND: usize = 176;

fn z_in_from_s11(s11: Complex64, z0: f64) -> Complex64 {
    Complex64::new(z0, 0.0) * (Complex64::new(1.0, 0.0) + s11) / (Complex64::new(1.0, 0.0) - s11)
}

/// mom-001 gate: half-wave dipole impedance at resonance.
///
/// Geometry: L = 1.0 m, radius = 5 mm, cylinder lateral surface (no end
/// caps), delta-gap at central edge. Resonance frequency f = c/(2L).
/// Reference: NEC-4 finite-radius wire MoM, `Z ≈ 87 + j41 Ω`.
#[test]
fn dipole_z_at_resonance() {
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, N_AXIAL, N_AROUND);
    let f0 = yee_core::units::C0 / 2.0;
    let freq = FreqRange::new(f0, f0 + 1.0, 1).unwrap();
    let solver = PlanarMoM::default();
    let s = solver.run(&mesh, freq).expect("solve must succeed");
    let s11 = s.data[0][0];
    let z_in = z_in_from_s11(s11, Z0_REF);

    let err_re = (z_in.re - REFERENCE_RE).abs() / REFERENCE_RE;
    let err_im = (z_in.im - REFERENCE_IM).abs() / REFERENCE_IM;
    assert!(
        err_re <= TOL_RE_REL,
        "Re(Z_in) = {:.3} vs NEC-4 reference {}; rel err {:.4} > {}",
        z_in.re,
        REFERENCE_RE,
        err_re,
        TOL_RE_REL
    );
    assert!(
        err_im <= TOL_IM_REL,
        "Im(Z_in) = {:.3} vs NEC-4 reference {}; rel err {:.4} > {}",
        z_in.im,
        REFERENCE_IM,
        err_im,
        TOL_IM_REL
    );
}

/// Phase 1 diagnostic: always-on, never asserts. Prints port edge count,
/// port edge lengths, total RWG count, Z_in, |Z_in|/arg, LU residual,
/// per-port-edge currents (to detect orientation flips), and a radius /
/// mesh sweep showing the convergence behaviour. Numbers are read
/// manually from `--nocapture` output to triangulate any future
/// `dipole_z_at_resonance` regression.
#[test]
fn dipole_z_diagnostics() {
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, N_AXIAL, N_AROUND);
    let basis = yee_mom::__internal::build_basis(&mesh).expect("basis");
    let f0 = yee_core::units::C0 / 2.0;

    let port_indices: Vec<usize> = basis.port_basis_indices(1).collect();
    let port_lengths: Vec<f64> = port_indices
        .iter()
        .map(|&k| basis.edges[k].length)
        .collect();
    let total_port_length: f64 = port_lengths.iter().sum();

    eprintln!("Port edge count       = {}", port_indices.len());
    eprintln!(
        "Port edge length min  = {:.6e} m",
        port_lengths.iter().copied().fold(f64::INFINITY, f64::min)
    );
    eprintln!(
        "Port edge length max  = {:.6e} m",
        port_lengths.iter().copied().fold(0.0, f64::max)
    );
    eprintln!(
        "Sum of port lengths   = {:.6e} m (expected 2π·r = {:.6e})",
        total_port_length,
        std::f64::consts::TAU * 0.005
    );
    eprintln!("Total RWG count       = {}", basis.n_basis());

    let (z_in, lu_residual_norm) =
        yee_mom::__internal::z_in_and_residual_at_freq(&mesh, 1, f0, 50.0).expect("solve");
    eprintln!("Z_in   = {:.4} + j{:.4} Ω", z_in.re, z_in.im);
    eprintln!(
        "|Z|    = {:.4} Ω, arg(Z) = {:.4} rad",
        z_in.norm(),
        z_in.arg()
    );
    eprintln!("Reference (NEC-4, a=5mm): 87 + j41 Ω");
    eprintln!("LU residual ||Zi - b|| / ||b|| = {lu_residual_norm:.3e}");

    // Per-port-edge currents — confirms +/- orientation consistency around
    // the symmetric cylinder ring and reveals partial cancellation, if any.
    // On the N_AXIAL × N_AROUND cylinder the port edges are related by
    // exact discrete rotational symmetry, so all `i_k` are expected to be
    // identical to machine precision; the printed first/last few lines
    // verify this.
    let currents = yee_mom::__internal::port_edge_currents(&mesh, 1, f0).expect("currents");
    let mut total = Complex64::new(0.0, 0.0);
    for (idx, (len, i_k)) in currents.iter().enumerate() {
        total += Complex64::new(*len, 0.0) * *i_k;
        if idx < 3 || idx >= currents.len() - 3 {
            eprintln!(
                "  port[{idx}]: length={:.4e} m, i_k = {:.4e} + j{:.4e} A",
                len, i_k.re, i_k.im
            );
        } else if idx == 3 {
            eprintln!("  ... (middle entries elided)");
        }
    }
    eprintln!(
        "Σ length_k · i_k = {:.4e} + j{:.4e} A (matches I_port from residual helper)",
        total.re, total.im
    );

    // Radius / mesh sweep: a true wire-limit MPIE Z_in approaches the
    // a → 0 sinusoidal-current value (~73 + j42 Ω) only as a/L → 0.
    // Finite a/L produces a smooth monotonic shift to the NEC-4
    // finite-radius value 87 + j41 Ω at a = 5 mm. This block lets a
    // reviewer eyeball the convergence trend without rerunning the
    // gate test.
    for &(name, radius, n_ax, n_ar) in &[
        ("a=5mm,   24×24", 0.005_f64, 24_usize, 24_usize),
        ("a=2mm,   24×24", 0.002_f64, 24_usize, 24_usize),
        ("a=1mm,   24×24", 0.001_f64, 24_usize, 24_usize),
        ("a=0.5mm, 24×24", 0.0005_f64, 24_usize, 24_usize),
        ("a=5mm,   48×24", 0.005_f64, 48_usize, 24_usize),
        ("a=5mm,   24×48", 0.005_f64, 24_usize, 48_usize),
    ] {
        let m = fixtures::cylinder::thin_cylinder(1.0, radius, n_ax, n_ar);
        let (z, _resid) =
            yee_mom::__internal::z_in_and_residual_at_freq(&m, 1, f0, 50.0).expect("solve");
        eprintln!(
            "  [{name}] Z_in = {:.3} + j{:.3} Ω  |Z|={:.3}",
            z.re,
            z.im,
            z.norm()
        );
    }
}

/// Mesh size used for the conditioning regression. Independent of the gate
/// mesh because `condition_number_at_freq` runs an O(N³) SVD that scales
/// catastrophically on the 24×176 gate mesh (≈ 13 632 RWGs, SVD takes
/// hours). On the 24×24 mesh the SVD completes in ~5 s and gives a
/// structurally-equivalent conditioning signal (Phase 1.0 Task 11
/// diagnostics showed cond ≈ 4e7 here; finer meshes inflate this further,
/// but the 1e8 budget exists to catch *regressions*, not to track the
/// genuine fine-mesh growth).
const N_AXIAL_COND: usize = 24;
const N_AROUND_COND: usize = 24;

#[test]
fn condition_number_within_bound() {
    use yee_mom::__internal::condition_number_at_freq;
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, N_AXIAL_COND, N_AROUND_COND);
    let f0 = yee_core::units::C0 / 2.0;
    let cond = condition_number_at_freq(&mesh, 1, f0).expect("cond must succeed");
    // NOTE: cond(Z) on sub-wavelength MPIE meshes is structurally
    // ill-conditioned; 1e8 budget acknowledges this and tracks Phase 1.1
    // loop-tree work to tighten. The 24×24 baseline shows ~4e7; the gate
    // mesh (24×176) pushes conditioning higher but stays solvable via
    // partial-pivoting LU because the matrix is symmetric — measuring it
    // here would require a several-hour SVD run.
    assert!(
        cond <= 1.0e8,
        "cond(Z) = {cond:.3e} exceeds 1e8 — mesh quality regression"
    );
}

#[test]
#[ignore]
fn dipole_full_sweep() {
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, N_AXIAL, N_AROUND);
    let freq = FreqRange::new(130.0e6, 170.0e6, 21).unwrap();
    let solver = PlanarMoM::default();
    let s = solver.run(&mesh, freq).expect("solve must succeed");

    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("results");
    std::fs::create_dir_all(&out_dir).unwrap();
    let path = out_dir.join("dipole.s1p");

    s.write_touchstone(&path, Z0_REF).expect("write_touchstone");

    let file = yee_io::touchstone::read(&path).expect("read back");
    let s2 = SParameters::from_touchstone(&file);
    assert_eq!(s.freq_hz.len(), s2.freq_hz.len());
    for (a, b) in s.freq_hz.iter().zip(s2.freq_hz.iter()) {
        assert!((a - b).abs() <= 1.0e-12 * a.abs().max(1.0));
    }
}
