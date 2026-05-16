//! Phase 1.0 dipole validation against Balanis Ch. 8 §8.2 reference.

#[path = "fixtures/mod.rs"]
mod fixtures;

use num_complex::Complex64;
use yee_core::{FreqRange, Solver};
use yee_mom::{PlanarMoM, SParameters};

const Z0_REF: f64 = 50.0;
const REFERENCE_RE: f64 = 73.0;
const REFERENCE_IM: f64 = 42.0;
const TOLERANCE_REL: f64 = 0.05;

fn reference_z_in() -> Complex64 {
    Complex64::new(REFERENCE_RE, REFERENCE_IM)
}

fn rel_diff(a: Complex64, b: Complex64) -> f64 {
    (a - b).norm() / b.norm()
}

fn z_in_from_s11(s11: Complex64, z0: f64) -> Complex64 {
    Complex64::new(z0, 0.0) * (Complex64::new(1.0, 0.0) + s11)
        / (Complex64::new(1.0, 0.0) - s11)
}

/// mom-001 gate: half-wave dipole impedance at resonance.
///
/// Geometry: L = 1.0 m, radius = 5 mm, cylinder lateral surface (no end
/// caps), delta-gap at central edge. Resonance frequency f = c/(2L) = 150 MHz.
/// Reference: Balanis, *Antenna Theory* (4th ed.) Ch. 8 §8.2, Z ≈ 73 + j42 Ω.
#[test]
fn dipole_z_at_resonance() {
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, 24, 24);
    let f0 = yee_core::units::C0 / 2.0; // exactly λ = 2 m
    // Single-point FreqRange requires n_points = 1 with start == stop allowed only
    // if start < stop. FreqRange::new rejects start >= stop. Use a tiny ε.
    let freq = FreqRange::new(f0, f0 + 1.0, 1).unwrap();
    let solver = PlanarMoM::default();
    let s = solver.run(&mesh, freq).expect("solve must succeed");
    let s11 = s.data[0][0];
    let z_in = z_in_from_s11(s11, Z0_REF);
    let err = rel_diff(z_in, reference_z_in());
    assert!(
        err <= TOLERANCE_REL,
        "Z_in = {z_in:.3} vs reference 73+j42 Ω; rel err {err:.4} > {TOLERANCE_REL}"
    );
}

#[test]
fn condition_number_within_bound() {
    use yee_mom::__internal::condition_number_at_freq;
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, 24, 24);
    let f0 = yee_core::units::C0 / 2.0;
    let cond = condition_number_at_freq(&mesh, 1, f0).expect("cond must succeed");
    assert!(
        cond <= 1.0e6,
        "cond(Z) = {cond:.3e} exceeds 1e6 — mesh quality regression"
    );
}

#[test]
#[ignore]
fn dipole_full_sweep() {
    let mesh = fixtures::cylinder::thin_cylinder(1.0, 0.005, 24, 24);
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
