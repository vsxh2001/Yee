//! Phase 1.1.0 multilayer Green's function integration tests.
//!
//! These exercise the trait-dispatch plumbing and the limit behaviour of
//! the one-image DCIM placeholder; they are NOT a physics validation gate
//! for microstrip or patch antennas (Phase 1.1.1+ will add those once a
//! real Sommerfeld / DCIM extraction lands).

#[path = "fixtures/mod.rs"]
mod fixtures;

use num_complex::Complex64;
use yee_mom::__internal::{MultilayerGreens, z_in_free_space, z_in_with_greens};

// 24×24 cylinder — same axial / circumferential counts the dipole's
// `condition_number_within_bound` test uses, which keeps the assemble +
// LU solve sub-30 s on the validation host.
const N_AXIAL: usize = 24;
const N_AROUND: usize = 24;
const DIPOLE_LEN_M: f64 = 1.0;
const DIPOLE_RADIUS_M: f64 = 0.005;

fn relative_complex_diff(a: Complex64, b: Complex64) -> f64 {
    (a - b).norm() / b.norm().max(1e-300)
}

/// Trait-dispatch sanity: `MultilayerGreens` with `eps_r = 1` cancels
/// every image contribution (the dielectric-contrast factor
/// `(ε_r − 1) / (ε_r + 1)` is zero), so the resulting Z_in must match the
/// free-space `FreeSpaceGreen` solve on the same mesh to machine
/// precision. This is the Phase 1.1 generic-impedance-matrix smoke test.
#[test]
fn multilayer_eps_r_one_matches_free_space() {
    let mesh = fixtures::cylinder::thin_cylinder(DIPOLE_LEN_M, DIPOLE_RADIUS_M, N_AXIAL, N_AROUND);
    let f0 = yee_core::units::C0 / 2.0; // dipole resonance, λ = 2 m

    // eps_r = 1.0, h = 100 mm — placeholder microstrip stack. With Γ = 0
    // the image weight is zero so the image distance is irrelevant.
    let mg = MultilayerGreens::new_microstrip(f0, 1.0, 100.0e-3);
    let z_in_ml = z_in_with_greens(&mesh, 1, &mg).expect("multilayer solve");
    let z_in_fs = z_in_free_space(&mesh, 1, f0).expect("free-space solve");

    let rel = relative_complex_diff(z_in_ml, z_in_fs);
    // The 1e-10 budget is well above LU residual noise and well below the
    // 10% sanity bound the spec allows for the walking-skeleton DCIM —
    // it catches genuine trait-dispatch bugs without flagging benign
    // numerical jitter.
    assert!(
        rel <= 1.0e-10,
        "MultilayerGreens(eps_r=1) Z_in = {:.4} + j{:.4} \
         differs from FreeSpace Z_in = {:.4} + j{:.4} by rel = {:.3e}",
        z_in_ml.re,
        z_in_ml.im,
        z_in_fs.re,
        z_in_fs.im,
        rel
    );
}

/// `h → ∞` limit: the image pushes to infinite distance and decays as
/// `1/R_image`, so its contribution to G vanishes regardless of the
/// dielectric contrast factor `Γ`. `MultilayerGreens` evaluated with a
/// large h on a non-unity ε_r must therefore still match the free-space
/// Z_in. This isolates the image-geometry path from the Γ-scaling path
/// of `multilayer_eps_r_one_matches_free_space`.
#[test]
fn multilayer_large_h_matches_free_space() {
    let mesh = fixtures::cylinder::thin_cylinder(DIPOLE_LEN_M, DIPOLE_RADIUS_M, N_AXIAL, N_AROUND);
    let f0 = yee_core::units::C0 / 2.0;

    // ε_r = 4.4 (typical FR-4), h = 1e10 m → image is ≈ 2e10 m below the
    // dipole; image distance ≈ 2e10 m → 1/R ≈ 5e-11. The image is
    // numerically zero, so Z_in must equal the free-space result well
    // within 10% — we tighten to 1e-3 to actually exercise the limit.
    let mg = MultilayerGreens::new_microstrip(f0, 4.4, 1.0e10);
    let z_in_ml = z_in_with_greens(&mesh, 1, &mg).expect("multilayer solve");
    let z_in_fs = z_in_free_space(&mesh, 1, f0).expect("free-space solve");

    let rel = relative_complex_diff(z_in_ml, z_in_fs);
    assert!(
        rel <= 1.0e-3,
        "MultilayerGreens(h → ∞) Z_in = {:.4} + j{:.4} \
         differs from FreeSpace Z_in = {:.4} + j{:.4} by rel = {:.3e}",
        z_in_ml.re,
        z_in_ml.im,
        z_in_fs.re,
        z_in_fs.im,
        rel
    );
}

/// Spec-mandated 10% bound: build the substrate stack at the exact
/// parameters called out in the Phase 1.1 plan
/// (`MultilayerGreens::new_microstrip(150 MHz, 1.0, 100 mm)`), run the
/// same dipole solver as the free-space gate, and verify that the
/// multilayer Z_in lands within 10% of the free-space Z_in. With ε_r = 1
/// the image weight is identically zero, so the bound is trivially
/// satisfied to machine precision in this Phase 1.1.0 implementation —
/// the explicit 10% number is preserved here as a tripwire for any
/// future DCIM extraction that re-introduces approximation noise.
#[test]
fn multilayer_spec_bound_at_150_mhz() {
    let mesh = fixtures::cylinder::thin_cylinder(DIPOLE_LEN_M, DIPOLE_RADIUS_M, N_AXIAL, N_AROUND);
    let freq_hz = 150.0e6;

    let mg = MultilayerGreens::new_microstrip(freq_hz, 1.0, 100.0e-3);
    let z_in_ml = z_in_with_greens(&mesh, 1, &mg).expect("multilayer solve");
    let z_in_fs = z_in_free_space(&mesh, 1, freq_hz).expect("free-space solve");

    let rel = relative_complex_diff(z_in_ml, z_in_fs);
    assert!(
        rel <= 0.10,
        "Phase 1.1.0 walking-skeleton DCIM exceeded 10% bound: rel = {:.3e} \
         (multilayer Z = {:.3} + j{:.3}, free-space Z = {:.3} + j{:.3})",
        rel,
        z_in_ml.re,
        z_in_ml.im,
        z_in_fs.re,
        z_in_fs.im
    );
}
