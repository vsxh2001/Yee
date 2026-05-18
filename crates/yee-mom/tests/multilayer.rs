//! Multilayer Green's function integration tests.
//!
//! These exercise the trait-dispatch plumbing and the limit behaviour of
//! both the Phase 1.1.0 one-image DCIM placeholder and the Phase 1.1.1.0
//! N-image DCIM extension. They are NOT a physics validation gate for
//! microstrip or patch antennas (Phase 1.1.1.1+ will tighten those once
//! real Sommerfeld-integral extraction with surface-wave pole subtraction
//! lands).

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

/// Diagnostic-only: PEC-mirror image at `b = -1`, `a = -2h` (the N=1
/// limit of the PEC-backed slab in the thin-substrate limit). Compare
/// against the GPOF-fitted N-image result; if the simple PEC mirror
/// gives a number close to the GPOF result, that is evidence the
/// substrate is "too thin" for the additional images to add
/// information at this geometry.
#[test]
#[ignore = "diagnostic-only: PEC mirror Z at the microstrip"]
fn dcim_pec_mirror_microstrip() {
    use nalgebra::Vector3;
    use num_complex::Complex64;
    use yee_mom::__internal::{MultilayerGreens, z_in_with_greens};

    let length_m = 30.0e-3;
    let width_m = 2.94e-3;
    let n_length = 30usize;
    let n_width = 2usize;
    let h = 1.6e-3;
    let f_hz = 1.0e9;

    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::new();
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;
    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..ny {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == 0 {
                1
            } else if i == 1 {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    let mesh = yee_mesh::TriMesh::new(vertices, triangles, tags).unwrap();

    // Build MultilayerGreens with manually-set N=1 PEC-mirror image.
    let mut mg = MultilayerGreens::new_microstrip(f_hz, 4.4, h);
    // Replace placeholder dielectric-mirror with PEC-mirror image.
    mg.vector_images = vec![(Complex64::new(-1.0, 0.0), Complex64::new(-2.0 * h, 0.0))];
    mg.scalar_images = vec![(Complex64::new(-1.0, 0.0), Complex64::new(-2.0 * h, 0.0))];

    let z_in = z_in_with_greens(&mesh, 1, &mg).unwrap();
    let s11 = (z_in - Complex64::new(50.0, 0.0)) / (z_in + Complex64::new(50.0, 0.0));
    println!(
        "  PEC-mirror b=-1 at a=-2h: Z_in = {:.3} + j{:.3} Ohm, |Z| = {:.3}, |S11| = {:.4}",
        z_in.re,
        z_in.im,
        z_in.norm(),
        s11.norm(),
    );
}

/// Diagnostic-only: sweep `n_images` for the FR-4 microstrip strip
/// from `mom-002` (30 mm × 2.94 mm at 1 GHz) and report the recovered
/// `|Z_in|`. Used to characterise the GPOF fit's stability under the
/// spec escape hatch ("fall back to N=3 if N=5 unstable"). Not part
/// of the gated CI surface; ignored so the regular fast tests stay
/// silent.
#[test]
#[ignore = "diagnostic-only: sweep n_images on the microstrip"]
fn dcim_n_sweep_microstrip() {
    use yee_core::{FreqRange, Solver};
    use yee_mom::{GreensSpec, PlanarMoM};

    let length_m = 30.0e-3;
    let width_m = 2.94e-3;
    let n_length = 30usize;
    let n_width = 2usize;
    let eps_r = 4.4;
    let h = 1.6e-3;
    let f_hz = 1.0e9;

    let nx = n_length + 1;
    let ny = n_width + 1;
    use nalgebra::Vector3;
    let mut vertices: Vec<Vector3<f64>> = Vec::new();
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;
    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..ny {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == 0 {
                1
            } else if i == 1 {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    let mesh = yee_mesh::TriMesh::new(vertices, triangles, tags).unwrap();

    let freq = FreqRange::new(f_hz, f_hz + 1.0, 1).unwrap();
    for n in [1usize, 2, 3, 5, 7, 10] {
        let solver = PlanarMoM::default().with_greens(GreensSpec::microstrip_dcim(eps_r, h, n));
        let s = solver.run(&mesh, freq).unwrap();
        let s11 = s.data[0][0];
        let z = Complex64::new(50.0, 0.0) * (Complex64::new(1.0, 0.0) + s11)
            / (Complex64::new(1.0, 0.0) - s11);
        println!(
            "  N={n}: Z_in = {:.3} + j{:.3} Ohm, |Z| = {:.3} Ohm",
            z.re,
            z.im,
            z.norm()
        );
    }
}

/// ADR-0020 tripwire at the `GreensSpec` enum level: building a
/// `PlanarMoM` with `GreensSpec::MicrostripSommerfeld { ..,
/// n_surface_wave_poles: 0 }` must produce the *same* S-matrix as
/// `GreensSpec::MicrostripDcim { .. }` with matching `eps_r`, `h_m`,
/// `n_images`. The constructor-level tripwire lives in
/// `multilayer::tests::sommerfeld_n_sw_poles_zero_matches_phase_1_1_1_0`
/// — this one closes the loop on the *public* enum dispatch path so a
/// later refactor of the `build` match arm cannot silently regress.
///
/// Uses the short FR-4 strip from `dcim_n_sweep_microstrip` minus the
/// frequency loop. The 1e-12 tolerance is the same element-wise budget
/// the constructor tripwire enforces; anything tighter is below LU
/// residual noise on the ~120-edge mesh.
#[test]
fn greens_spec_sommerfeld_zero_poles_matches_dcim() {
    use nalgebra::Vector3;
    use yee_core::{FreqRange, Solver};
    use yee_mom::{GreensSpec, PlanarMoM};

    let length_m = 30.0e-3;
    let width_m = 2.94e-3;
    let n_length = 10usize;
    let n_width = 2usize;
    let eps_r = 4.4;
    let h = 1.6e-3;
    let n_images = 5usize;
    let f_hz = 1.0e9;

    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::new();
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;
    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..ny {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            // Port edge convention: tag=1 / tag=2 across the first column
            // boundary creates a delta-gap port (basis::port_basis_indices
            // expects DIFFERENT non-zero tags on adjacent triangles).
            let tag = if i == 0 {
                1
            } else if i == 1 {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    let mesh = yee_mesh::TriMesh::new(vertices, triangles, tags).unwrap();
    let freq = FreqRange::new(f_hz, f_hz + 1.0, 1).unwrap();

    let dcim = PlanarMoM::default()
        .with_greens(GreensSpec::microstrip_dcim(eps_r, h, n_images))
        .run(&mesh, freq)
        .expect("dcim solve");
    let somm0 = PlanarMoM::default()
        .with_greens(GreensSpec::microstrip_sommerfeld(eps_r, h, n_images, 0))
        .run(&mesh, freq)
        .expect("sommerfeld n_sw=0 solve");

    let s_dcim = dcim.data[0][0];
    let s_somm = somm0.data[0][0];
    let diff = (s_dcim - s_somm).norm();
    assert!(
        diff <= 1.0e-12,
        "GreensSpec::MicrostripSommerfeld {{ n_surface_wave_poles: 0 }} \
         must reduce to MicrostripDcim bit-for-bit (ADR-0020 tripwire): \
         |ΔS11| = {diff:.3e}, S_dcim = {s_dcim}, S_somm0 = {s_somm}"
    );
}

/// Sanity check that the Sommerfeld pole-subtraction is actually
/// exercised: with `n_surface_wave_poles = 1` on FR-4 / 1.6 mm / 1 GHz
/// the TM₀ pole is found and the analytic Hankel residue is added back
/// to the kernel, so the S-matrix must differ from the
/// `n_surface_wave_poles = 0` (pure DCIM) path by a non-trivial margin.
/// Guards against a regression that silently routes the `n_sw > 0` arm
/// through the unsubtracted code path.
#[test]
fn greens_spec_sommerfeld_one_pole_differs_from_zero() {
    use nalgebra::Vector3;
    use yee_core::{FreqRange, Solver};
    use yee_mom::{GreensSpec, PlanarMoM};

    let length_m = 30.0e-3;
    let width_m = 2.94e-3;
    let n_length = 10usize;
    let n_width = 2usize;
    let eps_r = 4.4;
    let h = 1.6e-3;
    let n_images = 5usize;
    let f_hz = 1.0e9;

    let nx = n_length + 1;
    let ny = n_width + 1;
    let mut vertices: Vec<Vector3<f64>> = Vec::new();
    let dx = length_m / (n_length as f64);
    let dy = width_m / (n_width as f64);
    let y0 = -width_m / 2.0;
    for i in 0..nx {
        let x = (i as f64) * dx;
        for j in 0..ny {
            let y = y0 + (j as f64) * dy;
            vertices.push(Vector3::new(x, y, 0.0));
        }
    }
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();
    for i in 0..n_length {
        for j in 0..n_width {
            let a = (i * ny + j) as u32;
            let b = ((i + 1) * ny + j) as u32;
            let c = ((i + 1) * ny + (j + 1)) as u32;
            let d = (i * ny + (j + 1)) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            // Port edge convention: tag=1 / tag=2 across the first column
            // boundary creates a delta-gap port (basis::port_basis_indices
            // expects DIFFERENT non-zero tags on adjacent triangles).
            let tag = if i == 0 {
                1
            } else if i == 1 {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }
    let mesh = yee_mesh::TriMesh::new(vertices, triangles, tags).unwrap();
    let freq = FreqRange::new(f_hz, f_hz + 1.0, 1).unwrap();

    let somm0 = PlanarMoM::default()
        .with_greens(GreensSpec::microstrip_sommerfeld(eps_r, h, n_images, 0))
        .run(&mesh, freq)
        .expect("sommerfeld n_sw=0 solve");
    let somm1 = PlanarMoM::default()
        .with_greens(GreensSpec::microstrip_sommerfeld(eps_r, h, n_images, 1))
        .run(&mesh, freq)
        .expect("sommerfeld n_sw=1 solve");

    let s0 = somm0.data[0][0];
    let s1 = somm1.data[0][0];
    let diff = (s0 - s1).norm();
    assert!(
        diff > 1.0e-6,
        "n_surface_wave_poles=1 must differ from n_surface_wave_poles=0 \
         on FR-4 / 1 GHz (TM₀ pole subtraction is non-trivial): \
         |ΔS11| = {diff:.3e}, S_n0 = {s0}, S_n1 = {s1}"
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
