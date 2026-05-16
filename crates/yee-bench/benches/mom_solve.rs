//! Benchmark the full `PlanarMoM::run` on a small thin-cylinder dipole mesh.
//!
//! Mesh size is `n_axial = 8`, `n_around = 8` (128 triangles), and the sweep
//! is a single frequency near the half-wave dipole resonance. This is small
//! enough that one iteration of `b.iter` completes in well under a second
//! on a typical dev laptop, leaving Criterion's default warm-up / sample
//! configuration entirely within a ~1 min wall budget for the whole bench.

use criterion::{Criterion, criterion_group, criterion_main};
use nalgebra::Vector3;
use yee_core::{FreqRange, Solver};
use yee_mesh::TriMesh;
use yee_mom::PlanarMoM;

fn mom_solve(c: &mut Criterion) {
    let mesh = thin_cylinder_inline(1.0, 0.005, 8, 8);
    let freq = FreqRange::new(1.49e8, 1.50e8, 1).unwrap();
    c.bench_function("mom_solve_dipole_8x8_single_freq", |b| {
        b.iter(|| {
            let _ = PlanarMoM::default().run(&mesh, freq);
        })
    });
}

/// Inlined thin-cylinder mesh generator. Copied from
/// `crates/yee-mom/tests/fixtures/cylinder.rs` because the fixture lives in
/// the `yee-mom` test tree and is not exported on the crate's public surface.
/// Keep in sync with the source-of-truth fixture if the tagging convention
/// changes.
fn thin_cylinder_inline(length_m: f64, radius_m: f64, n_axial: usize, n_around: usize) -> TriMesh {
    assert!(
        n_axial >= 2 && n_axial.is_multiple_of(2),
        "n_axial must be even and >= 2"
    );
    assert!(n_around >= 3, "n_around must be >= 3");

    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity((n_axial + 1) * n_around);
    let dz = length_m / (n_axial as f64);
    let z0 = -length_m / 2.0;
    let dtheta = std::f64::consts::TAU / (n_around as f64);

    for i in 0..=n_axial {
        let z = z0 + (i as f64) * dz;
        for j in 0..n_around {
            let theta = (j as f64) * dtheta;
            vertices.push(Vector3::new(
                radius_m * theta.cos(),
                radius_m * theta.sin(),
                z,
            ));
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_axial * n_around);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_axial * n_around);
    let central_ring = n_axial / 2;

    for i in 0..n_axial {
        for j in 0..n_around {
            let j_next = (j + 1) % n_around;
            let a = (i * n_around + j) as u32;
            let b = (i * n_around + j_next) as u32;
            let c = ((i + 1) * n_around + j_next) as u32;
            let d = ((i + 1) * n_around + j) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == central_ring - 1 {
                1
            } else if i == central_ring {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }

    TriMesh::new(vertices, triangles, tags).expect("cylinder mesh invariants")
}

criterion_group!(benches, mom_solve);
criterion_main!(benches);
