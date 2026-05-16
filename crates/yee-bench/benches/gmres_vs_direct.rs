//! Compare iterative `gmres_jacobi` against a direct `PartialPivLu` solve on
//! a diagonally-dominant 128×128 Hermitian-positive-definite-like matrix.
//!
//! The test matrix is built as `A = 10·I + 0.01·sin((i+j))·(off-diagonal)`,
//! matching the diagonal-dominant fixture used in the
//! `yee_mom::iterative::gmres_diagonal_dominant_converges` unit test. At
//! n = 128 both paths should complete in single-digit milliseconds, and the
//! two timings are expected to land within an order of magnitude of each
//! other on a typical dev laptop.

use criterion::{Criterion, criterion_group, criterion_main};
use faer::{
    Mat,
    linalg::solvers::{PartialPivLu, Solve},
};
use num_complex::Complex64;
use yee_mom::{GmresParams, gmres_jacobi};

fn build_hpd(n: usize) -> Mat<Complex64> {
    let mut a = Mat::<Complex64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let v = if i == j {
                Complex64::new(10.0, 0.0)
            } else {
                Complex64::new(0.01 * ((i + j) as f64).sin(), 0.0)
            };
            a[(i, j)] = v;
        }
    }
    a
}

fn gmres_vs_direct(c: &mut Criterion) {
    let n = 128;
    let a = build_hpd(n);
    let mut b = Mat::<Complex64>::zeros(n, 1);
    for i in 0..n {
        b[(i, 0)] = Complex64::new(1.0, 0.0);
    }

    c.bench_function("direct_lu_128", |bb| {
        bb.iter(|| {
            let lu = PartialPivLu::new(a.as_ref());
            let _x = lu.solve(b.as_ref());
        })
    });

    c.bench_function("gmres_jacobi_128", |bb| {
        bb.iter(|| {
            let _ = gmres_jacobi(a.as_ref(), b.as_ref(), GmresParams::default());
        })
    });
}

criterion_group!(benches, gmres_vs_direct);
criterion_main!(benches);
