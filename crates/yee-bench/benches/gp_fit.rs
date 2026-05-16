//! Gaussian-process fit + fit_ml benchmarks.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use nalgebra::{DMatrix, DVector};
use yee_surrogate::{GaussianProcess, MlFitConfig};

fn build_sin_dataset(n: usize) -> (DMatrix<f64>, DVector<f64>) {
    let mut x = DMatrix::<f64>::zeros(n, 1);
    let mut y = DVector::<f64>::zeros(n);
    for i in 0..n {
        let xi = (i as f64) * 2.0 * std::f64::consts::PI / (n as f64 - 1.0);
        x[(i, 0)] = xi;
        y[i] = xi.sin();
    }
    (x, y)
}

fn gp_fit_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("gp_fit");
    for &n in &[10usize, 25, 50, 100] {
        let (x, y) = build_sin_dataset(n);
        group.bench_with_input(BenchmarkId::new("sin_1d", n), &n, |b, _| {
            b.iter(|| {
                let _ = GaussianProcess::fit(x.clone(), y.clone(), 0.5, 1.0, 1e-4);
            })
        });
    }
    group.finish();
}

fn gp_fit_ml_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("gp_fit_ml");
    for &n in &[10usize, 25, 50] {
        let (x, y) = build_sin_dataset(n);
        group.bench_with_input(BenchmarkId::new("sin_1d", n), &n, |b, _| {
            b.iter(|| {
                let cfg = MlFitConfig::default();
                let _ = GaussianProcess::fit_ml(x.clone(), y.clone(), cfg);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, gp_fit_bench, gp_fit_ml_bench);
criterion_main!(benches);
