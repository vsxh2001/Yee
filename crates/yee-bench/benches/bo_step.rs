//! Bayesian optimization end-to-end (cheap synthetic objective).

use criterion::{Criterion, criterion_group, criterion_main};
use nalgebra::DVector;
use yee_surrogate::{BoConfig, minimize as bo_minimize};

fn bo_full_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("bo");
    group.sample_size(10); // BO is expensive; 10 samples is enough for a coarse number
    group.bench_function("deceptive_1d_n_evals_25", |b| {
        b.iter(|| {
            let f = |x: &DVector<f64>| (x[0] - 3.0).powi(2) + (5.0 * x[0]).sin();
            let cfg = BoConfig {
                n_initial: 5,
                n_iters: 20,
                n_candidates: 1024,
                xi: 0.01,
                seed: 0xC0FFEE,
            };
            let _ = bo_minimize(f, vec![(0.0, 6.0)], cfg);
        })
    });
    group.finish();
}

criterion_group!(benches, bo_full_run);
criterion_main!(benches);
