//! `LumpedRlcPort::correct_e` benchmark.
//!
//! Measures the per-call cost of the lumped-resistor semi-implicit
//! E-field correction at a single Yee cell. The grid is 30³ vacuum
//! cells with `dx = 5 mm`; the port is a pure 50 Ω resistor at the
//! grid centre with no driving waveform.
//!
//! [`yee_fdtd::LumpedRlcPort`] is `Clone`, so each iteration starts
//! from a fresh port state via `iter_batched` — this keeps the
//! discrete `α`-saturation state and `e_z_prev` history clean between
//! samples. The grid clone in the setup closure is excluded from
//! Criterion's timing.

use criterion::{Criterion, criterion_group, criterion_main};
use yee_fdtd::{LumpedRlcPort, SourceWaveform, YeeGrid};

fn lumped_resistor_correct(c: &mut Criterion) {
    let grid = YeeGrid::vacuum(30, 30, 30, 5.0e-3);
    let dt = 0.99 * grid.courant_limit();
    let port = LumpedRlcPort::pure_resistor((15, 15, 15), 50.0, SourceWaveform::None);

    c.bench_function("lumped_resistor_correct_e_at_one_cell", |b| {
        b.iter_batched(
            || (grid.clone(), port.clone()),
            |(mut g, mut p)| {
                p.correct_e(&mut g, 0, dt);
            },
            criterion::BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, lumped_resistor_correct);
criterion_main!(benches);
