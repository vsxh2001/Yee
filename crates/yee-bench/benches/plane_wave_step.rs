//! TF/SF plane-wave source step benchmark.
//!
//! Measures the cost of one full 1-D auxiliary-grid advance of
//! [`yee_fdtd::PlaneWaveSource`] (i.e. one `step_incident_h` followed by
//! one `step_incident_e`, which together advance the incident field by
//! a single Yee leapfrog pair). The 3D-grid corrections (`correct_h` /
//! `correct_e`) are deliberately excluded so this bench isolates the
//! 1-D auxiliary-grid hot path.
//!
//! The grid is fixed at 50³ vacuum cells with `dx = 5 mm` and a TF box
//! spanning `i ∈ [15, 35]`, full `j, k`. `PlaneWaveSource` derives
//! `Clone`, so each iteration starts from a fresh source state via
//! `iter_batched` — the clone is paid in the setup closure and excluded
//! from Criterion's timing.
//!
//! The brief originally referenced `step_incident(n)`; the real API
//! exposes the leapfrog as two methods (`step_incident_h` then
//! `step_incident_e`), so the bench wraps both into one "advance" call.

use criterion::{Criterion, criterion_group, criterion_main};
use yee_fdtd::{PlaneWaveDirection, PlaneWaveSource, YeeGrid};

fn plane_wave_step(c: &mut Criterion) {
    let nx = 50;
    let ny = 50;
    let nz = 50;
    let dx = 5.0e-3;
    let grid = YeeGrid::vacuum(nx, ny, nz, dx);
    let dt = 0.99 * grid.courant_limit();

    // TF slab over the full transverse extent, propagating along +x.
    let src = PlaneWaveSource::new(
        15,
        35,
        0,
        ny - 1,
        0,
        nz - 1,
        PlaneWaveDirection::PlusX,
        1.0e9,
        80, // ramp_steps
        dx,
        dt,
        4, // pad
    );

    c.bench_function("plane_wave_step_one_inc_advance", |b| {
        b.iter_batched(
            || src.clone(),
            |mut s| {
                s.step_incident_h();
                s.step_incident_e();
            },
            criterion::BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, plane_wave_step);
criterion_main!(benches);
