//! Benchmark one `update_h + update_e` pair on a 50³ vacuum grid.
//!
//! Uses the `WalkingSkeletonSolver` so the boundary handling (PEC clamp)
//! is included in the per-step cost. The grid is cloned per-iteration via
//! `iter_with_setup` so that the measured region only contains the field
//! updates — the clone is paid in the setup closure and excluded from
//! Criterion's timing.

use criterion::{Criterion, criterion_group, criterion_main};
use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver, YeeGrid};

fn fdtd_step(c: &mut Criterion) {
    let grid = YeeGrid::vacuum(50, 50, 50, 1.0e-3);
    c.bench_function("fdtd_step_50cubed_vacuum", |b| {
        b.iter_batched(
            || WalkingSkeletonSolver::new(grid.clone()),
            |mut solver| {
                solver.step();
            },
            criterion::BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, fdtd_step);
criterion_main!(benches);
