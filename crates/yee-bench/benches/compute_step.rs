//! E.4 engine benchmarks: one leapfrog step at 64³ on the three execution
//! paths — `yee-fdtd` scalar reference, `yee-compute` rayon CPU, and (when
//! an adapter exists) `yee-compute` wgpu GPU.
//!
//! The GPU case measures `step_n(STEPS_PER_ITER)` plus a full field readback
//! per iteration and reports per-step time; the readback amortizes to
//! `1/STEPS_PER_ITER` of an iteration — real workloads read back far less
//! often, so the GPU number here is a *lower* bound on its advantage.

use criterion::{Criterion, criterion_group, criterion_main};
use yee_compute::{CpuFdtd, FdtdSpec, Fields, GpuFdtd};
use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver, YeeGrid};

const N: usize = 64;
const DX: f64 = 1.0e-3;
const GPU_STEPS_PER_ITER: usize = 64;

fn compute_step(c: &mut Criterion) {
    // --- yee-fdtd scalar reference (incl. PEC clamp, as fdtd_step.rs) ---
    let grid = YeeGrid::vacuum(N, N, N, DX);
    c.bench_function("step_64cubed_scalar_reference", |b| {
        b.iter_batched(
            || WalkingSkeletonSolver::new(grid.clone()),
            |mut solver| solver.step(),
            criterion::BatchSize::SmallInput,
        );
    });

    // --- yee-compute rayon CPU (raw kernels, Boundary::None) ---
    let spec = FdtdSpec::vacuum(N, N, N, DX);
    let init = Fields::with_gaussian_ez(&spec, (N / 2, N / 2, N / 2), 3.0);
    c.bench_function("step_64cubed_compute_cpu_rayon", |b| {
        b.iter_batched(
            || CpuFdtd::new(spec, init.clone()),
            |mut engine| engine.step_n(1),
            criterion::BatchSize::SmallInput,
        );
    });

    // --- yee-compute wgpu GPU (skipped when no adapter) ---
    match GpuFdtd::new(spec, init.clone()) {
        Err(_) => eprintln!("compute_step: no wgpu adapter — GPU bench skipped"),
        Ok(_probe) => {
            drop(_probe);
            c.bench_function(
                &format!("step_64cubed_compute_gpu_x{GPU_STEPS_PER_ITER}_incl_readback"),
                |b| {
                    b.iter_batched(
                        || GpuFdtd::new(spec, init.clone()).expect("adapter vanished"),
                        |mut engine| {
                            engine.step_n(GPU_STEPS_PER_ITER).expect("GPU step failed");
                            let _ = engine.read_fields().expect("readback failed");
                        },
                        criterion::BatchSize::SmallInput,
                    );
                },
            );
        }
    }
}

criterion_group!(benches, compute_step);
criterion_main!(benches);
