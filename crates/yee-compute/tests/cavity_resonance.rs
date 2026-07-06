//! Gate `compute-006`: rectangular PEC-cavity TE₁₀₁ resonance on the engine
//! vs the **analytic** Pozar §6.3 frequency — a strong closed-form reference:
//!
//! ```text
//! f₁₀₁ = (c/2)·√((1/a)² + (1/d)²)
//! ```
//!
//! Same method as `yee-fdtd`'s fdtd-201 gate (soft Gaussian E_y drive,
//! 30 000 steps, 400-bin single-DFT scan, ±2.5 % window that grid dispersion
//! on a ~28 cells/λ mesh comfortably meets). Runs on the CPU backend and —
//! when an adapter exists — on the GPU backend, additionally asserting the
//! two extract the same peak to < 0.5 %.
//!
//! `#[ignore]`'d like the reference (~5–15 s release per backend); runs in
//! the compute release gate and on the GPU nightly:
//!
//! ```bash
//! cargo test -p yee-compute --release --test cavity_resonance -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_compute::{
    Boundary, ComputeError, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, GpuFdtd, Materials,
    Probe, SoftSource, Waveform,
};

const NX: usize = 20;
const NY: usize = 10;
const NZ: usize = 20;
const DX: f64 = 0.010;
const N_STEPS: usize = 30_000;
const N_FREQ_BINS: usize = 400;

const SRC: (usize, usize, usize) = (NX / 4, NY / 2, NZ / 4);
const PRB: (usize, usize, usize) = (NX * 3 / 4, NY / 2, NZ * 3 / 4);

const C0: f64 = 299_792_458.0;

fn analytic_f101() -> f64 {
    let a = NX as f64 * DX;
    let d = NZ as f64 * DX;
    0.5 * C0 * ((1.0 / (a * a)) + (1.0 / (d * d))).sqrt()
}

/// Single-bin DFT scan over [0.65, 1.50]·f_ref; returns the peak frequency.
fn peak_frequency(series: &[f64], dt: f64, f_ref: f64) -> f64 {
    let f_lo = 0.65 * f_ref;
    let f_hi = 1.50 * f_ref;
    let df = (f_hi - f_lo) / (N_FREQ_BINS - 1) as f64;
    let mut peak_power = 0.0_f64;
    let mut peak_freq = f_lo;
    for bin in 0..N_FREQ_BINS {
        let f = f_lo + bin as f64 * df;
        let omega = 2.0 * PI * f;
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for (n, &x) in series.iter().enumerate() {
            let phase = omega * n as f64 * dt;
            re += x * phase.cos();
            im -= x * phase.sin();
        }
        let power = re * re + im * im;
        if power > peak_power {
            peak_power = power;
            peak_freq = f;
        }
    }
    peak_freq
}

fn cavity_drive(dt: f64) -> Drive {
    Drive {
        soft_sources: vec![SoftSource {
            component: EComponent::Ey,
            cell: SRC,
            waveform: Waveform::Gaussian {
                t0: 12.0 * dt,
                sigma: 4.0 * dt,
            },
        }],
        ports: vec![],
        probes: vec![Probe {
            component: EComponent::Ey,
            cell: PRB,
        }],
    }
}

#[test]
#[ignore = "slow: ~5-15 s release per backend; compute-006 TE101 cavity resonance gate (E.2)"]
fn te101_resonance_matches_analytic_on_both_backends() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let dt = spec.dt;
    let f_ref = analytic_f101();

    // ---- CPU backend ----
    let mut cpu = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        Boundary::PecBox,
        cavity_drive(dt),
    );
    cpu.step_n(N_STEPS);
    let f_cpu = peak_frequency(&cpu.probe_series()[0], dt, f_ref);
    let cpu_err = (f_cpu - f_ref) / f_ref;
    eprintln!(
        "compute-006 [cpu]: analytic f101 = {:.6} GHz, extracted = {:.6} GHz, err = {:.4} %",
        f_ref * 1e-9,
        f_cpu * 1e-9,
        cpu_err * 100.0
    );
    assert!(
        cpu_err.abs() < 0.025,
        "compute-006 [cpu]: extracted {f_cpu:.6e} vs analytic {f_ref:.6e} — err {:.4} % > 2.5 %",
        cpu_err * 100.0
    );

    // ---- GPU backend (skips without an adapter) ----
    let mut gpu = match GpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        Boundary::PecBox,
        cavity_drive(dt),
        N_STEPS,
    ) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!("compute-006 [gpu]: SKIPPED — no wgpu adapter");
            return;
        }
        Err(other) => panic!("compute-006 [gpu]: construction failed: {other}"),
    };
    eprintln!("compute-006 [gpu]: adapter '{}'", gpu.adapter_name());
    gpu.step_n(N_STEPS).expect("GPU stepping failed");
    let gpu_series = gpu.read_probes().expect("GPU probe readback failed");
    let f_gpu = peak_frequency(&gpu_series[0], dt, f_ref);
    let gpu_err = (f_gpu - f_ref) / f_ref;
    let cross_err = (f_gpu - f_cpu) / f_cpu;
    eprintln!(
        "compute-006 [gpu]: extracted = {:.6} GHz, err vs analytic = {:.4} %, vs CPU = {:.4} %",
        f_gpu * 1e-9,
        gpu_err * 100.0,
        cross_err * 100.0
    );
    assert!(
        gpu_err.abs() < 0.025,
        "compute-006 [gpu]: extracted {f_gpu:.6e} vs analytic {f_ref:.6e} — err {:.4} % > 2.5 %",
        gpu_err * 100.0
    );
    assert!(
        cross_err.abs() < 0.005,
        "compute-006: CPU/GPU peak disagreement {:.4} % > 0.5 %",
        cross_err * 100.0
    );
}
