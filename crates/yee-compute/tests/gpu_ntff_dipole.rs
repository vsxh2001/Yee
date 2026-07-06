//! Gate `compute-013` (E.5b): far-field from a **GPU-resident** run — the
//! full-field DFT phasor is accumulated on the GPU every step (no per-step
//! readback), read back once, and projected by the reference
//! `yee_fdtd::NtffState` through two synthetic samples (its accumulation is
//! linear: `sample(Ê_re, t: e^{-jωt}=1)` + `sample(Ê_im, t: e^{-jωt}=+j)`
//! reconstruct `Σ F·e^{-jωt}·Δt` exactly).
//!
//! Two assertions against strong references:
//! 1. the **analytic dipole pattern**: broadside/endfire ≥ 20 dB (the
//!    sin θ endfire null), as in the reference NTFF gate;
//! 2. cross-backend: the broadside magnitude matches the CPU host-adapter
//!    path (gate `compute-010`'s method) within 1 % — FP32 accumulation
//!    noise is O(1e-4).
//!
//! Self-skipping without a wgpu adapter.
//!
//! ```bash
//! cargo test -p yee-compute --release --test gpu_ntff_dipole -- --ignored --nocapture
//! ```

#![cfg(feature = "gpu")]

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use yee_compute::{
    Boundary, ComputeError, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, GpuFdtd,
    Materials, SoftSource, Waveform,
};
use yee_fdtd::{NtffParams, NtffState, YeeGrid};

const N: usize = 50;
const DX: f64 = 1.0e-3;
const NPML: usize = 10;
const N_STEPS: usize = 2000;
const F_PROBE: f64 = 15.0e9;
const SRC: (usize, usize, usize) = (25, 25, 25);
const BOX_MARGIN_CELLS: usize = NPML + 5;

fn copy_into_grid(fields: &Fields, grid: &mut YeeGrid) {
    grid.ex.as_slice_mut().unwrap().copy_from_slice(&fields.ex);
    grid.ey.as_slice_mut().unwrap().copy_from_slice(&fields.ey);
    grid.ez.as_slice_mut().unwrap().copy_from_slice(&fields.ez);
    grid.hx.as_slice_mut().unwrap().copy_from_slice(&fields.hx);
    grid.hy.as_slice_mut().unwrap().copy_from_slice(&fields.hy);
    grid.hz.as_slice_mut().unwrap().copy_from_slice(&fields.hz);
}

fn ntff_state(scratch: &YeeGrid) -> NtffState {
    NtffState::new(
        scratch,
        NtffParams {
            f_probe: F_PROBE,
            box_margin_cells: BOX_MARGIN_CELLS,
            theta_rad: FRAC_PI_2,
            phi_rad: 0.0,
        },
    )
}

#[test]
#[ignore = "slow: 2000-step 50^3 GPU run with per-step on-GPU DFT (release ~1-3 min on llvmpipe); compute-013 GPU NTFF gate (E.5b)"]
fn gpu_accumulated_ntff_recovers_dipole_pattern() {
    let mut scratch = YeeGrid::vacuum(N, N, N, DX);
    let dt = scratch.dt;
    let mut spec = FdtdSpec::vacuum(N, N, N, DX);
    spec.dt = dt;
    let drive = Drive {
        soft_sources: vec![SoftSource {
            component: EComponent::Ez,
            cell: SRC,
            waveform: Waveform::Gaussian {
                t0: 12.0 * dt,
                sigma: 4.0 * dt,
            },
        }],
        ports: vec![],
        aperture_ports: vec![],
        probes: vec![],
    };
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML));

    // ---- GPU run with on-GPU DFT accumulation ----
    let mut gpu = match GpuFdtd::with_ntff_dft(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary.clone(),
        drive.clone(),
        N_STEPS,
        F_PROBE,
    ) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!("SKIPPED compute-013: no wgpu adapter");
            return;
        }
        Err(other) => panic!("compute-013: GPU construction failed: {other}"),
    };
    eprintln!("compute-013: adapter '{}'", gpu.adapter_name());
    gpu.step_n(N_STEPS).expect("GPU stepping failed");
    let (dft_re, dft_im) = gpu.read_dft_fields().expect("DFT readback failed");

    // Feed the phasor pair to the reference NTFF through two synthetic
    // samples: e^{-jω·0} = 1 picks up Ê_re; e^{-jω·t₁} = +j at ωt₁ = 3π/2
    // picks up j·Ê_im. Together: (Ê_re + j·Ê_im)·Δt — the reference
    // accumulation.
    let omega = TAU * F_PROBE;
    let mut ntff_gpu = ntff_state(&scratch);
    copy_into_grid(&dft_re, &mut scratch);
    ntff_gpu.sample(&scratch, 0.0);
    copy_into_grid(&dft_im, &mut scratch);
    ntff_gpu.sample(&scratch, 3.0 * PI / (2.0 * omega));

    let gpu_broad = ntff_gpu.far_field_at(FRAC_PI_2, 0.0).norm();
    let gpu_end = ntff_gpu.far_field_at(0.0, 0.0).norm();
    assert!(gpu_broad > 0.0, "GPU broadside is zero — no radiation");
    let gpu_db = 20.0 * (gpu_broad / gpu_end.max(f64::MIN_POSITIVE)).log10();
    eprintln!(
        "compute-013 [gpu]: |E(broadside)| = {gpu_broad:.3e}, |E(endfire)| = {gpu_end:.3e}, \
         ratio = {gpu_db:.2} dB"
    );
    assert!(
        gpu_db >= 20.0,
        "compute-013: GPU broadside/endfire {gpu_db:.2} dB < 20 dB (analytic sin θ null)"
    );

    // ---- CPU cross-check (the compute-010 host-adapter method) ----
    let mut cpu = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive,
    );
    let mut ntff_cpu = ntff_state(&scratch);
    for _ in 0..N_STEPS {
        cpu.step_n(1);
        copy_into_grid(cpu.fields(), &mut scratch);
        ntff_cpu.sample(&scratch, cpu.current_time());
    }
    let cpu_broad = ntff_cpu.far_field_at(FRAC_PI_2, 0.0).norm();
    let rel = (gpu_broad - cpu_broad).abs() / cpu_broad;
    eprintln!(
        "compute-013: broadside GPU {gpu_broad:.6e} vs CPU {cpu_broad:.6e} → rel {:.3e}",
        rel
    );
    assert!(
        rel < 0.01,
        "compute-013: GPU-accumulated broadside deviates {rel:.3e} (≥ 1 %) from the CPU path"
    );
}
