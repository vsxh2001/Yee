//! Gate `compute-015` (R.3, ADR-0196): the wgpu FP32 backend runs the
//! **aperture-port** design-flow stack — per-cell ε_r substrate, PEC
//! trace/ground masks, CPML, a driven and a matched multi-cell aperture
//! port — and its probe series match the FP64 CPU backend (itself
//! bit-exact vs `yee_fdtd`, gate compute-014) within FP32 accumulation
//! tolerance, family-relative (the E.3 idiom).
//!
//! **Self-skipping** without a wgpu adapter; real on Mesa llvmpipe and the
//! GPU nightly runner.

#![cfg(feature = "gpu")]

use yee_compute::{
    AperturePort, Boundary, ComputeError, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields,
    GpuFdtd, Materials, Probe, Waveform,
};

const NX: usize = 40;
const NY: usize = 18;
const NZ: usize = 14;
const DX: f64 = 5e-4;
const NPML: usize = 6;
const STEPS: usize = 400;

const REL_L2_TOL: f64 = 1e-3;
const REL_LINF_TOL: f64 = 5e-3;

const J_LO: usize = 7; // trace band [7, 11) → 4 columns
const J_HI: usize = 11;
const K_SUB: usize = 4; // substrate cells k = 0..4, trace at k = 4
const PORT_I_DRIVE: usize = 10;
const PORT_I_LOAD: usize = 30;

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

/// Miniature S.10 board: FR-4-ish substrate filling k < K_SUB, PEC ground
/// sheet (E_x and E_y masked at k = 0) and a PEC trace strip (E_x masked at
/// k = K_SUB over the band).
fn materials() -> Materials {
    let celld = (NX + 1, NY + 1, NZ + 1);
    let mut eps_r = vec![1.0; celld.0 * celld.1 * celld.2];
    for i in 0..celld.0 {
        for j in 0..celld.1 {
            for k in 0..K_SUB {
                eps_r[(i * celld.1 + j) * celld.2 + k] = 4.4;
            }
        }
    }
    let exd = (NX, NY + 1, NZ + 1);
    let mut mask_ex = vec![false; exd.0 * exd.1 * exd.2];
    let eyd = (NX + 1, NY, NZ + 1);
    let mut mask_ey = vec![false; eyd.0 * eyd.1 * eyd.2];
    for i in 0..exd.0 {
        for j in 0..exd.1 {
            mask_ex[(i * exd.1 + j) * exd.2] = true; // ground, k = 0
            if (J_LO..J_HI).contains(&j) {
                mask_ex[(i * exd.1 + j) * exd.2 + K_SUB] = true; // trace
            }
        }
    }
    for i in 0..eyd.0 {
        for j in 0..eyd.1 {
            mask_ey[(i * eyd.1 + j) * eyd.2] = true; // ground, k = 0
        }
    }
    Materials {
        eps_r_cells: Some(eps_r),
        pec_mask_ex: Some(mask_ex),
        pec_mask_ey: Some(mask_ey),
        ..Materials::default()
    }
}

fn aperture(i: usize, resistance: f64, v0: f64, dt: f64) -> AperturePort {
    let cells: Vec<(usize, usize, usize)> = (J_LO..J_HI)
        .flat_map(|j| (0..K_SUB).map(move |k| (i, j, k)))
        .collect();
    let n_columns = J_HI - J_LO;
    let height = K_SUB as f64 * DX;
    let area = n_columns as f64 * DX * height;
    let f0 = 20.0e9;
    let bw = 0.8 * f0;
    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw))
        / dt)
        .ceil() as usize;
    AperturePort {
        cells,
        n_columns,
        area,
        height,
        resistance,
        waveform: Waveform::GaussianPulse {
            v0,
            f0,
            bw,
            t0_steps,
        },
        record: false,
    }
}

#[test]
fn gpu_aperture_ports_match_cpu_probe_series() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let dt = spec.dt;
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML).with_axes([true, true, false]));
    let drive = Drive {
        soft_sources: vec![],
        ports: vec![],
        aperture_ports: vec![
            aperture(PORT_I_DRIVE, 50.0, 1.0, dt),
            aperture(PORT_I_LOAD, 50.0, 0.0, dt),
        ],
        probes: vec![
            Probe {
                component: EComponent::Ez,
                cell: (16, (J_LO + J_HI) / 2, K_SUB - 1),
            },
            Probe {
                component: EComponent::Ez,
                cell: (24, (J_LO + J_HI) / 2, K_SUB - 1),
            },
            Probe {
                component: EComponent::Ez,
                cell: (24, J_HI + 2, K_SUB + 3),
            },
        ],
        h_probes: vec![],
    };

    let mut gpu = match GpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        materials(),
        boundary.clone(),
        drive.clone(),
        STEPS,
    ) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!(
                "SKIPPED compute-015: no wgpu adapter on this machine; \
                 llvmpipe / the GPU nightly runner exercise this gate for real"
            );
            return;
        }
        Err(other) => panic!("compute-015: GPU backend construction failed: {other}"),
    };
    eprintln!("compute-015: running on adapter '{}'", gpu.adapter_name());

    let mut cpu = CpuFdtd::with_drive(spec, Fields::zero(&spec), materials(), boundary, drive);
    cpu.step_n(STEPS);
    gpu.step_n(STEPS).expect("GPU stepping failed");
    let gpu_probes = gpu.read_probes().expect("GPU probe readback failed");
    let cpu_probes = cpu.probe_series();

    assert_eq!(gpu_probes.len(), cpu_probes.len());
    let family: Vec<f64> = cpu_probes.iter().flatten().copied().collect();
    let family_l2 = l2(&family);
    let family_linf = linf(&family);
    assert!(
        family_l2 > 0.0,
        "compute-015: the CPU reference never energized the probes — scenario broken"
    );

    for (q, (reference, candidate)) in cpu_probes.iter().zip(&gpu_probes).enumerate() {
        assert_eq!(
            reference.len(),
            candidate.len(),
            "probe {q} length mismatch"
        );
        let diff: Vec<f64> = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| a - b)
            .collect();
        let rel_l2 = l2(&diff) / family_l2;
        let rel_linf = linf(&diff) / family_linf;
        eprintln!(
            "compute-015: probe {q}: family-rel L2 = {rel_l2:.3e}, family-rel L∞ = {rel_linf:.3e}"
        );
        assert!(
            rel_l2 < REL_L2_TOL,
            "compute-015: probe {q} family-rel L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
        );
        assert!(
            rel_linf < REL_LINF_TOL,
            "compute-015: probe {q} family-rel L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
        );
    }
}
