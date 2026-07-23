//! Parity gate (FS.4.2a, deliverable 1): the GPU H-probe recording path
//! (`record_h_probes` WGSL kernel / `GpuFdtd::read_h_probes`) matches the
//! CPU H-probe series within FP32 accumulation tolerance — same posture and
//! same tolerances as `gpu_cpu_parity.rs` (compute-002).
//!
//! **Self-skipping:** when no wgpu adapter exists (typical hosted CI), the
//! test prints a SKIPPED notice and returns green. The GPU nightly runner is
//! where this gate actually bites — a green default-CI run is NOT evidence
//! the GPU H-probe path works (CLAUDE.md §10).

#![cfg(feature = "gpu")]

use yee_compute::{
    Boundary, ComputeError, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, GpuFdtd, HComponent,
    HProbe, Materials, Probe,
};

const NX: usize = 20;
const NY: usize = 18;
const NZ: usize = 16;
const DX: f64 = 1e-3;
const CENTER: (usize, usize, usize) = (NX / 2, NY / 2, NZ / 2);
const STEPS: usize = 100;

/// Same FP32 idiom as `gpu_cpu_parity.rs` (compute-002): relative L2/L∞
/// tolerance against each probe's own series norm.
const REL_L2_TOL: f64 = 1e-4;
const REL_LINF_TOL: f64 = 1e-3;

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

/// One E probe (regression: the pre-existing E-probe GPU path) plus two H
/// probes (FS.4.2a) one cell off the source in each transverse axis, on a
/// free (undriven) vacuum Gaussian evolution — same scenario shape as
/// `gpu_cpu_parity.rs`.
fn drive() -> Drive {
    Drive {
        soft_sources: vec![],
        ports: vec![],
        aperture_ports: vec![],
        probes: vec![Probe {
            component: EComponent::Ez,
            cell: (CENTER.0, CENTER.1, CENTER.2 - 1),
        }],
        h_probes: vec![
            HProbe {
                component: HComponent::Hx,
                cell: (CENTER.0, CENTER.1 - 1, CENTER.2),
            },
            HProbe {
                component: HComponent::Hy,
                cell: (CENTER.0 - 1, CENTER.1, CENTER.2),
            },
        ],
    }
}

/// Assert `candidate` matches `reference` within the FP32 tolerance,
/// normalized by `reference`'s own L2/L∞ norm (mirrors `gpu_cpu_parity.rs`'s
/// per-family normalization, here per-probe-series since there is no
/// separate E/H "family" grouping for a handful of scalar probes).
fn assert_series_matches(name: &str, reference: &[f64], candidate: &[f64]) {
    assert_eq!(reference.len(), candidate.len(), "{name}: length mismatch");
    let ref_l2 = l2(reference);
    let ref_linf = linf(reference);
    assert!(
        ref_l2 > 0.0,
        "{name}: the CPU reference series never left zero — scenario broken"
    );
    let diff: Vec<f64> = reference
        .iter()
        .zip(candidate)
        .map(|(a, b)| a - b)
        .collect();
    let rel_l2 = l2(&diff) / ref_l2;
    let rel_linf = linf(&diff) / ref_linf;
    eprintln!("FS.4.2a {name}: rel L2 = {rel_l2:.3e}, rel L∞ = {rel_linf:.3e}");
    assert!(
        rel_l2 < REL_L2_TOL,
        "FS.4.2a {name}: rel L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
    );
    assert!(
        rel_linf < REL_LINF_TOL,
        "FS.4.2a {name}: rel L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
    );
}

#[test]
fn gpu_h_probes_match_cpu_within_fp32_tolerance() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, CENTER, 2.0);

    let mut gpu = match GpuFdtd::with_drive(
        spec,
        init.clone(),
        Materials::default(),
        Boundary::None,
        drive(),
        STEPS,
    ) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!(
                "SKIPPED FS.4.2a h-probe parity: no wgpu adapter on this machine; \
                 the GPU nightly runner exercises this gate for real"
            );
            return;
        }
        Err(other) => panic!("FS.4.2a h-probe parity: GPU backend construction failed: {other}"),
    };
    eprintln!(
        "FS.4.2a h-probe parity: running on adapter '{}'",
        gpu.adapter_name()
    );

    let mut cpu = CpuFdtd::with_drive(spec, init, Materials::default(), Boundary::None, drive());
    cpu.step_n(STEPS);
    gpu.step_n(STEPS).expect("GPU stepping failed");

    let cpu_h = cpu.h_probe_series();
    let gpu_h = gpu.read_h_probes().expect("GPU H-probe readback failed");
    assert_eq!(cpu_h.len(), 2);
    assert_eq!(gpu_h.len(), 2);
    for (q, (reference, candidate)) in cpu_h.iter().zip(&gpu_h).enumerate() {
        assert_series_matches(&format!("h-probe[{q}]"), reference, candidate);
    }

    // The pre-existing E-probe GPU path stays green alongside the new H
    // probes on the SAME drive (regression: adding `h_probes` must not
    // disturb `probes`' recording).
    let cpu_e = cpu.probe_series();
    let gpu_e = gpu.read_probes().expect("GPU E-probe readback failed");
    assert_eq!(cpu_e.len(), 1);
    assert_eq!(gpu_e.len(), 1);
    assert_series_matches("e-probe[0]", &cpu_e[0], &gpu_e[0]);
}
