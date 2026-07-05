//! Gate `compute-002`: the wgpu FP32 GPU backend matches the FP64 CPU
//! backend within FP32 accumulation tolerance on the E.0 vacuum scenario.
//!
//! **Self-skipping:** when no wgpu adapter exists (typical hosted CI), the
//! test prints a SKIPPED notice and returns green. The GPU nightly runner is
//! where this gate actually bites — a green default-CI run is NOT evidence
//! the GPU path works (same posture as the CUDA lane, CLAUDE.md §10).

#![cfg(feature = "gpu")]

use yee_compute::{ComputeError, CpuFdtd, FdtdSpec, Fields, GpuFdtd};

const NX: usize = 20;
const NY: usize = 18;
const NZ: usize = 16;
const DX: f64 = 1e-3;
const STEPS: usize = 100;

/// Relative L2 tolerance per component (FP32 round-off accumulated over
/// `STEPS` leapfrog updates; measured ~4e-7 on llvmpipe, so ~250× headroom).
const REL_L2_TOL: f64 = 1e-4;
/// Normalized L∞ tolerance per component.
const REL_LINF_TOL: f64 = 1e-3;

/// L2 norm of a slice.
fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// L∞ norm of a slice.
fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

#[test]
fn gpu_matches_cpu_within_fp32_tolerance() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, NZ / 2), 2.0);

    let mut gpu = match GpuFdtd::new(spec, init.clone()) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!(
                "SKIPPED compute-002: no wgpu adapter on this machine; \
                 the GPU nightly runner exercises this gate for real"
            );
            return;
        }
        Err(other) => panic!("compute-002: GPU backend construction failed: {other}"),
    };
    eprintln!("compute-002: running on adapter '{}'", gpu.adapter_name());

    let mut cpu = CpuFdtd::new(spec, init);
    cpu.step_n(STEPS);
    gpu.step_n(STEPS).expect("GPU stepping failed");
    let gpu_fields = gpu.read_fields().expect("GPU readback failed");
    let cpu_fields = cpu.fields();

    // Normalize each component against its *family* (E or H) norm: with a
    // pure-E_z initial pulse H_z is excited only at second order (its norm
    // is ~10³ smaller than H_x's here), so a per-component relative test
    // would amplify ordinary FP32 round-off into a spurious failure. E and H
    // carry different units (factor η₀), so the families stay separate.
    let e_l2 = l2(&[&cpu_fields.ex[..], &cpu_fields.ey, &cpu_fields.ez].concat());
    let h_l2 = l2(&[&cpu_fields.hx[..], &cpu_fields.hy, &cpu_fields.hz].concat());
    let e_linf = linf(&cpu_fields.ex)
        .max(linf(&cpu_fields.ey))
        .max(linf(&cpu_fields.ez));
    let h_linf = linf(&cpu_fields.hx)
        .max(linf(&cpu_fields.hy))
        .max(linf(&cpu_fields.hz));
    assert!(
        e_l2 > 0.0 && h_l2 > 0.0,
        "compute-002: the CPU reference never energized E or H — scenario broken"
    );

    struct Check<'a> {
        reference: &'a [f64],
        candidate: &'a [f64],
        family_l2: f64,
        family_linf: f64,
        name: &'static str,
    }
    let checks = [
        ("ex", &cpu_fields.ex, &gpu_fields.ex, e_l2, e_linf),
        ("ey", &cpu_fields.ey, &gpu_fields.ey, e_l2, e_linf),
        ("ez", &cpu_fields.ez, &gpu_fields.ez, e_l2, e_linf),
        ("hx", &cpu_fields.hx, &gpu_fields.hx, h_l2, h_linf),
        ("hy", &cpu_fields.hy, &gpu_fields.hy, h_l2, h_linf),
        ("hz", &cpu_fields.hz, &gpu_fields.hz, h_l2, h_linf),
    ]
    .map(
        |(name, reference, candidate, family_l2, family_linf)| Check {
            reference,
            candidate,
            family_l2,
            family_linf,
            name,
        },
    );
    for Check {
        reference,
        candidate,
        family_l2,
        family_linf,
        name,
    } in checks
    {
        assert_eq!(reference.len(), candidate.len(), "{name} length mismatch");
        let diff: Vec<f64> = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| a - b)
            .collect();
        let rel_l2 = l2(&diff) / family_l2;
        let rel_linf = linf(&diff) / family_linf;
        eprintln!(
            "compute-002: {name}: family-rel L2 = {rel_l2:.3e}, family-rel L∞ = {rel_linf:.3e}"
        );
        assert!(
            rel_l2 < REL_L2_TOL,
            "compute-002: {name} family-rel L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
        );
        assert!(
            rel_linf < REL_LINF_TOL,
            "compute-002: {name} family-rel L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
        );
    }
}
