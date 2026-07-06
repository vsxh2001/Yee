//! Gate `compute-005`: the wgpu FP32 backend matches the FP64 CPU backend
//! on the full E.1 scenario — CPML absorbing boundaries, per-cell ε_r / σ /
//! μ_r, and an interior PEC sheet with a slot — within FP32 accumulation
//! tolerance (family-relative norms, as `compute-002`). Also asserts the
//! CPML actually absorbs on the GPU: after the pulse leaves the domain, the
//! CPML run holds far less field energy than the same run in a PEC box.
//!
//! **Self-skipping** without a wgpu adapter (see `compute-002`); real on
//! the GPU nightly runner and on Mesa llvmpipe.

#![cfg(feature = "gpu")]

use yee_compute::{
    Boundary, ComputeError, CpmlConfig, CpuFdtd, FdtdSpec, Fields, GpuFdtd, Materials,
};

const NX: usize = 24;
const NY: usize = 22;
const NZ: usize = 20;
const DX: f64 = 1e-3;
const NPML: usize = 6;
const STEPS: usize = 100;

const REL_L2_TOL: f64 = 1e-4;
const REL_LINF_TOL: f64 = 1e-3;

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

fn materials() -> Materials {
    let celld = (NX + 1, NY + 1, NZ + 1);
    let n_cells = celld.0 * celld.1 * celld.2;
    let cell = |i: usize, j: usize, k: usize| (i * celld.1 + j) * celld.2 + k;
    let mut eps_r = vec![1.0; n_cells];
    let mut mu_r = vec![1.0; n_cells];
    let mut sigma = vec![0.0; n_cells];
    for i in 0..celld.0 {
        for j in 0..celld.1 {
            for k in 0..celld.2 {
                if (3..8).contains(&k) {
                    eps_r[cell(i, j, k)] = 4.3;
                }
                if (8..14).contains(&i) && (8..14).contains(&j) && (9..13).contains(&k) {
                    sigma[cell(i, j, k)] = 0.5;
                }
                if (14..18).contains(&j) {
                    mu_r[cell(i, j, k)] = 2.0;
                }
            }
        }
    }
    // PEC sheet at k = 14 on E_x with a slot at j ∈ [9, 13).
    let exd = (NX, NY + 1, NZ + 1);
    let mut mask_ex = vec![false; exd.0 * exd.1 * exd.2];
    for i in 0..exd.0 {
        for j in 0..exd.1 {
            if !(9..13).contains(&j) {
                mask_ex[(i * exd.1 + j) * exd.2 + 14] = true;
            }
        }
    }
    Materials {
        eps_r_cells: Some(eps_r),
        mu_r_cells: Some(mu_r),
        sigma_cells: Some(sigma),
        pec_mask_ex: Some(mask_ex),
        pec_mask_ey: None,
        pec_mask_ez: None,
    }
}

#[test]
fn gpu_matches_cpu_on_cpml_heterogeneous_scenario() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, NZ / 2), 2.0);
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML));

    let mut gpu = match GpuFdtd::with_config(spec, init.clone(), materials(), boundary.clone()) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!(
                "SKIPPED compute-005: no wgpu adapter on this machine; \
                 the GPU nightly runner exercises this gate for real"
            );
            return;
        }
        Err(other) => panic!("compute-005: GPU backend construction failed: {other}"),
    };
    eprintln!("compute-005: running on adapter '{}'", gpu.adapter_name());

    let mut cpu = CpuFdtd::with_config(spec, init.clone(), materials(), boundary);
    cpu.step_n(STEPS);
    gpu.step_n(STEPS).expect("GPU stepping failed");
    let gpu_fields = gpu.read_fields().expect("GPU readback failed");
    let cpu_fields = cpu.fields();

    // Family-relative comparison; see compute-002 for the rationale.
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
        "compute-005: the CPU reference never energized E or H — scenario broken"
    );

    type Check<'a> = (&'a str, &'a [f64], &'a [f64], f64, f64);
    let checks: [Check; 6] = [
        ("ex", &cpu_fields.ex, &gpu_fields.ex, e_l2, e_linf),
        ("ey", &cpu_fields.ey, &gpu_fields.ey, e_l2, e_linf),
        ("ez", &cpu_fields.ez, &gpu_fields.ez, e_l2, e_linf),
        ("hx", &cpu_fields.hx, &gpu_fields.hx, h_l2, h_linf),
        ("hy", &cpu_fields.hy, &gpu_fields.hy, h_l2, h_linf),
        ("hz", &cpu_fields.hz, &gpu_fields.hz, h_l2, h_linf),
    ];
    for (name, reference, candidate, family_l2, family_linf) in checks {
        assert_eq!(reference.len(), candidate.len(), "{name} length mismatch");
        let diff: Vec<f64> = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| a - b)
            .collect();
        let rel_l2 = l2(&diff) / family_l2;
        let rel_linf = linf(&diff) / family_linf;
        eprintln!(
            "compute-005: {name}: family-rel L2 = {rel_l2:.3e}, family-rel L∞ = {rel_linf:.3e}"
        );
        assert!(
            rel_l2 < REL_L2_TOL,
            "compute-005: {name} family-rel L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
        );
        assert!(
            rel_linf < REL_LINF_TOL,
            "compute-005: {name} family-rel L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
        );
    }

    // Absorption evidence on the GPU, measured on the **H family only**:
    // the initial E_z ball carries a large electrostatic (curl-free)
    // residual that never propagates — no boundary can absorb it, so the
    // total-field norm barely drops. A static charge holds no H, though:
    // once the radiating part is absorbed, the CPML run's H collapses,
    // while the PEC box keeps the wave bouncing.
    let spec_vac = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init_vac = Fields::with_gaussian_ez(&spec_vac, (NX / 2, NY / 2, NZ / 2), 2.0);
    let h_norm = |f: &Fields| l2(&[&f.hx[..], &f.hy, &f.hz].concat());
    let run = |boundary: Boundary| -> Option<f64> {
        let mut engine =
            GpuFdtd::with_config(spec_vac, init_vac.clone(), Materials::default(), boundary)
                .ok()?;
        engine.step_n(STEPS).ok()?;
        Some(h_norm(&engine.read_fields().ok()?))
    };
    if let (Some(pec_h), Some(cpml_h)) = (
        run(Boundary::PecBox),
        run(Boundary::Cpml(CpmlConfig::for_spec(&spec_vac, NPML))),
    ) {
        eprintln!("compute-005: post-run ‖H‖₂ — PEC box {pec_h:.3e}, CPML {cpml_h:.3e}");
        assert!(
            cpml_h < 0.1 * pec_h,
            "compute-005: GPU CPML absorbed too little (‖H‖ CPML {cpml_h:e} vs PEC {pec_h:e})"
        );
    }
}
