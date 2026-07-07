//! Gate `compute-016` (R.3, ADR-0196): the wgpu backend honors **per-face
//! CPML masks** — the antenna-track boundary (A.2: open top over a PEC
//! ground, `faces = [[t,t],[t,t],[f,t]]`) the per-axis shader mask could
//! not express. Full final field state compared family-relative against
//! the FP64 CPU backend (the compute-005 idiom), plus evidence the
//! enabled z-max face actually absorbs (H-family norm far below a PEC-box
//! run) while the disabled z-min face stays reflective.
//!
//! **Self-skipping** without a wgpu adapter; real on Mesa llvmpipe and the
//! GPU nightly runner.

#![cfg(feature = "gpu")]

use yee_compute::{
    Boundary, ComputeError, CpmlConfig, CpuFdtd, FdtdSpec, Fields, GpuFdtd, Materials,
};

const NX: usize = 24;
const NY: usize = 22;
const NZ: usize = 20;
const DX: f64 = 1e-3;
const NPML: usize = 6;
const STEPS: usize = 150;

const REL_L2_TOL: f64 = 1e-4;
const REL_LINF_TOL: f64 = 1e-3;

/// Open top over a reflective floor: every face absorbs except z-min.
const FACES: [[bool; 2]; 3] = [[true, true], [true, true], [false, true]];

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

#[test]
fn gpu_matches_cpu_on_per_face_cpml() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    // Source ball low over the floor so both the reflective z-min face and
    // the absorbing z-max face see real energy.
    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, 6), 2.0);
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML).with_faces(FACES));

    let mut gpu =
        match GpuFdtd::with_config(spec, init.clone(), Materials::default(), boundary.clone()) {
            Ok(gpu) => gpu,
            Err(ComputeError::NoAdapter) => {
                eprintln!(
                    "SKIPPED compute-016: no wgpu adapter on this machine; \
                 llvmpipe / the GPU nightly runner exercise this gate for real"
                );
                return;
            }
            Err(other) => panic!("compute-016: GPU backend construction failed: {other}"),
        };
    eprintln!("compute-016: running on adapter '{}'", gpu.adapter_name());

    let mut cpu = CpuFdtd::with_config(spec, init.clone(), Materials::default(), boundary);
    cpu.step_n(STEPS);
    gpu.step_n(STEPS).expect("GPU stepping failed");
    let gpu_fields = gpu.read_fields().expect("GPU readback failed");
    let cpu_fields = cpu.fields();

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
        "compute-016: the CPU reference never energized E or H — scenario broken"
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
            "compute-016: {name}: family-rel L2 = {rel_l2:.3e}, family-rel L∞ = {rel_linf:.3e}"
        );
        assert!(
            rel_l2 < REL_L2_TOL,
            "compute-016: {name} family-rel L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
        );
        assert!(
            rel_linf < REL_LINF_TOL,
            "compute-016: {name} family-rel L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
        );
    }

    // Absorption evidence (H family, as compute-005): the per-face CPML run
    // must shed most of its radiating energy through the five open faces;
    // a PEC box keeps it bouncing.
    let h_norm = |f: &Fields| l2(&[&f.hx[..], &f.hy, &f.hz].concat());
    let run_pec = || -> Option<f64> {
        let mut engine =
            GpuFdtd::with_config(spec, init.clone(), Materials::default(), Boundary::PecBox)
                .ok()?;
        engine.step_n(STEPS).ok()?;
        Some(h_norm(&engine.read_fields().ok()?))
    };
    if let Some(pec_h) = run_pec() {
        let open_h = h_norm(&gpu_fields);
        eprintln!("compute-016: post-run ‖H‖₂ — PEC box {pec_h:.3e}, open-top CPML {open_h:.3e}");
        assert!(
            open_h < 0.2 * pec_h,
            "compute-016: per-face CPML absorbed too little (‖H‖ {open_h:e} vs PEC {pec_h:e})"
        );
    }
}
