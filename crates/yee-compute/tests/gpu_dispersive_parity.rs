//! Gate `compute-012` (E.5c): the GPU unified-ADE dispersive update matches
//! the FP64 CPU backend (itself bit-exact vs `yee_fdtd::dispersive`, gate
//! `compute-011`) on a four-arm scenario — Drude + Lorentz + Debye blocks in
//! vacuum, Gaussian E_z ball in a PEC box, 400 steps.
//!
//! **Differential gate.** This scenario's family-relative normalization is
//! structurally harsh on H: the dispersive blocks absorb the traveling wave
//! while the ball's electrostatic E residual persists, so ‖H‖ stays ~η₀
//! below ‖E‖ and ordinary E-scale FP32 noise (~1e-7 absolute — verified via
//! a standard-path control during E.5c bring-up) reads as ~2e-4 relative to
//! H no matter how few steps run. The gate therefore asserts the thing
//! E.5c is actually responsible for: the ADE pair's GPU↔CPU error must be
//! within a bounded factor of the *standard* pair's error on the identical
//! scenario (per component; measured ≤ 6×, gated at 20×), plus an absolute drift-class
//! backstop (1e-3 L2 / 5e-3 L∞, the compute-009 calibration).
//!
//! Self-skipping without a wgpu adapter; real on llvmpipe and the nightly.

#![cfg(feature = "gpu")]

use yee_compute::{
    Boundary, ComputeError, CpuFdtd, DispersiveMap, DispersiveMaterial, Drive, FdtdSpec, Fields,
    GpuFdtd, Materials,
};

const NX: usize = 18;
const NY: usize = 16;
const NZ: usize = 14;
const DX: f64 = 1e-3;
const STEPS: usize = 400;

const ABS_L2_TOL: f64 = 1e-3;
const ABS_LINF_TOL: f64 = 5e-3;
/// Measured on llvmpipe at E.5c bring-up: the ADE pair runs ≤ 6× the
/// standard pair (the f32 aux-state recursions — Lorentz α ≈ 1.99 —
/// amplify rounding relative to the CPU's f64 aux state; the all-vacuum
/// ADE path is bit-identical to the standard path, ruling out indexing
/// bugs, which would show O(1) errors thousands of times the control).
const DIFFERENTIAL_FACTOR: f64 = 20.0;
/// Floor for the standard-pair error so a fortuitously-exact standard run
/// cannot make the differential bound impossible to meet.
const STD_FLOOR: f64 = 1e-6;

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

fn dispersive_map(spec: &FdtdSpec) -> DispersiveMap {
    let mut map = DispersiveMap::vacuum(spec);
    map.set_box(
        spec,
        2,
        7,
        2,
        7,
        2,
        7,
        DispersiveMaterial::Drude {
            eps_inf: 1.0,
            omega_p: 2.0e10 * std::f64::consts::TAU,
            gamma: 1.0e9,
        },
    );
    map.set_box(
        spec,
        10,
        16,
        3,
        9,
        4,
        10,
        DispersiveMaterial::Lorentz {
            eps_inf: 2.0,
            delta_eps: 1.5,
            omega_0: 1.0e10 * std::f64::consts::TAU,
            delta: 2.0e8,
        },
    );
    map.set_box(
        spec,
        4,
        12,
        10,
        15,
        6,
        12,
        DispersiveMaterial::Debye {
            eps_inf: 1.5,
            delta_eps: 8.0,
            tau: 5.0e-11,
        },
    );
    map
}

/// Per-component (L2, L∞) errors of `candidate` vs `reference`, normalized
/// by the reference's E/H family norms. Order: ex, ey, ez, hx, hy, hz.
fn family_rel_errors(reference: &Fields, candidate: &Fields) -> [(f64, f64); 6] {
    let e_l2 = l2(&[&reference.ex[..], &reference.ey, &reference.ez].concat());
    let h_l2 = l2(&[&reference.hx[..], &reference.hy, &reference.hz].concat());
    let e_linf = linf(&reference.ex)
        .max(linf(&reference.ey))
        .max(linf(&reference.ez));
    let h_linf = linf(&reference.hx)
        .max(linf(&reference.hy))
        .max(linf(&reference.hz));
    assert!(e_l2 > 0.0 && h_l2 > 0.0, "reference never energized");

    let pairs: [(&[f64], &[f64], f64, f64); 6] = [
        (&reference.ex, &candidate.ex, e_l2, e_linf),
        (&reference.ey, &candidate.ey, e_l2, e_linf),
        (&reference.ez, &candidate.ez, e_l2, e_linf),
        (&reference.hx, &candidate.hx, h_l2, h_linf),
        (&reference.hy, &candidate.hy, h_l2, h_linf),
        (&reference.hz, &candidate.hz, h_l2, h_linf),
    ];
    pairs.map(|(r, c, fam_l2, fam_linf)| {
        let diff: Vec<f64> = r.iter().zip(c).map(|(a, b)| a - b).collect();
        (l2(&diff) / fam_l2, linf(&diff) / fam_linf)
    })
}

#[test]
fn gpu_dispersive_matches_cpu_within_fp32_tolerance() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, NZ / 2), 2.0);
    let map = dispersive_map(&spec);

    // ---- ADE pair ----
    let mut gpu_ade = match GpuFdtd::with_dispersive(
        spec,
        init.clone(),
        Materials::default(),
        Boundary::PecBox,
        Drive::default(),
        0,
        &map,
    ) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => {
            eprintln!("SKIPPED compute-012: no wgpu adapter");
            return;
        }
        Err(other) => panic!("compute-012: GPU construction failed: {other}"),
    };
    eprintln!("compute-012: adapter '{}'", gpu_ade.adapter_name());
    let mut cpu_ade =
        CpuFdtd::with_config(spec, init.clone(), Materials::default(), Boundary::PecBox);
    cpu_ade.set_dispersive(map);
    gpu_ade.step_n(STEPS).expect("GPU ADE stepping failed");
    cpu_ade.step_n(STEPS);
    let ade_err = family_rel_errors(cpu_ade.fields(), &gpu_ade.read_fields().unwrap());

    // ---- standard-pair control on the identical scenario ----
    let mut gpu_std =
        GpuFdtd::with_config(spec, init.clone(), Materials::default(), Boundary::PecBox)
            .expect("adapter vanished");
    let mut cpu_std = CpuFdtd::with_config(spec, init, Materials::default(), Boundary::PecBox);
    gpu_std.step_n(STEPS).expect("GPU std stepping failed");
    cpu_std.step_n(STEPS);
    let std_err = family_rel_errors(cpu_std.fields(), &gpu_std.read_fields().unwrap());

    let names = ["ex", "ey", "ez", "hx", "hy", "hz"];
    for ((name, (ade_l2, ade_linf)), (std_l2, std_linf)) in names.iter().zip(ade_err).zip(std_err) {
        eprintln!(
            "compute-012: {name}: ADE-pair L2 = {ade_l2:.3e} (std {std_l2:.3e}), \
             L∞ = {ade_linf:.3e} (std {std_linf:.3e})"
        );
        assert!(
            ade_l2 <= DIFFERENTIAL_FACTOR * std_l2.max(STD_FLOOR),
            "compute-012: {name} ADE L2 {ade_l2:e} exceeds {DIFFERENTIAL_FACTOR}× the \
             standard-pair control {std_l2:e}"
        );
        assert!(
            ade_linf <= DIFFERENTIAL_FACTOR * std_linf.max(STD_FLOOR),
            "compute-012: {name} ADE L∞ {ade_linf:e} exceeds {DIFFERENTIAL_FACTOR}× the \
             standard-pair control {std_linf:e}"
        );
        assert!(
            ade_l2 < ABS_L2_TOL && ade_linf < ABS_LINF_TOL,
            "compute-012: {name} ADE error exceeds the absolute backstop \
             ({ade_l2:e} / {ade_linf:e} vs {ABS_L2_TOL:e} / {ABS_LINF_TOL:e})"
        );
    }
}
