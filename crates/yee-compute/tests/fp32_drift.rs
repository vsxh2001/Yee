//! Gate `compute-009` (E.3 precision policy): FP32 accumulation drift of the
//! GPU backend over a **long** run — 10⁴ leapfrog steps in a closed PEC box
//! (energy-conserving, so round-off cannot hide in decay) — bounded against
//! the FP64 CPU backend.
//!
//! Round-off behaves as a random walk: per-step relative error is O(1e-7)
//! (FP32 ulp), so 10⁴ steps predicts O(√10⁴·1e-7) = O(1e-5) family-relative
//! drift; the gate allows 100× headroom over the E.0/E.1 short-run numbers
//! (~4e-7 at 100 steps) with 1e-3 L2 / 5e-3 L∞. Measured on llvmpipe when
//! this gate landed: see ADR-0177.
//!
//! Self-skipping without an adapter; runs on the GPU nightly.
//!
//! ```bash
//! cargo test -p yee-compute --release --test fp32_drift -- --ignored --nocapture
//! ```

#![cfg(feature = "gpu")]

use yee_compute::{Boundary, ComputeError, CpuFdtd, FdtdSpec, Fields, GpuFdtd, Materials};

const NX: usize = 20;
const NY: usize = 18;
const NZ: usize = 16;
const DX: f64 = 1e-3;
const STEPS: usize = 10_000;

const REL_L2_TOL: f64 = 1e-3;
const REL_LINF_TOL: f64 = 5e-3;

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

#[test]
#[ignore = "slow: 10^4-step GPU drift characterization (E.3); run with --release --ignored"]
fn fp32_drift_over_ten_thousand_steps_is_bounded() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, NZ / 2), 2.0);

    let mut gpu =
        match GpuFdtd::with_config(spec, init.clone(), Materials::default(), Boundary::PecBox) {
            Ok(gpu) => gpu,
            Err(ComputeError::NoAdapter) => {
                eprintln!("SKIPPED compute-009: no wgpu adapter");
                return;
            }
            Err(other) => panic!("compute-009: GPU construction failed: {other}"),
        };
    eprintln!("compute-009: adapter '{}'", gpu.adapter_name());

    let mut cpu = CpuFdtd::with_config(spec, init, Materials::default(), Boundary::PecBox);
    cpu.step_n(STEPS);
    gpu.step_n(STEPS).expect("GPU stepping failed");
    let g = gpu.read_fields().expect("readback failed");
    let c = cpu.fields();

    let e_l2 = l2(&[&c.ex[..], &c.ey, &c.ez].concat());
    let h_l2 = l2(&[&c.hx[..], &c.hy, &c.hz].concat());
    let e_linf = linf(&c.ex).max(linf(&c.ey)).max(linf(&c.ez));
    let h_linf = linf(&c.hx).max(linf(&c.hy)).max(linf(&c.hz));
    assert!(
        e_l2 > 0.0 && h_l2 > 0.0,
        "reference lost the field — broken"
    );

    type Check<'a> = (&'a str, &'a [f64], &'a [f64], f64, f64);
    let checks: [Check; 6] = [
        ("ex", &c.ex, &g.ex, e_l2, e_linf),
        ("ey", &c.ey, &g.ey, e_l2, e_linf),
        ("ez", &c.ez, &g.ez, e_l2, e_linf),
        ("hx", &c.hx, &g.hx, h_l2, h_linf),
        ("hy", &c.hy, &g.hy, h_l2, h_linf),
        ("hz", &c.hz, &g.hz, h_l2, h_linf),
    ];
    for (name, reference, candidate, fam_l2, fam_linf) in checks {
        let diff: Vec<f64> = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| a - b)
            .collect();
        let rel_l2 = l2(&diff) / fam_l2;
        let rel_linf = linf(&diff) / fam_linf;
        eprintln!(
            "compute-009: {name}: drift after {STEPS} steps — family-rel L2 = {rel_l2:.3e}, \
             L∞ = {rel_linf:.3e}"
        );
        assert!(
            rel_l2 < REL_L2_TOL,
            "compute-009: {name} drift L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
        );
        assert!(
            rel_linf < REL_LINF_TOL,
            "compute-009: {name} drift L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
        );
    }
}
