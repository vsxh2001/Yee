//! Gate `compute-004`: the E.1 CPML reproduces the `yee-fdtd` reflection
//! criterion — **≥ 30 dB** reduction of the reflected-wave peak versus a
//! reflecting PEC box — running entirely on `yee-compute`'s CPU backend.
//!
//! Same methodology as `crates/yee-fdtd/tests/cpml_reflection.rs`: 50³
//! vacuum grid, Gaussian-in-time soft source at the centre, `E_z` probed at
//! (38, 25, 25) for 300 steps. The PEC-reflected wave is isolated as
//! `max|pec(t) − cpml(t)|` (static source residual and the outgoing pulse
//! cancel exactly); the CPML's own residual is the late-window oscillation
//! around its static floor.

use yee_compute::{Boundary, CpmlConfig, CpuFdtd, FdtdSpec, Fields, Materials};

const N: usize = 50;
const DX: f64 = 1.0e-3;
const NPML: usize = 10;
const N_STEPS: usize = 300;
const SOURCE: (usize, usize, usize) = (25, 25, 25);
const PROBE: (usize, usize, usize) = (38, 25, 25);

fn run_trace(boundary: Boundary, t0: f64, sigma: f64) -> Vec<f64> {
    let spec = FdtdSpec::vacuum(N, N, N, DX);
    let mut engine =
        CpuFdtd::with_config(spec, Fields::zero(&spec), Materials::default(), boundary);
    let ezd = spec.ez_dims();
    let probe_idx = (PROBE.0 * ezd.1 + PROBE.1) * ezd.2 + PROBE.2;
    let mut trace = Vec::with_capacity(N_STEPS);
    for _ in 0..N_STEPS {
        engine.step_with_gaussian_ez(SOURCE, t0, sigma);
        trace.push(engine.fields().ez[probe_idx]);
    }
    trace
}

#[test]
fn cpml_attenuates_reflection_vs_pec() {
    let spec = FdtdSpec::vacuum(N, N, N, DX);
    let dt = spec.dt;
    let t0 = 20.0 * dt;
    let sigma = 6.0 * dt;

    let pec_trace = run_trace(Boundary::PecBox, t0, sigma);
    let cpml_trace = run_trace(Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML)), t0, sigma);

    assert!(
        pec_trace.iter().all(|x| x.is_finite()),
        "PEC trace went non-finite"
    );
    assert!(
        cpml_trace.iter().all(|x| x.is_finite()),
        "CPML trace went non-finite"
    );

    // PEC reflection: the traces are identical until a wall reflection
    // returns, so their difference isolates the reflected wave.
    let pec_reflection_peak = pec_trace
        .iter()
        .zip(cpml_trace.iter())
        .map(|(p, c)| (p - c).abs())
        .fold(0.0_f64, f64::max);

    // CPML residual: late-window oscillation around the (DC-like) static
    // source residual, estimated by the late-window mean.
    const REFLECTION_START: usize = 80;
    let late_cpml = &cpml_trace[REFLECTION_START..];
    let static_floor = late_cpml.iter().copied().sum::<f64>() / (late_cpml.len() as f64);
    let cpml_reflection_peak = late_cpml
        .iter()
        .map(|x| (x - static_floor).abs())
        .fold(0.0_f64, f64::max);

    let peak_outgoing = pec_trace
        .iter()
        .chain(cpml_trace.iter())
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);

    let pec_db = 20.0 * (pec_reflection_peak / peak_outgoing).log10();
    let cpml_db = 20.0 * (cpml_reflection_peak / peak_outgoing).log10();
    let reduction_db = pec_db - cpml_db;

    eprintln!("compute-004: outgoing peak        = {peak_outgoing:.3e}");
    eprintln!("compute-004: PEC reflection peak  = {pec_reflection_peak:.3e} ({pec_db:.2} dB)");
    eprintln!("compute-004: CPML reflection peak = {cpml_reflection_peak:.3e} ({cpml_db:.2} dB)");
    eprintln!("compute-004: reduction            = {reduction_db:.2} dB");

    assert!(peak_outgoing > 0.0, "no outgoing pulse seen");
    assert!(
        pec_reflection_peak > 0.0,
        "PEC and CPML are identical (no reflection seen)"
    );
    assert!(
        cpml_reflection_peak > 0.0,
        "CPML reflection peak is zero (test logic error)"
    );
    assert!(
        reduction_db >= 30.0,
        "compute-004: CPML reflection reduction {reduction_db:.2} dB is below the 30 dB target"
    );
}
