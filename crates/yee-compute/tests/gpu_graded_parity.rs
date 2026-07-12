//! Gates `compute-020` / `compute-021` (FS.0b.2, ADR-0214): the graded
//! (nonuniform) grid on the wgpu FP32 GPU backend.
//!
//! - `compute-020` (`gpu_graded_uniform_parity`, fast): a uniformly filled
//!   [`GradedSpacings`] must reproduce the GPU's own scalar-dx run
//!   **bit-for-bit** (the uniform inverse fill is bit-equal to the retired
//!   scalar `inv_*` uniforms, and multiplication is a pure function of its
//!   operand bit patterns — this is the off-by-one tripwire for the
//!   primal/dual divisor indexing), and must match the FP64 CPU backend
//!   within the compute-002 tolerances.
//! - `compute-021` (`gpu_graded_taper_parity`, `#[ignore]`, release): the
//!   compute-019 taper scenario run graded on the GPU matches the CPU
//!   probe time series within FP32 tolerance, and the GPU-side grading
//!   reflection (difference against a GPU uniform reference) sits under
//!   the same −48 dB pinned floor.
//!
//! **Self-skipping:** when no wgpu adapter exists, both tests print a
//! SKIPPED notice and return green (the `gpu_cpu_parity.rs` idiom). The GPU
//! nightly runner is the certifying environment; measured numbers below are
//! from llvmpipe (LLVM 20.1.2) in this workspace.

#![cfg(feature = "gpu")]

use yee_compute::{
    Boundary, ComputeError, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, GpuFdtd,
    GradedSpacings, Materials, Probe, ResistivePort, SoftSource, Waveform,
};

/// Relative L2 tolerance per component family (the compute-002 values).
const REL_L2_TOL: f64 = 1e-4;
/// Normalized L∞ tolerance per component family.
const REL_LINF_TOL: f64 = 1e-3;

fn l2(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn linf(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

// ===================== compute-020: uniform-fill parity =====================

const NX: usize = 20;
const NY: usize = 16;
const NZ: usize = 12;
const DX: f64 = 1.0e-3;
const N_STEPS: usize = 150;

/// The compute-018 drive: soft source + resistive port + aperture port +
/// probes, on asymmetric dims to catch any axis mix-up in the spacing
/// indexing. (`record: false` — aperture recording is not on the GPU.)
fn drive() -> Drive {
    let mut drive = Drive::default();
    drive.soft_sources.push(SoftSource {
        component: EComponent::Ez,
        cell: (10, 8, 6),
        waveform: Waveform::Gaussian {
            t0: 20.0 * 1.0e-12,
            sigma: 6.0 * 1.0e-12,
        },
    });
    drive.ports.push(ResistivePort {
        cell: (6, 8, 6),
        resistance: 50.0,
        waveform: Waveform::GaussianPulse {
            v0: 1.0,
            f0: 10.0e9,
            bw: 8.0e9,
            t0_steps: 30,
        },
    });
    drive.aperture_ports.push(yee_compute::AperturePort {
        cells: (7..10)
            .flat_map(|j| (0..3).map(move |k| (14, j, k)))
            .collect(),
        n_columns: 3,
        area: 3.0 * DX * (3.0 * DX),
        height: 3.0 * DX,
        resistance: 50.0,
        waveform: Waveform::GaussianPulse {
            v0: 0.5,
            f0: 12.0e9,
            bw: 10.0e9,
            t0_steps: 25,
        },
        record: false,
    });
    for cell in [(8, 8, 6), (12, 10, 4)] {
        drive.probes.push(Probe {
            component: EComponent::Ez,
            cell,
        });
    }
    drive
}

fn uniform_spacings(spec: &FdtdSpec) -> GradedSpacings {
    GradedSpacings {
        dx: vec![spec.dx; spec.nx],
        dy: vec![spec.dy; spec.ny],
        dz: vec![spec.dz; spec.nz],
    }
}

fn gpu_run(spec: FdtdSpec, boundary: Boundary, graded: bool) -> Option<GpuFdtd> {
    let mut gpu = match GpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive(),
        N_STEPS,
    ) {
        Ok(gpu) => gpu,
        Err(ComputeError::NoAdapter) => return None,
        Err(other) => panic!("compute-020: GPU backend construction failed: {other}"),
    };
    if graded {
        gpu.set_spacings(&uniform_spacings(&spec))
            .expect("uniform spacings must be accepted");
    }
    gpu.step_n(N_STEPS).expect("GPU stepping failed");
    Some(gpu)
}

#[test]
fn gpu_graded_uniform_parity() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, 4));

    let Some(mut gpu_graded) = gpu_run(spec, boundary.clone(), true) else {
        eprintln!(
            "SKIPPED compute-020: no wgpu adapter on this machine; \
             the GPU nightly runner exercises this gate for real"
        );
        return;
    };
    eprintln!(
        "compute-020: running on adapter '{}'",
        gpu_graded.adapter_name()
    );
    let mut gpu_scalar = gpu_run(spec, boundary.clone(), false).expect("adapter vanished mid-test");

    // --- (a) GPU graded-uniform vs GPU scalar: BIT-FOR-BIT. The widened
    // f64 values are equal iff the f32 bits are (widening is injective on
    // finite values), so exact f64 equality is the bit-exactness assert.
    let fg = gpu_graded.read_fields().expect("GPU readback failed");
    let fs = gpu_scalar.read_fields().expect("GPU readback failed");
    for (name, a, b) in [
        ("ex", &fg.ex, &fs.ex),
        ("ey", &fg.ey, &fs.ey),
        ("ez", &fg.ez, &fs.ez),
        ("hx", &fg.hx, &fs.hx),
        ("hy", &fg.hy, &fs.hy),
        ("hz", &fg.hz, &fs.hz),
    ] {
        assert!(
            a.iter().all(|v| v.is_finite()),
            "compute-020: {name} went non-finite on the graded path"
        );
        let n_diff = a.iter().zip(b.iter()).filter(|(u, v)| u != v).count();
        assert!(
            n_diff == 0,
            "compute-020: {name}: {n_diff} elements differ between the \
             uniform-filled graded run and the scalar run (must be bit-for-bit)"
        );
    }
    let pg = gpu_graded.read_probes().expect("GPU probe readback failed");
    let ps = gpu_scalar.read_probes().expect("GPU probe readback failed");
    for (q, (a, b)) in pg.iter().zip(&ps).enumerate() {
        assert_eq!(a.len(), N_STEPS, "probe {q}: wrong series length");
        assert!(
            a.iter().any(|v| *v != 0.0),
            "compute-020: probe {q} stayed silent (test logic error)"
        );
        for (n, (u, v)) in a.iter().zip(b.iter()).enumerate() {
            assert!(
                u == v,
                "compute-020: probe {q} step {n} diverged (graded {u:e}, scalar {v:e})"
            );
        }
    }
    eprintln!("compute-020: graded-uniform vs scalar GPU run: bit-for-bit PASS");

    // --- (b) GPU graded-uniform vs CPU FP64: compute-002 tolerances,
    // family-normalized (see gpu_cpu_parity.rs for the rationale).
    let mut cpu = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive(),
    );
    cpu.set_spacings(&uniform_spacings(&spec));
    cpu.step_n(N_STEPS);
    let fc = cpu.fields();
    let e_l2 = l2(&[&fc.ex[..], &fc.ey, &fc.ez].concat());
    let h_l2 = l2(&[&fc.hx[..], &fc.hy, &fc.hz].concat());
    let e_linf = linf(&fc.ex).max(linf(&fc.ey)).max(linf(&fc.ez));
    let h_linf = linf(&fc.hx).max(linf(&fc.hy)).max(linf(&fc.hz));
    assert!(
        e_l2 > 0.0 && h_l2 > 0.0,
        "compute-020: the CPU reference never energized E or H — scenario broken"
    );
    for (name, reference, candidate, family_l2, family_linf) in [
        ("ex", &fc.ex, &fg.ex, e_l2, e_linf),
        ("ey", &fc.ey, &fg.ey, e_l2, e_linf),
        ("ez", &fc.ez, &fg.ez, e_l2, e_linf),
        ("hx", &fc.hx, &fg.hx, h_l2, h_linf),
        ("hy", &fc.hy, &fg.hy, h_l2, h_linf),
        ("hz", &fc.hz, &fg.hz, h_l2, h_linf),
    ] {
        let diff: Vec<f64> = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| a - b)
            .collect();
        let rel_l2 = l2(&diff) / family_l2;
        let rel_linf = linf(&diff) / family_linf;
        eprintln!(
            "compute-020: {name}: family-rel L2 = {rel_l2:.3e}, family-rel L∞ = {rel_linf:.3e}"
        );
        assert!(
            rel_l2 < REL_L2_TOL,
            "compute-020: {name} family-rel L2 {rel_l2:e} ≥ {REL_L2_TOL:e}"
        );
        assert!(
            rel_linf < REL_LINF_TOL,
            "compute-020: {name} family-rel L∞ {rel_linf:e} ≥ {REL_LINF_TOL:e}"
        );
    }
}

// ===================== compute-021: taper-scenario parity =====================
// The compute-019 setup (graded_interface_reflection.rs), verbatim: a pulse
// crosses a 0.5 → 0.25 → 0.5 mm geometric taper (ratio 2^(1/6) ≈ 1.122 per
// cell) in free space, CPML on every face, dt shared from the graded
// Courant limit.

const COARSE: f64 = 0.5e-3;
const FINE: f64 = 0.25e-3;
const NPML: usize = 10;
const N_TAPER: usize = 6;
const TNY: usize = 40;
const TNZ: usize = 40;
const T_STEPS: usize = 560;
const SOURCE_I: usize = 16;
const PROBE_I: usize = 34;

/// CPU↔GPU probe-series tolerance, normalized by the CPU trace peak / L2.
/// Measured 2026-07-12 on llvmpipe (release): rel L∞ = 4.709e-6, rel
/// L2 = 5.853e-6 over the 560-step series (FP32 accumulation through
/// CPML). Pinned at 1e-4 (~20× headroom); do not weaken without
/// re-measuring.
const SERIES_REL_TOL: f64 = 1.0e-4;

fn graded_dx() -> Vec<f64> {
    let r = (FINE / COARSE).powf(1.0 / N_TAPER as f64);
    let mut dx = vec![COARSE; NPML + 30];
    for m in 1..=N_TAPER {
        dx.push(COARSE * r.powi(m as i32));
    }
    dx.extend(std::iter::repeat_n(FINE, 40));
    for m in (1..N_TAPER).rev() {
        dx.push(COARSE * r.powi(m as i32));
    }
    dx.extend(std::iter::repeat_n(COARSE, 60 + NPML));
    dx
}

fn taper_spec_drive(nx: usize, dt: f64) -> (FdtdSpec, Boundary, Drive) {
    let mut spec = FdtdSpec::vacuum(nx, TNY, TNZ, COARSE);
    spec.dt = dt;
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML));
    let mut drive = Drive::default();
    drive.soft_sources.push(SoftSource {
        component: EComponent::Ez,
        cell: (SOURCE_I, TNY / 2, TNZ / 2),
        waveform: Waveform::Gaussian {
            t0: 48.0 * dt,
            sigma: 12.0 * dt,
        },
    });
    drive.probes.push(Probe {
        component: EComponent::Ez,
        cell: (PROBE_I, TNY / 2, TNZ / 2),
    });
    (spec, boundary, drive)
}

fn cpu_taper_trace(nx: usize, spacings: &GradedSpacings, dt: f64) -> Vec<f64> {
    let (spec, boundary, drive) = taper_spec_drive(nx, dt);
    let mut engine = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive,
    );
    engine.set_spacings(spacings);
    engine.step_n(T_STEPS);
    engine.probe_series()[0].clone()
}

/// GPU probe trace; `None` when no adapter exists (skip).
fn gpu_taper_trace(nx: usize, spacings: Option<&GradedSpacings>, dt: f64) -> Option<Vec<f64>> {
    let (spec, boundary, drive) = taper_spec_drive(nx, dt);
    let mut engine = match GpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive,
        T_STEPS,
    ) {
        Ok(engine) => engine,
        Err(ComputeError::NoAdapter) => return None,
        Err(other) => panic!("compute-021: GPU backend construction failed: {other}"),
    };
    if let Some(g) = spacings {
        engine.set_spacings(g).expect("taper spacings rejected");
    }
    engine.step_n(T_STEPS).expect("GPU stepping failed");
    Some(engine.read_probes().expect("GPU probe readback failed")[0].clone())
}

#[test]
#[ignore = "slow: three 560-step ~260k-cell runs (one CPU, two GPU); compute-021 graded GPU taper gate, run in release"]
fn gpu_graded_taper_parity() {
    let dx = graded_dx();
    let graded = GradedSpacings {
        dx: dx.clone(),
        dy: vec![COARSE; TNY],
        dz: vec![COARSE; TNZ],
    };
    let dt = 0.9 * graded.courant_limit();
    let len_m: f64 = dx.iter().sum();
    let nx_ref = (len_m / COARSE).round() as usize;

    let Some(gpu_trace) = gpu_taper_trace(dx.len(), Some(&graded), dt) else {
        eprintln!(
            "SKIPPED compute-021: no wgpu adapter on this machine; \
             the GPU nightly runner exercises this gate for real"
        );
        return;
    };

    // --- (a) graded-GPU stability + CPU↔GPU probe-series parity.
    assert!(
        gpu_trace.iter().all(|x| x.is_finite()),
        "compute-021: GPU graded trace went non-finite (instability)"
    );
    let cpu_trace = cpu_taper_trace(dx.len(), &graded, dt);
    assert_eq!(cpu_trace.len(), gpu_trace.len());
    let cpu_peak = linf(&cpu_trace);
    assert!(
        cpu_peak > 0.0,
        "no pulse at the CPU probe (test logic error)"
    );
    let diff: Vec<f64> = cpu_trace
        .iter()
        .zip(&gpu_trace)
        .map(|(c, g)| c - g)
        .collect();
    let rel_linf = linf(&diff) / cpu_peak;
    let rel_l2 = l2(&diff) / l2(&cpu_trace);
    eprintln!("compute-021: CPU peak            = {cpu_peak:.3e}");
    eprintln!("compute-021: series rel L∞       = {rel_linf:.3e}");
    eprintln!("compute-021: series rel L2       = {rel_l2:.3e}");
    assert!(
        rel_linf < SERIES_REL_TOL,
        "compute-021: probe-series rel L∞ {rel_linf:e} ≥ {SERIES_REL_TOL:e}"
    );
    assert!(
        rel_l2 < SERIES_REL_TOL,
        "compute-021: probe-series rel L2 {rel_l2:e} ≥ {SERIES_REL_TOL:e}"
    );

    // --- (b) the compute-019 physics figure, on the GPU: grading reflection
    // via the uniform-reference difference method. Upstream of the grading
    // the two GPU runs execute bit-identical FP32 operations (identical
    // inverse-spacing values, shared dt), so the difference isolates the
    // grading reflection exactly as on the CPU; the early window pins that
    // isolation property. Same −48 dB pinned floor as compute-019.
    let gpu_uniform = gpu_taper_trace(nx_ref, None, dt).expect("adapter vanished mid-test");
    let incident_peak = linf(&gpu_uniform);
    assert!(incident_peak > 0.0, "no incident pulse at the GPU probe");
    let early_max = gpu_trace[..60]
        .iter()
        .zip(&gpu_uniform[..60])
        .map(|(g, u)| (g - u).abs())
        .fold(0.0_f64, f64::max);
    eprintln!("compute-021: early-window Δ      = {early_max:.3e}");
    assert!(
        early_max < 1.0e-9 * incident_peak,
        "compute-021: GPU runs diverged before the pulse reached the grading \
         (Δ = {early_max:e}) — the difference no longer isolates the reflection"
    );
    let reflected_peak = gpu_trace
        .iter()
        .zip(&gpu_uniform)
        .map(|(g, u)| (g - u).abs())
        .fold(0.0_f64, f64::max);
    let reflection_db = 20.0 * (reflected_peak / incident_peak).log10();
    eprintln!("compute-021: incident peak       = {incident_peak:.3e}");
    eprintln!("compute-021: grading reflection  = {reflected_peak:.3e}");
    eprintln!("compute-021: reflection level    = {reflection_db:.2} dB");
    assert!(
        reflection_db <= -48.0,
        "compute-021: GPU grading reflection {reflection_db:.2} dB is above \
         the pinned −48 dB compute-019 floor"
    );
}
