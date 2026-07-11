//! Gate `compute-019` (FS.0b.0, ADR-0208): the numerical reflection caused
//! by a smooth graded-mesh transition must sit below a measured-then-pinned
//! floor.
//!
//! Setup: a pulse propagates along a graded x-axis — 0.5 mm cells, a
//! geometric taper (ratio 1/r ≈ 1.122 per cell ≤ 1.3) down to 0.25 mm over
//! 6 cells, a 40-cell fine region, a taper back — in free space with CPML on
//! every face (uniform 0.5 mm spacing inside the x absorbers, per the
//! FS.0b.0 scope rule). The reference run is a uniform 0.5 mm grid with
//! **identical upstream geometry and the same dt** (the graded Courant
//! step): the two runs execute bit-identical FP operations at the upstream
//! probe until the pulse reaches the grading, so the probe **difference**
//! isolates the grading-caused reflection exactly (source residual,
//! backward pulse, and its near-CPML return all cancel). The measurement
//! window closes before the far-wall CPML returns (which differ between the
//! runs) can reach the probe.
//!
//! Measured 2026-07-11 (release, this gate's `--nocapture` output):
//! incident peak 2.721e-3, grading reflection 6.322e-6 → **−52.68 dB** for
//! the ratio-1.122 taper — comfortably under the −40 dB expectation for
//! ratio ≤ 1.3 grading. Pinned at −48 dB (~4.7 dB margin); see ADR-0208.

use yee_compute::{
    Boundary, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, GradedSpacings, Materials,
    Probe, SoftSource, Waveform,
};

/// Pinned gate floor (dB): measured −52.68 dB, pinned with ~4.7 dB margin
/// (ADR-0208). Do not weaken without re-measuring.
const REFLECTION_FLOOR_DB: f64 = -48.0;

const COARSE: f64 = 0.5e-3;
const FINE: f64 = 0.25e-3;
const NPML: usize = 10;
const N_TAPER: usize = 6; // coarse → fine cells; ratio 2^(1/6) ≈ 1.122
const NY: usize = 40;
const NZ: usize = 40;
const N_STEPS: usize = 560;

const SOURCE_I: usize = 16;
const PROBE_I: usize = 34;

/// Primal x widths: [10 CPML + 30 run-up]·0.5 mm, 6-cell geometric taper to
/// 0.25 mm, 40 fine cells, 5-cell taper back, [60 + 10 CPML]·0.5 mm (the
/// long tail keeps the far-wall CPML return — which does NOT cancel between
/// the runs — from reaching the probe inside the measurement window).
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

fn probe_trace(nx: usize, spacings: Option<GradedSpacings>, dt: f64) -> Vec<f64> {
    let mut spec = FdtdSpec::vacuum(nx, NY, NZ, COARSE);
    spec.dt = dt;
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML));
    let mut drive = Drive::default();
    drive.soft_sources.push(SoftSource {
        component: EComponent::Ez,
        cell: (SOURCE_I, NY / 2, NZ / 2),
        waveform: Waveform::Gaussian {
            t0: 48.0 * dt,
            sigma: 12.0 * dt,
        },
    });
    drive.probes.push(Probe {
        component: EComponent::Ez,
        cell: (PROBE_I, NY / 2, NZ / 2),
    });
    let mut engine = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive,
    );
    if let Some(g) = spacings {
        engine.set_spacings(&g);
    }
    engine.step_n(N_STEPS);
    engine.probe_series()[0].clone()
}

#[test]
#[ignore = "slow: two 560-step ~260k-cell runs; compute-019 graded reflection gate, run in release"]
fn graded_transition_reflection_is_below_floor() {
    let dx = graded_dx();
    let graded = GradedSpacings {
        dx: dx.clone(),
        dy: vec![COARSE; NY],
        dz: vec![COARSE; NZ],
    };
    // Both runs share the graded Courant dt so the upstream FP evolution is
    // bit-identical until the pulse reaches the grading.
    let dt = 0.9 * graded.courant_limit();

    // Reference: uniform coarse grid of (approximately) the same physical
    // length — only the upstream identity and far-return timing matter.
    let len_m: f64 = dx.iter().sum();
    let nx_ref = (len_m / COARSE).round() as usize;

    let graded_trace = probe_trace(dx.len(), Some(graded), dt);
    let uniform_trace = probe_trace(nx_ref, None, dt);

    assert!(
        graded_trace.iter().all(|x| x.is_finite()),
        "graded trace went non-finite (instability)"
    );
    let incident_peak = uniform_trace.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    assert!(incident_peak > 0.0, "no incident pulse at the probe");

    // The physical wavefront cannot complete source → grading → probe
    // (15 mm ≈ 82 steps) before step ~82; the early window must therefore
    // hold only the (super-exponentially small) outside-light-cone
    // precursor. This pins the isolation property the difference method
    // relies on: upstream, the two runs are the same computation.
    let early_max = graded_trace[..60]
        .iter()
        .zip(&uniform_trace[..60])
        .map(|(g, u)| (g - u).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        early_max < 1.0e-12 * incident_peak,
        "runs diverged before the pulse reached the grading (Δ = {early_max:e}) — \
         the difference no longer isolates the grading reflection"
    );

    let reflected_peak = graded_trace
        .iter()
        .zip(&uniform_trace)
        .map(|(g, u)| (g - u).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        reflected_peak > 0.0,
        "no grading reflection seen at all (test logic error)"
    );
    let reflection_db = 20.0 * (reflected_peak / incident_peak).log10();

    eprintln!("compute-019: incident peak       = {incident_peak:.3e}");
    eprintln!("compute-019: grading reflection  = {reflected_peak:.3e}");
    eprintln!("compute-019: reflection level    = {reflection_db:.2} dB");

    assert!(
        reflection_db <= REFLECTION_FLOOR_DB,
        "compute-019: grading reflection {reflection_db:.2} dB is above the \
         pinned {REFLECTION_FLOOR_DB} dB floor"
    );
}
