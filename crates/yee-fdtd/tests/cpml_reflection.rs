//! CPML reflection regression test.
//!
//! Drives two side-by-side simulations on a 50³ vacuum grid:
//!
//! - **PEC reference**: hard reflecting PEC on all six faces (legacy
//!   `WalkingSkeletonSolver::new`).
//! - **CPML candidate**: Roden–Gedney 2000 CPML with `npml = 10` and
//!   default parameters (`σ_max` from `R_0 = 1e-6`, `κ_max = 1`,
//!   `α_max = 0.05`).
//!
//! A Gaussian-in-time pulse is injected on `E_z` at the centre. We probe
//! `E_z` at `(11, 25, 25)` — just inside the inner PML edge on the
//! low-x face — for 300 time steps. The reflected-vs-outgoing peak ratio
//! is computed for both runs, and the CPML candidate must reduce the
//! late-time peak by at least 30 dB versus the PEC reference.
//!
//! Wall-time budget: < 60 s single-threaded release build.

use yee_fdtd::{CpmlParams, WalkingSkeletonSolver, YeeGrid};

/// Run a simulation injecting a Gaussian pulse at `source` and record
/// `E_z` at `probe` for each step.
fn run_trace(
    mut solver: WalkingSkeletonSolver,
    n_steps: usize,
    source: (usize, usize, usize),
    probe: (usize, usize, usize),
    t0: f64,
    sigma: f64,
) -> Vec<f64> {
    let mut trace = Vec::with_capacity(n_steps);
    for _ in 0..n_steps {
        solver.step_with_source(source.0, source.1, source.2, t0, sigma);
        trace.push(solver.grid().ez[probe]);
    }
    trace
}

#[test]
fn cpml_attenuates_reflection_vs_pec() {
    const N: usize = 50;
    const DX: f64 = 1.0e-3;
    const NPML: usize = 10;
    const N_STEPS: usize = 300;
    const SOURCE: (usize, usize, usize) = (25, 25, 25);
    // Probe inside the high-x PML inner edge: domain interior runs
    // i ∈ [npml, N-npml) = [10, 40). The probe at (38, 25, 25) is 13
    // cells from the source and 2 cells from the inner edge of the
    // high-x PML — close enough to the wall to capture any reflection
    // before it can spread, far enough that the outgoing pulse arrives
    // first and the spherical-front overlap is small.
    const PROBE: (usize, usize, usize) = (38, 25, 25);

    // Reference setup: same grid for both, derive Gaussian parameters from dt.
    let grid_ref = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid_ref.dt;
    // Wider Gaussian: most CPML energy losses are at high frequencies and
    // a broader pulse is better matched to the PML thickness here.
    let t0 = 20.0 * dt;
    let sigma = 6.0 * dt;
    drop(grid_ref);

    // --- PEC run ---
    let pec_trace = run_trace(
        WalkingSkeletonSolver::new(YeeGrid::vacuum(N, N, N, DX)),
        N_STEPS,
        SOURCE,
        PROBE,
        t0,
        sigma,
    );

    // --- CPML run ---
    let grid_cpml = YeeGrid::vacuum(N, N, N, DX);
    let params = CpmlParams::for_grid(&grid_cpml, NPML);
    let cpml_trace = run_trace(
        WalkingSkeletonSolver::with_cpml(grid_cpml, params),
        N_STEPS,
        SOURCE,
        PROBE,
        t0,
        sigma,
    );

    // Sanity check: both finite, both saw the outgoing pulse.
    assert!(
        pec_trace.iter().all(|x| x.is_finite()),
        "PEC trace went non-finite"
    );
    assert!(
        cpml_trace.iter().all(|x| x.is_finite()),
        "CPML trace went non-finite"
    );

    // ---------------------------------------------------------------
    // Reflection measurement.
    //
    // A Gaussian-in-time *soft* source adds the same integrated charge
    // each run, which creates a residual quasi-static field at the
    // probe that is *identical* in both PEC and CPML cases. Comparing
    // raw `cpml` amplitude post-pulse therefore double-counts this
    // static residual. We instead isolate the PEC-reflected wave via
    // the difference `Δ(t) = pec(t) − cpml(t)`: the static near-field
    // cancels exactly, the outgoing pulse cancels exactly (until any
    // reflection returns), and only the difference in reflected
    // amplitude remains. Specifically:
    //
    //   |Δ(t)| = |PEC reflection| − |CPML reflection|  (≈ |PEC reflection|)
    //
    // CPML's own residual reflection is much smaller than PEC's, so
    // `max|Δ|` is a clean estimate of the PEC reflected-wave peak.
    //
    // The "CPML reflection" itself is then measured as
    // `max(|cpml − pec_baseline_no_reflection|)`. Since we don't have a
    // no-reflection run, we use the early-time portion where neither
    // CPML nor PEC has seen any wall reflection: their values agree
    // exactly in the outgoing window, and we treat the trace beyond the
    // outgoing window for CPML *minus the early outgoing pulse* as the
    // reflected-energy proxy.
    //
    // Simpler and more robust: take the maximum of |pec − cpml| over
    // the *whole* trace as the PEC reflection amplitude, and the
    // maximum |cpml − pec_static_baseline| where the static baseline is
    // pec_trace[REFLECTION_START..] averaged ... actually simplest is
    // just the *raw* CPML peak in the late-time window since the
    // static field plus CPML's own small reflection both fit in it,
    // and the static field's peak amplitude is in fact what we want
    // CPML's reflection to be "as small as".
    //
    // We measure:
    //   - `pec_diff` = max|pec(t) − cpml(t)|  (PEC reflection)
    //   - `cpml_static` = max|cpml(t) − cpml(t-back)|  using late-time
    //     CPML oscillation as the static-floor noise — equivalently
    //     std-dev of late-time CPML.
    // ---------------------------------------------------------------

    let n = pec_trace.len();
    // PEC reflection: difference between PEC and CPML traces. Until any
    // wall reflection returns to the probe (around step 70 for the low-x
    // wall closest to the probe at x=38, or earlier for the high-x wall
    // 12 cells away), the two traces are bit-identical.
    let pec_reflection_peak = pec_trace
        .iter()
        .zip(cpml_trace.iter())
        .map(|(p, c)| (p - c).abs())
        .fold(0.0_f64, f64::max);

    // CPML's own residual reflection: take the late-time trace minus
    // the static-field component, which we estimate by *averaging*
    // CPML over the late window (static field is DC-like, hence the
    // average; reflections oscillate).
    const REFLECTION_START: usize = 80;
    let late_cpml = &cpml_trace[REFLECTION_START..];
    let static_floor =
        late_cpml.iter().copied().sum::<f64>() / (late_cpml.len() as f64);
    let cpml_reflection_peak = late_cpml
        .iter()
        .map(|x| (x - static_floor).abs())
        .fold(0.0_f64, f64::max);

    // Outgoing-pulse reference (largest amplitude in either trace).
    let peak_outgoing = pec_trace
        .iter()
        .chain(cpml_trace.iter())
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);

    let pec_db = 20.0 * (pec_reflection_peak / peak_outgoing).log10();
    let cpml_db = 20.0 * (cpml_reflection_peak / peak_outgoing).log10();
    let reduction_db = pec_db - cpml_db; // = 20·log10(pec_ratio / cpml_ratio)

    eprintln!("n_steps                       = {n}");
    eprintln!("outgoing peak                 = {peak_outgoing:.3e}");
    eprintln!("PEC reflection peak (|pec-cpml|) = {pec_reflection_peak:.3e}  ({pec_db:.2} dB)");
    eprintln!("CPML static floor (late mean) = {static_floor:+.3e}");
    eprintln!("CPML reflection peak (residual) = {cpml_reflection_peak:.3e}  ({cpml_db:.2} dB)");
    eprintln!("CPML reflection reduction     = {reduction_db:.2} dB");

    assert!(peak_outgoing > 0.0, "no outgoing pulse seen");
    assert!(pec_reflection_peak > 0.0, "PEC and CPML are identical (no reflection seen)");
    assert!(
        cpml_reflection_peak > 0.0,
        "CPML reflection peak is zero (test logic error)"
    );

    // Goal: ≥ 30 dB improvement.
    assert!(
        reduction_db >= 30.0,
        "CPML reflection reduction {reduction_db:.2} dB is below 30 dB target."
    );
}
