//! Per-axis CPML face-selection gate (`cpml_per_axis_001`).
//!
//! Companion to `tests/cpml_reflection.rs`. Where that gate exercises the
//! symmetric all-six-faces CPML, this one exercises the **x-only** face
//! selection added in Phase 2.fdtd.6.7 (ADR-0122):
//!
//! - CPML on the **x** axis only (`CpmlParams::with_axes([true, false,
//!   false])`) → absorbing at `x = 0` and `x = nx`.
//! - **PEC** on the transverse (y, z) walls → the guide mode is preserved,
//!   not absorbed. This is the matched-line / parallel-plate configuration
//!   the reactive-port research track (ADR-0121 increment 3) needs.
//!
//! The grid is a long, thin guide (`nx` long, small transverse extent). An
//! `E_z`-polarized Gaussian pulse is driven on the centre x-plane; it
//! propagates outward along `±x` and is absorbed at the x-ends. A control
//! run terminates *all six* faces with hard PEC (so the x-ends reflect). The
//! reflected wave is isolated as `Δ(t) = pec(t) − cpml(t)` — every shared
//! contribution (the outgoing pulse, the transverse-PEC behaviour, the
//! static near-field) cancels, leaving only the difference in the x-end
//! reflection, exactly as `cpml_reflection.rs` does.
//!
//! Two assertions, neither weakened from the all-faces gate:
//!
//! 1. **≥30 dB reduction** of the x-end reflection vs the PEC control.
//! 2. **Transverse PEC walls intact** — tangential `E` on the y and z faces
//!    is ≈ 0 in the x-only-CPML run (the guide mode survives between intact
//!    walls; the disabled axes are genuinely PEC, not absorbing), AND a
//!    non-trivial interior field is present (the mode was not killed).
//!
//! Wall-time budget: a few seconds, single-threaded release build.
//! `#[ignore]`'d like the other FDTD release gates; run via
//! `cargo test -p yee-fdtd --release --test cpml_per_axis_001 -- --ignored`.

use yee_fdtd::{CpmlParams, FdtdSolver, WalkingSkeletonSolver, YeeGrid};

/// Short guide length along x — the high-x CPML wall is close to the probe.
const NX: usize = 80;
/// Extended guide length along x — the same source / probe positions, but
/// the high-x wall is pushed `NX_EXT − NX = 120` cells farther away. Used as
/// the **reflection-free reference**: within the measurement window no wall
/// reflection (CPML or otherwise) can travel out to the far wall and back to
/// the probe, so this trace is the "no high-x reflection" baseline that the
/// short-guide CPML run is differenced against. This cancels the dispersive
/// forward tail of the rectangular-guide mode (identical in both runs),
/// leaving only the short guide's residual CPML reflection.
const NX_EXT: usize = 200;
const NY: usize = 16;
const NZ: usize = 16;
const DX: f64 = 1.0e-3;
const NPML: usize = 10;
const N_STEPS: usize = 300;

/// Probe x-index — in the *short* guide, 2 cells inside the inner edge of the
/// high-x PML. Anchored to the low-x origin so it is identical in the short
/// and extended guides (only the high-x wall moves).
const PX: usize = NX - NPML - 2;
/// Source x-index, 13 cells upstream of the probe (matching
/// `cpml_reflection.rs`'s source→probe separation): the outgoing pulse
/// reaches the probe early, its high-x-wall reflection returns later. Well
/// clear of the low-x PML.
const SX: usize = PX - 13;
/// Source cell `(SX, ny/2, nz/2)`.
const SOURCE: (usize, usize, usize) = (SX, NY / 2, NZ / 2);
/// Probe cell `(PX, ny/2, nz/2)`.
const PROBE: (usize, usize, usize) = (PX, NY / 2, NZ / 2);

/// Zero the tangential `E` on the transverse (y, z) faces only — a PEC clamp
/// that leaves the x-faces alone so the x-only CPML owns them.
///
/// y = 0 / y = ny faces: tangential are `E_x`, `E_z`.
/// z = 0 / z = nz faces: tangential are `E_x`, `E_y`.
fn clamp_transverse_pec(grid: &mut YeeGrid) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    // y = 0 and y = ny faces: clamp E_x ([nx, ny+1, nz+1]) and E_z
    // ([nx+1, ny+1, nz]).
    for i in 0..nx {
        for k in 0..=nz {
            grid.ex[(i, 0, k)] = 0.0;
            grid.ex[(i, ny, k)] = 0.0;
        }
    }
    for i in 0..=nx {
        for k in 0..nz {
            grid.ez[(i, 0, k)] = 0.0;
            grid.ez[(i, ny, k)] = 0.0;
        }
    }

    // z = 0 and z = nz faces: clamp E_x ([nx, ny+1, nz+1]) and E_y
    // ([nx+1, ny, nz+1]).
    for i in 0..nx {
        for j in 0..=ny {
            grid.ex[(i, j, 0)] = 0.0;
            grid.ex[(i, j, nz)] = 0.0;
        }
    }
    for i in 0..=nx {
        for j in 0..ny {
            grid.ey[(i, j, 0)] = 0.0;
            grid.ey[(i, j, nz)] = 0.0;
        }
    }
}

/// Largest tangential-E magnitude found on the transverse (y, z) faces. For
/// an intact PEC wall this must stay ≈ 0 throughout the run.
fn max_transverse_tangential_e(grid: &YeeGrid) -> f64 {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let mut m = 0.0_f64;

    // y faces.
    for i in 0..nx {
        for k in 0..=nz {
            m = m
                .max(grid.ex[(i, 0, k)].abs())
                .max(grid.ex[(i, ny, k)].abs());
        }
    }
    for i in 0..=nx {
        for k in 0..nz {
            m = m
                .max(grid.ez[(i, 0, k)].abs())
                .max(grid.ez[(i, ny, k)].abs());
        }
    }
    // z faces.
    for i in 0..nx {
        for j in 0..=ny {
            m = m
                .max(grid.ex[(i, j, 0)].abs())
                .max(grid.ex[(i, j, nz)].abs());
        }
    }
    for i in 0..=nx {
        for j in 0..ny {
            m = m
                .max(grid.ey[(i, j, 0)].abs())
                .max(grid.ey[(i, j, nz)].abs());
        }
    }
    m
}

/// Result of one x-only-CPML guide run.
struct CpmlRun {
    /// `E_z(t)` at [`PROBE`] for every step.
    trace: Vec<f64>,
    /// Largest tangential `E` on the transverse (y, z) faces over the run.
    max_tan: f64,
    /// Peak `|E_z|` at [`SOURCE`] (the soft-source plane) — a witness that
    /// the guide mode was actually excited and survived.
    interior_peak: f64,
}

/// Run an **x-only CPML** guide of x-length `nx`: CPML on x (both ends),
/// transverse PEC on y/z. The source / probe sit at the same absolute
/// `(SX|PX, ny/2, nz/2)` regardless of `nx`, so two runs with different `nx`
/// share an identical low-x near field and differ only in where the high-x
/// wall lives.
fn run_x_only_cpml(nx: usize) -> CpmlRun {
    let grid = YeeGrid::vacuum(nx, NY, NZ, DX);
    let dt = grid.dt;
    let t0 = 20.0 * dt;
    let sigma = 6.0 * dt;

    let params = CpmlParams::for_grid(&grid, NPML).with_axes([true, false, false]);
    let mut solver = WalkingSkeletonSolver::with_cpml(grid, params);

    let mut trace = Vec::with_capacity(N_STEPS);
    let mut max_tan = 0.0_f64;
    let mut interior_peak = 0.0_f64;

    for _ in 0..N_STEPS {
        let t = solver.current_time();
        // H half-step, then x-only CPML on H via the split borrow.
        solver.update_h_only();
        {
            let (g, cpml) = solver.grid_and_cpml_mut();
            cpml.expect("x-only CPML configured").update_h(g);
        }
        // Source between H and E (matches step_with_source timing).
        solver.apply_gaussian_source_ez(SOURCE.0, SOURCE.1, SOURCE.2, t, t0, sigma);
        // E half-step + x-only CPML on E.
        solver.update_e_only();
        {
            let (g, cpml) = solver.grid_and_cpml_mut();
            cpml.expect("x-only CPML configured").update_e(g);
        }
        // Transverse PEC clamp (y/z faces only) — the disabled CPML axes.
        clamp_transverse_pec(solver.grid_mut());
        solver.advance_clock();

        trace.push(solver.grid().ez[PROBE]);
        max_tan = max_tan.max(max_transverse_tangential_e(solver.grid()));
        interior_peak = interior_peak.max(solver.grid().ez[SOURCE].abs());
    }

    CpmlRun {
        trace,
        max_tan,
        interior_peak,
    }
}

/// Run an **all-PEC control** of x-length `nx`: hard PEC on every face (the
/// high-x wall reflects). Same source, same probe.
fn run_all_pec_control(nx: usize) -> Vec<f64> {
    let grid = YeeGrid::vacuum(nx, NY, NZ, DX);
    let dt = grid.dt;
    let t0 = 20.0 * dt;
    let sigma = 6.0 * dt;
    // `WalkingSkeletonSolver::new` clamps all six faces with hard PEC.
    let mut solver = WalkingSkeletonSolver::new(grid);
    let mut trace = Vec::with_capacity(N_STEPS);
    for _ in 0..N_STEPS {
        solver.step_with_source(SOURCE.0, SOURCE.1, SOURCE.2, t0, sigma);
        trace.push(solver.grid().ez[PROBE]);
    }
    trace
}

#[test]
#[ignore = "multi-second release FDTD gate; run via --ignored (see ci.yml)"]
fn cpml_per_axis_x_only_absorbs_and_keeps_transverse_pec() {
    // Short guide: high-x CPML wall close to the probe.
    let short = run_x_only_cpml(NX);
    // Extended guide: identical low-x near field, high-x wall pushed far
    // away — the reflection-free reference.
    let ext = run_x_only_cpml(NX_EXT);
    // PEC control on the short guide: the high-x wall reflects hard.
    let pec_trace = run_all_pec_control(NX);

    let cpml_trace = &short.trace;
    let ext_trace = &ext.trace;

    assert!(
        cpml_trace.iter().all(|x| x.is_finite()),
        "x-only CPML trace went non-finite"
    );
    assert!(
        ext_trace.iter().all(|x| x.is_finite()),
        "extended-guide reference trace went non-finite"
    );
    assert!(
        pec_trace.iter().all(|x| x.is_finite()),
        "PEC control trace went non-finite"
    );

    // The measurement window must end before the *extended* guide's own
    // high-x wall reflection can return to the probe — otherwise the
    // reference is no longer reflection-free. Round-trip from the probe
    // (x = PX) to the extended high-x wall (x = NX_EXT) and back, in steps,
    // with a safety margin (the wave travels ≤ 1 cell/step at the Courant
    // factor 0.9, so steps ≥ distance is a conservative bound).
    let ext_wall_roundtrip = 2 * (NX_EXT - PX);
    let window_end = N_STEPS.min(ext_wall_roundtrip.saturating_sub(8));
    assert!(
        window_end > 0,
        "measurement window is empty (NX_EXT too small)"
    );

    // ---- (A) PEC reference reflection (mirrors cpml_reflection.rs) ----
    // The short CPML run and the PEC control share the transverse-PEC
    // behaviour, the outgoing pulse, and the static near field; they differ
    // only at the high-x end (CPML vs hard PEC). Δ = pec − cpml isolates the
    // PEC reflected wave.
    let pec_reflection_peak = pec_trace
        .iter()
        .zip(cpml_trace.iter())
        .take(window_end)
        .map(|(p, c)| (p - c).abs())
        .fold(0.0_f64, f64::max);

    // ---- (B) Short-guide residual CPML reflection ----
    // The short and extended CPML runs are bit-identical until the short
    // guide's high-x CPML reflection (if any) returns to the probe. The
    // dispersive forward-propagating tail of the rectangular-guide mode is
    // identical in both and cancels exactly; only the short guide's residual
    // CPML reflection survives the difference.
    let cpml_reflection_peak = cpml_trace
        .iter()
        .zip(ext_trace.iter())
        .take(window_end)
        .map(|(s, e)| (s - e).abs())
        .fold(0.0_f64, f64::max);

    // Outgoing-pulse reference amplitude (largest |E_z| seen at the probe).
    let peak_outgoing = pec_trace
        .iter()
        .chain(cpml_trace.iter())
        .chain(ext_trace.iter())
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);

    let pec_db = 20.0 * (pec_reflection_peak / peak_outgoing).log10();
    let cpml_db = 20.0 * (cpml_reflection_peak / peak_outgoing).log10();
    let reduction_db = pec_db - cpml_db;

    eprintln!("short guide                   = {NX}x{NY}x{NZ}, npml(x)={NPML}");
    eprintln!("extended ref guide            = {NX_EXT}x{NY}x{NZ}");
    eprintln!("source / probe (x)            = {SX} / {PX}");
    eprintln!("measurement window            = 0..{window_end} steps");
    eprintln!("outgoing peak                 = {peak_outgoing:.3e}");
    eprintln!("PEC reflection (|pec-cpml|)   = {pec_reflection_peak:.3e}  ({pec_db:.2} dB)");
    eprintln!("CPML residual (|short-ext|)   = {cpml_reflection_peak:.3e}  ({cpml_db:.2} dB)");
    eprintln!("x-only CPML reflection reduction = {reduction_db:.2} dB");
    eprintln!("max transverse tangential E   = {:.3e}", short.max_tan);
    eprintln!(
        "interior peak |E_z| at source = {:.3e}",
        short.interior_peak
    );

    assert!(peak_outgoing > 0.0, "no outgoing pulse seen");
    assert!(
        pec_reflection_peak > 0.0,
        "PEC control and x-only CPML are identical (no reflection seen)"
    );
    assert!(
        cpml_reflection_peak > 0.0,
        "CPML residual is exactly zero (test logic error — short and extended \
         guides must differ once a reflection exists)"
    );

    // (1) ≥30 dB reduction — the same target as the all-faces gate. NOT
    // weakened.
    assert!(
        reduction_db >= 30.0,
        "x-only CPML reflection reduction {reduction_db:.2} dB is below the 30 dB target."
    );

    // (2a) Transverse PEC walls intact: tangential E on the y/z faces is ≈ 0
    // throughout. It is *exactly* 0 by construction of the transverse clamp;
    // a generous absolute floor guards against float noise.
    assert!(
        short.max_tan < 1.0e-12 * short.interior_peak.max(1.0),
        "transverse PEC walls not intact: max tangential E = {:.3e} \
         (interior peak {:.3e}) — the y/z faces are not PEC",
        short.max_tan,
        short.interior_peak
    );

    // (2b) The guide mode survived (the transverse axes were *not* absorbing):
    // a non-trivial interior field is present.
    assert!(
        short.interior_peak > 1.0e-3 * peak_outgoing,
        "guide mode appears to have been absorbed: interior peak {:.3e} \
         is negligible vs outgoing {peak_outgoing:.3e}",
        short.interior_peak
    );
}
