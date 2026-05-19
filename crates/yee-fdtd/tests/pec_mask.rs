//! Per-cell interior PEC-mask reflection / shielding regression test
//! (Phase 2.fdtd.7.z infrastructure).
//!
//! Verifies that a [`yee_fdtd::YeeGrid`] equipped with per-component PEC
//! masks via [`yee_fdtd::YeeGrid::with_pec_mask_ez`] (et al.) behaves
//! as an interior conducting sheet: an incident pulse driven from one
//! side of the sheet is reflected, and the field on the opposite side
//! is suppressed to < 1 % of the incident.
//!
//! ## Geometry choice (why an `x = const` mask plane)
//!
//! For a Gaussian `E_z` point source the source currents are vertical
//! ("z-oriented dipoles"). The corresponding radiated E field at a
//! field point in the equatorial (`z = const`) plane is also dominated
//! by `E_z` plus near-field `H_x` / `H_y`. To make the PEC sheet
//! actually intercept the radiated `E_z` we orient the sheet so
//! `E_z` is *tangential* to it — i.e. the sheet is a plane of constant
//! `x` (or constant `y`), with normal `x̂` (resp. `ŷ`). On such a
//! sheet the tangential E components are `E_y` and `E_z`; clamping
//! both to zero models a thin sheet of perfect conductor and reflects
//! the incident `E_z` cleanly.
//!
//! ## Setup
//!
//! - Grid: `60 × 60 × 60` cells, `dx = 1 mm`. Deprecated outer-face
//!   PEC on all six faces. The `mask − vac` difference cancels any
//!   spurious wall reflection that appears identically in both runs.
//! - PEC mask: the full `i = i_mask = 35` cross-section, i.e. the
//!   entire `(j, k) ∈ [0, N) × [0, N)` rectangle. The mask sets the
//!   `pec_mask_ey` and `pec_mask_ez` arrays to `true` across that
//!   plane. The full-cross-section sheet leaves the wave nowhere to
//!   diffract around — the grid is bisected into a source half-space
//!   (`i < I_MASK`) and a shielded half-space (`i > I_MASK`) and the
//!   only signal path between them is through the PEC sheet itself.
//! - Source: soft Gaussian-in-time pulse on `E_z` at `(25, 30, 30)`
//!   — 10 cells "behind" the sheet (lower-x side).
//! - Probe A (reflection): `E_z` at `(30, 30, 30)`, between source
//!   and sheet.
//! - Probe B (shielded): `E_z` at `(42, 30, 30)`, 7 cells beyond the
//!   sheet on the far side. The field at probe B in the masked run
//!   must be < 1 % of the incident peak (and is in practice exactly
//!   zero — the full-cross-section mask blocks all signal paths).

use ndarray::Array3;

use yee_fdtd::{WalkingSkeletonSolver, YeeGrid};

fn run_trace(
    mut solver: WalkingSkeletonSolver,
    n_steps: usize,
    source: (usize, usize, usize),
    probes: &[(usize, usize, usize)],
    t0: f64,
    sigma: f64,
) -> Vec<Vec<f64>> {
    let mut traces: Vec<Vec<f64>> = probes.iter().map(|_| Vec::with_capacity(n_steps)).collect();
    for _ in 0..n_steps {
        solver.step_with_source(source.0, source.1, source.2, t0, sigma);
        for (t, &p) in traces.iter_mut().zip(probes.iter()) {
            t.push(solver.grid().ez[p]);
        }
    }
    traces
}

#[test]
fn pec_mask_reflects_and_shields() {
    const N: usize = 60;
    const DX: f64 = 1.0e-3;
    const N_STEPS: usize = 200;

    // Source at i=25; sheet at i=35; probe-reflect at i=30 between
    // them; probe-shielded at i=42, well past the sheet.
    const SOURCE: (usize, usize, usize) = (25, 30, 30);
    const PROBE_REFLECT: (usize, usize, usize) = (30, 30, 30);
    const PROBE_SHIELDED: (usize, usize, usize) = (42, 30, 30);
    const I_MASK: usize = 35;
    // Sheet extent in the transverse (j, k) plane. We use the full
    // grid in j and k so the wave has nowhere to diffract around
    // the sheet — the outer-face PEC clamps the field along the
    // grid edges, and the interior PEC mask spans the entire i =
    // I_MASK cross-section. Effectively the grid is bisected into a
    // "source side" (i < I_MASK) and a "shielded side" (i > I_MASK)
    // with no signal path between them other than through the PEC
    // sheet itself.
    const MASK_J_LO: usize = 0;
    const MASK_J_HI: usize = N;
    const MASK_K_LO: usize = 0;
    const MASK_K_HI: usize = N;

    let dt = YeeGrid::vacuum(N, N, N, DX).dt;
    let t0 = 20.0 * dt;
    let sigma = 5.0 * dt;

    // ---- Build per-component PEC masks ----
    //
    // On the `i = I_MASK` plane (a constant-x sheet), the tangential
    // E components are `E_y` (shape [N+1, N, N+1]) and `E_z`
    // (shape [N+1, N+1, N]). The corresponding mask shapes match the
    // E component arrays exactly; we clamp both to true on the
    // rectangle `(j ∈ [MASK_J_LO, MASK_J_HI)) × (k ∈ [MASK_K_LO,
    // MASK_K_HI))` at i = I_MASK.
    let mut mask_ey = Array3::<bool>::from_elem((N + 1, N, N + 1), false);
    for j in MASK_J_LO..MASK_J_HI {
        for k in MASK_K_LO..(MASK_K_HI + 1) {
            mask_ey[(I_MASK, j, k)] = true;
        }
    }
    let mut mask_ez = Array3::<bool>::from_elem((N + 1, N + 1, N), false);
    for j in MASK_J_LO..(MASK_J_HI + 1) {
        for k in MASK_K_LO..MASK_K_HI {
            mask_ez[(I_MASK, j, k)] = true;
        }
    }

    // ---- Masked run ----
    let grid_mask = YeeGrid::vacuum(N, N, N, DX)
        .with_pec_mask_ey(mask_ey)
        .with_pec_mask_ez(mask_ez);
    let traces_mask = run_trace(
        WalkingSkeletonSolver::new(grid_mask),
        N_STEPS,
        SOURCE,
        &[PROBE_REFLECT, PROBE_SHIELDED],
        t0,
        sigma,
    );

    // ---- Vacuum reference (no mask) ----
    let grid_vac = YeeGrid::vacuum(N, N, N, DX);
    let traces_vac = run_trace(
        WalkingSkeletonSolver::new(grid_vac),
        N_STEPS,
        SOURCE,
        &[PROBE_REFLECT, PROBE_SHIELDED],
        t0,
        sigma,
    );

    // ---- Sanity ----
    for (label, traces) in [("mask", &traces_mask), ("vac", &traces_vac)] {
        for (i, tr) in traces.iter().enumerate() {
            assert!(
                tr.iter().all(|x| x.is_finite()),
                "{label} probe {i} trace went non-finite"
            );
        }
    }

    let inc_peak = traces_vac[0]
        .iter()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);
    assert!(inc_peak > 0.0, "no incident pulse observed at probe A");

    // ---- Reflection amplitude at probe A ----
    //
    // `mask − vac` at probe A: the incident pulse cancels (same in
    // both runs until the reflection arrives), leaving the reflected
    // pulse. For a planar PEC at normal incidence Γ = −1, so the
    // reflected pulse magnitude equals the incident magnitude. With
    // a near-field point source the measured magnitude is reduced
    // by 1/r spreading (image source farther from probe than direct
    // source) and by the Ez dipole's anisotropic radiation pattern
    // (most energy radiates broadside, not toward the sheet). We
    // test that the reflected amplitude is at least 10 % of the
    // incident peak — well above any "no reflection" noise floor
    // (~10⁻³ × incident) and easily satisfied by the analytical
    // Γ = −1 plus geometric attenuation (~0.20 in this setup).
    let diff_a: Vec<f64> = traces_mask[0]
        .iter()
        .zip(traces_vac[0].iter())
        .map(|(m, v)| m - v)
        .collect();
    let refl_peak = diff_a.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);

    // ---- Transmission amplitude at probe B (shielded side) ----
    //
    // The DoD is "transmitted field amplitude below mask < 1 % of
    // incident". We compute the peak |E_z| at probe B in the masked
    // run and compare it to the *incident peak at probe A* (the
    // strongest signal the source produces). A correct PEC mask
    // reduces the transmitted-side field to numerical noise.
    let trans_peak_mask = traces_mask[1]
        .iter()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);
    let trans_peak_vac = traces_vac[1]
        .iter()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);

    eprintln!("incident peak at probe A      = {inc_peak:.3e}");
    eprintln!("reflected peak |mask-vac| (A) = {refl_peak:.3e}");
    eprintln!(
        "reflection / incident         = {:.3}",
        refl_peak / inc_peak
    );
    eprintln!("|E_z| at probe B (mask)       = {trans_peak_mask:.3e}");
    eprintln!("|E_z| at probe B (vac)        = {trans_peak_vac:.3e}");
    eprintln!(
        "transmission / vac-amplitude  = {:.4}",
        trans_peak_mask / trans_peak_vac.max(1e-30)
    );
    eprintln!(
        "transmission / incident       = {:.4}",
        trans_peak_mask / inc_peak
    );

    // ---- Assertion thresholds ----
    //
    // Reflection: at normal incidence on a planar PEC, Γ = −1 (full
    // reflection). Because the source is a Gaussian E_z point pulse
    // and the probe sits in the near field, geometric 1/r spreading
    // reduces the measured "reflection / incident" ratio below
    // unity. With the chosen geometry (source 10 cells from sheet,
    // probe 5 cells from source on the source side), the image-source
    // distance to the probe is 15 cells vs the direct 5 cells, so the
    // ideal geometric reflection ratio is 5/15 ≈ 0.33 times Γ. We
    // measure ~0.20 (some additional reduction from the source's
    // anisotropic Ez-dipole radiation pattern, which puts most energy
    // broadside rather than along the symmetry axis toward the
    // sheet). The 0.10 threshold is a conservative
    // "is-clearly-reflecting" gate.
    //
    // Transmission: the brief's < 1 % DoD; in practice with a
    // full-cross-section mask the transmission is exactly 0 to
    // machine precision because the masked plane is the only signal
    // path between the source and shielded half-spaces.
    assert!(
        refl_peak / inc_peak > 0.10,
        "PEC mask reflection {:.3} (relative to incident peak) is \
         below 0.10 — the mask is not acting as a reflector",
        refl_peak / inc_peak
    );
    assert!(
        trans_peak_mask / inc_peak < 0.01,
        "PEC mask transmission {:.4} (peak |E_z| beyond the mask, \
         relative to incident peak) exceeds the 1 % DoD; mask may be \
         letting field leak through",
        trans_peak_mask / inc_peak
    );
}
