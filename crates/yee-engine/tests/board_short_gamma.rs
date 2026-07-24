//! Gate `engine-sparams-003` (R.2b, ADR-0200): **measured complex Γ of a
//! known one-port against transmission-line theory** — the calibration
//! step the ADR-0195 negative result called for. A via-shorted microstrip
//! line is the cleanest one-port there is: at the via plane Γ = −1
//! exactly, and at a reference plane a distance `d` before it,
//!
//! `Γ(f) = −e^{−2jβ(f)d}`,   `β = 2π f √ε_eff / c`  (Hammerstad–Jensen)
//!
//! — unit magnitude, and an unwrapped phase whose slope is the **round
//! trip** `dφ/df = −4π d √ε_eff / c` (twice the R.2 through-line slope).
//! The gate is **differential over two fence distances** (two release
//! solves): the phase-slope DIFFERENCE between d₂ and d₁ depends only on
//! the extra line length — the termination's own reactance cancels
//! exactly. The instrumented single-run iterations that forced each
//! design choice:
//!
//! - a single-centre-cell via read |Γ| ≈ 0.98 but a 3× slope — it is a
//!   partial shunt on a 10-cell-wide line and the passed wave reflected
//!   off the open line end (composite reflector) → the short became a
//!   full-width **via fence**;
//! - a fence with line continuing beyond still read ~1.8× — the fence
//!   leaks → the trace now ENDS at the fence plane (one reflection
//!   plane);
//! - the single-plane run read +14 % slope excess ≈ the fence's own
//!   inductive reactance (an inductive short's dφ/df adds apparent
//!   depth; ~0.3 nH at these dimensions) → the assert went differential,
//!   which cancels it.
//!
//! Bins are selected by the fit's own quality flags (ADR-0189 residual +
//! fitted β vs HJ): with a TOTAL reflector the standing-wave nulls sweep
//! across the probes and the 3-probe split degenerates where a null sits
//! near the middle probe (measured |Γ| → 0 fallbacks there).
//!
//! Unlike the ADR-0195 THRU case (where plane-A "Γ" was the far port's
//! residual reflection — a small, fixture-dominated quantity), the short
//! is a total reflector: the backward wave is as large as the forward
//! one, so the directional split works with full signal on both arms.
//!
//! `#[ignore]`'d (two multi-minute release FDTD runs):
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_short_gamma -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, voxelize_microstrip, with_via_at_cell};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const BW_HZ: f64 = 4.0e9;
const N_STEPS: usize = 9000;
const SPACING_CELLS: usize = 17;
/// The two reference-plane-to-short distances (metres) of the
/// differential measurement; ~0.73 / 1.46 λ_g at 5 GHz, so the phase
/// wraps and the unwrap is exercised.
const D1_M: f64 = 12.0e-3;
const D2_M: f64 = 24.0e-3;

/// One shorted-line run: returns (mean |Γ| over quality bins, fitted
/// dφ/df, cell-snapped d).
fn measure_short(d_m: f64) -> (f64, f64, f64) {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    // Drive at x≈0; probe triple starts one λ_g in (past the launch
    // near-field); the line ENDS at the via-fence plane D_M beyond the
    // first probe — fence and trace end coincide, so there is exactly ONE
    // reflection plane. (Earlier iterations left 5 mm of open line beyond
    // the fence; the fence leaks enough that the composite two-plane
    // reflector read a ~1.8× phase slope.)
    let x_a0 = lam_g;
    let l_m = x_a0 + d_m;
    let traces = vec![Polygon::rect(0.0, 0.0, l_m, W_M)];
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: Substrate {
            eps_r: EPS_R,
            height_m: H_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
            PortRef {
                at: Point2::new(l_m - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
        ],
        bbox,
    };

    let mut model = voxelize_microstrip(
        &layout,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_a0 = i_for(x_a0);
    // The last trace column: the fence lands exactly at the line's end.
    let i_via = i_for(l_m - dx / 2.0);
    // The distance the phase gate actually sees: cell-snapped.
    let d_snapped = (i_via - i_a0) as f64 * dx;

    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < W_M / 2.0 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);

    // The short: a via FENCE across the full trace width — a PEC wall
    // from ground to trace, so the whole quasi-TEM field is terminated.
    for j in j_lo..j_hi {
        with_via_at_cell(&mut model, i_via, j, k_top);
    }

    let materials = MaterialsSpec {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ex: model
            .grid
            .pec_mask_ex
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ey: model
            .grid
            .pec_mask_ey
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ez: model
            .grid
            .pec_mask_ez
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        ..MaterialsSpec::default()
    };

    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * BW_HZ)) / dt).ceil() as usize;

    let mk_probe = |i: usize| ProbeSpec {
        component: "ez".into(),
        cell: (i, j_strip, k_probe),
    };
    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: DX_M,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        // Drive port only: the via IS the termination under test.
        aperture_ports: vec![AperturePortSpec {
            i: model.port_cells[0].0,
            j_lo,
            j_hi,
            k_lo: 0,
            k_top,
            resistance_ohm: Z0_OHM,
            v0: 1.0,
            f0_hz: F0_HZ,
            bw_hz: BW_HZ,
            t0_steps,
            record: false,
        }],
        thin_wires: vec![],
        probes: vec![
            mk_probe(i_a0),
            mk_probe(i_a0 + SPACING_CELLS),
            mk_probe(i_a0 + 2 * SPACING_CELLS),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        spacings: None,
        backend: BackendChoice::Cpu,
    };

    let handle = yee_engine::submit(spec);
    let result = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");
    let p = &result.probes;
    let spacing_m = SPACING_CELLS as f64 * DX_M;

    // 4.0–6.0 GHz, 50 MHz raster. With a TOTAL reflector the line carries
    // a deep standing wave whose nulls sweep across the probe triple as
    // frequency sweeps; bins where a null sits near the middle probe make
    // the 3-probe fit ill-conditioned (the first fence run measured
    // |Γ| → 0 fallbacks and a garbage average slope there). Select bins by
    // the fit's own quality flags — the ADR-0189 residual |Im cos βd| and
    // the fitted β against Hammerstad–Jensen — and measure on those.
    let freqs: Vec<f64> = (0..=40).map(|n| 4.0e9 + n as f64 * 50.0e6).collect();
    let mut good_freqs: Vec<f64> = vec![];
    let mut gamma: Vec<(f64, f64)> = vec![];
    for &f in &freqs {
        let v: Vec<(f64, f64)> = (0..3)
            .map(|m| sparams::single_bin_dft(&p[m], dt, f))
            .collect();
        let split = sparams::fit_standing_wave(v[0], v[1], v[2], spacing_m);
        let beta_hj = 2.0 * PI * f * e_eff.sqrt() / C0_M_S;
        let beta_ok = (split.beta_rad_m - beta_hj).abs() / beta_hj < 0.10;
        let fwd_mag = split.fwd.0.hypot(split.fwd.1);
        if split.residual < 0.15 && beta_ok && fwd_mag > 0.0 {
            let n = fwd_mag * fwd_mag;
            let g = (
                (split.bwd.0 * split.fwd.0 + split.bwd.1 * split.fwd.1) / n,
                (split.bwd.1 * split.fwd.0 - split.bwd.0 * split.fwd.1) / n,
            );
            good_freqs.push(f);
            gamma.push(g);
        }
    }
    assert!(
        good_freqs.len() >= 20,
        "engine-sparams-003 FAILED: only {} / {} bins pass the fit-quality          selectors",
        good_freqs.len(),
        freqs.len()
    );
    let freqs = good_freqs;

    // ---- 1. Magnitude: a short reflects everything ----
    let mags: Vec<f64> = gamma.iter().map(|g| g.0.hypot(g.1)).collect();
    let mag_mean = mags.iter().sum::<f64>() / mags.len() as f64;

    // ---- 2. Unwrapped phase slope: the round trip to the short ----
    let mut phase: Vec<f64> = gamma.iter().map(|g| g.1.atan2(g.0)).collect();
    for n in 1..phase.len() {
        while phase[n] - phase[n - 1] > PI {
            phase[n] -= 2.0 * PI;
        }
        while phase[n] - phase[n - 1] < -PI {
            phase[n] += 2.0 * PI;
        }
    }
    let nf = freqs.len() as f64;
    let fm = freqs.iter().sum::<f64>() / nf;
    let pm = phase.iter().sum::<f64>() / nf;
    let slope = freqs
        .iter()
        .zip(&phase)
        .map(|(f, ph)| (f - fm) * (ph - pm))
        .sum::<f64>()
        / freqs.iter().map(|f| (f - fm).powi(2)).sum::<f64>();
    (mag_mean, slope, d_snapped)
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs; engine-sparams-003 gate (R.2b) — run with --release --ignored"]
fn shorted_line_gamma_matches_transmission_line_theory() {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let (mag1, slope1, d1) = measure_short(D1_M);
    let (mag2, slope2, d2) = measure_short(D2_M);

    // Differential: the termination reactance cancels; only the extra
    // line length remains.
    let dslope = slope2 - slope1;
    let dslope_ref = -4.0 * PI * (d2 - d1) * e_eff.sqrt() / C0_M_S;
    let dslope_err = (dslope - dslope_ref).abs() / dslope_ref.abs();

    // Single-run sanity: each slope must sit between the ideal-short value
    // and a bounded inductive excess (measured +14 % at d1).
    let ratio1 = slope1 / (-4.0 * PI * d1 * e_eff.sqrt() / C0_M_S);
    let ratio2 = slope2 / (-4.0 * PI * d2 * e_eff.sqrt() / C0_M_S);

    eprintln!(
        "engine-sparams-003: |Γ| mean {mag1:.3} / {mag2:.3} | single-run slope \
         ratios {ratio1:.3} / {ratio2:.3} (ideal 1, inductive excess) | \
         differential dφ/df = {dslope:.4e} vs TL {dslope_ref:.4e} → err {:.2} % \
         (Δd = {:.2} mm)",
        dslope_err * 100.0,
        (d2 - d1) * 1e3,
    );

    assert!(
        (mag1 - 1.0).abs() <= 0.15 && (mag2 - 1.0).abs() <= 0.15,
        "engine-sparams-003 FAILED: short must reflect near-unity \
         (|Γ| means {mag1:.3} / {mag2:.3}, ±15 %)"
    );
    assert!(
        dslope_err <= 0.05,
        "engine-sparams-003 FAILED: differential round-trip phase slope err \
         {:.2} % (> 5 %)",
        dslope_err * 100.0
    );
    assert!(
        (0.95..=1.35).contains(&ratio1) && (0.95..=1.35).contains(&ratio2),
        "engine-sparams-003 FAILED: single-run slope ratios {ratio1:.3} / \
         {ratio2:.3} outside [0.95, 1.35] (ideal short + bounded inductive \
         excess)"
    );
}
