//! Gate `engine-bpf-bo-001` (R.4 / F1.2.1, ADR-0197): the **BPF
//! end-to-end loop with EM-in-the-loop Bayesian optimization** — the
//! deferred F1.2.1 core, now closable because every piece exists:
//! synthesis (`synthesize` → coupling matrix), closed-form seeding
//! (`dimension_hairpin_with_fold`: fold-corrected arms, per-section gaps,
//! qe→tap), full-wave measurement over the job protocol (S.9/S.10 stack,
//! S.12 directional |S21|), the validated coupling-matrix reference
//! (`coupling_matrix_s_params`), and `yee_surrogate::bo::minimize`.
//!
//! **Why BO is load-bearing, not garnish** (the R.4a instrumented runs,
//! probe-dump forensics): the closed-form hairpin seed measures
//! |Γ_in| ≈ 1.0 across the band — the end resonator sits ~+17 % high
//! (corner/open-end effects the midline model can't see) and the
//! effective tap coupling lands far below the designed qe. The seed's
//! S21 peaks near −19 dB: a detuned, under-coupled response that no
//! closed form on this stack repairs. That measured seed→spec residual
//! is exactly what F1.2.1 scheduled surrogate refinement for.
//!
//! The loop: one straight-line reference solve (launch normalization,
//! shared grid via a fixed envelope bbox), then BO over three knobs —
//! `arm_scale` (retunes the resonators), `tap_scale` (retunes the
//! external Q), `gap_scale` (retunes the inter-resonator coupling) —
//! each objective call one full-wave DUT solve, minimizing the RMS
//! misfit (dB, floored) between measured directional |S21| and the
//! coupling-matrix response over 3.5–6.5 GHz. Gate: BO must improve the
//! misfit over the seed and land a real passband (assert numbers set
//! from the first honest converged run, the S.8 pattern).
//!
//! `#[ignore]`'d (~13 multi-minute release FDTD solves, ~1 h):
//!
//! ```bash
//! cargo test -p yee-filter --release --test engine_bpf_bo -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use nalgebra::DVector;
use yee_filter::{
    Approximation, FilterSpec, HairpinDimensions, Response, SpecMask, coupling_matrix_s_params,
    dimension_hairpin_with_fold, synthesize,
};
use yee_layout::{BBox, HairpinSectionParams, Layout, Polygon, Substrate, hairpin_bpf_sections};
use yee_surrogate::bo::{BoConfig, minimize};

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_voxel::{MicrostripModel, VoxelOptions, voxelize_microstrip};

const ORDER: usize = 3;
const F0_HZ: f64 = 5.0e9;
// Stack + FBW sized so the qe→tap lands on the (fold-corrected) arm with
// the feed half-width clear of the bend: on h = 0.8 mm the 50 Ω line is
// ~1.5 mm wide, so the 2-width fold consumes only ~3 mm of the half-wave
// (a 1.6 mm board's 3 mm line leaves no arm for the tap — the dims error
// out with TapNotRealizable there, by design).
const FBW: f64 = 0.22;
const FOLD_WIDTHS: f64 = 2.0;
const Z0_OHM: f64 = 50.0;
const EPS_R: f64 = 4.4;
const H_M: f64 = 0.8e-3;
const DX_M: f64 = 0.2e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const DRIVE_V0: f64 = 1.0;
const BW_HZ: f64 = 4.0e9; // drive envelope covers ~3–7 GHz
const N_STEPS: usize = 13000;
/// Feed length: long enough to carry a 3-probe triple with clearance.
const FEED_LEN_M: f64 = 12.0e-3;
/// Probe-triple spacing (βd ≈ 0.7 rad at f₀ — well-conditioned fit).
const SPACING_CELLS: usize = 12;

/// BO knob bounds: arm retune (the seed measured ~+17 % high, so the fix
/// is longer arms), tap retune, gap retune (the gap floor 2·dx is clamped
/// in `candidate_layout` — finer gaps are not honestly resolvable).
const BOUNDS: [(f64, f64); 3] = [(0.9, 1.4), (0.5, 1.3), (0.7, 1.1)];
/// Solve budget: n_initial + n_iters objective calls, one DUT solve each.
const BO_INITIAL: usize = 5;
const BO_ITERS: usize = 7;

fn fr4() -> Substrate {
    Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

fn filter_spec() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Butterworth,
        f0_hz: F0_HZ,
        fbw: FBW,
        order: Some(ORDER),
        z0_ohm: Z0_OHM,
        mask: SpecMask {
            passband_ripple_db: 3.0,
            return_loss_db: 10.0,
            stopband: vec![(6.5e9, 20.0)],
        },
    }
}

/// Build the knob-scaled hairpin layout. Tap is clamped onto the arm
/// (feed half-width clear of the bend) and gaps to the 2·dx grid floor,
/// so every candidate is a physical, resolvable geometry.
fn candidate_layout(dims: &HairpinDimensions, knobs: &[f64]) -> Layout {
    let (arm_scale, tap_scale, gap_scale) = (knobs[0], knobs[1], knobs[2]);
    let arm = dims.arm_length_m * arm_scale;
    let tap_max = arm - dims.line_width_m / 2.0 - DX_M;
    let tap = (dims.tap_offset_m * tap_scale).min(tap_max).max(DX_M);
    let gaps: Vec<f64> = dims
        .gaps_m
        .iter()
        .map(|g| (g * gap_scale).max(2.0 * DX_M))
        .collect();
    hairpin_bpf_sections(&HairpinSectionParams {
        substrate: fr4(),
        arm_length_m: arm,
        line_width_m: dims.line_width_m,
        fold_spacing_m: dims.fold_spacing_m,
        gaps_m: gaps,
        tap_offset_m: tap,
        feed_width_m: dims.line_width_m,
        feed_length_m: FEED_LEN_M,
    })
}

/// The reference: a straight Z₀ through line at the port height, spanning
/// port to port on the shared envelope bbox → the identical voxel grid.
fn reference_layout(seed: &Layout) -> Layout {
    let p0 = seed.ports[0].at;
    let p1 = seed.ports[1].at;
    let w = seed.ports[0].width_m;
    Layout {
        substrate: seed.substrate,
        traces: vec![Polygon::rect(p0.x, p0.y - w / 2.0, p1.x - p0.x, w)],
        ports: seed.ports.clone(),
        bbox: seed.bbox,
    }
}

/// Voxelize and express one run as a JobSpec; returns (spec, dt).
fn job_for(layout: &Layout) -> (JobSpec, f64) {
    let model: MicrostripModel = voxelize_microstrip(
        layout,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let k_top = model.port_cells[0].2;
    let load_cell = model.port_cells[1];
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;

    // Aperture / probe j band: the feed width centred on the tap height.
    let tap_y = layout.ports[0].at.y;
    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx - tap_y).abs() < w_feed / 2.0 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);
    assert!(j_hi > j_lo, "aperture band empty");
    let j_strip = (j_lo + j_hi) / 2;

    // Probe triples: A on the input feed, B on the output feed, both
    // ordered along +x with SPACING_CELLS spacing.
    let spacing_m = SPACING_CELLS as f64 * dx;
    let i_a0 = i_for(layout.ports[0].at.x + 2.4e-3);
    let i_b0 = i_for(layout.ports[1].at.x - 2.4e-3 - 2.0 * spacing_m);

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
        // Side-wall CPML, PEC ground/lid (S.9) — the board-level boundary.
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![
            AperturePortSpec {
                i: model.port_cells[0].0,
                j_lo,
                j_hi,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: DRIVE_V0,
                f0_hz: F0_HZ,
                bw_hz: BW_HZ,
                t0_steps,
            },
            AperturePortSpec {
                i: load_cell.0,
                j_lo,
                j_hi,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: 0.0,
                f0_hz: F0_HZ,
                bw_hz: BW_HZ,
                t0_steps,
            },
        ],
        probes: vec![
            mk_probe(i_a0),
            mk_probe(i_a0 + SPACING_CELLS),
            mk_probe(i_a0 + 2 * SPACING_CELLS),
            mk_probe(i_b0),
            mk_probe(i_b0 + SPACING_CELLS),
            mk_probe(i_b0 + 2 * SPACING_CELLS),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        backend: BackendChoice::Cpu,
    };
    (spec, dt)
}

fn run(spec: JobSpec) -> Vec<Vec<f64>> {
    let handle = yee_engine::submit(spec);
    let result = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");
    assert_eq!(result.steps_done, N_STEPS);
    result.probes
}

/// Directional |S21|(f) in dB for a DUT probe set against the shared
/// reference probe set.
fn s21_db(dut_p: &[Vec<f64>], ref_p: &[Vec<f64>], dt: f64, freqs: &[f64]) -> Vec<f64> {
    let spacing_m = SPACING_CELLS as f64 * DX_M;
    sparams::directional_transmission_db(
        [&dut_p[3], &dut_p[4], &dut_p[5]],
        [&ref_p[3], &ref_p[4], &ref_p[5]],
        dt,
        spacing_m,
        freqs,
    )
}

/// RMS misfit between two dB curves, both floored at −40 dB so the deep
/// stopband (measurement noise floor) does not dominate the fit.
fn misfit_db(measured: &[f64], designed: &[f64]) -> f64 {
    let floor = -40.0;
    let n = measured.len() as f64;
    (measured
        .iter()
        .zip(designed)
        .map(|(m, d)| (m.max(floor) - d.max(floor)).powi(2))
        .sum::<f64>()
        / n)
        .sqrt()
}

/// Passband features from a dB curve: (−3 dB centroid, −3 dB width, peak).
fn features(freqs: &[f64], db: &[f64]) -> (f64, f64, f64) {
    let peak = db.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let in_band: Vec<usize> = (0..db.len()).filter(|&n| db[n] >= peak - 3.0).collect();
    let lo = freqs[*in_band.first().unwrap()];
    let hi = freqs[*in_band.last().unwrap()];
    ((lo + hi) / 2.0, hi - lo, peak)
}

#[test]
#[ignore = "slow: ~13 multi-minute release FDTD solves (~1 h); engine-bpf-bo-001 gate (R.4) — run with --release --ignored"]
fn bo_closes_the_synthesized_hairpin_toward_its_coupling_matrix() {
    let project = synthesize(&filter_spec());
    let dims = dimension_hairpin_with_fold(&project, &fr4(), FOLD_WIDTHS)
        .expect("hairpin dimensioning failed");

    // Fixed envelope bbox at the knob extremes → one shared grid for the
    // reference and every BO candidate.
    let envelope = candidate_layout(&dims, &[BOUNDS[0].1, 1.0, BOUNDS[2].1]);
    let bbox = BBox {
        min: envelope.bbox.min,
        max: envelope.bbox.max,
    };
    let with_bbox = |mut l: Layout| {
        l.bbox = bbox;
        l
    };

    let seed = with_bbox(candidate_layout(&dims, &[1.0, 1.0, 1.0]));
    let reference = reference_layout(&seed);
    let (ref_spec, dt) = job_for(&reference);
    let ref_p = run(ref_spec);
    assert!(ref_p[3].iter().any(|v| *v != 0.0), "reference probe silent");

    // 3.5–6.5 GHz, 100 MHz raster (31 points — the BO objective's grid).
    let freqs: Vec<f64> = (0..=30).map(|n| 3.5e9 + n as f64 * 100.0e6).collect();
    let cm = coupling_matrix_s_params(&project.coupling, &freqs, F0_HZ, FBW);
    let cm_db: Vec<f64> = cm
        .iter()
        .map(|(_, s21)| 20.0 * s21.norm().log10())
        .collect();

    let measure = |knobs: &[f64]| -> Vec<f64> {
        let layout = with_bbox(candidate_layout(&dims, knobs));
        let (spec, dt2) = job_for(&layout);
        assert_eq!(dt, dt2, "candidate grid diverged from the reference");
        s21_db(&run(spec), &ref_p, dt, &freqs)
    };

    let seed_curve = measure(&[1.0, 1.0, 1.0]);
    let seed_misfit = misfit_db(&seed_curve, &cm_db);
    let (seed_fc, seed_bw, seed_peak) = features(&freqs, &seed_curve);
    eprintln!(
        "engine-bpf-bo-001: seed misfit {seed_misfit:.2} dB RMS | centre {:.3} GHz | \
         BW {:.0} MHz | peak {seed_peak:.2} dB",
        seed_fc / 1e9,
        seed_bw / 1e6
    );

    let evals = std::cell::RefCell::new(0usize);
    let result = minimize(
        |x: &DVector<f64>| {
            let knobs = [x[0], x[1], x[2]];
            let curve = measure(&knobs);
            let m = misfit_db(&curve, &cm_db);
            let mut n = evals.borrow_mut();
            *n += 1;
            eprintln!(
                "  BO eval {:>2}: arm {:.3} tap {:.3} gap {:.3} → misfit {m:.2} dB RMS",
                *n, knobs[0], knobs[1], knobs[2]
            );
            m
        },
        BOUNDS.to_vec(),
        BoConfig {
            n_initial: BO_INITIAL,
            n_iters: BO_ITERS,
            seed: 42,
            ..BoConfig::default()
        },
    );

    let best = [result.x_best[0], result.x_best[1], result.x_best[2]];
    let best_curve = measure(&best);
    let best_misfit = misfit_db(&best_curve, &cm_db);
    let (fc, bw, peak) = features(&freqs, &best_curve);
    eprintln!(
        "engine-bpf-bo-001: best knobs arm {:.3} tap {:.3} gap {:.3} → misfit \
         {best_misfit:.2} dB RMS (seed {seed_misfit:.2})",
        best[0], best[1], best[2]
    );
    for (f, (m, d)) in freqs.iter().zip(best_curve.iter().zip(&cm_db)) {
        if ((f / 1e8).round() as u64).is_multiple_of(5) {
            eprintln!(
                "  {:>5.2} GHz: measured {m:>7.2} dB | coupling-matrix {d:>7.2} dB",
                f / 1e9
            );
        }
    }
    eprintln!(
        "  best passband: centre {:.3} GHz (designed {:.2}) | BW {:.0} MHz (designed {:.0}) | \
         peak {peak:.2} dB",
        fc / 1e9,
        F0_HZ / 1e9,
        bw / 1e6,
        FBW * F0_HZ / 1e6
    );

    // Gate 1: BO must strictly improve the seed misfit — the EM loop earns
    // its cost.
    assert!(
        best_misfit < seed_misfit,
        "engine-bpf-bo-001 FAILED: BO did not improve the seed \
         ({best_misfit:.2} vs {seed_misfit:.2} dB RMS)"
    );
    // Gate 2: the optimized filter is a real band-pass at the right place —
    // a passband peak clear of the seed's buried −19 dB and a centre inside
    // the design neighbourhood. (Provisional walking-skeleton numbers from
    // the diagnosis; tightened by the first converged run, S.8 style.)
    assert!(
        peak >= -10.0,
        "engine-bpf-bo-001 FAILED: optimized passband peak {peak:.2} dB still buried"
    );
    let fc_err = (fc - F0_HZ).abs() / F0_HZ;
    assert!(
        fc_err <= 0.10,
        "engine-bpf-bo-001 FAILED: optimized centre {:.3} GHz vs {:.2} GHz (err {:.1} %)",
        fc / 1e9,
        F0_HZ / 1e9,
        fc_err * 100.0
    );
}
