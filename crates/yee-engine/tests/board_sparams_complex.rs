//! Gate `engine-sparams-002` (R.2, ADR-0195): **complex S-parameters +
//! Touchstone export** — the engine measures a complex plane-to-plane
//! transfer whose unwrapped phase matches transmission-line theory, and
//! the measured two-port round-trips through the project's Touchstone
//! writer.
//!
//! One release solve of the lossless 6 λ_g FR-4 line (the S.5-certified
//! scenario on the S.9/S.10 stack), two directional probe triples
//! d ≈ 3 λ_g apart. `sparams::forward_transfer` gives the complex
//! `T(f) = e^{−jβ(f)d}`:
//!
//! 1. **Phase gate**: the unwrapped-phase slope `dφ/df = −2π d √ε_eff/c`
//!    vs the Hammerstad–Jensen closed form, ±5 %.
//! 2. **Magnitude sanity**: |T| within ±1 dB of lossless unity across the
//!    band.
//! 3. **Touchstone**: assemble the measured two-port of the plane-A→B
//!    line segment — S21 = S12 = T measured, S11 = S22 = 0 (a uniform
//!    line in its own reference impedance has zero reflection; the Γ a
//!    probe triple reads on a THRU is the measurement fixture's
//!    load-port reflection, not the DUT's S11 — see ADR-0195, R.2b) —
//!    **enforce passivity** (the raw |T| ripples ≤ ~0.3 dB above unity,
//!    and `yee_io::touchstone::read` rejects λ_max(S†S) > 1 + 1e−9;
//!    scale offending samples by 1/σ_max, asserted ≤ 0.5 dB), write
//!    `.s2p`, read it back, and assert the data survives to the
//!    writer's precision.
//!
//! `#[ignore]`'d (one multi-minute release FDTD run):
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_sparams_complex -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use num_complex::Complex64;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
const N_STEPS: usize = 9000;
const SPACING_CELLS: usize = 17;

/// Scale a symmetric reciprocal 2-port sample `[[a, b], [b, a]]` to be
/// passive. Its singular values are `|a ± b|` in closed form, so if
/// σ_max > 1 divide the whole sample by σ_max. Returns the applied
/// σ_max so the caller can report the worst correction.
fn enforce_passivity_sym2(a: &mut Complex64, b: &mut Complex64) -> f64 {
    let sigma_max = (*a + *b).norm().max((*a - *b).norm());
    if sigma_max > 1.0 {
        *a /= sigma_max;
        *b /= sigma_max;
    }
    sigma_max
}

#[test]
fn passivity_enforcement_caps_sigma_max_at_unity() {
    // |a+b| = 1.3 dominates: expect both entries scaled by 1/1.3.
    let mut a = Complex64::new(0.3, 0.0);
    let mut b = Complex64::new(1.0, 0.0);
    let sigma = enforce_passivity_sym2(&mut a, &mut b);
    assert!((sigma - 1.3).abs() < 1e-12);
    assert!(((a + b).norm() - 1.0).abs() < 1e-12);
    // Already-passive sample is untouched.
    let mut a2 = Complex64::new(0.1, 0.05);
    let mut b2 = Complex64::new(0.6, -0.4);
    let (a0, b0) = (a2, b2);
    let sigma2 = enforce_passivity_sym2(&mut a2, &mut b2);
    assert!(sigma2 <= 1.0);
    assert_eq!((a2, b2), (a0, b0));
}

#[test]
#[ignore = "slow: one multi-minute release FDTD run; engine-sparams-002 gate (R.2) — run with --release --ignored"]
fn complex_transfer_phase_matches_theory_and_exports_to_touchstone() {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 6.0 * lam_g;
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

    let model = voxelize_microstrip(
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
    let load_cell = model.port_cells[1];
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_a0 = i_for(1.5 * lam_g);
    let i_b0 = i_for(4.5 * lam_g);
    let d_m = (i_b0 - i_a0) as f64 * dx;

    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < W_M / 2.0 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);

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

    let bw = 0.8 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

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
        aperture_ports: vec![
            AperturePortSpec {
                i: model.port_cells[0].0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: DRIVE_V0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
                record: false,
            },
            AperturePortSpec {
                i: load_cell.0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: 0.0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
                record: false,
            },
        ],
        thin_wires: vec![],
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

    // Offline-analysis escape hatch: dump the raw probe series once so fit
    // diagnostics can be iterated without re-running the FDTD solve.
    if let Ok(dir) = std::env::var("R2_DUMP_DIR") {
        let dump = serde_json::json!({
            "dt_s": dt,
            "spacing_m": spacing_m,
            "d_m": d_m,
            "e_eff": e_eff,
            "probes": p,
        });
        std::fs::write(
            std::path::Path::new(&dir).join("r2_probes.json"),
            serde_json::to_string(&dump).unwrap(),
        )
        .expect("probe dump failed");
    }

    // 4.0–6.0 GHz (well inside the drive band), 50 MHz raster.
    let freqs: Vec<f64> = (0..=40).map(|n| 4.0e9 + n as f64 * 50.0e6).collect();
    let transfer = sparams::forward_transfer(
        [&p[0], &p[1], &p[2]],
        [&p[3], &p[4], &p[5]],
        dt,
        spacing_m,
        &freqs,
    );

    // ---- 1. Unwrapped-phase slope vs HJ ----
    let mut phase: Vec<f64> = transfer.iter().map(|c| c.1.atan2(c.0)).collect();
    for n in 1..phase.len() {
        while phase[n] - phase[n - 1] > PI {
            phase[n] -= 2.0 * PI;
        }
        while phase[n] - phase[n - 1] < -PI {
            phase[n] += 2.0 * PI;
        }
    }
    // Least-squares slope dφ/df.
    let nf = freqs.len() as f64;
    let fm = freqs.iter().sum::<f64>() / nf;
    let pm = phase.iter().sum::<f64>() / nf;
    let slope = freqs
        .iter()
        .zip(&phase)
        .map(|(f, ph)| (f - fm) * (ph - pm))
        .sum::<f64>()
        / freqs.iter().map(|f| (f - fm).powi(2)).sum::<f64>();
    let slope_ref = -2.0 * PI * d_m * e_eff.sqrt() / C0_M_S;
    let slope_err = (slope - slope_ref).abs() / slope_ref.abs();

    // ---- 2. Magnitude sanity (lossless line) ----
    let worst_mag_db = transfer
        .iter()
        .map(|c| 20.0 * (c.0.hypot(c.1)).log10())
        .fold(
            0.0_f64,
            |acc, db| if db.abs() > acc.abs() { db } else { acc },
        );

    eprintln!(
        "engine-sparams-002: dφ/df = {slope:.4e} rad/Hz vs HJ {slope_ref:.4e} → err {:.2} % \
         | worst |T| deviation {worst_mag_db:+.2} dB over 4–6 GHz",
        slope_err * 100.0,
    );

    assert!(
        slope_err <= 0.05,
        "engine-sparams-002 FAILED: phase slope err {:.2} % (> 5 %)",
        slope_err * 100.0
    );
    assert!(
        worst_mag_db.abs() <= 1.0,
        "engine-sparams-002 FAILED: |T| deviates {worst_mag_db:+.2} dB from lossless unity"
    );

    // ---- 3. Touchstone round-trip of the measured two-port ----
    // The DUT is the plane-A→B segment of a uniform line: in its own
    // reference impedance S11 = S22 = 0 exactly. (The complex Γ a probe
    // triple reads on a THRU is the fixture's load-port reflection — and
    // the finite window truncates the port-to-port multi-bounce ring-down,
    // pushing raw |Γ| above 1 at the fwd-ripple minima. Exporting a
    // *measured* S11 needs de-embedding/calibration: R.2b, ADR-0195.)
    let mut worst_sigma = 0.0_f64;
    let data: Vec<Vec<Complex64>> = transfer
        .iter()
        .map(|t| {
            let mut s21 = Complex64::new(t.0, t.1);
            let mut s11 = Complex64::new(0.0, 0.0);
            // Raw measured |T| sits fractionally above unity (the ±1 dB
            // ripple asserted above); yee-io's read() rejects non-passive
            // matrices, so enforce passivity before export.
            let sigma = enforce_passivity_sym2(&mut s11, &mut s21);
            worst_sigma = worst_sigma.max(sigma);
            // Reciprocal, symmetric passive line: row-major [S11 S12; S21 S22].
            vec![s11, s21, s21, s11]
        })
        .collect();
    eprintln!(
        "  passivity enforcement: worst raw sigma_max = {worst_sigma:.6} \
         ({:+.3} dB)",
        20.0 * worst_sigma.log10()
    );
    assert!(
        worst_sigma < 1.06,
        "engine-sparams-002 FAILED: passivity correction {:.3} dB (> 0.5 dB) — \
         the export would silently reshape the measurement",
        20.0 * worst_sigma.log10()
    );
    let file = yee_io::touchstone::File {
        n_ports: 2,
        z0: Z0_OHM,
        freq_unit: yee_io::touchstone::FreqUnit::GHz,
        format: yee_io::touchstone::Format::RealImag,
        freq_hz: freqs.clone(),
        data,
        comments: vec!["engine-sparams-002: engine-measured 6 lambda_g FR-4 line".into()],
    };
    let path = std::env::temp_dir().join("engine_sparams_002.s2p");
    yee_io::touchstone::write(&path, &file).expect("touchstone write failed");
    let back = yee_io::touchstone::read(&path).expect("touchstone read failed");
    assert_eq!(back.n_ports, 2);
    assert_eq!(back.freq_hz.len(), file.freq_hz.len());
    for (k, (a, b)) in file.data.iter().zip(&back.data).enumerate() {
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            let err = (x - y).norm();
            assert!(
                err < 1e-6 * x.norm().max(1e-12) + 1e-9,
                "touchstone round-trip diverged at point {k} entry {i}: {x} vs {y}"
            );
        }
    }
    eprintln!(
        "  .s2p round-trip OK ({} points) at {}",
        back.freq_hz.len(),
        path.display()
    );
}
