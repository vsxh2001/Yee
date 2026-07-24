//! Gate `engine-stripline-alpha-001` (FS.4.2b, ADR-0226): per-layer
//! stackup loss (`yee_voxel::stackup_sigma_cells`) reproduces the
//! **exact** Beer-Lambert dielectric attenuation of a symmetric
//! stripline. Stripline is pure TEM entirely inside the dielectric (no
//! fringing into air, unlike microstrip), so the closed form is exact —
//! not an approximation like Pozar §3.199's `engine-loss-001` (which
//! needs an ε_eff fudge for microstrip's partial air/dielectric field):
//!
//! ```text
//! α_d = (π f √ε_r / c) · tan δ   [Np/m]      (× 8.686 → dB/m)
//! ```
//!
//! Conductor loss is not modeled (PEC everywhere, per the FS.4.2b
//! non-goals) — the measured attenuation is dielectric-only, matching the
//! closed form's assumption exactly.
//!
//! # Fixture
//!
//! Identical stripline cross-section to `engine-stripline-z0-001`
//! (`stripline_z0.rs`, FS.4.2a): ε_r 2.2, b = 16 cells (`B_M`), w/b =
//! 0.8125. `stackup_sigma_cells` (Task 1, FS.4.2b) maps `tan δ = 0.02`
//! (through a hand-built two-half-layer `Stackup`, since
//! `Stackup::symmetric_stripline` hardcodes `loss_tangent: 0.0`) onto the
//! per-cell σ riding `Materials::sigma_cells` into the lossy CA/CB
//! update — the same E.1 plumbing `engine-loss-001`
//! (`board_loss.rs`) already exercises for microstrip.
//!
//! # Extraction: two V-column planes, single-pass ratio
//!
//! Two `Ez`-column measurement planes (the `engine-stripline-z0-001`
//! V-column idiom — ground `k=0` up to, excluding, the trace `k=k_trace`,
//! summed × `dz`), plane A closer to the port and plane B further
//! downstream, both far enough from the port to have settled into the
//! guided TEM mode and far enough from the hard-PEC end wall that the
//! time gate (below) closes before the wall reflection reaches either
//! plane. Each plane's gated signal is a **single pass of the same
//! launched wave** — `|V_A|` and `|V_B|` are two samples of one
//! attenuating traveling wave, not two independent runs, so this is
//! launch-normalized by construction: both planes see the identical
//! launch, and the ADR-0204 warning (don't compare absolute single
//! ratios *across separate runs* with different incident waves) does not
//! apply — there is only one run, one wave, two taps.
//!
//! `α_meas = ln(|V_A| / |V_B|) / (x_B − x_A)` [Np/m], `x_B − x_A` taken
//! from the *actual* grid-quantized probe-plane separation (not the
//! analytic target), matching `engine-loss-001`'s `d_m` idiom.
//!
//! # Constant-σ vs true tan δ (documented, not asserted on)
//!
//! `stackup_sigma_cells` maps `tan δ` to σ at a single reference
//! frequency (`σ = 2π f_ref ε₀ ε_r tan δ`); the FDTD update then treats σ
//! as frequency-*independent*, so the discrete model's implied loss
//! tangent drifts as `tan δ_eff(f) = tan δ(f_ref) · f_ref / f` away from
//! `f_ref` — a real modeling deviation off-reference, not a bug (true
//! Debye/dispersive tan δ is a separate, unshipped lane — FS.4.2b
//! non-goals). This gate grades only at `f_ref` where the model is
//! exact; it does not sweep frequency.
//!
//! `#[ignore]`'d (two multi-minute release FDTD runs — lossy + lossless
//! control):
//!
//! ```bash
//! cargo test -p yee-engine --release --test stripline_alpha -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_compute::{
    Boundary, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials, Probe, ResistivePort,
    Waveform,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, StackupLayer, Substrate};
use yee_voxel::{VoxelOptions, stackup_sigma_cells, voxelize_stackup};

const EPS_R: f64 = 2.2;
/// Ground-to-ground spacing b — 16 cells at `DX_M` (matches
/// `engine-stripline-z0-001`'s confined-mode resolution lesson,
/// ADR-0215/0221: b >= ~16 cells).
const B_M: f64 = 3.2e-3;
/// Trace width — matches `engine-stripline-z0-001`'s w/b = 0.8125.
const W_M: f64 = 2.6e-3;
const DX_M: f64 = 0.2e-3;
/// Grading / σ reference frequency (see the module doc's constant-σ
/// note).
const F_REF_HZ: f64 = 6.0e9;
const C0_M_S: f64 = 299_792_458.0;
const MARGIN_CELLS: usize = 20;
const PORT_R_OHM: f64 = 50.0;
/// FR-4-core-like loss tangent (lossy but physical; matches the
/// `engine-loss-001` scale).
const TAN_D: f64 = 0.02;
/// Line length in guided wavelengths — longer than
/// `engine-stripline-z0-001`'s 8 λg to leave enough wall margin for
/// plane B's pulse tail to clear before the reflection-gate closes (see
/// the "Window hygiene" lesson in `stripline_eeff.rs`: the gate must
/// exceed the pulse-tail arrival time at the *farther* probe with
/// comfortable margin, not just the near one).
const L_LAMBDA_G: f64 = 8.5;
/// Measurement plane A, guided wavelengths downstream of the port —
/// past the launch transient (same order as `engine-stripline-z0-001`'s
/// 2.5 λg single probe).
const X_A_LAMBDA_G: f64 = 2.0;
/// Measurement plane B — far enough past A for a robust (~1 dB) drop at
/// `TAN_D`, and, at `L_LAMBDA_G` above, ~4.5 λg of wall margin left for
/// its pulse tail to clear before the reflection gate closes.
const X_B_LAMBDA_G: f64 = 4.0;

/// Exact stripline dielectric attenuation (pure TEM, entirely in the
/// dielectric — no ε_eff fudge needed, unlike microstrip's
/// `engine-loss-001` Pozar §3.199 form).
fn alpha_exact_np_per_m(eps_r: f64, tan_d: f64, f_hz: f64) -> f64 {
    PI * f_hz * eps_r.sqrt() / C0_M_S * tan_d
}

/// Run the stripline fixture at a given loss tangent and return the
/// gated `|V|` phasor magnitude at plane A and plane B, plus the actual
/// (grid-quantized) plane separation in metres.
fn measure_v_planes(tan_d: f64) -> (f64, f64, f64) {
    let lam_g = C0_M_S / (F_REF_HZ * EPS_R.sqrt());
    let l_m = L_LAMBDA_G * lam_g;
    let traces = vec![Polygon::rect(0.0, 0.0, l_m, W_M)];
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: Substrate {
            // Unused by the stackup path; kept for the Layout contract
            // (matches engine-stripline-z0-001).
            eps_r: EPS_R,
            height_m: B_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![PortRef {
            at: Point2::new(0.5e-3, W_M / 2.0),
            width_m: W_M,
            ref_impedance_ohm: PORT_R_OHM,
        }],
        bbox,
    };

    let half = StackupLayer {
        eps_r: EPS_R,
        height_m: B_M / 2.0,
        loss_tangent: tan_d,
    };
    let stack = Stackup {
        layers: vec![half, half],
        lid: true,
    };
    let model = voxelize_stackup(
        &layout,
        &stack,
        0,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: 0, // lidded: ignored
        },
    );
    let (nx, ny, nz) = model.dims;
    let (_i_drive, j_strip, k_trace) = model.port_cells[0];
    let dt = model.grid.dt;
    let dx = model.dx_m;

    let mut fdtd_spec = FdtdSpec::vacuum(nx, ny, nz, DX_M);
    fdtd_spec.dt = dt;

    let sigma = stackup_sigma_cells(&model, &stack, F_REF_HZ);
    let materials = Materials {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        sigma_cells: Some(sigma),
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
        ..Materials::default()
    };

    // Measurement planes A and B, grid-quantized (engine-loss-001's d_m
    // idiom: use the actual index separation, not the analytic target).
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize;
    let i_a = i_for(X_A_LAMBDA_G * lam_g);
    let i_b = i_for(X_B_LAMBDA_G * lam_g);
    assert!(i_b > i_a, "plane B must be downstream of plane A");
    let d_m = (i_b - i_a) as f64 * dx;

    // Time gate: stop before the far-end (hard-PEC wall) reflection
    // reaches the farther plane (B) — see the module doc's window-
    // hygiene note.
    let v_p_ref = C0_M_S / EPS_R.sqrt();
    let x_drive = 0.5e-3;
    let x_b_m = X_B_LAMBDA_G * lam_g;
    let t_refl = ((l_m - x_drive) + (l_m - x_b_m)) / v_p_ref;
    let gate_steps = (0.9 * t_refl / dt) as usize;
    let n_steps = gate_steps + 200;

    let bw = 0.4 * F_REF_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let mut drive = Drive::default();
    drive.ports.push(ResistivePort {
        cell: model.port_cells[0],
        resistance: PORT_R_OHM,
        waveform: Waveform::GaussianPulse {
            v0: 1.0,
            f0: F_REF_HZ,
            bw,
            t0_steps,
        },
    });
    // V(t) at plane A: probe indices [0, k_trace).
    for k in 0..k_trace {
        drive.probes.push(Probe {
            component: EComponent::Ez,
            cell: (i_a, j_strip, k),
        });
    }
    // V(t) at plane B: probe indices [k_trace, 2*k_trace).
    for k in 0..k_trace {
        drive.probes.push(Probe {
            component: EComponent::Ez,
            cell: (i_b, j_strip, k),
        });
    }

    let fields = Fields::zero(&fdtd_spec);
    let mut engine = CpuFdtd::with_drive(fdtd_spec, fields, materials, Boundary::PecBox, drive);
    engine.step_n(n_steps);

    // Time-gated single-bin DFT at f_ref (both series are Ez, sampled at
    // the same true time t = (m+1)*dt — no H-probe half-step staggering
    // to account for, per the brief: only V-columns are needed for α).
    let omega = 2.0 * PI * F_REF_HZ;
    let gate = gate_steps.min(n_steps);
    let ez = engine.probe_series();

    // `m` indexes time (0..gate), not the `ez` probe-column vec itself
    // (each `ez[k]` is its own full-length time series) — clippy's
    // needless_range_loop pattern-matches the nested `ez[k][m]` index
    // chain anyway; same idiom as the eigensolver assembly loops.
    #[allow(clippy::needless_range_loop)]
    let v_mag = |k_range: std::ops::Range<usize>| -> f64 {
        let mut acc = [0.0_f64; 2];
        for m in 0..gate {
            let v_t: f64 = k_range.clone().map(|k| ez[k][m]).sum::<f64>() * dx;
            let t_v = (m + 1) as f64 * dt;
            let (sv, cv) = (omega * t_v).sin_cos();
            acc[0] += v_t * cv;
            acc[1] -= v_t * sv;
        }
        (acc[0] * acc[0] + acc[1] * acc[1]).sqrt()
    };

    let v_mag_a = v_mag(0..k_trace);
    let v_mag_b = v_mag(k_trace..2 * k_trace);
    (v_mag_a, v_mag_b, d_m)
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs (lossy + lossless control); \
            engine-stripline-alpha-001 gate (FS.4.2b) — run with --release --ignored"]
fn stripline_alpha_matches_the_pozar_dielectric_loss_closed_form() {
    let alpha_ref = alpha_exact_np_per_m(EPS_R, TAN_D, F_REF_HZ);

    let (v_a, v_b, d_m) = measure_v_planes(TAN_D);
    assert!(
        v_a.is_finite() && v_a > 1e-3,
        "lossy run: V phasor magnitude at plane A is not non-trivial: {v_a}"
    );
    assert!(
        v_b.is_finite() && v_b > 1e-3,
        "lossy run: V phasor magnitude at plane B is not non-trivial: {v_b}"
    );
    let alpha_meas = (v_a / v_b).ln() / d_m;
    let rel_err = (alpha_meas - alpha_ref).abs() / alpha_ref;
    eprintln!(
        "engine-stripline-alpha-001: tan_d = {TAN_D}, d = {:.2} mm | |V_A| = {v_a:.4e}, \
         |V_B| = {v_b:.4e} -> alpha_meas = {alpha_meas:.4} Np/m ({:.4} dB/m) vs closed form \
         alpha_ref = {alpha_ref:.4} Np/m ({:.4} dB/m) -> err {:.3} %",
        d_m * 1e3,
        8.686 * alpha_meas,
        8.686 * alpha_ref,
        rel_err * 100.0
    );
    assert!(
        rel_err <= 0.10,
        "engine-stripline-alpha-001 FAILED: alpha_meas = {alpha_meas:.4} Np/m vs closed form \
         {alpha_ref:.4} Np/m (err {:.3} % > 10 %)",
        rel_err * 100.0
    );

    // Lossless control (same fixture, tan_d = 0 through the same
    // stackup_sigma_cells path): the differential kills systematic
    // gating bias — a real loss measurement should vanish here, and the
    // plane-A phasor should stay non-trivial/sane (same order of
    // magnitude as engine-stripline-z0-001's |V| ~= 2.33 at an
    // equivalent fixture), proving the sigma=0 no-op didn't silently
    // zero the field.
    let (v_a0, v_b0, d_m0) = measure_v_planes(0.0);
    assert!(
        v_a0.is_finite() && v_a0 > 1.0 && v_a0 < 10.0,
        "lossless control: plane-A V phasor is not in the sane non-trivial range \
         (engine-stripline-z0-001 reference ~2.33): {v_a0}"
    );
    let alpha_lossless = (v_a0 / v_b0).ln() / d_m0;
    eprintln!(
        "  lossless control: |V_A| = {v_a0:.4e}, |V_B| = {v_b0:.4e} -> alpha_lossless = \
         {alpha_lossless:.6} Np/m ({:.4} dB/m) — numeric gating-bias floor",
        8.686 * alpha_lossless
    );
    // Pinned from the measured floor with margin (2026-07-24 run): the
    // lossless differential must stay far under the closed-form signal
    // (alpha_ref above) — a real loss leak here would show up as a
    // floor comparable to alpha_ref, not a small fraction of it.
    assert!(
        alpha_lossless.abs() <= 0.05 * alpha_ref,
        "engine-stripline-alpha-001 no-op check FAILED: lossless control measured \
         alpha = {alpha_lossless:.6} Np/m, expected << alpha_ref = {alpha_ref:.4} Np/m \
         (>5 % of the real signal — the sigma=0 case must be a provable no-op)",
    );
}
