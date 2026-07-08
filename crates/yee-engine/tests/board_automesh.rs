//! Gate `engine-automesh-001` (FS.0a, ADR-0204): **push-button meshing**
//! — the S.6 λ/4 open-stub notch board solved with *no hand-set dx
//! anywhere*: [`yee_engine::automesh::auto_dx`] seeds the grid from the
//! layout + drive band, and [`converge_two_port`] refines dx by 1/√2 per
//! pass until the |S21| curve stops moving. The gate then holds the
//! converged notch frequency against transmission-line theory
//! (f_notch = c/(4·(l_stub + ΔL)·√ε_eff), the S.6 reference) — the
//! "novice gets a trustworthy answer without knowing the λ/20 rules"
//! criterion, machine-checked.
//!
//! `#[ignore]`'d (2 solves per pass × up to 3 passes, release):
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_automesh -- --ignored --nocapture
//! ```

use yee_engine::automesh::{auto_dx, converge_two_port};
use yee_engine::board::{TwoPortBoardOptions, reference_through_line, two_port_board_job};
use yee_engine::sparams::{fit_standing_wave, single_bin_dft};
use yee_engine::{JobEvent, submit};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff, open_end_delta_l};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const Z0_OHM: f64 = 50.0;

/// The S.6 scenario: a through line with a λ/4 open stub at mid-length.
fn stub_layout() -> Layout {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - open_end_delta_l(W_M, H_M, e_eff);
    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    let traces = vec![line, stub];
    let bbox = BBox::from_polygons(&traces);
    Layout {
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
    }
}

/// Manual diagnostic (not a gate): one pass of the automesh scenario at
/// `YEE_AUTOMESH_DX_MM` (default 0.267), dumping the raw three-probe
/// wave-split (|fwd|, |bwd|, fitted β vs closed-form, fit residual) for
/// reference and DUT at every bin — the forensic view of the pass-2
/// measurement blowup.
#[test]
#[ignore = "manual diagnostic: 2 release solves, dumps raw wave-splits"]
fn automesh_pass_fit_diagnostics() {
    let dx = std::env::var("YEE_AUTOMESH_DX_MM")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.267)
        * 1e-3;
    let layout = stub_layout();
    let reference = reference_through_line(&layout);

    // The same settings the convergence loop uses at this dx (constant
    // physical window/margins/absorber relative to the auto_dx pass 0).
    let dx0 = auto_dx(&layout, 6.0e9);
    let mut opts = TwoPortBoardOptions::for_band(F0_HZ, 0.8 * F0_HZ);
    opts.dx_m = dx;
    opts.n_steps = (9000.0 * 0.3e-3 / dx).round() as usize;
    opts.margin_cells = (34.0 * dx0 / dx).round() as usize;
    opts.air_above_cells = (34.0 * dx0 / dx).round() as usize;
    opts.npml = (10.0 * dx0 / dx).round() as usize;
    opts.spacing_cells = (12.0 * dx0 / dx).round() as usize;

    let run = |l: &Layout| -> (Vec<Vec<f64>>, f64, f64) {
        let job = two_port_board_job(l, &opts).expect("job build failed");
        let (dt, spacing) = (job.dt_s, job.spacing_m);
        let handle = submit(job.spec);
        for event in handle.events() {
            match event {
                JobEvent::Done { result } => return (result.probes, dt, spacing),
                JobEvent::Error { message } => panic!("job failed: {message}"),
                _ => {}
            }
        }
        panic!("engine stream ended without a result");
    };
    let (ref_p, dt, spacing) = run(&reference);
    let (dut_p, _, _) = run(&layout);

    let e_eff = eps_eff(W_M, H_M, EPS_R);
    eprintln!(
        "automesh-diag: dx = {:.3} mm, spacing = {:.3} mm, npml = {}, margin = {}",
        dx * 1e3,
        spacing * 1e3,
        opts.npml,
        opts.margin_cells
    );
    eprintln!(
        "  f/GHz | refA|fwd| dutA|fwd| dutA/refA_dB | refB|fwd| dutB|fwd| B-ratio_dB | double-ratio S21 dB | beta_hj refB_beta dutB_resid"
    );
    for n in 0..=25 {
        let f = 3.5e9 + n as f64 * 100.0e6;
        let mag = |c: (f64, f64)| (c.0 * c.0 + c.1 * c.1).sqrt();
        let split = |p: &[Vec<f64>], k0: usize| {
            let v: Vec<(f64, f64)> = (k0..k0 + 3).map(|k| single_bin_dft(&p[k], dt, f)).collect();
            fit_standing_wave(v[0], v[1], v[2], spacing)
        };
        let (ra, rb) = (split(&ref_p, 0), split(&ref_p, 3));
        let (da, db) = (split(&dut_p, 0), split(&dut_p, 3));
        let beta_hj = 2.0 * std::f64::consts::PI * f * e_eff.sqrt() / C0_M_S;
        let t_ref = mag(rb.fwd) / mag(ra.fwd);
        let t_dut = mag(db.fwd) / mag(da.fwd);
        eprintln!(
            "  {:5.2} | {:.3e} {:.3e} {:6.2} | {:.3e} {:.3e} {:6.2} | {:6.2} | {:7.1} {:7.1} {:.3}",
            f / 1e9,
            mag(ra.fwd),
            mag(da.fwd),
            20.0 * (mag(da.fwd) / mag(ra.fwd)).log10(),
            mag(rb.fwd),
            mag(db.fwd),
            20.0 * (mag(db.fwd) / mag(rb.fwd)).log10(),
            20.0 * (t_dut / t_ref).log10(),
            beta_hj,
            rb.beta_rad_m,
            db.residual,
        );
    }
}

#[test]
#[ignore = "slow: up to 6 release FDTD solves; engine-automesh-001 gate (FS.0a) — run with --release --ignored"]
fn push_button_meshing_converges_to_the_tl_notch() {
    let layout = stub_layout();
    let reference = reference_through_line(&layout);

    // No hand-set dx: the rulebook picks it from the layout + band.
    let f_max = 6.0e9;
    let dx0 = auto_dx(&layout, f_max);
    eprintln!(
        "engine-automesh-001: auto_dx = {:.3} mm (λ/20 = {:.3}, h/3 = {:.3}, feature/2 = {:.3})",
        dx0 * 1e3,
        C0_M_S / (f_max * EPS_R.sqrt()) / 20.0 * 1e3,
        H_M / 3.0 * 1e3,
        yee_engine::automesh::min_feature_m(&layout) / 2.0 * 1e3,
    );

    let mut opts = TwoPortBoardOptions::for_band(F0_HZ, 0.8 * F0_HZ);
    opts.dx_m = dx0;
    // Base time window: the R.0-family scenarios ring down in ~5 ns; the
    // loop rescales steps as dx shrinks.
    opts.n_steps = (9000.0 * 0.3e-3 / dx0).round() as usize;

    let freqs: Vec<f64> = (0..=50).map(|n| 3.5e9 + n as f64 * 50.0e6).collect();
    // Linear ΔS tolerance: HFSS's reference point is ~0.02; staircased
    // uniform FDTD gets 0.10 at walking-skeleton fidelity (FS.0b's graded
    // grid tightens this).
    let result = converge_two_port(&layout, &reference, opts, &freqs, 0.10, 3)
        .expect("convergence loop failed");

    for (n, pass) in result.passes.iter().enumerate() {
        let (i_min, db_min) = pass
            .s21_db
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .unwrap();
        eprintln!(
            "  pass {n}: dx = {:.3} mm → notch {:.3} GHz @ {:.1} dB",
            pass.dx_m * 1e3,
            freqs[i_min] / 1e9,
            db_min
        );
    }
    eprintln!(
        "  final Δ|S| = {:.4} (linear), converged = {}",
        result.final_delta, result.converged
    );
    // Per-bin dump of the last two passes (diagnostic).
    if result.passes.len() >= 2 {
        let prev = &result.passes[result.passes.len() - 2];
        let last = result.passes.last().unwrap();
        let lin = |db: f64| 10.0_f64.powf(db / 20.0);
        for (i, f) in freqs.iter().enumerate() {
            let (a, b) = (prev.s21_db[i], last.s21_db[i]);
            eprintln!(
                "    {:.2} GHz: {:8.2} dB → {:8.2} dB (Δlin {:.4})",
                f / 1e9,
                a,
                b,
                (lin(a) - lin(b)).abs()
            );
        }
    }

    // The converged notch frequency vs TL theory (the S.6 reference).
    let last = result.passes.last().unwrap();
    let (i_min, &db_min) = last
        .s21_db
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .unwrap();
    let f_notch = freqs[i_min];
    let err = (f_notch - F0_HZ).abs() / F0_HZ;
    eprintln!(
        "  converged notch {:.3} GHz vs designed {:.1} GHz → err {:.1} % (depth {:.1} dB)",
        f_notch / 1e9,
        F0_HZ / 1e9,
        err * 100.0,
        db_min
    );

    assert!(
        err <= 0.05,
        "engine-automesh-001 FAILED: converged notch {:.3} GHz vs {:.1} GHz (err {:.1} % > 5 %)",
        f_notch / 1e9,
        F0_HZ / 1e9,
        err * 100.0
    );
    assert!(
        db_min <= -20.0,
        "engine-automesh-001 FAILED: converged notch only {db_min:.1} dB deep (need ≤ −20)"
    );
    // The loop's own verdict is part of the contract.
    assert!(
        result.converged,
        "engine-automesh-001 FAILED: not converged (final Δ|S| = {:.4})",
        result.final_delta
    );
}
