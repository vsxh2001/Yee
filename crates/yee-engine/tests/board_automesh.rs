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
use yee_engine::board::{TwoPortBoardOptions, reference_through_line};
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
