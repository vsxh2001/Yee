//! Gate `engine-automesh-002` (FS.0b.2b, ADR-0216): **push-button meshing
//! on graded grids** — the S.6 λ/4 open-stub notch board with no hand-set
//! dx anywhere, run through [`converge_two_port_graded`]: each pass builds
//! its grid from the FS.0b.1 rulebook and refines the one `scale` knob by
//! 1/√2 until the |S21| curve stops moving (the identical linear-ΔS
//! criterion to `engine-automesh-001`). The gate holds the converged notch
//! against TL theory AND the graded economics: every pass must come in
//! under 0.35× the cells of the equivalent-resolution uniform grid.
//!
//! `#[ignore]`'d (2 solves per pass, release):
//!
//! ```bash
//! cargo test -p yee-engine --release --test automesh_graded -- --ignored --nocapture
//! ```

use std::f64::consts::SQRT_2;

use yee_engine::automesh::{auto_dx_bulk, auto_spacings, converge_two_port_graded};
use yee_engine::board::GradedBoardOptions;
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff, open_end_delta_l};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const Z0_OHM: f64 = 50.0;
const F_MAX_HZ: f64 = 6.0e9;

/// The S.6 scenario (same board as `engine-automesh-001`): a through line
/// with a λ/4 open stub at mid-length.
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
#[ignore = "slow: 2 release FDTD solves per pass; engine-automesh-002 gate (FS.0b.2b) — run with --release --ignored"]
fn graded_push_button_meshing_converges_to_the_tl_notch() {
    let layout = stub_layout();
    let opts = GradedBoardOptions::for_board(&layout, F_MAX_HZ, F0_HZ, 0.8 * F0_HZ);
    eprintln!(
        "engine-automesh-002: coarse = {:.3} mm (rulebook), growth {}, guard {:.1} mm",
        auto_dx_bulk(&layout, F_MAX_HZ) * 1e3,
        opts.mesh.growth,
        opts.mesh.guard_m * 1e3
    );

    let freqs: Vec<f64> = (0..=50).map(|n| 3.5e9 + n as f64 * 50.0e6).collect();
    // Same linear ΔS tolerance as engine-automesh-001 (0.20, the measured
    // walking-skeleton value; ADR-0204's rationale applies verbatim).
    let result = converge_two_port_graded(&layout, F_MAX_HZ, opts.clone(), &freqs, 0.20, 3)
        .expect("graded convergence loop failed");

    // Per-pass report + the graded economics assert: each pass vs the
    // uniform grid at the same resolution (its own extents at fine_m).
    for (n, pass) in result.passes.iter().enumerate() {
        let mut mesh = opts.mesh;
        mesh.scale = 1.0 / SQRT_2.powi(n as i32);
        // npml in coarse cells, as the loop set it for this pass.
        let coarse0 = auto_dx_bulk(&layout, F_MAX_HZ);
        mesh.npml = ((opts.mesh.npml as f64 * coarse0) / pass.dx_m)
            .round()
            .max(1.0) as usize;
        let spac = auto_spacings(&layout, F_MAX_HZ, &mesh).expect("rulebook failed");
        assert_eq!(
            spac.cell_count(),
            pass.cells,
            "pass {n}: reconstructed spacings disagree with the loop's grid"
        );
        let ext = |v: &[f64]| v.iter().sum::<f64>();
        let uniform_eq = (ext(&spac.dx) / spac.fine_m).ceil()
            * (ext(&spac.dy) / spac.fine_m).ceil()
            * (ext(&spac.dz) / spac.fine_m).ceil();
        let ratio = pass.cells as f64 / uniform_eq;
        let (i_min, db_min) = pass
            .s21_db
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .unwrap();
        eprintln!(
            "  pass {n}: coarse = {:.3} mm, fine = {:.3} mm, {} cells \
             ({:.3}x uniform-eq {:.2e}) → notch {:.3} GHz @ {:.1} dB",
            pass.dx_m * 1e3,
            spac.fine_m * 1e3,
            pass.cells,
            ratio,
            uniform_eq,
            freqs[i_min] / 1e9,
            db_min
        );
        assert!(
            ratio <= 0.35,
            "engine-automesh-002 FAILED: pass {n} used {:.3}x the equivalent-resolution \
             uniform cells (need ≤ 0.35)",
            ratio
        );
    }
    eprintln!(
        "  final Δ|S| = {:.4} (linear), converged = {}",
        result.final_delta, result.converged
    );

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
        "engine-automesh-002 FAILED: converged notch {:.3} GHz vs {:.1} GHz (err {:.1} % > 5 %)",
        f_notch / 1e9,
        F0_HZ / 1e9,
        err * 100.0
    );
    assert!(
        db_min <= -20.0,
        "engine-automesh-002 FAILED: converged notch only {db_min:.1} dB deep (need ≤ −20)"
    );
    assert!(
        result.converged,
        "engine-automesh-002 FAILED: not converged (final Δ|S| = {:.4})",
        result.final_delta
    );
}
