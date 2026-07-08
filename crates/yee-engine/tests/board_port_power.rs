//! Gate `engine-power-001` (FS.2a, ADR-0207): **energy bookkeeping from
//! port records** — the accepted-power observable that FS.2b's gain in
//! dBi normalizes by. Both aperture ports of a lossless matched through
//! line record their per-step `(v_src, v_terminal, i_branch)`. The
//! honest identity is on CIRCUIT-side quantities (**measured closure:
//! 0.9917**, with accepted-by-field = 51.3 % of the EMF supply — the
//! textbook matched-source halving): the EMF's supply
//! (Σ −v_src·i·dt) is dissipated in the two branch resistors (Σ i²R·dt
//! each) plus what the field keeps/loses — so the resistor share closes
//! to a large-but-<1 fraction, and P_accepted for FS.2b's gain is
//! E_emf − E_R(A). The naive field-side v_term·i is printed but NOT
//! asserted: it mixes real transfer with the reversible β-term storage
//! (β = dt·h/2ε₀A ≈ 14.5 Ω here, comparable to R — the first run's
//! measured lesson: it read a non-physical B/A ratio of 1.596).
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_port_power -- --ignored --nocapture
//! ```

use yee_engine::JobEvent;
use yee_engine::board::{TwoPortBoardOptions, two_port_board_job};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};

const F0_HZ: f64 = 5.0e9;
const BW_HZ: f64 = 4.0e9;
const W_M: f64 = 3.0e-3;
const L_M: f64 = 40.0e-3;

fn line_layout() -> Layout {
    let traces = vec![Polygon::rect(0.0, 0.0, L_M, W_M)];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: 50.0,
            },
            PortRef {
                at: Point2::new(L_M - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: 50.0,
            },
        ],
        bbox,
    }
}

#[test]
#[ignore = "slow: one ~1 min release FDTD run; engine-power-001 gate (FS.2a) — run with --release --ignored"]
fn port_energy_bookkeeping_closes_on_a_lossless_line() {
    let layout = line_layout();
    let mut opts = TwoPortBoardOptions::for_band(F0_HZ, BW_HZ);
    opts.dx_m = 0.4e-3;
    opts.n_steps = 7000;
    let mut job = two_port_board_job(&layout, &opts).expect("job build failed");
    for p in &mut job.spec.aperture_ports {
        p.record = true;
    }
    let dt = job.dt_s;

    let handle = yee_engine::submit(job.spec);
    let result = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");
    let records = result.port_records.expect("no port records returned");
    assert_eq!(records.len(), 2, "both ports recorded");
    assert_eq!(records[0].len(), opts.n_steps, "one sample per step");

    let r_ohm = 50.0;
    // Circuit-side quantities, unambiguous per branch:
    //   E_emf  = Σ −v_src·i·dt   (energy the EMF supplies)
    //   E_r    = Σ i²R·dt        (energy the branch resistor dissipates)
    // Field-side naive: E_ap = Σ v_term·i·dt (aperture → branch), which
    // mixes real transfer with the reversible β-term storage.
    let breakdown = |rec: &[(f64, f64, f64)]| {
        let (mut e_emf, mut e_r, mut e_ap) = (0.0, 0.0, 0.0);
        for &(v_src, v_term, i) in rec {
            e_emf += -v_src * i * dt;
            e_r += i * i * r_ohm * dt;
            e_ap += v_term * i * dt;
        }
        (e_emf, e_r, e_ap)
    };
    let (a_emf, a_r, a_ap) = breakdown(&records[0]);
    let (b_emf, b_r, b_ap) = breakdown(&records[1]);
    eprintln!("engine-power-001 energy breakdown (J):");
    eprintln!("  A: E_emf = {a_emf:.4e}, E_R = {a_r:.4e}, E_ap(v·i) = {a_ap:.4e}");
    eprintln!("  B: E_emf = {b_emf:.4e}, E_R = {b_r:.4e}, E_ap(v·i) = {b_ap:.4e}");

    // The passive branch never generates: v·i = i²(R+β) ≥ 0 sample-wise.
    for &(_, v, i) in &records[1] {
        assert!(
            v * i >= -1e-18,
            "engine-power-001 FAILED: passive port produced power (v·i = {:.3e})",
            v * i
        );
    }
    assert!(b_emf.abs() < 1e-30, "passive port has no EMF");

    // The conservation identity on circuit-side quantities only: the EMF's
    // supply is dissipated in the two branch resistors plus what the field
    // keeps/loses (CPML leakage, residual ring) — so the resistor share
    // must be a large-but-<1 fraction. P_accepted-by-field for FS.2b's
    // gain is E_emf − E_R(A).
    let closure = (a_r + b_r) / a_emf;
    let p_acc = a_emf - a_r;
    eprintln!(
        "  closure (E_R(A)+E_R(B))/E_emf = {closure:.4}; accepted-by-field = {p_acc:.4e} J          ({:.1} % of EMF supply)",
        100.0 * p_acc / a_emf
    );

    assert!(
        a_emf > 0.0 && a_emf.is_finite(),
        "engine-power-001 FAILED: no EMF supply measured"
    );
    // Measured-then-pinned: the first identity-closing run read 0.9917
    // (0.8 % to CPML leakage + residual ring on 7000 steps).
    assert!(
        (0.95..=1.0).contains(&closure),
        "engine-power-001 FAILED: resistor closure {closure:.4} outside [0.95, 1.0] \
         (measured 0.9917 at ship time) — energy bookkeeping does not close"
    );
}
