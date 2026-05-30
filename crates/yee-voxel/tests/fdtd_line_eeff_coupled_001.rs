//! Validation gate **fdtd-line-eeff-coupled-001** — FDTD-extracted even/odd
//! effective permittivities of a *coupled* microstrip pair, from a
//! **propagation** (phase-velocity) measurement, cross-checked against the
//! analytic Kirschning-Jansen coupled-line `ε_eff,e` / `ε_eff,o` (Filter Phase
//! F1.1b.1 coupled follow-on, ADR-0108).
//!
//! # What this gate proves
//!
//! It is the coupled-line extension of `fdtd-line-eeff-001`: the first
//! full-wave EM solve of a *coupled* pair in the filter-design pipeline. It
//! drives two parallel edge-coupled strips into their even and odd supermodes
//! (in-phase / anti-phase) and measures each mode's phase velocity →
//! `ε_eff,e` / `ε_eff,o`, comparing to the analytic even/odd effective
//! permittivities of the same geometry. This is the quantity a coupled-resonator
//! filter's coupling derives from.
//!
//! # Why propagation, not a resonant split
//!
//! Same reason as `fdtd-line-eeff-001` (PR #1, 7 iterations): a resonant split
//! is unworkable for an open microstrip (PEC box confines or resonates; open
//! CPML kills the Q). A driven, time-gated propagation measurement on a long
//! line is robust — no Q, no cavity modes, no peak-picking.
//!
//! # Geometry
//!
//! Two parallel edge-coupled microstrip lines on FR-4 (`εr = 4.4`,
//! `h = 1.6 mm`), each of width `W = 3 mm`, separated by an edge-to-edge gap
//! `S = 1 mm`, several guided wavelengths long at the 5 GHz drive. The analytic
//! even/odd effective permittivities come from
//! [`yee_layout::coupled_microstrip`] `(W, S, h, εr)`.
//!
//! # Tolerance
//!
//! The gate asserts the FDTD `ε_eff,e` and `ε_eff,o` each match the analytic
//! value within **≤ 15 % relative** — the same loose walking-skeleton band as
//! the single-line gate.
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-voxel --release -- --ignored fdtd_line_eeff_coupled_001 --nocapture
//! ```

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, coupled_microstrip};
use yee_voxel::{LineRunConfig, run_coupled_line_eeff};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const S_M: f64 = 1.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;

/// Yee cell size (metres). 0.3 mm → ~5 cells through the FR-4 substrate.
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;

/// `fdtd-line-eeff-coupled-001`: FDTD propagation-extracted even/odd `ε_eff` vs
/// analytic coupled-line `ε_eff,e` / `ε_eff,o`, each ≤ 15 % rel.
#[test]
#[ignore = "slow: multi-minute FDTD (two solves); fdtd-line-eeff-coupled-001 (F1.1b.1, ADR-0108); run with --release --ignored"]
fn fdtd_line_eeff_coupled_001_matches_analytic_within_fifteen_percent() {
    // Analytic even/odd effective permittivities (Kirschning-Jansen).
    let model = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    let eps_e_ref = model.eps_eff_e;
    let eps_o_ref = model.eps_eff_o;
    let k_ref = (eps_e_ref - eps_o_ref) / (eps_e_ref + eps_o_ref);

    // Length ≈ 6 λ_g (use the even-mode εeff for λ_g — the longer of the two).
    let lam_g = C0_M_S / (F0_HZ * eps_e_ref.sqrt());
    let l_m = 6.0 * lam_g;

    // Strip 1 spans x = [0, L] at y = [0, W]; strip 2 is offset by W + S in y.
    let y1 = 0.0;
    let y2 = W_M + S_M;
    let strip1 = Polygon::rect(0.0, y1, l_m, W_M);
    let strip2 = Polygon::rect(0.0, y2, l_m, W_M);
    let traces = vec![strip1, strip2];
    let bbox = BBox::from_polygons(&traces);

    // One drive port at the −x end of each strip (centre of its width).
    let ports = vec![
        PortRef {
            at: Point2::new(0.5e-3, y1 + W_M / 2.0),
            width_m: W_M,
            ref_impedance_ohm: 50.0,
        },
        PortRef {
            at: Point2::new(0.5e-3, y2 + W_M / 2.0),
            width_m: W_M,
            ref_impedance_ohm: 50.0,
        },
    ];

    let layout = Layout {
        substrate: Substrate {
            eps_r: EPS_R,
            height_m: H_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports,
        bbox,
    };

    // Upstream probe A at 2.5 λ_g; downstream probe B a third of a wavelength on.
    let x_a = 2.5 * lam_g;
    let probe_b_offset = lam_g / 3.0;
    let x_b = x_a + probe_b_offset;

    // Time-gate to the forward pulse, before the far-PEC reflection returns to B.
    let v_p = C0_M_S / eps_e_ref.sqrt();
    let dt = 0.9 * DX_M / (C0_M_S * 3.0_f64.sqrt());
    let x_drive = 0.5e-3;
    let t_refl_b = ((l_m - x_drive) + (l_m - x_b)) / v_p;
    let gate_steps = (0.9 * t_refl_b / dt) as usize;

    let cfg = LineRunConfig {
        f0_hz: F0_HZ,
        dx_m: DX_M,
        xy_margin_cells: MARGIN_CELLS,
        air_above_cells: AIR_ABOVE_CELLS,
        probe_a_x_m: x_a,
        probe_b_x_m: x_b,
        gate_steps: Some(gate_steps),
        n_steps: gate_steps + 200,
        ..LineRunConfig::default()
    };

    let result = run_coupled_line_eeff(&layout, &cfg, probe_b_offset);

    let err_e = (result.eps_eff_e - eps_e_ref).abs() / eps_e_ref.abs();
    let err_o = (result.eps_eff_o - eps_o_ref).abs() / eps_o_ref.abs();

    eprintln!(
        "\nfdtd-line-eeff-coupled-001 even/odd ε_eff gate (F1.1b.1, ADR-0108)
  geometry:      W = {:.3} mm, S = {:.3} mm, h = {:.3} mm, εr = {EPS_R}, L = {:.2} mm
  drive:         f0 = {:.3} GHz, λ_g = {:.2} mm
  analytic:      ε_eff,e = {:.4}, ε_eff,o = {:.4}  (split k_ref = {:.4})
  FDTD:          ε_eff,e = {:.4}, ε_eff,o = {:.4}  (split k     = {:.4})
  rel err:       even {:.3} %, odd {:.3} %  (threshold ≤ 15 %)
",
        W_M * 1e3,
        S_M * 1e3,
        H_M * 1e3,
        l_m * 1e3,
        F0_HZ * 1e-9,
        lam_g * 1e3,
        eps_e_ref,
        eps_o_ref,
        k_ref,
        result.eps_eff_e,
        result.eps_eff_o,
        result.k_split,
        err_e * 100.0,
        err_o * 100.0,
    );

    assert!(
        err_e <= 0.15 && err_o <= 0.15,
        "fdtd-line-eeff-coupled-001 FAILED: even ε_eff = {:.4} (ref {:.4}, err {:.3} %), \
         odd ε_eff = {:.4} (ref {:.4}, err {:.3} %); threshold ≤ 15 %",
        result.eps_eff_e,
        eps_e_ref,
        err_e * 100.0,
        result.eps_eff_o,
        eps_o_ref,
        err_o * 100.0,
    );
}
