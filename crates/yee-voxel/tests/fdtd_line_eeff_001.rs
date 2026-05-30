//! Validation gate **fdtd-line-eeff-001** — FDTD-extracted effective
//! permittivity `ε_eff` of a single straight microstrip line, from a
//! **propagation** (phase-velocity) measurement, cross-checked against the
//! analytic Hammerstad-Jensen / Pozar single-line `ε_eff` (Filter Phase
//! F1.1b.1, ADR-0108).
//!
//! # What this gate proves
//!
//! It is the first *full-wave* EM solve in the filter-design pipeline. Every
//! step before it is closed-form (synthesis → dimensional synthesis → layout →
//! manufacturing files); this gate confirms that a *dimensioned microstrip
//! line* actually propagates a wave at the right phase velocity by running an
//! FDTD simulation ([`yee_voxel::run_line_eeff`]) and comparing the extracted
//! `ε_eff = (c / v_p)²` to the analytic single-line value. A validated
//! full-wave `ε_eff` of a real dimensioned line *is* "full filter simulation
//! works" at the walking-skeleton level.
//!
//! # Why a propagation measurement (not a resonant split)
//!
//! PR #1 (7 CI/local iterations) proved the resonant-split method unworkable
//! for a microstrip on an open domain: a small hard-PEC box confines the
//! fringing fields (wrong `ε_eff`), a large PEC box becomes a resonant cavity
//! (box modes swamp the spectrum), and open CPML walls kill the resonator Q
//! (zero peaks). There is no box that is simultaneously high-Q and
//! non-confining/non-resonant. The robust, textbook alternative is a
//! propagation measurement: drive a long *non-resonant* line, terminate both
//! ends in matched CPML loads, and read the phase velocity off two probe planes
//! a known distance apart — no Q, no cavity modes, no peak-picking.
//!
//! # Geometry
//!
//! One straight microstrip line on FR-4 (`εr = 4.4`, `h = 1.6 mm`) of width
//! `W = 3 mm` (the FR-4 50 Ω width) and a length of several guided wavelengths
//! at the 5 GHz drive frequency. The analytic single-line `ε_eff` comes from
//! [`yee_layout::eps_eff`] `(W, h, εr)`.
//!
//! # Tolerance
//!
//! The gate asserts the FDTD `ε_eff` matches the analytic value within
//! **≤ 15 % relative** — a deliberately loose walking-skeleton band: coarse-grid
//! staircased FDTD is approximate (a 0.4 mm cell stairsteps the strip edges and
//! under-resolves the substrate). Tightening is a follow-on once the skeleton is
//! green.
//!
//! # Why `#[ignore]`'d + CI-routed
//!
//! The FDTD run is multi-second to a minute and the dev box is
//! memory-/CPU-constrained, so this gate never runs in the default
//! `cargo test`. It runs in a dedicated CI `--release` job
//! (`fdtd-coupling-gate` in `.github/workflows/ci.yml`), mirroring the
//! `mom-001` / GPU-nightly `--release -- --ignored` idiom. Per CLAUDE.md §4 the
//! feature must not merge until this gate is GREEN somewhere.
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-voxel --release -- --ignored fdtd_line_eeff_001 --nocapture
//! ```

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{LineRunConfig, run_line_eeff};

/// FR-4 substrate.
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;

/// Strip width `W` — the FR-4 50 Ω width (≈ 3 mm).
const W_M: f64 = 3.0e-3;

/// Drive centre frequency (Hz). 5 GHz → short guided wavelength → fewer cells
/// per λ_g → a cheap grid for a given probe spacing.
const F0_HZ: f64 = 5.0e9;

const C0_M_S: f64 = 299_792_458.0;

/// Yee cell size (metres). 0.3 mm → ~5 cells through the 1.6 mm FR-4 substrate.
const DX_M: f64 = 0.3e-3;

/// Lateral air margin and air-above clearance, in cells. Microstrip fringing
/// fields extend ~several substrate heights into the air; too-close PEC box
/// walls confine them and bias the measured ε_eff. ~10 mm of clearance each
/// way (≈ 6 h) keeps the box walls clear of the fringing region.
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;

/// `fdtd-line-eeff-001`: FDTD propagation-extracted `ε_eff` vs analytic
/// single-line `ε_eff`, ≤ 15 % rel.
#[test]
#[ignore = "slow: multi-second FDTD; fdtd-line-eeff-001 propagation ε_eff gate (F1.1b.1, ADR-0108); run with --release --ignored"]
fn fdtd_line_eeff_001_matches_analytic_within_fifteen_percent() {
    // ------------------------------------------------------------------
    // Analytic reference: single-line ε_eff (Pozar §3.8 closed form).
    // ------------------------------------------------------------------
    let eps_eff_ref = eps_eff(W_M, H_M, EPS_R);

    // ------------------------------------------------------------------
    // Build a single straight microstrip line several guided wavelengths long
    // at the drive frequency, fed at one end. Length L ≈ 6·λ_g leaves a long
    // run downstream of the probes so the reflection off the far hard-PEC wall
    // returns *after* the forward pulse has cleared the probes (the time gate
    // below cuts the DFT off before then). A long line + a short look is the
    // standard reflection-free FDTD line characterization.
    // ------------------------------------------------------------------
    let lam_g = C0_M_S / (F0_HZ * eps_eff_ref.sqrt());
    let l_m = 6.0 * lam_g;

    // Strip spans x = [0, L] at y = [0, W].
    let trace = Polygon::rect(0.0, 0.0, l_m, W_M);
    let traces = vec![trace];
    let bbox = BBox::from_polygons(&traces);

    // One drive port near the −x end, centred on the strip width.
    let ports = vec![PortRef {
        at: Point2::new(0.5e-3, W_M / 2.0),
        width_m: W_M,
        ref_impedance_ohm: 50.0,
    }];

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

    // ------------------------------------------------------------------
    // Run the propagation-based ε_eff driver: drive one end, CPML-terminate
    // both ends, read the phase velocity off two probe planes (open domain,
    // matched terminations — the determined PR #1 / ADR-0108 method).
    // ------------------------------------------------------------------
    // Probe planes near the middle of the line: A at 2.5 λ_g (well past the
    // feed transient), B a third of a wavelength downstream → Δx ≈ λ_g/3, a
    // ~120° phase advance: unambiguous (< 2π) and well-resolved.
    let x_a = 2.5 * lam_g;
    let x_b = x_a + lam_g / 3.0;

    // Time-gate the DFT to the forward pulse, before the far-PEC reflection
    // returns to the downstream probe B. The reflection path is
    // drive → far end (≈ L) → back to B, i.e. distance (L − x_drive) + (L − x_B)
    // at the phase velocity v_p = c/√ε_eff. Stop ~10 % short of that, as a
    // margin against the finite pulse width. The forward pulse (≈ 800 steps at
    // 80 % bandwidth) reaches B in a few hundred steps, so the gate still
    // captures its full passage.
    let v_p = C0_M_S / eps_eff_ref.sqrt();
    // dt is 0.9× the cubic Courant limit at the grid `dx`; recompute it here to
    // size the gate in steps.
    let dx = DX_M;
    let dt = 0.9 * dx / (C0_M_S * 3.0_f64.sqrt());
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
        // Run a little past the gate so the forward pulse is fully integrated.
        n_steps: gate_steps + 200,
        ..LineRunConfig::default()
    };
    let result = run_line_eeff(&layout, &cfg);

    let rel_err = (result.eps_eff - eps_eff_ref).abs() / eps_eff_ref.abs();

    eprintln!(
        "\nfdtd-line-eeff-001 propagation ε_eff gate (F1.1b.1, ADR-0108)
  geometry:      W = {:.3} mm, h = {:.3} mm, εr = {EPS_R}, L = {:.2} mm ({:.2} λ_g)
  drive:         f0 = {:.3} GHz, λ_g = {:.2} mm
  probes:        Δx = {:.3} mm, Δφ = {:.4} rad, v_p = {:.4e} m/s
  analytic:      ε_eff = {:.4}
  FDTD:          ε_eff = {:.4}
  relative err:  {:.3} %  (threshold ≤ 15 %)
",
        W_M * 1e3,
        H_M * 1e3,
        l_m * 1e3,
        l_m / lam_g,
        F0_HZ * 1e-9,
        lam_g * 1e3,
        result.delta_x * 1e3,
        result.delta_phi,
        result.v_p,
        eps_eff_ref,
        result.eps_eff,
        rel_err * 100.0,
    );

    assert!(
        rel_err <= 0.15,
        "fdtd-line-eeff-001 FAILED: FDTD ε_eff = {:.4}, analytic ε_eff = {:.4}, \
         relative error = {:.3} % (threshold ≤ 15 %)",
        result.eps_eff,
        eps_eff_ref,
        rel_err * 100.0,
    );
}
