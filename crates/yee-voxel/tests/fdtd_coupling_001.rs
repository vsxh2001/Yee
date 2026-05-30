//! Validation gate **fdtd-coupling-001** — FDTD-extracted inter-resonator
//! coupling `k` of a coupled microstrip-resonator pair, cross-checked against
//! the shipped analytic Kirschning-Jansen coupled-line reference
//! (Filter Phase F1.1b.1, ADR-0108).
//!
//! # What this gate proves
//!
//! It is the first *full-wave* EM solve in the filter-design pipeline. Every
//! step before it is closed-form (synthesis → dimensional synthesis → layout →
//! manufacturing files); this gate confirms that a *dimensioned coupled pair*
//! actually realizes its *target coupling* by running an FDTD simulation
//! ([`yee_voxel::run_coupled_pair`]) and comparing the extracted `k` to the
//! analytic [`yee_layout::coupling_coefficient`] of the same geometry.
//!
//! # Geometry
//!
//! Two parallel, edge-coupled half-wave microstrip resonators on FR-4
//! (`εr = 4.4`, `h = 1.6 mm`), each of width `W` and length `L ≈ λ_g/2` at the
//! synchronous centre `f0 = 2.4 GHz`, separated by an edge-to-edge gap `S`.
//! The analytic even/odd impedances come from [`yee_layout::coupled_microstrip`]
//! `(W, S, h, εr)`, and the analytic coupling from
//! [`yee_layout::coupling_coefficient`].
//!
//! The half-wave length is `L = c / (2·f0·√εeff)`. Using the even-mode
//! effective permittivity as a representative `εeff` keeps the resonator near
//! the synchronous frequency the driver scans around. The exact length is not
//! load-bearing for the *coupling* comparison (the driver scans a ±35 % window
//! and reads the split, not the absolute resonance), so a closed-form estimate
//! is used directly — no FDTD tuning (ADR-0108).
//!
//! # Tolerance
//!
//! The gate asserts the FDTD `k` matches the analytic `k` within **≤ 15 %
//! relative** — a deliberately loose walking-skeleton band: coarse-grid FDTD is
//! approximate, the two extraction models differ (split-peak inversion vs
//! coupled-line `(Z₀ₑ − Z₀ₒ)/(Z₀ₑ + Z₀ₒ)`), and this is the first end-to-end
//! run. Tightening is a follow-on once the skeleton is green.
//!
//! # Why `#[ignore]`'d + CI-routed
//!
//! The FDTD run is multi-minute and the dev box is memory-/CPU-constrained, so
//! this gate never runs in the default `cargo test`. It runs in a dedicated CI
//! `--release` job (`fdtd-coupling-gate` in `.github/workflows/ci.yml`),
//! mirroring the `mom-001` / GPU-nightly `--release -- --ignored` idiom. Per
//! CLAUDE.md §4 the feature must not merge until this gate is GREEN somewhere.
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-voxel --release -- --ignored fdtd_coupling_001 --nocapture
//! ```

use yee_layout::{
    BBox, Layout, Point2, Polygon, PortRef, Substrate, coupled_microstrip, coupling_coefficient,
};
use yee_voxel::{CoupledRunConfig, run_coupled_pair};

const C0_M_S: f64 = 299_792_458.0;

/// FR-4 substrate.
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;

/// Strip width `W` and edge-to-edge coupling gap `S` (metres). A moderate gap
/// gives a measurable-but-not-extreme split; `W ≈ 3 mm` is the FR-4 50 Ω width.
const W_M: f64 = 3.0e-3;
const S_M: f64 = 1.0e-3;

/// Synchronous resonator centre frequency (Hz).
const F0_HZ: f64 = 2.4e9;

/// `fdtd-coupling-001`: FDTD `k` vs analytic coupled-line `k`, ≤ 15 % relative.
#[test]
#[ignore = "slow: multi-minute FDTD; fdtd-coupling-001 coupled-resonator k gate (F1.1b.1, ADR-0108); run with --release --ignored"]
fn fdtd_coupling_001_matches_analytic_within_fifteen_percent() {
    // ------------------------------------------------------------------
    // Analytic reference: Kirschning-Jansen even/odd model -> coupling k.
    // ------------------------------------------------------------------
    let model = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    let k_analytic = coupling_coefficient(&model);

    // ------------------------------------------------------------------
    // Build the coupled-pair Layout: two parallel half-wave resonators
    // separated by the gap S, each fed at one end.
    //
    // Half-wave length L = c / (2·f0·√εeff). Use the even-mode εeff as a
    // representative effective permittivity (the driver scans a wide window
    // around f0 and reads the split, not the absolute resonance).
    // ------------------------------------------------------------------
    let eps_eff = model.eps_eff_e;
    let l_m = C0_M_S / (2.0 * F0_HZ * eps_eff.sqrt());

    // Strip 1 spans x = [0, L] at y = [0, W]; strip 2 is offset by W + S in y.
    let y1 = 0.0;
    let y2 = W_M + S_M;
    let strip1 = Polygon::rect(0.0, y1, l_m, W_M);
    let strip2 = Polygon::rect(0.0, y2, l_m, W_M);
    let traces = vec![strip1, strip2];
    let bbox = BBox::from_polygons(&traces);

    // One port near the fed end of each resonator (centre of the strip width).
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

    // ------------------------------------------------------------------
    // Run the FDTD coupled-resonator driver (walking-skeleton defaults,
    // synchronous centre = F0_HZ).
    // ------------------------------------------------------------------
    let cfg = CoupledRunConfig {
        f0_hz: F0_HZ,
        ..CoupledRunConfig::default()
    };
    let result = run_coupled_pair(&layout, &cfg);

    let rel_err = (result.k - k_analytic).abs() / k_analytic.abs();

    eprintln!(
        "\nfdtd-coupling-001 coupled-resonator k gate (F1.1b.1, ADR-0108)
  geometry:      W = {:.3} mm, S = {:.3} mm, h = {:.3} mm, εr = {EPS_R}
  half-wave L:   {:.3} mm  (εeff,e = {:.4}, f0 = {:.3} GHz)
  analytic:      Z0e = {:.2} Ω, Z0o = {:.2} Ω  ->  k = {:.5}
  FDTD:          f_even = {:.5} GHz, f_odd = {:.5} GHz  ->  k = {:.5}
  relative err:  {:.3} %  (threshold ≤ 15 %)
",
        W_M * 1e3,
        S_M * 1e3,
        H_M * 1e3,
        l_m * 1e3,
        eps_eff,
        F0_HZ * 1e-9,
        model.z0e_ohm,
        model.z0o_ohm,
        k_analytic,
        result.f_even * 1e-9,
        result.f_odd * 1e-9,
        result.k,
        rel_err * 100.0,
    );

    assert!(
        rel_err <= 0.15,
        "fdtd-coupling-001 FAILED: FDTD k = {:.5}, analytic k = {:.5}, \
         relative error = {:.3} % (threshold ≤ 15 %)",
        result.k,
        k_analytic,
        rel_err * 100.0,
    );
}
