//! Validation gate **fdtd-coupling-001** — FDTD-extracted even/odd modal
//! resonant split of a coupled microstrip-resonator pair, cross-checked against
//! the analytic **εeff-split** of the shipped Kirschning-Jansen coupled-line
//! model (Filter Phase F1.1b.1, ADR-0108).
//!
//! # What this gate proves
//!
//! It is the first *full-wave* EM solve in the filter-design pipeline. Every
//! step before it is closed-form (synthesis → dimensional synthesis → layout →
//! manufacturing files); this gate confirms that a *dimensioned coupled pair*
//! actually realizes its even/odd modal split by running an FDTD simulation
//! ([`yee_voxel::run_coupled_pair`]) with even/odd modal excitation and
//! comparing the extracted `k = (f_odd² − f_even²)/(f_odd² + f_even²)` to the
//! analytic εeff-split of the same geometry.
//!
//! # Geometry
//!
//! Two parallel, edge-coupled half-wave microstrip resonators on FR-4
//! (`εr = 4.4`, `h = 1.6 mm`), each of width `W` and length `L ≈ λ_g/2` at the
//! synchronous centre `f0 = 2.4 GHz`, separated by an edge-to-edge gap `S`.
//! The analytic even/odd effective permittivities come from
//! [`yee_layout::coupled_microstrip`] `(W, S, h, εr)`.
//!
//! The half-wave length is `L = c / (2·f0·√εeff)`. Using the even-mode
//! effective permittivity as a representative `εeff` keeps the resonator near
//! the synchronous frequency the driver scans around. The exact length is not
//! load-bearing for the *coupling* comparison (the driver scans a ±35 % window
//! and reads each mode's dominant peak, not the absolute resonance), so a
//! closed-form estimate is used directly — no FDTD tuning (ADR-0108).
//!
//! # Reference: the εeff-split, NOT the impedance coupling
//!
//! This geometry is two *full-length* coupled λ/2 resonators (coupled over
//! their entire length). Their even/odd resonant frequencies split by the
//! even/odd *phase-velocity* difference — `√εeff,e` vs `√εeff,o` — so the
//! physically-correct frequency-split reference is the **εeff-split**
//!
//! ```text
//! k_ref = (εeff_e − εeff_o) / (εeff_e + εeff_o).
//! ```
//!
//! The PR #1 root-cause analysis (5 CI iterations) established that the prior
//! reference [`yee_layout::coupling_coefficient`] =
//! `(z0e − z0o)/(z0e + z0o)` is the **impedance** coupling of a λ/4-overlap
//! parallel-coupled-line BPF *section* — the wrong quantity for the resonant
//! split of two full-length coupled lines. The FDTD even/odd modal split is
//! `(f_odd² − f_even²)/(f_odd² + f_even²)`, which is an εeff-split by
//! construction, so the apples-to-apples reference is `k_ref` above.
//!
//! # Tolerance
//!
//! The gate asserts the FDTD `k` matches the analytic εeff-split within
//! **≤ 15 % relative** — a deliberately loose walking-skeleton band: coarse-grid
//! FDTD is approximate (under-resolved εeff, finite air box) and this is the
//! first end-to-end run. Tightening is a follow-on once the skeleton is green.
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

use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, coupled_microstrip};
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

/// `fdtd-coupling-001`: FDTD modal-split `k` vs analytic εeff-split, ≤ 15 % rel.
#[test]
#[ignore = "slow: multi-minute FDTD; fdtd-coupling-001 coupled-resonator k gate (F1.1b.1, ADR-0108); run with --release --ignored"]
fn fdtd_coupling_001_matches_analytic_within_fifteen_percent() {
    // ------------------------------------------------------------------
    // Analytic reference: the εeff-split of the Kirschning-Jansen even/odd
    // model. For two FULL-LENGTH coupled λ/2 resonators the even/odd resonant
    // frequencies split by the even/odd phase-velocity difference (√εeff,e vs
    // √εeff,o), so the physically-correct frequency-split reference is
    //   k_ref = (εeff_e − εeff_o)/(εeff_e + εeff_o),
    // NOT the impedance coupling (z0e − z0o)/(z0e + z0o) — which is the
    // λ/4-overlap coupled-line-SECTION quantity and the wrong reference for this
    // geometry (PR #1 5-iteration root-cause analysis, ADR-0108).
    // ------------------------------------------------------------------
    let model = coupled_microstrip(W_M, S_M, H_M, EPS_R);
    let k_ref = (model.eps_eff_e - model.eps_eff_o) / (model.eps_eff_e + model.eps_eff_o);

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
    // Run the FDTD coupled-resonator driver with even/odd modal excitation
    // (synchronous centre = F0_HZ). The driver does two sub-runs: in-phase
    // (even supermode) and anti-phase (odd supermode), each giving a single
    // dominant resonance, then forms k = (f_odd² − f_even²)/(f_odd² + f_even²).
    // ------------------------------------------------------------------
    // A large air box is retained from the PR #1 iter#5 finding: hard-PEC
    // outer walls (WalkingSkeletonSolver::new) confine the microstrip fringing
    // + air-gap fields that SET the even/odd εeff difference, suppressing the
    // split — a box-size (mm) effect, hence grid-independent. dx = 0.4 mm with
    // air_above 48 / xy margin 24 cells affords a much larger open region.
    let cfg = CoupledRunConfig {
        f0_hz: F0_HZ,
        dx_m: 0.4e-3,
        xy_margin_cells: 24,
        air_above_cells: 48,
        n_steps: 40_000,
        ..CoupledRunConfig::default()
    };
    let result = run_coupled_pair(&layout, &cfg);

    let rel_err = (result.k - k_ref).abs() / k_ref.abs();

    eprintln!(
        "\nfdtd-coupling-001 coupled-resonator k gate (F1.1b.1, ADR-0108)
  geometry:      W = {:.3} mm, S = {:.3} mm, h = {:.3} mm, εr = {EPS_R}
  half-wave L:   {:.3} mm  (εeff,e = {:.4}, f0 = {:.3} GHz)
  analytic:      εeff,e = {:.4}, εeff,o = {:.4}  ->  εeff-split k_ref = {:.5}
  FDTD modal:    f_even = {:.5} GHz, f_odd = {:.5} GHz  ->  k = {:.5}
  relative err:  {:.3} %  (threshold ≤ 15 %)
",
        W_M * 1e3,
        S_M * 1e3,
        H_M * 1e3,
        l_m * 1e3,
        eps_eff,
        F0_HZ * 1e-9,
        model.eps_eff_e,
        model.eps_eff_o,
        k_ref,
        result.f_even * 1e-9,
        result.f_odd * 1e-9,
        result.k,
        rel_err * 100.0,
    );

    assert!(
        rel_err <= 0.15,
        "fdtd-coupling-001 FAILED: FDTD modal-split k = {:.5}, analytic \
         εeff-split k_ref = {:.5}, relative error = {:.3} % (threshold ≤ 15 %)",
        result.k,
        k_ref,
        rel_err * 100.0,
    );
}
