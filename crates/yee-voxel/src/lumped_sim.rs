//! Lumped-LC **FDTD EM simulation** of a synthesized filter board (Filter
//! Phase F2.3, ADR-0115).
//!
//! [`simulate_lumped_board`] takes a synthesized [`yee_filter::LumpedLadder`]
//! (F2.0), places it as SMD footprints on a microstrip board with
//! [`yee_filter::lumped_board`] (F2.2), voxelizes that board onto a Yee grid
//! with [`crate::voxelize_microstrip`] (F1.1a), drops a [`LumpedRlcPort`] at
//! each placed L/C, drives one port / matches the other, time-steps the FDTD,
//! and extracts the forward transmission `|S21|(f)` over a frequency sweep —
//! the lumped analogue of the distributed [`crate::run_line_eeff`] gate.
//!
//! It is the goal's named "EM simulation" of the lumped track, cross-validated
//! against the analytic circuit response [`yee_filter::ladder_s21`] (the
//! `fdtd_lumped_001` gate).
//!
//! # Series / shunt decomposition (the key subtlety)
//!
//! [`LumpedRlcPort::series_rlc`] is a **series** R-L-C at one `E_z` edge. The
//! ladder's two resonator topologies map onto it differently:
//!
//! - A **series-branch** resonator (a series L–C in the through arm) is exactly
//!   one [`LumpedRlcPort::series_rlc`]`(cell, r, L, C)` placed on the in-line
//!   gap cell of that footprint: its impedance `Z = jωL + 1/(jωC)` is in series
//!   with the signal path, as the ABCD model wants.
//! - A **shunt-branch** resonator (a *parallel* L‖C from the line to ground)
//!   cannot be one series element, so it becomes **two** ports at the *same*
//!   shunt cell — a pure inductor [`LumpedRlcPort::series_rlc`]`(cell, r, L,
//!   ∞)` (`c = ∞` shorts the capacitor) **in parallel with** a pure capacitor
//!   [`LumpedRlcPort::series_rlc`]`(cell, r, 0, C)` (`l = 0` removes the
//!   inductor). Both correct the *same* `E_z` edge each step, so the lattice
//!   sees their currents summed → a parallel `L‖C` admittance
//!   `Y = jωC + 1/(jωL)`, the correct shunt-resonator topology.
//!
//! A small finite series resistance ([`SERIES_ESR_OHM`]) is used on every
//! element because [`LumpedRlcPort::series_rlc`] requires `r > 0`.
//!
//! # S21 extraction (thru-normalized)
//!
//! The single drive port (a modulated-Gaussian series EMF through a `Z0`
//! resistor) launches a pulse down the signal line; the far port is a matched
//! [`LumpedRlcPort::pure_resistor`]`(Z0)` load. The voltage `V = E_z · dz`
//! sensed at the load cell is single-bin-DFT'd at each sweep frequency. To
//! divide out the (frequency-dependent, coarse-grid-dependent) feed + line +
//! port coupling, the **same** board is run a second time with the filter
//! elements removed (a bare through line) and its load voltage `V_thru(f)`
//! recorded; then
//!
//! ```text
//! S21(f) = V_dut(f) / V_thru(f)
//! ```
//!
//! which is ≈ 1 (0 dB) in the passband and rolls off in the stopband — the
//! transmission *relative to the matched thru*, exactly the quantity
//! [`yee_filter::ladder_s21`] computes for the ideal circuit. This thru
//! calibration is robust against the lumped-port `E_z` voltage convention and
//! the one-way circuit→field coupling (see [`yee_fdtd::LumpedRlcPort`]).

use std::f64::consts::PI;

use yee_fdtd::{LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver};
use yee_filter::{BranchKind, Footprint, LcBranch, LumpedLadder, Placement, lumped_board};
use yee_layout::Substrate;

use crate::{VoxelOptions, voxelize_microstrip};

/// Series ESR (Ω) used on every lumped element.
///
/// [`LumpedRlcPort::series_rlc`] requires `r > 0`; a tiny value approximates an
/// ideal (lossless) L/C while keeping the constructor happy. Matches the
/// `r ≈ 1e-3` the F2.3 spec calls for.
pub const SERIES_ESR_OHM: f64 = 1.0e-3;

/// Configuration for [`simulate_lumped_board`].
///
/// The defaults target a ~2 GHz lumped band-pass filter on FR-4: a coarse but
/// tractable cubic grid, a modulated-Gaussian drive centred in-band, enough
/// time steps for the pulse to transit the board, and a frequency sweep that
/// spans the passband and the stopband cross-check point. Heavy (multi-minute
/// FDTD) → iterate in the bounded dev container (CLAUDE.md §10).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LumpedSimConfig {
    /// Isotropic Yee cell size `dx = dy = dz`, metres.
    pub dx_m: f64,
    /// Air margin (in cells) around the board bounding box in `x`/`y`
    /// (forwarded to [`VoxelOptions::xy_margin_cells`]).
    pub xy_margin_cells: usize,
    /// Air layers above the top metal (forwarded to
    /// [`VoxelOptions::air_above_cells`]).
    pub air_above_cells: usize,
    /// SMD footprint the board is built with (forwarded to
    /// [`yee_filter::lumped_board`]).
    pub footprint: Footprint,
    /// Total number of FDTD time steps per solve. Must be long enough for the
    /// launched pulse to fully transit the board and be integrated by the DFT.
    pub n_steps: usize,
    /// Drive centre frequency `f0` (Hz) of the modulated-Gaussian pulse. Set
    /// near the filter centre so the launched spectrum covers the passband.
    pub drive_f0_hz: f64,
    /// Fractional FWHM bandwidth of the modulated-Gaussian drive, as a fraction
    /// of `drive_f0_hz` (`bw = drive_bw_frac · drive_f0_hz`). Wide enough to
    /// cover the passband *and* the stopband cross-check point in one launch.
    pub drive_bw_frac: f64,
    /// Peak drive voltage (V) of the Gaussian-modulated pulse.
    pub drive_v0: f64,
    /// Sweep lower frequency bound (Hz).
    pub f_lo_hz: f64,
    /// Sweep upper frequency bound (Hz).
    pub f_hi_hz: f64,
    /// Number of (linearly spaced) sweep points in `[f_lo_hz, f_hi_hz]`.
    pub n_freq: usize,
}

impl Default for LumpedSimConfig {
    /// Walking-skeleton defaults for a ~2 GHz lumped band-pass on FR-4. See the
    /// per-field docs; container-validated to pass `fdtd_lumped_001` (ADR-0115).
    fn default() -> Self {
        Self {
            dx_m: 0.4e-3,
            xy_margin_cells: 8,
            air_above_cells: 8,
            footprint: Footprint::Smd0603,
            n_steps: 4_000,
            drive_f0_hz: 2.0e9,
            drive_bw_frac: 1.2,
            drive_v0: 1.0,
            f_lo_hz: 1.0e9,
            f_hi_hz: 3.0e9,
            n_freq: 41,
        }
    }
}

/// Single-bin DFT accumulator at one frequency: real/imag of `Σ x[n]·e^{-jωt}`.
#[derive(Clone, Copy)]
struct Bin {
    omega: f64,
    re: f64,
    im: f64,
}

impl Bin {
    fn new(omega: f64) -> Self {
        Self {
            omega,
            re: 0.0,
            im: 0.0,
        }
    }
    fn accumulate(&mut self, x: f64, t: f64) {
        let phase = self.omega * t;
        self.re += x * phase.cos();
        self.im -= x * phase.sin();
    }
    fn mag(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

/// EM-simulate a synthesized lumped-LC filter board and return its forward
/// transmission `|S21|(f)` over a frequency sweep.
///
/// Builds the board ([`yee_filter::lumped_board`]), voxelizes it
/// ([`crate::voxelize_microstrip`]), places each ladder element as
/// [`LumpedRlcPort`](s) on the grid (series branch → one series-RLC; shunt
/// branch → pure-L ‖ pure-C, see the [module docs](self)), drives the input
/// port with a modulated-Gaussian source and matches the output port with a
/// `Z0` resistor, time-steps the FDTD, and single-bin-DFTs the load voltage.
/// A second, element-free *thru* solve normalizes the response, so the returned
/// `|S21|` is the transmission relative to the matched thru (≈ 0 dB in-band).
///
/// Returns `(freq_hz, |S21|)` for each of [`LumpedSimConfig::n_freq`] linearly
/// spaced points in `[f_lo_hz, f_hi_hz]`.
///
/// # Panics
///
/// Panics if `ladder.resonators` is empty, if the board does not voxelize to at
/// least two ports (a drive + a load), or if [`LumpedSimConfig`] has a
/// non-positive `dx_m` / `n_freq`.
pub fn simulate_lumped_board(
    ladder: &LumpedLadder,
    substrate: &Substrate,
    cfg: &LumpedSimConfig,
) -> Vec<(f64, f64)> {
    assert!(
        !ladder.resonators.is_empty(),
        "simulate_lumped_board: ladder has no resonators"
    );
    assert!(
        cfg.dx_m.is_finite() && cfg.dx_m > 0.0,
        "simulate_lumped_board: dx_m must be positive and finite"
    );
    assert!(
        cfg.n_freq >= 1,
        "simulate_lumped_board: n_freq must be >= 1"
    );

    // Sweep frequencies (linear).
    let freqs: Vec<f64> = if cfg.n_freq == 1 {
        vec![0.5 * (cfg.f_lo_hz + cfg.f_hi_hz)]
    } else {
        (0..cfg.n_freq)
            .map(|i| {
                cfg.f_lo_hz + (cfg.f_hi_hz - cfg.f_lo_hz) * (i as f64) / (cfg.n_freq as f64 - 1.0)
            })
            .collect()
    };

    // DUT: filter elements present. THRU: elements removed (bare line) for the
    // normalization. Both run the identical board geometry / grid / drive.
    let dut = run_board_solve(ladder, substrate, cfg, &freqs, true);
    let thru = run_board_solve(ladder, substrate, cfg, &freqs, false);

    freqs
        .iter()
        .enumerate()
        .map(|(fi, &f)| {
            let v_dut = dut[fi].mag();
            let v_thru = thru[fi].mag();
            // S21 = transmission relative to the matched thru. A vanishing thru
            // (outside the launched spectrum) maps to 0 transmission rather than
            // a divide-by-zero blow-up.
            let s21 = if v_thru > 0.0 { v_dut / v_thru } else { 0.0 };
            (f, s21)
        })
        .collect()
}

/// One FDTD solve of the board (DUT if `place_elements`, else a bare thru line)
/// returning the single-bin DFT of the *load* port voltage at each sweep
/// frequency. Factored out so the DUT and thru runs share identical geometry,
/// grid, drive, and step body.
fn run_board_solve(
    ladder: &LumpedLadder,
    substrate: &Substrate,
    cfg: &LumpedSimConfig,
    freqs: &[f64],
    place_elements: bool,
) -> Vec<Bin> {
    // --- 1. Board + voxelize. ----------------------------------------------
    let board = lumped_board(ladder, substrate, cfg.footprint);
    let opts = VoxelOptions {
        dx_m: cfg.dx_m,
        xy_margin_cells: cfg.xy_margin_cells,
        air_above_cells: cfg.air_above_cells,
    };
    let model = voxelize_microstrip(&board.layout, &opts);
    assert!(
        model.port_cells.len() >= 2,
        "simulate_lumped_board: board must voxelize to >= 2 ports (drive + load); got {}",
        model.port_cells.len()
    );
    let dx = model.dx_m;
    let dt = model.grid.dt;
    let k_top = model.port_cells[0].2;
    // Place lumped elements on the substrate `E_z` edge directly under the top
    // metal (`k_top − 1`) — that is where the quasi-TEM vertical field lives,
    // the same column the propagation driver probes. The drive / load ports sit
    // at the trace plane `k_top`.
    let k_elem = k_top.saturating_sub(1).max(1);

    // Map a board (x, y) centre to the grid `(i, j)` column with the SAME
    // origin / floor convention `voxelize_microstrip` uses for ports.
    let x0 = board.layout.bbox.min.x - cfg.xy_margin_cells as f64 * dx;
    let y0 = board.layout.bbox.min.y - cfg.xy_margin_cells as f64 * dx;
    let (nx, ny, _nz) = model.dims;
    let cell_for = |cx: f64, cy: f64, k: usize| -> (usize, usize, usize) {
        let i = (((cx - x0) / dx).floor() as isize).clamp(0, nx as isize - 1) as usize;
        let j = (((cy - y0) / dx).floor() as isize).clamp(0, ny as isize - 1) as usize;
        (i, j, k)
    };

    let z0 = ladder.z0_ohm;

    // --- 2. Drive (input port) + matched load (output port). ----------------
    // Input port is the −x end (`port_cells[0]`), output the +x end (last).
    let drive_cell = model.port_cells[0];
    let load_cell = *model.port_cells.last().unwrap();

    let mut solver = WalkingSkeletonSolver::new(model.grid);

    // Modulated-Gaussian drive centred in-band. Centre the pulse a few time
    // constants in so its t = 0 tail is negligible.
    let bw = cfg.drive_bw_frac * cfg.drive_f0_hz;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;
    let wave = SourceWaveform::GaussianPulse {
        v0: cfg.drive_v0,
        f0: cfg.drive_f0_hz,
        bw,
        t0_steps,
    };
    let mut drive_port = LumpedRlcPort::pure_resistor(drive_cell, z0, wave);
    let mut load_port = LumpedRlcPort::pure_resistor(load_cell, z0, SourceWaveform::None);

    // --- 3. Per-branch lumped elements (DUT only). --------------------------
    let mut elements: Vec<LumpedRlcPort> = Vec::new();
    if place_elements {
        // `board.placements` is `L1, C1, L2, C2, …` — two footprints per
        // resonator (the inductor then the capacitor), in ladder order. The
        // resonator's L/C values live on the ladder; the footprint centres on
        // the placements. Pair them up by resonator index.
        for (ri, res) in ladder.resonators.iter().enumerate() {
            // The two placements for resonator `ri` (L then C).
            let l_pl: &Placement = &board.placements[2 * ri];
            let c_pl: &Placement = &board.placements[2 * ri + 1];
            debug_assert_eq!(
                l_pl.kind == BranchKind::Series,
                res.branch == LcBranch::Series
            );

            match res.branch {
                LcBranch::Series => {
                    // One series R-L-C in the through arm, on the in-line gap
                    // cell. Use the midpoint of the two footprint centres so the
                    // single element sits in the middle of the broken line span.
                    let cx = 0.5 * (l_pl.center_m.0 + c_pl.center_m.0);
                    let cy = 0.5 * (l_pl.center_m.1 + c_pl.center_m.1);
                    let cell = cell_for(cx, cy, k_elem);
                    elements.push(LumpedRlcPort::series_rlc(
                        cell,
                        SERIES_ESR_OHM,
                        res.l_henry,
                        res.c_farad,
                        SourceWaveform::None,
                    ));
                }
                LcBranch::Shunt => {
                    // Parallel L‖C from line to ground: two elements at the SAME
                    // shunt cell — a pure inductor (c = ∞) and a pure capacitor
                    // (l = 0). Their summed currents form the parallel admittance.
                    let cx = 0.5 * (l_pl.center_m.0 + c_pl.center_m.0);
                    let cy = 0.5 * (l_pl.center_m.1 + c_pl.center_m.1);
                    let cell = cell_for(cx, cy, k_elem);
                    elements.push(LumpedRlcPort::series_rlc(
                        cell,
                        SERIES_ESR_OHM,
                        res.l_henry,
                        f64::INFINITY,
                        SourceWaveform::None,
                    ));
                    elements.push(LumpedRlcPort::series_rlc(
                        cell,
                        SERIES_ESR_OHM,
                        0.0,
                        res.c_farad,
                        SourceWaveform::None,
                    ));
                }
            }
        }
    }

    // Sense the load-port voltage `V = E_z · dz`; DFT at every sweep frequency.
    let mut bins: Vec<Bin> = freqs.iter().map(|&f| Bin::new(2.0 * PI * f)).collect();

    // --- 4. Step loop. Custom body mirroring `run_line_eeff`: H + boundary,
    //        E + boundary, then every lumped-port `correct_e` (drive, load, and
    //        the filter elements), then advance the clock, then record. -------
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        solver.apply_cpml_e();

        drive_port.correct_e(solver.grid_mut(), n, dt);
        load_port.correct_e(solver.grid_mut(), n, dt);
        for el in elements.iter_mut() {
            el.correct_e(solver.grid_mut(), n, dt);
        }

        solver.advance_clock();

        let grid = solver.grid();
        let v_load = grid.ez[load_cell] * grid.dz;
        let t = n as f64 * dt;
        for b in bins.iter_mut() {
            b.accumulate(v_load, t);
        }
    }

    bins
}
