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
//! # Aperture-port placement (Filter Phase F2.3-c, ADR-0126)
//!
//! The earlier F2.3-b approach placed each element as a value-distributed
//! full-width *sheet* of single-edge [`LumpedRlcPort::series_rlc`] ports
//! (`C/N`, `N·L` per cell). The dx-sweep investigation behind ADR-0125 showed
//! that loads the line but **cannot resonate**: a single-edge port references
//! its field coupling to one Yee cell, so under grid refinement the inductor's
//! two-way back-action collapses as **O(dx²)** while the capacitor freezes at a
//! fixed per-cell short. A parallel L‖C needs both arms balanced; the inductor
//! goes inert while the cap stays a short, so the sheet degenerates to a
//! frequency-flat shunt capacitance that the DUT/thru normalization divides out
//! to `|S21| ≈ 1.0` — no notch ever forms (ADR-0124/0125 Outcome).
//!
//! So each ladder element is now placed with the **multi-cell aperture lumped
//! port** ([`LumpedRlcPort::aperture`], Phase 2.fdtd.6.9, ADR-0125), one
//! aggregate-`R/L/C` branch per ladder element. The aperture port references
//! its modal terminal voltage `V = ∫E_z·dz` to the **full substrate height**
//! and its field back-action to the **physical port-face area `A = w·h`**
//! (trace width × substrate height), NOT a single `dx²` cell — which removes the
//! `O(dx²)` inductor collapse and presents a dx-stable reactance, so the L‖C
//! tanks resonate.
//!
//! The `(y, z)` port-face aperture at an element's x-column is built from:
//!
//! - the transverse trace band `[j_lo, j_hi)` (the `y` rows the signal line
//!   spans), read directly off the voxel model's top-metal PEC mask at the
//!   drive-port column (the contiguous copper run containing the port row) so it
//!   tracks the Hammerstad-Jensen trace width at any `dx` / `Z0`; and
//! - the substrate height cells `k = 0 .. k_top` (= `0 .. n_sub`), the `E_z`
//!   edges spanning the ground-to-trace gap the quasi-TEM mode lives in (the
//!   [`crate::voxelize_microstrip`] z-stack: ground at `k = 0`, trace at
//!   `k_top = n_sub`).
//!
//! The [`ApertureSpec`] tiles every `(j, k)` `E_z` edge in that band at the
//! element's x-plane, with `n_columns = j_hi − j_lo` (so a wider aperture does
//! not multiply the modal `V`), physical area `A = w·h` and substrate height
//! `h = k_top·dx`. **No per-cell value-splitting** — the aperture port carries
//! the aggregate element value and handles the modal `V` + area-`A` back-action
//! internally.
//!
//! Each resonator topology maps onto the aperture port as:
//!
//! - A **shunt-branch** resonator (a *parallel* L‖C from line to ground) becomes
//!   **two** aperture ports over the same face: a pure-inductor aperture
//!   (`aperture(spec, ESR, L, ∞, None)`, `c = ∞` shorts the cap) **in parallel
//!   with** a pure-capacitor aperture (`aperture(spec, ESR, 0, C, None)`,
//!   `l = 0` removes the inductor). Both correct the same `E_z` edges each step,
//!   so the lattice sees their currents summed → the parallel `L‖C` admittance
//!   `Y = jωC + 1/(jωL)`, the correct shunt topology.
//! - A **series-branch** resonator (a series L–C in the through arm) becomes one
//!   series-RLC aperture (`aperture(spec, ESR, L, C, None)`) at the in-line gap
//!   column. The dominant **shunt** resonators set the band-pass selectivity.
//!
//! A small finite series resistance ([`SERIES_ESR_OHM`]) is used on every
//! element because [`LumpedRlcPort::aperture`] requires `r > 0`.
//!
//! An aperture port is **always two-way** coupled (the lumped branch current is
//! solved implicitly with the field and feeds back into `E_z`), so a
//! source-free inductor is not inert and the L‖C resonates. The aperture ports
//! are stepped with [`LumpedRlcPort::correct_e_aperture`]; the matched
//! drive/load resistors stay on the single-edge `pure_resistor` /
//! [`LumpedRlcPort::correct_e`] path (a pure resistor reflects identically).
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
//! the residual feed/line coupling (see [`yee_fdtd::LumpedRlcPort`]).

use std::f64::consts::PI;

use yee_fdtd::{ApertureSpec, LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver};
use yee_filter::{BranchKind, Footprint, LcBranch, LumpedLadder, Placement, lumped_board};
use yee_layout::Substrate;

use crate::{MicrostripModel, VoxelOptions, voxelize_microstrip};

/// Series ESR (Ω) used on every lumped element.
///
/// [`LumpedRlcPort::aperture`] requires `r > 0`; a tiny value approximates an
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
    /// launched pulse to fully transit the board *and* for the lumped
    /// capacitor's slow integrator tail to reach steady state — a band-pass
    /// shaped by L‖C tanks needs the steady-state reactance, and ADR-0125
    /// flagged that a single short pulse reads the cap as a near-short. The
    /// linear-system DFT-of-pulse equals the transfer function only if the
    /// record captures the full (slow) response, so the default is generous.
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
            n_steps: 24_000,
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

/// Find the transverse signal-line band `[j_lo, j_hi)` (a half-open `y`-row
/// range, width `N = j_hi − j_lo`) at the drive-port column, over which a
/// lumped element is spread as a full-width sheet (ADR-0124).
///
/// The signal line is the contiguous run of top-metal copper in `y` containing
/// the drive port's row, read off the voxel model's `Ex` PEC mask at the
/// top-metal plane `k_top` (`pec_mask_ex` is `true` where a trace covers the
/// `Ex` node). Scanning the mask makes `N` track the Hammerstad-Jensen trace
/// width at any `dx` / `Z0` without re-deriving it. If the grid carries no
/// `Ex` PEC mask, fall back to the single port row (`[j, j+1)`), so the driver
/// still produces a valid (if narrow) placement rather than panicking.
fn line_band_at(model: &MicrostripModel, port_cell: (usize, usize, usize)) -> (usize, usize) {
    let (i, j_port, k_top) = port_cell;
    let (_nx, ny, _nz) = model.dims;
    let Some(mask) = model.grid.pec_mask_ex.as_ref() else {
        return (j_port, j_port + 1);
    };
    // `pec_mask_ex` shape is `(nx, ny+1, nz+1)`; the Ex node row index runs
    // `0..=ny`. Clamp the port column / row into range defensively.
    let (mnx, mny, mnz) = mask.dim();
    let ic = i.min(mnx.saturating_sub(1));
    let kc = k_top.min(mnz.saturating_sub(1));
    let jp = j_port.min(mny.saturating_sub(1));
    let covered = |j: usize| j < mny && mask[(ic, j, kc)];
    if !covered(jp) {
        // Port row not under copper (unexpected) — single-row fallback.
        return (j_port, (j_port + 1).min(ny));
    }
    // Walk down then up from the port row to the contiguous-copper extent.
    let mut lo = jp;
    while lo > 0 && covered(lo - 1) {
        lo -= 1;
    }
    let mut hi = jp + 1; // exclusive
    while hi <= ny && covered(hi) {
        hi += 1;
    }
    (lo, hi)
}

/// EM-simulate a synthesized lumped-LC filter board and return its forward
/// transmission `|S21|(f)` over a frequency sweep.
///
/// Builds the board ([`yee_filter::lumped_board`]), voxelizes it
/// ([`crate::voxelize_microstrip`]), places each ladder element with a
/// **multi-cell aperture lumped port** ([`LumpedRlcPort::aperture`]) over the
/// `(y, z)` port-face aperture (trace width × substrate height) at the
/// element's x-column (series branch → one series-RLC aperture; shunt branch →
/// pure-L ‖ pure-C apertures, see the [module docs](self)), drives the input
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
    // The drive / load ports sit at the trace plane `k_top = n_sub` (ground at
    // `k = 0`, dielectric `E_z` edges `k = 0 .. n_sub`). The aperture lumped
    // ports span the full substrate height — the `E_z` edges `k = 0 .. k_top`
    // (= `0 .. n_sub`) that the quasi-TEM vertical field lives in.
    let k_top = model.port_cells[0].2;

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

    // Transverse signal-line band `[j_lo, j_hi)` (width `N` cells) — the `y`
    // extent of the `(y, z)` port-face aperture each lumped element bridges (so
    // it couples to the *whole* line admittance via the modal port face, not a
    // single ≈inert edge — ADR-0125/0126). Read it from the top-metal PEC mask
    // at the drive-port column: the contiguous copper run in `y` containing the
    // port row IS the signal line. Falls back to the single port row if the mask
    // is unavailable (so a mask-less grid still produces a valid placement).
    let (j_lo, j_hi) = line_band_at(&model, model.port_cells[0]);
    // Physical aperture geometry, held in metres (dx-stable): substrate height
    // `h = k_top·dx` (the `E_z` edges `k = 0 .. k_top`), trace width
    // `w = (j_hi − j_lo)·dx`, port-face area `A = w·h`. `n_columns = j_hi − j_lo`
    // (the modal voltage averages over the width columns, so a wider aperture
    // does not multiply the modal `V`).
    let n_columns = (j_hi - j_lo).max(1);
    let ap_height = k_top as f64 * dx;
    let ap_width = n_columns as f64 * dx;
    let ap_area = ap_width * ap_height;

    // Build the `(y, z)` aperture cells at a given x-plane `i_col`: every `E_z`
    // edge `(i_col, j, k)` over the trace band `j ∈ [j_lo, j_hi)` × the
    // substrate height `k ∈ [0, k_top)`.
    let aperture_cells = |i_col: usize| -> Vec<(usize, usize, usize)> {
        let mut v = Vec::with_capacity(n_columns * k_top.max(1));
        for j in j_lo..j_hi {
            for k in 0..k_top {
                v.push((i_col, j, k));
            }
        }
        v
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

            // x-column of this resonator: the midpoint of its two footprint
            // centres. The element bridges the `(y, z)` port-face aperture (the
            // trace band `[j_lo, j_hi)` × the substrate height `k = 0 .. k_top`)
            // at this column via ONE aggregate-`R/L/C` aperture port per branch
            // — no per-cell value-splitting; the aperture port (Phase 2.fdtd.6.9,
            // ADR-0125) handles the modal `V = ∫E_z·dz` (full height) and the
            // physical-area `A = w·h` back-action internally, removing the
            // O(dx²) inductor collapse so the L‖C tanks resonate. An aperture
            // port is always two-way.
            let cx = 0.5 * (l_pl.center_m.0 + c_pl.center_m.0);
            let i_col = cell_for(cx, 0.0, k_top).0;
            let spec = ApertureSpec {
                cells: aperture_cells(i_col),
                n_columns,
                area: ap_area,
                height: ap_height,
            };

            match res.branch {
                LcBranch::Series => {
                    // Series R-L-C arm in the through path at the in-line gap
                    // column — one aperture port carrying the aggregate
                    // `Z = R + jωL + 1/(jωC)` over the whole port face.
                    elements.push(LumpedRlcPort::aperture(
                        spec,
                        SERIES_ESR_OHM,
                        res.l_henry,
                        res.c_farad,
                        SourceWaveform::None,
                    ));
                }
                LcBranch::Shunt => {
                    // Parallel L‖C from line to ground, as TWO aperture ports
                    // over the SAME port face: a pure-inductor aperture
                    // (c = ∞ shorts the cap) in parallel with a pure-capacitor
                    // aperture (l = 0 removes the inductor). Their summed
                    // currents form the parallel admittance
                    // `Y = jωC + 1/(jωL)` loading the line.
                    elements.push(LumpedRlcPort::aperture(
                        spec.clone(),
                        SERIES_ESR_OHM,
                        res.l_henry,
                        f64::INFINITY,
                        SourceWaveform::None,
                    ));
                    elements.push(LumpedRlcPort::aperture(
                        spec,
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
    //        E + boundary, then the drive/load single-edge `correct_e` and the
    //        filter elements' multi-cell `correct_e_aperture`, then advance the
    //        clock, then record. ---------------------------------------------
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        solver.apply_cpml_e();

        drive_port.correct_e(solver.grid_mut(), n, dt);
        load_port.correct_e(solver.grid_mut(), n, dt);
        for el in elements.iter_mut() {
            el.correct_e_aperture(solver.grid_mut(), n, dt);
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
