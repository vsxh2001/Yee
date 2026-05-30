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
//! # CW per-frequency S21 extraction (thru-normalized, F2.3-d, ADR-0128)
//!
//! The transmission is measured **one frequency at a time** with a CW
//! steady-state drive — NOT a single broadband pulse. The earlier F2.3-c driver
//! launched one modulated-Gaussian pulse and broadband-DFT'd the load voltage;
//! ADR-0127 proved the aperture lumped port is correct (under CW the capacitor
//! presents `1/(jωC)` and the shunt `L‖C` tank resonates), but a pulse + DFT
//! measures an *unsettled transient* on a short standing-wave line, so the
//! high-Q (Q≈10) tanks never reach steady state and no band-pass forms.
//!
//! So for each measured frequency `f`:
//!
//! 1. The drive port (a [`SourceWaveform::HannSine`] series EMF through a `Z0`
//!    resistor) drives a pure CW sinusoid at `f`, **Hann-ramped** over the first
//!    [`LumpedSimConfig::cw_ramp_cycles`] cycles to suppress the turn-on
//!    transient; the far port is a matched [`LumpedRlcPort::pure_resistor`]`(Z0)`
//!    load.
//! 2. The solve runs [`LumpedSimConfig::cw_ramp_cycles`] +
//!    [`LumpedSimConfig::cw_settle_cycles`] cycles so the highest-Q tank's
//!    ring-up **and** the source→load line transit settle into a
//!    single-frequency steady state.
//! 3. The load voltage `V = E_z · dz` is single-bin-DFT'd at `f` over the final
//!    [`LumpedSimConfig::cw_measure_cycles`] cycles (the **settled window only**,
//!    not the whole record) → the steady-state amplitude `|V_ss(f)|`.
//!
//! To divide out the (frequency-dependent, coarse-grid-dependent) feed + line +
//! port coupling and the line standing wave, the **same** board is run a second
//! time with the filter elements removed (a bare through line); then
//!
//! ```text
//! S21(f) = |V_dut,ss(f)| / |V_thru,ss(f)|
//! ```
//!
//! which is ≈ 1 (0 dB) in the passband and rolls off in the stopband — the
//! transmission *relative to the matched thru*, exactly the quantity
//! [`yee_filter::ladder_s21`] computes for the ideal circuit. This thru
//! calibration is robust against the lumped-port `E_z` voltage convention and
//! the residual feed/line coupling (see [`yee_fdtd::LumpedRlcPort`]). Each
//! frequency costs two full FDTD solves (DUT + thru), so the frequency set
//! ([`LumpedSimConfig::cw_freqs_hz`]) is deliberately small — the gate-check
//! points plus a handful for the sweep shape, NOT a fine sweep.

use std::f64::consts::PI;

use yee_fdtd::{ApertureSpec, LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver, YeeGrid};
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
/// tractable cubic grid and a **CW per-frequency steady-state** drive (F2.3-d,
/// ADR-0128) at the gate-check points plus a handful for the sweep shape. Heavy
/// (multiple multi-minute FDTD solves) → iterate in the bounded dev container
/// (CLAUDE.md §10).
///
/// # Why CW per frequency (F2.3-d, ADR-0128)
///
/// The earlier F2.3-c driver launched a *single* modulated-Gaussian pulse and
/// broadband-DFT'd the load voltage. ADR-0127 proved the aperture lumped port is
/// **correct** — under a CW drive the capacitor presents `1/(jωC)` and the shunt
/// `L‖C` tank **resonates** — but a pulse + broadband DFT measures an *unsettled
/// transient* on a short standing-wave line: the high-Q (Q≈10) tanks never reach
/// the steady-state reactance the band-pass needs, so the response stays flat.
///
/// The fix is a CW measurement: at each measured frequency `f` drive a pure
/// sinusoid (Hann-ramped over the first cycles to suppress the turn-on
/// transient), run enough cycles for the highest-Q tank plus the source→load
/// line transit to settle, then measure the **steady-state** load-voltage
/// amplitude over the last few settled cycles (single-bin DFT over the settled
/// window only). The same DUT/thru normalization divides out the line standing
/// wave (see the [module docs](self#s21-extraction-thru-normalized)).
#[derive(Debug, Clone, PartialEq)]
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
    /// Peak drive voltage (V) of the CW sinusoid.
    pub drive_v0: f64,
    /// Hann (raised-cosine) ramp length of the CW drive, in **carrier cycles**.
    /// The first `cw_ramp_cycles` cycles ramp the amplitude from 0 to full to
    /// suppress the turn-on transient; the source is a clean single tone
    /// thereafter (matches the `cap_cw_001` diagnostic's ramp idiom).
    pub cw_ramp_cycles: f64,
    /// Settling time *after* the ramp, in **carrier cycles**, before the
    /// steady-state measurement window opens. Must cover the highest-Q tank's
    /// ring-up (Q≈10 ⇒ ~10–30 cycles) **plus** the source→load line transit, so
    /// the measured amplitude is the true steady state rather than a transient.
    pub cw_settle_cycles: f64,
    /// Length of the steady-state measurement window, in **carrier cycles**.
    /// The load voltage is single-bin-DFT'd at `f` over this settled window
    /// only (a few clean cycles → one sharp single-tone bin).
    pub cw_measure_cycles: f64,
    /// The explicit set of CW measurement frequencies (Hz). Deliberately small
    /// (the gate-check points + a handful for the sweep shape) because each
    /// frequency costs two full FDTD solves (DUT + thru) — NOT a fine sweep.
    pub cw_freqs_hz: Vec<f64>,
}

impl Default for LumpedSimConfig {
    /// Walking-skeleton defaults for a ~2 GHz lumped band-pass on FR-4. See the
    /// per-field docs; container-validated to pass `fdtd_lumped_001` (ADR-0115,
    /// F2.3-d CW drive ADR-0128). The default frequency set spans 1.6–2.6 GHz
    /// and includes both gate-check points (2.0 GHz passband, 2.4 GHz stopband).
    fn default() -> Self {
        Self {
            dx_m: 0.4e-3,
            xy_margin_cells: 8,
            air_above_cells: 8,
            footprint: Footprint::Smd0603,
            drive_v0: 1.0,
            cw_ramp_cycles: 12.0,
            // The passband amplitude is settle-converged by ~60 cycles, but the
            // high-Q stopband notch keeps deepening with more settling (the tank
            // is still ringing up); 140 cycles measurably deepens the notch
            // (2.4 GHz rejection 2.7 → 5.1 dB on the 0.4 mm grid). See ADR-0128.
            cw_settle_cycles: 140.0,
            cw_measure_cycles: 16.0,
            cw_freqs_hz: vec![1.6e9, 1.8e9, 2.0e9, 2.2e9, 2.4e9, 2.6e9],
        }
    }
}

/// Single-bin DFT accumulator at one frequency: real/imag of `Σ x[n]·e^{-jωt}`.
///
/// Used to measure the **steady-state** load-voltage amplitude over the final
/// settled cycles of a CW solve (the settled window only — NOT the whole
/// record): `mag()` is the single-tone amplitude of `x[n]` at `ω` over the
/// accumulated samples (F2.3-d, ADR-0128).
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
/// pure-L ‖ pure-C apertures, see the [module docs](self)), then measures the
/// transmission **one frequency at a time** with a CW steady-state drive (F2.3-d,
/// ADR-0128): a Hann-ramped CW sinusoid at the input port through a `Z0` resistor
/// and a matched `Z0` load at the output, run until the high-Q tanks + line
/// transit settle, with the load voltage single-bin-DFT'd over the settled window
/// only. A second, element-free *thru* solve per frequency normalizes the
/// response, so the returned `|S21|` is the steady-state transmission relative to
/// the matched thru (≈ 0 dB in-band, rolling off in the stopband).
///
/// Returns `(freq_hz, |S21|)` for each frequency in
/// [`LumpedSimConfig::cw_freqs_hz`], in the order given.
///
/// # Panics
///
/// Panics if `ladder.resonators` is empty, if the board does not voxelize to at
/// least two ports (a drive + a load), or if [`LumpedSimConfig`] has a
/// non-positive `dx_m` or an empty `cw_freqs_hz`.
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
        !cfg.cw_freqs_hz.is_empty(),
        "simulate_lumped_board: cw_freqs_hz must be non-empty"
    );

    let freqs = &cfg.cw_freqs_hz;

    // DUT: filter elements present. THRU: elements removed (bare line) for the
    // normalization. Both run the identical board geometry / grid / CW drive at
    // each measured frequency. Each entry is the steady-state load-voltage
    // amplitude `|V_ss(f)|` measured over the settled window (F2.3-d, ADR-0128).
    let dut = run_board_solve(ladder, substrate, cfg, freqs, true);
    let thru = run_board_solve(ladder, substrate, cfg, freqs, false);

    freqs
        .iter()
        .enumerate()
        .map(|(fi, &f)| {
            let v_dut = dut[fi];
            let v_thru = thru[fi];
            // S21 = steady-state transmission relative to the matched thru. A
            // vanishing thru maps to 0 transmission rather than a
            // divide-by-zero blow-up.
            let s21 = if v_thru > 0.0 { v_dut / v_thru } else { 0.0 };
            (f, s21)
        })
        .collect()
}

/// The lumped elements for one resonator, captured as reusable aperture-port
/// build recipes (fixed geometry / values across all CW frequencies). One
/// series branch → one [`Self::Series`] aperture; one shunt branch → a
/// pure-inductor + pure-capacitor [`Self::Shunt`] aperture pair over the same
/// face (the parallel `L‖C`). Materialized into fresh passive
/// [`LumpedRlcPort`]s for every per-frequency CW solve (so each solve starts
/// from a clean port state).
enum ElementRecipe {
    /// Series R-L-C arm: `(spec, l, c)`.
    Series(ApertureSpec, f64, f64),
    /// Shunt parallel `L‖C`: `(spec, l, c)` → a pure-L aperture (`c = ∞`) ‖ a
    /// pure-C aperture (`l = 0`) over the same face.
    Shunt(ApertureSpec, f64, f64),
}

impl ElementRecipe {
    /// Materialize this recipe into fresh passive aperture ports.
    fn build(&self) -> Vec<LumpedRlcPort> {
        match self {
            ElementRecipe::Series(spec, l, c) => vec![LumpedRlcPort::aperture(
                spec.clone(),
                SERIES_ESR_OHM,
                *l,
                *c,
                SourceWaveform::None,
            )],
            ElementRecipe::Shunt(spec, l, c) => vec![
                LumpedRlcPort::aperture(
                    spec.clone(),
                    SERIES_ESR_OHM,
                    *l,
                    f64::INFINITY,
                    SourceWaveform::None,
                ),
                LumpedRlcPort::aperture(
                    spec.clone(),
                    SERIES_ESR_OHM,
                    0.0,
                    *c,
                    SourceWaveform::None,
                ),
            ],
        }
    }
}

/// CW per-frequency steady-state solve of the board (DUT if `place_elements`,
/// else a bare thru line), returning the **steady-state load-voltage amplitude**
/// `|V_ss(f)|` at each measured frequency (F2.3-d, ADR-0128).
///
/// The board is voxelized once; then for every frequency `f` a fresh solver
/// (cloned grid, zeroed fields) is driven with a Hann-ramped CW sinusoid at `f`,
/// run for `ramp + settle + measure` cycles, and the load voltage is single-bin
/// DFT'd over the **settled measurement window only** (the last
/// `cw_measure_cycles` cycles) — so the returned amplitude is the steady-state
/// response, not an unsettled transient. Factored out so the DUT and thru runs
/// share identical geometry, grid, drive, and step body.
fn run_board_solve(
    ladder: &LumpedLadder,
    substrate: &Substrate,
    cfg: &LumpedSimConfig,
    freqs: &[f64],
    place_elements: bool,
) -> Vec<f64> {
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

    // --- 3. Per-branch lumped-element recipes (DUT only). -------------------
    // Captured once (fixed geometry/values across frequencies) and materialized
    // into fresh passive aperture ports for each per-frequency CW solve.
    let mut recipes: Vec<ElementRecipe> = Vec::new();
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
                // Series R-L-C arm in the through path at the in-line gap column
                // — one aperture port carrying `Z = R + jωL + 1/(jωC)` over the
                // whole port face.
                LcBranch::Series => {
                    recipes.push(ElementRecipe::Series(spec, res.l_henry, res.c_farad))
                }
                // Parallel L‖C from line to ground (built as a pure-L ‖ pure-C
                // aperture pair over the same face — see [`ElementRecipe`]). The
                // dominant shunt tanks set the band-pass selectivity.
                LcBranch::Shunt => {
                    recipes.push(ElementRecipe::Shunt(spec, res.l_henry, res.c_farad))
                }
            }
        }
    }

    // --- 4. CW per-frequency steady-state loop (F2.3-d, ADR-0128). ----------
    // For each measured frequency, drive a Hann-ramped CW sinusoid, let the
    // high-Q tanks + the line transit settle, then single-bin DFT the load
    // voltage over the settled measurement window ONLY → the steady-state
    // amplitude `|V_ss(f)|`.
    freqs
        .iter()
        .map(|&f| {
            cw_steady_state_amplitude(
                model.grid.clone(),
                cfg,
                f,
                z0,
                drive_cell,
                load_cell,
                &recipes,
            )
        })
        .collect()
}

/// Run one CW steady-state FDTD solve at frequency `f` on a fresh `grid` (zeroed
/// fields) and return the steady-state load-voltage amplitude `|V_ss(f)|`.
///
/// Mirrors the `cap_cw_001` CW idiom: a Hann-ramped sinusoidal series EMF at the
/// drive port suppresses the turn-on transient; after `cw_ramp_cycles +
/// cw_settle_cycles` the fields reach a single-frequency steady state (the
/// highest-Q tank ring-up + the source→load line transit); the load voltage is
/// single-bin DFT'd over the final `cw_measure_cycles` (the settled window) into
/// the steady-state amplitude. The step body matches `run_line_eeff`: H +
/// boundary, E + boundary, the drive/load single-edge `correct_e`, the filter
/// elements' multi-cell `correct_e_aperture`, then advance the clock and record.
#[allow(clippy::too_many_arguments)]
fn cw_steady_state_amplitude(
    grid: YeeGrid,
    cfg: &LumpedSimConfig,
    f: f64,
    z0: f64,
    drive_cell: (usize, usize, usize),
    load_cell: (usize, usize, usize),
    recipes: &[ElementRecipe],
) -> f64 {
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    // Cycle counts → step counts at THIS frequency.
    let steps_per_cycle = (1.0 / (f * dt)).round().max(1.0) as usize;
    let ramp_steps = (cfg.cw_ramp_cycles * steps_per_cycle as f64).round() as usize;
    let settle_steps = (cfg.cw_settle_cycles * steps_per_cycle as f64).round() as usize;
    let measure_steps = (cfg.cw_measure_cycles * steps_per_cycle as f64)
        .round()
        .max(1.0) as usize;
    // The steady-state window opens after the ramp + settle and runs to the end.
    let measure_start = ramp_steps + settle_steps;
    let n_steps = measure_start + measure_steps;

    // Hann-ramped CW sinusoid at `f` driving the input port through `Z0`; matched
    // `Z0` load at the output port.
    let wave = SourceWaveform::HannSine {
        v0: cfg.drive_v0,
        frequency: f,
        ramp_steps,
    };
    let mut drive_port = LumpedRlcPort::pure_resistor(drive_cell, z0, wave);
    let mut load_port = LumpedRlcPort::pure_resistor(load_cell, z0, SourceWaveform::None);

    // Fresh passive aperture ports for this solve.
    let mut elements: Vec<LumpedRlcPort> = recipes.iter().flat_map(ElementRecipe::build).collect();

    // Single-bin DFT of the load voltage, accumulated over the SETTLED window
    // only (not the whole record) → the steady-state single-tone amplitude.
    let mut bin = Bin::new(2.0 * PI * f);

    for n in 0..n_steps {
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

        if n >= measure_start {
            let grid = solver.grid();
            let v_load = grid.ez[load_cell] * grid.dz;
            let t = n as f64 * dt;
            bin.accumulate(v_load, t);
        }
    }

    bin.mag()
}
