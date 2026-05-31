//! Lumped-LC **FDTD EM simulation** of a synthesized filter board (Filter
//! Phase F2.3, ADR-0115).
//!
//! [`simulate_lumped_board`] takes a synthesized [`yee_filter::LumpedLadder`]
//! (F2.0), places it as SMD footprints on a microstrip board with
//! [`yee_filter::lumped_board`] (F2.2), **lengthens the signal line** and
//! voxelizes that board onto a Yee grid with [`crate::voxelize_microstrip`]
//! (F1.1a), drops a [`LumpedRlcPort`] at each placed L/C, terminates the line
//! with an **x-only CPML matched termination** at both ends, drives a CW soft
//! source at one end, time-steps the FDTD, and extracts the forward
//! transmission `|S21|(f)` over a frequency sweep — the lumped analogue of the
//! distributed [`crate::run_line_eeff`] gate.
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
//! # CW per-frequency S21 extraction (PEC-box 2-point standing-wave de-embed, F2.3-g, ADR-0132)
//!
//! The transmission is measured **one frequency at a time** with a CW
//! steady-state drive — NOT a single broadband pulse. The earlier F2.3-c driver
//! launched one modulated-Gaussian pulse and broadband-DFT'd the load voltage;
//! ADR-0127 proved the aperture lumped port is correct (under CW the capacitor
//! presents `1/(jωC)` and the shunt `L‖C` tank resonates), but a pulse + DFT
//! measures an *unsettled transient* on a short standing-wave line, so the
//! high-Q (Q≈10) tanks never reach steady state and no band-pass forms.
//!
//! ## Why a PEC box + a 2-point standing-wave de-embed (F2.3-g, ADR-0132)
//!
//! F2.3-d's CW drive cured the transient, but **every** prior de-embed
//! (short-board F2.3-c/-d, finer-grid F2.3-e, matched-CPML F2.3-f) gave a
//! **monotone + over-unity** `|S21|` — never a band-pass. Two distinct failures
//! compounded:
//!
//! 1. **Over-unity is a bad-de-embed signature.** All the prior measurements
//!    took a *single* reference-plane voltage and never separated the forward
//!    (incident / transmitted) wave from the backward (reflected) wave. On a
//!    short, mismatched line the standing wave's peaks/nulls move with frequency,
//!    so the single-point DUT/thru ratio reads `> 1` in places — unphysical for a
//!    passive filter.
//! 2. **CPML-into-substrate is unstable.** F2.3-f's x-only CPML termination hit
//!    the *documented* instability ([`crate::run_line_eeff`] / ADR-0108): CPML on
//!    a microstrip whose PEC ground + high-ε substrate run into the boundary is
//!    late-time unstable, and it still read monotone + over-unity.
//!
//! F2.3-g fixes both. It mirrors the **stable** [`crate::run_line_eeff`] pattern:
//! a **plain PEC box** ([`WalkingSkeletonSolver::new`] — no CPML; the
//! `apply_cpml_*` calls fall back to a hard-PEC clamp), with the line lengthened
//! so the elements clear the ends and a steady standing wave develops. Then it
//! does a proper **2-point (here 3-point) standing-wave de-embed** to separate
//! the forward and backward travelling waves.
//!
//! ### The standing-wave decomposition
//!
//! On a uniform section of the guided line the steady-state CW voltage phasor is
//!
//! ```text
//! V(x) = a · e^{−jβx} + b · e^{+jβx}
//! ```
//!
//! a forward (`+x`) travelling amplitude `a` plus a backward (`−x`) amplitude
//! `b`. We single-bin-DFT the steady-state `E_z` (×`dz` → modal voltage) phasor
//! at **three equally-spaced columns** `x₀, x₀+d, x₀+2d` in each port reference
//! region. The standing wave obeys the exact recurrence
//!
//! ```text
//! V₀ + V₂ = 2·cos(βd)·V₁   ⟹   cos(βd) = (V₀ + V₂) / (2·V₁)
//! ```
//!
//! so `β` is read **directly from the data** (`β = acos(Re[(V₀+V₂)/(2V₁)]) / d`)
//! — self-consistent with the actual FDTD grid's numerical dispersion, needing
//! **no separate ε_eff calibration run**. (`Im[(V₀+V₂)/(2V₁)] → 0` as the fit is
//! consistent; it is reported as a quality check.) With `β` known, the 2×2 system
//! from `V₀` and `V₂` (local origin `x₀ = 0`) solves for the forward `a` and
//! backward `b`:
//!
//! ```text
//! V₀ = a + b
//! V₂ = a·e^{−j2βd} + b·e^{+j2βd}
//! ```
//!
//! ## The clean forward-wave launch (F2.3-h, ADR-0133)
//!
//! F2.3-g made the de-embed **physical** (no over-unity) but exposed a new
//! limiter: a single-column **soft `E_z` source** in a PEC box launches
//! *symmetrically* (equal `+x` and `−x` halves), so the `−x` half reflects off
//! the input wall and the input region becomes a **near-pure standing wave**
//! (`|fwd a₁| ≈ |bwd refl|`), while only a small, cavity-resonance-dependent
//! fraction of forward power reaches the output region — the bare thru read
//! `β_out = 0` at some gate freqs and `|b₂|` at the floor (`~0.02`). So the
//! `S21 = (b₂/a₁)_dut/(b₂/a₁)_thru` ratio divided tiny, partly-degenerate output
//! readings and the "notch at f0" was a likely floor artifact, not a real result.
//!
//! F2.3-h gives the launch two fixes so `a₁` and `b₂` become **trustworthy**
//! (`β > 0` at all gate freqs, `b₂` well above the floor):
//!
//! 1. **Directional two-column phased source.** Instead of one soft-`E_z` column,
//!    the launcher drives **two adjacent columns** `src_i` and `src_i + 1` with
//!    the downstream column *retarded* by the one-cell wave-transit phase
//!    `Δφ = β·dx` (β from the calibration pre-pass below). A pair of soft sources
//!    one cell apart, the downstream one delayed by the inter-cell transit time,
//!    **adds constructively in `+x` and destructively in `−x`** — a poor-man's
//!    TF/SF / Huygens launcher (cf. ADR-0014/0021/0026). It injects
//!    *predominantly forward*, so the input is no longer a near-pure standing wave
//!    and real forward power reaches the output region (raising `|b₂|`, fixing
//!    `β_out = 0`).
//! 2. **Time-gated incident `a₁` reference (`run_line_eeff` pattern, ADR-0108).**
//!    A separate, short **pulse** pre-pass (`calibrate_launch`) launches a
//!    modulated-Gaussian from the same directional source and **time-gates** the
//!    DFT to the forward passage at the input reference region, *before* the first
//!    far-wall reflection returns. With no reflected wave in the window the 3-point
//!    fit reads a **pure forward** incident amplitude `a₁_gated` (`β > 0`,
//!    `|bwd| ≈ 0`) — a trustworthy launch reference and the numerical `β` used to
//!    phase the directional source. Because `a₁_gated` is a property of the launch
//!    + lead-in line only (identical for the DUT and thru runs), it cancels in the
//!    thru-normalization, but it is kept explicit so the scheme is verifiable.
//!
//! ### S21 = forward-out / forward-in, thru-normalized (hybrid)
//!
//! For each measured frequency `f` (the **hybrid**: a time-gated incident `a₁`
//! reference + a CW-settled `b₂`):
//!
//! 1. A **directional two-column phased `HannSine` soft `E_z` source sheet**
//!    across the strip's `(y, z)` face at the input launches a predominantly `+x`
//!    travelling quasi-TEM wave, **Hann-ramped** over the first
//!    [`LumpedSimConfig::cw_ramp_cycles`] cycles to suppress the turn-on
//!    transient.
//! 2. The CW solve runs [`LumpedSimConfig::cw_ramp_cycles`] +
//!    [`LumpedSimConfig::cw_settle_cycles`] cycles so the highest-Q (Q≈10) tank's
//!    ring-up **and** the line transit settle into a single-frequency steady
//!    state on the PEC line (the tanks **must** ring up to show the band-pass, so
//!    the DUT response stays CW).
//! 3. The 3-point standing-wave fit at the **output** reference region (downstream
//!    of the last element, on the lengthened lead-in clear of the far wall) gives
//!    the transmitted forward amplitude `b₂` (the `+x` part). The incident forward
//!    amplitude `a₁` is the **time-gated** `a₁_gated` from `calibrate_launch` (the
//!    trustworthy clean-launch reference) rather than the CW-input standing-wave
//!    fit.
//!
//! The raw transmission is `(b₂/a₁)`. To divide out the (frequency-dependent,
//! coarse-grid-dependent) feed + line + port coupling, the **same** board is run
//! a second time with the filter elements removed (a bare through line); then
//!
//! ```text
//! S21(f) = (b₂/a₁)_dut / (b₂/a₁)_thru
//! ```
//!
//! which is ≈ 1 (0 dB) in the passband and rolls off in the stopband — the
//! transmission *relative to the thru*, exactly the quantity
//! [`yee_filter::ladder_s21`] computes for the ideal circuit. Because the forward
//! and backward waves are separated, the standing wave no longer corrupts the
//! ratio (no over-unity); the directional launch + time-gated `a₁` keep both
//! amplitudes well above the floor; and the PEC box is unconditionally stable (no
//! CPML divergence). Each frequency costs two full FDTD solves (DUT + thru) plus a
//! shared one-off calibration pulse, so the frequency set
//! ([`LumpedSimConfig::cw_freqs_hz`]) is deliberately small — the gate-check
//! points plus a handful for the sweep shape, NOT a fine sweep.

use std::f64::consts::PI;

use yee_fdtd::{ApertureSpec, LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver};
use yee_filter::{BranchKind, Footprint, LcBranch, LumpedLadder, Placement, lumped_board};
use yee_layout::{Point2, Polygon, Substrate};

use crate::{C0_M_S, MicrostripModel, VoxelOptions, voxelize_microstrip};

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
    /// PEC end-wall **guard band** in cells (F2.3-g, ADR-0132). The line runs in
    /// a plain PEC box (NO CPML — the stable [`crate::run_line_eeff`] pattern,
    /// per ADR-0108: CPML on a microstrip whose PEC ground / high-ε substrate run
    /// into the boundary is late-time unstable). `npml` cells of clear line are
    /// kept between each `x`-end PEC wall and the nearest source / probe column,
    /// so the launcher and the standing-wave reference regions sit on undisturbed
    /// line rather than right against the reflecting wall. (The field name is
    /// retained from the F2.3-f CPML driver to keep the public config API stable;
    /// it now sizes a PEC guard band, not an absorber thickness.)
    pub npml: usize,
    /// Length (in cells) of the straight microstrip **lead-in** added to each
    /// end of the synthesized board before voxelizing (F2.3-g, ADR-0132). The
    /// board's signal trace is extended `lead_in_cells` cells past each port so
    /// the source and the two standing-wave reference regions sit on a uniform
    /// section of line, clear of the element + original-port discontinuities, and
    /// the line is long enough past each port that a steady standing wave
    /// develops on the PEC box. The input reference region sits just downstream
    /// of the source; the output reference region just upstream of the far PEC
    /// wall, downstream of the last element.
    pub lead_in_cells: usize,
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
            // PEC box (F2.3-g, ADR-0132): the standing wave must fully establish
            // — multiple end-to-end transits + the high-Q (Q≈10) tank ring-up
            // settle into a single-frequency steady state before the phasors are
            // read. 120 cycles is a generous margin on a ~few-λ_g PEC line while
            // staying bounded (a single solve stays in the minutes, not hours).
            cw_settle_cycles: 120.0,
            cw_measure_cycles: 16.0,
            cw_freqs_hz: vec![1.6e9, 1.8e9, 2.0e9, 2.2e9, 2.4e9, 2.6e9],
            // PEC end-wall guard band (NOT a CPML thickness — F2.3-g, ADR-0132):
            // cells of clear line kept between each x-end PEC wall and the
            // nearest source / probe column.
            npml: 6,
            // Straight-line lead-in past each port (F2.3-g, ADR-0132): long
            // enough that the source + the two 3-point standing-wave reference
            // regions (each spanning 2·d ≈ 2·λ_g/12 ≈ 36 cells past the npml=6
            // guard) sit on a uniform section of line clear of the element /
            // original-port discontinuities, with a steady standing wave. 52
            // cells (~21 mm at 0.4 mm) fits npml(6)+guard(3)+2d(~36)+margin and
            // keeps the run tractable (a single solve stays in the minutes).
            lead_in_cells: 52,
        }
    }
}

/// A minimal complex number `(re, im)` for the standing-wave de-embed math
/// (F2.3-g, ADR-0132). `yee-voxel` carries no `num-complex` dependency, so —
/// matching the by-hand phasor arithmetic already used in [`Bin`] and
/// [`crate::run_line_eeff`] — the few complex operations the 3-point fit needs
/// are implemented inline.
#[derive(Clone, Copy, Debug)]
struct Cplx {
    re: f64,
    im: f64,
}

impl Cplx {
    const ZERO: Cplx = Cplx { re: 0.0, im: 0.0 };

    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// `e^{jθ}` (unit phasor).
    fn expj(theta: f64) -> Self {
        Self {
            re: theta.cos(),
            im: theta.sin(),
        }
    }

    fn add(self, o: Cplx) -> Cplx {
        Cplx::new(self.re + o.re, self.im + o.im)
    }

    fn sub(self, o: Cplx) -> Cplx {
        Cplx::new(self.re - o.re, self.im - o.im)
    }

    fn mul(self, o: Cplx) -> Cplx {
        Cplx::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }

    fn scale(self, s: f64) -> Cplx {
        Cplx::new(self.re * s, self.im * s)
    }

    /// `self / o`.
    fn div(self, o: Cplx) -> Cplx {
        let d = o.re * o.re + o.im * o.im;
        Cplx::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }

    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
}

/// Single-bin DFT accumulator at one frequency: real/imag of `Σ x[n]·e^{-jωt}`.
///
/// Used to measure the **steady-state** line-voltage phasor over the final
/// settled cycles of a CW solve (the settled window only — NOT the whole
/// record): [`Self::phasor`] is the complex single-tone amplitude of `x[n]` at
/// `ω` over the accumulated samples (F2.3-d/-g, ADR-0128/0132).
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
    /// The accumulated complex phasor `Σ x[n]·e^{-jωt}` (un-normalized; only
    /// phasor *ratios* are used downstream, so the common `1/N` cancels).
    fn phasor(&self) -> Cplx {
        Cplx::new(self.re, self.im)
    }
}

/// The travelling-wave decomposition `V(x) = a·e^{−jβx} + b·e^{+jβx}` of a
/// standing-wave region, fitted from 3 equally-spaced phasors (F2.3-g,
/// ADR-0132).
#[derive(Clone, Copy, Debug)]
struct StandingWaveFit {
    /// Forward (`+x`-travelling) complex amplitude `a` at the region origin.
    fwd: Cplx,
    /// Backward (`−x`-travelling) complex amplitude `b` at the region origin.
    bwd: Cplx,
    /// Extracted phase constant `β` (rad/m) — self-consistent with the FDTD
    /// grid's numerical dispersion (read from the data, not a calibration).
    beta: f64,
    /// `|Im[(V₀+V₂)/(2V₁)]|` — the imaginary residual of the `cos(βd)`
    /// recurrence; → 0 for a consistent fit (a quality diagnostic).
    cos_imag_residual: f64,
}

/// Fit `V(x) = a·e^{−jβx} + b·e^{+jβx}` from three equally-spaced steady-state
/// phasors `v0 = V(x₀)`, `v1 = V(x₀+d)`, `v2 = V(x₀+2d)` of spacing `d` metres
/// (F2.3-g, ADR-0132). Returns the forward `a` / backward `b` amplitudes at the
/// region origin `x₀`, the data-derived `β`, and the recurrence residual.
///
/// `β` comes from the exact 3-point recurrence `V₀ + V₂ = 2·cos(βd)·V₁`, so it
/// is self-consistent with the FDTD grid (no separate ε_eff calibration). With
/// `β` known, the 2×2 system from `V₀, V₂` (local origin `x₀ = 0`) gives `a, b`:
/// `V₀ = a + b`, `V₂ = a·e^{−j2βd} + b·e^{+j2βd}`.
fn fit_standing_wave(v0: Cplx, v1: Cplx, v2: Cplx, d: f64) -> StandingWaveFit {
    // cos(βd) = (V₀ + V₂) / (2·V₁). The real part is the physical cosine; the
    // imaginary part is a consistency residual (→ 0 for a clean standing wave).
    let cos_bd = v0.add(v2).div(v1.scale(2.0));
    let cos_imag_residual = cos_bd.im.abs();
    // acos needs a real argument in [-1, 1]; clamp so grid noise that pushes the
    // measured cosine slightly outside the unit interval does not produce a NaN.
    let c = cos_bd.re.clamp(-1.0, 1.0);
    let beta_d = c.acos();
    // d > 0 always (caller passes a positive cell spacing); β in rad/m.
    let beta = beta_d / d;

    // Solve V₀ = a + b, V₂ = a·p + b·q with p = e^{−j2βd}, q = e^{+j2βd}.
    // a = (V₀·q − V₂) / (q − p),  b = V₀ − a.
    let p = Cplx::expj(-2.0 * beta_d);
    let q = Cplx::expj(2.0 * beta_d);
    let denom = q.sub(p);
    let (fwd, bwd) = if denom.abs() > 1e-12 {
        let a = v0.mul(q).sub(v2).div(denom);
        let b = v0.sub(a);
        (a, b)
    } else {
        // Degenerate (βd ≈ 0 or π → q ≈ p): cannot separate the two waves.
        // Fall back to "all forward" so the caller still gets a finite number;
        // flagged by the residual / β being at a band edge.
        (v0, Cplx::ZERO)
    };

    StandingWaveFit {
        fwd,
        bwd,
        beta,
        cos_imag_residual,
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
/// Builds the board ([`yee_filter::lumped_board`]), **lengthens its signal
/// line**, voxelizes it ([`crate::voxelize_microstrip`]), places each ladder
/// element with a **multi-cell aperture lumped port** ([`LumpedRlcPort::aperture`])
/// over the `(y, z)` port-face aperture (trace width × substrate height) at the
/// element's x-column (series branch → one series-RLC aperture; shunt branch →
/// pure-L ‖ pure-C apertures, see the [module docs](self)), then measures the
/// transmission **one frequency at a time** with a CW steady-state drive in a
/// **plain PEC box** + a **3-point standing-wave de-embed** (F2.3-g, ADR-0132):
/// a Hann-ramped CW soft `E_z` source sheet at the input launches a clean
/// traveling wave on a PEC-bounded line (the stable [`crate::run_line_eeff`]
/// pattern — no CPML), run until the high-Q tanks + the standing wave settle,
/// then the steady-state `E_z` phasor is sampled at three equally-spaced columns
/// in the **input** and **output** reference regions and fitted to
/// `V(x) = a·e^{−jβx} + b·e^{+jβx}` ([`fit_standing_wave`]). This separates the
/// forward (incident / transmitted) from the backward (reflected) travelling
/// wave, so the raw transmission is the **complex ratio** `b₂/a₁` (transmitted
/// forward over incident forward) — free of the standing-wave / over-unity
/// artifact that corrupted the single-point de-embeds (ADR-0129/0131). A second,
/// element-free *thru* solve normalizes it:
/// `S21(f) = (b₂/a₁)_dut / (b₂/a₁)_thru` (≈ 0 dB in-band, rolling off in the
/// stopband).
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
    // each measured frequency. Each entry is the **complex** raw transmission
    // `b₂/a₁` (transmitted-forward over incident-forward) from the 3-point
    // standing-wave de-embed (F2.3-g, ADR-0132).
    let dut = run_board_solve(ladder, substrate, cfg, freqs, true);
    let thru = run_board_solve(ladder, substrate, cfg, freqs, false);

    freqs
        .iter()
        .enumerate()
        .map(|(fi, &f)| {
            let t_dut = dut[fi];
            let t_thru = thru[fi];
            // S21 = (b₂/a₁)_dut / (b₂/a₁)_thru — the transmission relative to the
            // bare-line thru, with the forward/backward waves already separated
            // so the standing wave cannot corrupt the ratio (no over-unity). A
            // vanishing thru maps to 0 transmission rather than a divide-by-zero
            // blow-up. We return the magnitude (the gate is a |S21| shape check).
            let s21 = if t_thru.abs() > 0.0 {
                t_dut.div(t_thru).abs()
            } else {
                0.0
            };
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

/// Extend the synthesized board's signal trace by `lead_m` metres of straight
/// microstrip past **each** port, and move the ports to the new line ends
/// (F2.3-g, ADR-0132).
///
/// The lead-in gives the PEC-box standing-wave de-embed uniform line on each
/// side: the source and the two 3-point reference regions sit clear of the
/// element + original port discontinuities, and the line is long enough past
/// each port that a steady standing wave develops. The added segments are
/// full-width (`w_line`) copper rectangles at the signal-line `y`-band, read off
/// the existing ports (centre `y_sig`, width `w_line`). The substrate / ground
/// extent follows automatically from the widened bounding box.
fn lengthen_board(board: &mut yee_filter::LumpedBoard, lead_m: f64) {
    debug_assert!(board.layout.ports.len() >= 2, "board needs ≥ 2 ports");
    let y_sig = board.layout.ports[0].at.y;
    let w_line = board.layout.ports[0].width_m;
    let x_min = board.layout.bbox.min.x;
    let x_max = board.layout.bbox.max.x;
    let line_y0 = y_sig - w_line / 2.0;

    // Lead-in copper at each end, abutting the existing trace ends.
    board
        .layout
        .traces
        .push(Polygon::rect(x_min - lead_m, line_y0, lead_m, w_line));
    board
        .layout
        .traces
        .push(Polygon::rect(x_max, line_y0, lead_m, w_line));

    // Widen the bbox in x to cover the new copper.
    board.layout.bbox.min.x = x_min - lead_m;
    board.layout.bbox.max.x = x_max + lead_m;

    // Move the ports to the new line ends (the source / reference plane will be
    // placed relative to these).
    board.layout.ports[0].at = Point2::new(x_min - lead_m, y_sig);
    let last = board.layout.ports.len() - 1;
    board.layout.ports[last].at = Point2::new(x_max + lead_m, y_sig);
}

/// CW per-frequency steady-state solve of the board (DUT if `place_elements`,
/// else a bare thru line), returning the **complex raw transmission** `b₂/a₁`
/// (transmitted-forward over incident-forward) at each measured frequency from
/// the 3-point standing-wave de-embed (F2.3-g, ADR-0132).
///
/// The board's signal line is lengthened by [`LumpedSimConfig::lead_in_cells`]
/// cells past each port ([`lengthen_board`]) and voxelized once. Then for every
/// frequency `f` a fresh solver (cloned grid, zeroed fields) runs in a **plain
/// PEC box** ([`WalkingSkeletonSolver::new`] — no CPML; the stable
/// [`crate::run_line_eeff`] pattern, ADR-0108), driven by a Hann-ramped CW soft
/// `E_z` source sheet at the input. After the high-Q tanks + the standing wave
/// settle, the steady-state `E_z` phasor is single-bin DFT'd at **three
/// equally-spaced columns** in the input reference region and three in the
/// output reference region (downstream of the last element). Each triple is
/// fitted to `V(x) = a·e^{−jβx} + b·e^{+jβx}` ([`fit_standing_wave`]): the
/// input fit's forward amplitude is the incident `a₁`; the output fit's forward
/// amplitude is the transmitted `b₂`. The returned per-frequency value is the
/// complex ratio `b₂/a₁`; [`simulate_lumped_board`] then thru-normalizes it.
/// Factored out so the DUT and thru runs share identical geometry, grid, drive,
/// and step body.
fn run_board_solve(
    ladder: &LumpedLadder,
    substrate: &Substrate,
    cfg: &LumpedSimConfig,
    freqs: &[f64],
    place_elements: bool,
) -> Vec<Cplx> {
    // --- 1. Board (lengthened) + voxelize. ----------------------------------
    let mut board = lumped_board(ladder, substrate, cfg.footprint);
    // Lengthen the line so the source + the two 3-point standing-wave reference
    // regions sit on uniform line clear of the element / original-port
    // discontinuities, and a steady standing wave develops on the PEC box
    // (F2.3-g, ADR-0132).
    let lead_m = cfg.lead_in_cells as f64 * cfg.dx_m;
    lengthen_board(&mut board, lead_m);
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

    // --- 2. Source column + two 3-point reference regions (F2.3-g, ADR-0132).
    // PEC box (no CPML): the source sits a guard band in from the input PEC
    // wall, and the standing-wave reference regions live on the uniform lead-in
    // line on each side. `npml` is repurposed as the PEC end-wall guard band
    // (see `LumpedSimConfig::npml`); a small extra `guard` keeps the source +
    // the input reference region clear of the launch transient.
    let npml = cfg.npml;
    let guard = 3usize; // extra clear cells past the wall-guard band
    // Strip centre row (the signal-line band midpoint) and the substrate probe
    // depth (the `E_z` node just below the trace, where the quasi-TEM vertical
    // field is strongest — matches `run_line_eeff`'s `k_probe`).
    let j_strip = (j_lo + j_hi) / 2;
    let k_probe = k_top.saturating_sub(1).max(1);

    // Probe spacing `d` (cells) for the 3-point standing-wave fit. It only needs
    // to land `βd` in a well-conditioned range (≈ π/6 … π/2) — `β` itself is
    // extracted from the data — so an *approximate* λ_g from a coarse FR-4 ε_eff
    // guess is fine. λ_g_guess = c / (f̄·√ε_eff_guess) at the band-centre f̄;
    // d ≈ λ_g_guess / 12.
    let eps_eff_guess = 0.5 * (substrate.eps_r + 1.0); // ≈ 2.7 for FR-4 (ε_r 4.4)
    let f_centre = {
        let s: f64 = freqs.iter().copied().sum();
        s / freqs.len() as f64
    };
    let lambda_g_guess = C0_M_S / (f_centre * eps_eff_guess.sqrt());
    let probe_d = ((lambda_g_guess / 12.0 / dx).round() as usize).max(3);

    // Input source column: a guard band in from the input PEC wall.
    let src_i = (npml + guard).min(nx.saturating_sub(1));
    // Input reference region: three columns `in0, in0+d, in0+2d`, a few cells
    // downstream of the source (clear of the launch transient).
    let in0 = (src_i + guard).min(nx.saturating_sub(1));
    // Output reference region: three columns ending a guard band short of the
    // far PEC wall (downstream of the last element), `out0, out0+d, out0+2d`.
    let out_end = nx.saturating_sub(npml + guard + 1);
    let out0 = out_end.saturating_sub(2 * probe_d);
    let in_cols = [in0, in0 + probe_d, in0 + 2 * probe_d];
    let out_cols = [out0, out0 + probe_d, out0 + 2 * probe_d];
    let drive_cell = (src_i, j_strip, k_top);
    assert!(
        out_cols[0] > in_cols[2],
        "run_board_solve: input and output reference regions overlap \
         (in_cols={in_cols:?}, out_cols={out_cols:?}, d={probe_d}); lengthen the \
         board (lead_in_cells) or shrink npml/guard/d"
    );
    assert!(
        out_cols[2] < nx,
        "run_board_solve: output reference region runs past the grid \
         (out_cols={out_cols:?}, nx={nx}); shrink d or npml/guard"
    );

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

    // --- 4. Clean-launch calibration pre-pass (F2.3-h, ADR-0133). -----------
    // One short, time-gated Gaussian-pulse run on this geometry (no elements in
    // the path that the gated forward incident wave reaches) measures, at each
    // frequency, the numerical phase constant `β(f)` (to phase the directional
    // source) and the trustworthy time-gated incident forward amplitude
    // `a₁_gated(f)` — the `run_line_eeff` pattern (ADR-0108): the DFT window
    // closes before the first far-wall reflection returns, so the input region
    // sees a pure forward wave. The calibration is launch/line-only, so it is run
    // once per board geometry (DUT and thru share the same lead-in line).
    let cal = calibrate_launch(
        model.grid.clone(),
        cfg,
        freqs,
        drive_cell,
        in_cols,
        j_strip,
        k_probe,
        (j_lo, j_hi),
        k_top,
        probe_d,
        dx,
    );

    // --- 5. CW per-frequency PEC-box de-embed loop (F2.3-g/-h, ADR-0132/0133).
    // For each measured frequency, drive a Hann-ramped CW sinusoid on the PEC
    // line through the DIRECTIONAL two-column phased source (β from `cal`), let
    // the high-Q tanks + the standing wave settle, then read the steady-state
    // `E_z` phasor at the three output reference columns, fit them to a
    // forward/backward travelling-wave pair for the transmitted forward `b₂`, and
    // return the complex raw transmission `b₂/a₁_gated` (transmitted-forward over
    // the trustworthy time-gated incident-forward). (`ladder.z0_ohm` is unused —
    // there is no lumped load.)
    freqs
        .iter()
        .enumerate()
        .map(|(fi, &f)| {
            cw_deembed_b2_over_a1(
                model.grid.clone(),
                cfg,
                f,
                cal[fi].beta,
                cal[fi].a1_gated,
                drive_cell,
                in_cols,
                out_cols,
                j_strip,
                k_probe,
                (j_lo, j_hi),
                k_top,
                probe_d,
                dx,
                &recipes,
            )
        })
        .collect()
}

/// Per-frequency clean-launch calibration (F2.3-h, ADR-0133): the numerical
/// phase constant `β` (to phase the directional source) and the trustworthy
/// time-gated incident forward amplitude `a₁_gated`.
#[derive(Clone, Copy, Debug)]
struct LaunchCal {
    /// Numerical phase constant `β` (rad/m) at this frequency, read from the
    /// time-gated forward incident wave (self-consistent with the FDTD grid's
    /// numerical dispersion). Used to retard the directional source's downstream
    /// column by `Δφ = β·dx`.
    beta: f64,
    /// Time-gated incident **forward** amplitude `a₁` at the input reference
    /// region — a pure forward wave (the gate closes before the first far-wall
    /// reflection returns), so it is a trustworthy launch reference.
    a1_gated: Cplx,
}

/// Run one short, time-gated **Gaussian-pulse** solve per frequency on a fresh
/// PEC-box `grid` and return, for each, the numerical `β` and the time-gated
/// incident forward amplitude `a₁_gated` (F2.3-h, ADR-0133; the
/// [`crate::run_line_eeff`] time-gated incident-wave pattern, ADR-0108).
///
/// The directional two-column phased source (see [`inject_directional_source`])
/// launches a modulated-Gaussian forward pulse; the DFT at each frequency is
/// **time-gated** to the forward passage at the input reference region, closing
/// before the first far-wall reflection returns. With no reflected wave in the
/// window, the 3-point fit at the input reference region reads a **pure forward**
/// incident amplitude (`β > 0`, `|bwd| ≈ 0`). The same physical phasing as the CW
/// directional source is used, but the gate's `β` is bootstrapped from the
/// quasi-TEM ε_eff guess (the forward pulse is insensitive to the exact phasing —
/// any residual backward leak is gated out anyway), then refined from the data.
///
/// Each pulse run is short (a few hundred steps — the forward transit, not the
/// high-Q ring-up), so the calibration is cheap relative to the CW solves.
#[allow(clippy::too_many_arguments)]
fn calibrate_launch(
    grid: yee_fdtd::YeeGrid,
    cfg: &LumpedSimConfig,
    freqs: &[f64],
    drive_cell: (usize, usize, usize),
    in_cols: [usize; 3],
    j_strip: usize,
    k_probe: usize,
    strip_band: (usize, usize),
    k_top: usize,
    probe_d: usize,
    dx: f64,
) -> Vec<LaunchCal> {
    let dt = grid.dt;
    let dz = grid.dz;
    let (nx, _ny, _nz) = grid.ez.dim();
    let (j_lo, j_hi) = strip_band;
    let src_i = drive_cell.0;

    // ε_eff guess (quasi-TEM upper bound) → an initial β/v_p for the directional
    // source phasing and for sizing the time gate. The forward pulse is robust to
    // imperfect phasing (any backward leak is gated out), and the *measured* β is
    // returned for the CW source; this guess only bootstraps the gate length.
    // Use the highest swept frequency's λ_g as the conservative (shortest)
    // wavelength for the gate-length transit estimate.
    let eps_eff_guess = 4.4_f64; // FR-4 ε_r — upper bound on quasi-TEM ε_eff (β large → safe gate)
    let v_p_guess = C0_M_S / eps_eff_guess.sqrt();

    freqs
        .iter()
        .map(|&f| {
            // β guess for THIS frequency (rad/m): the source phasing target.
            let beta_guess = 2.0 * PI * f / v_p_guess;

            // Modulated-Gaussian pulse centred at `f`, ~80 % fractional bandwidth
            // (matches `run_line_eeff`): a clean forward launch whose tail is
            // negligible at t = 0.
            let bw = 0.8 * f;
            let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt()
                / (std::f64::consts::PI * bw))
                / dt)
                .ceil() as usize;
            let wave = SourceWaveform::GaussianPulse {
                v0: cfg.drive_v0,
                f0: f,
                bw,
                t0_steps,
            };

            // Time gate: stop ~10 % before the first far-wall reflection returns
            // to the downstream input-reference column. Reflection path =
            // source → far wall (≈ (nx − src_i)·dx) → back to in_cols[2]
            // (≈ (nx − in_cols[2])·dx), at v_p_guess. Run a little past the gate so
            // the forward pulse is fully integrated.
            let refl_dist = ((nx - src_i) as f64 + (nx - in_cols[2]) as f64) * dx;
            let t_refl = refl_dist / v_p_guess;
            let gate = (0.9 * t_refl / dt) as usize;
            // Ensure the pulse has time to reach the probes (t0 + a transit).
            let n_steps = gate + 200;

            let mut solver = WalkingSkeletonSolver::new(grid.clone());
            let omega = 2.0 * PI * f;
            let mut in_bins = [Bin::new(omega), Bin::new(omega), Bin::new(omega)];

            for n in 0..n_steps {
                solver.update_h_only();
                solver.apply_cpml_h(); // no CPML → PEC clamp

                inject_directional_source(
                    solver.grid_mut(),
                    &wave,
                    n,
                    dt,
                    src_i,
                    j_lo,
                    j_hi,
                    k_top,
                    beta_guess,
                    dx,
                    f,
                );

                solver.update_e_only();
                solver.apply_cpml_e();
                solver.advance_clock();

                if n < gate {
                    let g = solver.grid();
                    let t = n as f64 * dt;
                    for (idx, &i) in in_cols.iter().enumerate() {
                        let v = g.ez[(i, j_strip, k_probe)] * dz;
                        in_bins[idx].accumulate(v, t);
                    }
                }
            }

            let d = probe_d as f64 * dx;
            let fit = fit_standing_wave(
                in_bins[0].phasor(),
                in_bins[1].phasor(),
                in_bins[2].phasor(),
                d,
            );

            eprintln!(
                "[F2.3-h CAL] f={:.3} GHz | gate={gate} steps | β={:.2} rad/m \
                 (guess {:.2}) | |a₁_gated(fwd)|={:.3e} |bwd|={:.3e} (ratio {:.3}) | \
                 cos-resid={:.2e}",
                f * 1e-9,
                fit.beta,
                beta_guess,
                fit.fwd.abs(),
                fit.bwd.abs(),
                if fit.fwd.abs() > 0.0 {
                    fit.bwd.abs() / fit.fwd.abs()
                } else {
                    f64::INFINITY
                },
                fit.cos_imag_residual,
            );

            LaunchCal {
                beta: fit.beta,
                a1_gated: fit.fwd,
            }
        })
        .collect()
}

/// Inject the **directional two-column phased soft `E_z` source** (F2.3-h,
/// ADR-0133) into the strip's `(y, z)` face at the input.
///
/// Two soft `E_z` source sheets one cell apart in `x` (`src_i` and `src_i + 1`),
/// the downstream one **retarded** by the one-cell wave-transit phase
/// `Δφ = β·dx` (= `ω·dx/v_p`), add constructively in `+x` and destructively in
/// `−x` — a poor-man's Huygens / TF-SF launcher (ADR-0014/0021/0026) that injects
/// predominantly forward, so the PEC-box input is no longer a near-pure standing
/// wave and real forward power reaches the output region. The retardation is a
/// time shift of the same waveform: `s₂(t) = s(t − Δφ/ω)` for a CW tone; for the
/// pulse the same time shift retards the envelope and carrier together. `f` is
/// the carrier frequency (for the CW-tone time-shift); a non-positive `β` (the
/// degenerate calibration fallback) collapses to a single-column launch.
#[allow(clippy::too_many_arguments)]
fn inject_directional_source(
    grid: &mut yee_fdtd::YeeGrid,
    wave: &SourceWaveform,
    n: usize,
    dt: f64,
    src_i: usize,
    j_lo: usize,
    j_hi: usize,
    k_top: usize,
    beta: f64,
    dx: f64,
    f: f64,
) {
    let (nx, _ny, _nz) = grid.ez.dim();
    // Upstream column: the bare waveform.
    let s0 = wave.value(n, dt);
    for j in j_lo..j_hi {
        for k in 0..k_top {
            grid.ez[(src_i, j, k)] += s0;
        }
    }
    // Downstream column (src_i + 1): retard by the one-cell transit phase
    // Δφ = β·dx, i.e. a time delay Δt = Δφ/ω = β·dx/(2πf). A non-positive β
    // (degenerate fit) → single-column launch (no second column).
    let i2 = src_i + 1;
    if beta > 0.0 && i2 < nx && f > 0.0 {
        let dt_shift = beta * dx / (2.0 * PI * f);
        let n_shift = dt_shift / dt; // in (fractional) steps
        // Evaluate the same waveform at the retarded continuous time t − Δt by
        // using a fractional step index (SourceWaveform::value samples a
        // continuous t = n·dt). For n below the shift the retarded source is zero
        // (causal: the second sheet has not "seen" the wave yet).
        let n_eff = n as f64 - n_shift;
        let s1 = if n_eff >= 0.0 {
            wave_value_frac(wave, n_eff, dt)
        } else {
            0.0
        };
        for j in j_lo..j_hi {
            for k in 0..k_top {
                grid.ez[(i2, j, k)] += s1;
            }
        }
    }
}

/// Evaluate a [`SourceWaveform`] at a **fractional** step index `n_frac` (the
/// continuous time `t = n_frac·dt`), for the directional source's sub-step
/// retardation (F2.3-h, ADR-0133).
///
/// [`SourceWaveform::value`] takes an integer step; the directional launcher
/// needs the same waveform at `t − Δt` where `Δt = β·dx/ω` is generally not an
/// integer number of steps. Both `HannSine` and `GaussianPulse` are closed-form
/// in continuous `t`, so a fractional `n_frac` simply samples that closed form
/// (the ramp / envelope use `n_frac` in place of the integer step).
fn wave_value_frac(wave: &SourceWaveform, n_frac: f64, dt: f64) -> f64 {
    let t = n_frac * dt;
    match *wave {
        SourceWaveform::None => 0.0,
        SourceWaveform::HannSine {
            v0,
            frequency,
            ramp_steps,
        } => {
            let ramp = if ramp_steps == 0 || n_frac >= ramp_steps as f64 {
                1.0
            } else {
                0.5 * (1.0 - (PI * n_frac / ramp_steps as f64).cos())
            };
            v0 * ramp * (2.0 * PI * frequency * t).sin()
        }
        SourceWaveform::GaussianPulse {
            v0,
            f0,
            bw,
            t0_steps,
        } => {
            let t0 = t0_steps as f64 * dt;
            let tau = if bw > 0.0 {
                (2.0 * std::f64::consts::LN_2).sqrt() / (PI * bw)
            } else {
                f64::INFINITY
            };
            let env = if tau.is_infinite() {
                1.0
            } else {
                let arg = (t - t0) / tau;
                (-arg * arg).exp()
            };
            v0 * env * (2.0 * PI * f0 * (t - t0)).sin()
        }
    }
}

/// Run one CW steady-state FDTD solve at frequency `f` on a fresh `grid` (zeroed
/// fields) in a **plain PEC box** and return the complex raw transmission
/// `b₂/a₁` (F2.3-g/-h, ADR-0132/0133).
///
/// The line runs in a hard-PEC box ([`WalkingSkeletonSolver::new`] — the
/// `apply_cpml_*` calls fall back to a PEC clamp when no CPML is configured;
/// this is the *stable* [`crate::run_line_eeff`] pattern, ADR-0108: CPML into a
/// microstrip's PEC-ground / high-ε substrate is late-time unstable). A
/// Hann-ramped CW **soft** `E_z` source sheet, launched through the
/// **directional two-column phased** launcher ([`inject_directional_source`],
/// `β` from the calibration pre-pass) at `drive_cell`'s column, sends a
/// predominantly `+x` travelling quasi-TEM wave; after
/// `cw_ramp_cycles + cw_settle_cycles` the fields settle into a single-frequency
/// steady standing wave (the highest-Q tank ring-up + the box transits).
///
/// Over the final `cw_measure_cycles` (the settled window) the steady-state
/// `E_z·dz` voltage phasor is single-bin DFT'd at the three `in_cols` (a
/// diagnostic) and three `out_cols` columns (substrate depth `k_probe`, strip
/// centre `j_strip`). The output triple of spacing `probe_d·dx` is fitted to
/// `V(x) = a·e^{−jβx} + b·e^{+jβx}` ([`fit_standing_wave`]) → the transmitted
/// forward amplitude `b₂`. The incident forward `a₁` is the **time-gated**
/// `a1_gated` from the calibration pre-pass (F2.3-h, ADR-0133) — a trustworthy
/// pure-forward launch reference, NOT the CW input fit (which the input-wall
/// reflection contaminates). The returned complex ratio `b₂/a₁_gated` is the
/// transmitted-forward over incident-forward, free of the standing-wave /
/// over-unity artifact (ADR-0129/0131) and well above the floor (ADR-0132).
///
/// The per-step body mirrors [`crate::run_line_eeff`]: `update_h_only` →
/// `apply_cpml_h` (= PEC clamp) → directional soft CW source → `update_e_only` →
/// `apply_cpml_e` (= PEC clamp + interior mask) → the filter elements'
/// multi-cell `correct_e_aperture` → advance the clock and record.
#[allow(clippy::too_many_arguments)]
fn cw_deembed_b2_over_a1(
    grid: yee_fdtd::YeeGrid,
    cfg: &LumpedSimConfig,
    f: f64,
    beta: f64,
    a1_gated: Cplx,
    drive_cell: (usize, usize, usize),
    in_cols: [usize; 3],
    out_cols: [usize; 3],
    j_strip: usize,
    k_probe: usize,
    strip_band: (usize, usize),
    k_top: usize,
    probe_d: usize,
    dx: f64,
    recipes: &[ElementRecipe],
) -> Cplx {
    let dt = grid.dt;
    let dz = grid.dz;

    // Plain PEC box (no CPML): the stable run_line_eeff configuration. The
    // `apply_cpml_*` calls below fall back to a hard-PEC tangential-E clamp when
    // no CPML state is configured (ADR-0108).
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

    // Hann-ramped CW sinusoid at `f`, injected as a SOFT `E_z` source over the
    // strip's `(y, z)` face at the source column (the quasi-TEM launcher).
    let wave = SourceWaveform::HannSine {
        v0: cfg.drive_v0,
        frequency: f,
        ramp_steps,
    };
    let (j_lo, j_hi) = strip_band;
    let src_i = drive_cell.0;

    // Fresh passive aperture ports for this solve.
    let mut elements: Vec<LumpedRlcPort> = recipes.iter().flat_map(ElementRecipe::build).collect();

    // One single-bin DFT accumulator per reference column (3 input + 3 output),
    // accumulated over the SETTLED window only → the steady-state phasors.
    let omega = 2.0 * PI * f;
    let mut in_bins = [Bin::new(omega), Bin::new(omega), Bin::new(omega)];
    let mut out_bins = [Bin::new(omega), Bin::new(omega), Bin::new(omega)];

    for n in 0..n_steps {
        solver.update_h_only();
        solver.apply_cpml_h(); // no CPML → PEC clamp

        // Directional two-column phased soft `E_z` source sheet across the strip
        // face at the input (F2.3-h, ADR-0133): the downstream column is retarded
        // by the one-cell transit phase `Δφ = β·dx` so the launch adds forward and
        // cancels backward — predominantly forward injection (no near-pure input
        // standing wave). `β` from the time-gated calibration pre-pass.
        inject_directional_source(
            solver.grid_mut(),
            &wave,
            n,
            dt,
            src_i,
            j_lo,
            j_hi,
            k_top,
            beta,
            dx,
            f,
        );

        solver.update_e_only();
        solver.apply_cpml_e(); // no CPML → PEC clamp + interior PEC mask

        for el in elements.iter_mut() {
            el.correct_e_aperture(solver.grid_mut(), n, dt);
        }

        solver.advance_clock();

        if n >= measure_start {
            let grid = solver.grid();
            let t = n as f64 * dt;
            for (idx, &i) in in_cols.iter().enumerate() {
                let v = grid.ez[(i, j_strip, k_probe)] * dz;
                in_bins[idx].accumulate(v, t);
            }
            for (idx, &i) in out_cols.iter().enumerate() {
                let v = grid.ez[(i, j_strip, k_probe)] * dz;
                out_bins[idx].accumulate(v, t);
            }
        }
    }

    // Fit each triple to V(x) = a·e^{−jβx} + b·e^{+jβx}. The probe columns are
    // equally spaced by `probe_d` cells → physical spacing `d = probe_d·dx`.
    let d = probe_d as f64 * dx;
    let in_fit = fit_standing_wave(
        in_bins[0].phasor(),
        in_bins[1].phasor(),
        in_bins[2].phasor(),
        d,
    );
    let out_fit = fit_standing_wave(
        out_bins[0].phasor(),
        out_bins[1].phasor(),
        out_bins[2].phasor(),
        d,
    );

    // Incident forward `a₁`: the TRUSTWORTHY time-gated reference from the
    // calibration pre-pass (F2.3-h, ADR-0133), NOT the CW input standing-wave fit
    // (which the input-wall reflection contaminates). `a₁_gated` is a launch/line
    // property — identical for the DUT and thru runs — so it cancels in the
    // thru-normalization, but referencing `b₂` to it makes the raw transmission a
    // clean transmitted-forward-over-incident-forward ratio at a well-resolved
    // amplitude. The CW input fit (`in_fit`) is kept only as a diagnostic
    // (β_in / the residual reflection) to confirm the directional launch worked.
    let a1 = a1_gated;
    // b₂ = transmitted forward (output region, CW-settled so the high-Q tanks
    // have rung up). Raw transmission is the complex ratio b₂/a₁.
    let b2 = out_fit.fwd;
    let t_raw = if a1.abs() > 0.0 {
        b2.div(a1)
    } else {
        Cplx::ZERO
    };

    eprintln!(
        "[F2.3-h DIAG] f={:.3} GHz | β_in(cw)={:.2} β_out={:.2} rad/m (d={:.2} mm) | \
         in(cw): |fwd|={:.3e} |bwd(refl)|={:.3e} (refl/fwd {:.2}) | a₁_gated={:.3e} | \
         out: |fwd b₂|={:.3e} |bwd(wall)|={:.3e} | |b₂/a₁|={:.4} | \
         cos-resid in={:.2e} out={:.2e}",
        f * 1e-9,
        in_fit.beta,
        out_fit.beta,
        d * 1e3,
        in_fit.fwd.abs(),
        in_fit.bwd.abs(),
        if in_fit.fwd.abs() > 0.0 {
            in_fit.bwd.abs() / in_fit.fwd.abs()
        } else {
            f64::INFINITY
        },
        a1.abs(),
        b2.abs(),
        out_fit.bwd.abs(),
        t_raw.abs(),
        in_fit.cos_imag_residual,
        out_fit.cos_imag_residual,
    );

    t_raw
}
