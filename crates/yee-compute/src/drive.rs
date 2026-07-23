//! Sources, lumped ports, and probes (E.2).
//!
//! The waveforms and the resistive-port update are verbatim ports of
//! `yee_fdtd::sources::gaussian_pulse_ez` and
//! `yee_fdtd::lumped::{SourceWaveform::GaussianPulse, update_pure_resistor}`
//! — gate `compute-007` asserts the driven CPU step is bit-exact against
//! `WalkingSkeletonSolver` + `LumpedRlcPort`, so keep the arithmetic
//! byte-identical to the reference.

use std::f64::consts::{PI, TAU};

use crate::spec::{FdtdSpec, idx3, len3};

/// Electric-field component selector for sources and probes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EComponent {
    /// `E_x`
    Ex,
    /// `E_y`
    Ey,
    /// `E_z`
    Ez,
}

impl EComponent {
    /// Dims of this component's array under `spec`.
    pub(crate) fn dims(self, spec: &FdtdSpec) -> (usize, usize, usize) {
        match self {
            Self::Ex => spec.ex_dims(),
            Self::Ey => spec.ey_dims(),
            Self::Ez => spec.ez_dims(),
        }
    }

    /// Flat index of `cell` within this component's array.
    pub(crate) fn flat(self, spec: &FdtdSpec, cell: (usize, usize, usize)) -> usize {
        let dims = self.dims(spec);
        assert!(
            cell.0 < dims.0 && cell.1 < dims.1 && cell.2 < dims.2,
            "cell {cell:?} out of bounds for {self:?} dims {dims:?}"
        );
        idx3(dims, cell.0, cell.1, cell.2)
    }

    /// Offset of this component within the packed E/H field arena
    /// (ex | ey | ez | …), matching the GPU packing.
    #[cfg_attr(not(feature = "gpu"), allow(dead_code))]
    pub(crate) fn arena_offset(self, spec: &FdtdSpec) -> usize {
        match self {
            Self::Ex => 0,
            Self::Ey => len3(spec.ex_dims()),
            Self::Ez => len3(spec.ex_dims()) + len3(spec.ey_dims()),
        }
    }
}

/// Time-domain waveform, sampled at step boundaries (`t = n·dt`).
#[derive(Debug, Clone, Copy)]
pub enum Waveform {
    /// Plain Gaussian `exp(−((t−t0)/σ)²)` — the
    /// `sources::gaussian_pulse_ez` shape used by the cavity and CPML gates.
    Gaussian {
        /// Pulse centre (seconds).
        t0: f64,
        /// Pulse width (seconds).
        sigma: f64,
    },
    /// Modulated Gaussian `v0·exp(−((t−t0)/τ)²)·sin(2πf0(t−t0))` with
    /// `τ = √(2 ln 2)/(π·bw)` and `t0 = t0_steps·dt` — verbatim
    /// `SourceWaveform::GaussianPulse` (the line-eeff launch).
    GaussianPulse {
        /// Peak EMF (volts).
        v0: f64,
        /// Carrier frequency (Hz).
        f0: f64,
        /// FWHM spectral bandwidth (Hz); `0` degenerates to pure CW.
        bw: f64,
        /// Pulse centre in steps.
        t0_steps: usize,
    },
}

impl Waveform {
    /// Sample the waveform at step `n_step` (`t = n_step·dt`), matching the
    /// reference implementations bit-for-bit.
    pub fn value(&self, n_step: usize, dt: f64) -> f64 {
        let t = n_step as f64 * dt;
        match *self {
            Self::Gaussian { t0, sigma } => {
                let arg = (t - t0) / sigma;
                (-arg * arg).exp()
            }
            Self::GaussianPulse {
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
                let arg = (t - t0) / tau;
                let env = if tau.is_infinite() {
                    1.0
                } else {
                    (-arg * arg).exp()
                };
                v0 * env * (TAU * f0 * (t - t0)).sin()
            }
        }
    }
}

/// Additive ("soft") point source: `field[cell] += waveform(t)` between the
/// H and E half-steps, exactly where the reference solvers inject.
#[derive(Debug, Clone, Copy)]
pub struct SoftSource {
    /// Driven component.
    pub component: EComponent,
    /// Driven cell `(i, j, k)` in the component's own staggered indexing.
    pub cell: (usize, usize, usize),
    /// Amplitude waveform.
    pub waveform: Waveform,
}

/// Lumped resistive drive port on a single `E_z` cell — the pure-resistor
/// arm of `yee_fdtd::lumped::LumpedRlcPort` (the validated launch used by
/// `fdtd-line-eeff-001`). Applied after the E half-step and boundary phase:
///
/// ```text
/// α = Δt·dz / (2·ε₀·R·dx·dy),  γ = Δt / (ε₀·R·dx·dy)
/// E¹ = (E¹* − α·E⁰ + γ·V_src(n)) / (1 + α)
/// ```
#[derive(Debug, Clone, Copy)]
pub struct ResistivePort {
    /// Port `E_z` cell.
    pub cell: (usize, usize, usize),
    /// Port resistance (Ω). Must be positive and finite.
    pub resistance: f64,
    /// Thevenin EMF waveform.
    pub waveform: Waveform,
}

/// Multi-cell **aperture** resistive port (S.10, ADR-0187) — a verbatim
/// port of the pure-R arm of `yee_fdtd::LumpedRlcPort::aperture` /
/// `correct_e_aperture` (Phase 2.fdtd.6.9, ADR-0125): one aggregate R
/// branch bridging the modal port face. The branch voltage is the modal
/// `V = ∫E_z·dz` over the full substrate height averaged across the width
/// columns, the semi-implicit two-way solve uses the aperture back-action
/// impedance `β = dt·h/(2·ε₀·A)`, and the branch current is distributed
/// back onto every aperture cell as a sheet current referenced to the
/// **physical** area `A` — the dx-stable lumped port a single-cell
/// [`ResistivePort`] cannot approximate on a multi-cell substrate.
#[derive(Debug, Clone)]
pub struct AperturePort {
    /// The `E_z` cells `(i, j, k)` spanning the `(y, z)` aperture face
    /// (trace width × substrate height), all sharing the port-plane `i`.
    pub cells: Vec<(usize, usize, usize)>,
    /// Number of width-direction (`y`) columns; the modal voltage averages
    /// the per-column height integral over these.
    pub n_columns: usize,
    /// Physical aperture area `A = w·h` (m²).
    pub area: f64,
    /// Physical substrate height `h` (m) the modal voltage integrates over.
    pub height: f64,
    /// Aggregate branch resistance (Ω); `f64::INFINITY` = open.
    pub resistance: f64,
    /// Series EMF waveform (`v0 = 0` for a passive matched load).
    pub waveform: Waveform,
    /// Record the per-step `(v_src, v_terminal, i_branch)` triple
    /// (FS.2a, ADR-0207) — the accepted-power observables. CPU-only for
    /// now; the GPU backend rejects recording ports (`Unsupported`).
    pub record: bool,
}

/// Per-step field sample recorded after each full step.
#[derive(Debug, Clone, Copy)]
pub struct Probe {
    /// Sampled component.
    pub component: EComponent,
    /// Sampled cell.
    pub cell: (usize, usize, usize),
}

/// Magnetic-field component selector for [`HProbe`] (FS.4.2a).
///
/// Kept as its own type — not folded into [`EComponent`] — so no existing
/// `EComponent`/`Probe` call site (source, port, or E-probe) has to change:
/// H sampling lives entirely on the parallel [`Drive::h_probes`] field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HComponent {
    /// `H_x`
    Hx,
    /// `H_y`
    Hy,
    /// `H_z`
    Hz,
}

impl HComponent {
    /// Dims of this component's array under `spec`.
    fn dims(self, spec: &FdtdSpec) -> (usize, usize, usize) {
        match self {
            Self::Hx => spec.hx_dims(),
            Self::Hy => spec.hy_dims(),
            Self::Hz => spec.hz_dims(),
        }
    }

    /// Flat index of `cell` within this component's array.
    pub(crate) fn flat(self, spec: &FdtdSpec, cell: (usize, usize, usize)) -> usize {
        let dims = self.dims(spec);
        assert!(
            cell.0 < dims.0 && cell.1 < dims.1 && cell.2 < dims.2,
            "cell {cell:?} out of bounds for {self:?} dims {dims:?}"
        );
        idx3(dims, cell.0, cell.1, cell.2)
    }

    /// Offset of this component within the packed field arena, counting
    /// past the E block (`ex | ey | ez |` **`hx | hy | hz`**) — matches the
    /// GPU packing order in `gpu.rs` / `off_hx`/`off_hy`/`off_hz` in
    /// `fdtd.wgsl`.
    #[cfg_attr(not(feature = "gpu"), allow(dead_code))]
    pub(crate) fn arena_offset(self, spec: &FdtdSpec) -> usize {
        let e_len = len3(spec.ex_dims()) + len3(spec.ey_dims()) + len3(spec.ez_dims());
        match self {
            Self::Hx => e_len,
            Self::Hy => e_len + len3(spec.hx_dims()),
            Self::Hz => e_len + len3(spec.hx_dims()) + len3(spec.hy_dims()),
        }
    }
}

/// Per-step H-field sample recorded after each full step (FS.4.2a) —
/// parallel to [`Probe`] on [`Drive::h_probes`]. Recorded alongside the E
/// probes (after the port/aperture-port correction phase, same step
/// iteration): the H state read there is `H` at `t = (n + ½)·Δt` — the
/// half-step **before** the co-recorded `E` sample at `t = (n + 1)·Δt`,
/// because `update_h` runs once at the top of the leapfrog iteration and
/// nothing after it touches H. Callers doing Ampère-loop current extraction
/// (e.g. the stripline Z₀ gate) must account for this half-step offset from
/// the co-recorded V(E) sample explicitly — it is not corrected here.
#[derive(Debug, Clone, Copy)]
pub struct HProbe {
    /// Sampled component.
    pub component: HComponent,
    /// Sampled cell.
    pub cell: (usize, usize, usize),
}

/// Drive configuration: what excites the grid and what gets recorded.
#[derive(Debug, Clone, Default)]
pub struct Drive {
    /// Additive point sources (injected between the H and E half-steps).
    pub soft_sources: Vec<SoftSource>,
    /// Resistive drive ports (applied after the E boundary phase).
    pub ports: Vec<ResistivePort>,
    /// Multi-cell aperture ports (applied after the single-cell ports;
    /// CPU backend only — the GPU backend rejects a drive that carries
    /// any, see `ComputeError::Unsupported`).
    pub aperture_ports: Vec<AperturePort>,
    /// Probes recorded once per step, in order, after the ports.
    pub probes: Vec<Probe>,
    /// H-field probes (FS.4.2a) recorded once per step, alongside `probes`
    /// — see [`HProbe`] for the staggering caveat.
    pub h_probes: Vec<HProbe>,
}

impl Drive {
    /// Panic unless every source/port/probe cell is in bounds.
    pub(crate) fn validate(&self, spec: &FdtdSpec) {
        for s in &self.soft_sources {
            let _ = s.component.flat(spec, s.cell);
        }
        for p in &self.ports {
            assert!(
                p.resistance > 0.0 && p.resistance.is_finite(),
                "ResistivePort: resistance must be positive and finite (got {})",
                p.resistance
            );
            let _ = EComponent::Ez.flat(spec, p.cell);
        }
        for p in &self.aperture_ports {
            assert!(
                (p.resistance > 0.0 && p.resistance.is_finite()) || p.resistance.is_infinite(),
                "AperturePort: resistance must be positive (got {}); use f64::INFINITY for open",
                p.resistance
            );
            assert!(!p.cells.is_empty(), "AperturePort: no cells");
            assert!(p.n_columns >= 1, "AperturePort: n_columns must be >= 1");
            assert!(
                p.area.is_finite() && p.area > 0.0,
                "AperturePort: area must be finite and positive (got {})",
                p.area
            );
            assert!(
                p.height.is_finite() && p.height > 0.0,
                "AperturePort: height must be finite and positive (got {})",
                p.height
            );
            for &cell in &p.cells {
                let _ = EComponent::Ez.flat(spec, cell);
            }
        }
        for p in &self.probes {
            let _ = p.component.flat(spec, p.cell);
        }
        for p in &self.h_probes {
            let _ = p.component.flat(spec, p.cell);
        }
    }

    /// True when there is nothing to inject, correct, or record.
    #[cfg_attr(not(feature = "gpu"), allow(dead_code))]
    pub(crate) fn is_empty(&self) -> bool {
        self.soft_sources.is_empty()
            && self.ports.is_empty()
            && self.aperture_ports.is_empty()
            && self.probes.is_empty()
            && self.h_probes.is_empty()
    }
}
