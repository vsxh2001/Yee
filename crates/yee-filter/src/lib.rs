//! # yee-filter
//!
//! Filter-domain data model, ideal-response evaluation, and spec-mask gating
//! for the Yee electromagnetic-simulation studio (Filter Phase F0).
//!
//! This crate is the **data model every later filter phase plugs into**. It
//! consumes [`yee_synth`] (prototype g-values, coupling matrix, external Q) and
//! exposes the end-to-end pipe:
//! `FilterSpec → synthesize → FilterProject → ideal_response / check_mask`.
//!
//! The ideal response uses the **closed-form** lowpass transfer function
//! applied to the bandpass-mapped lowpass variable `Ω` (spec §2.4):
//! Chebyshev `|S21|² = 1/(1 + ε²·T_N²(Ω))` with `ε = √(10^{L_Ar/10} − 1)`, and
//! Butterworth `|S21|² = 1/(1 + Ω^{2N})`. Reflection follows from losslessness
//! (`|S11|² = 1 − |S21|²`). Driving S-parameters *from* the coupling matrix
//! (Hong-Lancaster `[A] = [q] + pU − jM`) is later-phase work — see ADR-0084.
//!
//! ## Example
//!
//! ```
//! use yee_filter::{Approximation, FilterSpec, Response, SpecMask, synthesize, check_mask};
//!
//! let spec = FilterSpec {
//!     response: Response::Bandpass,
//!     approximation: Approximation::Chebyshev { ripple_db: 0.5 },
//!     f0_hz: 2.0e9,
//!     fbw: 0.10,
//!     order: Some(5),
//!     z0_ohm: 50.0,
//!     mask: SpecMask {
//!         passband_ripple_db: 0.5,
//!         return_loss_db: 10.0,
//!         stopband: vec![(2.4e9, 30.0)],
//!     },
//! };
//! let proj = synthesize(&spec);
//! assert_eq!(proj.prototype.order(), 5);
//! ```

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Coupling/Qe extraction (Filter Phase F1.1b.0): measured response → `k`/`Q`.
pub mod extract;
pub use extract::{CouplingExtraction, extract_coupling, extract_q_ringdown};

/// Closed-form dimensional synthesis (Filter Phases F1.2.0 / F1.2.2 / F1.2.3 /
/// F1.2.5 / F1.2.7): coupling matrix → physical microstrip dimensions
/// (edge-coupled, hairpin, combline, interdigital) and low-pass prototype →
/// stepped-impedance line sections.
pub mod dimension;
pub use dimension::{
    ComblineDimensions, DimError, EdgeCoupledDimensions, HairpinDimensions, InterdigitalDimensions,
    SteppedImpedanceDimensions, SteppedSection, dimension_combline, dimension_combline_layout,
    dimension_edge_coupled, dimension_edge_coupled_layout, dimension_hairpin,
    dimension_hairpin_layout, dimension_interdigital, dimension_interdigital_layout,
    dimension_stepped_impedance, dimension_stepped_impedance_layout,
};

/// Closed-form lumped-element LC ladder synthesis (Filter Phase F2.0):
/// prototype g-values → ideal series/shunt LC resonators for a band-pass
/// filter.
pub mod lumped;
pub use lumped::{
    LcBranch, LcResonator, LumpedError, LumpedLadder, MaskVerdict, mask_verdict, synthesize_lumped,
};
// `ladder_s21` / `ladder_s21_lossy` / `ladder_s_params_lossy` are
// `#[doc(hidden)] pub`: the realized-response ABCD helpers (lossless and
// finite-Q S21, plus the full lossy `(S11, S21)` pair), kept out of the
// documented API surface but reachable by the `lumped_001` / `lumped_q_001`
// gates and the CLI's finite-Q Touchstone export.
#[doc(hidden)]
pub use lumped::{ladder_s_params_lossy, ladder_s21, ladder_s21_lossy};

/// Top-C-coupled (capacitively-coupled) band-pass synthesis (JLCPCB narrow-band
/// track, ADR-0165 brick T1): low-pass prototype → `N` shunt LC resonators +
/// `N+1` series coupling capacitors (admittance-inverter coupled resonators).
pub mod top_c;
pub use top_c::{ShuntResonator, TopCNetwork, synthesize_top_c_coupled};
// `top_c_s21` is `#[doc(hidden)] pub`: the realized-response ABCD helper, kept
// out of the documented API surface but reachable by the `top_c_coupled_001`
// gate (a separate crate) for the non-circular S21-mask validation.
#[doc(hidden)]
pub use top_c::top_c_s21;

/// Monte-Carlo tolerance / yield analysis (Filter Phase F2.4): snap each L/C to
/// an E-series value, perturb within tolerance over many seeded trials, and
/// report the fraction of realized ladders that meet the spec mask.
pub mod tolerance;
pub use tolerance::{YieldResult, monte_carlo_yield};

/// E-series component selection + bill of materials (Filter Phase F2.1): ideal
/// LC ladder values → nearest IEC 60063 standard parts + a grouped [`Bom`].
pub mod parts;
pub use parts::{Bom, BomLine, CompKind, ESeries, select_components};

/// Lumped-LC PCB board generator (Filter Phase F2.2): place an LC ladder's
/// resonators as SMD footprints + pads + traces on a [`yee_layout::Layout`].
pub mod board;
pub use board::{
    BranchKind, Footprint, LumpedBoard, PadSpec, Placement, lumped_board, top_c_board,
};

/// LCSC part autopick + bundled real-parts table (JLCPCB production track,
/// ADR-0164 brick J1): map an E-series [`BomLine`] to a real, orderable JLCPCB
/// **Basic** LCSC part by kind + [`Footprint`] + value, preferring Basic.
pub mod jlcpcb;
pub use jlcpcb::{
    DEFAULT_TOLERANCE_PCT, LCSC_PARTS, LcscPart, autopick, autopick_bom, autopick_within,
};

/// JLCPCB assembly upload CSV export (JLCPCB production track, ADR-0164 bricks
/// J2 + J3): a placed [`LumpedBoard`]'s [`Placement`]s + its [`LumpedLadder`]
/// values → the JLCPCB BOM CSV (`Comment,Designator,Footprint,LCSC Part #`) and
/// CPL/centroid CSV (`Designator,Mid X,Mid Y,Layer,Rotation`).
pub mod jlcpcb_export;
pub use jlcpcb_export::{
    BOM_HEADER, CPL_HEADER, JlcpcbFiles, PlacedPart, jlcpcb_bom_csv, jlcpcb_cpl_csv, jlcpcb_files,
    jlcpcb_footprint_name, join_placed_parts, join_top_c_parts, value_comment,
};

/// Guided technique-recommender (App.2.0, ADR-0136): a deterministic decision
/// tree mapping a [`FilterSpec`] to a recommended physical realization
/// technique with a plain-language rationale + ranked alternatives.
pub mod recommend;
pub use recommend::{RealizationTechnique, TechniqueRecommendation, recommend_technique};

/// Orderable board-topology auto-selector (JLCPCB production track, ADR-0167
/// brick T3): for a [`FilterProject`], pick the lumped board topology that
/// yields a fully-orderable JLCPCB board (alternating ladder for wideband,
/// top-C-coupled for the sub-GHz / moderate-band corner), or honestly report
/// that neither lumped topology can (the distributed/planar track). The board-
/// realization [`BoardTopology`] is distinct from the synthesis-realization
/// [`Topology`] enum — see the module docs.
pub mod topology;
pub use topology::{BoardTopology, OrderableBoard, synthesize_orderable, synthesize_orderable_on};

pub use yee_synth::Approximation;
use yee_synth::{Prototype, coupling_design, lowpass_to_bandpass, min_order, prototype};

/// Filter frequency-response class.
///
/// Phase F0 synthesizes and evaluates [`Response::Bandpass`]; the other
/// variants are reserved for later phases.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Response {
    /// Lowpass.
    Lowpass,
    /// Highpass.
    Highpass,
    /// Bandpass.
    Bandpass,
    /// Bandstop.
    Bandstop,
}

/// Passband + stopband requirements the synthesized response must satisfy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpecMask {
    /// Maximum allowed passband insertion-loss ripple, dB (e.g. `0.5`).
    pub passband_ripple_db: f64,
    /// Minimum required in-band return loss, dB (e.g. `10.0`). A larger value
    /// is a stricter match requirement.
    pub return_loss_db: f64,
    /// Stopband points: `(frequency_hz, minimum_rejection_db)`. Each point
    /// requires `|S21|` at that frequency to be at least `min_reject` dB down.
    pub stopband: Vec<(f64, f64)>,
}

/// A filter specification: the design intent the synthesis consumes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FilterSpec {
    /// Response class (Phase F0 synthesizes [`Response::Bandpass`]).
    pub response: Response,
    /// Approximation (response shape), re-exported from [`yee_synth`].
    pub approximation: Approximation,
    /// Centre frequency, Hz.
    pub f0_hz: f64,
    /// Fractional bandwidth `(f2 − f1) / f0`.
    pub fbw: f64,
    /// Explicit filter order. `None` → estimate the minimum order from the
    /// stopband mask via [`yee_synth::min_order`].
    pub order: Option<usize>,
    /// System reference impedance, Ω (written into the Touchstone option line).
    pub z0_ohm: f64,
    /// The spec mask the response is graded against by [`check_mask`].
    pub mask: SpecMask,
}

/// The synthesized normalized coupling matrix plus external Q values.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CouplingMatrix {
    /// Normalized `N × N` coupling matrix (synchronous → zero diagonal).
    pub m: Vec<Vec<f64>>,
    /// Input external quality factor.
    pub qe_in: f64,
    /// Output external quality factor.
    pub qe_out: f64,
}

/// Realization topology of the synthesized filter.
///
/// Phase F0 synthesizes a [`Topology::CoupledResonator`] all-pole network; the
/// `non_exhaustive` marker reserves room for `Ladder`, `Iris`, … in later
/// phases.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Topology {
    /// All-pole coupled-resonator topology (coupling matrix + external Q).
    CoupledResonator,
}

/// The persisted design document: spec + synthesized prototype, coupling
/// matrix, and topology.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FilterProject {
    /// The originating specification.
    pub spec: FilterSpec,
    /// The synthesized lowpass prototype (g-values).
    pub prototype: Prototype,
    /// The synthesized coupling matrix + external Q.
    pub coupling: CouplingMatrix,
    /// The realization topology.
    pub topology: Topology,
}

/// Synthesize a [`FilterProject`] from a [`FilterSpec`] (spec §3, stages 1–2).
///
/// Resolves the order (explicit, or estimated from the worst-case stopband
/// point via [`yee_synth::min_order`]), builds the lowpass prototype, and
/// synthesizes the all-pole coupling matrix + external Q at the spec's
/// fractional bandwidth.
///
/// # Panics
///
/// Panics if neither `spec.order` is set nor a stopband point is available to
/// estimate the order from, or if `spec.fbw <= 0.0`.
pub fn synthesize(spec: &FilterSpec) -> FilterProject {
    let order = match spec.order {
        Some(n) => n,
        None => estimate_order(spec),
    };
    let proto = prototype(spec.approximation, order);
    let design = coupling_design(&proto, spec.fbw);
    let coupling = CouplingMatrix {
        m: design.m,
        qe_in: design.qe_in,
        qe_out: design.qe_out,
    };
    FilterProject {
        spec: spec.clone(),
        prototype: proto,
        coupling,
        topology: Topology::CoupledResonator,
    }
}

/// Estimate the minimum order from the spec mask: take the stopband point with
/// the largest required rejection mapped through the bandpass transform, and
/// ask [`yee_synth::min_order`] for the order meeting it.
fn estimate_order(spec: &FilterSpec) -> usize {
    assert!(
        !spec.mask.stopband.is_empty(),
        "cannot estimate order: spec.order is None and the mask has no stopband points"
    );
    let omega0 = spec.f0_hz;
    let mut best = 1usize;
    for &(f_hz, reject_db) in &spec.mask.stopband {
        // Map the stopband frequency to the lowpass variable Ω; the lowpass
        // stopband ratio is |Ω| (the passband edge is Ω = ±1).
        let omega_s = lowpass_to_bandpass(f_hz, omega0, spec.fbw).abs();
        if omega_s <= 1.0 {
            // Inside the (mapped) passband — cannot constrain the order; skip.
            continue;
        }
        let n = min_order(spec.approximation, reject_db, omega_s);
        best = best.max(n);
    }
    best
}

/// Evaluate the ideal forward transmission `S21` over `freqs_hz` (spec §3).
///
/// Uses the closed-form lowpass transfer function evaluated at the
/// bandpass-mapped lowpass variable `Ω`:
///
/// - Chebyshev: `|S21|² = 1 / (1 + ε²·T_N²(Ω))`, `ε = √(10^{L_Ar/10} − 1)`,
///   `T_N` the order-`N` Chebyshev polynomial.
/// - Butterworth: `|S21|² = 1 / (1 + Ω^{2N})`.
///
/// The returned values are the (real, zero-phase) magnitude of `S21`; the
/// closed-form response models magnitude only.
pub fn ideal_response(proj: &FilterProject, freqs_hz: &[f64]) -> Vec<Complex64> {
    let n = proj.prototype.order();
    let omega0 = proj.spec.f0_hz;
    let fbw = proj.spec.fbw;
    freqs_hz
        .iter()
        .map(|&f| {
            let s21_sq = if f <= 0.0 {
                0.0
            } else {
                let omega = lowpass_to_bandpass(f, omega0, fbw);
                lowpass_s21_squared(proj.spec.approximation, n, omega)
            };
            Complex64::new(s21_sq.sqrt(), 0.0)
        })
        .collect()
}

/// Evaluate the ideal forward transmission `S21` of a **low-pass** filter over
/// `freqs_hz` (Filter Phase App.2.2, ADR-0139).
///
/// The low-pass analogue of [`ideal_response`]: it reuses the *same* closed-form
/// magnitude response [`lowpass_s21_squared`] but evaluates it at the bare
/// low-pass frequency variable `Ω = f / f_c` — there is **no** band-pass
/// frequency transform (a low-pass filter is already expressed in the prototype
/// `Ω` domain, scaled only by the cutoff `f_c`):
///
/// - Butterworth: `|S21|² = 1 / (1 + Ω^{2N})`.
/// - Chebyshev: `|S21|² = 1 / (1 + ε²·T_N²(Ω))`, `ε = √(10^{L_Ar/10} − 1)`.
///
/// At `Ω = 1` (i.e. `f = f_c`) Butterworth gives the defining `−3.01 dB`
/// half-power edge and Chebyshev the `−ripple_db` equi-ripple edge; past `f_c`
/// the response rolls off toward the `−20·N·log10(f/f_c)` asymptote. Frequencies
/// `f ≤ 0` map to a fully-rejected `0`. The returned values are the (real,
/// zero-phase) magnitude of `S21`; the closed-form model is magnitude only.
///
/// Unlike [`ideal_response`] this takes the [`Approximation`] / `order` / cutoff
/// directly rather than a [`FilterProject`], because a low-pass design has no
/// band-pass coupling matrix — the magnitude response is fully determined by
/// `(approx, order, f_c)`.
pub fn ideal_response_lowpass(
    approx: Approximation,
    order: usize,
    cutoff_hz: f64,
    freqs_hz: &[f64],
) -> Vec<Complex64> {
    freqs_hz
        .iter()
        .map(|&f| {
            let s21_sq = if f <= 0.0 {
                0.0
            } else {
                let omega = f / cutoff_hz;
                lowpass_s21_squared(approx, order, omega)
            };
            Complex64::new(s21_sq.sqrt(), 0.0)
        })
        .collect()
}

/// Closed-form lowpass `|S21|²(Ω)` for the given approximation and order.
fn lowpass_s21_squared(approx: Approximation, n: usize, omega: f64) -> f64 {
    match approx {
        Approximation::Butterworth => 1.0 / (1.0 + omega.powi(2 * n as i32)),
        Approximation::Chebyshev { ripple_db } => {
            let eps_sq = 10f64.powf(ripple_db / 10.0) - 1.0;
            let t = chebyshev_t(n, omega);
            1.0 / (1.0 + eps_sq * t * t)
        }
    }
}

/// Chebyshev polynomial of the first kind, `T_N(x)`.
///
/// `T_N(x) = cos(N·acos(x))` for `|x| ≤ 1` and `cosh(N·acosh(|x|))` (with the
/// correct sign for odd `N`) for `|x| > 1`.
fn chebyshev_t(n: usize, x: f64) -> f64 {
    if x.abs() <= 1.0 {
        ((n as f64) * x.acos()).cos()
    } else {
        // For |x| > 1: T_N(x) = sign·cosh(N·acosh(|x|)); sign = sgn(x)^N.
        let mag = ((n as f64) * x.abs().acosh()).cosh();
        if x < 0.0 && n % 2 == 1 { -mag } else { mag }
    }
}

/// Per-point detail of a mask check.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaskPoint {
    /// Frequency, Hz.
    pub freq_hz: f64,
    /// Insertion loss `−20·log10(|S21|)`, dB (positive = loss).
    pub insertion_loss_db: f64,
    /// Return loss `−20·log10(|S11|)`, dB (positive = better match).
    pub return_loss_db: f64,
}

/// Result of grading a synthesized response against its [`SpecMask`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaskReport {
    /// Overall verdict: `true` iff every mask requirement is met.
    pub pass: bool,
    /// Worst-case passband insertion-loss ripple observed, dB.
    pub worst_passband_ripple_db: f64,
    /// Worst-case (smallest) in-band return loss observed, dB.
    pub worst_return_loss_db: f64,
    /// Per stopband point: `(freq_hz, achieved_rejection_db, required_db, met)`.
    pub stopband: Vec<(f64, f64, f64, bool)>,
    /// Human-readable reasons for any failure (empty when `pass`).
    pub failures: Vec<String>,
}

/// Grade the synthesized ideal response against the spec mask (spec §3,
/// stage-3 gate).
///
/// Passes iff, over the swept `freqs_hz`:
/// - in-band insertion-loss ripple ≤ `mask.passband_ripple_db`,
/// - in-band return loss ≥ `mask.return_loss_db`, and
/// - at every stopband point, the rejection (insertion loss) ≥ the point's
///   required minimum.
///
/// "In-band" is `|Ω| ≤ 1` under the bandpass map (i.e. between the band edges
/// `f1`, `f2`). Stopband rejection is evaluated at each mask point's frequency
/// directly (interpolation is not required — the closed-form response is
/// evaluated exactly at that frequency).
pub fn check_mask(proj: &FilterProject, freqs_hz: &[f64]) -> MaskReport {
    let mask = &proj.spec.mask;
    let omega0 = proj.spec.f0_hz;
    let fbw = proj.spec.fbw;

    // ---- passband sweep ---------------------------------------------------
    let mut min_il = f64::INFINITY; // best (smallest) in-band insertion loss
    let mut max_il = f64::NEG_INFINITY; // worst (largest) in-band insertion loss
    let mut worst_rl = f64::INFINITY; // smallest in-band return loss
    let mut saw_passband = false;

    let s21 = ideal_response(proj, freqs_hz);
    for (&f, s21f) in freqs_hz.iter().zip(s21.iter()) {
        if f <= 0.0 {
            continue;
        }
        let omega = lowpass_to_bandpass(f, omega0, fbw);
        if omega.abs() > 1.0 {
            continue; // out of band; graded by the stopband points instead
        }
        saw_passband = true;
        let s21_mag = s21f.norm();
        let s11_sq = (1.0 - s21_mag * s21_mag).max(0.0);
        let il_db = -20.0 * s21_mag.max(1e-300).log10();
        let rl_db = if s11_sq <= 0.0 {
            f64::INFINITY
        } else {
            -10.0 * s11_sq.log10()
        };
        min_il = min_il.min(il_db);
        max_il = max_il.max(il_db);
        worst_rl = worst_rl.min(rl_db);
    }

    let ripple = if saw_passband {
        (max_il - min_il).max(0.0)
    } else {
        0.0
    };

    let mut failures = Vec::new();
    if saw_passband {
        if ripple > mask.passband_ripple_db + 1e-9 {
            failures.push(format!(
                "passband ripple {ripple:.3} dB exceeds spec {:.3} dB",
                mask.passband_ripple_db
            ));
        }
        if worst_rl + 1e-9 < mask.return_loss_db {
            failures.push(format!(
                "in-band return loss {worst_rl:.3} dB below spec {:.3} dB",
                mask.return_loss_db
            ));
        }
    } else {
        failures.push("no swept frequency fell inside the passband".to_string());
    }

    // ---- stopband points --------------------------------------------------
    let mut stopband = Vec::with_capacity(mask.stopband.len());
    for &(f_hz, required_db) in &mask.stopband {
        let s21_sq = if f_hz <= 0.0 {
            0.0
        } else {
            let omega = lowpass_to_bandpass(f_hz, omega0, fbw);
            lowpass_s21_squared(proj.spec.approximation, proj.prototype.order(), omega)
        };
        let s21_mag = s21_sq.sqrt();
        let rejection_db = -20.0 * s21_mag.max(1e-300).log10();
        let met = rejection_db + 1e-9 >= required_db;
        if !met {
            failures.push(format!(
                "stopband {f_hz:.3e} Hz rejection {rejection_db:.3} dB below required {required_db:.3} dB"
            ));
        }
        stopband.push((f_hz, rejection_db, required_db, met));
    }

    MaskReport {
        pass: failures.is_empty(),
        worst_passband_ripple_db: ripple,
        worst_return_loss_db: if worst_rl.is_finite() { worst_rl } else { 0.0 },
        stopband,
        failures,
    }
}
