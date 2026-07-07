//! Closed-form microstrip dimensional synthesis (Filter Phases F1.2.0 / F1.2.2 /
//! F1.2.3 / F1.2.5 / F1.2.7 / F1.2.8).
//!
//! Turns an abstract synthesized [`crate::CouplingMatrix`] (or low-pass
//! prototype) into **physical microstrip dimensions** for four coupled-resonator
//! band-pass topologies — the edge-coupled half-wave filter
//! ([`dimension_edge_coupled`], F1.2.0), the U-folded **hairpin** filter
//! ([`dimension_hairpin`], F1.2.2), the capacitively-loaded **combline**
//! filter ([`dimension_combline`], F1.2.5), and the **interdigital** filter
//! ([`dimension_interdigital`], F1.2.7; board [`dimension_interdigital_layout`],
//! F1.2.8) — plus the **stepped-impedance**
//! low-pass ([`dimension_stepped_impedance`], F1.2.3), by inverting the
//! already-validated `yee-layout` closed-form models. Pure `f64`, WASM-safe, NO
//! FDTD, NO surrogate — this is the *initial* dimensioning that seeds the later
//! EM-in-the-loop refinement (F1.2.1).
//!
//! The four coupled-resonator band-pass topologies share the **same
//! inter-resonator coupling mechanism**: adjacent resonators couple through the
//! edge gap between their lines — the edge-coupled gap→`k` inversion. They
//! therefore reuse the identical [`solve_gap`] bisection and the
//! `target_k = FBW · m_{i,i+1}` derivation below; only the resonator geometry
//! differs — a straight λ/2 strip (edge-coupled), a folded half-wave = two ≈λ/4
//! arms (hairpin, see [`dimension_hairpin`]), a short-circuited θ0 < π/2 line
//! capacitively loaded to resonance (combline, see [`dimension_combline`]), or a
//! full λ/4 line short-circuited at *alternating* ends with no cap (interdigital,
//! the θ = π/2 limit of combline, see [`dimension_interdigital`]).
//!
//! # Method (Hong & Lancaster ch. 8 / Pozar §8.7)
//!
//! For an `N`-pole edge-coupled half-wave filter:
//!
//! - **Line width** — the spec-`Z0` Hammerstad-Jensen synthesis width
//!   ([`yee_layout::microstrip_width`]).
//! - **Resonator length** — a half guided wavelength at `f0`,
//!   `ℓ = λ_g/2 = c / (2·f0·√ε_eff)`, with `ε_eff` from
//!   [`yee_layout::eps_eff`] at the synthesized width (`c = 299_792_458` m/s).
//! - **Inter-resonator gaps** — for each adjacent resonator pair `(i, i+1)` the
//!   coupling coefficient `k_{i,i+1} = FBW · m_{i,i+1}` is realized by a coupled
//!   section whose voltage coupling
//!   `(Z0e − Z0o)/(Z0e + Z0o)` ([`yee_layout::coupling_coefficient`]) equals
//!   `k_{i,i+1}`. Because that coupling is **strictly decreasing in the gap `s`**
//!   (`yee-layout`'s `coupled_002` gate), the inverse "gap that realizes a target
//!   `k`" is found exactly by **bisection** — no optimizer, no FDTD.
//!
//! ## Cross-check: `target_k = FBW · m_{i,i+1}` equals `yee-synth`'s `k`
//!
//! `yee-synth::coupling_design` builds the normalized matrix with
//! `m[i][i+1] = 1/√(g_i g_{i+1})` and the inter-resonator coupling
//! `k_{i,i+1} = FBW / √(g_i g_{i+1})`. Hence
//! `FBW · m[i][i+1] = FBW / √(g_i g_{i+1}) = k_{i,i+1}` exactly — so multiplying
//! the off-diagonal of [`crate::CouplingMatrix::m`] by `spec.fbw` reproduces the
//! synthesized `k` vector, which is the target each gap is solved for.

use serde::{Deserialize, Serialize};

use yee_layout::{
    BBox, EdgeCoupledParams, EdgeCoupledSection, HairpinSectionParams, Layout, Point2, Polygon,
    PortRef, Substrate, coupled_microstrip, coupling_coefficient, edge_coupled_bpf, eps_eff,
    hairpin_bpf_sections, microstrip_width,
};

use crate::{FilterProject, Topology};

/// Speed of light in vacuum, m/s (exact, SI definition).
const C: f64 = 299_792_458.0;

/// Gap-bisection bracket lower bound, metres (5 µm — tightest realizable etch
/// gap; tighter gaps over-couple and are unmanufacturable).
const GAP_MIN_M: f64 = 5.0e-6;
/// Gap-bisection bracket upper bound, metres (5 mm — beyond this the strips are
/// effectively uncoupled and `k → 0`).
const GAP_MAX_M: f64 = 5.0e-3;
/// Relative tolerance on the realized coupling for the gap bisection.
const GAP_REL_TOL: f64 = 1.0e-4;
/// Hard cap on bisection iterations (≈ log2((5e-3 − 5e-6)/(5e-6·1e-4)) ≈ 33, so
/// 200 is comfortably above the worst case and guards against non-convergence).
const GAP_MAX_ITERS: usize = 200;

/// First-order physical dimensions of an edge-coupled half-wave microstrip
/// band-pass filter, synthesized from a [`crate::CouplingMatrix`].
///
/// All lengths are in metres. `gaps_m` and `target_k` are both length `N − 1`
/// (one per adjacent resonator pair) and index-aligned: `gaps_m[i]` is the gap
/// that realizes `target_k[i]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EdgeCoupledDimensions {
    /// Resonator / feed line width for the spec `Z0`, metres (Hammerstad-Jensen).
    pub line_width_m: f64,
    /// Resonator length `≈ λ_g/2` at `f0`, metres (via `ε_eff`).
    pub resonator_length_m: f64,
    /// Inter-resonator coupled-section gaps, metres (length `N − 1`).
    pub gaps_m: Vec<f64>,
    /// The `FBW · m_{i,i+1}` coupling each gap was solved for (length `N − 1`).
    pub target_k: Vec<f64>,
}

/// Errors from [`dimension_edge_coupled`] / [`dimension_edge_coupled_layout`].
#[derive(Debug, Clone, PartialEq)]
pub enum DimError {
    /// The project's [`Topology`] is not [`Topology::CoupledResonator`], which is
    /// the only topology edge-coupled dimensioning supports.
    UnsupportedTopology,
    /// The filter order `N < 2`: there is no inter-resonator coupling to realize.
    OrderTooSmall,
    /// A `target_k` could not be realized by any gap in the bisection bracket.
    /// Carries the resonator-pair index, the unreachable target, and the
    /// achievable coupling range `(k_at_max_gap, k_at_min_gap)` over the bracket
    /// (`k` decreases with the gap, so `k_at_min_gap` is the largest realizable).
    GapNotBracketed {
        /// Adjacent-resonator-pair index `i` (the `i`-th of `N − 1` gaps).
        index: usize,
        /// The `FBW · m_{i,i+1}` target that fell outside the achievable range.
        target_k: f64,
        /// Smallest realizable coupling (at the maximum bracket gap).
        k_min: f64,
        /// Largest realizable coupling (at the minimum bracket gap).
        k_max: f64,
    },
    /// The synthesized external Q cannot be realized by a tapped feed on this
    /// resonator (R.4/F1.2.1): a tap can only realize
    /// `qe ∈ [(π/2)(Z0/Zr), (π/2)(Z0/Zr)/cos²(π·arm/L)]` — the lower bound is
    /// the voltage antinode (open end), the upper bound the last physical
    /// point on the arm before the bend. Below it the coupling is stronger
    /// than any tap can deliver; above it the tap lands on the bend.
    TapNotRealizable {
        /// The synthesized singly-loaded external Q (`g0·g1/FBW`).
        qe: f64,
        /// Smallest realizable external Q (tap at the antinode).
        qe_min: f64,
        /// Largest realizable external Q (tap at the end of the arm);
        /// `+∞` when no arm bound applies (bare [`tap_offset_from_qe`]).
        qe_max: f64,
    },
    /// A stepped-impedance input was non-physical: the prototype order is `0`,
    /// or the cut-off frequency / impedances (`f_c`, `Z₀`, `Z_high`, `Z_low`) are
    /// not strictly positive. Carries a human-readable description.
    NonPhysicalInput(&'static str),
    /// A combline resonator electrical length `θ0` was not in the open interval
    /// `(0, π/2)`. The combline loading capacitor is `C_L = cot(θ0)/(2π·f0·Z0)`,
    /// which is only positive (physical) for `θ0 ∈ (0, π/2)`; at `θ0 = π/2` the
    /// line is already self-resonant (`cot = 0` → `C_L = 0`) and beyond it
    /// `cot < 0` would demand a non-physical negative capacitance. Carries the
    /// offending `θ0` in radians.
    InvalidTheta0(f64),
}

impl std::fmt::Display for DimError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DimError::UnsupportedTopology => write!(
                f,
                "dimensional synthesis supports only Topology::CoupledResonator"
            ),
            DimError::OrderTooSmall => write!(
                f,
                "filter order N must be >= 2 for inter-resonator coupling to realize"
            ),
            DimError::GapNotBracketed {
                index,
                target_k,
                k_min,
                k_max,
            } => write!(
                f,
                "target_k[{index}] = {target_k:.6} is unreachable in the gap bracket \
                 [{GAP_MIN_M:.1e}, {GAP_MAX_M:.1e}] m; achievable coupling range is \
                 [{k_min:.6}, {k_max:.6}]"
            ),
            DimError::TapNotRealizable { qe, qe_min, qe_max } => write!(
                f,
                "external Q qe = {qe:.4} is outside the tapped-feed realizable range \
                 [{qe_min:.4}, {qe_max:.4}] (antinode tap .. end-of-arm tap)"
            ),
            DimError::NonPhysicalInput(why) => {
                write!(f, "non-physical stepped-impedance input: {why}")
            }
            DimError::InvalidTheta0(theta0) => write!(
                f,
                "combline resonator electrical length theta0 = {theta0:.6} rad must be in \
                 (0, pi/2); cot(theta0) <= 0 outside it gives a non-physical loading cap"
            ),
        }
    }
}

impl std::error::Error for DimError {}

/// Tap offset (metres, from the **open end** of the resonator) realizing a
/// singly-loaded external Q on a tapped half-wave resonator (R.4/F1.2.1;
/// Hong & Lancaster ch. 6 tapped-line coupling).
///
/// For a half-wave resonator of characteristic impedance `zr_ohm` and total
/// electrical length `L = halfwave_m` (an unfolded hairpin: `2 * arm`),
/// tapped by a `z0_ohm` feed at distance `t` from the open end, the standing
/// wave `V(x) = V0 cos(pi x / L)` (antinode at the open ends, null at the
/// fold) stores `W = C_l V0^2 L / 4` and delivers `P = V(t)^2 / (2 Z0)` to
/// the feed, so
///
/// `Qe(t) = w0 W / P = (pi/2) (Z0/Zr) / cos^2(pi t / L)`
///
/// which inverts to `t = (L/pi) acos( sqrt( (pi/2)(Z0/Zr) / qe ) )`. `Qe` is
/// smallest (strongest coupling) at the antinode `t = 0` and diverges toward
/// the fold.
///
/// # Errors
///
/// [`DimError::TapNotRealizable`] when `qe < (pi/2)(Z0/Zr)` — the requested
/// coupling is stronger than a tap at the antinode can deliver.
pub fn tap_offset_from_qe(
    qe: f64,
    z0_ohm: f64,
    zr_ohm: f64,
    halfwave_m: f64,
) -> Result<f64, DimError> {
    let qe_min = std::f64::consts::FRAC_PI_2 * z0_ohm / zr_ohm;
    if qe < qe_min {
        return Err(DimError::TapNotRealizable {
            qe,
            qe_min,
            qe_max: f64::INFINITY,
        });
    }
    Ok(halfwave_m / std::f64::consts::PI * (qe_min / qe).sqrt().acos())
}

/// Solve for the gap `s` (metres) whose edge-coupled-line coupling coefficient
/// equals `target_k`, by bisection over `[GAP_MIN_M, GAP_MAX_M]`.
///
/// `coupling_coefficient(coupled_microstrip(w, s, h, εr))` is strictly
/// decreasing in `s`, so the coupling is largest at `GAP_MIN_M` and smallest at
/// `GAP_MAX_M`. If `target_k` falls outside `[k(GAP_MAX_M), k(GAP_MIN_M)]` the
/// target is not bracketed and [`DimError::GapNotBracketed`] is returned (no
/// silent clamping). Converges to a relative tolerance of [`GAP_REL_TOL`] on the
/// realized coupling, capped at [`GAP_MAX_ITERS`] iterations.
fn solve_gap(index: usize, target_k: f64, w_m: f64, h_m: f64, eps_r: f64) -> Result<f64, DimError> {
    let k_of = |s: f64| coupling_coefficient(&coupled_microstrip(w_m, s, h_m, eps_r));

    // k decreases with gap: k_max at the smallest gap, k_min at the largest.
    let k_max = k_of(GAP_MIN_M);
    let k_min = k_of(GAP_MAX_M);
    if !(k_min..=k_max).contains(&target_k) {
        return Err(DimError::GapNotBracketed {
            index,
            target_k,
            k_min,
            k_max,
        });
    }

    // Bisect: `lo` is the small-gap (high-k) end, `hi` the large-gap (low-k) end.
    let mut lo = GAP_MIN_M;
    let mut hi = GAP_MAX_M;
    let mut mid = 0.5 * (lo + hi);
    for _ in 0..GAP_MAX_ITERS {
        mid = 0.5 * (lo + hi);
        let k_mid = k_of(mid);
        // Relative-tolerance convergence on the realized coupling.
        if (k_mid - target_k).abs() <= GAP_REL_TOL * target_k.abs().max(f64::MIN_POSITIVE) {
            return Ok(mid);
        }
        if k_mid > target_k {
            // Coupling too strong → widen the gap.
            lo = mid;
        } else {
            // Coupling too weak → narrow the gap.
            hi = mid;
        }
    }
    // Loop exhausted without hitting GAP_REL_TOL. For a strictly-monotone
    // `coupling_coefficient` over [GAP_MIN_M, GAP_MAX_M] (proven by `coupled_002`)
    // bisection converges far inside GAP_MAX_ITERS, so reaching here signals the
    // monotonicity assumption broke (e.g. a future coupled-microstrip model
    // change). Trip it in debug/test builds; release returns the best estimate.
    debug_assert!(
        (k_of(mid) - target_k).abs() <= 100.0 * GAP_REL_TOL * target_k.abs().max(f64::MIN_POSITIVE),
        "solve_gap: bisection did not converge in {GAP_MAX_ITERS} iters (index {index}, \
         target_k {target_k}) — coupling_coefficient may be non-monotone over the gap bracket"
    );
    Ok(mid)
}

/// Invert the validated coupled-microstrip model to size an edge-coupled
/// half-wave BPF from a synthesized [`FilterProject`] and a [`Substrate`].
///
/// Closed-form throughout: the line width is the spec-`Z0` Hammerstad-Jensen
/// width, the resonator length is `λ_g/2` at `f0` (via `ε_eff`), and each
/// inter-resonator gap is found by bisecting the (monotonic) coupled-line
/// coupling coefficient onto `FBW · m_{i,i+1}`. See the [module docs](self) for
/// the method and the `target_k = FBW · m` cross-check.
///
/// # Errors
///
/// - [`DimError::UnsupportedTopology`] if the project is not
///   [`Topology::CoupledResonator`].
/// - [`DimError::OrderTooSmall`] if the order `N < 2` (no inter-resonator
///   coupling to realize).
/// - [`DimError::GapNotBracketed`] if a `target_k` is unreachable for any gap in
///   the `[5 µm, 5 mm]` bracket at the synthesized width (no silent clamping).
pub fn dimension_edge_coupled(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<EdgeCoupledDimensions, DimError> {
    if project.topology != Topology::CoupledResonator {
        return Err(DimError::UnsupportedTopology);
    }

    let n = project.coupling.m.len();
    if n < 2 {
        return Err(DimError::OrderTooSmall);
    }

    let eps_r = substrate.eps_r;
    let h_m = substrate.height_m;
    let f0 = project.spec.f0_hz;
    let fbw = project.spec.fbw;
    let z0 = project.spec.z0_ohm;

    // 1. Line width from the Hammerstad-Jensen Z0 synthesis.
    let line_width_m = microstrip_width(z0, eps_r, h_m);

    // 2. Resonator length = λ_g/2 = c / (2·f0·√ε_eff).
    let e_eff = eps_eff(line_width_m, h_m, eps_r);
    let resonator_length_m = C / (2.0 * f0 * e_eff.sqrt());

    // 3. Inter-resonator gaps: target_k[i] = FBW · m[i][i+1] (= yee-synth's
    //    k_{i,i+1} = FBW/√(g_i g_{i+1}); see module docs), solved by bisection.
    let mut target_k = Vec::with_capacity(n - 1);
    let mut gaps_m = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let k_i = fbw * project.coupling.m[i][i + 1];
        let gap = solve_gap(i, k_i, line_width_m, h_m, eps_r)?;
        target_k.push(k_i);
        gaps_m.push(gap);
    }

    Ok(EdgeCoupledDimensions {
        line_width_m,
        resonator_length_m,
        gaps_m,
        target_k,
    })
}

/// Convenience: assemble a [`yee_layout::Layout`] from the synthesized
/// dimensions via the existing [`yee_layout::edge_coupled_bpf`].
///
/// Builds the `N` coupled half-wave sections (all of width `line_width_m` and
/// length `resonator_length_m`) with the `N − 1` solved inter-resonator gaps.
/// `edge_coupled_bpf` reads each section's `gap_m` as the gap *to the next*
/// section, so the last section has no real successor; its `gap_m` is set to a
/// **documented placeholder** — the first inter-resonator gap — purely so the
/// struct is well-formed. The feed-line width is `line_width_m` and the feed
/// length is one resonator length (a neutral default). Mapping the external Q
/// (`qe_in`/`qe_out`) to a feed/tap geometry is **deferred to F1.2.1**; this
/// function does **not** invent a `qe`→gap formula.
///
/// # Errors
///
/// Propagates every [`DimError`] from [`dimension_edge_coupled`].
pub fn dimension_edge_coupled_layout(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<Layout, DimError> {
    let dims = dimension_edge_coupled(project, substrate)?;

    let n = dims.gaps_m.len() + 1; // N resonators, N−1 gaps.
    // Placeholder gap for the trailing section (no successor strip); documented
    // above — qe→feed dimensioning is F1.2.1.
    let placeholder_gap = dims.gaps_m[0];

    let sections: Vec<EdgeCoupledSection> = (0..n)
        .map(|i| EdgeCoupledSection {
            length_m: dims.resonator_length_m,
            width_m: dims.line_width_m,
            gap_m: if i < n - 1 {
                dims.gaps_m[i]
            } else {
                placeholder_gap
            },
        })
        .collect();

    let params = EdgeCoupledParams {
        substrate: *substrate,
        sections,
        feed_width_m: dims.line_width_m,
        feed_length_m: dims.resonator_length_m,
    };

    Ok(edge_coupled_bpf(&params))
}

// ---------------------------------------------------------------------------
// Hairpin (U-folded half-wave) dimensional synthesis — Filter Phase F1.2.2.
// ---------------------------------------------------------------------------

/// Multiple of the line width used for the intra-hairpin arm spacing.
///
/// `fold_spacing_m` is the centre-to-centre distance between the *two arms of
/// one* hairpin — a weak self-coupling internal to a single resonator, **not**
/// the inter-resonator coupling that sets the filter response (that is the edge
/// gap, solved by [`solve_gap`]). A fixed sensible value is therefore adequate
/// for the walking skeleton; two line widths keeps the arms close enough to fold
/// a compact U without the arms shorting, and is the conventional Hong &
/// Lancaster ch. 6 starting choice. F1.2.1 BO refines it against EM.
const HAIRPIN_FOLD_SPACING_WIDTHS: f64 = 2.0;

/// First-order physical dimensions of a **hairpin** (U-folded half-wave)
/// microstrip band-pass filter, synthesized from a [`crate::CouplingMatrix`].
///
/// All lengths are in metres. `gaps_m` and `target_k` are both length `N − 1`
/// (one per adjacent resonator pair) and index-aligned: `gaps_m[i]` is the edge
/// gap between resonators `i` and `i + 1` that realizes `target_k[i]`.
///
/// Mirrors [`EdgeCoupledDimensions`]; the difference is the resonator geometry —
/// `arm_length_m` is `λ_g/4` (a folded half-wave is two ≈λ/4 arms) rather than
/// the edge-coupled `λ_g/2` straight length, plus `fold_spacing_m` for the
/// intra-hairpin arm pitch.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HairpinDimensions {
    /// Resonator-arm / feed line width for the spec `Z0`, metres
    /// (Hammerstad-Jensen, via [`yee_layout::microstrip_width`]).
    pub line_width_m: f64,
    /// Length of each resonator arm `≈ λ_g/4` at `f0`, metres (via `ε_eff`). The
    /// U-folded half-wave resonator is two of these arms joined by a bend.
    pub arm_length_m: f64,
    /// Centre-to-centre spacing of the two arms of one hairpin, metres (a fixed
    /// closed-form choice — [`HAIRPIN_FOLD_SPACING_WIDTHS`] line widths — since it
    /// is intra-resonator self-coupling, not the inter-resonator coupling).
    pub fold_spacing_m: f64,
    /// Inter-resonator edge-coupling gaps, metres (length `N − 1`).
    pub gaps_m: Vec<f64>,
    /// The `FBW · m_{i,i+1}` coupling each gap was solved for (length `N − 1`).
    pub target_k: Vec<f64>,
    /// Tapped-feed offset up the end-resonator arms realizing the synthesized
    /// external Q (R.4/F1.2.1; [`tap_offset_from_qe`] with `L = 2·arm`),
    /// metres from the open end.
    pub tap_offset_m: f64,
}

/// Invert the validated coupled-microstrip model to size a **hairpin**
/// (U-folded half-wave) BPF from a synthesized [`FilterProject`] and a
/// [`Substrate`].
///
/// Closed-form throughout and a direct mirror of [`dimension_edge_coupled`]; the
/// only physical difference is the resonator geometry:
///
/// - **Line width** — the spec-`Z0` Hammerstad-Jensen width
///   ([`yee_layout::microstrip_width`]).
/// - **Arm length** — `arm_length_m = λ_g/4 = c / (4·f0·√ε_eff)`. A hairpin
///   resonator is a half-wave (`λ_g/2`) line *folded into a U*, i.e. two arms of
///   `≈ λ_g/4` joined by a bend — hence the **factor-4** here versus the
///   edge-coupled straight half-wave's **factor-2** (`λ_g/2`). `ε_eff` is
///   evaluated at the synthesized width via [`yee_layout::eps_eff`]. (Hong &
///   Lancaster, *Microstrip Filters for RF/Microwave Applications*, ch. 6.)
/// - **Fold spacing** — a fixed closed-form choice
///   ([`HAIRPIN_FOLD_SPACING_WIDTHS`] line widths); the two arms of one hairpin
///   are weakly self-coupled (intra-resonator), *not* the inter-resonator
///   coupling, so a sensible constant suffices for the walking skeleton.
/// - **Inter-resonator gaps** — identical to edge-coupled: for each adjacent
///   pair `(i, i+1)`, `target_k[i] = FBW · m_{i,i+1}` is realized by bisecting the
///   monotone coupled-line coupling coefficient with the shared [`solve_gap`]
///   helper (no optimizer, no FDTD). See the [module docs](self) for the
///   `target_k = FBW · m` cross-check.
///
/// # Errors
///
/// - [`DimError::UnsupportedTopology`] if the project is not
///   [`Topology::CoupledResonator`] (the only synthesized topology today; the
///   hairpin is a *realization* of that coupling network).
/// - [`DimError::OrderTooSmall`] if the order `N < 2` (no inter-resonator
///   coupling to realize).
/// - [`DimError::GapNotBracketed`] if a `target_k` is unreachable for any gap in
///   the `[5 µm, 5 mm]` bracket at the synthesized width (no silent clamping).
pub fn dimension_hairpin(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<HairpinDimensions, DimError> {
    dimension_hairpin_with_fold(project, substrate, HAIRPIN_FOLD_SPACING_WIDTHS)
}

/// [`dimension_hairpin`] with an explicit fold spacing in **line widths**
/// (R.4a). The fold pitch trades layout compactness against intra-resonator
/// self-coupling; it also consumes resonator length (see the arm formula), so
/// a tighter fold leaves more arm for the tap. The default wrapper uses
/// [`HAIRPIN_FOLD_SPACING_WIDTHS`].
///
/// **Fold-corrected arm length (R.4a).** The hairpin is a half-wave resonator
/// *folded into a U*: its total midline length — up one arm, across the bend,
/// down the other arm — must be `λ_g/2`, so
///
/// `arm = (λ_g/2 − fold_spacing) / 2`
///
/// with the bend's midline run ≈ the centre-to-centre `fold_spacing`. The
/// F1.2.2 walking skeleton used `arm = λ_g/4` (fold length *extra*), which the
/// first full-wave BPF run (gate `engine-bpf-verify-001`) measured as a
/// low-shifted, wrecked passband — every resonator was electrically long by
/// the full bend path (~37 % here). First-order correction; the corner
/// mitres are R.4b's BO residual.
///
/// # Errors
///
/// All of [`dimension_hairpin`]'s errors, plus
/// [`DimError::NonPhysicalInput`] when the fold consumes the whole half-wave
/// (`arm ≤ 0`) and [`DimError::TapNotRealizable`] when the qe→tap offset
/// falls beyond the (fold-corrected) arm.
pub fn dimension_hairpin_with_fold(
    project: &FilterProject,
    substrate: &Substrate,
    fold_widths: f64,
) -> Result<HairpinDimensions, DimError> {
    if project.topology != Topology::CoupledResonator {
        return Err(DimError::UnsupportedTopology);
    }

    let n = project.coupling.m.len();
    if n < 2 {
        return Err(DimError::OrderTooSmall);
    }

    let eps_r = substrate.eps_r;
    let h_m = substrate.height_m;
    let f0 = project.spec.f0_hz;
    let fbw = project.spec.fbw;
    let z0 = project.spec.z0_ohm;

    // 1. Line width from the Hammerstad-Jensen Z0 synthesis.
    let line_width_m = microstrip_width(z0, eps_r, h_m);

    // 2. Fold spacing: intra-hairpin arm pitch (centre-to-centre). At
    //    fold_widths <= 1 the two arms of the U touch and the resonator
    //    degenerates into a solid block (measured as a wrecked response by
    //    an early engine-bpf-verify-001 iteration) — reject it.
    if fold_widths <= 1.0 {
        return Err(DimError::NonPhysicalInput(
            "hairpin fold spacing is centre-to-centre: fold_widths <= 1 merges the two arms",
        ));
    }
    let fold_spacing_m = fold_widths * line_width_m;

    // 3. Fold-corrected arm length: the U's midline — arm + fold + arm —
    //    is the half-wave, so arm = (lambda_g/2 - fold)/2 (doc above).
    let e_eff = eps_eff(line_width_m, h_m, eps_r);
    let halfwave_m = C / (2.0 * f0 * e_eff.sqrt());
    let arm_length_m = (halfwave_m - fold_spacing_m) / 2.0;
    if arm_length_m <= 0.0 {
        return Err(DimError::NonPhysicalInput(
            "hairpin fold spacing consumes the entire half-wave resonator",
        ));
    }

    // 4. Inter-resonator gaps: target_k[i] = FBW * m[i][i+1] (= yee-synth's
    //    k_{i,i+1}; see module docs), solved by the SAME bisection as
    //    edge-coupled because adjacent hairpins couple through the edge gap
    //    between their arms.
    let mut target_k = Vec::with_capacity(n - 1);
    let mut gaps_m = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let k_i = fbw * project.coupling.m[i][i + 1];
        let gap = solve_gap(i, k_i, line_width_m, h_m, eps_r)?;
        target_k.push(k_i);
        gaps_m.push(gap);
    }

    // 5. Tapped-feed offset from the synthesized external Q (R.4/F1.2.1).
    //    The hairpin resonator line is the spec-Z0 width, so Zr = Z0; qe_in
    //    and qe_out are equal for the symmetric prototypes synthesis emits.
    //    t runs along the resonator from the open end, so it must land on
    //    the physical arm (the region [arm, halfwave/2] is the bend).
    let tap_offset_m = tap_offset_from_qe(project.coupling.qe_in, z0, z0, halfwave_m)?;
    if tap_offset_m > arm_length_m {
        let qe_min = std::f64::consts::FRAC_PI_2;
        let qe_max = qe_min
            / (std::f64::consts::PI * arm_length_m / halfwave_m)
                .cos()
                .powi(2);
        return Err(DimError::TapNotRealizable {
            qe: project.coupling.qe_in,
            qe_min,
            qe_max,
        });
    }

    Ok(HairpinDimensions {
        line_width_m,
        arm_length_m,
        fold_spacing_m,
        gaps_m,
        target_k,
        tap_offset_m,
    })
}

/// Assemble a [`yee_layout::Layout`] from the synthesized hairpin dimensions
/// via the **per-section** generator [`yee_layout::hairpin_bpf_sections`]
/// (R.4/F1.2.1 — this closed the two gaps the F1.2.2 walking skeleton
/// documented):
///
/// - **Per-section gaps (gap option (a)).** Each adjacent resonator pair sits
///   at its own solved [`HairpinDimensions::gaps_m`] — the mean-gap collapse
///   through the single-`coupling_gap_m` [`HairpinParams`] is gone. The
///   uniform-gap [`yee_layout::hairpin_bpf`] (and its committed `geo-003`
///   gate) is untouched.
/// - **qe→tap.** The feed taps the end-resonator arms at
///   [`HairpinDimensions::tap_offset_m`], the [`tap_offset_from_qe`] solution
///   for the synthesized external Q — the `arm_length/3` placeholder is gone.
///
/// `feed_width_m = line_width_m` and `feed_length_m = arm_length_m` remain
/// simple closed-form choices (the feed is a plain spec-`Z0` line).
///
/// # Errors
///
/// Propagates every [`DimError`] from [`dimension_hairpin`] (including
/// [`DimError::TapNotRealizable`] from the qe→tap solve).
pub fn dimension_hairpin_layout(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<Layout, DimError> {
    let dims = dimension_hairpin(project, substrate)?;
    let params = HairpinSectionParams {
        substrate: *substrate,
        arm_length_m: dims.arm_length_m,
        line_width_m: dims.line_width_m,
        fold_spacing_m: dims.fold_spacing_m,
        gaps_m: dims.gaps_m.clone(),
        tap_offset_m: dims.tap_offset_m,
        feed_width_m: dims.line_width_m,
        feed_length_m: dims.arm_length_m,
    };
    Ok(hairpin_bpf_sections(&params))
}

// ---------------------------------------------------------------------------
// Stepped-impedance low-pass (alternating high-Z / low-Z lines) — F1.2.3.
// ---------------------------------------------------------------------------

/// One transmission-line section of a stepped-impedance low-pass filter.
///
/// Each low-pass-prototype reactive element `g_k` (k = 1..N) becomes one short
/// microstrip section, alternating shunt-capacitor (low-Z) / series-inductor
/// (high-Z) **starting with a shunt capacitor (low-Z)**. All lengths in metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SteppedSection {
    /// `true` for a series-inductor **high-Z** line (realizes an inductor);
    /// `false` for a shunt-capacitor **low-Z** line (realizes a capacitor).
    pub high_z: bool,
    /// The section's characteristic impedance, ohms — `Z_high` when
    /// [`high_z`](Self::high_z) is `true`, else `Z_low`.
    pub z_ohm: f64,
    /// Electrical length `βl` of the section, radians (Pozar §8.6:
    /// `g_k·Z_low/Z₀` for a low-Z line, `g_k·Z₀/Z_high` for a high-Z line).
    pub electrical_length_rad: f64,
    /// Physical microstrip width for `z_ohm` (Hammerstad-Jensen), metres.
    pub width_m: f64,
    /// Physical section length `l = (βl / 2π)·λ_g` at the section width, metres.
    pub length_m: f64,
}

/// First-order physical dimensions of a stepped-impedance low-pass microstrip
/// filter, synthesized from a low-pass [`yee_synth::Prototype`].
///
/// The `sections` are in physical order, **source → load**, one per reactive
/// prototype element `g_k` (k = 1..N), alternating low-Z / high-Z starting with
/// a low-Z (shunt-capacitor) line. Mirrors [`EdgeCoupledDimensions`] /
/// [`HairpinDimensions`] in shape.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SteppedImpedanceDimensions {
    /// The line sections in order, source → load (length `N`).
    pub sections: Vec<SteppedSection>,
    /// Substrate relative permittivity `ε_r` (carried for the layout step).
    pub eps_r: f64,
    /// Substrate height `h`, metres (carried for the layout step).
    pub h_m: f64,
}

/// Synthesize the alternating line sections of a **stepped-impedance low-pass
/// filter** from a low-pass prototype and the line-impedance choices, on a
/// [`Substrate`] (Filter Phase F1.2.3, Pozar §8.6).
///
/// Closed-form throughout, mirroring [`dimension_edge_coupled`]. For each
/// reactive prototype element `g_k` (k = 1..N), one short microstrip line
/// section is produced, alternating shunt-capacitor (low-Z) / series-inductor
/// (high-Z) **starting with a shunt capacitor (low-Z)** — the standard low-pass
/// prototype begins with a shunt element, so `sections[0].high_z == false`:
///
/// - **Shunt capacitor → low-Z line** (`z_low`): electrical length
///   `βl = g_k · Z_low / Z₀`.
/// - **Series inductor → high-Z line** (`z_high`): electrical length
///   `βl = g_k · Z₀ / Z_high`.
///
/// (Derivation: a high-Z line of electrical length `βl` looks inductive with
/// `L = (Z_high/ω)·βl`; matching the prototype inductance `L = g_k·Z₀/ω_c` at
/// `ω = ω_c` gives `βl = g_k·Z₀/Z_high`. Dually for the capacitive low-Z line.)
///
/// Physical dimensions per section: the width is the Hammerstad-Jensen
/// synthesis width for that section's impedance ([`yee_layout::microstrip_width`]);
/// the guided wavelength is `λ_g = c / (f_c · √ε_eff)` with `ε_eff` from
/// [`yee_layout::eps_eff`] at that section's width; the physical length is
/// `l = (βl / 2π) · λ_g`.
///
/// # Errors
///
/// Returns [`DimError::NonPhysicalInput`] if the prototype order is `0` or if
/// `f_c_hz`, `z0`, `z_high`, or `z_low` is not strictly positive.
pub fn dimension_stepped_impedance(
    proto: &yee_synth::Prototype,
    f_c_hz: f64,
    z0: f64,
    z_high: f64,
    z_low: f64,
    sub: &Substrate,
) -> Result<SteppedImpedanceDimensions, DimError> {
    let n = proto.order();
    if n == 0 {
        return Err(DimError::NonPhysicalInput("prototype order N must be >= 1"));
    }
    if f_c_hz <= 0.0 {
        return Err(DimError::NonPhysicalInput("f_c must be > 0"));
    }
    if z0 <= 0.0 || z_high <= 0.0 || z_low <= 0.0 {
        return Err(DimError::NonPhysicalInput(
            "Z0, Z_high and Z_low must all be > 0",
        ));
    }

    let eps_r = sub.eps_r;
    let h_m = sub.height_m;
    // `proto.g` is `[g0, g1, …, gN, g_{N+1}]`; `g[1..=N]` are the reactive
    // elements. Iterate those by `enumerate()` so the 1-based prototype index
    // `k = idx + 1` drives the low-Z-first alternation.
    let reactive = &proto.g[1..=n];

    let mut sections = Vec::with_capacity(n);
    for (idx, &g_k) in reactive.iter().enumerate() {
        let k = idx + 1; // 1-based prototype element index.
        // Section 1 (k = 1) is the shunt capacitor → low-Z; alternate from there.
        let high_z = k % 2 == 0;

        // Pozar §8.6 electrical length βl (radians).
        let (z_ohm, electrical_length_rad) = if high_z {
            // Series inductor → high-Z line: βl = g_k·Z₀/Z_high.
            (z_high, g_k * z0 / z_high)
        } else {
            // Shunt capacitor → low-Z line: βl = g_k·Z_low/Z₀.
            (z_low, g_k * z_low / z0)
        };

        // Physical width for this section's impedance (Hammerstad-Jensen).
        let width_m = microstrip_width(z_ohm, eps_r, h_m);
        // Guided wavelength at the section width: λ_g = c / (f_c·√ε_eff).
        let e_eff = eps_eff(width_m, h_m, eps_r);
        let lambda_g = C / (f_c_hz * e_eff.sqrt());
        // Physical length: l = (βl / 2π)·λ_g.
        let length_m = electrical_length_rad / (2.0 * std::f64::consts::PI) * lambda_g;

        sections.push(SteppedSection {
            high_z,
            z_ohm,
            electrical_length_rad,
            width_m,
            length_m,
        });
    }

    Ok(SteppedImpedanceDimensions {
        sections,
        eps_r,
        h_m,
    })
}

/// Convenience: assemble a [`yee_layout::Layout`] placing the synthesized
/// stepped-impedance sections **in-line** along `x`, source → load.
///
/// Each section is a width-`width_m` × length-`length_m` rectangle laid end to
/// end along `x`, centred on the `y = 0` axis (so the abrupt width steps are
/// symmetric about the line centre, as in a real stepped-impedance line). A
/// `feed_length` feed stub of the `Z₀` synthesis width attaches at each end,
/// with a `Z₀`-referenced [`yee_layout::PortRef`] at the two outer feed ends.
///
/// There is no dedicated in-line generator in `yee-layout` (the existing
/// `edge_coupled_bpf` / `hairpin_bpf` generators lay strips offset in `y` with a
/// single uniform width), so this composes the [`yee_layout`] primitives
/// directly rather than inventing a new generator. The feed length is one `Z₀`
/// guided quarter-wave at `f_c` (a neutral default); port → feed de-embedding is
/// out of scope for this increment.
///
/// # Errors
///
/// Propagates every [`DimError`] from [`dimension_stepped_impedance`].
pub fn dimension_stepped_impedance_layout(
    proto: &yee_synth::Prototype,
    f_c_hz: f64,
    z0: f64,
    z_high: f64,
    z_low: f64,
    sub: &Substrate,
) -> Result<Layout, DimError> {
    let dims = dimension_stepped_impedance(proto, f_c_hz, z0, z_high, z_low, sub)?;

    // Z0 feed line: synthesis width, a quarter guided wavelength long (neutral).
    let feed_width_m = microstrip_width(z0, sub.eps_r, sub.height_m);
    let feed_e_eff = eps_eff(feed_width_m, sub.height_m, sub.eps_r);
    let feed_length_m = C / (4.0 * f_c_hz * feed_e_eff.sqrt());

    let mut traces: Vec<Polygon> = Vec::with_capacity(dims.sections.len() + 2);

    // Input feed: extends leftward (−x) from the line start at x = 0.
    traces.push(Polygon::rect(
        -feed_length_m,
        -feed_width_m / 2.0,
        feed_length_m,
        feed_width_m,
    ));
    let in_port = PortRef {
        at: Point2::new(-feed_length_m, 0.0),
        width_m: feed_width_m,
        ref_impedance_ohm: z0,
    };

    // Lay the sections in-line along +x, each centred on y = 0.
    let mut x = 0.0_f64;
    for sec in &dims.sections {
        traces.push(Polygon::rect(
            x,
            -sec.width_m / 2.0,
            sec.length_m,
            sec.width_m,
        ));
        x += sec.length_m;
    }

    // Output feed: extends rightward (+x) from the line end at x.
    traces.push(Polygon::rect(
        x,
        -feed_width_m / 2.0,
        feed_length_m,
        feed_width_m,
    ));
    let out_port = PortRef {
        at: Point2::new(x + feed_length_m, 0.0),
        width_m: feed_width_m,
        ref_impedance_ohm: z0,
    };

    let bbox = BBox::from_polygons(&traces);
    Ok(Layout {
        substrate: *sub,
        traces,
        ports: vec![in_port, out_port],
        bbox,
    })
}

// ---------------------------------------------------------------------------
// Combline (capacitively-loaded short-circuited coupled lines) — F1.2.5.
// ---------------------------------------------------------------------------

/// First-order physical dimensions of a **combline** microstrip band-pass
/// filter, synthesized from a [`crate::CouplingMatrix`] (Filter Phase F1.2.5,
/// Hong & Lancaster §5.2.5).
///
/// A combline resonator is a short-circuited microstrip line of characteristic
/// impedance `Z0` and electrical length `θ0 < π/2` at `f0`, **capacitively
/// loaded** by a shunt capacitor `C_L` at its open end; the shorter-than-λ/4
/// line is brought to resonance at `f0` by that load. Adjacent resonators couple
/// through the line-to-line edge gap (coupled even/odd modes) — *exactly* the
/// edge-coupled / hairpin mechanism — so the coupling realization **reuses** the
/// validated [`solve_gap`] bisection and the `target_k = FBW · m_{i,i+1}`
/// derivation (see the [module docs](self)). The combline-**distinct** pieces
/// are the short-circuited `θ0` resonator and its loading cap; this struct
/// carries both alongside the shared `gaps_m` / `target_k`.
///
/// All lengths are in metres; `loading_cap_f` is in farads. `gaps_m` and
/// `target_k` are both length `N − 1` (one per adjacent resonator pair) and
/// index-aligned: `gaps_m[i]` is the edge gap that realizes `target_k[i]`.
///
/// This first-order engine reuses the proven `solve_gap` coupling realization
/// (like hairpin) rather than the rigorous Getsinger/Cristal self-/mutual-
/// capacitance coupled-bar synthesis (H&L eq 5.44); that, the discrete E-series
/// selection of `C_L`, and the via/short-circuit 3-D modelling are out of scope
/// for this increment (ADR-0144).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ComblineDimensions {
    /// Resonator / feed line width for the spec `Z0`, metres (Hammerstad-Jensen,
    /// via [`yee_layout::microstrip_width`]).
    pub line_width_m: f64,
    /// Chosen resonator electrical length `θ0` at `f0`, radians (must be in
    /// `(0, π/2)`; default design choice 45° = `π/4` = λ_g/8 for compactness).
    pub theta0_rad: f64,
    /// Physical resonator length `L = θ0 / β(f0)`, metres, with
    /// `β(f0) = 2π·f0·√ε_eff/c` at the synthesized width.
    pub resonator_length_m: f64,
    /// Loading capacitance `C_L = cot(θ0)/(2π·f0·Z0)`, farads, placed at the
    /// resonator's **open** end. The opposite end is **short-circuited** — a via
    /// to the ground plane (the short-circuit / via 3-D model is out of scope;
    /// here it is the ideal `Y → ∞` boundary the synthesis assumes).
    pub loading_cap_f: f64,
    /// Inter-resonator edge-coupling gaps, metres (length `N − 1`).
    pub gaps_m: Vec<f64>,
    /// The `FBW · m_{i,i+1}` coupling each gap was solved for (length `N − 1`).
    pub target_k: Vec<f64>,
}

/// Synthesize the physical dimensions of a **combline** microstrip band-pass
/// filter from a synthesized [`FilterProject`], a chosen resonator electrical
/// length `θ0`, and a [`Substrate`] (Filter Phase F1.2.5, Hong & Lancaster
/// §5.2.5).
///
/// Closed-form throughout and a direct mirror of [`dimension_hairpin`]; the
/// combline-distinct pieces are the short-circuited `θ0` resonator and its
/// loading cap:
///
/// - **Line width** — the spec-`Z0` Hammerstad-Jensen width
///   ([`yee_layout::microstrip_width`]).
/// - **Resonator length** — `L = θ0 / β(f0)` with `β(f0) = 2π·f0·√ε_eff/c`
///   (`ε_eff` from [`yee_layout::eps_eff`] at the synthesized width). A combline
///   resonator is a *short* (`θ0 < π/2`) short-circuited line, not the
///   edge-coupled λ_g/2 strip nor the hairpin's λ_g/4 arm.
/// - **Loading cap** — `C_L = cot(θ0)/(2π·f0·Z0)` (H&L eq 5.43): the shunt cap
///   at the open end that resonates the short-circuited `θ0` line at `f0`. The
///   short-circuited stub has input susceptance `B_stub = −(1/Z0)·cot(θ0·f/f0)`;
///   adding the cap's `2π·f·C_L` and forcing the sum to zero at `f = f0` gives
///   exactly this `C_L` (the `dim_combline_001` gate re-derives that resonance
///   independently rather than inverting this formula).
/// - **Inter-resonator gaps** — identical to edge-coupled / hairpin: for each
///   adjacent pair `(i, i+1)`, `target_k[i] = FBW · m_{i,i+1}` is realized by
///   bisecting the monotone coupled-line coupling coefficient with the shared
///   [`solve_gap`] helper (no optimizer, no FDTD). See the [module docs](self)
///   for the `target_k = FBW · m` cross-check.
///
/// # Errors
///
/// - [`DimError::InvalidTheta0`] if `θ0_rad` is not in the open interval
///   `(0, π/2)` (outside it `cot(θ0) ≤ 0` → a non-physical loading cap).
/// - [`DimError::UnsupportedTopology`] if the project is not
///   [`Topology::CoupledResonator`] (the only synthesized topology today; the
///   combline is a *realization* of that coupling network).
/// - [`DimError::OrderTooSmall`] if the order `N < 2` (no inter-resonator
///   coupling to realize).
/// - [`DimError::GapNotBracketed`] if a `target_k` is unreachable for any gap in
///   the `[5 µm, 5 mm]` bracket at the synthesized width (no silent clamping).
pub fn dimension_combline(
    project: &FilterProject,
    theta0_rad: f64,
    substrate: &Substrate,
) -> Result<ComblineDimensions, DimError> {
    // θ0 must be strictly inside (0, π/2): cot(θ0) ≤ 0 outside it yields a
    // non-physical (zero / negative) loading capacitance.
    if !(theta0_rad.is_finite() && theta0_rad > 0.0 && theta0_rad < std::f64::consts::FRAC_PI_2) {
        return Err(DimError::InvalidTheta0(theta0_rad));
    }

    if project.topology != Topology::CoupledResonator {
        return Err(DimError::UnsupportedTopology);
    }

    let n = project.coupling.m.len();
    if n < 2 {
        return Err(DimError::OrderTooSmall);
    }

    let eps_r = substrate.eps_r;
    let h_m = substrate.height_m;
    let f0 = project.spec.f0_hz;
    let fbw = project.spec.fbw;
    let z0 = project.spec.z0_ohm;

    // 1. Line width from the Hammerstad-Jensen Z0 synthesis.
    let line_width_m = microstrip_width(z0, eps_r, h_m);

    // 2. Resonator length L = θ0 / β(f0), β(f0) = 2π·f0·√ε_eff/c. The combline
    //    resonator is a short (θ0 < π/2) short-circuited line.
    let e_eff = eps_eff(line_width_m, h_m, eps_r);
    let beta0 = 2.0 * std::f64::consts::PI * f0 * e_eff.sqrt() / C;
    let resonator_length_m = theta0_rad / beta0;

    // 3. Loading cap C_L = cot(θ0)/(2π·f0·Z0) (H&L eq 5.43): the shunt cap at the
    //    open end that resonates the short-circuited θ0 line at f0.
    let cot_theta0 = theta0_rad.cos() / theta0_rad.sin();
    let loading_cap_f = cot_theta0 / (2.0 * std::f64::consts::PI * f0 * z0);

    // 4. Inter-resonator gaps: target_k[i] = FBW · m[i][i+1] (= yee-synth's
    //    k_{i,i+1}; see module docs), solved by the SAME bisection as edge-coupled
    //    / hairpin because adjacent combline resonators couple through the edge
    //    gap between their lines (coupled even/odd modes).
    let mut target_k = Vec::with_capacity(n - 1);
    let mut gaps_m = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let k_i = fbw * project.coupling.m[i][i + 1];
        let gap = solve_gap(i, k_i, line_width_m, h_m, eps_r)?;
        target_k.push(k_i);
        gaps_m.push(gap);
    }

    Ok(ComblineDimensions {
        line_width_m,
        theta0_rad,
        resonator_length_m,
        loading_cap_f,
        gaps_m,
        target_k,
    })
}

/// Compose a [`yee_layout::Layout`] for a **combline** band-pass filter (Filter
/// Phase F1.2.6), placing the synthesized [`dimension_combline`] dimensions as an
/// honest **comb**: `N` aligned, short-circuited resonator lines on a common
/// ground spine, capacitively loaded at the open ends, with tapped input/output
/// feeds.
///
/// This is the board-layout companion of [`dimension_combline`] — the
/// prerequisite for lighting combline in the studio (which renders a `Layout`).
/// It calls [`dimension_combline`] for the physics (no recompute) and composes
/// the comb from [`yee_layout`] primitives directly (`Polygon::rect` / `PortRef`
/// / `BBox::from_polygons` / `Layout`), the same approach
/// [`dimension_stepped_impedance_layout`] used — there is no first-class
/// `yee-layout::combline_bpf` generator (a later refactor; ADR-0145).
///
/// A combline is deliberately **not** drawn like [`edge_coupled_bpf`], which lays
/// staggered *open* half-wave lines: combline resonators are aligned, all
/// short-circuited at a common spine and capacitively loaded at the open ends.
/// The geometry composed here is:
///
/// - **Resonator lines** — `N` (= `gaps_m.len() + 1`) parallel vertical
///   rectangles, each `line_width_m` wide (x) × `resonator_length_m` long (y).
///   The short-circuit end is at `y = 0`, the open end at
///   `y = resonator_length_m`. Resonator `i`'s left edge sits at
///   `x_i = Σ_{j<i}(line_width_m + gaps_m[j])`, so the centre-to-centre pitch is
///   `line_width_m + gaps_m[i]` — the **real solved per-section gaps**, not a
///   uniform placeholder.
/// - **Ground spine** — a horizontal bar at the short-circuit end (`y` just ≤ 0,
///   `line_width_m` tall) spanning all `N` lines' x-range: the comb spine
///   (grounded via vias; vias are a fabrication annotation, not separate copper).
/// - **Loading-cap pads** — a small `line_width_m`-square pad at each open end
///   (`y = resonator_length_m`) where the SMD loading cap `C_L` mounts (the cap
///   value lives in the dimensions / studio table, not the copper).
/// - **Feeds + ports** — a tapped feed line to the first and last resonator
///   (neutral defaults, mirroring hairpin / stepped-Z: feed width =
///   `line_width_m`, feed length = one resonator length), each ending in a
///   [`PortRef`] referenced to the spec `Z0`. Mapping the external Q to a tap
///   position is deferred (as in [`dimension_hairpin_layout`]).
///
/// `bbox = BBox::from_polygons(&traces)`. No physics is recomputed and
/// `yee-layout` is not edited.
///
/// # Errors
///
/// Propagates every [`DimError`] from [`dimension_combline`].
pub fn dimension_combline_layout(
    project: &FilterProject,
    theta0_rad: f64,
    substrate: &Substrate,
) -> Result<Layout, DimError> {
    let dims = dimension_combline(project, theta0_rad, substrate)?;

    let n = dims.gaps_m.len() + 1; // N resonators, N−1 gaps.
    let w = dims.line_width_m;
    let l = dims.resonator_length_m;

    // Left-edge x of each resonator: x_0 = 0, x_i = x_{i-1} + w + gaps_m[i-1].
    // The centre-to-centre pitch is therefore `w + gaps_m[i-1]` — the solved
    // per-section gap, not a uniform placeholder.
    let mut resonator_x = Vec::with_capacity(n);
    let mut x = 0.0_f64;
    for i in 0..n {
        resonator_x.push(x);
        if i < dims.gaps_m.len() {
            x += w + dims.gaps_m[i];
        }
    }
    // Total x-extent spanned by the N resonator lines (left edge of #0 to right
    // edge of #(N−1)).
    let comb_right = resonator_x[n - 1] + w;

    // traces: N resonator lines + 1 ground spine + N cap pads + 2 feeds.
    let mut traces: Vec<Polygon> = Vec::with_capacity(2 * n + 3);

    // N resonator lines: short-circuit end at y = 0, open end at y = l.
    for &xi in &resonator_x {
        traces.push(Polygon::rect(xi, 0.0, w, l));
    }

    // Ground spine: a w-tall horizontal bar just below the short-circuit ends
    // (y ∈ [−w, 0]), spanning the full x-range of the resonators — the comb spine.
    traces.push(Polygon::rect(0.0, -w, comb_right, w));

    // Loading-cap pads: a w × w square at each open end (y ∈ [l, l + w]).
    for &xi in &resonator_x {
        traces.push(Polygon::rect(xi, l, w, w));
    }

    // Tapped input/output feeds + ports (neutral defaults, mirroring hairpin /
    // stepped-Z): feed width = line_width_m, feed length = one resonator length.
    // The input feed taps the first resonator's side edge partway up (at tap_y, a
    // neutral tap height) and extends in −x; the output feed taps the last
    // resonator's side edge and extends in +x. (A combline tap is up the resonator
    // from the grounded spine, not at the open / cap end.)
    let feed_width_m = w;
    let feed_length_m = l;
    // The tap height up the resonator (a neutral default; qe→tap dimensioning is
    // deferred, as in dimension_hairpin_layout).
    let tap_y = l / 3.0;

    // Input feed: extends leftward (−x) from the first resonator's left edge.
    let in_x0 = resonator_x[0] - feed_length_m;
    traces.push(Polygon::rect(
        in_x0,
        tap_y - feed_width_m / 2.0,
        feed_length_m,
        feed_width_m,
    ));
    let in_port = PortRef {
        at: Point2::new(in_x0, tap_y),
        width_m: feed_width_m,
        ref_impedance_ohm: project.spec.z0_ohm,
    };

    // Output feed: extends rightward (+x) from the last resonator's right edge.
    let out_x0 = resonator_x[n - 1] + w;
    traces.push(Polygon::rect(
        out_x0,
        tap_y - feed_width_m / 2.0,
        feed_length_m,
        feed_width_m,
    ));
    let out_port = PortRef {
        at: Point2::new(out_x0 + feed_length_m, tap_y),
        width_m: feed_width_m,
        ref_impedance_ohm: project.spec.z0_ohm,
    };

    let bbox = BBox::from_polygons(&traces);
    Ok(Layout {
        substrate: *substrate,
        traces,
        ports: vec![in_port, out_port],
        bbox,
    })
}

// ---------------------------------------------------------------------------
// Interdigital (short-circuited λg/4 coupled lines, no loading cap) — F1.2.7.
// ---------------------------------------------------------------------------

/// First-order physical dimensions of an **interdigital** microstrip band-pass
/// filter, synthesized from a [`crate::CouplingMatrix`] (Filter Phase F1.2.7,
/// Hong & Lancaster §5).
///
/// An interdigital resonator is a straight microstrip line of characteristic
/// impedance `Z0` that is **short-circuited at one end** and a **full quarter
/// guided wavelength (`λ_g/4`) long** — i.e. its electrical length at `f0` is
/// fixed at `θ = π/2`. Adjacent resonators are grounded at *alternating* ends
/// (the interdigital finger structure) and couple through the line-to-line edge
/// gap (coupled even/odd modes) — *exactly* the edge-coupled / hairpin / combline
/// mechanism — so the coupling realization **reuses** the validated [`solve_gap`]
/// bisection and the `target_k = FBW · m_{i,i+1}` derivation (see the
/// [module docs](self)).
///
/// The interdigital-**distinct** point is the resonator: it is the `θ = π/2`
/// limit of the combline line, where the short-circuited stub's input
/// susceptance `B(f0) = −(1/Z0)·cot(π/2) = 0` is **already zero** because the
/// full `λ_g/4` line is self-resonant. Consequently there is **no loading
/// capacitor** at all — the structural contrast with [`ComblineDimensions`],
/// which shortens the line to `θ0 < π/2` and adds `C_L = cot(θ0)/(2π·f0·Z0)`
/// precisely to make up the missing susceptance. This struct therefore carries
/// **neither** a `loading_cap_f` **nor** a `theta0_rad` field (θ is fixed at
/// π/2 by definition).
///
/// All lengths are in metres. `gaps_m` and `target_k` are both length `N − 1`
/// (one per adjacent resonator pair) and index-aligned: `gaps_m[i]` is the edge
/// gap between resonators `i` and `i + 1` that realizes `target_k[i]`.
///
/// This first-order engine reuses the proven `solve_gap` coupling realization
/// (like hairpin / combline) rather than the rigorous alternating-ground
/// even/odd-mode coupled-bar synthesis; that EM coupling refinement and the
/// via/short-circuit 3-D modelling are out of scope for this increment
/// (ADR-0148), exactly the scope boundary combline drew around its loading cap's
/// effect on coupling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InterdigitalDimensions {
    /// Resonator / feed line width for the spec `Z0`, metres (Hammerstad-Jensen,
    /// via [`yee_layout::microstrip_width`]).
    pub line_width_m: f64,
    /// Physical resonator length `L = λ_g/4 = (π/2)/β(f0)`, metres, with
    /// `β(f0) = 2π·f0·√ε_eff/c` at the synthesized width. The line is
    /// short-circuited at one end; the full quarter-wave length makes it
    /// self-resonant at `f0` with **no** loading cap.
    pub resonator_length_m: f64,
    /// Inter-resonator edge-coupling gaps, metres (length `N − 1`).
    pub gaps_m: Vec<f64>,
    /// The `FBW · m_{i,i+1}` coupling each gap was solved for (length `N − 1`).
    pub target_k: Vec<f64>,
}

/// Synthesize the physical dimensions of an **interdigital** microstrip
/// band-pass filter from a synthesized [`FilterProject`] and a [`Substrate`]
/// (Filter Phase F1.2.7, Hong & Lancaster §5).
///
/// Closed-form throughout and a direct mirror of [`dimension_hairpin`] /
/// [`dimension_combline`]; the interdigital-distinct point is the resonator — a
/// **short-circuited `λ_g/4` line, `θ = π/2` fixed, with no loading cap**. Unlike
/// [`dimension_combline`] this takes **no** `θ0` parameter (θ is fixed at π/2 by
/// definition):
///
/// - **Line width** — the spec-`Z0` Hammerstad-Jensen width
///   ([`yee_layout::microstrip_width`]).
/// - **Resonator length** — `resonator_length_m = λ_g/4 = (π/2)/β(f0)` with
///   `β(f0) = 2π·f0·√ε_eff/c` (`ε_eff` from [`yee_layout::eps_eff`] at the
///   synthesized width). This is the **factor-4** quarter-wave (like the
///   hairpin's λ_g/4 *arm*) but a *straight* line (not folded into a U), with
///   the far end short-circuited. Equivalently `c / (4·f0·√ε_eff)`.
/// - **No loading cap** — the `θ = π/2` line is self-resonant at `f0`
///   (`cot(π/2) = 0` → `B(f0) = 0`), so [`InterdigitalDimensions`] carries no
///   loading-cap field. This is the `θ = π/2` limit of combline, which instead
///   shortens the line and adds `C_L`; [`dimension_combline`] deliberately
///   *errors* at `θ0 = π/2`, so this is a genuinely distinct function.
/// - **Inter-resonator gaps** — identical to edge-coupled / hairpin / combline:
///   for each adjacent pair `(i, i+1)`, `target_k[i] = FBW · m_{i,i+1}` is
///   realized by bisecting the monotone coupled-line coupling coefficient with
///   the shared [`solve_gap`] helper (no optimizer, no FDTD). See the
///   [module docs](self) for the `target_k = FBW · m` cross-check.
///
/// # Errors
///
/// - [`DimError::UnsupportedTopology`] if the project is not
///   [`Topology::CoupledResonator`] (the only synthesized topology today; the
///   interdigital is a *realization* of that coupling network).
/// - [`DimError::OrderTooSmall`] if the order `N < 2` (no inter-resonator
///   coupling to realize).
/// - [`DimError::GapNotBracketed`] if a `target_k` is unreachable for any gap in
///   the `[5 µm, 5 mm]` bracket at the synthesized width (no silent clamping).
pub fn dimension_interdigital(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<InterdigitalDimensions, DimError> {
    if project.topology != Topology::CoupledResonator {
        return Err(DimError::UnsupportedTopology);
    }

    let n = project.coupling.m.len();
    if n < 2 {
        return Err(DimError::OrderTooSmall);
    }

    let eps_r = substrate.eps_r;
    let h_m = substrate.height_m;
    let f0 = project.spec.f0_hz;
    let fbw = project.spec.fbw;
    let z0 = project.spec.z0_ohm;

    // 1. Line width from the Hammerstad-Jensen Z0 synthesis.
    let line_width_m = microstrip_width(z0, eps_r, h_m);

    // 2. Resonator length L = λ_g/4 = (π/2)/β(f0), β(f0) = 2π·f0·√ε_eff/c. The
    //    interdigital resonator is a short-circuited *full* quarter-wave line
    //    (θ = π/2 fixed) — self-resonant at f0, so no loading cap. (Factor-4 like
    //    the hairpin arm, but a straight line.)
    let e_eff = eps_eff(line_width_m, h_m, eps_r);
    let beta0 = 2.0 * std::f64::consts::PI * f0 * e_eff.sqrt() / C;
    let resonator_length_m = std::f64::consts::FRAC_PI_2 / beta0;

    // 3. Inter-resonator gaps: target_k[i] = FBW · m[i][i+1] (= yee-synth's
    //    k_{i,i+1}; see module docs), solved by the SAME bisection as edge-coupled
    //    / hairpin / combline because adjacent interdigital resonators couple
    //    through the edge gap between their lines (coupled even/odd modes).
    let mut target_k = Vec::with_capacity(n - 1);
    let mut gaps_m = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let k_i = fbw * project.coupling.m[i][i + 1];
        let gap = solve_gap(i, k_i, line_width_m, h_m, eps_r)?;
        target_k.push(k_i);
        gaps_m.push(gap);
    }

    Ok(InterdigitalDimensions {
        line_width_m,
        resonator_length_m,
        gaps_m,
        target_k,
    })
}

/// Compose a [`yee_layout::Layout`] for an **interdigital** band-pass filter
/// (Filter Phase F1.2.8, Hong & Lancaster §5), placing the synthesized
/// [`dimension_interdigital`] dimensions as an honest interdigital **comb**: `N`
/// aligned, short-circuited `λ_g/4` resonator lines grounded at **alternating**
/// ends between **two** ground rails, with tapped input/output feeds.
///
/// This is the board-layout companion of [`dimension_interdigital`] — the
/// prerequisite for lighting interdigital in the studio (which renders a
/// `Layout`). It calls [`dimension_interdigital`] for the physics (no recompute)
/// and composes the comb from [`yee_layout`] primitives directly (`Polygon::rect`
/// / `PortRef` / `BBox::from_polygons` / `Layout`), exactly as
/// [`dimension_combline_layout`] does — there is no first-class
/// `yee-layout::interdigital_bpf` generator (ADR-0149). Unlike
/// [`dimension_combline_layout`] it takes **no** `θ0` parameter (interdigital is
/// `θ = π/2` fixed by definition).
///
/// An interdigital comb is the same aligned coupled-line comb at the solved
/// per-section gaps as combline, but differs in three concrete, drawable ways
/// (combline → interdigital):
///
/// 1. **Alternating-end grounding (the "finger" structure).** Combline shorts
///    *all* resonators at one common ground spine (`y = 0`). Interdigital shorts
///    adjacent resonators at **alternating** ends, so there are **two** ground
///    rails (bottom + top): the **bottom rail** (`y ∈ [−w, 0]`) grounds the
///    **even**-index resonators and the **top rail** (`y ∈ [l + g_open,
///    l + g_open + w]`) grounds the **odd**-index resonators. Each resonator's
///    open end is gapped `g_open` from the opposite rail; **no** resonator
///    touches both rails (an accidental short would make a cavity).
/// 2. **No loading-cap pads.** Combline draws a `w × w` cap pad at each open end
///    (the SMD `C_L` mounts there). The interdigital `λ_g/4` line is
///    self-resonant and needs **no** cap, so there are **no pads** — the trace
///    count is `N + 2 + 2` (N lines + 2 rails + 2 feeds), not combline's
///    `2N + 3`.
/// 3. **Full `λ_g/4` lines** — `resonator_length_m` from the engine (the
///    `θ = π/2` quarter-wave; combline's were the shortened `θ0 < π/2` line).
///
/// The geometry composed here is:
///
/// - **N resonator lines** — `N` (= `gaps_m.len() + 1`) parallel vertical
///   rectangles, each `line_width_m` wide (x) × `resonator_length_m` long (y),
///   **alternately offset** in `y` so the grounded end touches its rail and the
///   open end is gapped `g_open` from the opposite rail:
///   - **even `i`** (grounded bottom): `Polygon::rect(x_i, 0, w, l)` — shares the
///     `y = 0` edge with the bottom rail; open top at `y = l`, gapped `g_open`
///     below the top rail.
///   - **odd `i`** (grounded top): `Polygon::rect(x_i, g_open, w, l)` — top at
///     `y = l + g_open` shares the top rail's edge; open bottom at `y = g_open`,
///     gapped `g_open` above the bottom rail.
///
///   Resonator `i`'s left edge sits at `x_i = Σ_{j<i}(line_width_m + gaps_m[j])`,
///   so the centre-to-centre pitch is `line_width_m + gaps_m[i]` — the **real
///   solved per-section gaps**, not a uniform placeholder.
/// - **Bottom ground rail** — `Polygon::rect(0, −w, comb_right, w)`
///   (`y ∈ [−w, 0]`), spanning the comb x-range; grounds the even resonators
///   (vias are a fabrication annotation, not separate copper, as in combline).
/// - **Top ground rail** — `Polygon::rect(0, l + g_open, comb_right, w)`
///   (`y ∈ [l + g_open, l + g_open + w]`); grounds the odd resonators.
/// - **Feeds + ports** — a tapped feed line to the first (`i = 0`, grounded
///   bottom) and last (`i = N−1`) resonators (neutral defaults, mirroring
///   combline / hairpin: feed width = `line_width_m`, feed length = one resonator
///   length), tapped up the resonator from its grounded end at a neutral
///   `tap_y`, each ending in a [`PortRef`] referenced to the spec `Z0`. Mapping
///   the external Q to a tap position is deferred (F1.2.1, as in
///   [`dimension_hairpin_layout`]).
///
/// The open-end coupling gap `g_open = line_width_m` is a neutral fixed default
/// (like the hairpin fold spacing); mapping it to a precise end-coupling is an EM
/// follow-on, not first-order (ADR-0149).
///
/// `bbox = BBox::from_polygons(&traces)`. No physics is recomputed and
/// `yee-layout` is not edited.
///
/// # Errors
///
/// Propagates every [`DimError`] from [`dimension_interdigital`].
pub fn dimension_interdigital_layout(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<Layout, DimError> {
    let dims = dimension_interdigital(project, substrate)?;

    let n = dims.gaps_m.len() + 1; // N resonators, N−1 gaps.
    let w = dims.line_width_m;
    let l = dims.resonator_length_m;
    // Open-end coupling gap = a neutral fixed default (the hairpin fold spacing
    // analog); the precise end-coupling is an EM follow-on (ADR-0149).
    let g_open = w;

    // Left-edge x of each resonator: x_0 = 0, x_i = x_{i-1} + w + gaps_m[i-1].
    // The centre-to-centre pitch is therefore `w + gaps_m[i-1]` — the solved
    // per-section gap, not a uniform placeholder (identical to combline).
    let mut resonator_x = Vec::with_capacity(n);
    let mut x = 0.0_f64;
    for i in 0..n {
        resonator_x.push(x);
        if i < dims.gaps_m.len() {
            x += w + dims.gaps_m[i];
        }
    }
    // Total x-extent spanned by the N resonator lines (left edge of #0 to right
    // edge of #(N−1)).
    let comb_right = resonator_x[n - 1] + w;

    // traces: N resonator lines + 2 ground rails + 2 feeds (NO cap pads).
    let mut traces: Vec<Polygon> = Vec::with_capacity(n + 4);

    // N resonator lines, alternately offset (the interdigital finger structure):
    //   even i grounded at the bottom rail → y0 = 0 (open top at y = l),
    //   odd  i grounded at the top    rail → y0 = g_open (top at y = l + g_open).
    // Each line is the full λ_g/4 length l; no line touches both rails.
    for (i, &xi) in resonator_x.iter().enumerate() {
        let y0 = if i % 2 == 0 { 0.0 } else { g_open };
        traces.push(Polygon::rect(xi, y0, w, l));
    }

    // Bottom ground rail: a w-tall horizontal bar just below the bottom-grounded
    // ends (y ∈ [−w, 0]), spanning the full x-range — grounds the even resonators.
    traces.push(Polygon::rect(0.0, -w, comb_right, w));

    // Top ground rail: a w-tall horizontal bar just above the top-grounded ends
    // (y ∈ [l + g_open, l + g_open + w]), spanning the full x-range — grounds the
    // odd resonators. There is NO cap pad anywhere (the λ_g/4 line self-resonates).
    traces.push(Polygon::rect(0.0, l + g_open, comb_right, w));

    // Tapped input/output feeds + ports (neutral defaults, mirroring combline /
    // hairpin): feed width = line_width_m, feed length = one resonator length. The
    // tap is up the resonator from its grounded end at tap_y (a neutral tap
    // height; qe→tap dimensioning is deferred, F1.2.1). The first (i = 0) and last
    // (i = N−1) resonators are both grounded at the bottom rail when N is odd; tap
    // up from y = 0 in both cases (a neutral default — the interdigital tap is up
    // the line from the short, not at an open / cap end).
    let feed_width_m = w;
    let feed_length_m = l;
    let tap_y = l / 3.0;

    // Input feed: extends leftward (−x) from the first resonator's left edge.
    let in_x0 = resonator_x[0] - feed_length_m;
    traces.push(Polygon::rect(
        in_x0,
        tap_y - feed_width_m / 2.0,
        feed_length_m,
        feed_width_m,
    ));
    let in_port = PortRef {
        at: Point2::new(in_x0, tap_y),
        width_m: feed_width_m,
        ref_impedance_ohm: project.spec.z0_ohm,
    };

    // Output feed: extends rightward (+x) from the last resonator's right edge.
    let out_x0 = resonator_x[n - 1] + w;
    traces.push(Polygon::rect(
        out_x0,
        tap_y - feed_width_m / 2.0,
        feed_length_m,
        feed_width_m,
    ));
    let out_port = PortRef {
        at: Point2::new(out_x0 + feed_length_m, tap_y),
        width_m: feed_width_m,
        ref_impedance_ohm: project.spec.z0_ohm,
    };

    let bbox = BBox::from_polygons(&traces);
    Ok(Layout {
        substrate: *substrate,
        traces,
        ports: vec![in_port, out_port],
        bbox,
    })
}
