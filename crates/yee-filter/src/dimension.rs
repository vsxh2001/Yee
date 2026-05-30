//! Closed-form microstrip dimensional synthesis (Filter Phases F1.2.0 / F1.2.2).
//!
//! Turns an abstract synthesized [`crate::CouplingMatrix`] into **physical
//! microstrip dimensions** for two band-pass topologies — the edge-coupled
//! half-wave filter ([`dimension_edge_coupled`], F1.2.0) and the U-folded
//! **hairpin** filter ([`dimension_hairpin`], F1.2.2) — by inverting the
//! already-validated `yee-layout` closed-form models. Pure `f64`, WASM-safe, NO
//! FDTD, NO surrogate — this is the *initial* dimensioning that seeds the later
//! EM-in-the-loop refinement (F1.2.1).
//!
//! Both topologies share the **same inter-resonator coupling mechanism**: a
//! hairpin is a half-wave line folded into a U, and adjacent hairpins couple
//! through the edge gap between their adjacent arms — exactly the edge-coupled
//! gap→`k` inversion. The two paths therefore reuse the identical
//! [`solve_gap`] bisection and the `target_k = FBW · m_{i,i+1}` derivation
//! below; only the resonator geometry differs (a folded half-wave = two ≈λ/4
//! arms vs a single λ/2 straight strip — see [`dimension_hairpin`]).
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
    EdgeCoupledParams, EdgeCoupledSection, HairpinParams, Layout, Substrate, coupled_microstrip,
    coupling_coefficient, edge_coupled_bpf, eps_eff, hairpin_bpf, microstrip_width,
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
        }
    }
}

impl std::error::Error for DimError {}

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

    // 2. Arm length = λ_g/4 = c / (4·f0·√ε_eff). The hairpin half-wave resonator
    //    is folded into a U, so it is two ≈λ/4 arms (factor-4) vs the edge-coupled
    //    straight λ/2 strip (factor-2). (Hong & Lancaster ch. 6.)
    let e_eff = eps_eff(line_width_m, h_m, eps_r);
    let arm_length_m = C / (4.0 * f0 * e_eff.sqrt());

    // 3. Fold spacing: intra-hairpin arm pitch — a fixed closed-form choice (not
    //    the inter-resonator coupling), refined later by F1.2.1 BO.
    let fold_spacing_m = HAIRPIN_FOLD_SPACING_WIDTHS * line_width_m;

    // 4. Inter-resonator gaps: target_k[i] = FBW · m[i][i+1] (= yee-synth's
    //    k_{i,i+1}; see module docs), solved by the SAME bisection as edge-coupled
    //    because adjacent hairpins couple through the edge gap between their arms.
    let mut target_k = Vec::with_capacity(n - 1);
    let mut gaps_m = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let k_i = fbw * project.coupling.m[i][i + 1];
        let gap = solve_gap(i, k_i, line_width_m, h_m, eps_r)?;
        target_k.push(k_i);
        gaps_m.push(gap);
    }

    Ok(HairpinDimensions {
        line_width_m,
        arm_length_m,
        fold_spacing_m,
        gaps_m,
        target_k,
    })
}

/// Convenience: assemble a [`yee_layout::Layout`] from the synthesized hairpin
/// dimensions via the existing [`yee_layout::hairpin_bpf`].
///
/// Builds the `N` U-folded resonators (all of width `line_width_m`, arm length
/// `arm_length_m`, arm pitch `fold_spacing_m`) with a tapped feed.
///
/// **Uniform-gap walking-skeleton limitation (gap option (b)).**
/// [`HairpinParams`] today carries a *single* `coupling_gap_m`, which
/// [`yee_layout::hairpin_bpf`] bakes into a uniform resonator pitch — it has no
/// per-section gap field. Synthesis, however, produces `N − 1` *distinct* gaps
/// (one per coupling `k_{i,i+1}`). Extending `hairpin_bpf` to per-section gaps
/// would rework the generator's coordinate math and perturb the committed
/// `geo-003` geometry gate, so this skeleton instead passes a **single
/// representative gap** — the mean of the solved [`HairpinDimensions::gaps_m`] —
/// and documents the limitation. The per-section `gaps_m` are still returned in
/// full by [`dimension_hairpin`] (and round-trip-validated by `hairpin_dim_001`),
/// so the synthesis fidelity is unaffected; only the convenience `Layout` here is
/// uniform-gap. A per-section `hairpin_bpf` (gap option (a)) and the `qe`→tap
/// dimensioning are both deferred to F1.2.1.
///
/// The tapped-feed geometry uses neutral defaults: `tap_offset_m` is a third of
/// the arm length, `feed_width_m = line_width_m`, and `feed_length_m` is one arm
/// length. Mapping the external Q (`qe_in`/`qe_out`) to a tap position is
/// **deferred to F1.2.1**; this function does **not** invent a `qe`→tap formula.
///
/// # Errors
///
/// Propagates every [`DimError`] from [`dimension_hairpin`].
pub fn dimension_hairpin_layout(
    project: &FilterProject,
    substrate: &Substrate,
) -> Result<Layout, DimError> {
    let dims = dimension_hairpin(project, substrate)?;

    let n = dims.gaps_m.len() + 1; // N resonators, N−1 gaps.

    // Uniform-gap walking skeleton (gap option (b)): hairpin_bpf takes a single
    // coupling_gap_m, so collapse the N−1 distinct solved gaps to their mean. The
    // full per-section gaps stay in `dims.gaps_m`; see the doc-comment above.
    let representative_gap_m = dims.gaps_m.iter().sum::<f64>() / dims.gaps_m.len() as f64;

    let params = HairpinParams {
        substrate: *substrate,
        n,
        arm_length_m: dims.arm_length_m,
        line_width_m: dims.line_width_m,
        fold_spacing_m: dims.fold_spacing_m,
        coupling_gap_m: representative_gap_m,
        // Neutral tapped-feed defaults; qe→tap dimensioning is F1.2.1.
        tap_offset_m: dims.arm_length_m / 3.0,
        feed_width_m: dims.line_width_m,
        feed_length_m: dims.arm_length_m,
    };

    Ok(hairpin_bpf(&params))
}
