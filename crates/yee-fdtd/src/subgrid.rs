//! Subgridding scaffold — Phase 2.fdtd.7.0 (Q2).
//!
//! Ships the type surface for a single axis-aligned, cuboidal **fine
//! sub-region** nested at 2× resolution inside a parent [`YeeGrid`], plus a
//! solver wrapper that owns it. This is the *scaffold* step: the fine grid is
//! allocated, sized, and inherits the parent's scalar `eps_r` / `mu_r`, but
//! the coupling logic between coarse and fine grids (spatial / temporal
//! interpolation of `E_t`, area-averaging of `H_t`) lands in subsequent
//! steps Q3, Q4, Q5.
//!
//! At this step [`SubgriddedSolver::step`] is a placeholder that delegates
//! straight to the wrapped [`WalkingSkeletonSolver`]'s helper sequence —
//! the fine grid is dormant. That keeps the type surface stable so
//! downstream tracks can be developed against it in parallel while the
//! actual interleave is added under the same module.
//!
//! ## Phase 2.fdtd.7.0 scope (this module)
//!
//! - Single nested region.
//! - 2× isotropic refinement: `dx_fine = dx_coarse / 2`, ditto `dy`, `dz`,
//!   `dt`.
//! - Axis-aligned, cuboidal.
//! - Non-dispersive, isotropic, lossless materials inherited from the parent.
//! - **No co-location with CPML thickness or a TF/SF box face.** Co-location
//!   is a documented runtime error in 7.0 (see [`SubgridContext`]).
//!
//! ## References
//!
//! - Spec: `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`
//! - Plan: `docs/superpowers/plans/2026-05-18-phase-2-fdtd-7-subgridding.md`,
//!   Step Q2.
//! - Chevalier, M. W., Luebbers, R. J., Cable, V. P., "FDTD local grid with
//!   material traverse", *IEEE Trans. Antennas Propag.* 45(3), 1997.

use ndarray::{Array2, Array3};
use yee_core::Error;
use yee_core::units::{C0, EPS0, MU0};

use crate::FdtdSolver;
use crate::WalkingSkeletonSolver;
use crate::grid::YeeGrid;
use crate::update;

/// Coarse-cell `(lo, hi)` extent of an axis-aligned box, inclusive-low /
/// exclusive-high. Used for the TF/SF-box placement check.
pub type CoarseBox = ((usize, usize, usize), (usize, usize, usize));

/// Which bracketing snapshot — start- or end-of-coarse-step — a copy
/// targets. Internal helper for the Q3 snapshot pair.
#[derive(Debug, Clone, Copy)]
enum SnapshotKind {
    /// Parent `E_t` just before the coarse E-update.
    Start,
    /// Parent `E_t` just after the coarse E-update.
    End,
}

/// 2D linear interpolation on the unit square with weights
/// `w_a, w_b ∈ [0, 1]` along the two axes.
///
/// The four samples are `f00, f01, f10, f11` where the first index is the
/// `w_a` axis and the second is the `w_b` axis (so `f01` sits at
/// `(w_a, w_b) = (0, 1)`).
#[inline]
fn bilerp(f00: f64, f01: f64, f10: f64, f11: f64, w_a: f64, w_b: f64) -> f64 {
    let lo = (1.0 - w_b) * f00 + w_b * f01;
    let hi = (1.0 - w_b) * f10 + w_b * f11;
    (1.0 - w_a) * lo + w_a * hi
}

/// Optional context for validating the placement of a [`SubgridRegion`]
/// against parent-domain features (CPML thickness, TF/SF box) that are not
/// stored on the bare [`YeeGrid`].
///
/// `Default` is all-`None`, which disables the corresponding checks. Callers
/// that *do* know the parent's CPML thickness or TF/SF box should populate
/// the matching field so [`SubgridRegion::new_with_context`] can enforce the
/// "fine region is interior to CPML and TF/SF" invariant from spec §6.
#[derive(Debug, Default, Clone, Copy)]
pub struct SubgridContext {
    /// CPML layer thickness in coarse cells on every outer face. If `Some`,
    /// the constructor rejects subgrid regions whose `lo`/`hi` overlap any
    /// CPML cell on any axis.
    pub cpml_thickness: Option<usize>,
    /// Inclusive-low / exclusive-high coarse-cell coordinates of a TF/SF
    /// box. If `Some`, the constructor rejects subgrid regions that
    /// straddle (cross) any face of the TF/SF box. A region entirely
    /// inside or entirely outside the TF/SF box is permitted; the failure
    /// mode this guards against is a face *intersecting* the nest, which
    /// breaks the TF/SF reciprocity argument.
    pub tfsf_box: Option<CoarseBox>,
}

/// Identifier for one of the six Huygens-surface faces of a cuboidal
/// [`SubgridRegion`], used by the Berenger 2006 fine → coarse closure to
/// enumerate per-face equivalent-current injections.
///
/// The naming convention is `<axis><Low|High>` where the axis is the
/// direction of the outward unit normal `n̂`:
///
/// - `XLow`  — outward normal `−x̂`, coarse cells with `i_c = lo.0 − 1`
///   (the cell layer of the parent grid just outside the fine box).
/// - `XHigh` — outward normal `+x̂`, coarse cells with `i_c = hi.0`.
/// - `YLow`  — outward normal `−ŷ`, `j_c = lo.1 − 1`.
/// - `YHigh` — outward normal `+ŷ`, `j_c = hi.1`.
/// - `ZLow`  — outward normal `−ẑ`, `k_c = lo.2 − 1`.
/// - `ZHigh` — outward normal `+ẑ`, `k_c = hi.2`.
///
/// The axis index (`X = 0`, `Y = 1`, `Z = 2`) is used by the
/// edge-ownership rule (lower-numbered axis wins the shared edge) — see
/// [`SubgridRegion::face_index_table`].
///
/// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`.
/// ADR: `docs/src/decisions/0035-berenger-huygens-subgridding.md`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BerengerHuygensFace {
    /// `−x̂` face (outward normal points in the −x direction).
    XLow,
    /// `+x̂` face.
    XHigh,
    /// `−ŷ` face.
    YLow,
    /// `+ŷ` face.
    YHigh,
    /// `−ẑ` face.
    ZLow,
    /// `+ẑ` face.
    ZHigh,
}

impl BerengerHuygensFace {
    /// Axis index (`X = 0`, `Y = 1`, `Z = 2`) that this face's outward
    /// normal is parallel to. Used by the lower-axis-wins edge-ownership
    /// rule — see [`SubgridRegion::face_index_table`] and
    /// [`assign_edge_to_face`].
    pub fn axis(self) -> usize {
        match self {
            Self::XLow | Self::XHigh => 0,
            Self::YLow | Self::YHigh => 1,
            Self::ZLow | Self::ZHigh => 2,
        }
    }

    /// All six face identifiers, ordered `[XLow, XHigh, YLow, YHigh, ZLow,
    /// ZHigh]`. Stable across calls; matches the layout returned by
    /// [`SubgridRegion::face_index_table`].
    pub fn all() -> [Self; 6] {
        [
            Self::XLow,
            Self::XHigh,
            Self::YLow,
            Self::YHigh,
            Self::ZLow,
            Self::ZHigh,
        ]
    }
}

/// Given two faces of the cuboidal Huygens surface that share an edge,
/// return the face that **owns** the edge under the spec's
/// "lower-numbered-axis-wins" rule (see
/// `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`
/// §3 risks register, item 2).
///
/// Two faces share an edge iff their axes differ; if they have the same
/// axis they are either identical (`XLow` ≡ `XLow`) or opposite (`XLow`,
/// `XHigh`) and do not share an edge. In those degenerate cases this
/// helper returns the first argument unchanged — callers that care
/// should not invoke it on a same-axis pair (the unit tests enumerate
/// only the 12 axis-distinct cuboid edges).
///
/// Example:
///
/// ```text
/// assign_edge_to_face(XLow, YLow)  == XLow   (X axis 0 < Y axis 1)
/// assign_edge_to_face(YLow, ZLow)  == YLow   (Y axis 1 < Z axis 2)
/// assign_edge_to_face(XHigh, ZLow) == XHigh  (X axis 0 < Z axis 2)
/// ```
pub fn assign_edge_to_face(
    face_a: BerengerHuygensFace,
    face_b: BerengerHuygensFace,
) -> BerengerHuygensFace {
    if face_a.axis() <= face_b.axis() {
        face_a
    } else {
        face_b
    }
}

/// Per-face start/end snapshots of the parent grid's tangential E field on
/// the six interface planes of a [`SubgridRegion`].
///
/// Each face owns two 2D arrays (one per tangential `E` component). The
/// arrays are sized to *match the parent's natural-grid samples* on the
/// interface plane (one coarse cell per coarse edge), not the fine-grid
/// resolution — spatial interpolation onto the fine grid happens inside
/// [`SubgridRegion::interpolate_coarse_e_to_fine`].
///
/// "Start" is the parent `E_t` at the beginning of the coarse step (before
/// `parent.update_e_only`); "end" is the parent `E_t` after the coarse
/// E-update. Linear-in-time blending between them gives the fine-grid
/// boundary `E_t` at fractional time `frac ∈ (0, 1)`.
///
/// Cf. Chevalier 1997 §III — spatial part — and Okoniewski 1997 — the
/// 2× temporal-subcycling pattern that motivates the two-sample blend.
#[derive(Debug, Clone, Default)]
pub struct InterfaceSnapshots {
    /// `(start, end)` snapshots of `E_y` on the ±x faces.
    /// Shape `[hi.1 - lo.1, hi.2 - lo.2 + 1]` — one entry per coarse
    /// `E_y` edge incident on the face plane.
    pub xmin_ey: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_z` on the −x face.
    /// Shape `[hi.1 - lo.1 + 1, hi.2 - lo.2]`.
    pub xmin_ez: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_y` on the +x face.
    pub xmax_ey: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_z` on the +x face.
    pub xmax_ez: (Array2<f64>, Array2<f64>),

    /// `(start, end)` snapshots of `E_x` on the −y face.
    /// Shape `[hi.0 - lo.0, hi.2 - lo.2 + 1]`.
    pub ymin_ex: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_z` on the −y face.
    /// Shape `[hi.0 - lo.0 + 1, hi.2 - lo.2]`.
    pub ymin_ez: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_x` on the +y face.
    pub ymax_ex: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_z` on the +y face.
    pub ymax_ez: (Array2<f64>, Array2<f64>),

    /// `(start, end)` snapshots of `E_x` on the −z face.
    /// Shape `[hi.0 - lo.0, hi.1 - lo.1 + 1]`.
    pub zmin_ex: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_y` on the −z face.
    /// Shape `[hi.0 - lo.0 + 1, hi.1 - lo.1]`.
    pub zmin_ey: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_x` on the +z face.
    pub zmax_ex: (Array2<f64>, Array2<f64>),
    /// `(start, end)` snapshots of `E_y` on the +z face.
    pub zmax_ey: (Array2<f64>, Array2<f64>),
}

/// Axis-aligned, cuboidal sub-region nested at 2× resolution inside a parent
/// [`YeeGrid`].
///
/// Owns its own fine `YeeGrid` instance whose cell sizes (`dx`, `dy`, `dz`)
/// and time step (`dt`) are half the parent's, sized
/// `(2·(hi.0 − lo.0), 2·(hi.1 − lo.1), 2·(hi.2 − lo.2))` cells. The fine
/// grid inherits the parent's scalar `eps_r` and `mu_r`.
///
/// Carries `InterfaceSnapshots` (Q3) for coarse → fine `E_t` interpolation
/// during fine sub-steps. The fine → coarse `H_t` area-averaging closure
/// (Q4) and the seven-stage `step` (Q5) plug into this same state.
#[derive(Debug, Clone)]
pub struct SubgridRegion {
    /// Coarse-cell index of the nest corner (inclusive lower bound).
    pub lo: (usize, usize, usize),
    /// Coarse-cell index of the nest corner (exclusive upper bound).
    pub hi: (usize, usize, usize),
    /// The fine grid backing this region. `dx_fine = dx_coarse / 2`;
    /// `dt_fine = dt_coarse / 2`.
    fine: YeeGrid,
    /// Bracketing parent `E_t` snapshots on the six interface faces,
    /// blended linearly in time during each fine sub-step.
    snapshots: InterfaceSnapshots,
    /// Mid-coarse-step snapshot of the fine `H` field, taken between the
    /// two fine sub-steps in the Q5 seven-stage `step`. Populated by
    /// [`SubgridRegion::snapshot_fine_h_mid_step`] at fine wall-clock time
    /// `t = n + 1/4` (i.e. right after sub-step 1's `update_fine_h`); the
    /// post-sub-step-2 fine H (at `t = n + 3/4`) is averaged against this
    /// snapshot in [`SubgridRegion::average_fine_h_to_coarse`] to recover
    /// the time-centered value `t = n + 1/2` the coarse `H_t` slot
    /// represents. `None` until the first snapshot of a coarse step;
    /// callers that invoke `average_fine_h_to_coarse` directly (without
    /// stepping) get the legacy single-sample behaviour.
    fine_h_snapshot: Option<FineHSnapshot>,
    /// End-of-coarse-step snapshot of the fine `E` field, captured at
    /// `t = n + 1` (after sub-step 2's `update_fine_e`) and consumed at
    /// the **start of the next coarse step** to inject the
    /// `M = -n̂ × E_tot` magnetic-current source onto the coarse `H`
    /// arrays just before that step's `update_h_only`.
    ///
    /// `None` until the first call to
    /// [`SubgridRegion::snapshot_fine_e_end_of_step`]; the first coarse
    /// step's `inject_m_to_coarse_h` call therefore treats the source
    /// as identically zero (correct: no fine fields existed at `t = 0`
    /// before sub-step 2 ran).
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md` §3.
    /// ADR:  `docs/src/decisions/0035-berenger-huygens-subgridding.md`.
    fine_e_snapshot: Option<FineESnapshot>,
    /// Phase 2.fdtd.7.y Option β `E_pre` snapshot — fine `E_t` captured in
    /// sub-step 2 immediately after [`Self::interpolate_coarse_e_to_fine`]
    /// at `frac = 0.75` writes the Q3 Dirichlet value and **before** that
    /// sub-step's [`Self::update_fine_e`] runs. Paired with
    /// [`Self::fine_e_post_snapshot`] to compute the compensating M source
    /// `M = -n̂ × (E_post − E_pre)` (consumed by
    /// [`Self::inject_m_to_coarse_h`] starting in Phase 2.fdtd.7.y Step C2;
    /// populated-but-unread in C1).
    ///
    /// `None` until [`Self::snapshot_fine_e_pre_update`] has been called at
    /// least once.
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3.
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
    fine_e_pre_snapshot: Option<FineESnapshot>,
    /// Phase 2.fdtd.7.y Option β `E_post` snapshot — fine `E_t` captured
    /// immediately after sub-step 2's [`Self::update_fine_e`] completes.
    /// Paired with [`Self::fine_e_pre_snapshot`]; see that field's docs
    /// for the compensating-source rationale.
    ///
    /// `None` until [`Self::snapshot_fine_e_post_update`] has been called
    /// at least once.
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3.
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
    fine_e_post_snapshot: Option<FineESnapshot>,
    /// Phase 2.fdtd.7.y Step C5 (Option α — Mur ABC) per-face cache of
    /// the fine-grid tangential `E` values at the **boundary** and the
    /// **first-cell-inside** layers, captured **immediately before** a
    /// fine sub-step's [`Self::update_fine_e`] runs. Consumed
    /// **immediately after** that `update_fine_e` by
    /// [`Self::apply_mur_abc_to_fine_outer_e`] to evaluate the Mur
    /// 1st-order absorbing boundary update
    /// `E_t^{n+1}(boundary) = E_t^n(inside) + α · (E_t^{n+1}(inside) −
    /// E_t^n(boundary))` with `α = (c·dt_fine − dx_fine) /
    /// (c·dt_fine + dx_fine)`.
    ///
    /// Boundary-vs-inside layer indices per face (in fine-cell space):
    ///
    /// - `−x` face: boundary `i_f = 0`, inside `i_f = 1`.
    /// - `+x` face: boundary `i_f = fine_nx`, inside `i_f = fine_nx − 1`.
    /// - Similarly on `±y`, `±z`.
    ///
    /// `None` until [`Self::snapshot_fine_e_for_mur`] has been called at
    /// least once. This snapshot replaces the Q3 Dirichlet
    /// [`Self::interpolate_coarse_e_to_fine`] writes that previously fed
    /// the fine outer `E_t` layer; see ADR-0038 Option α (escape hatch)
    /// for the rationale.
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3 Option α.
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md` "Consequences".
    /// Reference: Mur, G., "Absorbing Boundary Conditions for the
    /// Finite-Difference Approximation of the Time-Domain
    /// Electromagnetic-Field Equations", IEEE Trans. EMC EMC-23(4),
    /// 1981, pp. 377–382, eq. 5 — 1st-order Mur on a Cartesian Yee
    /// boundary.
    fine_e_mur_snapshot: Option<FineEMurSnapshot>,
}

/// Mid-coarse-step snapshot of the fine `H` field components used by the
/// Q4.1 time-centered fine → coarse closure. Holds full-fine-grid copies
/// of `H_x`, `H_y`, `H_z`; the closure only reads cells adjacent to the
/// six interface faces, so the full-grid clone is a small constant-factor
/// memory cost (≈ same size as the fine grid's own H arrays) that keeps
/// the snapshot self-contained and avoids per-face slice plumbing.
#[derive(Debug, Clone)]
struct FineHSnapshot {
    hx: Array3<f64>,
    hy: Array3<f64>,
    hz: Array3<f64>,
}

/// End-of-coarse-step snapshot of the fine `E` field components used by
/// the Phase 2.fdtd.7.x B2.1 split Berenger closure to defer the
/// `M = -n̂ × E_tot` magnetic-current injection to the **start of the
/// next coarse step**, just before the next coarse `update_h_only`.
///
/// Captured at `t = n + 1` (immediately after sub-step 2's
/// `update_fine_e`) by [`SubgridRegion::snapshot_fine_e_end_of_step`].
/// Consumed at the top of the next coarse step by
/// [`SubgridRegion::inject_m_to_coarse_h`]. Decouples the `M` injection
/// from the J injection so each source enters its respective leapfrog
/// update before that update runs, per the Phase 2.fdtd.7.x B2.1 fix
/// for the spec §6 risk 3 time-centering bug diagnosed by Track
/// HHHHHHH (commit `997e706` divergence; the J/M source values were
/// being stacked on top of already-updated coarse `E`/`H` slots,
/// closing a unit-magnitude feedback loop through the Q3 coarse → fine
/// Dirichlet boundary).
#[derive(Debug, Clone)]
struct FineESnapshot {
    ex: Array3<f64>,
    ey: Array3<f64>,
    ez: Array3<f64>,
}

/// Phase 2.fdtd.7.y Step C5 (Option α — Mur ABC) per-face cache used to
/// evaluate the 1st-order Mur absorbing boundary update on the fine
/// outer `E_t` plane. Each per-face `(Array2, Array2)` tuple holds, in
/// order:
///
/// - `.0` — the fine `E_t` value at the **boundary** cell (i.e. on the
///   outer Huygens surface) at fine wall-clock time `t = n` (sampled
///   right before [`SubgridRegion::update_fine_e`] runs).
/// - `.1` — the fine `E_t` value at the **first-cell-inside** layer
///   (one fine cell inward of the boundary along the face normal) at
///   the same fine wall-clock time `t = n`.
///
/// After `update_fine_e` advances the inside layer to `t = n + 1` (the
/// boundary cell is untouched — `crate::update::update_e` skips outer
/// tangential faces, see the loop ranges in that module), the Mur
/// update reads `E_t^n(inside) = .1`, `E_t^n(boundary) = .0`,
/// `E_t^{n+1}(inside) = grid.e?[(boundary +/- 1, …)]` (live), and writes
/// the new boundary value back into `grid.e?[(boundary, …)]`.
///
/// Component layout per face mirrors the existing
/// [`InterfaceSnapshots`]-derived shapes but in **fine-cell** index
/// space (every coarse cell is doubled in each tangential direction):
///
/// - `±x` faces — tangential components are `E_y` and `E_z`.
/// - `±y` faces — tangential components are `E_x` and `E_z`.
/// - `±z` faces — tangential components are `E_x` and `E_y`.
///
/// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3 Option α.
/// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
#[derive(Debug, Clone)]
struct FineEMurSnapshot {
    /// −x face: `(boundary, inside)` pairs for `E_y` and `E_z`.
    xmin_ey: (Array2<f64>, Array2<f64>),
    xmin_ez: (Array2<f64>, Array2<f64>),
    /// +x face.
    xmax_ey: (Array2<f64>, Array2<f64>),
    xmax_ez: (Array2<f64>, Array2<f64>),
    /// −y face: `(boundary, inside)` pairs for `E_x` and `E_z`.
    ymin_ex: (Array2<f64>, Array2<f64>),
    ymin_ez: (Array2<f64>, Array2<f64>),
    /// +y face.
    ymax_ex: (Array2<f64>, Array2<f64>),
    ymax_ez: (Array2<f64>, Array2<f64>),
    /// −z face: `(boundary, inside)` pairs for `E_x` and `E_y`.
    zmin_ex: (Array2<f64>, Array2<f64>),
    zmin_ey: (Array2<f64>, Array2<f64>),
    /// +z face.
    zmax_ex: (Array2<f64>, Array2<f64>),
    zmax_ey: (Array2<f64>, Array2<f64>),
}

impl SubgridRegion {
    /// Build a 2× nest covering coarse cells `lo..hi` of `parent`.
    ///
    /// Performs the *base* validity checks — `lo < hi` on every axis and
    /// `hi` inside the parent grid's cell count. Callers that need to
    /// additionally check the region against CPML thickness or a TF/SF
    /// box should call [`Self::new_with_context`] with a populated
    /// [`SubgridContext`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] if `lo >= hi` on any axis or if any
    /// component of `hi` exceeds the matching parent dimension.
    pub fn new(
        parent: &YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
    ) -> Result<Self, Error> {
        Self::new_with_context(parent, lo, hi, SubgridContext::default())
    }

    /// Like [`Self::new`], but additionally enforces the documented
    /// runtime-error cases from spec §6: the region must not overlap the
    /// CPML thickness on any face, and must not cross any face of the
    /// supplied TF/SF box.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Invalid`] for the base checks (see [`Self::new`])
    /// or for either of the co-location checks when the corresponding
    /// field of `ctx` is populated.
    pub fn new_with_context(
        parent: &YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        ctx: SubgridContext,
    ) -> Result<Self, Error> {
        Self::check_bounds(parent, lo, hi)?;

        if let Some(npml) = ctx.cpml_thickness {
            Self::check_cpml_disjoint(parent, lo, hi, npml)?;
        }
        if let Some((tlo, thi)) = ctx.tfsf_box {
            Self::check_tfsf_disjoint(lo, hi, tlo, thi)?;
        }

        let fine_nx = 2 * (hi.0 - lo.0);
        let fine_ny = 2 * (hi.1 - lo.1);
        let fine_nz = 2 * (hi.2 - lo.2);

        // Parent is uniform-cell (dy = dz = dx per YeeGrid::vacuum); inherit
        // the parent's dx halved and let YeeGrid::vacuum derive its own dt
        // from the resulting Courant limit, then overwrite with exactly
        // parent.dt / 2 so the 2× temporal-subcycle invariant is exact to
        // f64::EPSILON.
        let dx_fine = parent.dx / 2.0;
        let mut fine = YeeGrid::vacuum(fine_nx, fine_ny, fine_nz, dx_fine);
        fine.dy = parent.dy / 2.0;
        fine.dz = parent.dz / 2.0;
        fine.dt = parent.dt / 2.0;
        fine.eps_r = parent.eps_r;
        fine.mu_r = parent.mu_r;

        let snapshots = Self::allocate_snapshots(lo, hi);
        Ok(Self {
            lo,
            hi,
            fine,
            snapshots,
            fine_h_snapshot: None,
            fine_e_snapshot: None,
            fine_e_pre_snapshot: None,
            fine_e_post_snapshot: None,
            fine_e_mur_snapshot: None,
        })
    }

    /// Allocate zero-initialised parent `E_t` snapshot buffers sized to
    /// match the coarse-grid sample count on each of the six interface
    /// faces. Shapes follow the parent's natural Yee staggering, not the
    /// fine resolution.
    fn allocate_snapshots(
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
    ) -> InterfaceSnapshots {
        let nx_c = hi.0 - lo.0;
        let ny_c = hi.1 - lo.1;
        let nz_c = hi.2 - lo.2;

        let face_x_ey = (
            Array2::<f64>::zeros((ny_c, nz_c + 1)),
            Array2::<f64>::zeros((ny_c, nz_c + 1)),
        );
        let face_x_ez = (
            Array2::<f64>::zeros((ny_c + 1, nz_c)),
            Array2::<f64>::zeros((ny_c + 1, nz_c)),
        );

        let face_y_ex = (
            Array2::<f64>::zeros((nx_c, nz_c + 1)),
            Array2::<f64>::zeros((nx_c, nz_c + 1)),
        );
        let face_y_ez = (
            Array2::<f64>::zeros((nx_c + 1, nz_c)),
            Array2::<f64>::zeros((nx_c + 1, nz_c)),
        );

        let face_z_ex = (
            Array2::<f64>::zeros((nx_c, ny_c + 1)),
            Array2::<f64>::zeros((nx_c, ny_c + 1)),
        );
        let face_z_ey = (
            Array2::<f64>::zeros((nx_c + 1, ny_c)),
            Array2::<f64>::zeros((nx_c + 1, ny_c)),
        );

        InterfaceSnapshots {
            xmin_ey: face_x_ey.clone(),
            xmin_ez: face_x_ez.clone(),
            xmax_ey: face_x_ey,
            xmax_ez: face_x_ez,
            ymin_ex: face_y_ex.clone(),
            ymin_ez: face_y_ez.clone(),
            ymax_ex: face_y_ex,
            ymax_ez: face_y_ez,
            zmin_ex: face_z_ex.clone(),
            zmin_ey: face_z_ey.clone(),
            zmax_ex: face_z_ex,
            zmax_ey: face_z_ey,
        }
    }

    /// Immutable borrow of the fine [`YeeGrid`] backing this region.
    pub fn fine_grid(&self) -> &YeeGrid {
        &self.fine
    }

    /// Mutable borrow of the fine [`YeeGrid`] backing this region.
    ///
    /// Escape hatch for callers that need to write into the fine grid's
    /// material state or seed initial fields before stepping.
    pub fn fine_grid_mut(&mut self) -> &mut YeeGrid {
        &mut self.fine
    }

    /// Apply the bulk H-field Yee update to the fine grid.
    ///
    /// Companion to [`WalkingSkeletonSolver::update_h_only`], specialised
    /// to the fine sub-grid. Pure spatial-curl update with no boundary,
    /// no source, no clock advance — the fine grid's "boundary" is the
    /// Dirichlet `E_t` written by [`Self::interpolate_coarse_e_to_fine`]
    /// just before this call.
    pub fn update_fine_h(&mut self) {
        update::update_h(&mut self.fine);
    }

    /// Apply the bulk E-field Yee update to the fine grid.
    ///
    /// Companion to [`Self::update_fine_h`]; closes the fine half-step.
    /// Phase 2.fdtd.7.0 carries no CPML or PEC closure on the fine grid
    /// (the fine region sits strictly interior to any coarse CPML); the
    /// Dirichlet boundary `E_t` is fixed by the prior interpolation.
    pub fn update_fine_e(&mut self) {
        update::update_e(&mut self.fine);
    }

    /// Snapshot the current fine `H` field for the Q4.1 time-centered
    /// fine → coarse closure.
    ///
    /// Call once **between** the two fine sub-steps in the Q5 seven-stage
    /// `step` — specifically, after sub-step 1's
    /// [`Self::update_fine_h`] (which lands fine H at wall-clock time
    /// `t = n + 1/4`) and before sub-step 2's `update_fine_h`. The snapshot
    /// is then averaged against the post-sub-step-2 fine H (at `t = n + 3/4`)
    /// inside [`Self::average_fine_h_to_coarse`] to recover the
    /// time-centered value `t = n + 1/2` that the coarse `H_t` slot
    /// represents, eliminating the quarter-step phase error otherwise
    /// fed into the coarse H closure each coarse step.
    pub fn snapshot_fine_h_mid_step(&mut self) {
        self.fine_h_snapshot = Some(FineHSnapshot {
            hx: self.fine.hx.clone(),
            hy: self.fine.hy.clone(),
            hz: self.fine.hz.clone(),
        });
    }

    /// Snapshot the post-sub-step-2 fine `E` field for the Phase
    /// 2.fdtd.7.x B2.1 split-injection closure.
    ///
    /// Call once at the **end** of each coarse step — after sub-step 2's
    /// [`Self::update_fine_e`] has advanced fine `E` to wall-clock time
    /// `t = n + 1`, and before the next coarse step begins. The snapshot
    /// is then consumed at the **top** of the next coarse step by
    /// [`Self::inject_m_to_coarse_h`] to apply the
    /// `M = -n̂ × E_tot^{n+1}` magnetic-current source onto the coarse
    /// `H_t` arrays *before* the next coarse `update_h_only` consumes
    /// them. This defers the M injection by one coarse step so the
    /// source enters its respective leapfrog update before the update
    /// runs — the spec §6 risk 3 time-centering fix per Phase 2.fdtd.7.x
    /// B2.1.
    pub fn snapshot_fine_e_end_of_step(&mut self) {
        self.fine_e_snapshot = Some(FineESnapshot {
            ex: self.fine.ex.clone(),
            ey: self.fine.ey.clone(),
            ez: self.fine.ez.clone(),
        });
    }

    /// Snapshot the fine `E_t` field **before** sub-step 2's
    /// [`Self::update_fine_e`] runs — i.e. immediately after sub-step 2's
    /// [`Self::interpolate_coarse_e_to_fine`] has written the Q3 Dirichlet
    /// value onto the outer fine layer. Companion to
    /// [`Self::snapshot_fine_e_post_update`]; the pair forms the Phase
    /// 2.fdtd.7.y Option β compensating M source
    /// `M = -n̂ × (E_post − E_pre)`.
    ///
    /// Time-level diagram (sub-step 2 within one coarse step `n`):
    ///
    /// ```text
    ///   coarse E^n                          coarse E^{n+1}
    ///       │  Q3 Dirichlet                       │
    ///       └─► interpolate_coarse_e_to_fine(0.75)
    ///              │
    ///              ▼   E_pre captured here (post-Q3, pre-update_fine_e)
    ///       snapshot_fine_e_pre_update   ◄── this method
    ///              │
    ///              ▼
    ///         update_fine_h
    ///              │
    ///              ▼
    ///         update_fine_e   (fine E advances by one fine dt curl-of-H step)
    ///              │
    ///              ▼   E_post captured here (post-update_fine_e)
    ///       snapshot_fine_e_post_update
    /// ```
    ///
    /// By construction `E_pre ≈ time-interpolated coarse E_surface` (the
    /// Q3 Dirichlet tie); `E_post − E_pre` then isolates the
    /// Maxwell-evolved correction added by the fine sub-step's
    /// `update_fine_e`, which is the physically non-zero quantity the
    /// compensating source draws on.
    ///
    /// Phase 2.fdtd.7.y Step C1 — snapshot is captured but not yet
    /// consumed by [`Self::inject_m_to_coarse_h`] (that wires up in
    /// Step C2).
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3.
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
    pub fn snapshot_fine_e_pre_update(&mut self) {
        self.fine_e_pre_snapshot = Some(FineESnapshot {
            ex: self.fine.ex.clone(),
            ey: self.fine.ey.clone(),
            ez: self.fine.ez.clone(),
        });
    }

    /// Snapshot the fine `E_t` field **after** sub-step 2's
    /// [`Self::update_fine_e`] completes. Companion to
    /// [`Self::snapshot_fine_e_pre_update`]; see that method's docstring
    /// for the time-level diagram.
    ///
    /// Time-level diagram (sub-step 2 within one coarse step `n`):
    ///
    /// ```text
    ///         snapshot_fine_e_pre_update          (E_pre = post-Q3 Dirichlet)
    ///              │
    ///              ▼
    ///         update_fine_h
    ///              │
    ///              ▼
    ///         update_fine_e   (fine E advances by one fine dt curl-of-H step)
    ///              │
    ///              ▼   E_post captured here
    ///       snapshot_fine_e_post_update   ◄── this method
    /// ```
    ///
    /// Together with `E_pre`, this pair feeds the Phase 2.fdtd.7.y
    /// Option β compensating M source
    /// `M = -n̂ × (E_post − E_pre)`; the difference isolates the
    /// Maxwell-evolved part of fine E that escapes the Q3 Dirichlet tie
    /// and discards the Dirichlet-tied part that would otherwise
    /// nullify against the coarse ghost in Berenger's canonical form
    /// (see ADR-0038 / OOOOOOO's empirical regression).
    ///
    /// Phase 2.fdtd.7.y Step C1 — snapshot is captured but not yet
    /// consumed by [`Self::inject_m_to_coarse_h`] (that wires up in
    /// Step C2). Distinct from [`Self::snapshot_fine_e_end_of_step`],
    /// which targets the deferred B2.1 M source `M = -n̂ × E_fine^{n+1}`
    /// at the *end* of the coarse step rather than the compensating
    /// per-sub-step delta.
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3.
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
    pub fn snapshot_fine_e_post_update(&mut self) {
        self.fine_e_post_snapshot = Some(FineESnapshot {
            ex: self.fine.ex.clone(),
            ey: self.fine.ey.clone(),
            ez: self.fine.ez.clone(),
        });
    }

    /// Immutable borrow of the Phase 2.fdtd.7.y Option β `E_pre` snapshot
    /// (`ex`, `ey`, `ez`), populated by
    /// [`Self::snapshot_fine_e_pre_update`]. Returns `None` if no
    /// snapshot has been taken yet.
    ///
    /// Diagnostic/test accessor for Phase 2.fdtd.7.y Step C1; the
    /// production consumer ([`Self::inject_m_to_coarse_h`] in Step C2)
    /// reaches into the field directly.
    pub fn fine_e_pre_snapshot(&self) -> Option<(&Array3<f64>, &Array3<f64>, &Array3<f64>)> {
        self.fine_e_pre_snapshot
            .as_ref()
            .map(|s| (&s.ex, &s.ey, &s.ez))
    }

    /// Immutable borrow of the Phase 2.fdtd.7.y Option β `E_post`
    /// snapshot (`ex`, `ey`, `ez`), populated by
    /// [`Self::snapshot_fine_e_post_update`]. Returns `None` if no
    /// snapshot has been taken yet.
    ///
    /// Diagnostic/test accessor for Phase 2.fdtd.7.y Step C1.
    pub fn fine_e_post_snapshot(&self) -> Option<(&Array3<f64>, &Array3<f64>, &Array3<f64>)> {
        self.fine_e_post_snapshot
            .as_ref()
            .map(|s| (&s.ex, &s.ey, &s.ez))
    }

    // ----------------------------------------------------------------
    // Phase 2.fdtd.7.y Step C5 (Option α) — 1st-order Mur ABC on the
    // fine outer E_t plane. Replaces the Q3 coarse → fine Dirichlet
    // tie that previously held the outer fine `E_t` slaved to the
    // coarse-side interpolation. The Q3 helper
    // [`Self::interpolate_coarse_e_to_fine`] is retained for tests and
    // as a `#[doc(hidden)]` rollback path; production callers
    // ([`SubgriddedSolver::step`] / `step_with_gaussian_source_ez`)
    // route through the Mur helpers below instead.
    //
    // Reference: Mur, G., "Absorbing Boundary Conditions for the
    // Finite-Difference Approximation of the Time-Domain
    // Electromagnetic-Field Equations", IEEE Trans. EMC EMC-23(4),
    // 1981, pp. 377–382, eq. 5. The 1st-order Mur formula for an
    // outward-propagating plane wave hitting a Yee grid boundary
    // at `i = i_max` (analogous on the other 5 faces) reads
    //
    //   E_t^{n+1}(i_max, j, k) = E_t^n(i_max − 1, j, k)
    //                          + α · (E_t^{n+1}(i_max − 1, j, k)
    //                                 − E_t^n(i_max, j, k))
    //
    // with `α = (c·dt − dx) / (c·dt + dx)`. `α → 0` as `c·dt → dx`
    // (the lattice CFL limit) — the boundary then becomes a perfect
    // 1-cell upwind extrapolation. For the standard 0.9 × Courant
    // step used in this crate `α ≈ −0.0526`, which yields the
    // documented ~−40 dB reflection floor on plane waves at normal
    // incidence (spec §3 Option α trade-off).
    // ----------------------------------------------------------------

    /// Snapshot the fine `E_t` field on the **boundary** and on the
    /// **first-cell-inside** layer of each of the six outer Huygens
    /// faces, at fine wall-clock time `t = n`. Call **immediately
    /// before** [`Self::update_fine_e`] runs.
    ///
    /// The companion [`Self::apply_mur_abc_to_fine_outer_e`] consumes
    /// this snapshot **immediately after** that same `update_fine_e`
    /// completes (when the inside layer has advanced to `t = n + 1`
    /// while the boundary cell, which `crate::update::update_e` skips,
    /// still holds `t = n`). The Mur update then writes the new
    /// `t = n + 1` boundary value back into `self.fine.e?`.
    ///
    /// Pairs of snapshots are kept per face per tangential component:
    /// each pair `(boundary, inside)` is a 2-D fine-cell array sized to
    /// match the corresponding face's tangential extent.
    ///
    /// Phase 2.fdtd.7.y Step C5 (Option α) — replaces the Q3
    /// coarse → fine Dirichlet write that previously tied the fine
    /// outer `E_t` to the coarse-side interpolation. See ADR-0038
    /// "Consequences" and spec §3 Option α for the rationale (the Q3
    /// tie made Berenger's canonical M source
    /// `M = -n̂ × (E_TF_fine − E_SF_coarse_ghost)` degenerate to zero
    /// because `update_fine_e` skips boundary cells; replacing Q3 with
    /// a Mur ABC lets the fine outer `E_t` evolve as a function of
    /// adjacent-inside fine `E_t` instead, recovering non-zero
    /// `E_post − E_pre` differencing for the compensating-source M
    /// path).
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3 Option α.
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
    pub fn snapshot_fine_e_for_mur(&mut self) {
        let fine_nx = 2 * (self.hi.0 - self.lo.0);
        let fine_ny = 2 * (self.hi.1 - self.lo.1);
        let fine_nz = 2 * (self.hi.2 - self.lo.2);
        let f = &self.fine;

        // ±x faces — tangential `E_y` (shape `[*, fine_ny, fine_nz+1]`)
        // and `E_z` (shape `[*, fine_ny+1, fine_nz]`).
        let mut xmin_ey_b = Array2::<f64>::zeros((fine_ny, fine_nz + 1));
        let mut xmin_ey_i = Array2::<f64>::zeros((fine_ny, fine_nz + 1));
        let mut xmax_ey_b = Array2::<f64>::zeros((fine_ny, fine_nz + 1));
        let mut xmax_ey_i = Array2::<f64>::zeros((fine_ny, fine_nz + 1));
        for j in 0..fine_ny {
            for k in 0..=fine_nz {
                xmin_ey_b[(j, k)] = f.ey[(0, j, k)];
                xmin_ey_i[(j, k)] = f.ey[(1, j, k)];
                xmax_ey_b[(j, k)] = f.ey[(fine_nx, j, k)];
                xmax_ey_i[(j, k)] = f.ey[(fine_nx - 1, j, k)];
            }
        }
        let mut xmin_ez_b = Array2::<f64>::zeros((fine_ny + 1, fine_nz));
        let mut xmin_ez_i = Array2::<f64>::zeros((fine_ny + 1, fine_nz));
        let mut xmax_ez_b = Array2::<f64>::zeros((fine_ny + 1, fine_nz));
        let mut xmax_ez_i = Array2::<f64>::zeros((fine_ny + 1, fine_nz));
        for j in 0..=fine_ny {
            for k in 0..fine_nz {
                xmin_ez_b[(j, k)] = f.ez[(0, j, k)];
                xmin_ez_i[(j, k)] = f.ez[(1, j, k)];
                xmax_ez_b[(j, k)] = f.ez[(fine_nx, j, k)];
                xmax_ez_i[(j, k)] = f.ez[(fine_nx - 1, j, k)];
            }
        }

        // ±y faces — tangential `E_x` (shape `[fine_nx, *, fine_nz+1]`)
        // and `E_z` (shape `[fine_nx+1, *, fine_nz]`).
        let mut ymin_ex_b = Array2::<f64>::zeros((fine_nx, fine_nz + 1));
        let mut ymin_ex_i = Array2::<f64>::zeros((fine_nx, fine_nz + 1));
        let mut ymax_ex_b = Array2::<f64>::zeros((fine_nx, fine_nz + 1));
        let mut ymax_ex_i = Array2::<f64>::zeros((fine_nx, fine_nz + 1));
        for i in 0..fine_nx {
            for k in 0..=fine_nz {
                ymin_ex_b[(i, k)] = f.ex[(i, 0, k)];
                ymin_ex_i[(i, k)] = f.ex[(i, 1, k)];
                ymax_ex_b[(i, k)] = f.ex[(i, fine_ny, k)];
                ymax_ex_i[(i, k)] = f.ex[(i, fine_ny - 1, k)];
            }
        }
        let mut ymin_ez_b = Array2::<f64>::zeros((fine_nx + 1, fine_nz));
        let mut ymin_ez_i = Array2::<f64>::zeros((fine_nx + 1, fine_nz));
        let mut ymax_ez_b = Array2::<f64>::zeros((fine_nx + 1, fine_nz));
        let mut ymax_ez_i = Array2::<f64>::zeros((fine_nx + 1, fine_nz));
        for i in 0..=fine_nx {
            for k in 0..fine_nz {
                ymin_ez_b[(i, k)] = f.ez[(i, 0, k)];
                ymin_ez_i[(i, k)] = f.ez[(i, 1, k)];
                ymax_ez_b[(i, k)] = f.ez[(i, fine_ny, k)];
                ymax_ez_i[(i, k)] = f.ez[(i, fine_ny - 1, k)];
            }
        }

        // ±z faces — tangential `E_x` (shape `[fine_nx, fine_ny+1, *]`)
        // and `E_y` (shape `[fine_nx+1, fine_ny, *]`).
        let mut zmin_ex_b = Array2::<f64>::zeros((fine_nx, fine_ny + 1));
        let mut zmin_ex_i = Array2::<f64>::zeros((fine_nx, fine_ny + 1));
        let mut zmax_ex_b = Array2::<f64>::zeros((fine_nx, fine_ny + 1));
        let mut zmax_ex_i = Array2::<f64>::zeros((fine_nx, fine_ny + 1));
        for i in 0..fine_nx {
            for j in 0..=fine_ny {
                zmin_ex_b[(i, j)] = f.ex[(i, j, 0)];
                zmin_ex_i[(i, j)] = f.ex[(i, j, 1)];
                zmax_ex_b[(i, j)] = f.ex[(i, j, fine_nz)];
                zmax_ex_i[(i, j)] = f.ex[(i, j, fine_nz - 1)];
            }
        }
        let mut zmin_ey_b = Array2::<f64>::zeros((fine_nx + 1, fine_ny));
        let mut zmin_ey_i = Array2::<f64>::zeros((fine_nx + 1, fine_ny));
        let mut zmax_ey_b = Array2::<f64>::zeros((fine_nx + 1, fine_ny));
        let mut zmax_ey_i = Array2::<f64>::zeros((fine_nx + 1, fine_ny));
        for i in 0..=fine_nx {
            for j in 0..fine_ny {
                zmin_ey_b[(i, j)] = f.ey[(i, j, 0)];
                zmin_ey_i[(i, j)] = f.ey[(i, j, 1)];
                zmax_ey_b[(i, j)] = f.ey[(i, j, fine_nz)];
                zmax_ey_i[(i, j)] = f.ey[(i, j, fine_nz - 1)];
            }
        }

        self.fine_e_mur_snapshot = Some(FineEMurSnapshot {
            xmin_ey: (xmin_ey_b, xmin_ey_i),
            xmin_ez: (xmin_ez_b, xmin_ez_i),
            xmax_ey: (xmax_ey_b, xmax_ey_i),
            xmax_ez: (xmax_ez_b, xmax_ez_i),
            ymin_ex: (ymin_ex_b, ymin_ex_i),
            ymin_ez: (ymin_ez_b, ymin_ez_i),
            ymax_ex: (ymax_ex_b, ymax_ex_i),
            ymax_ez: (ymax_ez_b, ymax_ez_i),
            zmin_ex: (zmin_ex_b, zmin_ex_i),
            zmin_ey: (zmin_ey_b, zmin_ey_i),
            zmax_ex: (zmax_ex_b, zmax_ex_i),
            zmax_ey: (zmax_ey_b, zmax_ey_i),
        });
    }

    /// Apply the 1st-order Mur absorbing boundary update to the fine
    /// outer `E_t` plane on every one of the six Huygens faces. Call
    /// **immediately after** [`Self::update_fine_e`] completes, with
    /// the matching [`Self::snapshot_fine_e_for_mur`] having been
    /// captured **immediately before** that same `update_fine_e`.
    ///
    /// Mur (1981) eq. 5 for the `+x` boundary (analogous on the other
    /// five faces, with the appropriate "inside" index offset of `−1`
    /// on `−` faces and `+1` on `+` faces — i.e. always one cell
    /// *inward* of the boundary node):
    ///
    /// ```text
    /// E_t^{n+1}(i_max, j, k) = E_t^n(i_max−1, j, k)
    ///                        + α · (E_t^{n+1}(i_max−1, j, k)
    ///                               − E_t^n(i_max,   j, k))
    /// ```
    ///
    /// with `α = (c·dt_fine − dx_fine) / (c·dt_fine + dx_fine)`. The
    /// inside cell's new value at `t = n + 1` is read live from the
    /// fine grid (it was advanced by the preceding `update_fine_e`);
    /// the two `E_t^n` reads come from the snapshot.
    ///
    /// `c` is the vacuum speed of light; the fine grid carries the
    /// parent's `eps_r`, `mu_r` so the in-medium phase velocity is
    /// `c / sqrt(eps_r · mu_r)`. For the Phase 2.fdtd.7 walking
    /// skeleton both relative parameters are uniform vacuum (= 1), so
    /// `c` is the bare vacuum value [`yee_core::units::C0`].
    ///
    /// No-op if no snapshot has been taken yet (defensive guard for
    /// the first call after construction or for test paths that drive
    /// `update_fine_e` directly without first taking the Mur
    /// snapshot). The pre-Mur step in [`SubgriddedSolver::step`] /
    /// `step_with_gaussian_source_ez` always pairs snapshot + apply
    /// around `update_fine_e`.
    ///
    /// Phase 2.fdtd.7.y Step C5 (Option α). Spec / ADR references on
    /// [`Self::snapshot_fine_e_for_mur`].
    pub fn apply_mur_abc_to_fine_outer_e(&mut self) {
        let Some(snap) = self.fine_e_mur_snapshot.as_ref() else {
            return;
        };

        let fine_nx = 2 * (self.hi.0 - self.lo.0);
        let fine_ny = 2 * (self.hi.1 - self.lo.1);
        let fine_nz = 2 * (self.hi.2 - self.lo.2);

        // Mur 1st-order coefficient. The fine grid is uniform-cubic
        // (dx_fine = dy_fine = dz_fine) per `SubgridRegion::new`, so a
        // single `α` per axis is correct. Compute per-axis anyway for
        // robustness if anisotropic fine cells are ever wired (each
        // face uses its own normal-axis `dx`).
        let dt_f = self.fine.dt;
        let c_eff = C0 / (self.fine.eps_r * self.fine.mu_r).sqrt();
        let alpha_x = (c_eff * dt_f - self.fine.dx) / (c_eff * dt_f + self.fine.dx);
        let alpha_y = (c_eff * dt_f - self.fine.dy) / (c_eff * dt_f + self.fine.dy);
        let alpha_z = (c_eff * dt_f - self.fine.dz) / (c_eff * dt_f + self.fine.dz);

        let f = &mut self.fine;

        // ±x faces (Mur along x): outer `E_y` and `E_z`.
        // `−x`: boundary `i = 0`, inside `i = 1`.
        // `+x`: boundary `i = fine_nx`, inside `i = fine_nx − 1`.
        for j in 0..fine_ny {
            for k in 0..=fine_nz {
                let e_inside_new = f.ey[(1, j, k)];
                let e_inside_old = snap.xmin_ey.1[(j, k)];
                let e_bdry_old = snap.xmin_ey.0[(j, k)];
                f.ey[(0, j, k)] = e_inside_old + alpha_x * (e_inside_new - e_bdry_old);

                let e_inside_new = f.ey[(fine_nx - 1, j, k)];
                let e_inside_old = snap.xmax_ey.1[(j, k)];
                let e_bdry_old = snap.xmax_ey.0[(j, k)];
                f.ey[(fine_nx, j, k)] = e_inside_old + alpha_x * (e_inside_new - e_bdry_old);
            }
        }
        for j in 0..=fine_ny {
            for k in 0..fine_nz {
                let e_inside_new = f.ez[(1, j, k)];
                let e_inside_old = snap.xmin_ez.1[(j, k)];
                let e_bdry_old = snap.xmin_ez.0[(j, k)];
                f.ez[(0, j, k)] = e_inside_old + alpha_x * (e_inside_new - e_bdry_old);

                let e_inside_new = f.ez[(fine_nx - 1, j, k)];
                let e_inside_old = snap.xmax_ez.1[(j, k)];
                let e_bdry_old = snap.xmax_ez.0[(j, k)];
                f.ez[(fine_nx, j, k)] = e_inside_old + alpha_x * (e_inside_new - e_bdry_old);
            }
        }

        // ±y faces (Mur along y): outer `E_x` and `E_z`.
        for i in 0..fine_nx {
            for k in 0..=fine_nz {
                let e_inside_new = f.ex[(i, 1, k)];
                let e_inside_old = snap.ymin_ex.1[(i, k)];
                let e_bdry_old = snap.ymin_ex.0[(i, k)];
                f.ex[(i, 0, k)] = e_inside_old + alpha_y * (e_inside_new - e_bdry_old);

                let e_inside_new = f.ex[(i, fine_ny - 1, k)];
                let e_inside_old = snap.ymax_ex.1[(i, k)];
                let e_bdry_old = snap.ymax_ex.0[(i, k)];
                f.ex[(i, fine_ny, k)] = e_inside_old + alpha_y * (e_inside_new - e_bdry_old);
            }
        }
        for i in 0..=fine_nx {
            for k in 0..fine_nz {
                let e_inside_new = f.ez[(i, 1, k)];
                let e_inside_old = snap.ymin_ez.1[(i, k)];
                let e_bdry_old = snap.ymin_ez.0[(i, k)];
                f.ez[(i, 0, k)] = e_inside_old + alpha_y * (e_inside_new - e_bdry_old);

                let e_inside_new = f.ez[(i, fine_ny - 1, k)];
                let e_inside_old = snap.ymax_ez.1[(i, k)];
                let e_bdry_old = snap.ymax_ez.0[(i, k)];
                f.ez[(i, fine_ny, k)] = e_inside_old + alpha_y * (e_inside_new - e_bdry_old);
            }
        }

        // ±z faces (Mur along z): outer `E_x` and `E_y`.
        for i in 0..fine_nx {
            for j in 0..=fine_ny {
                let e_inside_new = f.ex[(i, j, 1)];
                let e_inside_old = snap.zmin_ex.1[(i, j)];
                let e_bdry_old = snap.zmin_ex.0[(i, j)];
                f.ex[(i, j, 0)] = e_inside_old + alpha_z * (e_inside_new - e_bdry_old);

                let e_inside_new = f.ex[(i, j, fine_nz - 1)];
                let e_inside_old = snap.zmax_ex.1[(i, j)];
                let e_bdry_old = snap.zmax_ex.0[(i, j)];
                f.ex[(i, j, fine_nz)] = e_inside_old + alpha_z * (e_inside_new - e_bdry_old);
            }
        }
        for i in 0..=fine_nx {
            for j in 0..fine_ny {
                let e_inside_new = f.ey[(i, j, 1)];
                let e_inside_old = snap.zmin_ey.1[(i, j)];
                let e_bdry_old = snap.zmin_ey.0[(i, j)];
                f.ey[(i, j, 0)] = e_inside_old + alpha_z * (e_inside_new - e_bdry_old);

                let e_inside_new = f.ey[(i, j, fine_nz - 1)];
                let e_inside_old = snap.zmax_ey.1[(i, j)];
                let e_bdry_old = snap.zmax_ey.0[(i, j)];
                f.ey[(i, j, fine_nz)] = e_inside_old + alpha_z * (e_inside_new - e_bdry_old);
            }
        }
    }

    /// Base bounds validation: `lo < hi` per axis, `hi` inside the parent.
    fn check_bounds(
        parent: &YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
    ) -> Result<(), Error> {
        if lo.0 >= hi.0 || lo.1 >= hi.1 || lo.2 >= hi.2 {
            return Err(Error::Invalid(format!(
                "SubgridRegion: lo must be strictly less than hi on every axis, got lo={lo:?}, hi={hi:?}"
            )));
        }
        if hi.0 > parent.nx || hi.1 > parent.ny || hi.2 > parent.nz {
            return Err(Error::Invalid(format!(
                "SubgridRegion: hi={:?} out of parent bounds (nx={}, ny={}, nz={})",
                hi, parent.nx, parent.ny, parent.nz
            )));
        }
        Ok(())
    }

    /// Read-only borrow of the cached interface E-field snapshots. Useful
    /// for diagnostics, energy-balance probes, and the Q5 step driver.
    pub fn snapshots(&self) -> &InterfaceSnapshots {
        &self.snapshots
    }

    // ----------------------------------------------------------------
    // Q3 — coarse → fine E_t spatial + temporal interpolation
    // ----------------------------------------------------------------

    /// Cache the parent grid's tangential `E` field on the six interface
    /// faces as the **start-of-coarse-step** snapshot.
    ///
    /// Call once at the top of each coarse step before `parent.update_e_only`.
    /// Pair with [`Self::snapshot_coarse_e_t_end`] after the coarse E-update
    /// so the two snapshots bracket the time interval over which the fine
    /// sub-steps interpolate.
    pub fn snapshot_coarse_e_t(&mut self, parent: &YeeGrid) {
        Self::copy_face_e_t(
            &mut self.snapshots,
            parent,
            self.lo,
            self.hi,
            SnapshotKind::Start,
        );
    }

    /// Cache the parent grid's tangential `E` field on the six interface
    /// faces as the **end-of-coarse-step** snapshot.
    ///
    /// Call once after `parent.update_e_only` (and any source / CPML / PEC
    /// closure on the coarse grid) so the cached pair brackets the coarse
    /// E-update interval. Linear blending between the two snapshots at the
    /// fine sub-step fractions `frac ∈ {0.25, 0.75}` yields the Dirichlet
    /// fine boundary `E_t` per Okoniewski 1997.
    pub fn snapshot_coarse_e_t_end(&mut self, parent: &YeeGrid) {
        Self::copy_face_e_t(
            &mut self.snapshots,
            parent,
            self.lo,
            self.hi,
            SnapshotKind::End,
        );
    }

    /// Write Dirichlet `E_t` values on the six outer fine-grid faces by
    /// blending the start/end coarse snapshots in time at fraction `frac`
    /// (typically `0.25` for the first fine sub-step and `0.75` for the
    /// second) and linearly interpolating in space between bracketing
    /// coarse edges per Chevalier 1997 §III.
    ///
    /// `frac` is clamped to `[0, 1]`; values outside the unit interval
    /// would imply extrapolation across a coarse interval, which is
    /// outside the Phase 2.fdtd.7.0 scope.
    ///
    /// **Phase 2.fdtd.7.y Step C5 (Option α) note:** as of Step C5 the
    /// production pipeline ([`SubgriddedSolver::step`] /
    /// `step_with_gaussian_source_ez`) no longer calls this helper —
    /// the fine outer `E_t` plane is governed by the 1st-order Mur
    /// absorbing BC ([`Self::snapshot_fine_e_for_mur`] +
    /// [`Self::apply_mur_abc_to_fine_outer_e`]) instead. The Q3
    /// Dirichlet helper is retained as a `#[doc(hidden)]` rollback
    /// path per ADR-0038 "Consequences" and continues to be exercised
    /// by the standalone Q3 / Q4 unit tests under
    /// `crates/yee-fdtd/tests/subgrid_e_interp.rs` and
    /// `crates/yee-fdtd/tests/subgrid_h_average.rs`.
    #[doc(hidden)]
    pub fn interpolate_coarse_e_to_fine(&mut self, frac: f64) {
        let frac = frac.clamp(0.0, 1.0);
        let lo = self.lo;
        let hi = self.hi;
        let nx_c = hi.0 - lo.0;
        let ny_c = hi.1 - lo.1;
        let nz_c = hi.2 - lo.2;
        let fine_nx = 2 * nx_c;
        let fine_ny = 2 * ny_c;
        let fine_nz = 2 * nz_c;

        // ±x faces — write fine E_y and E_z at fine_i ∈ {0, fine_nx}.
        Self::interp_face_x(
            &mut self.fine,
            &self.snapshots.xmin_ey,
            &self.snapshots.xmin_ez,
            0,
            fine_ny,
            fine_nz,
            frac,
        );
        Self::interp_face_x(
            &mut self.fine,
            &self.snapshots.xmax_ey,
            &self.snapshots.xmax_ez,
            fine_nx,
            fine_ny,
            fine_nz,
            frac,
        );

        // ±y faces — write fine E_x and E_z at fine_j ∈ {0, fine_ny}.
        Self::interp_face_y(
            &mut self.fine,
            &self.snapshots.ymin_ex,
            &self.snapshots.ymin_ez,
            0,
            fine_nx,
            fine_nz,
            frac,
        );
        Self::interp_face_y(
            &mut self.fine,
            &self.snapshots.ymax_ex,
            &self.snapshots.ymax_ez,
            fine_ny,
            fine_nx,
            fine_nz,
            frac,
        );

        // ±z faces — write fine E_x and E_y at fine_k ∈ {0, fine_nz}.
        Self::interp_face_z(
            &mut self.fine,
            &self.snapshots.zmin_ex,
            &self.snapshots.zmin_ey,
            0,
            fine_nx,
            fine_ny,
            frac,
        );
        Self::interp_face_z(
            &mut self.fine,
            &self.snapshots.zmax_ex,
            &self.snapshots.zmax_ey,
            fine_nz,
            fine_nx,
            fine_ny,
            frac,
        );
    }

    /// Copy the coarse `E_t` on the six interface faces into the
    /// matching `start` or `end` snapshot buffer.
    fn copy_face_e_t(
        snap: &mut InterfaceSnapshots,
        parent: &YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        which: SnapshotKind,
    ) {
        let nx_c = hi.0 - lo.0;
        let ny_c = hi.1 - lo.1;
        let nz_c = hi.2 - lo.2;

        // ±x faces (E_y, E_z; i_c ∈ {lo.0, hi.0}).
        for (face_i_c, ey_slot, ez_slot) in [
            (lo.0, &mut snap.xmin_ey, &mut snap.xmin_ez),
            (hi.0, &mut snap.xmax_ey, &mut snap.xmax_ez),
        ] {
            let ey_buf = match which {
                SnapshotKind::Start => &mut ey_slot.0,
                SnapshotKind::End => &mut ey_slot.1,
            };
            for j_c in 0..ny_c {
                for k_c in 0..=nz_c {
                    ey_buf[(j_c, k_c)] = parent.ey[(face_i_c, lo.1 + j_c, lo.2 + k_c)];
                }
            }
            let ez_buf = match which {
                SnapshotKind::Start => &mut ez_slot.0,
                SnapshotKind::End => &mut ez_slot.1,
            };
            for j_c in 0..=ny_c {
                for k_c in 0..nz_c {
                    ez_buf[(j_c, k_c)] = parent.ez[(face_i_c, lo.1 + j_c, lo.2 + k_c)];
                }
            }
        }

        // ±y faces (E_x, E_z; j_c ∈ {lo.1, hi.1}).
        for (face_j_c, ex_slot, ez_slot) in [
            (lo.1, &mut snap.ymin_ex, &mut snap.ymin_ez),
            (hi.1, &mut snap.ymax_ex, &mut snap.ymax_ez),
        ] {
            let ex_buf = match which {
                SnapshotKind::Start => &mut ex_slot.0,
                SnapshotKind::End => &mut ex_slot.1,
            };
            for i_c in 0..nx_c {
                for k_c in 0..=nz_c {
                    ex_buf[(i_c, k_c)] = parent.ex[(lo.0 + i_c, face_j_c, lo.2 + k_c)];
                }
            }
            let ez_buf = match which {
                SnapshotKind::Start => &mut ez_slot.0,
                SnapshotKind::End => &mut ez_slot.1,
            };
            for i_c in 0..=nx_c {
                for k_c in 0..nz_c {
                    ez_buf[(i_c, k_c)] = parent.ez[(lo.0 + i_c, face_j_c, lo.2 + k_c)];
                }
            }
        }

        // ±z faces (E_x, E_y; k_c ∈ {lo.2, hi.2}).
        for (face_k_c, ex_slot, ey_slot) in [
            (lo.2, &mut snap.zmin_ex, &mut snap.zmin_ey),
            (hi.2, &mut snap.zmax_ex, &mut snap.zmax_ey),
        ] {
            let ex_buf = match which {
                SnapshotKind::Start => &mut ex_slot.0,
                SnapshotKind::End => &mut ex_slot.1,
            };
            for i_c in 0..nx_c {
                for j_c in 0..=ny_c {
                    ex_buf[(i_c, j_c)] = parent.ex[(lo.0 + i_c, lo.1 + j_c, face_k_c)];
                }
            }
            let ey_buf = match which {
                SnapshotKind::Start => &mut ey_slot.0,
                SnapshotKind::End => &mut ey_slot.1,
            };
            for i_c in 0..=nx_c {
                for j_c in 0..ny_c {
                    ey_buf[(i_c, j_c)] = parent.ey[(lo.0 + i_c, lo.1 + j_c, face_k_c)];
                }
            }
        }
    }

    /// Linear blend in time between the start and end snapshot pair at
    /// fractional offset `frac` along the coarse interval.
    #[inline]
    fn blend_time(start: f64, end: f64, frac: f64) -> f64 {
        (1.0 - frac) * start + frac * end
    }

    /// Linear-interpolation bracket for a half-integer fine-edge position
    /// against the coarse half-integer edge grid.
    ///
    /// Returns `(lo_idx, hi_idx, w_hi)` where the interpolated value is
    /// `snap[lo_idx] * (1 - w_hi) + snap[hi_idx] * w_hi`. Inside the
    /// subgrid domain `(lo_idx, hi_idx)` is the natural floor/ceil
    /// bracket and `w_hi ∈ [0, 1]`. Outside the domain the bracket is
    /// pinned to the boundary pair `(0, 1)` or `(n_coarse - 2, n_coarse - 1)`
    /// and `w_hi` is the **linear-extrapolation** weight against that
    /// pair, so a linear field in the parent reproduces exactly through
    /// the boundary cells (Chevalier 1997 §III). `n_coarse` must be ≥ 2.
    ///
    /// `t` is the position in coarse-cell units measured from the
    /// snapshot's first half-integer edge (`snap[0]` lives at `t = 0`).
    #[inline]
    fn bracket_half(t: f64, n_coarse: usize) -> (usize, usize, f64) {
        if n_coarse <= 1 {
            // Degenerate: collapse to the single sample. The caller's
            // bilerp then ignores the second tap.
            return (0, 0, 0.0);
        }
        let max_idx = n_coarse - 1;
        let lo_f = t.floor();
        if lo_f < 0.0 {
            // Linear extrapolation off the low end against (0, 1).
            return (0, 1, t);
        }
        let lo_i = lo_f as usize;
        if lo_i + 1 > max_idx {
            // Linear extrapolation off the high end against (max-1, max).
            let lo_i = max_idx - 1;
            return (lo_i, max_idx, t - (lo_i as f64));
        }
        let w_hi = t - lo_f;
        (lo_i, lo_i + 1, w_hi)
    }

    /// Linear-interpolation bracket for an integer fine-edge position
    /// (lives on a coarse node) against the coarse-node grid.
    ///
    /// Returns `(lo_idx, hi_idx, w_hi)` analogously to [`Self::bracket_half`].
    /// At even fine indices the bracket collapses (`w_hi = 0`); at odd
    /// fine indices the bracket spans one coarse cell with `w_hi = 0.5`.
    /// At the high end (fine_idx = 2·(n_coarse_nodes - 1), the last
    /// coarse node) it pins to (`max-1`, `max`) with `w_hi = 1` so the
    /// caller still gets a valid pair.
    #[inline]
    fn bracket_int(fine_idx: usize, n_coarse_nodes: usize) -> (usize, usize, f64) {
        if n_coarse_nodes <= 1 {
            return (0, 0, 0.0);
        }
        let max_idx = n_coarse_nodes - 1;
        let lo_i_natural = fine_idx / 2;
        if lo_i_natural >= max_idx {
            return (max_idx - 1, max_idx, 1.0);
        }
        let w_hi = if fine_idx.is_multiple_of(2) { 0.0 } else { 0.5 };
        (lo_i_natural, lo_i_natural + 1, w_hi)
    }

    /// Interpolate ±x face: write fine `E_y` and `E_z` at column
    /// `fine_i_face` using the snapshot pair on that face. The fine
    /// boundary edges along this face vary in `(j_f, k_f)`; spatial
    /// interpolation is in `j_f` (and `k_f`) against the coarse-edge grid
    /// stored in the snapshot, temporal blending is at `frac`.
    fn interp_face_x(
        fine: &mut YeeGrid,
        ey_snap: &(Array2<f64>, Array2<f64>),
        ez_snap: &(Array2<f64>, Array2<f64>),
        fine_i_face: usize,
        fine_ny: usize,
        fine_nz: usize,
        frac: f64,
    ) {
        let (ny_c, nz_c_p1) = ey_snap.0.dim();
        let (ny_c_p1, nz_c) = ez_snap.0.dim();

        // E_y on the face: fine_i = fine_i_face, j_f in [0, fine_ny),
        // k_f in [0, fine_nz + 1). Spatial interp: half-integer in y,
        // integer (on node) in z.
        for j_f in 0..fine_ny {
            // Half-integer fine-y position (in coarse units relative to
            // the snapshot origin): (j_f + 0.5)/2; coarse E_y edges sit
            // at coarse-y = j_c + 0.5, so the bracket target is
            // (j_f + 0.5)/2 - 0.5 = (j_f as f64 - 0.5) / 2.0.
            let t_y = ((j_f as f64) - 0.5) / 2.0;
            let (j_lo, j_hi, w_jy) = Self::bracket_half(t_y, ny_c);
            for k_f in 0..=fine_nz {
                let (k_lo, k_hi, w_kz) = Self::bracket_int(k_f, nz_c_p1);
                let s00 = ey_snap.0[(j_lo, k_lo)];
                let s01 = ey_snap.0[(j_lo, k_hi)];
                let s10 = ey_snap.0[(j_hi, k_lo)];
                let s11 = ey_snap.0[(j_hi, k_hi)];
                let e00 = ey_snap.1[(j_lo, k_lo)];
                let e01 = ey_snap.1[(j_lo, k_hi)];
                let e10 = ey_snap.1[(j_hi, k_lo)];
                let e11 = ey_snap.1[(j_hi, k_hi)];
                let start = bilerp(s00, s01, s10, s11, w_jy, w_kz);
                let end = bilerp(e00, e01, e10, e11, w_jy, w_kz);
                fine.ey[(fine_i_face, j_f, k_f)] = Self::blend_time(start, end, frac);
            }
        }

        // E_z on the face: fine_i = fine_i_face, j_f in [0, fine_ny + 1),
        // k_f in [0, fine_nz). Integer in y, half-integer in z.
        for j_f in 0..=fine_ny {
            let (j_lo, j_hi, w_jy) = Self::bracket_int(j_f, ny_c_p1);
            for k_f in 0..fine_nz {
                let t_z = ((k_f as f64) - 0.5) / 2.0;
                let (k_lo, k_hi, w_kz) = Self::bracket_half(t_z, nz_c);
                let s00 = ez_snap.0[(j_lo, k_lo)];
                let s01 = ez_snap.0[(j_lo, k_hi)];
                let s10 = ez_snap.0[(j_hi, k_lo)];
                let s11 = ez_snap.0[(j_hi, k_hi)];
                let e00 = ez_snap.1[(j_lo, k_lo)];
                let e01 = ez_snap.1[(j_lo, k_hi)];
                let e10 = ez_snap.1[(j_hi, k_lo)];
                let e11 = ez_snap.1[(j_hi, k_hi)];
                let start = bilerp(s00, s01, s10, s11, w_jy, w_kz);
                let end = bilerp(e00, e01, e10, e11, w_jy, w_kz);
                fine.ez[(fine_i_face, j_f, k_f)] = Self::blend_time(start, end, frac);
            }
        }
    }

    /// Interpolate ±y face: write fine `E_x` and `E_z` at row
    /// `fine_j_face`. Spatial interp varies in `(i_f, k_f)`.
    fn interp_face_y(
        fine: &mut YeeGrid,
        ex_snap: &(Array2<f64>, Array2<f64>),
        ez_snap: &(Array2<f64>, Array2<f64>),
        fine_j_face: usize,
        fine_nx: usize,
        fine_nz: usize,
        frac: f64,
    ) {
        let (nx_c, nz_c_p1) = ex_snap.0.dim();
        let (nx_c_p1, nz_c) = ez_snap.0.dim();

        // E_x: half-integer in x, integer in z.
        for i_f in 0..fine_nx {
            let t_x = ((i_f as f64) - 0.5) / 2.0;
            let (i_lo, i_hi, w_ix) = Self::bracket_half(t_x, nx_c);
            for k_f in 0..=fine_nz {
                let (k_lo, k_hi, w_kz) = Self::bracket_int(k_f, nz_c_p1);
                let s00 = ex_snap.0[(i_lo, k_lo)];
                let s01 = ex_snap.0[(i_lo, k_hi)];
                let s10 = ex_snap.0[(i_hi, k_lo)];
                let s11 = ex_snap.0[(i_hi, k_hi)];
                let e00 = ex_snap.1[(i_lo, k_lo)];
                let e01 = ex_snap.1[(i_lo, k_hi)];
                let e10 = ex_snap.1[(i_hi, k_lo)];
                let e11 = ex_snap.1[(i_hi, k_hi)];
                let start = bilerp(s00, s01, s10, s11, w_ix, w_kz);
                let end = bilerp(e00, e01, e10, e11, w_ix, w_kz);
                fine.ex[(i_f, fine_j_face, k_f)] = Self::blend_time(start, end, frac);
            }
        }

        // E_z: integer in x, half-integer in z.
        for i_f in 0..=fine_nx {
            let (i_lo, i_hi, w_ix) = Self::bracket_int(i_f, nx_c_p1);
            for k_f in 0..fine_nz {
                let t_z = ((k_f as f64) - 0.5) / 2.0;
                let (k_lo, k_hi, w_kz) = Self::bracket_half(t_z, nz_c);
                let s00 = ez_snap.0[(i_lo, k_lo)];
                let s01 = ez_snap.0[(i_lo, k_hi)];
                let s10 = ez_snap.0[(i_hi, k_lo)];
                let s11 = ez_snap.0[(i_hi, k_hi)];
                let e00 = ez_snap.1[(i_lo, k_lo)];
                let e01 = ez_snap.1[(i_lo, k_hi)];
                let e10 = ez_snap.1[(i_hi, k_lo)];
                let e11 = ez_snap.1[(i_hi, k_hi)];
                let start = bilerp(s00, s01, s10, s11, w_ix, w_kz);
                let end = bilerp(e00, e01, e10, e11, w_ix, w_kz);
                fine.ez[(i_f, fine_j_face, k_f)] = Self::blend_time(start, end, frac);
            }
        }
    }

    /// Interpolate ±z face: write fine `E_x` and `E_y` at plane
    /// `fine_k_face`. Spatial interp varies in `(i_f, j_f)`.
    fn interp_face_z(
        fine: &mut YeeGrid,
        ex_snap: &(Array2<f64>, Array2<f64>),
        ey_snap: &(Array2<f64>, Array2<f64>),
        fine_k_face: usize,
        fine_nx: usize,
        fine_ny: usize,
        frac: f64,
    ) {
        let (nx_c, ny_c_p1) = ex_snap.0.dim();
        let (nx_c_p1, ny_c) = ey_snap.0.dim();

        // E_x: half-integer in x, integer in y.
        for i_f in 0..fine_nx {
            let t_x = ((i_f as f64) - 0.5) / 2.0;
            let (i_lo, i_hi, w_ix) = Self::bracket_half(t_x, nx_c);
            for j_f in 0..=fine_ny {
                let (j_lo, j_hi, w_jy) = Self::bracket_int(j_f, ny_c_p1);
                let s00 = ex_snap.0[(i_lo, j_lo)];
                let s01 = ex_snap.0[(i_lo, j_hi)];
                let s10 = ex_snap.0[(i_hi, j_lo)];
                let s11 = ex_snap.0[(i_hi, j_hi)];
                let e00 = ex_snap.1[(i_lo, j_lo)];
                let e01 = ex_snap.1[(i_lo, j_hi)];
                let e10 = ex_snap.1[(i_hi, j_lo)];
                let e11 = ex_snap.1[(i_hi, j_hi)];
                let start = bilerp(s00, s01, s10, s11, w_ix, w_jy);
                let end = bilerp(e00, e01, e10, e11, w_ix, w_jy);
                fine.ex[(i_f, j_f, fine_k_face)] = Self::blend_time(start, end, frac);
            }
        }

        // E_y: integer in x, half-integer in y.
        for i_f in 0..=fine_nx {
            let (i_lo, i_hi, w_ix) = Self::bracket_int(i_f, nx_c_p1);
            for j_f in 0..fine_ny {
                let t_y = ((j_f as f64) - 0.5) / 2.0;
                let (j_lo, j_hi, w_jy) = Self::bracket_half(t_y, ny_c);
                let s00 = ey_snap.0[(i_lo, j_lo)];
                let s01 = ey_snap.0[(i_lo, j_hi)];
                let s10 = ey_snap.0[(i_hi, j_lo)];
                let s11 = ey_snap.0[(i_hi, j_hi)];
                let e00 = ey_snap.1[(i_lo, j_lo)];
                let e01 = ey_snap.1[(i_lo, j_hi)];
                let e10 = ey_snap.1[(i_hi, j_lo)];
                let e11 = ey_snap.1[(i_hi, j_hi)];
                let start = bilerp(s00, s01, s10, s11, w_ix, w_jy);
                let end = bilerp(e00, e01, e10, e11, w_ix, w_jy);
                fine.ey[(i_f, j_f, fine_k_face)] = Self::blend_time(start, end, frac);
            }
        }
    }

    // ----------------------------------------------------------------
    // Q4 — fine → coarse H_t area-average and E_t edge-average closure
    // ----------------------------------------------------------------

    /// Overwrite the parent grid's tangential `H` on the six interface
    /// faces with the area-weighted average of the four fine-grid
    /// `H_t` cells covering each coarse face.
    ///
    /// This is the closure step from Chevalier 1997 §IV: after both fine
    /// sub-steps have completed for the current coarse interval, the
    /// coarse `H_t` cells **just inside the subgrid boundary** (one coarse
    /// cell layer adjacent to each face) inherit the fine-grid solution.
    /// The next coarse `E` update outside the subgrid then sees a
    /// consistent `H_t` on the interface — which is what closes the
    /// discrete energy balance to second order in `dx_coarse`.
    ///
    /// Layer overwritten per face (in coarse-cell indices):
    /// - `−x` face: `i_c = lo.0`
    /// - `+x` face: `i_c = hi.0 − 1`
    /// - `−y` face: `j_c = lo.1`
    /// - `+y` face: `j_c = hi.1 − 1`
    /// - `−z` face: `k_c = lo.2`
    /// - `+z` face: `k_c = hi.2 − 1`
    ///
    /// On 2× refinement each coarse face covers a 2×2 tile of fine faces
    /// (uniform refinement → equal area weighting → arithmetic mean of
    /// four fine samples).
    ///
    /// **Time centering (Q4.1).** When a mid-coarse-step snapshot has
    /// been taken via [`Self::snapshot_fine_h_mid_step`] (the Q5 step
    /// driver does so between the two fine sub-steps), the values fed
    /// into the area-average are the **arithmetic mean of the snapshot
    /// and the current fine H**, i.e. `(H_f^{n+1/4} + H_f^{n+3/4}) / 2`.
    /// That recovers the time-centered fine-H value `t = n + 1/2` that
    /// the coarse slot represents. Absent a snapshot the closure falls
    /// back to a single-sample area-average of the current fine H — the
    /// pre-Q4.1 behaviour, preserved so the Q4 unit tests
    /// (which only exercise the spatial average) keep working without
    /// having to fake a snapshot.
    ///
    /// **Phase 2.fdtd.7.x B2 status — retired from the step pipeline.**
    /// [`SubgriddedSolver::step`] no longer calls this method; the
    /// fine → coarse coupling is now done by
    /// [`Self::inject_equivalent_currents_to_coarse`] (Berenger 2006
    /// equivalent-current re-radiation). This method is retained
    /// `#[doc(hidden)]` for posterity per ADR-0035 — see
    /// `docs/src/decisions/0035-berenger-huygens-subgridding.md`. The
    /// existing Q4 unit tests pin its bit-for-bit behaviour; do not
    /// modify them without a spec amendment.
    #[doc(hidden)]
    pub fn average_fine_h_to_coarse(&self, parent: &mut YeeGrid) {
        let lo = self.lo;
        let hi = self.hi;
        let fine_nx = 2 * (hi.0 - lo.0);
        let fine_ny = 2 * (hi.1 - lo.1);
        let fine_nz = 2 * (hi.2 - lo.2);

        // Q4.1: time-center fine H against the mid-step snapshot if one
        // was taken. Compute owned, time-averaged H_x/H_y/H_z arrays once
        // and read them through the spatial-average helpers. If no
        // snapshot is present, the helpers read the live fine H instead
        // (pre-Q4.1 behaviour).
        let snap = self.fine_h_snapshot.as_ref();
        let hx_view = match snap {
            Some(s) => Self::time_avg(&self.fine.hx, &s.hx),
            None => self.fine.hx.clone(),
        };
        let hy_view = match snap {
            Some(s) => Self::time_avg(&self.fine.hy, &s.hy),
            None => self.fine.hy.clone(),
        };
        let hz_view = match snap {
            Some(s) => Self::time_avg(&self.fine.hz, &s.hz),
            None => self.fine.hz.clone(),
        };

        // ±x faces — overwrite coarse H_y, H_z on the layer i_c ∈ {lo.0, hi.0 − 1}.
        Self::avg_face_x(&hy_view, &hz_view, parent, lo, hi, lo.0, 0);
        Self::avg_face_x(&hy_view, &hz_view, parent, lo, hi, hi.0 - 1, fine_nx - 2);

        // ±y faces — overwrite coarse H_x, H_z on the layer j_c ∈ {lo.1, hi.1 − 1}.
        Self::avg_face_y(&hx_view, &hz_view, parent, lo, hi, lo.1, 0);
        Self::avg_face_y(&hx_view, &hz_view, parent, lo, hi, hi.1 - 1, fine_ny - 2);

        // ±z faces — overwrite coarse H_x, H_y on the layer k_c ∈ {lo.2, hi.2 − 1}.
        Self::avg_face_z(&hx_view, &hy_view, parent, lo, hi, lo.2, 0);
        Self::avg_face_z(&hx_view, &hy_view, parent, lo, hi, hi.2 - 1, fine_nz - 2);
    }

    /// Elementwise arithmetic mean of two equal-shape arrays. Allocates
    /// a fresh owning array (Q4.1 fine-H time-centering helper).
    fn time_avg(a: &Array3<f64>, b: &Array3<f64>) -> Array3<f64> {
        let mut out = a.clone();
        out.zip_mut_with(b, |x, y| *x = 0.5 * (*x + *y));
        out
    }

    /// Overwrite the parent grid's tangential `E` on the six interface
    /// faces with the edge-averaged fine `E_t` (two fine edges per coarse
    /// edge under 2× refinement).
    ///
    /// Symmetric closure to [`Self::average_fine_h_to_coarse`] for stage 7
    /// of the spec §3 time-step pattern. Coarse `E_t` edges lie *on* the
    /// interface plane (unlike coarse `H_t` which is cell-centered in the
    /// normal direction), so the affected coarse indices are the boundary
    /// nodes `i ∈ {lo.0, hi.0}`, `j ∈ {lo.1, hi.1}`, `k ∈ {lo.2, hi.2}`
    /// for each respective face.
    ///
    /// **Phase 2.fdtd.7.x B2 status — retired from the step pipeline.**
    /// [`SubgriddedSolver::step`] no longer calls this method; the
    /// fine → coarse coupling is now done by
    /// [`Self::inject_equivalent_currents_to_coarse`] (Berenger 2006
    /// equivalent-current re-radiation). Retained `#[doc(hidden)]` for
    /// posterity per ADR-0035 — see
    /// `docs/src/decisions/0035-berenger-huygens-subgridding.md`.
    #[doc(hidden)]
    pub fn overwrite_coarse_e_from_fine(&self, parent: &mut YeeGrid) {
        let lo = self.lo;
        let hi = self.hi;
        let fine_nx = 2 * (hi.0 - lo.0);
        let fine_ny = 2 * (hi.1 - lo.1);
        let fine_nz = 2 * (hi.2 - lo.2);

        // ±x faces — overwrite coarse E_y, E_z on the planes i_c ∈ {lo.0, hi.0}.
        Self::overwrite_face_x(&self.fine, parent, lo, hi, lo.0, 0);
        Self::overwrite_face_x(&self.fine, parent, lo, hi, hi.0, fine_nx);

        // ±y faces — overwrite coarse E_x, E_z on the planes j_c ∈ {lo.1, hi.1}.
        Self::overwrite_face_y(&self.fine, parent, lo, hi, lo.1, 0);
        Self::overwrite_face_y(&self.fine, parent, lo, hi, hi.1, fine_ny);

        // ±z faces — overwrite coarse E_x, E_y on the planes k_c ∈ {lo.2, hi.2}.
        Self::overwrite_face_z(&self.fine, parent, lo, hi, lo.2, 0);
        Self::overwrite_face_z(&self.fine, parent, lo, hi, hi.2, fine_nz);
    }

    /// Area-average fine `H_y`, `H_z` onto a single coarse-cell layer
    /// `i_c_face` adjacent to a ±x interface face.
    ///
    /// `fine_i_lo` is the first fine x-index covered by the coarse layer
    /// (`0` for the −x face, `fine_nx − 2` for the +x face). `fine_hy` /
    /// `fine_hz` are the (already time-averaged, Q4.1) fine `H_y` / `H_z`
    /// arrays.
    fn avg_face_x(
        fine_hy: &Array3<f64>,
        fine_hz: &Array3<f64>,
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        i_c_face: usize,
        fine_i_lo: usize,
    ) {
        // H_y on the layer: coarse hy[(i_c_face, j_c, k_c)], j_c ∈ [lo.1, hi.1],
        // k_c ∈ [lo.2, hi.2). 4 fine H_y cells per coarse: fine_i ∈ {fine_i_lo,
        // fine_i_lo+1}, fine_j_node = 2*(j_c − lo.1), fine_k ∈ {2*(k_c − lo.2),
        // 2*(k_c − lo.2) + 1}.
        for j_c in lo.1..=hi.1 {
            let j_f = 2 * (j_c - lo.1);
            for k_c in lo.2..hi.2 {
                let k_f0 = 2 * (k_c - lo.2);
                let s = fine_hy[(fine_i_lo, j_f, k_f0)]
                    + fine_hy[(fine_i_lo + 1, j_f, k_f0)]
                    + fine_hy[(fine_i_lo, j_f, k_f0 + 1)]
                    + fine_hy[(fine_i_lo + 1, j_f, k_f0 + 1)];
                parent.hy[(i_c_face, j_c, k_c)] = 0.25 * s;
            }
        }
        // H_z on the layer: coarse hz[(i_c_face, j_c, k_c)], j_c ∈ [lo.1, hi.1),
        // k_c ∈ [lo.2, hi.2]. 4 fine H_z cells: fine_i ∈ {fine_i_lo, fine_i_lo+1},
        // fine_j ∈ {2*(j_c − lo.1), 2*(j_c − lo.1) + 1}, fine_k_node = 2*(k_c − lo.2).
        for j_c in lo.1..hi.1 {
            let j_f0 = 2 * (j_c - lo.1);
            for k_c in lo.2..=hi.2 {
                let k_f = 2 * (k_c - lo.2);
                let s = fine_hz[(fine_i_lo, j_f0, k_f)]
                    + fine_hz[(fine_i_lo + 1, j_f0, k_f)]
                    + fine_hz[(fine_i_lo, j_f0 + 1, k_f)]
                    + fine_hz[(fine_i_lo + 1, j_f0 + 1, k_f)];
                parent.hz[(i_c_face, j_c, k_c)] = 0.25 * s;
            }
        }
    }

    /// Area-average fine `H_x`, `H_z` onto a coarse layer adjacent to a
    /// ±y interface face.
    fn avg_face_y(
        fine_hx: &Array3<f64>,
        fine_hz: &Array3<f64>,
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        j_c_face: usize,
        fine_j_lo: usize,
    ) {
        // H_x: coarse hx[(i_c, j_c_face, k_c)], i_c ∈ [lo.0, hi.0], k_c ∈ [lo.2, hi.2).
        // 4 fine H_x cells: fine_i_node = 2*(i_c − lo.0), fine_j ∈ {fine_j_lo,
        // fine_j_lo+1}, fine_k ∈ {2*(k_c − lo.2), 2*(k_c − lo.2)+1}.
        for i_c in lo.0..=hi.0 {
            let i_f = 2 * (i_c - lo.0);
            for k_c in lo.2..hi.2 {
                let k_f0 = 2 * (k_c - lo.2);
                let s = fine_hx[(i_f, fine_j_lo, k_f0)]
                    + fine_hx[(i_f, fine_j_lo + 1, k_f0)]
                    + fine_hx[(i_f, fine_j_lo, k_f0 + 1)]
                    + fine_hx[(i_f, fine_j_lo + 1, k_f0 + 1)];
                parent.hx[(i_c, j_c_face, k_c)] = 0.25 * s;
            }
        }
        // H_z: coarse hz[(i_c, j_c_face, k_c)], i_c ∈ [lo.0, hi.0), k_c ∈ [lo.2, hi.2].
        // 4 fine H_z cells: fine_i ∈ {2*(i_c − lo.0), 2*(i_c − lo.0)+1}, fine_j ∈
        // {fine_j_lo, fine_j_lo+1}, fine_k_node = 2*(k_c − lo.2).
        for i_c in lo.0..hi.0 {
            let i_f0 = 2 * (i_c - lo.0);
            for k_c in lo.2..=hi.2 {
                let k_f = 2 * (k_c - lo.2);
                let s = fine_hz[(i_f0, fine_j_lo, k_f)]
                    + fine_hz[(i_f0 + 1, fine_j_lo, k_f)]
                    + fine_hz[(i_f0, fine_j_lo + 1, k_f)]
                    + fine_hz[(i_f0 + 1, fine_j_lo + 1, k_f)];
                parent.hz[(i_c, j_c_face, k_c)] = 0.25 * s;
            }
        }
    }

    /// Area-average fine `H_x`, `H_y` onto a coarse layer adjacent to a
    /// ±z interface face.
    fn avg_face_z(
        fine_hx: &Array3<f64>,
        fine_hy: &Array3<f64>,
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        k_c_face: usize,
        fine_k_lo: usize,
    ) {
        // H_x: coarse hx[(i_c, j_c, k_c_face)], i_c ∈ [lo.0, hi.0], j_c ∈ [lo.1, hi.1).
        // 4 fine H_x cells: fine_i_node = 2*(i_c − lo.0), fine_j ∈ {2*(j_c − lo.1),
        // 2*(j_c − lo.1)+1}, fine_k ∈ {fine_k_lo, fine_k_lo+1}.
        for i_c in lo.0..=hi.0 {
            let i_f = 2 * (i_c - lo.0);
            for j_c in lo.1..hi.1 {
                let j_f0 = 2 * (j_c - lo.1);
                let s = fine_hx[(i_f, j_f0, fine_k_lo)]
                    + fine_hx[(i_f, j_f0 + 1, fine_k_lo)]
                    + fine_hx[(i_f, j_f0, fine_k_lo + 1)]
                    + fine_hx[(i_f, j_f0 + 1, fine_k_lo + 1)];
                parent.hx[(i_c, j_c, k_c_face)] = 0.25 * s;
            }
        }
        // H_y: coarse hy[(i_c, j_c, k_c_face)], i_c ∈ [lo.0, hi.0), j_c ∈ [lo.1, hi.1].
        // 4 fine H_y cells: fine_i ∈ {2*(i_c − lo.0), 2*(i_c − lo.0)+1}, fine_j_node
        // = 2*(j_c − lo.1), fine_k ∈ {fine_k_lo, fine_k_lo+1}.
        for i_c in lo.0..hi.0 {
            let i_f0 = 2 * (i_c - lo.0);
            for j_c in lo.1..=hi.1 {
                let j_f = 2 * (j_c - lo.1);
                let s = fine_hy[(i_f0, j_f, fine_k_lo)]
                    + fine_hy[(i_f0 + 1, j_f, fine_k_lo)]
                    + fine_hy[(i_f0, j_f, fine_k_lo + 1)]
                    + fine_hy[(i_f0 + 1, j_f, fine_k_lo + 1)];
                parent.hy[(i_c, j_c, k_c_face)] = 0.25 * s;
            }
        }
    }

    /// Edge-average fine `E_y`, `E_z` onto the coarse `E_t` plane at
    /// `i_c_face` (a ±x interface face).
    ///
    /// `fine_i_face` is the fine x-node index that coincides with the
    /// coarse boundary plane (`0` for −x, `fine_nx` for +x).
    fn overwrite_face_x(
        fine: &YeeGrid,
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        i_c_face: usize,
        fine_i_face: usize,
    ) {
        // E_y on the face: coarse ey[(i_c_face, j_c, k_c)], j_c ∈ [lo.1, hi.1),
        // k_c ∈ [lo.2, hi.2]. 2 fine E_y edges per coarse: fine_i_node = fine_i_face,
        // fine_j ∈ {2*(j_c − lo.1), 2*(j_c − lo.1)+1}, fine_k_node = 2*(k_c − lo.2).
        for j_c in lo.1..hi.1 {
            let j_f0 = 2 * (j_c - lo.1);
            for k_c in lo.2..=hi.2 {
                let k_f = 2 * (k_c - lo.2);
                let s = fine.ey[(fine_i_face, j_f0, k_f)] + fine.ey[(fine_i_face, j_f0 + 1, k_f)];
                parent.ey[(i_c_face, j_c, k_c)] = 0.5 * s;
            }
        }
        // E_z on the face: coarse ez[(i_c_face, j_c, k_c)], j_c ∈ [lo.1, hi.1],
        // k_c ∈ [lo.2, hi.2). 2 fine E_z edges: fine_i_node = fine_i_face,
        // fine_j_node = 2*(j_c − lo.1), fine_k ∈ {2*(k_c − lo.2), 2*(k_c − lo.2)+1}.
        for j_c in lo.1..=hi.1 {
            let j_f = 2 * (j_c - lo.1);
            for k_c in lo.2..hi.2 {
                let k_f0 = 2 * (k_c - lo.2);
                let s = fine.ez[(fine_i_face, j_f, k_f0)] + fine.ez[(fine_i_face, j_f, k_f0 + 1)];
                parent.ez[(i_c_face, j_c, k_c)] = 0.5 * s;
            }
        }
    }

    /// Edge-average fine `E_x`, `E_z` onto the coarse `E_t` plane at
    /// `j_c_face` (a ±y interface face).
    fn overwrite_face_y(
        fine: &YeeGrid,
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        j_c_face: usize,
        fine_j_face: usize,
    ) {
        // E_x: coarse ex[(i_c, j_c_face, k_c)], i_c ∈ [lo.0, hi.0), k_c ∈ [lo.2, hi.2].
        for i_c in lo.0..hi.0 {
            let i_f0 = 2 * (i_c - lo.0);
            for k_c in lo.2..=hi.2 {
                let k_f = 2 * (k_c - lo.2);
                let s = fine.ex[(i_f0, fine_j_face, k_f)] + fine.ex[(i_f0 + 1, fine_j_face, k_f)];
                parent.ex[(i_c, j_c_face, k_c)] = 0.5 * s;
            }
        }
        // E_z: coarse ez[(i_c, j_c_face, k_c)], i_c ∈ [lo.0, hi.0], k_c ∈ [lo.2, hi.2).
        for i_c in lo.0..=hi.0 {
            let i_f = 2 * (i_c - lo.0);
            for k_c in lo.2..hi.2 {
                let k_f0 = 2 * (k_c - lo.2);
                let s = fine.ez[(i_f, fine_j_face, k_f0)] + fine.ez[(i_f, fine_j_face, k_f0 + 1)];
                parent.ez[(i_c, j_c_face, k_c)] = 0.5 * s;
            }
        }
    }

    /// Edge-average fine `E_x`, `E_y` onto the coarse `E_t` plane at
    /// `k_c_face` (a ±z interface face).
    fn overwrite_face_z(
        fine: &YeeGrid,
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        k_c_face: usize,
        fine_k_face: usize,
    ) {
        // E_x: coarse ex[(i_c, j_c, k_c_face)], i_c ∈ [lo.0, hi.0), j_c ∈ [lo.1, hi.1].
        for i_c in lo.0..hi.0 {
            let i_f0 = 2 * (i_c - lo.0);
            for j_c in lo.1..=hi.1 {
                let j_f = 2 * (j_c - lo.1);
                let s = fine.ex[(i_f0, j_f, fine_k_face)] + fine.ex[(i_f0 + 1, j_f, fine_k_face)];
                parent.ex[(i_c, j_c, k_c_face)] = 0.5 * s;
            }
        }
        // E_y: coarse ey[(i_c, j_c, k_c_face)], i_c ∈ [lo.0, hi.0], j_c ∈ [lo.1, hi.1).
        for i_c in lo.0..=hi.0 {
            let i_f = 2 * (i_c - lo.0);
            for j_c in lo.1..hi.1 {
                let j_f0 = 2 * (j_c - lo.1);
                let s = fine.ey[(i_f, j_f0, fine_k_face)] + fine.ey[(i_f, j_f0 + 1, fine_k_face)];
                parent.ey[(i_c, j_c, k_c_face)] = 0.5 * s;
            }
        }
    }

    /// Reject regions that touch a CPML cell on any face. The interior
    /// (non-CPML) coarse cells are `[npml, n - npml)` per axis; the
    /// subgrid `lo..hi` interval must lie inside that half-open range.
    fn check_cpml_disjoint(
        parent: &YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        npml: usize,
    ) -> Result<(), Error> {
        let dims = [parent.nx, parent.ny, parent.nz];
        let lows = [lo.0, lo.1, lo.2];
        let highs = [hi.0, hi.1, hi.2];
        for axis in 0..3 {
            let n = dims[axis];
            if n < 2 * npml {
                return Err(Error::Invalid(format!(
                    "SubgridRegion: parent axis {axis} of length {n} is too small for CPML thickness {npml}"
                )));
            }
            let inner_lo = npml;
            let inner_hi = n - npml;
            if lows[axis] < inner_lo || highs[axis] > inner_hi {
                return Err(Error::Invalid(format!(
                    "SubgridRegion: region [{}, {}) on axis {axis} overlaps CPML thickness (inner range [{}, {}), npml={})",
                    lows[axis], highs[axis], inner_lo, inner_hi, npml
                )));
            }
        }
        Ok(())
    }

    /// Reject regions whose `lo..hi` interval *crosses* (straddles) any
    /// face of the supplied TF/SF box. A region wholly inside or wholly
    /// outside the box is permitted; only a face-crossing nest breaks
    /// the TF/SF reciprocity argument (spec §6).
    fn check_tfsf_disjoint(
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        tlo: (usize, usize, usize),
        thi: (usize, usize, usize),
    ) -> Result<(), Error> {
        let lows = [lo.0, lo.1, lo.2];
        let highs = [hi.0, hi.1, hi.2];
        let tlows = [tlo.0, tlo.1, tlo.2];
        let thighs = [thi.0, thi.1, thi.2];
        // Per-axis, the region "crosses" a TF/SF face if its interval
        // strictly straddles either tlo or thi.
        for axis in 0..3 {
            let l = lows[axis];
            let h = highs[axis];
            let tl = tlows[axis];
            let th = thighs[axis];
            let crosses_lo_face = l < tl && h > tl;
            let crosses_hi_face = l < th && h > th;
            if crosses_lo_face || crosses_hi_face {
                return Err(Error::Invalid(format!(
                    "SubgridRegion: region [{l}, {h}) on axis {axis} crosses TF/SF box face (tfsf [{tl}, {th}))"
                )));
            }
        }
        Ok(())
    }

    // ----------------------------------------------------------------
    // Phase 2.fdtd.7.x B1 — Berenger Huygens-surface skeleton
    // ----------------------------------------------------------------

    /// Inject equivalent surface currents `J = +n̂ × H_tot` and
    /// `M = −n̂ × E_tot` from the fine subdomain onto the six Huygens
    /// faces of the parent coarse grid (Berenger 2006, *IEEE T-AP*
    /// 54(12), pp. 3797–3804, §III).
    ///
    /// One-directional, post-coarse-update closure: the fine grid's
    /// storage is read only; the coarse grid's `E` and `H` arrays are
    /// mutated in-place at the Huygens surface only.
    ///
    /// ## Time centering
    ///
    /// - The `J = +n̂ × H_tot` source is sampled at fine wall-clock time
    ///   `t = n + 1/2` by averaging the Q4.1 mid-step snapshot
    ///   ([`Self::snapshot_fine_h_mid_step`], taken at `t = n + 1/4`)
    ///   with the current fine `H` (at `t = n + 3/4` after sub-step 2);
    ///   absent a snapshot the closure falls back to the single-sample
    ///   current fine `H` (the bit-exact-fresh-construction behaviour
    ///   that keeps the B1 no-op pin valid on a zero-field fine grid).
    /// - The `M = -n̂ × E_tot` source is sampled at `t = n + 1` from the
    ///   post-`update_fine_e` outer-layer fine `E_t` (sub-step 2's
    ///   final `update_fine_e`).
    ///
    /// ## TF/SF convention (sign discipline)
    ///
    /// Following spec §3.1: the fine box's coarse-cell footprint is
    /// the **scattered-field** (SF) side; everything strictly outside
    /// is **total-field** (TF). The six Huygens faces sit at the SF/TF
    /// boundary. The per-face corrections are the equivalent-current
    /// counterpart of the existing TF/SF plane-wave injection
    /// (`crates/yee-fdtd/src/sources.rs`) with the SF and TF roles
    /// inverted — see the per-face derivation in
    /// `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`
    /// §3 and ADR-0035.
    ///
    /// Concretely the additions per face are
    ///
    /// ```text
    /// E_t_coarse += s · (-dt / (ε₀·ε_r·dx_n)) · H_t_fine_just_inside    (J term)
    /// H_t_coarse += s · (-dt / (μ₀·μ_r·dx_n)) · E_t_fine_on_surface     (M term)
    /// ```
    ///
    /// where `dx_n` is the coarse cell size normal to the face and the
    /// per-component sign `s` comes from the `n̂ × ·` cross product
    /// (see implementation per-face sign table in the source body).
    /// The global `-1` factor on `coeff_e` / `coeff_h` encodes the
    /// inverted-TF/SF convention: standard TF/SF (TF inside) subtracts
    /// an `H_inc` / `E_inc` excess from an SF stencil reading a TF
    /// neighbour; Berenger inverts the role of the SF and TF sides,
    /// which flips the correction sign relative to the natural "add
    /// J·dt/ε" form.
    ///
    /// The cuboid-edge cells (where two adjacent faces meet) are
    /// assigned to a single owning face under the lower-numbered-axis
    /// wins rule ([`assign_edge_to_face`]). Concretely the tangential
    /// index range on each face is half-open on the higher-axis edges:
    /// `±x` faces use the full `[lo.1, hi.1)` × `[lo.2, hi.2)` extent
    /// (X axis is lowest, wins every shared edge); `±y` faces use
    /// `[lo.0, hi.0)` for the X tangent (no contribution at the X-edge
    /// because X already owns it) but the full `[lo.2, hi.2)` for the Z
    /// tangent; `±z` faces use `[lo.0, hi.0)` × `[lo.1, hi.1)`
    /// (everything is owned by X or Y). This convention is what the B1
    /// `cuboid_edge_owned_by_one_face_only` unit test (sibling to
    /// `assign_edge_to_face`) verifies independently.
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`.
    /// ADR: `docs/src/decisions/0035-berenger-huygens-subgridding.md`.
    pub fn inject_equivalent_currents_to_coarse(&self, parent: &mut YeeGrid) {
        // Legacy monolithic entry point — preserved without B2.2 coarse-
        // ghost subtraction so the B1 skeleton no-op test
        // (`berenger_skeleton::inject_equivalent_currents_to_coarse_is_currently_noop`)
        // continues to pass bit-exactly: when the fine grid is zero,
        // the un-ghosted `J = sign·coef·H_fine_avg = 0` and `M =
        // sign·coef·E_fine_avg = 0`, leaving the parent grid untouched
        // even if the parent carries seeded non-zero fields. The ghost-
        // subtracted (B2.2) form is reserved for the split-injection
        // pipeline used by [`SubgriddedSolver::step`].
        self.inject_currents_inner(parent, true, true, false);
    }

    /// Inject the electric equivalent current
    /// `J = +n̂ × (H_TF_fine − H_SF_coarse_ghost)` onto the coarse `E_t`
    /// arrays at the Huygens surface.
    ///
    /// Phase 2.fdtd.7.x B2.2 coarse-ghost-subtraction form (Berenger
    /// 2006 §III canonical equivalent-source convention). LLLLLLL's
    /// B2.1 diagnosis traced the residual B2 divergence (peak coarse
    /// `|E_z|` ≈ 1e30 at 500 steps even with the split-injection
    /// reorder) to the J source representing `+n̂ × H_fine_inside`
    /// **without** subtracting the coarse `H` slot the coarse natural
    /// curl already reads for its leapfrog update. Without that
    /// subtraction the J term carries gain > 1 around the loop
    /// `fine E → coarse E → coarse H → fine boundary E → fine H →
    /// larger J`. Subtracting the same-spatial-location coarse-H ghost
    /// converts the injection from full TF to the canonical
    /// scattered-correction form (Berenger 2006 §III), which closes
    /// the equivalence principle at the discrete level.
    ///
    /// `H_tot` is sampled at `t = n + 1/2` via the Q4.1 mid-step snapshot
    /// plus post-sub-step-2 average (same time-centering as the legacy
    /// monolithic [`Self::inject_equivalent_currents_to_coarse`]).
    pub fn inject_j_to_coarse_e(&self, parent: &mut YeeGrid) {
        self.inject_currents_inner(parent, true, false, true);
    }

    /// Inject the magnetic equivalent current
    /// `M = -n̂ × (E_post − E_pre)` (Phase 2.fdtd.7.y Option β
    /// compensating-source form) onto the coarse `H_t` arrays at the
    /// Huygens surface. The fine-`E` pair is sourced from the
    /// [`Self::snapshot_fine_e_pre_update`] and
    /// [`Self::snapshot_fine_e_post_update`] snapshots captured by
    /// [`SubgriddedSolver::step`] either side of sub-step 2's
    /// [`Self::update_fine_e`].
    ///
    /// **Why the compensating form, not the canonical `-n̂ × E_fine`
    /// or `-n̂ × (E_fine − E_coarse_ghost)`.** Berenger 2006 §III gives
    /// the canonical equivalent magnetic current
    /// `M = -n̂ × (E_TF_fine − E_SF_coarse_ghost)`. On this pipeline
    /// the "coarse-E ghost" is `parent.ey[(i_c_surface, …)]` /
    /// `parent.ez[(…)]` at the surface plane, i.e. the same coarse-E
    /// slot Q3 [`Self::interpolate_coarse_e_to_fine`] used to write
    /// the fine outer-layer `E_t` Dirichlet value just before
    /// `update_fine_e`. By Q3 construction `fine_E_surface ≈
    /// coarse_E_surface` to interpolation order, so the canonical
    /// `M_canonical ≈ 0` and is dominated by round-off noise; Track
    /// OOOOOOO confirmed this empirically (100-step canary peak
    /// `|E_z|_fine` regressing from 2.75 V/m with J-only ghost to
    /// ≈ 1.0e3 V/m once M ghost subtraction was enabled).
    ///
    /// The compensating form `M = -n̂ × (E_post − E_pre)` captures
    /// the **per-fine-sub-step Maxwell-evolved increment** to the
    /// outer-layer fine `E_t` — the part of `E_fine` that escapes the
    /// Q3 Dirichlet tie precisely because `update_fine_e` runs after
    /// Q3 writes. The Dirichlet-tied part nullifies inside the
    /// difference. The B2.1 J-side coarse-ghost subtraction
    /// ([`Self::inject_j_to_coarse_e`]) is **kept unchanged** — it
    /// remains the load-bearing improvement of Phase 2.fdtd.7.x and
    /// the differencing happens on a genuinely non-Q3-coupled `H`
    /// pair there.
    ///
    /// **Do not also enable coarse-ghost subtraction on M.** The
    /// `E_post − E_pre` sample already performs the differencing;
    /// re-subtracting a coarse `E` ghost would double-count and
    /// reintroduce OOOOOOO's failure mode. The per-face helpers carry
    /// a `debug_assert!(!(do_ghost && use_compensating_source))`
    /// guard.
    ///
    /// On the **first** coarse step after construction (or after a
    /// region is freshly attached), one or both snapshots are `None`
    /// and this helper is a no-op. The fine grid has zero fields at
    /// `t = 0`, so the correct compensating `M` is identically zero
    /// anyway; no special handling is required.
    ///
    /// Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3 (Option β).
    /// ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.
    pub fn inject_m_to_coarse_h(&self, parent: &mut YeeGrid) {
        let (Some(pre), Some(post)) = (
            self.fine_e_pre_snapshot.as_ref(),
            self.fine_e_post_snapshot.as_ref(),
        ) else {
            // First coarse step (or test path that never wired the
            // Step C1 snapshots): the compensating source is
            // identically zero. Return without touching the parent
            // grid to preserve the B1 no-op pin behaviour.
            return;
        };
        // Compensating-source path: pass both pre and post arrays
        // through; the per-face helpers read `E_post − E_pre`. The
        // `do_ghost = false` on the M path is enforced by the
        // `use_compensating_source = true` flag's debug-assert guard
        // inside `inject_*_face`; we wire it explicitly here for
        // grep-ability.
        self.inject_currents_inner_with_e_pair(
            parent, false, true, &pre.ex, &pre.ey, &pre.ez, &post.ex, &post.ey, &post.ez, false,
            true,
        );
    }

    /// Shared body for the monolithic and split-injection closures.
    /// Computes the time-centered fine `H` average (Q4.1) and dispatches
    /// to the per-face helpers with the requested `do_j` / `do_m` flags.
    /// When `do_m` is `true`, M is sourced from the current fine `E`
    /// (legacy monolithic behaviour preserved for the B1 no-op test).
    ///
    /// `do_ghost` enables the Phase 2.fdtd.7.x B2.2 coarse-ghost-
    /// subtraction form on the J term: when `true`, `J = +n̂ × (H_fine
    /// − H_coarse_ghost)` per Berenger 2006 §III; when `false`, `J =
    /// +n̂ × H_fine` (legacy un-ghosted form, kept for the monolithic
    /// entry point so the B1 skeleton no-op test continues to pass).
    fn inject_currents_inner(&self, parent: &mut YeeGrid, do_j: bool, do_m: bool, do_ghost: bool) {
        self.inject_currents_inner_with_e(
            parent,
            do_j,
            do_m,
            &self.fine.ex,
            &self.fine.ey,
            &self.fine.ez,
            do_ghost,
        );
    }

    /// Phase 2.fdtd.7.y Option β compensating-source dispatcher.
    /// Wraps [`Self::inject_currents_inner_with_e_pair`] for the
    /// no-pre-array case (legacy `M = -n̂ × E_fine`) by routing the
    /// post array as the sole M source — kept as a thin alias so the
    /// J-only / monolithic call paths above do not have to plumb the
    /// extra pre arrays. The compensating path (`E_post − E_pre`) is
    /// reached only through [`Self::inject_m_to_coarse_h`] →
    /// `inject_currents_inner_with_e_pair` with
    /// `use_compensating_source = true`.
    #[allow(clippy::too_many_arguments)]
    fn inject_currents_inner_with_e(
        &self,
        parent: &mut YeeGrid,
        do_j: bool,
        do_m: bool,
        fine_ex: &Array3<f64>,
        fine_ey: &Array3<f64>,
        fine_ez: &Array3<f64>,
        do_ghost: bool,
    ) {
        // Route through the pair helper with the same array for both
        // pre and post and `use_compensating_source = false`; the
        // per-face helpers then read the single array directly and
        // ignore the pre slot.
        self.inject_currents_inner_with_e_pair(
            parent, do_j, do_m, fine_ex, fine_ey, fine_ez, fine_ex, fine_ey, fine_ez, do_ghost,
            false,
        );
    }

    /// Phase 2.fdtd.7.y Option β-aware injection helper. Like
    /// [`Self::inject_currents_inner_with_e`] but takes both a `_pre`
    /// and a `_post` fine-E array triple plus a
    /// `use_compensating_source` flag that selects between two M-source
    /// sampling regimes:
    ///
    /// - `use_compensating_source = false` (legacy paths) — the M
    ///   source reads `fine_e*_post` directly as `M = -n̂ × E_fine`.
    ///   The `_pre` arrays are passed through but unread on the M side.
    /// - `use_compensating_source = true` (Phase 2.fdtd.7.y Option β) —
    ///   the M source reads `M = -n̂ × (E_post − E_pre)` per ADR-0038.
    ///   `do_ghost` MUST be `false` on the M side in this regime
    ///   (otherwise the per-face debug-assert fires); see
    ///   [`Self::inject_m_to_coarse_h`]'s docstring for the
    ///   double-counting rationale.
    #[allow(clippy::too_many_arguments)]
    fn inject_currents_inner_with_e_pair(
        &self,
        parent: &mut YeeGrid,
        do_j: bool,
        do_m: bool,
        fine_ex_pre: &Array3<f64>,
        fine_ey_pre: &Array3<f64>,
        fine_ez_pre: &Array3<f64>,
        fine_ex: &Array3<f64>,
        fine_ey: &Array3<f64>,
        fine_ez: &Array3<f64>,
        do_ghost: bool,
        use_compensating_source: bool,
    ) {
        let lo = self.lo;
        let hi = self.hi;
        let fine_nx = 2 * (hi.0 - lo.0);
        let fine_ny = 2 * (hi.1 - lo.1);
        let fine_nz = 2 * (hi.2 - lo.2);

        // Time-center fine H against the mid-step snapshot if one was
        // taken: `(H_f^{n+1/4} + H_f^{n+3/4}) / 2 = H_f^{n+1/2}`. Absent
        // a snapshot, fall through to the current fine H. The B1 no-op
        // pin relies on the fresh-construction path: a freshly allocated
        // SubgridRegion has zero fields everywhere and no snapshot, so
        // J = +n̂ × 0 = 0 and M = -n̂ × 0 = 0 — the injection is the
        // identity on the parent grid as the B1 test requires.
        let snap = self.fine_h_snapshot.as_ref();
        let hx_t = match snap {
            Some(s) => Self::time_avg(&self.fine.hx, &s.hx),
            None => self.fine.hx.clone(),
        };
        let hy_t = match snap {
            Some(s) => Self::time_avg(&self.fine.hy, &s.hy),
            None => self.fine.hy.clone(),
        };
        let hz_t = match snap {
            Some(s) => Self::time_avg(&self.fine.hz, &s.hz),
            None => self.fine.hz.clone(),
        };

        // Coefficient signs: the inverted TF/SF convention (SF inside the
        // fine box footprint, TF outside) flips both the J → E and the
        // M → H contributions relative to the natural "add J·dt/ε" /
        // "add M·dt/μ" forms — the standard TF/SF accumulation in
        // `sources.rs` subtracts an excess `H_inc` / `E_inc` baked into
        // the SF stencil reading a TF neighbour, and the Berenger
        // counterpart flips the role of SF and TF, which net-flips the
        // sign once more. Both effects compose to the global `-1`
        // pre-factor below; the per-face sign tables in
        // [`Self::inject_x_face`] / `inject_y_face` / `inject_z_face`
        // multiply by the outward-normal sign on top.
        let coeff_e = -parent.dt / (EPS0 * parent.eps_r);
        let coeff_h = -parent.dt / (MU0 * parent.mu_r);

        // ±x faces (normal = ±x̂). `fine_i_h` is the fine cell-centered
        // x-index of the fine H_y / H_z layer adjacent to the Huygens
        // surface from the SF (fine-interior) side. For the +x face
        // that's `fine_i = fine_nx - 1`; for the −x face `fine_i = 0`.
        // `fine_i_e` is the fine x-node index of the surface plane
        // (where fine E_y / E_z live).
        Self::inject_x_face(
            parent,
            lo,
            hi,
            fine_ey_pre,
            fine_ez_pre,
            fine_ey,
            fine_ez,
            &hy_t,
            &hz_t,
            1.0,
            hi.0,
            hi.0 - 1,
            fine_nx - 1,
            fine_nx,
            coeff_e,
            coeff_h,
            do_j,
            do_m,
            do_ghost,
            use_compensating_source,
        );
        Self::inject_x_face(
            parent,
            lo,
            hi,
            fine_ey_pre,
            fine_ez_pre,
            fine_ey,
            fine_ez,
            &hy_t,
            &hz_t,
            -1.0,
            lo.0,
            lo.0,
            0,
            0,
            coeff_e,
            coeff_h,
            do_j,
            do_m,
            do_ghost,
            use_compensating_source,
        );

        // ±y faces.
        Self::inject_y_face(
            parent,
            lo,
            hi,
            fine_ex_pre,
            fine_ez_pre,
            fine_ex,
            fine_ez,
            &hx_t,
            &hz_t,
            1.0,
            hi.1,
            hi.1 - 1,
            fine_ny - 1,
            fine_ny,
            coeff_e,
            coeff_h,
            do_j,
            do_m,
            do_ghost,
            use_compensating_source,
        );
        Self::inject_y_face(
            parent,
            lo,
            hi,
            fine_ex_pre,
            fine_ez_pre,
            fine_ex,
            fine_ez,
            &hx_t,
            &hz_t,
            -1.0,
            lo.1,
            lo.1,
            0,
            0,
            coeff_e,
            coeff_h,
            do_j,
            do_m,
            do_ghost,
            use_compensating_source,
        );

        // ±z faces.
        Self::inject_z_face(
            parent,
            lo,
            hi,
            fine_ex_pre,
            fine_ey_pre,
            fine_ex,
            fine_ey,
            &hx_t,
            &hy_t,
            1.0,
            hi.2,
            hi.2 - 1,
            fine_nz - 1,
            fine_nz,
            coeff_e,
            coeff_h,
            do_j,
            do_m,
            do_ghost,
            use_compensating_source,
        );
        Self::inject_z_face(
            parent,
            lo,
            hi,
            fine_ex_pre,
            fine_ey_pre,
            fine_ex,
            fine_ey,
            &hx_t,
            &hy_t,
            -1.0,
            lo.2,
            lo.2,
            0,
            0,
            coeff_e,
            coeff_h,
            do_j,
            do_m,
            do_ghost,
            use_compensating_source,
        );
    }

    /// Inject the equivalent-current corrections on one `±x` Huygens
    /// face. `sign = +1.0` for the `+x` face (outward `+x̂`), `-1.0`
    /// for the `-x` face. `i_c_surface` is the coarse `i` index of
    /// the surface plane (where coarse `E_y`, `E_z` live); `i_c_inside_h`
    /// is the coarse `i` index of the coarse `H_y`, `H_z` layer just
    /// **inside** the fine box (SF storage). `fine_i_h` is the fine-x
    /// cell-centered index of the fine `H_y`, `H_z` layer adjacent to
    /// the surface (SF side); `fine_i_e` is the fine-x node index of
    /// the surface plane (where fine `E_y`, `E_z` live).
    ///
    /// `do_j` and `do_m` independently enable the J → coarse E and
    /// M → coarse H paths. The Phase 2.fdtd.7.x B2.1 split-injection
    /// closure ([`Self::inject_j_to_coarse_e`] /
    /// [`Self::inject_m_to_coarse_h`]) calls this helper twice per
    /// face per coarse step — once with `do_j=true, do_m=false` before
    /// `update_e_only`, and once with `do_j=false, do_m=true` at the
    /// top of the next coarse step before `update_h_only` — so each
    /// source enters its leapfrog update before the update runs.
    ///
    /// `fine_e*_for_m` supplies the fine `E_x`, `E_y`, `E_z` arrays
    /// the `M = -n̂ × E` source samples. Behaviour depends on
    /// `use_compensating_source`:
    ///
    /// - `use_compensating_source = false` — legacy paths. `fine_e*_post`
    ///   is read directly (`M = -n̂ × E_fine`); the `_pre` slot is
    ///   ignored on the M side. Used by the monolithic
    ///   [`Self::inject_equivalent_currents_to_coarse`] entry point.
    /// - `use_compensating_source = true` — Phase 2.fdtd.7.y Option β
    ///   compensating-source form. The per-cell M source reads
    ///   `(E_post − E_pre)` cell-by-cell from the two snapshot arrays
    ///   captured by [`SubgridRegion::snapshot_fine_e_pre_update`] and
    ///   [`SubgridRegion::snapshot_fine_e_post_update`]; this isolates
    ///   the Maxwell-evolved fine-E increment per sub-step and
    ///   nullifies the Q3-Dirichlet-tied part that would otherwise
    ///   cancel against the coarse-E ghost. **`do_ghost = false` is
    ///   required in this regime** (see
    ///   [`SubgridRegion::inject_m_to_coarse_h`]'s docstring for the
    ///   double-counting rationale); a `debug_assert!` below catches
    ///   the misuse.
    ///
    /// The four signs per face come from `J = sign · x̂ × H` and
    /// `M = -sign · x̂ × E`:
    ///
    /// ```text
    /// J_y = -sign · H_z        ⇒ E_y += +sign · (dt/(ε·dx)) · H_z_fine
    /// J_z = +sign · H_y        ⇒ E_z += -sign · (dt/(ε·dx)) · H_y_fine
    /// M_y = +sign · E_z        ⇒ H_y += -sign · (dt/(μ·dx)) · E_z_fine
    /// M_z = -sign · E_y        ⇒ H_z += +sign · (dt/(μ·dx)) · E_y_fine
    /// ```
    ///
    /// (Sign of the E-equation contribution: `ε ∂E/∂t = ∇×H - J_s`,
    /// surface delta integrated over a coarse cell gives `J_s/dx_n`
    /// as the effective volumetric current, hence `E_t += -dt/ε ·
    /// J_t/dx`.)
    ///
    /// Spatial averaging: each coarse face cell `(i_c_*, j_c, k_c)`
    /// receives the arithmetic mean of the two fine `H_y` (or `H_z`)
    /// cells covering it along the tangential half-cell — analogous
    /// to [`Self::avg_face_x`]. The fine `E_y`, `E_z` are averaged
    /// over two fine edges per coarse edge — analogous to
    /// [`Self::overwrite_face_x`].
    #[allow(clippy::too_many_arguments)]
    fn inject_x_face(
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        fine_ey_pre_for_m: &Array3<f64>,
        fine_ez_pre_for_m: &Array3<f64>,
        fine_ey_for_m: &Array3<f64>,
        fine_ez_for_m: &Array3<f64>,
        fine_hy_t: &Array3<f64>,
        fine_hz_t: &Array3<f64>,
        sign: f64,
        i_c_surface: usize,
        i_c_inside_h: usize,
        fine_i_h: usize,
        fine_i_e: usize,
        coeff_e: f64,
        coeff_h: f64,
        do_j: bool,
        do_m: bool,
        do_ghost: bool,
        use_compensating_source: bool,
    ) {
        // ADR-0038: M-side coarse-ghost subtraction is incompatible
        // with the compensating-source sampling, because the
        // `E_post − E_pre` difference already does the differencing
        // that ghost subtraction would do again.
        debug_assert!(
            !(do_ghost && use_compensating_source),
            "inject_x_face: do_ghost and use_compensating_source are mutually exclusive on M"
        );
        let dx = parent.dx;
        let ce = coeff_e / dx;
        let ch = coeff_h / dx;

        if do_j {
            // -------- J term: write to E_y, E_z on the surface plane. --------
            //
            // Berenger 2006 §III equivalent-source form (B2.2,
            // `do_ghost = true`):
            //   J = +n̂ × (H_TF_fine − H_SF_coarse_ghost)
            // The "coarse ghost" is the coarse `H` slot at the same
            // spatial location that the coarse natural-curl stencil for
            // the target `E` cell already reads from inside the fine box
            // footprint (SF side). For `ey[(i_c_surface, j, k)]` the
            // natural curl reads `hz[(i_c_inside_h, j, k)]` — that IS the
            // ghost. Subtracting it converts the injection from full TF
            // ("J = n̂ × H_fine") to the canonical SF-subtracted form,
            // closing Berenger's discrete equivalence principle.
            //
            // `do_ghost = false` reverts to the legacy un-ghosted form
            // `J = +n̂ × H_TF_fine` so the monolithic
            // [`Self::inject_equivalent_currents_to_coarse`] entry point
            // continues to pass the B1 skeleton no-op test (which seeds
            // non-zero coarse fields and expects fine = 0 to yield a
            // bit-exact identity injection).

            // E_y on the surface: coarse ey[(i_c_surface, j_c, k_c)],
            // j_c ∈ [lo.1, hi.1), k_c ∈ [lo.2, hi.2].
            for j_c in lo.1..hi.1 {
                let j_f0 = 2 * (j_c - lo.1);
                for k_c in lo.2..=hi.2 {
                    let k_f = 2 * (k_c - lo.2);
                    let hz_avg = 0.5
                        * (fine_hz_t[(fine_i_h, j_f0, k_f)] + fine_hz_t[(fine_i_h, j_f0 + 1, k_f)]);
                    let hz_eff = if do_ghost {
                        hz_avg - parent.hz[(i_c_inside_h, j_c, k_c)]
                    } else {
                        hz_avg
                    };
                    parent.ey[(i_c_surface, j_c, k_c)] += sign * ce * hz_eff;
                }
            }

            // E_z on the surface: coarse ez[(i_c_surface, j_c, k_c)],
            // j_c ∈ [lo.1, hi.1], k_c ∈ [lo.2, hi.2).
            for j_c in lo.1..=hi.1 {
                let j_f = 2 * (j_c - lo.1);
                for k_c in lo.2..hi.2 {
                    let k_f0 = 2 * (k_c - lo.2);
                    let hy_avg = 0.5
                        * (fine_hy_t[(fine_i_h, j_f, k_f0)] + fine_hy_t[(fine_i_h, j_f, k_f0 + 1)]);
                    let hy_eff = if do_ghost {
                        hy_avg - parent.hy[(i_c_inside_h, j_c, k_c)]
                    } else {
                        hy_avg
                    };
                    parent.ez[(i_c_surface, j_c, k_c)] -= sign * ce * hy_eff;
                }
            }
        }

        if do_m {
            // -------- M term: write to H_y, H_z on the layer just inside. --------
            //
            // Phase 2.fdtd.7.y Option β (`use_compensating_source =
            // true`): `M = -n̂ × (E_post − E_pre)`. The per-cell fine-E
            // sample is `E_post − E_pre` cell-by-cell. The
            // Dirichlet-tied part of `E_post` cancels against `E_pre`;
            // the residual is the per-fine-sub-step Maxwell-evolved
            // increment to fine `E_t` on the surface plane, which is
            // the part Berenger's canonical M source needs.
            //
            // Legacy un-ghosted path (`use_compensating_source =
            // false`): reads only `E_post` (`M = -n̂ × E_fine`). Kept
            // for the monolithic
            // [`SubgridRegion::inject_equivalent_currents_to_coarse`]
            // entry point so the B1 skeleton no-op test continues to
            // pass.
            //
            // Coarse-E ghost subtraction on M is **not** supported in
            // either regime; see ADR-0038 and the docstring on
            // [`SubgridRegion::inject_m_to_coarse_h`] for the
            // empirical rationale.

            // H_y on the layer: coarse hy[(i_c_inside_h, j_c, k_c)],
            // j_c ∈ [lo.1, hi.1], k_c ∈ [lo.2, hi.2).
            // Source: fine E_z on the surface (fine_i = fine_i_e).
            for j_c in lo.1..=hi.1 {
                let j_f = 2 * (j_c - lo.1);
                for k_c in lo.2..hi.2 {
                    let k_f0 = 2 * (k_c - lo.2);
                    let ez0 = fine_ez_for_m[(fine_i_e, j_f, k_f0)];
                    let ez1 = fine_ez_for_m[(fine_i_e, j_f, k_f0 + 1)];
                    let ez_eff0 = if use_compensating_source {
                        ez0 - fine_ez_pre_for_m[(fine_i_e, j_f, k_f0)]
                    } else {
                        ez0
                    };
                    let ez_eff1 = if use_compensating_source {
                        ez1 - fine_ez_pre_for_m[(fine_i_e, j_f, k_f0 + 1)]
                    } else {
                        ez1
                    };
                    let ez_avg = 0.5 * (ez_eff0 + ez_eff1);
                    parent.hy[(i_c_inside_h, j_c, k_c)] -= sign * ch * ez_avg;
                }
            }

            // H_z on the layer: coarse hz[(i_c_inside_h, j_c, k_c)],
            // j_c ∈ [lo.1, hi.1), k_c ∈ [lo.2, hi.2].
            // Source: fine E_y on the surface (fine_i = fine_i_e).
            for j_c in lo.1..hi.1 {
                let j_f0 = 2 * (j_c - lo.1);
                for k_c in lo.2..=hi.2 {
                    let k_f = 2 * (k_c - lo.2);
                    let ey0 = fine_ey_for_m[(fine_i_e, j_f0, k_f)];
                    let ey1 = fine_ey_for_m[(fine_i_e, j_f0 + 1, k_f)];
                    let ey_eff0 = if use_compensating_source {
                        ey0 - fine_ey_pre_for_m[(fine_i_e, j_f0, k_f)]
                    } else {
                        ey0
                    };
                    let ey_eff1 = if use_compensating_source {
                        ey1 - fine_ey_pre_for_m[(fine_i_e, j_f0 + 1, k_f)]
                    } else {
                        ey1
                    };
                    let ey_avg = 0.5 * (ey_eff0 + ey_eff1);
                    parent.hz[(i_c_inside_h, j_c, k_c)] += sign * ch * ey_avg;
                }
            }
        }
    }

    /// Inject corrections on one `±y` face. Mirror of
    /// [`Self::inject_x_face`] with cyclic axis permutation
    /// `(x, y, z) → (y, z, x)`:
    ///
    /// ```text
    /// J_z = -sign · H_x        ⇒ E_z += +sign · (dt/(ε·dy)) · H_x_fine
    /// J_x = +sign · H_z        ⇒ E_x += -sign · (dt/(ε·dy)) · H_z_fine
    /// M_z = +sign · E_x        ⇒ H_z += -sign · (dt/(μ·dy)) · E_x_fine
    /// M_x = -sign · E_z        ⇒ H_x += +sign · (dt/(μ·dy)) · E_z_fine
    /// ```
    ///
    /// (Same derivation: `J = sign · ŷ × H = sign · (Hz, 0, -Hx)`,
    /// `M = -sign · ŷ × E = -sign · (Ez, 0, -Ex)`.)
    #[allow(clippy::too_many_arguments)]
    fn inject_y_face(
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        fine_ex_pre_for_m: &Array3<f64>,
        fine_ez_pre_for_m: &Array3<f64>,
        fine_ex_for_m: &Array3<f64>,
        fine_ez_for_m: &Array3<f64>,
        fine_hx_t: &Array3<f64>,
        fine_hz_t: &Array3<f64>,
        sign: f64,
        j_c_surface: usize,
        j_c_inside_h: usize,
        fine_j_h: usize,
        fine_j_e: usize,
        coeff_e: f64,
        coeff_h: f64,
        do_j: bool,
        do_m: bool,
        do_ghost: bool,
        use_compensating_source: bool,
    ) {
        // ADR-0038: M-side coarse-ghost subtraction is incompatible
        // with compensating-source sampling. See `inject_x_face` for
        // the rationale.
        debug_assert!(
            !(do_ghost && use_compensating_source),
            "inject_y_face: do_ghost and use_compensating_source are mutually exclusive on M"
        );
        let dy = parent.dy;
        let ce = coeff_e / dy;
        let ch = coeff_h / dy;

        if do_j {
            // -------- J term: write to E_x, E_z on the surface plane. --------
            //
            // Berenger 2006 §III equivalent-source form (B2.2,
            // `do_ghost = true`):
            //   J = +n̂ × (H_TF_fine − H_SF_coarse_ghost)
            // For target `ez[(i, j_c_surface, k)]` the coarse natural
            // curl reads `hx[(i, j_c_inside_h, k)]` from inside the
            // footprint — that's the ghost. See [`Self::inject_x_face`]
            // for the rationale on `do_ghost`.

            // E_z on the surface: coarse ez[(i_c, j_c_surface, k_c)],
            // i_c ∈ [lo.0, hi.0], k_c ∈ [lo.2, hi.2). (Cuboid edges at
            // i = lo.0 and i = hi.0 are owned by ±x faces — but the
            // surface E_z[(lo.0, j_c_surface, k_c)] / E_z[(hi.0, ...)] is
            // a different array cell than the ±x face's targets, so no
            // double-count. We use the full range here.)
            for i_c in lo.0..=hi.0 {
                let i_f = 2 * (i_c - lo.0);
                for k_c in lo.2..hi.2 {
                    let k_f0 = 2 * (k_c - lo.2);
                    let hx_avg = 0.5
                        * (fine_hx_t[(i_f, fine_j_h, k_f0)] + fine_hx_t[(i_f, fine_j_h, k_f0 + 1)]);
                    let hx_eff = if do_ghost {
                        hx_avg - parent.hx[(i_c, j_c_inside_h, k_c)]
                    } else {
                        hx_avg
                    };
                    parent.ez[(i_c, j_c_surface, k_c)] += sign * ce * hx_eff;
                }
            }

            // E_x on the surface: coarse ex[(i_c, j_c_surface, k_c)],
            // i_c ∈ [lo.0, hi.0), k_c ∈ [lo.2, hi.2].
            for i_c in lo.0..hi.0 {
                let i_f0 = 2 * (i_c - lo.0);
                for k_c in lo.2..=hi.2 {
                    let k_f = 2 * (k_c - lo.2);
                    let hz_avg = 0.5
                        * (fine_hz_t[(i_f0, fine_j_h, k_f)] + fine_hz_t[(i_f0 + 1, fine_j_h, k_f)]);
                    let hz_eff = if do_ghost {
                        hz_avg - parent.hz[(i_c, j_c_inside_h, k_c)]
                    } else {
                        hz_avg
                    };
                    parent.ex[(i_c, j_c_surface, k_c)] -= sign * ce * hz_eff;
                }
            }
        }

        if do_m {
            // -------- M term: write to H_z, H_x on the layer just inside. --------
            //
            // Phase 2.fdtd.7.y Option β: see [`Self::inject_x_face`]
            // for the compensating-source rationale. When
            // `use_compensating_source = true`, each fine-E sample is
            // `E_post − E_pre` cell-by-cell.

            // H_z on the layer: coarse hz[(i_c, j_c_inside_h, k_c)],
            // i_c ∈ [lo.0, hi.0), k_c ∈ [lo.2, hi.2].
            // Source: fine E_x on the surface (fine_j = fine_j_e).
            for i_c in lo.0..hi.0 {
                let i_f0 = 2 * (i_c - lo.0);
                for k_c in lo.2..=hi.2 {
                    let k_f = 2 * (k_c - lo.2);
                    let ex0 = fine_ex_for_m[(i_f0, fine_j_e, k_f)];
                    let ex1 = fine_ex_for_m[(i_f0 + 1, fine_j_e, k_f)];
                    let ex_eff0 = if use_compensating_source {
                        ex0 - fine_ex_pre_for_m[(i_f0, fine_j_e, k_f)]
                    } else {
                        ex0
                    };
                    let ex_eff1 = if use_compensating_source {
                        ex1 - fine_ex_pre_for_m[(i_f0 + 1, fine_j_e, k_f)]
                    } else {
                        ex1
                    };
                    let ex_avg = 0.5 * (ex_eff0 + ex_eff1);
                    parent.hz[(i_c, j_c_inside_h, k_c)] -= sign * ch * ex_avg;
                }
            }

            // H_x on the layer: coarse hx[(i_c, j_c_inside_h, k_c)],
            // i_c ∈ [lo.0, hi.0], k_c ∈ [lo.2, hi.2).
            // Source: fine E_z on the surface (fine_j = fine_j_e).
            for i_c in lo.0..=hi.0 {
                let i_f = 2 * (i_c - lo.0);
                for k_c in lo.2..hi.2 {
                    let k_f0 = 2 * (k_c - lo.2);
                    let ez0 = fine_ez_for_m[(i_f, fine_j_e, k_f0)];
                    let ez1 = fine_ez_for_m[(i_f, fine_j_e, k_f0 + 1)];
                    let ez_eff0 = if use_compensating_source {
                        ez0 - fine_ez_pre_for_m[(i_f, fine_j_e, k_f0)]
                    } else {
                        ez0
                    };
                    let ez_eff1 = if use_compensating_source {
                        ez1 - fine_ez_pre_for_m[(i_f, fine_j_e, k_f0 + 1)]
                    } else {
                        ez1
                    };
                    let ez_avg = 0.5 * (ez_eff0 + ez_eff1);
                    parent.hx[(i_c, j_c_inside_h, k_c)] += sign * ch * ez_avg;
                }
            }
        }
    }

    /// Inject corrections on one `±z` face. Mirror of
    /// [`Self::inject_x_face`] with cyclic axis permutation
    /// `(x, y, z) → (z, x, y)`:
    ///
    /// ```text
    /// J_x = -sign · H_y        ⇒ E_x += +sign · (dt/(ε·dz)) · H_y_fine
    /// J_y = +sign · H_x        ⇒ E_y += -sign · (dt/(ε·dz)) · H_x_fine
    /// M_x = +sign · E_y        ⇒ H_x += -sign · (dt/(μ·dz)) · E_y_fine
    /// M_y = -sign · E_x        ⇒ H_y += +sign · (dt/(μ·dz)) · E_x_fine
    /// ```
    ///
    /// (`J = sign · ẑ × H = sign · (-Hy, Hx, 0)`,
    /// `M = -sign · ẑ × E = -sign · (-Ey, Ex, 0)`.)
    #[allow(clippy::too_many_arguments)]
    fn inject_z_face(
        parent: &mut YeeGrid,
        lo: (usize, usize, usize),
        hi: (usize, usize, usize),
        fine_ex_pre_for_m: &Array3<f64>,
        fine_ey_pre_for_m: &Array3<f64>,
        fine_ex_for_m: &Array3<f64>,
        fine_ey_for_m: &Array3<f64>,
        fine_hx_t: &Array3<f64>,
        fine_hy_t: &Array3<f64>,
        sign: f64,
        k_c_surface: usize,
        k_c_inside_h: usize,
        fine_k_h: usize,
        fine_k_e: usize,
        coeff_e: f64,
        coeff_h: f64,
        do_j: bool,
        do_m: bool,
        do_ghost: bool,
        use_compensating_source: bool,
    ) {
        // ADR-0038: M-side coarse-ghost subtraction is incompatible
        // with compensating-source sampling. See `inject_x_face` for
        // the rationale.
        debug_assert!(
            !(do_ghost && use_compensating_source),
            "inject_z_face: do_ghost and use_compensating_source are mutually exclusive on M"
        );
        let dz = parent.dz;
        let ce = coeff_e / dz;
        let ch = coeff_h / dz;

        if do_j {
            // -------- J term: write to E_x, E_y on the surface plane. --------
            //
            // Berenger 2006 §III equivalent-source form (B2.2,
            // `do_ghost = true`):
            //   J = +n̂ × (H_TF_fine − H_SF_coarse_ghost)
            // For target `ex[(i, j, k_c_surface)]` the coarse natural
            // curl reads `hy[(i, j, k_c_inside_h)]` from inside the
            // footprint — that's the ghost. See [`Self::inject_x_face`]
            // for the rationale on `do_ghost`.

            // E_x on the surface: coarse ex[(i_c, j_c, k_c_surface)],
            // i_c ∈ [lo.0, hi.0), j_c ∈ [lo.1, hi.1].
            for i_c in lo.0..hi.0 {
                let i_f0 = 2 * (i_c - lo.0);
                for j_c in lo.1..=hi.1 {
                    let j_f = 2 * (j_c - lo.1);
                    let hy_avg = 0.5
                        * (fine_hy_t[(i_f0, j_f, fine_k_h)] + fine_hy_t[(i_f0 + 1, j_f, fine_k_h)]);
                    let hy_eff = if do_ghost {
                        hy_avg - parent.hy[(i_c, j_c, k_c_inside_h)]
                    } else {
                        hy_avg
                    };
                    parent.ex[(i_c, j_c, k_c_surface)] += sign * ce * hy_eff;
                }
            }

            // E_y on the surface: coarse ey[(i_c, j_c, k_c_surface)],
            // i_c ∈ [lo.0, hi.0], j_c ∈ [lo.1, hi.1).
            for i_c in lo.0..=hi.0 {
                let i_f = 2 * (i_c - lo.0);
                for j_c in lo.1..hi.1 {
                    let j_f0 = 2 * (j_c - lo.1);
                    let hx_avg = 0.5
                        * (fine_hx_t[(i_f, j_f0, fine_k_h)] + fine_hx_t[(i_f, j_f0 + 1, fine_k_h)]);
                    let hx_eff = if do_ghost {
                        hx_avg - parent.hx[(i_c, j_c, k_c_inside_h)]
                    } else {
                        hx_avg
                    };
                    parent.ey[(i_c, j_c, k_c_surface)] -= sign * ce * hx_eff;
                }
            }
        }

        if do_m {
            // -------- M term: write to H_x, H_y on the layer just inside. --------
            //
            // Phase 2.fdtd.7.y Option β: see [`Self::inject_x_face`]
            // for the compensating-source rationale. When
            // `use_compensating_source = true`, each fine-E sample is
            // `E_post − E_pre` cell-by-cell.

            // H_x on the layer: coarse hx[(i_c, j_c, k_c_inside_h)],
            // i_c ∈ [lo.0, hi.0], j_c ∈ [lo.1, hi.1).
            // Source: fine E_y on the surface (fine_k = fine_k_e).
            for i_c in lo.0..=hi.0 {
                let i_f = 2 * (i_c - lo.0);
                for j_c in lo.1..hi.1 {
                    let j_f0 = 2 * (j_c - lo.1);
                    let ey0 = fine_ey_for_m[(i_f, j_f0, fine_k_e)];
                    let ey1 = fine_ey_for_m[(i_f, j_f0 + 1, fine_k_e)];
                    let ey_eff0 = if use_compensating_source {
                        ey0 - fine_ey_pre_for_m[(i_f, j_f0, fine_k_e)]
                    } else {
                        ey0
                    };
                    let ey_eff1 = if use_compensating_source {
                        ey1 - fine_ey_pre_for_m[(i_f, j_f0 + 1, fine_k_e)]
                    } else {
                        ey1
                    };
                    let ey_avg = 0.5 * (ey_eff0 + ey_eff1);
                    parent.hx[(i_c, j_c, k_c_inside_h)] -= sign * ch * ey_avg;
                }
            }

            // H_y on the layer: coarse hy[(i_c, j_c, k_c_inside_h)],
            // i_c ∈ [lo.0, hi.0), j_c ∈ [lo.1, hi.1].
            // Source: fine E_x on the surface (fine_k = fine_k_e).
            for i_c in lo.0..hi.0 {
                let i_f0 = 2 * (i_c - lo.0);
                for j_c in lo.1..=hi.1 {
                    let j_f = 2 * (j_c - lo.1);
                    let ex0 = fine_ex_for_m[(i_f0, j_f, fine_k_e)];
                    let ex1 = fine_ex_for_m[(i_f0 + 1, j_f, fine_k_e)];
                    let ex_eff0 = if use_compensating_source {
                        ex0 - fine_ex_pre_for_m[(i_f0, j_f, fine_k_e)]
                    } else {
                        ex0
                    };
                    let ex_eff1 = if use_compensating_source {
                        ex1 - fine_ex_pre_for_m[(i_f0 + 1, j_f, fine_k_e)]
                    } else {
                        ex1
                    };
                    let ex_avg = 0.5 * (ex_eff0 + ex_eff1);
                    parent.hy[(i_c, j_c, k_c_inside_h)] += sign * ch * ex_avg;
                }
            }
        }
    }

    /// Enumerate the coarse-cell `(i, j, k)` triples on each of the six
    /// Huygens faces of this subgrid region, returned in the order
    /// `[XLow, XHigh, YLow, YHigh, ZLow, ZHigh]` matching
    /// [`BerengerHuygensFace::all`].
    ///
    /// Each face is a 2-D slice of coarse cells in the parent grid: the
    /// cells immediately *outside* the fine box in the direction of the
    /// outward normal. Concretely:
    ///
    /// - `XLow`  — `i_c = lo.0 − 1`, `j_c ∈ [lo.1, hi.1)`, `k_c ∈ [lo.2, hi.2)`.
    /// - `XHigh` — `i_c = hi.0`,     `j_c ∈ [lo.1, hi.1)`, `k_c ∈ [lo.2, hi.2)`.
    /// - `YLow`  — `j_c = lo.1 − 1`, `i_c ∈ [lo.0, hi.0)`, `k_c ∈ [lo.2, hi.2)`.
    /// - `YHigh` — `j_c = hi.1`,     `i_c ∈ [lo.0, hi.0)`, `k_c ∈ [lo.2, hi.2)`.
    /// - `ZLow`  — `k_c = lo.2 − 1`, `i_c ∈ [lo.0, hi.0)`, `j_c ∈ [lo.1, hi.1)`.
    /// - `ZHigh` — `k_c = hi.2`,     `i_c ∈ [lo.0, hi.0)`, `j_c ∈ [lo.1, hi.1)`.
    ///
    /// Each face owns
    /// `(hi.a − lo.a) · (hi.b − lo.b)` cells, where `a`, `b` are the two
    /// tangential axes for that face. For a `3 × 3 × 3` coarse subgrid
    /// (e.g. `lo = (2, 2, 2)`, `hi = (5, 5, 5)`) every face has 9 cells.
    ///
    /// Cuboid edges (where two faces meet) are owned by exactly one of
    /// the two adjacent faces under the lower-numbered-axis-wins rule
    /// — see [`assign_edge_to_face`]. This per-face index table
    /// enumerates each face's full tangential extent independently;
    /// the edge-ownership tie-break is consumed at the J/M
    /// accumulation site in B2 to avoid double-counting at the 12
    /// cuboid edges.
    ///
    /// **Note on `XLow`/`YLow`/`ZLow` indices** — `lo.a − 1` is the
    /// outward-side coarse cell. If `lo.a == 0` the subgrid would touch
    /// the parent grid boundary, which the constructor already rejects
    /// under CPML co-location (see [`SubgridContext::cpml_thickness`])
    /// for any physically meaningful run. This helper does not guard
    /// against that case and would underflow if `lo.a == 0`; callers
    /// invoking it on a degenerate region get an arithmetic-underflow
    /// panic, which is the desired loud-failure mode in debug builds.
    pub fn face_index_table(&self) -> [Vec<(usize, usize, usize)>; 6] {
        let lo = self.lo;
        let hi = self.hi;

        let mut x_low = Vec::with_capacity((hi.1 - lo.1) * (hi.2 - lo.2));
        let mut x_high = Vec::with_capacity((hi.1 - lo.1) * (hi.2 - lo.2));
        let mut y_low = Vec::with_capacity((hi.0 - lo.0) * (hi.2 - lo.2));
        let mut y_high = Vec::with_capacity((hi.0 - lo.0) * (hi.2 - lo.2));
        let mut z_low = Vec::with_capacity((hi.0 - lo.0) * (hi.1 - lo.1));
        let mut z_high = Vec::with_capacity((hi.0 - lo.0) * (hi.1 - lo.1));

        // XLow / XHigh — outward-side i_c, tangential (j, k).
        let i_low = lo.0 - 1;
        let i_high = hi.0;
        for j_c in lo.1..hi.1 {
            for k_c in lo.2..hi.2 {
                x_low.push((i_low, j_c, k_c));
                x_high.push((i_high, j_c, k_c));
            }
        }

        // YLow / YHigh — outward-side j_c, tangential (i, k).
        let j_low = lo.1 - 1;
        let j_high = hi.1;
        for i_c in lo.0..hi.0 {
            for k_c in lo.2..hi.2 {
                y_low.push((i_c, j_low, k_c));
                y_high.push((i_c, j_high, k_c));
            }
        }

        // ZLow / ZHigh — outward-side k_c, tangential (i, j).
        let k_low = lo.2 - 1;
        let k_high = hi.2;
        for i_c in lo.0..hi.0 {
            for j_c in lo.1..hi.1 {
                z_low.push((i_c, j_c, k_low));
                z_high.push((i_c, j_c, k_high));
            }
        }

        [x_low, x_high, y_low, y_high, z_low, z_high]
    }
}

/// Subgridded FDTD driver wrapping a [`WalkingSkeletonSolver`] and at most
/// one [`SubgridRegion`].
///
/// At Phase 2.fdtd.7.0 Q2 (this step) [`Self::step`] is a *placeholder*
/// that delegates straight to the wrapped solver's [`FdtdSolver::step`]
/// implementation — the fine grid does not influence the coarse fields.
/// This keeps the type surface stable so the Q3, Q4, Q5 tracks can wire
/// coarse ↔ fine coupling without re-jigging the call site.
///
/// [`FdtdSolver::step`]: crate::FdtdSolver::step
///
/// The Q1 helper-sequence refactor (`update_h_only`, `apply_cpml_h`,
/// `update_e_only`, `apply_cpml_e`, `advance_clock`) is the seam those
/// future tracks will inject the seven-stage spec §3 sequence into.
pub struct SubgriddedSolver {
    inner: WalkingSkeletonSolver,
    region: Option<SubgridRegion>,
}

impl SubgriddedSolver {
    /// Wrap a [`WalkingSkeletonSolver`] with no subgrid region attached.
    ///
    /// In this configuration [`Self::step`] is bit-for-bit identical to
    /// the wrapped solver's own [`FdtdSolver::step`] implementation.
    ///
    /// [`FdtdSolver::step`]: crate::FdtdSolver::step
    pub fn new(solver: WalkingSkeletonSolver) -> Self {
        Self {
            inner: solver,
            region: None,
        }
    }

    /// Attach a [`SubgridRegion`] to this solver, consuming `self`.
    ///
    /// The region is held but, at Q2, is not yet stepped. Q5 lands the
    /// seven-stage interleave that activates it.
    #[must_use]
    pub fn with_region(mut self, region: SubgridRegion) -> Self {
        self.region = Some(region);
        self
    }

    /// Advance the simulation by one coarse step.
    ///
    /// Phase 2.fdtd.7.x B2.2 — Berenger 2006 Huygens-surface closure
    /// with **split J / M injection** and **coarse-ghost-subtracted J
    /// source** (canonical Berenger §III equivalent-current form). The
    /// J injection now reads `J = +n̂ × (H_fine − H_coarse_ghost)`
    /// where `H_coarse_ghost` is the coarse `H` slot inside the fine-
    /// box footprint that the coarse natural-curl stencil already reads
    /// for its leapfrog update — converting the source from full TF
    /// to the canonical scattered-correction form (LLLLLLL's B2.1
    /// diagnosis identified the missing-ghost-subtraction as the root
    /// cause of the residual divergence after the split-injection
    /// reorder). The M side remains un-ghosted because its candidate
    /// ghost slot (coarse `E` at the surface plane) is Dirichlet-tied
    /// to the fine grid via Q3, so subtracting it would nullify M; see
    /// [`SubgridRegion::inject_m_to_coarse_h`] for the empirical
    /// rationale.
    ///
    /// The split-injection pipeline interleaves the J and M injections
    /// with the coarse leapfrog so each source enters its respective
    /// `update_*_only` stencil *before* the update runs, rather than
    /// being stacked on top of an already-updated field:
    ///
    /// 1. **M-source injection (deferred from previous step).**
    ///    [`SubgridRegion::inject_m_to_coarse_h`] applies
    ///    `M^n = -n̂ × E_fine^n` (sampled from the previous step's
    ///    [`SubgridRegion::snapshot_fine_e_end_of_step`]) onto coarse
    ///    `H_t`. No-op on the first coarse step (no snapshot yet).
    /// 2. Coarse `update_h_only` consumes the M source as it advances
    ///    coarse `H^{n−1/2}` → `H^{n+1/2}`.
    /// 3. Coarse `apply_cpml_h`.
    /// 4. Snapshot the coarse `E_t` (start of step).
    /// 5. Fine sub-step `k = 1`: `interpolate_coarse_e_to_fine(0.25)`,
    ///    `update_fine_h`, **`snapshot_fine_h_mid_step`** (Q4.1; captures
    ///    fine `H^{n+1/4}` for the J time-centering average),
    ///    `update_fine_e`.
    /// 6. Coarse `update_e_only` is **not yet called** — the J source
    ///    must enter `update_e_only` so we run fine sub-step 2 to a
    ///    state where `H_fine^{n+1/2}` is well-defined first.
    /// 7. Snapshot the coarse `E_t` (end-of-step proxy: equals start
    ///    because coarse `update_e_only` has not run yet). The fine
    ///    sub-step 2 reads the same value at `frac = 0.75` that
    ///    sub-step 1 read at `frac = 0.25`; the Berenger closure does
    ///    not depend on temporal blending here (it depends on the
    ///    equivalent-current injection instead).
    /// 8. Fine sub-step `k = 2`: `interpolate_coarse_e_to_fine(0.75)`,
    ///    `update_fine_h`, `update_fine_e`. Fine `H` now at `t = n + 3/4`,
    ///    fine `E` at `t = n + 1`.
    /// 9. **J-source injection** —
    ///    [`SubgridRegion::inject_j_to_coarse_e`] applies
    ///    `J^{n+1/2} = +n̂ × ((H_fine^{n+1/4} + H_fine^{n+3/4}) / 2)`
    ///    onto coarse `E_t` *before* coarse `update_e_only`.
    /// 10. Coarse `update_e_only` consumes the J source as it advances
    ///     coarse `E^n` → `E^{n+1}`.
    /// 11. Coarse `apply_cpml_e`.
    /// 12. [`SubgridRegion::snapshot_fine_e_end_of_step`] captures fine
    ///     `E^{n+1}` for the next coarse step's M injection.
    /// 13. Advance the coarse clock by one `dt_coarse`.
    ///
    /// Replaces the prior B2 monolithic
    /// `inject_equivalent_currents_to_coarse` call site that injected
    /// both J and M **after** both `update_*_only` stages had completed.
    /// Per HHHHHHH's diagnosis, that ordering closed a unit-magnitude
    /// feedback loop: fine E → coarse E (via J injected post-`update_e`)
    /// → coarse H (via next step's `update_h`) → fine E_t Dirichlet
    /// (via Q3 interpolation) → larger fine E updated by larger fine H
    /// → larger J. Decoupling J/M into the pre-update stencil breaks
    /// the loop.
    ///
    /// When no [`SubgridRegion`] is attached this collapses to the
    /// bare `WalkingSkeletonSolver` helper sequence, matching
    /// [`WalkingSkeletonSolver::step`] bit-for-bit.
    pub fn step(&mut self) {
        let Some(region) = self.region.as_mut() else {
            self.inner.update_h_only();
            self.inner.apply_cpml_h();
            self.inner.update_e_only();
            self.inner.apply_cpml_e();
            self.inner.advance_clock();
            return;
        };

        // Phase 2.fdtd.7.x B2.1 split-injection pipeline. The stage
        // ordering breaks HHHHHHH's feedback loop diagnosis (commit
        // `997e706`) — the J source is decoupled from the Q3
        // coarse → fine Dirichlet propagation by deferring its
        // application until **after** both fine sub-steps have read
        // coarse `E` and after coarse `update_e_only` has advanced its
        // natural-curl-only step. The M source is similarly deferred
        // by one coarse step so it lands on coarse `H` before the
        // next `update_h_only` consumes it.

        // 1. Deferred M injection from previous step's end-of-fine-E
        //    snapshot (no-op on the first coarse step). The M source
        //    enters coarse `H` *before* `update_h_only` runs so the
        //    next H half-step consumes it inside the leapfrog stencil
        //    rather than stacking on top of an already-updated H.
        {
            let (parent_grid, _) = self.inner.grid_and_cpml_mut();
            region.inject_m_to_coarse_h(parent_grid);
        }

        // 2–3. Coarse H half-step (consumes M^n) and outer-boundary
        //      closure.
        self.inner.update_h_only();
        self.inner.apply_cpml_h();

        // 4. Snapshot the parent E_t at the start of the coarse step.
        //    Taken before the coarse-E update so fine sub-step 1
        //    Dirichlet reads coarse E^n.
        region.snapshot_coarse_e_t(self.inner.grid());

        // 5. Fine sub-step k = 1. Phase 2.fdtd.7.y Step C5 (Option α):
        //    the Q3 coarse → fine Dirichlet `interpolate_coarse_e_to_fine`
        //    is replaced by a 1st-order Mur ABC on the fine outer `E_t`
        //    plane. The Mur sandwich is:
        //      (a) snapshot adjacent-inside + boundary fine E_t at t = n
        //          (`snapshot_fine_e_for_mur`),
        //      (b) advance fine H then fine E via the usual updates
        //          (`update_fine_h` / `update_fine_e`); the latter
        //          deliberately skips the outer boundary cells per
        //          `crate::update::update_e`,
        //      (c) write the new boundary value via the Mur formula
        //          using the snapshot + the post-update adjacent-inside
        //          value (`apply_mur_abc_to_fine_outer_e`).
        //    The mid-step fine-H snapshot stays where the spec §3
        //    seven-stage sequence placed it (between `update_fine_h` and
        //    `update_fine_e` for the Q4.1 time-centred J source).
        region.snapshot_fine_e_for_mur();
        region.update_fine_h();
        region.snapshot_fine_h_mid_step();
        region.update_fine_e();
        region.apply_mur_abc_to_fine_outer_e();

        // 6–7. Coarse E half-step *without* the J source. update_e_only
        //      advances coarse E using only the natural curl(H^{n+1/2}).
        //      The fine outer `E_t` is now absorbing (Mur) rather than
        //      Dirichlet-coupled to the coarse boundary, so the
        //      sub-step 2 reads of the coarse E_t snapshot remain in
        //      place only as a no-op pair for the snapshot accessor
        //      contract (Q3 helper kept available for tests and for
        //      rollback per ADR-0038).
        self.inner.update_e_only();
        self.inner.apply_cpml_e();

        // 8. End-of-coarse-E snapshot. Retained for accessor-test
        //    compatibility (the C5 Mur path does not read it; the
        //    bracketing snapshot pair is harmless and cheap).
        region.snapshot_coarse_e_t_end(self.inner.grid());

        // 9. Fine sub-step k = 2. Same Mur sandwich as sub-step 1, with
        //    the additional Phase 2.fdtd.7.y compensating-source pair
        //    `snapshot_fine_e_pre_update` / `snapshot_fine_e_post_update`
        //    around the inner `update_fine_e` + `apply_mur_abc_*` block.
        //    Crucially `E_pre` is captured *before* `update_fine_e` and
        //    `E_post` *after* `apply_mur_abc_to_fine_outer_e`, so the
        //    Mur-written boundary value is included in `E_post`. With
        //    Q3 retired, `E_post − E_pre` on the surface plane is now
        //    governed entirely by the Mur ABC (which depends on the
        //    adjacent-inside fine `E_t`), recovering non-zero
        //    differencing for the canonical M source per ADR-0038
        //    Option α.
        region.snapshot_fine_e_for_mur();
        region.snapshot_fine_e_pre_update();
        region.update_fine_h();
        region.update_fine_e();
        region.apply_mur_abc_to_fine_outer_e();
        region.snapshot_fine_e_post_update();

        // 10. J-source injection — applies J^{n+1/2} to coarse E^{n+1}
        //     *after* the Q3 snapshots and both fine sub-steps have
        //     read coarse E. The J modifies coarse E^{n+1} for the
        //     coarse leapfrog (so the next coarse `update_h_only` sees
        //     a J-augmented coarse `E`, as Maxwell's equations require),
        //     but the fine grid's next-step Dirichlet feed will
        //     re-snapshot post-`update_e` E *before* re-applying J at
        //     the same point in the next step's pipeline.
        {
            let (parent_grid, _) = self.inner.grid_and_cpml_mut();
            region.inject_j_to_coarse_e(parent_grid);
        }

        // 11. Snapshot fine E for next step's M injection.
        region.snapshot_fine_e_end_of_step();

        // 12. Advance the coarse clock.
        self.inner.advance_clock();
    }

    /// Advance the simulation by one coarse step while injecting a
    /// Gaussian-in-time pulse on the coarse `E_z` array at cell
    /// `(i, j, k)`.
    ///
    /// The source is added **only to the coarse grid**, between the
    /// coarse H-update / CPML-H closure and the coarse E-update — the
    /// same leapfrog timing that [`WalkingSkeletonSolver::step_with_source`]
    /// uses. The fine grid in Phase 2.fdtd.7.0 is sourceless; energy
    /// crosses into the fine region exclusively through the Dirichlet
    /// boundary `E_t` interpolation, which is the property exercised by
    /// the `subgrid_plane_wave_traversal` integration gate.
    ///
    /// The Gaussian is sampled at the current simulation time (before
    /// this step advances the coarse clock).
    pub fn step_with_gaussian_source_ez(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        t0: f64,
        sigma: f64,
    ) {
        let t = self.inner.current_time();
        let Some(region) = self.region.as_mut() else {
            self.inner.update_h_only();
            self.inner.apply_cpml_h();
            self.inner.apply_gaussian_source_ez(i, j, k, t, t0, sigma);
            self.inner.update_e_only();
            self.inner.apply_cpml_e();
            self.inner.advance_clock();
            return;
        };

        // Phase 2.fdtd.7.x B2.1 split-injection pipeline — see
        // [`Self::step`] for the stage-by-stage rationale.

        // 1. Deferred M injection from previous step.
        {
            let (parent_grid, _) = self.inner.grid_and_cpml_mut();
            region.inject_m_to_coarse_h(parent_grid);
        }

        // 2–3. Coarse H half-step (consumes M) and outer-boundary
        //      closure.
        self.inner.update_h_only();
        self.inner.apply_cpml_h();

        // 4. Start-of-step snapshot (taken before the coarse-E update).
        region.snapshot_coarse_e_t(self.inner.grid());

        // 5. Fine sub-step k = 1. Phase 2.fdtd.7.y Step C5 (Option α):
        //    Mur ABC replaces the Q3 coarse → fine Dirichlet write. See
        //    [`Self::step`] for the stage-by-stage rationale.
        region.snapshot_fine_e_for_mur();
        region.update_fine_h();
        region.snapshot_fine_h_mid_step();
        region.update_fine_e();
        region.apply_mur_abc_to_fine_outer_e();

        // 6. Coarse Gaussian source injection (between fine sub-step 1
        //    and the coarse E half-step — standard leapfrog timing).
        self.inner.apply_gaussian_source_ez(i, j, k, t, t0, sigma);

        // 7–8. Coarse E half-step *without* the Berenger J source, and
        //      CPML-E closure. update_e_only consumes the Gaussian
        //      source plus the natural curl(H^{n+1/2}); J^{n+1/2} is
        //      injected later (stage 11).
        self.inner.update_e_only();
        self.inner.apply_cpml_e();

        // 9. End-of-coarse-E snapshot retained for accessor contract
        //    only (Mur path does not consume it).
        region.snapshot_coarse_e_t_end(self.inner.grid());

        // 10. Fine sub-step k = 2 with the Mur sandwich plus the
        //     compensating-source `E_pre` / `E_post` capture. See
        //     [`Self::step`] for the rationale.
        region.snapshot_fine_e_for_mur();
        region.snapshot_fine_e_pre_update();
        region.update_fine_h();
        region.update_fine_e();
        region.apply_mur_abc_to_fine_outer_e();
        region.snapshot_fine_e_post_update();

        // 11. J-source injection — applies J^{n+1/2} to coarse E^{n+1}
        //     after both fine sub-steps have read coarse E (so the
        //     Q3 propagation path is not contaminated).
        {
            let (parent_grid, _) = self.inner.grid_and_cpml_mut();
            region.inject_j_to_coarse_e(parent_grid);
        }

        // 12. Snapshot fine E for next step's M injection.
        region.snapshot_fine_e_end_of_step();

        // 13. Advance the coarse clock.
        self.inner.advance_clock();
    }

    /// Immutable borrow of the wrapped coarse solver.
    pub fn inner(&self) -> &WalkingSkeletonSolver {
        &self.inner
    }

    /// Mutable borrow of the wrapped coarse solver (e.g. to inject
    /// per-step coarse-grid sources between [`Self::step`] calls).
    pub fn inner_mut(&mut self) -> &mut WalkingSkeletonSolver {
        &mut self.inner
    }

    /// Immutable borrow of the attached [`SubgridRegion`], if any.
    pub fn region(&self) -> Option<&SubgridRegion> {
        self.region.as_ref()
    }

    /// Mutable borrow of the attached [`SubgridRegion`], if any.
    pub fn region_mut(&mut self) -> Option<&mut SubgridRegion> {
        self.region.as_mut()
    }
}
