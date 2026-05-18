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

use yee_core::Error;

use crate::WalkingSkeletonSolver;
use crate::grid::YeeGrid;

/// Coarse-cell `(lo, hi)` extent of an axis-aligned box, inclusive-low /
/// exclusive-high. Used for the TF/SF-box placement check.
pub type CoarseBox = ((usize, usize, usize), (usize, usize, usize));

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

/// Axis-aligned, cuboidal sub-region nested at 2× resolution inside a parent
/// [`YeeGrid`].
///
/// Owns its own fine `YeeGrid` instance whose cell sizes (`dx`, `dy`, `dz`)
/// and time step (`dt`) are half the parent's, sized
/// `(2·(hi.0 − lo.0), 2·(hi.1 − lo.1), 2·(hi.2 − lo.2))` cells. The fine
/// grid inherits the parent's scalar `eps_r` and `mu_r`.
///
/// At Phase 2.fdtd.7.0 Q2 (this step) the region carries no coupling state;
/// coarse ↔ fine `E_t` interpolation and `H_t` area-averaging land in Q3,
/// Q4. The fine grid is therefore inert until Q5 wires it into
/// [`SubgriddedSolver::step`].
#[derive(Debug, Clone)]
pub struct SubgridRegion {
    /// Coarse-cell index of the nest corner (inclusive lower bound).
    pub lo: (usize, usize, usize),
    /// Coarse-cell index of the nest corner (exclusive upper bound).
    pub hi: (usize, usize, usize),
    /// The fine grid backing this region. `dx_fine = dx_coarse / 2`;
    /// `dt_fine = dt_coarse / 2`.
    fine: YeeGrid,
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

        Ok(Self { lo, hi, fine })
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
    /// **Q2 placeholder:** this delegates to the wrapped
    /// [`WalkingSkeletonSolver`]'s helper sequence; the fine grid is
    /// dormant. Q5 will replace this body with the seven-stage spec §3
    /// sequence (update_h_only → apply_cpml_h → fine k=1 → update_e_only
    /// → apply_cpml_e → fine k=2 → average fine→coarse → advance_clock).
    pub fn step(&mut self) {
        // Region present but dormant at Q2; Q5 wires the fine sub-steps.
        self.inner.update_h_only();
        self.inner.apply_cpml_h();
        self.inner.update_e_only();
        self.inner.apply_cpml_e();
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
