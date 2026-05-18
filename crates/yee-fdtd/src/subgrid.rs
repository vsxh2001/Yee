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

use ndarray::Array2;
use yee_core::Error;

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
    pub fn average_fine_h_to_coarse(&self, parent: &mut YeeGrid) {
        let lo = self.lo;
        let hi = self.hi;
        let fine_nx = 2 * (hi.0 - lo.0);
        let fine_ny = 2 * (hi.1 - lo.1);
        let fine_nz = 2 * (hi.2 - lo.2);

        // ±x faces — overwrite coarse H_y, H_z on the layer i_c ∈ {lo.0, hi.0 − 1}.
        Self::avg_face_x(&self.fine, parent, lo, hi, lo.0, 0);
        Self::avg_face_x(&self.fine, parent, lo, hi, hi.0 - 1, fine_nx - 2);

        // ±y faces — overwrite coarse H_x, H_z on the layer j_c ∈ {lo.1, hi.1 − 1}.
        Self::avg_face_y(&self.fine, parent, lo, hi, lo.1, 0);
        Self::avg_face_y(&self.fine, parent, lo, hi, hi.1 - 1, fine_ny - 2);

        // ±z faces — overwrite coarse H_x, H_y on the layer k_c ∈ {lo.2, hi.2 − 1}.
        Self::avg_face_z(&self.fine, parent, lo, hi, lo.2, 0);
        Self::avg_face_z(&self.fine, parent, lo, hi, hi.2 - 1, fine_nz - 2);
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
    /// (`0` for the −x face, `fine_nx − 2` for the +x face).
    fn avg_face_x(
        fine: &YeeGrid,
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
                let s = fine.hy[(fine_i_lo, j_f, k_f0)]
                    + fine.hy[(fine_i_lo + 1, j_f, k_f0)]
                    + fine.hy[(fine_i_lo, j_f, k_f0 + 1)]
                    + fine.hy[(fine_i_lo + 1, j_f, k_f0 + 1)];
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
                let s = fine.hz[(fine_i_lo, j_f0, k_f)]
                    + fine.hz[(fine_i_lo + 1, j_f0, k_f)]
                    + fine.hz[(fine_i_lo, j_f0 + 1, k_f)]
                    + fine.hz[(fine_i_lo + 1, j_f0 + 1, k_f)];
                parent.hz[(i_c_face, j_c, k_c)] = 0.25 * s;
            }
        }
    }

    /// Area-average fine `H_x`, `H_z` onto a coarse layer adjacent to a
    /// ±y interface face.
    fn avg_face_y(
        fine: &YeeGrid,
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
                let s = fine.hx[(i_f, fine_j_lo, k_f0)]
                    + fine.hx[(i_f, fine_j_lo + 1, k_f0)]
                    + fine.hx[(i_f, fine_j_lo, k_f0 + 1)]
                    + fine.hx[(i_f, fine_j_lo + 1, k_f0 + 1)];
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
                let s = fine.hz[(i_f0, fine_j_lo, k_f)]
                    + fine.hz[(i_f0 + 1, fine_j_lo, k_f)]
                    + fine.hz[(i_f0, fine_j_lo + 1, k_f)]
                    + fine.hz[(i_f0 + 1, fine_j_lo + 1, k_f)];
                parent.hz[(i_c, j_c_face, k_c)] = 0.25 * s;
            }
        }
    }

    /// Area-average fine `H_x`, `H_y` onto a coarse layer adjacent to a
    /// ±z interface face.
    fn avg_face_z(
        fine: &YeeGrid,
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
                let s = fine.hx[(i_f, j_f0, fine_k_lo)]
                    + fine.hx[(i_f, j_f0 + 1, fine_k_lo)]
                    + fine.hx[(i_f, j_f0, fine_k_lo + 1)]
                    + fine.hx[(i_f, j_f0 + 1, fine_k_lo + 1)];
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
                let s = fine.hy[(i_f0, j_f, fine_k_lo)]
                    + fine.hy[(i_f0 + 1, j_f, fine_k_lo)]
                    + fine.hy[(i_f0, j_f, fine_k_lo + 1)]
                    + fine.hy[(i_f0 + 1, j_f, fine_k_lo + 1)];
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
    /// Implements the seven-stage spec §3 time-step sequence:
    ///
    /// 1. Coarse `update_h_only`.
    /// 2. Coarse `apply_cpml_h` (no-op outside any configured CPML face).
    /// 3. Snapshot the coarse `E_t` on the six interface faces (start).
    /// 4. Fine sub-step `k = 1`: interpolate coarse `E_t` at `frac = 0.25`
    ///    onto the fine boundary, then bulk `update_h` and `update_e` on
    ///    the fine grid.
    /// 5. Coarse `update_e_only`, then `apply_cpml_e`.
    /// 6. Snapshot the coarse `E_t` again (end-of-coarse-step) so the
    ///    second fine sub-step blends against the post-E-update parent.
    /// 7. Fine sub-step `k = 2`: interpolate at `frac = 0.75`, then bulk
    ///    `update_h` and `update_e` on the fine grid.
    /// 8. Average the fine `H_t` onto the coarse interface layer
    ///    (`average_fine_h_to_coarse`) and overwrite the coarse `E_t` on
    ///    the interface planes (`overwrite_coarse_e_from_fine`) —
    ///    Chevalier 1997 §IV closure of the discrete energy balance.
    /// 9. Advance the coarse clock by one `dt_coarse`.
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

        // Per the brief, the two fine sub-steps bracket the coarse
        // E-update with frac = 0.25 / 0.75 against the start/end
        // snapshots. The start snapshot is taken at t = n (after the
        // coarse H-half-step but before the coarse E-half-step) and the
        // end snapshot at t = n + 1 (after the coarse E-half-step). On
        // the first sub-step the not-yet-written end buffer carries the
        // previous coarse interval's end value, which by construction
        // equals the current start value — so the 0.25 blend collapses
        // to the start value, leaving the fine sub-step momentarily
        // first-order in time at the boundary while sub-step 2 recovers
        // second order. See spec §3.

        // 1–2. Coarse H half-step and outer-boundary closure.
        self.inner.update_h_only();
        self.inner.apply_cpml_h();

        // 3. Snapshot the parent E_t at the start of the coarse step.
        region.snapshot_coarse_e_t(self.inner.grid());

        // 4. Fine sub-step k = 1 (interpolate at t = n + 1/4).
        region.interpolate_coarse_e_to_fine(0.25);
        region.update_fine_h();
        region.update_fine_e();

        // 5. Coarse E half-step and outer-boundary closure.
        self.inner.update_e_only();
        self.inner.apply_cpml_e();

        // 6. Snapshot the parent E_t at the end of the coarse step.
        region.snapshot_coarse_e_t_end(self.inner.grid());

        // 7. Fine sub-step k = 2 (interpolate at t = n + 3/4).
        region.interpolate_coarse_e_to_fine(0.75);
        region.update_fine_h();
        region.update_fine_e();

        // 8. Close the energy balance: fine → coarse H area-average plus
        //    coarse E_t overwrite on the interface planes.
        {
            let (parent_grid, _) = self.inner.grid_and_cpml_mut();
            region.average_fine_h_to_coarse(parent_grid);
            region.overwrite_coarse_e_from_fine(parent_grid);
        }

        // 9. Advance the coarse clock.
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

        // 1–2. Coarse H half-step.
        self.inner.update_h_only();
        self.inner.apply_cpml_h();

        // 3. Start-of-step snapshot. Done before the source contribution
        //    so the snapshot reflects the parent's pre-source E_t. The
        //    source landing on a cell strictly interior to the subgrid
        //    region would in principle want to advance the snapshot, but
        //    the integration gate places the source upstream of the
        //    region so this ordering is correct for the v7.0 scope.
        region.snapshot_coarse_e_t(self.inner.grid());

        // 4. Fine sub-step k = 1.
        region.interpolate_coarse_e_to_fine(0.25);
        region.update_fine_h();
        region.update_fine_e();

        // 5. Coarse source injection and E half-step.
        self.inner.apply_gaussian_source_ez(i, j, k, t, t0, sigma);
        self.inner.update_e_only();
        self.inner.apply_cpml_e();

        // 6. End-of-step snapshot (after the source-modulated coarse E).
        region.snapshot_coarse_e_t_end(self.inner.grid());

        // 7. Fine sub-step k = 2.
        region.interpolate_coarse_e_to_fine(0.75);
        region.update_fine_h();
        region.update_fine_e();

        // 8. Close the energy balance.
        {
            let (parent_grid, _) = self.inner.grid_and_cpml_mut();
            region.average_fine_h_to_coarse(parent_grid);
            region.overwrite_coarse_e_from_fine(parent_grid);
        }

        // 9. Advance the coarse clock.
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
