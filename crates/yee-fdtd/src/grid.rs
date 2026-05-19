//! Yee staggered grid for 3D FDTD.
//!
//! The grid stores the six electromagnetic field components on a staggered
//! lattice as described in Taflove & Hagness, *Computational Electrodynamics*,
//! 3rd ed., §3.6. Each component lives at a half-cell offset:
//!
//! - `E_x` on edges parallel to x: shape `[nx,   ny+1, nz+1]`
//! - `E_y` on edges parallel to y: shape `[nx+1, ny,   nz+1]`
//! - `E_z` on edges parallel to z: shape `[nx+1, ny+1, nz  ]`
//! - `H_x` on faces normal to x:   shape `[nx+1, ny,   nz  ]`
//! - `H_y` on faces normal to y:   shape `[nx,   ny+1, nz  ]`
//! - `H_z` on faces normal to z:   shape `[nx,   ny,   nz+1]`
//!
//! ## Per-cell heterogeneous materials (Phase 2.fdtd.7.z infrastructure)
//!
//! The grid carries optional per-cell relative-permittivity and
//! relative-permeability arrays ([`YeeGrid::eps_r_cells`],
//! [`YeeGrid::mu_r_cells`]) and three optional per-component PEC masks
//! ([`YeeGrid::pec_mask_ex`], [`YeeGrid::pec_mask_ey`],
//! [`YeeGrid::pec_mask_ez`]). When all five are `None`, the update stencils
//! and boundary helpers fall through to the scalar `eps_r` / `mu_r` and
//! deprecated outer-face PEC clamp respectively, so any existing call site
//! that builds a grid via [`YeeGrid::vacuum`] sees bit-exact-identical
//! behaviour. Per-cell maps are opt-in via the
//! [`YeeGrid::with_eps_r_cells`] / [`YeeGrid::with_mu_r_cells`] /
//! [`YeeGrid::with_pec_mask_ex`] / [`YeeGrid::with_pec_mask_ey`] /
//! [`YeeGrid::with_pec_mask_ez`] fluent builders.
//!
//! This infrastructure unblocks fdtd-007 (the Maloney-Smith
//! dielectric-loaded slot antenna gate, see
//! `yee-validation::run_fdtd_007_maloney_smith_slot`) by providing the two
//! missing physical primitives — a substrate dielectric slab and a PEC
//! ground plane with a slot aperture — without rebuilding the rest of the
//! solver surface.

use ndarray::Array3;

use yee_core::units::C0;

/// Yee staggered grid.
///
/// The six field arrays use the standard Taflove staggering. All fields are
/// zero-initialized; sources are injected by mutating cells in place.
///
/// By default the grid is homogeneous (scalar `eps_r` / `mu_r`); per-cell
/// heterogeneity is opt-in via the `with_eps_r_cells` / `with_mu_r_cells` /
/// `with_pec_mask_e{x,y,z}` builders. See the module-level docs for the
/// back-compat contract.
#[derive(Debug, Clone)]
pub struct YeeGrid {
    /// Number of primary cells along x.
    pub nx: usize,
    /// Number of primary cells along y.
    pub ny: usize,
    /// Number of primary cells along z.
    pub nz: usize,

    /// Cell size along x (meters).
    pub dx: f64,
    /// Cell size along y (meters).
    pub dy: f64,
    /// Cell size along z (meters).
    pub dz: f64,
    /// Time step (seconds). Must be ≤ [`YeeGrid::courant_limit`].
    pub dt: f64,

    /// Electric field x-component (staggered: `[nx, ny+1, nz+1]`).
    pub ex: Array3<f64>,
    /// Electric field y-component (staggered: `[nx+1, ny, nz+1]`).
    pub ey: Array3<f64>,
    /// Electric field z-component (staggered: `[nx+1, ny+1, nz]`).
    pub ez: Array3<f64>,
    /// Magnetic field x-component (staggered: `[nx+1, ny, nz]`).
    pub hx: Array3<f64>,
    /// Magnetic field y-component (staggered: `[nx, ny+1, nz]`).
    pub hy: Array3<f64>,
    /// Magnetic field z-component (staggered: `[nx, ny, nz+1]`).
    pub hz: Array3<f64>,

    /// Relative permittivity (uniform fallback; vacuum = 1.0).
    ///
    /// When [`YeeGrid::eps_r_cells`] is `None`, this scalar drives the E
    /// update everywhere. When `Some`, the per-cell array takes precedence
    /// and this scalar is ignored by the E update.
    pub eps_r: f64,
    /// Relative permeability (uniform fallback; vacuum = 1.0). Mirror of
    /// [`YeeGrid::eps_r`] for the H update.
    pub mu_r: f64,

    /// Optional per-cell relative permittivity.
    ///
    /// Shape `[nx+1, ny+1, nz+1]`: oversized to address any staggered E
    /// component by its primary `(i, j, k)` cell index. All three E
    /// components in a single primary cell see the same ε_r, matching the
    /// convention used by [`crate::dispersive::DispersiveState`].
    pub eps_r_cells: Option<Array3<f64>>,

    /// Optional per-cell relative permeability, sibling of [`Self::eps_r_cells`]
    /// for the H update.
    pub mu_r_cells: Option<Array3<f64>>,

    /// Optional per-component PEC mask for `E_x`. Shape matches `ex`:
    /// `[nx, ny+1, nz+1]`. Any cell with `mask[i, j, k] == true` has its
    /// `E_x` clamped to zero by [`YeeGrid::apply_pec_mask`] after the
    /// standard E update.
    pub pec_mask_ex: Option<Array3<bool>>,

    /// Optional per-component PEC mask for `E_y`. Shape matches `ey`:
    /// `[nx+1, ny, nz+1]`.
    pub pec_mask_ey: Option<Array3<bool>>,

    /// Optional per-component PEC mask for `E_z`. Shape matches `ez`:
    /// `[nx+1, ny+1, nz]`.
    pub pec_mask_ez: Option<Array3<bool>>,
}

impl YeeGrid {
    /// Build a vacuum-filled cubic-cell grid of `nx × ny × nz` cells with
    /// `dx = dy = dz` and `dt` set to 0.9× the Courant limit.
    ///
    /// # Panics
    ///
    /// Panics if any dimension is zero or `dx` is non-positive.
    pub fn vacuum(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0, "grid dimensions must be > 0");
        assert!(
            dx.is_finite() && dx > 0.0,
            "cell size must be positive and finite"
        );

        let dy = dx;
        let dz = dx;

        let ex = Array3::<f64>::zeros((nx, ny + 1, nz + 1));
        let ey = Array3::<f64>::zeros((nx + 1, ny, nz + 1));
        let ez = Array3::<f64>::zeros((nx + 1, ny + 1, nz));
        let hx = Array3::<f64>::zeros((nx + 1, ny, nz));
        let hy = Array3::<f64>::zeros((nx, ny + 1, nz));
        let hz = Array3::<f64>::zeros((nx, ny, nz + 1));

        let mut grid = Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            dt: 0.0, // filled in below from the Courant limit
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
            eps_r: 1.0,
            mu_r: 1.0,
            eps_r_cells: None,
            mu_r_cells: None,
            pec_mask_ex: None,
            pec_mask_ey: None,
            pec_mask_ez: None,
        };
        grid.dt = 0.9 * grid.courant_limit();
        grid
    }

    /// Courant-Friedrichs-Lewy stability limit for the 3D Yee scheme in vacuum:
    ///
    /// ```text
    /// dt_max = 1 / (c · sqrt(1/dx² + 1/dy² + 1/dz²))
    /// ```
    ///
    /// (Taflove & Hagness eq. 4.60.)
    pub fn courant_limit(&self) -> f64 {
        let inv_sq =
            1.0 / (self.dx * self.dx) + 1.0 / (self.dy * self.dy) + 1.0 / (self.dz * self.dz);
        1.0 / (C0 * inv_sq.sqrt())
    }

    /// Attach a per-cell relative-permittivity map.
    ///
    /// The map must have shape `[nx+1, ny+1, nz+1]`; any cell value `< 1.0`
    /// would correspond to a phase velocity above `c`, so the builder
    /// silently accepts any positive value but panics on non-finite or
    /// non-positive entries.
    ///
    /// # Panics
    ///
    /// Panics if `cells.dim()` ≠ `(nx + 1, ny + 1, nz + 1)` or any entry is
    /// non-finite / ≤ 0.
    pub fn with_eps_r_cells(mut self, cells: Array3<f64>) -> Self {
        let expected = (self.nx + 1, self.ny + 1, self.nz + 1);
        assert_eq!(
            cells.dim(),
            expected,
            "eps_r_cells shape {:?} must equal grid extent {:?}",
            cells.dim(),
            expected
        );
        assert!(
            cells.iter().all(|&v| v.is_finite() && v > 0.0),
            "eps_r_cells must be finite and > 0 everywhere"
        );
        self.eps_r_cells = Some(cells);
        self
    }

    /// Attach a per-cell relative-permeability map. See
    /// [`Self::with_eps_r_cells`] for shape and validity contracts.
    pub fn with_mu_r_cells(mut self, cells: Array3<f64>) -> Self {
        let expected = (self.nx + 1, self.ny + 1, self.nz + 1);
        assert_eq!(
            cells.dim(),
            expected,
            "mu_r_cells shape {:?} must equal grid extent {:?}",
            cells.dim(),
            expected
        );
        assert!(
            cells.iter().all(|&v| v.is_finite() && v > 0.0),
            "mu_r_cells must be finite and > 0 everywhere"
        );
        self.mu_r_cells = Some(cells);
        self
    }

    /// Attach a per-component PEC mask for `E_x`. Shape must be
    /// `[nx, ny+1, nz+1]` (the `ex` array shape).
    ///
    /// # Panics
    ///
    /// Panics if `mask.dim()` ≠ `(nx, ny + 1, nz + 1)`.
    pub fn with_pec_mask_ex(mut self, mask: Array3<bool>) -> Self {
        let expected = (self.nx, self.ny + 1, self.nz + 1);
        assert_eq!(
            mask.dim(),
            expected,
            "pec_mask_ex shape {:?} must equal ex extent {:?}",
            mask.dim(),
            expected
        );
        self.pec_mask_ex = Some(mask);
        self
    }

    /// Attach a per-component PEC mask for `E_y`. Shape must be
    /// `[nx+1, ny, nz+1]`.
    pub fn with_pec_mask_ey(mut self, mask: Array3<bool>) -> Self {
        let expected = (self.nx + 1, self.ny, self.nz + 1);
        assert_eq!(
            mask.dim(),
            expected,
            "pec_mask_ey shape {:?} must equal ey extent {:?}",
            mask.dim(),
            expected
        );
        self.pec_mask_ey = Some(mask);
        self
    }

    /// Attach a per-component PEC mask for `E_z`. Shape must be
    /// `[nx+1, ny+1, nz]`.
    pub fn with_pec_mask_ez(mut self, mask: Array3<bool>) -> Self {
        let expected = (self.nx + 1, self.ny + 1, self.nz);
        assert_eq!(
            mask.dim(),
            expected,
            "pec_mask_ez shape {:?} must equal ez extent {:?}",
            mask.dim(),
            expected
        );
        self.pec_mask_ez = Some(mask);
        self
    }

    /// Clamp `E` components flagged by the per-component PEC masks to zero.
    ///
    /// Called after each E half-step by the solver. If no masks are
    /// attached this is a no-op; otherwise every cell with `mask[i, j, k]
    /// == true` has the corresponding component reset to zero, modelling
    /// a perfect electric conductor occupying that edge.
    ///
    /// This is *additive* to (and orthogonal from) the outer-face
    /// [`crate::boundary::apply_pec`] / CPML boundary update — masks
    /// describe interior PEC structures (slot edges, ground planes,
    /// striplines) that the deprecated outer-face clamp cannot reach.
    pub fn apply_pec_mask(&mut self) {
        if let Some(mask) = &self.pec_mask_ex {
            for (e, &m) in self.ex.iter_mut().zip(mask.iter()) {
                if m {
                    *e = 0.0;
                }
            }
        }
        if let Some(mask) = &self.pec_mask_ey {
            for (e, &m) in self.ey.iter_mut().zip(mask.iter()) {
                if m {
                    *e = 0.0;
                }
            }
        }
        if let Some(mask) = &self.pec_mask_ez {
            for (e, &m) in self.ez.iter_mut().zip(mask.iter()) {
                if m {
                    *e = 0.0;
                }
            }
        }
    }
}
