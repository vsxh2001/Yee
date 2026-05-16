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
//! Phase 2 walking skeleton: vacuum only (uniform `eps_r`, `mu_r`). Heterogeneous
//! materials, dispersion, and PML are explicitly out of scope here.

use ndarray::Array3;

use yee_core::units::C0;

/// Yee staggered grid with vacuum (or uniform) material.
///
/// The six field arrays use the standard Taflove staggering. All fields are
/// zero-initialized; sources are injected by mutating cells in place.
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

    /// Relative permittivity (uniform; vacuum = 1.0).
    pub eps_r: f64,
    /// Relative permeability (uniform; vacuum = 1.0).
    pub mu_r: f64,
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
}
