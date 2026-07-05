//! Uniform-vacuum FDTD problem description (E.0 scope).

use yee_core::units::C0;

/// Description of a uniform, lossless FDTD problem on a Yee grid.
///
/// This is the E.0 walking-skeleton contract (ADR-0175): scalar `eps_r` /
/// `mu_r`, σ = 0, PEC outer box. The staggered component shapes follow
/// `yee_fdtd::grid::YeeGrid` exactly so parity is testable index-for-index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FdtdSpec {
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
    /// Time step (seconds). Must be ≤ [`FdtdSpec::courant_limit`].
    pub dt: f64,
    /// Relative permittivity (uniform; vacuum = 1.0).
    pub eps_r: f64,
    /// Relative permeability (uniform; vacuum = 1.0).
    pub mu_r: f64,
}

impl FdtdSpec {
    /// Build a vacuum cubic-cell spec with `dt` at 0.9× the Courant limit,
    /// mirroring `YeeGrid::vacuum`.
    ///
    /// # Panics
    ///
    /// Panics if any dimension is zero or `dx` is non-positive/non-finite.
    pub fn vacuum(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0, "grid dimensions must be > 0");
        assert!(
            dx.is_finite() && dx > 0.0,
            "cell size must be positive and finite"
        );
        let mut spec = Self {
            nx,
            ny,
            nz,
            dx,
            dy: dx,
            dz: dx,
            dt: 0.0,
            eps_r: 1.0,
            mu_r: 1.0,
        };
        spec.dt = 0.9 * spec.courant_limit();
        spec
    }

    /// CFL stability limit `1 / (c₀ √(1/dx² + 1/dy² + 1/dz²))`.
    pub fn courant_limit(&self) -> f64 {
        1.0 / (C0
            * (1.0 / (self.dx * self.dx) + 1.0 / (self.dy * self.dy) + 1.0 / (self.dz * self.dz))
                .sqrt())
    }

    /// Shape of the `E_x` array: `[nx, ny+1, nz+1]`.
    pub fn ex_dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny + 1, self.nz + 1)
    }

    /// Shape of the `E_y` array: `[nx+1, ny, nz+1]`.
    pub fn ey_dims(&self) -> (usize, usize, usize) {
        (self.nx + 1, self.ny, self.nz + 1)
    }

    /// Shape of the `E_z` array: `[nx+1, ny+1, nz]`.
    pub fn ez_dims(&self) -> (usize, usize, usize) {
        (self.nx + 1, self.ny + 1, self.nz)
    }

    /// Shape of the `H_x` array: `[nx+1, ny, nz]`.
    pub fn hx_dims(&self) -> (usize, usize, usize) {
        (self.nx + 1, self.ny, self.nz)
    }

    /// Shape of the `H_y` array: `[nx, ny+1, nz]`.
    pub fn hy_dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny + 1, self.nz)
    }

    /// Shape of the `H_z` array: `[nx, ny, nz+1]`.
    pub fn hz_dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny, self.nz + 1)
    }
}

/// Flat row-major index into an array of shape `dims`, matching `ndarray`'s
/// default (C-order) layout: `(i * dim_j + j) * dim_k + k`.
#[inline]
pub(crate) fn idx3(dims: (usize, usize, usize), i: usize, j: usize, k: usize) -> usize {
    debug_assert!(i < dims.0 && j < dims.1 && k < dims.2);
    (i * dims.1 + j) * dims.2 + k
}

/// Element count of an array of shape `dims`.
#[inline]
pub(crate) fn len3(dims: (usize, usize, usize)) -> usize {
    dims.0 * dims.1 * dims.2
}
