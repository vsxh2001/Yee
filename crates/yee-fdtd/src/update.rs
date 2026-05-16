//! Scalar Yee FDTD update kernels (vacuum, FP64, single-threaded).
//!
//! Implements the textbook leapfrog update from Taflove & Hagness,
//! *Computational Electrodynamics*, 3rd ed., §3.6. Each kernel walks the
//! interior of its target field array with a triply-nested loop.
//!
//! No SIMD, no `rayon`, no GPU. This module is the **walking skeleton**;
//! optimized kernels (CUDA, multi-threaded CPU) will replace it in later
//! sub-phases.
//!
//! ## Update equations (vacuum)
//!
//! ```text
//! H_x(i,j,k) += (dt / (μ₀ μ_r)) · ((E_y(i,j,k+1) - E_y(i,j,k)) / dz
//!                                 - (E_z(i,j+1,k) - E_z(i,j,k)) / dy)
//! H_y(i,j,k) += (dt / (μ₀ μ_r)) · ((E_z(i+1,j,k) - E_z(i,j,k)) / dx
//!                                 - (E_x(i,j,k+1) - E_x(i,j,k)) / dz)
//! H_z(i,j,k) += (dt / (μ₀ μ_r)) · ((E_x(i,j+1,k) - E_x(i,j,k)) / dy
//!                                 - (E_y(i+1,j,k) - E_y(i,j,k)) / dx)
//!
//! E_x(i,j,k) += (dt / (ε₀ ε_r)) · ((H_z(i,j,k) - H_z(i,j-1,k)) / dy
//!                                 - (H_y(i,j,k) - H_y(i,j,k-1)) / dz)
//! E_y(i,j,k) += (dt / (ε₀ ε_r)) · ((H_x(i,j,k) - H_x(i,j,k-1)) / dz
//!                                 - (H_z(i,j,k) - H_z(i-1,j,k)) / dx)
//! E_z(i,j,k) += (dt / (ε₀ ε_r)) · ((H_y(i,j,k) - H_y(i-1,j,k)) / dx
//!                                 - (H_x(i,j,k) - H_x(i,j-1,k)) / dy)
//! ```
//!
//! E-field updates skip the outer tangential faces; those are reset to zero
//! by [`crate::boundary::apply_pec`] each step.

use yee_core::units::{EPS0, MU0};

use crate::grid::YeeGrid;

/// Advance every magnetic-field component by one time step.
///
/// Reads the current `E` field and writes the new `H` field in place.
/// All six staggered field arrays are walked over their full extent for
/// the curl-of-E that produces them.
pub fn update_h(grid: &mut YeeGrid) {
    let dt = grid.dt;
    let dx = grid.dx;
    let dy = grid.dy;
    let dz = grid.dz;
    let coeff = dt / (MU0 * grid.mu_r);

    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    // ---- H_x: shape [nx+1, ny, nz] ----
    for i in 0..=nx {
        for j in 0..ny {
            for k in 0..nz {
                let dey_dz = (grid.ey[(i, j, k + 1)] - grid.ey[(i, j, k)]) / dz;
                let dez_dy = (grid.ez[(i, j + 1, k)] - grid.ez[(i, j, k)]) / dy;
                grid.hx[(i, j, k)] += coeff * (dey_dz - dez_dy);
            }
        }
    }

    // ---- H_y: shape [nx, ny+1, nz] ----
    for i in 0..nx {
        for j in 0..=ny {
            for k in 0..nz {
                let dez_dx = (grid.ez[(i + 1, j, k)] - grid.ez[(i, j, k)]) / dx;
                let dex_dz = (grid.ex[(i, j, k + 1)] - grid.ex[(i, j, k)]) / dz;
                grid.hy[(i, j, k)] += coeff * (dez_dx - dex_dz);
            }
        }
    }

    // ---- H_z: shape [nx, ny, nz+1] ----
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..=nz {
                let dex_dy = (grid.ex[(i, j + 1, k)] - grid.ex[(i, j, k)]) / dy;
                let dey_dx = (grid.ey[(i + 1, j, k)] - grid.ey[(i, j, k)]) / dx;
                grid.hz[(i, j, k)] += coeff * (dex_dy - dey_dx);
            }
        }
    }
}

/// Advance every electric-field component by one time step.
///
/// Reads the current `H` field and writes the new `E` field in place.
/// Tangential E components on the outer faces are deliberately skipped
/// — those cells are managed by the PEC boundary in
/// [`crate::boundary::apply_pec`].
pub fn update_e(grid: &mut YeeGrid) {
    let dt = grid.dt;
    let dx = grid.dx;
    let dy = grid.dy;
    let dz = grid.dz;
    let coeff = dt / (EPS0 * grid.eps_r);

    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    // ---- E_x: shape [nx, ny+1, nz+1] ----
    // Interior j ∈ [1, ny), k ∈ [1, nz); j == 0, ny and k == 0, nz are PEC faces.
    for i in 0..nx {
        for j in 1..ny {
            for k in 1..nz {
                let dhz_dy = (grid.hz[(i, j, k)] - grid.hz[(i, j - 1, k)]) / dy;
                let dhy_dz = (grid.hy[(i, j, k)] - grid.hy[(i, j, k - 1)]) / dz;
                grid.ex[(i, j, k)] += coeff * (dhz_dy - dhy_dz);
            }
        }
    }

    // ---- E_y: shape [nx+1, ny, nz+1] ----
    // Interior i ∈ [1, nx), k ∈ [1, nz).
    for i in 1..nx {
        for j in 0..ny {
            for k in 1..nz {
                let dhx_dz = (grid.hx[(i, j, k)] - grid.hx[(i, j, k - 1)]) / dz;
                let dhz_dx = (grid.hz[(i, j, k)] - grid.hz[(i - 1, j, k)]) / dx;
                grid.ey[(i, j, k)] += coeff * (dhx_dz - dhz_dx);
            }
        }
    }

    // ---- E_z: shape [nx+1, ny+1, nz] ----
    // Interior i ∈ [1, nx), j ∈ [1, ny).
    for i in 1..nx {
        for j in 1..ny {
            for k in 0..nz {
                let dhy_dx = (grid.hy[(i, j, k)] - grid.hy[(i - 1, j, k)]) / dx;
                let dhx_dy = (grid.hx[(i, j, k)] - grid.hx[(i, j - 1, k)]) / dy;
                grid.ez[(i, j, k)] += coeff * (dhy_dx - dhx_dy);
            }
        }
    }
}
