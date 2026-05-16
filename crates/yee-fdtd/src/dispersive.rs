//! Auxiliary fields for the Auxiliary Differential Equation (ADE) method.
//!
//! This module hosts the per-cell auxiliary state needed to integrate
//! dispersive materials in time alongside the Yee `E` and `H` fields. The
//! ADE update kernels themselves arrive in a follow-up commit; here we just
//! allocate the auxiliary state and provide a stub
//! [`DispersiveState::update_e_with_dispersion`] that handles the
//! all-vacuum case (i.e. equivalent to [`crate::update::update_e`]).
//!
//! Field layout (Taflove & Hagness §9.4 conventions):
//!
//! - `jp_{x,y,z}` — current-step polarization current density J_p (A/m²)
//!   for Drude / Lorentz cells; unused for Debye.
//! - `jp_{x,y,z}_prev` — previous-step J_p, needed by the two-time-level
//!   Lorentz recursion.
//! - `p_{x,y,z}` — polarization P (C/m²) used by the Debye ADE.
//!
//! ## Walking-skeleton scope (Phase 2.fdtd.3)
//!
//! - Auxiliary fields are full-grid (`[nx+1, ny+1, nz+1]`), not sparse.
//!   A truly sparse "dispersive cells only" representation is a future
//!   optimisation.
//! - All arrays start zero-initialized (medium unpolarized at `t = 0`).

use ndarray::Array3;

use yee_core::units::EPS0;

use crate::grid::YeeGrid;
use crate::material::MaterialMap;

/// Auxiliary fields for ADE dispersive updates.
///
/// Allocated full-grid; staggering matches the corresponding E component.
#[derive(Debug, Clone)]
pub struct DispersiveState {
    /// Current-step polarization-current J_p (Drude/Lorentz), x component.
    pub jp_x: Array3<f64>,
    /// Current-step polarization-current J_p (Drude/Lorentz), y component.
    pub jp_y: Array3<f64>,
    /// Current-step polarization-current J_p (Drude/Lorentz), z component.
    pub jp_z: Array3<f64>,
    /// Previous-step J_p (needed by Lorentz two-time-level recursion).
    pub jp_x_prev: Array3<f64>,
    /// Previous-step J_p, y component.
    pub jp_y_prev: Array3<f64>,
    /// Previous-step J_p, z component.
    pub jp_z_prev: Array3<f64>,
    /// Debye polarization P_x.
    pub p_x: Array3<f64>,
    /// Debye polarization P_y.
    pub p_y: Array3<f64>,
    /// Debye polarization P_z.
    pub p_z: Array3<f64>,
}

impl DispersiveState {
    /// Allocate auxiliary fields sized to the [`MaterialMap`].
    ///
    /// All fields start zero — the medium is assumed unpolarized at `t = 0`.
    pub fn new(materials: &MaterialMap) -> Self {
        let dim = (materials.nx + 1, materials.ny + 1, materials.nz + 1);
        let z = || Array3::<f64>::zeros(dim);
        Self {
            jp_x: z(),
            jp_y: z(),
            jp_z: z(),
            jp_x_prev: z(),
            jp_y_prev: z(),
            jp_z_prev: z(),
            p_x: z(),
            p_y: z(),
            p_z: z(),
        }
    }

    /// Fused ADE + E-field update (stub: all-vacuum path only).
    ///
    /// This walking-skeleton version performs the standard vacuum `E` update
    /// without any dispersive contribution. The Drude / Lorentz / Debye
    /// branches will be added in a follow-up commit.
    pub fn update_e_with_dispersion(
        &mut self,
        grid: &mut YeeGrid,
        _materials: &MaterialMap,
    ) {
        let dt = grid.dt;
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;
        let coeff = dt / EPS0;

        for i in 0..nx {
            for j in 1..ny {
                for k in 1..nz {
                    let dhz_dy = (grid.hz[(i, j, k)] - grid.hz[(i, j - 1, k)]) / dy;
                    let dhy_dz = (grid.hy[(i, j, k)] - grid.hy[(i, j, k - 1)]) / dz;
                    grid.ex[(i, j, k)] += coeff * (dhz_dy - dhy_dz);
                }
            }
        }
        for i in 1..nx {
            for j in 0..ny {
                for k in 1..nz {
                    let dhx_dz = (grid.hx[(i, j, k)] - grid.hx[(i, j, k - 1)]) / dz;
                    let dhz_dx = (grid.hz[(i, j, k)] - grid.hz[(i - 1, j, k)]) / dx;
                    grid.ey[(i, j, k)] += coeff * (dhx_dz - dhz_dx);
                }
            }
        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispersive_state_alloc_zeroes_all_fields() {
        let materials = MaterialMap::vacuum(4, 4, 4);
        let state = DispersiveState::new(&materials);
        assert_eq!(state.jp_x.dim(), (5, 5, 5));
        assert!(state.jp_x.iter().all(|&v| v == 0.0));
        assert!(state.p_z.iter().all(|&v| v == 0.0));
    }
}
