//! Outer-boundary handling for the walking skeleton.
//!
//! # Status
//!
//! **Deprecated** in favor of [`crate::cpml::CpmlState`].
//!
//! This module ships a *reflecting* hard PEC boundary: tangential
//! electric-field components on the six outer faces are clamped to zero each
//! step, so any energy that reaches the walls bounces back into the
//! computational domain. It is kept around for the regression test in
//! `tests/fdtd_propagation.rs` (which exercises propagation in a closed
//! cavity) and for diagnostic comparisons against CPML.
//!
//! For an absorbing boundary (correct long-term behavior for open-domain
//! antenna or radar-cross-section runs) use the Roden–Gedney 2000 CPML
//! implementation in [`crate::cpml`].

use crate::grid::YeeGrid;

/// Zero out the tangential `E` field on all six outer faces (perfect electric
/// conductor, PEC).
///
/// On the `x = 0` and `x = nx` faces the tangential components are `E_y` and
/// `E_z`. On the `y = 0` and `y = ny` faces they are `E_x` and `E_z`. On the
/// `z = 0` and `z = nz` faces they are `E_x` and `E_y`. Normal `E` and all
/// `H` components are left untouched.
#[deprecated(
    note = "use CpmlState for production; PEC is reflecting and only suitable for cavities"
)]
pub fn apply_pec(grid: &mut YeeGrid) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    // ----- x = 0 and x = nx faces: clamp E_y and E_z -----
    // E_y shape: [nx+1, ny, nz+1]
    for j in 0..ny {
        for k in 0..=nz {
            grid.ey[(0, j, k)] = 0.0;
            grid.ey[(nx, j, k)] = 0.0;
        }
    }
    // E_z shape: [nx+1, ny+1, nz]
    for j in 0..=ny {
        for k in 0..nz {
            grid.ez[(0, j, k)] = 0.0;
            grid.ez[(nx, j, k)] = 0.0;
        }
    }

    // ----- y = 0 and y = ny faces: clamp E_x and E_z -----
    // E_x shape: [nx, ny+1, nz+1]
    for i in 0..nx {
        for k in 0..=nz {
            grid.ex[(i, 0, k)] = 0.0;
            grid.ex[(i, ny, k)] = 0.0;
        }
    }
    // E_z shape: [nx+1, ny+1, nz]
    for i in 0..=nx {
        for k in 0..nz {
            grid.ez[(i, 0, k)] = 0.0;
            grid.ez[(i, ny, k)] = 0.0;
        }
    }

    // ----- z = 0 and z = nz faces: clamp E_x and E_y -----
    // E_x shape: [nx, ny+1, nz+1]
    for i in 0..nx {
        for j in 0..=ny {
            grid.ex[(i, j, 0)] = 0.0;
            grid.ex[(i, j, nz)] = 0.0;
        }
    }
    // E_y shape: [nx+1, ny, nz+1]
    for i in 0..=nx {
        for j in 0..ny {
            grid.ey[(i, j, 0)] = 0.0;
            grid.ey[(i, j, nz)] = 0.0;
        }
    }
}
