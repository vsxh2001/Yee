//! Outer-boundary handling for the walking skeleton.
//!
//! # Status
//!
//! **This is a *reflecting* hard PEC boundary**, not an absorbing one. Tangential
//! electric-field components on the six outer faces are clamped to zero each
//! step, so any energy that reaches the walls bounces back into the
//! computational domain.
//!
//! That is the wrong long-term behavior for an open-domain antenna or
//! radar-cross-section simulation, where the outer faces must *absorb*
//! outgoing waves. A real **CPML** (convolutional perfectly matched layer,
//! Roden & Gedney 2000) implementation is **Phase 2.1+ work** and is
//! deliberately out of scope for the walking skeleton.
//!
//! For now: pick a domain large enough that the wavefront does not reach the
//! walls during the simulation window, or accept the reflections for
//! resonant / cavity-style problems.

use crate::grid::YeeGrid;

/// Zero out the tangential `E` field on all six outer faces (perfect electric
/// conductor, PEC).
///
/// On the `x = 0` and `x = nx` faces the tangential components are `E_y` and
/// `E_z`. On the `y = 0` and `y = ny` faces they are `E_x` and `E_z`. On the
/// `z = 0` and `z = nz` faces they are `E_x` and `E_y`. Normal `E` and all
/// `H` components are left untouched.
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
