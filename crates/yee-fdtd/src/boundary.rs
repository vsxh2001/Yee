//! Outer-boundary handling for the walking skeleton.
//!
//! # Status
//!
//! Two boundary primitives live here:
//!
//! - [`apply_pec`] (**deprecated**) — reflecting outer-face PEC clamp.
//!   Replaced by [`crate::cpml::CpmlState`] for production work; kept for
//!   `tests/fdtd_propagation.rs` (closed-cavity regression) and for
//!   diagnostic comparisons against CPML.
//! - [`apply_pec_mask`] — per-component **interior** PEC clamp driven by
//!   the optional `pec_mask_e{x,y,z}` arrays on [`YeeGrid`]. This is the
//!   Phase 2.fdtd.7.z infrastructure that lets fdtd-007 (Maloney-Smith
//!   dielectric-loaded slot) model the ground-plane-with-slot geometry as
//!   an interior PEC sheet rather than an outer-face clamp.
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

/// Apply the per-component interior PEC mask (if any) to the grid's E
/// arrays.
///
/// Thin wrapper around [`YeeGrid::apply_pec_mask`] kept in `boundary` so
/// drivers that call `boundary::apply_pec` can mirror the call shape with
/// the per-cell variant. If no masks are attached this is a no-op.
///
/// Intended call site: **after** the E half-step (after `update_e` and
/// after the CPML auxiliary E update) so the clamp is the final word for
/// the step. Calling it before `update_e` would let the next E update
/// reintroduce non-zero values from the curl-of-H stencil.
pub fn apply_pec_mask(grid: &mut YeeGrid) {
    grid.apply_pec_mask();
}
