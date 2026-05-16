//! Near-to-far-field (NTFF) transformation for FDTD radiation patterns.
//!
//! Phase 2.fdtd.2 walking-skeleton stub: this commit only ships the
//! parameter type, the per-face DFT accumulator allocation, and the
//! public API surface; `sample` and `far_field` are filled in by the
//! follow-up commits in this same series.
//!
//! References:
//!
//! > A. Taflove and S. C. Hagness, *Computational Electrodynamics: The
//! > Finite-Difference Time-Domain Method*, 3rd ed., Artech House, 2005,
//! > Chapter 8.
//! >
//! > K. S. Yee, "Time-domain near-to-far-field transformation for radar
//! > cross-section computation", *IEEE Trans. Antennas Propag.* **39**(2),
//! > 1992.

use ndarray::Array2;
use num_complex::Complex64;

use crate::grid::YeeGrid;

/// NTFF user-supplied parameters.
///
/// All angles in radians. `box_margin_cells` is measured **inward** from
/// the outer face: the integration surface sits `box_margin_cells` cells
/// inside the outer face on every side.
#[derive(Debug, Clone, Copy)]
pub struct NtffParams {
    /// Probe frequency (Hz). The DFT accumulates exactly one bin at this
    /// frequency.
    pub f_probe: f64,
    /// Standoff of the integration surface from the outer face,
    /// measured in cells.
    pub box_margin_cells: usize,
    /// Polar angle θ (rad) of the single observation direction.
    pub theta_rad: f64,
    /// Azimuthal angle φ (rad) of the single observation direction.
    pub phi_rad: f64,
}

/// Outward-pointing face index for the box integration surface.
///
/// Order is `[0] = −x̂, [1] = +x̂, [2] = −ŷ, [3] = +ŷ, [4] = −ẑ, [5] = +ẑ`.
const NUM_FACES: usize = 6;

/// Per-face tangential axis pairs `(u_axis, v_axis)` chosen so that
/// `ê_u × ê_v = n̂` (right-handed).
const FACE_TANGENT_AXES: [(usize, usize); NUM_FACES] = [
    (2, 1), // −x: (z, y)
    (1, 2), // +x: (y, z)
    (0, 2), // −y: (x, z)
    (2, 0), // +y: (z, x)
    (1, 0), // −z: (y, x)
    (0, 1), // +z: (x, y)
];

/// Accumulator state for the NTFF DFT.
///
/// Each of the six faces holds two complex-valued tangential surface
/// current arrays (electric `J` and magnetic `M`). The DFT will run over
/// the FDTD time loop in a follow-up commit.
pub struct NtffState {
    /// Cached parameters.
    params: NtffParams,
    /// Cached grid geometry.
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    dy: f64,
    dz: f64,
    dt: f64,
    /// Integer cell index of the box face on each axis: `(i_min, i_max,
    /// j_min, j_max, k_min, k_max)`.
    bounds: (usize, usize, usize, usize, usize, usize),
    /// Number of `sample()` calls so far.
    n_steps: u64,
    /// DFT-accumulated tangential electric current per face. Each entry
    /// is shape `(n_u * n_v, 2)` — two tangential components flattened.
    j_face: [Array2<Complex64>; NUM_FACES],
    /// DFT-accumulated tangential magnetic current per face.
    m_face: [Array2<Complex64>; NUM_FACES],
}

/// Extract the two tangential axes for face `f` along with the number of
/// integer-E cells along each tangential axis.
fn face_axes(
    f: usize,
    bounds: (usize, usize, usize, usize, usize, usize),
) -> (
    /* u_axis */ usize,
    /* v_axis */ usize,
    /* n_u */ usize,
    /* n_v */ usize,
    /* plane index */ usize,
) {
    let (i_min, i_max, j_min, j_max, k_min, k_max) = bounds;
    let plane = match f {
        0 => i_min,
        1 => i_max,
        2 => j_min,
        3 => j_max,
        4 => k_min,
        5 => k_max,
        _ => unreachable!(),
    };
    let (u_axis, v_axis) = FACE_TANGENT_AXES[f];
    let extents = [
        (i_max - i_min) + 1,
        (j_max - j_min) + 1,
        (k_max - k_min) + 1,
    ];
    let n_u = extents[u_axis];
    let n_v = extents[v_axis];
    (u_axis, v_axis, n_u, n_v, plane)
}

impl NtffState {
    /// Build a fresh NTFF accumulator sized to `grid` with parameters
    /// `params`.
    ///
    /// The integration surface sits at `params.box_margin_cells` cells
    /// inside every outer face. If using CPML, choose
    /// `box_margin_cells ≥ npml + 1`.
    ///
    /// # Panics
    ///
    /// Panics if the box would have non-positive extent on any axis.
    pub fn new(grid: &YeeGrid, params: NtffParams) -> Self {
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;
        let m = params.box_margin_cells;

        assert!(
            nx > 2 * m + 1 && ny > 2 * m + 1 && nz > 2 * m + 1,
            "NTFF box margin {m} too large for grid {nx}×{ny}×{nz}",
        );
        assert!(
            params.f_probe.is_finite() && params.f_probe > 0.0,
            "NTFF f_probe must be positive and finite",
        );

        let bounds = (m, nx - m, m, ny - m, m, nz - m);

        let j_face: [Array2<Complex64>; NUM_FACES] = std::array::from_fn(|f| {
            let (_, _, n_u, n_v, _) = face_axes(f, bounds);
            Array2::<Complex64>::zeros((n_u * n_v, 2))
        });
        let m_face: [Array2<Complex64>; NUM_FACES] = std::array::from_fn(|f| {
            let (_, _, n_u, n_v, _) = face_axes(f, bounds);
            Array2::<Complex64>::zeros((n_u * n_v, 2))
        });

        Self {
            params,
            nx,
            ny,
            nz,
            dx: grid.dx,
            dy: grid.dy,
            dz: grid.dz,
            dt: grid.dt,
            bounds,
            n_steps: 0,
            j_face,
            m_face,
        }
    }

    /// Borrow the parameters used to build this accumulator.
    pub fn params(&self) -> &NtffParams {
        &self.params
    }

    /// Number of [`Self::sample`] calls so far.
    pub fn n_samples(&self) -> u64 {
        self.n_steps
    }

    /// Box-face indices `(i_min, i_max, j_min, j_max, k_min, k_max)`.
    pub fn bounds(&self) -> (usize, usize, usize, usize, usize, usize) {
        self.bounds
    }

    /// Sample tangential `E` and `H` on each face of the integration
    /// surface and add the contribution at simulation time `t` to the DFT
    /// bin at `params.f_probe`.
    ///
    /// **Stub for this commit**: implemented in the follow-up
    /// "yee-fdtd: sample() — surface-current sampling on 6 faces"
    /// commit. The signature is locked here so downstream code can
    /// reference it.
    pub fn sample(&mut self, _grid: &YeeGrid, _t: f64) {
        // Filled in by the next commit in this series.
        self.n_steps += 1;
        let _ = (&self.dx, &self.dy, &self.dz, &self.dt, &self.nx, &self.ny, &self.nz);
    }

    /// Project the accumulated currents to the single observation
    /// direction `(params.theta_rad, params.phi_rad)`.
    ///
    /// **Stub for this commit**: implemented in the follow-up
    /// "yee-fdtd: far_field() — Stratton-Chu single-direction projection"
    /// commit.
    pub fn far_field(&self) -> Complex64 {
        // Filled in by the next commit. Touch the buffers so clippy
        // doesn't flag them as unread.
        let _ = (&self.j_face, &self.m_face);
        Complex64::new(0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ntff_state_alloc_sizes_match_box() {
        let grid = YeeGrid::vacuum(40, 40, 40, 1.0e-3);
        let params = NtffParams {
            f_probe: 15.0e9,
            box_margin_cells: 12,
            theta_rad: std::f64::consts::FRAC_PI_2,
            phi_rad: 0.0,
        };
        let state = NtffState::new(&grid, params);
        let (i_min, i_max, j_min, j_max, k_min, k_max) = state.bounds();
        assert_eq!(i_min, 12);
        assert_eq!(i_max, 28);
        assert_eq!(j_min, 12);
        assert_eq!(j_max, 28);
        assert_eq!(k_min, 12);
        assert_eq!(k_max, 28);
    }
}
