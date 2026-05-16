//! Near-to-far-field (NTFF) transformation for FDTD radiation patterns.
//!
//! Implements the standard Stratton–Chu surface-equivalence transformation
//! used in FDTD to obtain a frequency-domain far-field amplitude from a
//! finite-domain time-marching solve. References:
//!
//! > A. Taflove and S. C. Hagness, *Computational Electrodynamics: The
//! > Finite-Difference Time-Domain Method*, 3rd ed., Artech House, 2005,
//! > Chapter 8.
//! >
//! > K. S. Yee, "Time-domain near-to-far-field transformation for radar
//! > cross-section computation", *IEEE Trans. Antennas Propag.* **39**(2),
//! > 1992.
//!
//! ## What this module does
//!
//! 1. Closes a **box-shaped integration surface** inside the FDTD domain,
//!    sitting `box_margin_cells` cells inside every outer face. The
//!    surface is the union of six planar patches normal to ±x̂, ±ŷ, ±ẑ.
//! 2. On every time step the user calls [`NtffState::sample`], which
//!    reads tangential E and H field components on each face,
//!    interpolates them from the Yee-staggered locations to a common
//!    face-centred sample point, computes the equivalent surface
//!    currents `J = n̂ × H`, `M = −n̂ × E`, and accumulates one
//!    **discrete-time Fourier transform bin** at `f_probe`:
//!
//!    ```text
//!    Ĵ(ω) += J(t) · exp(−jωt) · Δt
//!    M̂(ω) += M(t) · exp(−jωt) · Δt
//!    ```
//! 3. After the time loop, [`NtffState::far_field`] projects the
//!    accumulated currents to one (θ, φ) direction using the far-zone
//!    Stratton–Chu kernel:
//!
//!    ```text
//!    E_far(r̂) ∝ (jk / 4π) · r̂ × [ L(r̂) + η₀ · (r̂ × N(r̂)) ]
//!    N(r̂)     = ∫∫_S J(r') · e^{+jk r̂·r'} dS
//!    L(r̂)     = ∫∫_S M(r') · e^{+jk r̂·r'} dS
//!    ```
//!
//!    (Taflove eq. 8.35 with the standard far-zone simplification; the
//!    distance-independent `e^{−jkr}/r` envelope is dropped — re-apply
//!    it externally for a specific observation radius).
//!
//! ## Out of scope for Phase 2.fdtd.2
//!
//! - Multi-frequency probe (one `f_probe` only).
//! - Full θ/φ sweep (single observation direction;
//!   [`NtffState::far_field_at`] does allow extra directions after one
//!   solve so callers can sweep cheaply post-hoc).
//! - Stored-then-projected time-domain currents (we use the running-DFT
//!   approach, which is the standard FDTD trick for narrow-band patterns).

use std::f64::consts::{PI, TAU};

use ndarray::Array2;
use num_complex::Complex64;

use yee_core::units::{C0, ETA0};

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
    /// measured in cells. Typical value `npml + 3` to keep the surface
    /// clear of PML evanescent fields.
    pub box_margin_cells: usize,
    /// Polar angle θ (rad) of the single observation direction returned
    /// by [`NtffState::far_field`].
    pub theta_rad: f64,
    /// Azimuthal angle φ (rad) of the single observation direction.
    pub phi_rad: f64,
}

/// Outward-pointing face index for the box integration surface.
///
/// Order is `[0] = −x̂, [1] = +x̂, [2] = −ŷ, [3] = +ŷ, [4] = −ẑ, [5] = +ẑ`.
const NUM_FACES: usize = 6;

/// Outward face normal as a unit vector.
const FACE_NORMALS: [[f64; 3]; NUM_FACES] = [
    [-1.0, 0.0, 0.0], // −x
    [1.0, 0.0, 0.0],  // +x
    [0.0, -1.0, 0.0], // −y
    [0.0, 1.0, 0.0],  // +y
    [0.0, 0.0, -1.0], // −z
    [0.0, 0.0, 1.0],  // +z
];

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
/// current arrays (electric `J` and magnetic `M`). The DFT runs over the
/// caller's FDTD time loop via [`NtffState::sample`].
pub struct NtffState {
    /// Cached parameters.
    params: NtffParams,
    /// Cached grid geometry.
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
/// integer-E cells along each tangential axis and the fixed normal-axis
/// plane index.
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

/// Normal axis index (0 = x, 1 = y, 2 = z) for face `f`.
#[inline]
fn normal_axis(f: usize) -> usize {
    match f {
        0 | 1 => 0,
        2 | 3 => 1,
        4 | 5 => 2,
        _ => unreachable!(),
    }
}

impl NtffState {
    /// Build a fresh NTFF accumulator sized to `grid` with parameters
    /// `params`.
    ///
    /// The integration surface sits at `params.box_margin_cells` cells
    /// inside every outer face. If using CPML, choose
    /// `box_margin_cells ≥ npml + 1` so the surface stays clear of the
    /// PML.
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
    /// Call once per FDTD time step, after the matching `update_e` and
    /// `update_h` calls (see
    /// [`crate::WalkingSkeletonSolver::step_with_source_and_ntff`]).
    pub fn sample(&mut self, grid: &YeeGrid, t: f64) {
        let omega = TAU * self.params.f_probe;
        // DFT kernel × dt (so the accumulator approximates the
        // continuous-time Fourier integral after the time loop).
        //
        // E and H are half-step staggered in time; for a narrow-band
        // probe the half-step shift only multiplies the result by
        // exp(±jωΔt/2). For Phase 2.fdtd.2 walking skeleton we accept
        // that small shared phase factor; it does not affect magnitudes
        // or pattern ratios.
        let phase = Complex64::from_polar(1.0, -omega * t) * self.dt;

        for f in 0..NUM_FACES {
            self.sample_face(grid, f, phase);
        }

        self.n_steps += 1;
    }

    /// Sample one face. Split out from [`Self::sample`] for clarity.
    fn sample_face(&mut self, grid: &YeeGrid, f: usize, phase: Complex64) {
        let (u_axis, v_axis, n_u, n_v, plane) = face_axes(f, self.bounds);
        let normal = FACE_NORMALS[f];
        let n_axis = normal_axis(f);
        let (i_min, _i_max, j_min, _j_max, k_min, _k_max) = self.bounds;
        let mins = [i_min, j_min, k_min];

        for vi in 0..n_v {
            for ui in 0..n_u {
                // Construct the (i, j, k) integer index of the sample
                // point on this face.
                let mut ijk = [0usize; 3];
                ijk[n_axis] = plane;
                ijk[u_axis] = mins[u_axis] + ui;
                ijk[v_axis] = mins[v_axis] + vi;
                let (i, j, k) = (ijk[0], ijk[1], ijk[2]);

                // Sample E and H at this face point. Simple averaging
                // brings the staggered Yee components to the integer-E
                // node.
                let (ex, ey, ez) = sample_e_at(grid, i, j, k);
                let (hx, hy, hz) = sample_h_at(grid, i, j, k);

                let e = [ex, ey, ez];
                let h = [hx, hy, hz];

                // J = n̂ × H, M = −n̂ × E.
                let j_cross = cross(normal, h);
                let m_cross_pos = cross(normal, e);
                let m_cross = [-m_cross_pos[0], -m_cross_pos[1], -m_cross_pos[2]];

                // Project onto the two tangential axes.
                let j_u = j_cross[u_axis];
                let j_v = j_cross[v_axis];
                let m_u = m_cross[u_axis];
                let m_v = m_cross[v_axis];

                let flat = vi * n_u + ui;
                self.j_face[f][(flat, 0)] += phase * j_u;
                self.j_face[f][(flat, 1)] += phase * j_v;
                self.m_face[f][(flat, 0)] += phase * m_u;
                self.m_face[f][(flat, 1)] += phase * m_v;
            }
        }
    }

    /// Project the accumulated currents to the single observation
    /// direction `(params.theta_rad, params.phi_rad)` using the
    /// far-zone Stratton–Chu kernel. Returns the complex E-field
    /// pattern amplitude.
    pub fn far_field(&self) -> Complex64 {
        self.far_field_at(self.params.theta_rad, self.params.phi_rad)
    }

    /// Project the accumulated currents to an arbitrary observation
    /// direction `(theta, phi)`. Useful for sweeping after one solve.
    ///
    /// Computes the far-zone radiation integrals (Taflove eq. 8.34):
    ///
    /// ```text
    /// N(r̂) = ∫∫_S J(r') · e^{+jk r̂·r'} dS
    /// L(r̂) = ∫∫_S M(r') · e^{+jk r̂·r'} dS
    /// ```
    ///
    /// and forms the far-field electric vector (Taflove eq. 8.35a/b
    /// with the `e^{−jkr}/r` envelope dropped):
    ///
    /// ```text
    /// E_far = (jk / 4π) · r̂ × ( L + η₀ · (r̂ × N) )
    /// ```
    ///
    /// The returned `Complex64` is scalar with magnitude
    /// `|E_far| = √(|E_far,x|² + |E_far,y|² + |E_far,z|²)` and phase
    /// taken from the dominant Cartesian component (largest by
    /// magnitude). For pattern-ratio tests only the magnitude matters,
    /// but exposing a meaningful phase lets future code do coherent
    /// post-processing.
    pub fn far_field_at(&self, theta: f64, phi: f64) -> Complex64 {
        let k_wave = TAU * self.params.f_probe / C0;
        let r_hat = [
            theta.sin() * phi.cos(),
            theta.sin() * phi.sin(),
            theta.cos(),
        ];

        // Accumulate N and L over all six faces.
        let mut n_vec = [Complex64::new(0.0, 0.0); 3];
        let mut l_vec = [Complex64::new(0.0, 0.0); 3];

        for f in 0..NUM_FACES {
            let (u_axis, v_axis, n_u, n_v, plane) = face_axes(f, self.bounds);
            let n_axis = normal_axis(f);
            let (i_min, _i_max, j_min, _j_max, k_min, _k_max) = self.bounds;
            let mins = [i_min, j_min, k_min];
            let ds = self.face_ds(f);

            for vi in 0..n_v {
                for ui in 0..n_u {
                    let mut ijk = [0usize; 3];
                    ijk[n_axis] = plane;
                    ijk[u_axis] = mins[u_axis] + ui;
                    ijk[v_axis] = mins[v_axis] + vi;
                    let (i, j, k) = (ijk[0], ijk[1], ijk[2]);

                    // Position of this sample point in metres relative
                    // to the grid origin. Only differences in the
                    // phase factor between faces matter for the
                    // pattern shape; an overall e^{+jk r̂·r₀}
                    // multiplies the whole result uniformly.
                    let r_prime = [
                        (i as f64) * self.dx,
                        (j as f64) * self.dy,
                        (k as f64) * self.dz,
                    ];
                    let r_dot = r_hat[0] * r_prime[0]
                        + r_hat[1] * r_prime[1]
                        + r_hat[2] * r_prime[2];
                    let phase = Complex64::from_polar(1.0, k_wave * r_dot);

                    let flat = vi * n_u + ui;
                    let j_u = self.j_face[f][(flat, 0)];
                    let j_v = self.j_face[f][(flat, 1)];
                    let m_u = self.m_face[f][(flat, 0)];
                    let m_v = self.m_face[f][(flat, 1)];

                    // Expand the (u, v) tangential J and M back into
                    // 3-D Cartesian components for the integrals.
                    let mut j_xyz = [Complex64::new(0.0, 0.0); 3];
                    let mut m_xyz = [Complex64::new(0.0, 0.0); 3];
                    j_xyz[u_axis] = j_u;
                    j_xyz[v_axis] = j_v;
                    m_xyz[u_axis] = m_u;
                    m_xyz[v_axis] = m_v;

                    let weight = phase * ds;
                    for a in 0..3 {
                        n_vec[a] += weight * j_xyz[a];
                        l_vec[a] += weight * m_xyz[a];
                    }
                }
            }
        }

        // E_far = (jk / 4π) · r̂ × ( L + η₀ · (r̂ × N) ).
        let r_hat_c = [
            Complex64::new(r_hat[0], 0.0),
            Complex64::new(r_hat[1], 0.0),
            Complex64::new(r_hat[2], 0.0),
        ];
        let eta = Complex64::new(ETA0, 0.0);
        let rxn = cross_c(r_hat_c, n_vec);
        let inner = [
            l_vec[0] + eta * rxn[0],
            l_vec[1] + eta * rxn[1],
            l_vec[2] + eta * rxn[2],
        ];
        let outer = cross_c(r_hat_c, inner);
        let prefactor = Complex64::new(0.0, k_wave / (4.0 * PI));

        let e_far_vec = [
            prefactor * outer[0],
            prefactor * outer[1],
            prefactor * outer[2],
        ];

        let mag_sq =
            e_far_vec[0].norm_sqr() + e_far_vec[1].norm_sqr() + e_far_vec[2].norm_sqr();
        if mag_sq == 0.0 {
            return Complex64::new(0.0, 0.0);
        }
        // Pick the dominant Cartesian component for the returned phase.
        let (rep_idx, _) = e_far_vec
            .iter()
            .enumerate()
            .map(|(i, c)| (i, c.norm_sqr()))
            .fold((0usize, 0.0f64), |(bi, bm), (i, m)| {
                if m > bm { (i, m) } else { (bi, bm) }
            });
        let rep = e_far_vec[rep_idx];
        // Scale `rep` so its magnitude equals the full vector
        // magnitude — this gives a single Complex64 whose |.| is the
        // correct |E_far|.
        let scale = mag_sq.sqrt() / rep.norm();
        rep * scale
    }

    /// Area of a single face cell (dA) for face `f`.
    fn face_ds(&self, f: usize) -> f64 {
        match f {
            0 | 1 => self.dy * self.dz, // ±x face
            2 | 3 => self.dx * self.dz, // ±y face
            4 | 5 => self.dx * self.dy, // ±z face
            _ => unreachable!(),
        }
    }
}

/// Real 3-D vector cross product.
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Complex 3-D vector cross product (used for r̂ × ⟨complex⟩).
#[inline]
fn cross_c(a: [Complex64; 3], b: [Complex64; 3]) -> [Complex64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Interpolate the electric-field vector to the integer-E node at
/// `(i, j, k)`.
///
/// The Yee staggering puts each E component on the edge that points
/// along its own direction. To get a "collocated" estimate at the
/// integer node we average the two edges that share the node along the
/// component's axis. Out-of-range neighbours fall back to the available
/// value (one-sided sample at the boundary).
fn sample_e_at(grid: &YeeGrid, i: usize, j: usize, k: usize) -> (f64, f64, f64) {
    // E_x: shape [nx, ny+1, nz+1]; average E_x[(i-1, j, k)] and E_x[(i, j, k)].
    let ex = average_pair_axis(&grid.ex, 0, i, j, k);
    // E_y: shape [nx+1, ny, nz+1]; average E_y[(i, j-1, k)] and E_y[(i, j, k)].
    let ey = average_pair_axis(&grid.ey, 1, i, j, k);
    // E_z: shape [nx+1, ny+1, nz]; average E_z[(i, j, k-1)] and E_z[(i, j, k)].
    let ez = average_pair_axis(&grid.ez, 2, i, j, k);
    (ex, ey, ez)
}

/// Interpolate the magnetic-field vector to the integer-E node at
/// `(i, j, k)`. H lives on Yee faces (half-cell offsets along *two*
/// axes), so we average up to four neighbouring values to reach the
/// integer-E node.
fn sample_h_at(grid: &YeeGrid, i: usize, j: usize, k: usize) -> (f64, f64, f64) {
    // H_x: shape [nx+1, ny, nz], lives at (i, j+½, k+½); average over
    // j ∈ {j-1, j} and k ∈ {k-1, k}.
    let hx = average_face_two_axes(&grid.hx, [1, 2], i, j, k);
    // H_y: shape [nx, ny+1, nz], lives at (i+½, j, k+½); average over
    // i ∈ {i-1, i} and k ∈ {k-1, k}.
    let hy = average_face_two_axes(&grid.hy, [0, 2], i, j, k);
    // H_z: shape [nx, ny, nz+1], lives at (i+½, j+½, k); average over
    // i ∈ {i-1, i} and j ∈ {j-1, j}.
    let hz = average_face_two_axes(&grid.hz, [0, 1], i, j, k);
    (hx, hy, hz)
}

/// Average two values of an edge-staggered array along `axis`
/// (`axis = 0, 1, 2`). The pair is `index − 1` and `index` on `axis`;
/// other indices are clamped into range so boundary nodes degrade to
/// one-sided samples without panicking.
fn average_pair_axis(
    arr: &ndarray::Array3<f64>,
    axis: usize,
    i: usize,
    j: usize,
    k: usize,
) -> f64 {
    let dims = arr.shape();
    let cap = [dims[0], dims[1], dims[2]];
    let ijk = [i, j, k];
    let lo_axis = ijk[axis].checked_sub(1);
    let hi_axis = Some(ijk[axis]);
    let take = |a: Option<usize>| -> Option<f64> {
        let val = a?;
        if val >= cap[axis] {
            return None;
        }
        let mut idx = ijk;
        idx[axis] = val;
        // Clamp the off-axis indices into range. This only matters at
        // the very edge of the grid and never inside the NTFF box.
        for ax in 0..3 {
            if ax != axis && idx[ax] >= cap[ax] {
                idx[ax] = cap[ax] - 1;
            }
        }
        Some(arr[(idx[0], idx[1], idx[2])])
    };
    match (take(lo_axis), take(hi_axis)) {
        (Some(a), Some(b)) => 0.5 * (a + b),
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => 0.0,
    }
}

/// Average up to four face-staggered values around the integer-E node
/// `(i, j, k)` along the two given axes. The non-listed axis is kept
/// fixed (and clamped if it overruns the array bounds).
fn average_face_two_axes(
    arr: &ndarray::Array3<f64>,
    axes: [usize; 2],
    i: usize,
    j: usize,
    k: usize,
) -> f64 {
    let dims = arr.shape();
    let cap = [dims[0], dims[1], dims[2]];
    let ijk = [i, j, k];

    // The fixed axis is the one not in `axes`.
    let fixed_axis = 3 - axes[0] - axes[1];
    let mut base = ijk;
    if base[fixed_axis] >= cap[fixed_axis] {
        base[fixed_axis] = cap[fixed_axis] - 1;
    }

    let mut sum = 0.0f64;
    let mut count = 0u32;
    let a0_lo = ijk[axes[0]].checked_sub(1);
    let a0_hi = Some(ijk[axes[0]]);
    let a1_lo = ijk[axes[1]].checked_sub(1);
    let a1_hi = Some(ijk[axes[1]]);
    for &p0 in &[a0_lo, a0_hi] {
        let Some(p0) = p0 else { continue };
        if p0 >= cap[axes[0]] {
            continue;
        }
        for &p1 in &[a1_lo, a1_hi] {
            let Some(p1) = p1 else { continue };
            if p1 >= cap[axes[1]] {
                continue;
            }
            let mut idx = base;
            idx[axes[0]] = p0;
            idx[axes[1]] = p1;
            sum += arr[(idx[0], idx[1], idx[2])];
            count += 1;
        }
    }
    if count == 0 { 0.0 } else { sum / count as f64 }
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

    #[test]
    fn cross_basic() {
        assert_eq!(cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), [0.0, 0.0, 1.0]);
        assert_eq!(cross([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]), [1.0, 0.0, 0.0]);
        assert_eq!(cross([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]), [0.0, 1.0, 0.0]);
    }

    #[test]
    fn empty_far_field_is_zero() {
        let grid = YeeGrid::vacuum(40, 40, 40, 1.0e-3);
        let params = NtffParams {
            f_probe: 15.0e9,
            box_margin_cells: 12,
            theta_rad: std::f64::consts::FRAC_PI_2,
            phi_rad: 0.0,
        };
        let state = NtffState::new(&grid, params);
        let e = state.far_field();
        assert_eq!(e.norm(), 0.0);
    }

    #[test]
    fn sample_increments_step_counter() {
        let mut grid = YeeGrid::vacuum(40, 40, 40, 1.0e-3);
        // Drop a non-zero E_z and H_x so the DFT bins are non-trivial.
        grid.ez[(20, 20, 20)] = 1.0;
        grid.hx[(20, 20, 20)] = 1.0;
        let params = NtffParams {
            f_probe: 15.0e9,
            box_margin_cells: 12,
            theta_rad: std::f64::consts::FRAC_PI_2,
            phi_rad: 0.0,
        };
        let mut state = NtffState::new(&grid, params);
        state.sample(&grid, 0.0);
        assert_eq!(state.n_samples(), 1);
        // The single-cell source is at (20,20,20), well *inside* the
        // integration box [(12,28)]^3, so it should not appear on any
        // face. All face currents stay zero.
        for f in 0..6 {
            for v in state.j_face[f].iter() {
                assert_eq!(v.norm(), 0.0);
            }
        }
    }
}
