//! Convolutional Perfectly Matched Layer (CPML) absorbing boundary.
//!
//! Implements the Roden & Gedney (2000) CPML formulation on all six outer
//! faces of a [`crate::grid::YeeGrid`]. References below are to:
//!
//! > J. A. Roden and S. D. Gedney, "Convolutional PML (CPML): An efficient
//! > FDTD implementation of the CFS-PML for arbitrary media",
//! > *Microwave Opt. Technol. Lett.* **27** (5) 334–339, 2000.
//!
//! ## Theory (R&G 2000, §III)
//!
//! The CFS-PML stretching variables are
//!
//! ```text
//! s_w = κ_w + σ_w / (α_w + jωε₀)        (R&G eq. 3)
//! ```
//!
//! In the time domain this becomes a convolution. The convolution is updated
//! recursively via the auxiliary variable ψ (R&G eq. 14):
//!
//! ```text
//! ψ_E_x(n+1) = b_x · ψ_E_x(n) + c_x · ∂H/∂x|^(n+1/2)
//! ```
//!
//! with coefficients (R&G eq. 25):
//!
//! ```text
//! b_w = exp[-(σ_w/κ_w + α_w) · Δt / ε₀]
//! c_w = (σ_w / (σ_w·κ_w + κ_w²·α_w)) · (b_w − 1)
//! ```
//!
//! The Maxwell curl update is then
//!
//! ```text
//! ∂E_z/∂t|cpml = (1/ε) · (1/κ_x · ∂H_y/∂x − 1/κ_y · ∂H_x/∂y + ψ_Ezx − ψ_Ezy)
//! ```
//!
//! ## Walking-skeleton scope
//!
//! - Symmetric `npml`-layer thickness on all six faces (default 10 cells).
//! - Polynomial grading of order `m = 3` (R&G eq. 17 / Taflove §7.5).
//! - Standard parameter set: `σ_max = -(m+1) · ln(R_0) / (2·η₀·npml·dx)`,
//!   `κ_max = 1`, `α_max = 0.05`, reflection target `R_0 = 1e-6`.
//! - Auxiliary ψ arrays allocated full-grid for simple indexing; refactor to
//!   per-face slabs is a future optimization.

use ndarray::Array3;

use yee_core::units::{EPS0, ETA0, MU0};

use crate::grid::YeeGrid;

/// CPML configuration parameters.
///
/// Defaults follow the standard Roden–Gedney 2000 / Taflove §7.5 recipe:
/// 10-cell layer, third-order polynomial grading, `R_0 = 1e-6`,
/// `κ_max = 1`, `α_max = 0.05`. `σ_max` is set so that the theoretical
/// reflection of a normally-incident wave at the inner PML edge is `R_0`.
#[derive(Debug, Clone, Copy)]
pub struct CpmlParams {
    /// Number of PML layers on each face (symmetric). Standard: 10.
    pub npml: usize,
    /// Polynomial grading order. Standard: 3.
    pub m: i32,
    /// Peak conductivity inside the PML. Populated from `R_0 = 1e-6` by
    /// [`CpmlParams::for_grid`]; this raw field is exposed for advanced
    /// callers who want to override the standard choice.
    pub sigma_max: f64,
    /// Peak coordinate-stretching factor. Standard: 1.0 (no stretching).
    pub kappa_max: f64,
    /// Peak CFS shift parameter. Standard: 0.05. Larger `α` improves
    /// low-frequency / evanescent-wave absorption at the cost of more
    /// reflection of propagating waves.
    pub alpha_max: f64,
}

impl Default for CpmlParams {
    fn default() -> Self {
        // sigma_max is grid-dependent; build a placeholder here and let
        // `for_grid` populate it. Default uses npml=10, m=3, dx=1mm.
        let npml = 10;
        let m = 3i32;
        let dx = 1.0e-3;
        Self {
            npml,
            m,
            sigma_max: sigma_max_optimal(m, npml, dx),
            kappa_max: 1.0,
            alpha_max: 0.05,
        }
    }
}

impl CpmlParams {
    /// Standard parameter set sized to the given grid.
    ///
    /// `σ_max = -(m+1) · ln(R_0) / (2·η₀·npml·dx)` with `R_0 = 1e-6`.
    pub fn for_grid(grid: &YeeGrid, npml: usize) -> Self {
        let m = 3i32;
        Self {
            npml,
            m,
            sigma_max: sigma_max_optimal(m, npml, grid.dx),
            kappa_max: 1.0,
            alpha_max: 0.05,
        }
    }
}

/// `σ_max = -(m+1) · ln(R_0) / (2·η₀·npml·dx)` with the standard
/// `R_0 = 1e-6` reflection target (Taflove eq. 7.66).
fn sigma_max_optimal(m: i32, npml: usize, dx: f64) -> f64 {
    let r0: f64 = 1.0e-6;
    -(f64::from(m) + 1.0) * r0.ln() / (2.0 * ETA0 * (npml as f64) * dx)
}

/// CPML auxiliary state.
///
/// The six ψ_E arrays are indexed in the order
/// `(xy, xz, yx, yz, zx, zy)`. Each array is sized to match the
/// corresponding E (or H) field component, so indexing into the field
/// arrays and the ψ arrays uses the same `(i, j, k)` triple.
///
/// `b`, `c`, `kappa` are 1-D vectors of length `npml`, indexed by depth
/// into the PML (0 = innermost, `npml − 1` = outermost). The same
/// coefficients are used at both faces of each axis (the PML is
/// symmetric).
pub struct CpmlState {
    /// Auxiliary ψ arrays for E-field updates: order `(xy, xz, yx, yz, zx, zy)`.
    ///
    /// Indexing convention: `psi_E_uv[i,j,k]` is the running convolution of
    /// `∂H_w/∂v` (where `w ≠ u, v`) that enters the `E_u` update, sized to
    /// the `E_u` component array. The interior of the grid contains zeros
    /// and is never touched.
    pub psi_e: [Array3<f64>; 6],

    /// Auxiliary ψ arrays for H-field updates, same `(xy, xz, yx, yz, zx, zy)`
    /// ordering, each sized to its respective `H_u` array.
    pub psi_h: [Array3<f64>; 6],

    /// Convolutional decay factor `b` evaluated at integer cell positions
    /// (E-field grading), one entry per PML depth (R&G eq. 25).
    pub b: Vec<f64>,
    /// Convolutional source factor `c` evaluated at integer cell positions
    /// (E-field grading), one entry per PML depth (R&G eq. 25).
    pub c: Vec<f64>,
    /// Coordinate-stretching factor κ evaluated at integer cell positions
    /// (E-field grading), one entry per PML depth (R&G eq. 3).
    pub kappa: Vec<f64>,

    /// Convolutional decay factor `b` evaluated at *half-cell* positions
    /// (H-field grading). The H field is on a half-cell-shifted Yee grid;
    /// using a separate profile shifted by `+0.5` cell improves CPML
    /// performance versus reusing the E profile.
    pub b_h: Vec<f64>,
    /// `c` evaluated at half-cell positions (H-field grading).
    pub c_h: Vec<f64>,
    /// κ evaluated at half-cell positions (H-field grading).
    pub kappa_h: Vec<f64>,

    /// PML thickness in cells (cached from [`CpmlParams::npml`]).
    npml: usize,
}

/// Sample σ, κ, α at depth `rho_over_d` ∈ (0, 1] using the standard
/// polynomial grading (R&G eq. 17 / Taflove eq. 7.79):
///
/// ```text
/// σ(ρ) = σ_max · (ρ/d)^m
/// κ(ρ) = 1 + (κ_max − 1) · (ρ/d)^m
/// α(ρ) = α_max · (1 − ρ/d)
/// ```
fn grading_sample(params: &CpmlParams, rho_over_d: f64) -> (f64, f64, f64) {
    let rho_m = rho_over_d.powi(params.m);
    let sigma = params.sigma_max * rho_m;
    let kappa = 1.0 + (params.kappa_max - 1.0) * rho_m;
    let alpha = params.alpha_max * (1.0 - rho_over_d);
    (sigma, kappa, alpha)
}

/// Finalize `(b, c, κ)` from `(σ, κ, α)` and `dt` via R&G eq. 25:
///
/// ```text
/// b = exp(-(σ/κ + α) · dt / ε₀)
/// c = σ / (σ·κ + κ²·α) · (b − 1)
/// ```
fn finalize_coeffs(sigma: f64, kappa: f64, alpha: f64, dt: f64) -> (f64, f64) {
    let exponent = -(sigma / kappa + alpha) * dt / EPS0;
    let b = exponent.exp();
    let denom = sigma * kappa + kappa * kappa * alpha;
    let c = if denom.abs() > 1.0e-30 {
        sigma * (b - 1.0) / denom
    } else {
        0.0
    };
    (b, c)
}

/// One CPML coefficient profile triple: `(b, c, κ)` each of length `npml`.
type ProfileTriple = (Vec<f64>, Vec<f64>, Vec<f64>);

/// Compute the (b, c, κ) profiles at both E-cell and H-cell positions.
///
/// The PML profile vectors are indexed by `d` ∈ `[0, npml)`, where `d = 0`
/// is the **innermost** PML cell (smallest σ) and `d = npml − 1` is the
/// **outermost** PML cell (largest σ).
///
/// E nodes sit at integer cell centres along the PML axis. The outermost
/// E node coincides with the PEC outer face, so the depth fraction is
/// `ρ/d = (d + 1) / npml` — i.e. `d = npml − 1` gives ρ/d = 1 (maximum σ
/// on the wall) and `d = 0` gives ρ/d = 1/npml (innermost cell, small σ).
///
/// H nodes are half-cell-shifted toward the outer face along the PML
/// axis: the outermost H node sits half a cell *inside* the PEC face. So
/// the H depth fraction is `ρ/d = (d + 0.5) / npml`.
fn make_profiles(params: &CpmlParams, dt: f64) -> (ProfileTriple, ProfileTriple) {
    let n = params.npml;
    let mut b_e = vec![1.0; n];
    let mut c_e = vec![0.0; n];
    let mut kappa_e = vec![1.0; n];
    let mut b_h = vec![1.0; n];
    let mut c_h = vec![0.0; n];
    let mut kappa_h = vec![1.0; n];

    for d in 0..n {
        // E grading: depth fraction (d + 1)/npml.
        let rho_e = (d as f64 + 1.0) / (n as f64);
        let (sigma_e, kappa_ev, alpha_e) = grading_sample(params, rho_e);
        let (b, c) = finalize_coeffs(sigma_e, kappa_ev, alpha_e, dt);
        b_e[d] = b;
        c_e[d] = c;
        kappa_e[d] = kappa_ev;

        // H grading: half-cell shifted inward by -0.5 cell relative to E.
        let rho_h = (d as f64 + 0.5) / (n as f64);
        let (sigma_h, kappa_hv, alpha_h) = grading_sample(params, rho_h);
        let (b, c) = finalize_coeffs(sigma_h, kappa_hv, alpha_h, dt);
        b_h[d] = b;
        c_h[d] = c;
        kappa_h[d] = kappa_hv;
    }

    ((b_e, c_e, kappa_e), (b_h, c_h, kappa_h))
}

impl CpmlState {
    /// Build a fresh CPML state sized to `grid` with parameters `params`.
    ///
    /// The convolutional coefficients `b` and `c` are computed using the
    /// grid's `dt`; if the caller subsequently changes `grid.dt` they must
    /// rebuild the [`CpmlState`].
    pub fn new(grid: &YeeGrid, params: CpmlParams) -> Self {
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;

        // ψ_E arrays: same shape as the corresponding E-field component.
        // Order: xy, xz, yx, yz, zx, zy (Psi::idx).
        let psi_e = [
            Array3::<f64>::zeros((nx, ny + 1, nz + 1)), // E_x, ∂/∂y
            Array3::<f64>::zeros((nx, ny + 1, nz + 1)), // E_x, ∂/∂z
            Array3::<f64>::zeros((nx + 1, ny, nz + 1)), // E_y, ∂/∂x
            Array3::<f64>::zeros((nx + 1, ny, nz + 1)), // E_y, ∂/∂z
            Array3::<f64>::zeros((nx + 1, ny + 1, nz)), // E_z, ∂/∂x
            Array3::<f64>::zeros((nx + 1, ny + 1, nz)), // E_z, ∂/∂y
        ];

        // ψ_H arrays: same shape as the corresponding H-field component.
        // Order: Hxy, Hxz, Hyx, Hyz, Hzx, Hzy (Psi::idx).
        let psi_h = [
            Array3::<f64>::zeros((nx + 1, ny, nz)), // H_x, ∂E_z/∂y
            Array3::<f64>::zeros((nx + 1, ny, nz)), // H_x, ∂E_y/∂z
            Array3::<f64>::zeros((nx, ny + 1, nz)), // H_y, ∂E_z/∂x
            Array3::<f64>::zeros((nx, ny + 1, nz)), // H_y, ∂E_x/∂z
            Array3::<f64>::zeros((nx, ny, nz + 1)), // H_z, ∂E_y/∂x
            Array3::<f64>::zeros((nx, ny, nz + 1)), // H_z, ∂E_x/∂y
        ];

        // Build spatial profiles for both E-grid and H-grid (half-cell-shifted).
        let ((b, c, kappa), (b_h, c_h, kappa_h)) = make_profiles(&params, grid.dt);

        Self {
            psi_e,
            psi_h,
            b,
            c,
            kappa,
            b_h,
            c_h,
            kappa_h,
            npml: params.npml,
        }
    }

    /// PML thickness in cells.
    pub fn npml(&self) -> usize {
        self.npml
    }

    /// Layer-depth lookup: given an absolute grid index `i` on an axis of
    /// length `n`, return the PML profile index in `0..npml` if the cell
    /// is inside the low- or high-side PML, else `None`. Also returns the
    /// "side" (`false` = low/origin side, `true` = high/far side).
    #[inline]
    fn pml_depth(&self, i: usize, n: usize) -> Option<(usize, bool)> {
        if i < self.npml {
            // Low-side PML: i=0 is the outermost cell, i=npml-1 is innermost.
            // Profile index is "depth from inner edge", so depth = npml-1-i.
            Some((self.npml - 1 - i, false))
        } else if i >= n.saturating_sub(self.npml) && n >= self.npml {
            // High-side PML: i = n-npml is innermost, i = n-1 is outermost.
            // depth = i - (n - npml).
            let depth = i - (n - self.npml);
            if depth < self.npml {
                Some((depth, true))
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Apply the CPML correction during the E-field update.
    ///
    /// Call this **instead of** the interior portion of `update_e` for the
    /// PML region, *after* `update_e` has already done the standard
    /// "1/ε · curl H · Δt" interior update. This function adds the CPML
    /// correction term `(Δt/ε) · (ψ + (1/κ − 1) · ∂H/∂ξ)` to every E-field
    /// cell whose curl term came from a derivative *along* a PML axis.
    ///
    /// The convention used here (R&G eq. 14, written for an E-update):
    ///
    /// ```text
    /// ψ^{n+1}_{E,uv} = b_v · ψ^n_{E,uv} + c_v · ∂H_w/∂v|^{n+1/2}
    /// E_u^{n+1} += (Δt/ε) · ψ^{n+1}_{E,uv} − (Δt/ε)·(1−1/κ_v)·∂H_w/∂v
    /// ```
    ///
    /// The `(1 − 1/κ)` term collapses to zero for the standard choice
    /// `κ_max = 1` but is included for completeness.
    ///
    /// When `grid.eps_r_cells` is `Some`, the coefficient `Δt/(ε₀ ε_r)` is
    /// recomputed per cell from the per-cell map (mirroring the convention
    /// used by [`crate::update::update_e`]); otherwise the scalar
    /// `grid.eps_r` is used and the loop body matches the pre-percell
    /// behaviour bit-for-bit.
    pub fn update_e(&mut self, grid: &mut YeeGrid) {
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let dt = grid.dt;
        let coeff_scalar = dt / (EPS0 * grid.eps_r);
        let eps_r_cells = grid.eps_r_cells.as_ref();

        // ---- E_x: shape [nx, ny+1, nz+1] ----
        // ∂H_z/∂y term -> ψ_E_xy[i,j,k]; PML active when j-1 or j is in y-PML.
        // ∂H_y/∂z term -> ψ_E_xz[i,j,k]; PML active when k-1 or k is in z-PML.
        // E_x update region: i ∈ [0, nx), j ∈ [1, ny), k ∈ [1, nz).
        for i in 0..nx {
            for j in 1..ny {
                let dep_y = self.pml_depth(j, ny + 1);
                for k in 1..nz {
                    let dep_z = self.pml_depth(k, nz + 1);
                    if dep_y.is_none() && dep_z.is_none() {
                        continue;
                    }
                    let dhz_dy = (grid.hz[(i, j, k)] - grid.hz[(i, j - 1, k)]) / dy;
                    let dhy_dz = (grid.hy[(i, j, k)] - grid.hy[(i, j, k - 1)]) / dz;
                    let coeff = match eps_r_cells {
                        None => coeff_scalar,
                        Some(e) => dt / (EPS0 * e[(i, j, k)]),
                    };

                    if let Some((d, _)) = dep_y {
                        let b = self.b[d];
                        let c = self.c[d];
                        let k_y = self.kappa[d];
                        let p = b * self.psi_e[Psi::Exy.idx()][(i, j, k)] + c * dhz_dy;
                        self.psi_e[Psi::Exy.idx()][(i, j, k)] = p;
                        // The standard scalar update already added coeff*dhz_dy.
                        // CPML correction:  coeff * (p − (1 − 1/κ) · dhz_dy)
                        grid.ex[(i, j, k)] += coeff * (p - (1.0 - 1.0 / k_y) * dhz_dy);
                    }
                    if let Some((d, _)) = dep_z {
                        let b = self.b[d];
                        let c = self.c[d];
                        let k_z = self.kappa[d];
                        let p = b * self.psi_e[Psi::Exz.idx()][(i, j, k)] + c * dhy_dz;
                        self.psi_e[Psi::Exz.idx()][(i, j, k)] = p;
                        // Sign convention: standard update was coeff*(dhz_dy - dhy_dz).
                        // So the dhy_dz contribution enters with a minus sign;
                        // the CPML correction matches that sign:
                        grid.ex[(i, j, k)] -= coeff * (p - (1.0 - 1.0 / k_z) * dhy_dz);
                    }
                }
            }
        }

        // ---- E_y: shape [nx+1, ny, nz+1] ----
        // ∂H_x/∂z term -> ψ_E_yz, ∂H_z/∂x term -> ψ_E_yx.
        // E_y update region: i ∈ [1, nx), j ∈ [0, ny), k ∈ [1, nz).
        for i in 1..nx {
            let dep_x = self.pml_depth(i, nx + 1);
            for j in 0..ny {
                for k in 1..nz {
                    let dep_z = self.pml_depth(k, nz + 1);
                    if dep_x.is_none() && dep_z.is_none() {
                        continue;
                    }
                    let dhx_dz = (grid.hx[(i, j, k)] - grid.hx[(i, j, k - 1)]) / dz;
                    let dhz_dx = (grid.hz[(i, j, k)] - grid.hz[(i - 1, j, k)]) / dx;
                    let coeff = match eps_r_cells {
                        None => coeff_scalar,
                        Some(e) => dt / (EPS0 * e[(i, j, k)]),
                    };

                    if let Some((d, _)) = dep_z {
                        let b = self.b[d];
                        let c = self.c[d];
                        let k_z = self.kappa[d];
                        let p = b * self.psi_e[Psi::Eyz.idx()][(i, j, k)] + c * dhx_dz;
                        self.psi_e[Psi::Eyz.idx()][(i, j, k)] = p;
                        // Standard E_y curl is +dhx_dz − dhz_dx.
                        grid.ey[(i, j, k)] += coeff * (p - (1.0 - 1.0 / k_z) * dhx_dz);
                    }
                    if let Some((d, _)) = dep_x {
                        let b = self.b[d];
                        let c = self.c[d];
                        let k_x = self.kappa[d];
                        let p = b * self.psi_e[Psi::Eyx.idx()][(i, j, k)] + c * dhz_dx;
                        self.psi_e[Psi::Eyx.idx()][(i, j, k)] = p;
                        grid.ey[(i, j, k)] -= coeff * (p - (1.0 - 1.0 / k_x) * dhz_dx);
                    }
                }
            }
        }

        // ---- E_z: shape [nx+1, ny+1, nz] ----
        // ∂H_y/∂x term -> ψ_E_zx, ∂H_x/∂y term -> ψ_E_zy.
        // E_z update region: i ∈ [1, nx), j ∈ [1, ny), k ∈ [0, nz).
        for i in 1..nx {
            let dep_x = self.pml_depth(i, nx + 1);
            for j in 1..ny {
                let dep_y = self.pml_depth(j, ny + 1);
                if dep_x.is_none() && dep_y.is_none() {
                    continue;
                }
                for k in 0..nz {
                    let dhy_dx = (grid.hy[(i, j, k)] - grid.hy[(i - 1, j, k)]) / dx;
                    let dhx_dy = (grid.hx[(i, j, k)] - grid.hx[(i, j - 1, k)]) / dy;
                    let coeff = match eps_r_cells {
                        None => coeff_scalar,
                        Some(e) => dt / (EPS0 * e[(i, j, k)]),
                    };

                    if let Some((d, _)) = dep_x {
                        let b = self.b[d];
                        let c = self.c[d];
                        let k_x = self.kappa[d];
                        let p = b * self.psi_e[Psi::Ezx.idx()][(i, j, k)] + c * dhy_dx;
                        self.psi_e[Psi::Ezx.idx()][(i, j, k)] = p;
                        // Standard E_z curl is +dhy_dx − dhx_dy.
                        grid.ez[(i, j, k)] += coeff * (p - (1.0 - 1.0 / k_x) * dhy_dx);
                    }
                    if let Some((d, _)) = dep_y {
                        let b = self.b[d];
                        let c = self.c[d];
                        let k_y = self.kappa[d];
                        let p = b * self.psi_e[Psi::Ezy.idx()][(i, j, k)] + c * dhx_dy;
                        self.psi_e[Psi::Ezy.idx()][(i, j, k)] = p;
                        grid.ez[(i, j, k)] -= coeff * (p - (1.0 - 1.0 / k_y) * dhx_dy);
                    }
                }
            }
        }
    }

    /// Apply the CPML correction during the H-field update.
    ///
    /// Same structure as [`Self::update_e`] but for the H curl. The auxiliary
    /// variables in `psi_h` are updated with the appropriate `∂E/∂ξ`
    /// derivative and then added to the H field with the matching sign from
    /// the standard Yee curl.
    ///
    /// When `grid.mu_r_cells` is `Some`, the coefficient `Δt/(μ₀ μ_r)` is
    /// recomputed per cell from the per-cell map (mirroring the convention
    /// used by [`crate::update::update_h`]); otherwise the scalar
    /// `grid.mu_r` is used and the loop body matches the pre-percell
    /// behaviour bit-for-bit.
    pub fn update_h(&mut self, grid: &mut YeeGrid) {
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let dt = grid.dt;
        let coeff_scalar = dt / (MU0 * grid.mu_r);
        let mu_r_cells = grid.mu_r_cells.as_ref();

        // ---- H_x: shape [nx+1, ny, nz] ----
        // H_x sits at (x_int, y+0.5, z+0.5); along y and z it uses the
        // half-cell-shifted H profile.
        // Standard curl: +dey_dz − dez_dy.
        for i in 0..=nx {
            for j in 0..ny {
                let dep_y = self.pml_depth(j, ny);
                for k in 0..nz {
                    let dep_z = self.pml_depth(k, nz);
                    if dep_y.is_none() && dep_z.is_none() {
                        continue;
                    }
                    let dey_dz = (grid.ey[(i, j, k + 1)] - grid.ey[(i, j, k)]) / dz;
                    let dez_dy = (grid.ez[(i, j + 1, k)] - grid.ez[(i, j, k)]) / dy;
                    let coeff = match mu_r_cells {
                        None => coeff_scalar,
                        Some(m) => dt / (MU0 * m[(i, j, k)]),
                    };

                    if let Some((d, _)) = dep_z {
                        let b = self.b_h[d];
                        let c = self.c_h[d];
                        let k_z = self.kappa_h[d];
                        let p = b * self.psi_h[Psi::Hxz.idx()][(i, j, k)] + c * dey_dz;
                        self.psi_h[Psi::Hxz.idx()][(i, j, k)] = p;
                        grid.hx[(i, j, k)] += coeff * (p - (1.0 - 1.0 / k_z) * dey_dz);
                    }
                    if let Some((d, _)) = dep_y {
                        let b = self.b_h[d];
                        let c = self.c_h[d];
                        let k_y = self.kappa_h[d];
                        let p = b * self.psi_h[Psi::Hxy.idx()][(i, j, k)] + c * dez_dy;
                        self.psi_h[Psi::Hxy.idx()][(i, j, k)] = p;
                        grid.hx[(i, j, k)] -= coeff * (p - (1.0 - 1.0 / k_y) * dez_dy);
                    }
                }
            }
        }

        // ---- H_y: shape [nx, ny+1, nz] ----
        // H_y sits at (x+0.5, y_int, z+0.5); along x and z it uses the
        // half-cell-shifted H profile.
        // Standard curl: +dez_dx − dex_dz.
        for i in 0..nx {
            let dep_x = self.pml_depth(i, nx);
            for j in 0..=ny {
                for k in 0..nz {
                    let dep_z = self.pml_depth(k, nz);
                    if dep_x.is_none() && dep_z.is_none() {
                        continue;
                    }
                    let dez_dx = (grid.ez[(i + 1, j, k)] - grid.ez[(i, j, k)]) / dx;
                    let dex_dz = (grid.ex[(i, j, k + 1)] - grid.ex[(i, j, k)]) / dz;
                    let coeff = match mu_r_cells {
                        None => coeff_scalar,
                        Some(m) => dt / (MU0 * m[(i, j, k)]),
                    };

                    if let Some((d, _)) = dep_x {
                        let b = self.b_h[d];
                        let c = self.c_h[d];
                        let k_x = self.kappa_h[d];
                        let p = b * self.psi_h[Psi::Hyx.idx()][(i, j, k)] + c * dez_dx;
                        self.psi_h[Psi::Hyx.idx()][(i, j, k)] = p;
                        grid.hy[(i, j, k)] += coeff * (p - (1.0 - 1.0 / k_x) * dez_dx);
                    }
                    if let Some((d, _)) = dep_z {
                        let b = self.b_h[d];
                        let c = self.c_h[d];
                        let k_z = self.kappa_h[d];
                        let p = b * self.psi_h[Psi::Hyz.idx()][(i, j, k)] + c * dex_dz;
                        self.psi_h[Psi::Hyz.idx()][(i, j, k)] = p;
                        grid.hy[(i, j, k)] -= coeff * (p - (1.0 - 1.0 / k_z) * dex_dz);
                    }
                }
            }
        }

        // ---- H_z: shape [nx, ny, nz+1] ----
        // H_z sits at (x+0.5, y+0.5, z_int); along x and y it uses the
        // half-cell-shifted H profile.
        // Standard curl: +dex_dy − dey_dx.
        for i in 0..nx {
            let dep_x = self.pml_depth(i, nx);
            for j in 0..ny {
                let dep_y = self.pml_depth(j, ny);
                if dep_x.is_none() && dep_y.is_none() {
                    continue;
                }
                for k in 0..=nz {
                    let dex_dy = (grid.ex[(i, j + 1, k)] - grid.ex[(i, j, k)]) / dy;
                    let dey_dx = (grid.ey[(i + 1, j, k)] - grid.ey[(i, j, k)]) / dx;
                    let coeff = match mu_r_cells {
                        None => coeff_scalar,
                        Some(m) => dt / (MU0 * m[(i, j, k)]),
                    };

                    if let Some((d, _)) = dep_y {
                        let b = self.b_h[d];
                        let c = self.c_h[d];
                        let k_y = self.kappa_h[d];
                        let p = b * self.psi_h[Psi::Hzy.idx()][(i, j, k)] + c * dex_dy;
                        self.psi_h[Psi::Hzy.idx()][(i, j, k)] = p;
                        grid.hz[(i, j, k)] += coeff * (p - (1.0 - 1.0 / k_y) * dex_dy);
                    }
                    if let Some((d, _)) = dep_x {
                        let b = self.b_h[d];
                        let c = self.c_h[d];
                        let k_x = self.kappa_h[d];
                        let p = b * self.psi_h[Psi::Hzx.idx()][(i, j, k)] + c * dey_dx;
                        self.psi_h[Psi::Hzx.idx()][(i, j, k)] = p;
                        grid.hz[(i, j, k)] -= coeff * (p - (1.0 - 1.0 / k_x) * dey_dx);
                    }
                }
            }
        }
    }
}

/// Index labels for the six ψ auxiliary arrays.
///
/// The first letter selects the field component the ψ contributes to; the
/// second letter selects the derivative axis whose stretched-coordinate
/// convolution this ψ tracks. The same six indices are reused for `psi_e`
/// and `psi_h` — interpret `Exy` as "Hxy" when accessing `psi_h`.
#[derive(Debug, Clone, Copy)]
enum Psi {
    /// E_x update, ∂/∂y derivative (or H_x, ∂/∂y when used for psi_h).
    Exy = 0,
    /// E_x update, ∂/∂z derivative (or H_x, ∂/∂z).
    Exz = 1,
    /// E_y update, ∂/∂x derivative (or H_y, ∂/∂x).
    Eyx = 2,
    /// E_y update, ∂/∂z derivative (or H_y, ∂/∂z).
    Eyz = 3,
    /// E_z update, ∂/∂x derivative (or H_z, ∂/∂x).
    Ezx = 4,
    /// E_z update, ∂/∂y derivative (or H_z, ∂/∂y).
    Ezy = 5,
}

#[allow(non_upper_case_globals)]
impl Psi {
    /// H_x update, ∂/∂y derivative.
    const Hxy: Self = Self::Exy;
    /// H_x update, ∂/∂z derivative.
    const Hxz: Self = Self::Exz;
    /// H_y update, ∂/∂x derivative.
    const Hyx: Self = Self::Eyx;
    /// H_y update, ∂/∂z derivative.
    const Hyz: Self = Self::Eyz;
    /// H_z update, ∂/∂x derivative.
    const Hzx: Self = Self::Ezx;
    /// H_z update, ∂/∂y derivative.
    const Hzy: Self = Self::Ezy;

    #[inline]
    fn idx(self) -> usize {
        self as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_have_expected_shape() {
        let p = CpmlParams::default();
        assert_eq!(p.npml, 10);
        assert_eq!(p.m, 3);
        assert!(p.sigma_max > 0.0);
        assert_eq!(p.kappa_max, 1.0);
        assert!((p.alpha_max - 0.05).abs() < 1e-12);
    }

    #[test]
    fn sigma_max_optimal_matches_roden_gedney() {
        // For npml=10, dx=1mm, m=3, R_0=1e-6:
        // σ_max = -(4) · ln(1e-6) / (2 · 376.73 · 10 · 1e-3)
        //       ≈ -(4) · (-13.8155) / (7.5346) ≈ 7.334
        let s = sigma_max_optimal(3, 10, 1.0e-3);
        assert!(s > 7.0 && s < 8.0, "σ_max ≈ 7.33 expected, got {s}");
    }

    #[test]
    fn b_coeffs_are_in_unit_interval() {
        let grid = YeeGrid::vacuum(20, 20, 20, 1.0e-3);
        let p = CpmlParams::for_grid(&grid, 10);
        let state = CpmlState::new(&grid, p);
        for &b in &state.b {
            assert!(b > 0.0 && b <= 1.0, "b out of (0,1]: {b}");
        }
        // c is negative for σ > 0 (since b < 1 ⇒ b-1 < 0).
        for &c in &state.c {
            assert!(c <= 0.0, "c should be ≤ 0, got {c}");
        }
    }

    #[test]
    fn pml_depth_lookup_low_high_sides() {
        let grid = YeeGrid::vacuum(50, 50, 50, 1.0e-3);
        let state = CpmlState::new(&grid, CpmlParams::for_grid(&grid, 10));
        // n = 50: low PML covers i ∈ [0, 10), high PML covers i ∈ [40, 50).
        assert_eq!(state.pml_depth(0, 50), Some((9, false))); // outermost on low side
        assert_eq!(state.pml_depth(9, 50), Some((0, false))); // innermost on low side
        assert_eq!(state.pml_depth(10, 50), None); // interior
        assert_eq!(state.pml_depth(39, 50), None); // interior
        assert_eq!(state.pml_depth(40, 50), Some((0, true))); // innermost on high side
        assert_eq!(state.pml_depth(49, 50), Some((9, true))); // outermost on high side
    }
}
