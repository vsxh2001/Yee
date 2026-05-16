//! Auxiliary Differential Equation (ADE) update for dispersive materials.
//!
//! Implements ADE-FDTD updates for the three single-pole dispersion models in
//! [`crate::material::Material`]. References below follow the chapter and
//! equation numbering of Taflove & Hagness, *Computational Electrodynamics*,
//! 3rd ed., chapter 9 ("FDTD modeling of frequency-dependent media").
//!
//! ## Algorithm summary
//!
//! Each Yee step performs:
//!
//! 1. `update_h` (unchanged vacuum kernel — magnetic field has no
//!    dispersion in this walking skeleton).
//! 2. **[`DispersiveState::update_e_with_dispersion`]** which fuses the
//!    standard E-curl update with a per-cell branch on
//!    [`crate::material::Material`]:
//!    - **Vacuum**: textbook E update (`E += (dt/ε₀) · curl_H`).
//!    - **Drude**: Taflove §9.4.3 — auxiliary polarization current
//!      `J_p^{n+1/2} = α · J_p^{n-1/2} + β · E^n`, then
//!      `E^{n+1} = E^n + (dt/(ε∞·ε₀)) · (curl_H − (J_p^{n+1/2}+J_p^{n-1/2})/2)`.
//!      (We use the J average for second-order accuracy; the spec lists the
//!      simpler `−J^{n+1/2}` form which is first-order in the source term.)
//!    - **Lorentz**: Taflove §9.5.2 — three-time-level recursion on the
//!      polarization
//!      `P^{n+1} = α_L · P^n + β_L · P^{n-1} + γ_L · E^n`, then
//!      `E^{n+1} = E^n + (dt/(ε∞·ε₀)) · (curl_H − (P^{n+1} − P^n)/dt)`.
//!      For convenience we store P in the `jp_*` arrays (current) and
//!      `jp_*_prev` (previous step); this differs from the strict J-form
//!      of Taflove eq. 9.34a but is mathematically equivalent and far more
//!      numerically stable (no implicit E coupling).
//!    - **Debye**: Taflove §9.6.2 — first-order recursion on the
//!      polarization
//!      `P^{n+1} = P^n + (dt/τ) · (ε₀·Δε·E^n − P^n)`, then
//!      `E^{n+1} = E^n + (dt/(ε∞·ε₀)) · (curl_H − (P^{n+1} − P^n)/dt)`.
//!
//! ## Coefficient pre-computation
//!
//! Per-cell coefficient evaluation inside the hot E-loop would be wasteful, but
//! the walking skeleton allocates the auxiliary state full-grid (see scope
//! note below) so a per-cell `match` is the simplest correct implementation.
//! Coefficient memoisation and a sparse "dispersive cells only" list are
//! Phase 2.fdtd.3.1 work.
//!
//! ## Walking-skeleton scope (Phase 2.fdtd.3)
//!
//! - Auxiliary fields are full-grid (`[nx+1, ny+1, nz+1]`), not sparse.
//! - Same staggering as `E` (each component lives at the same Yee location
//!   as the corresponding E component); we evaluate material tags at the
//!   primary cell `(i, j, k)` for *all three* components — i.e. all three
//!   components in cell `(i, j, k)` see the same material. Sub-cell or
//!   tensorial materials are out of scope.
//! - Per-step branching on `Material` per cell. Vacuum cells fall through
//!   to the standard curl-H update with zero overhead beyond the branch.

use ndarray::Array3;

use yee_core::units::EPS0;

use crate::grid::YeeGrid;
use crate::material::{Material, MaterialMap};

/// Auxiliary fields for ADE dispersive updates.
///
/// Holds the per-cell polarization currents (Drude) or polarizations
/// (Lorentz, Debye) needed to advance the E field through a dispersive
/// medium. All arrays are shape `[nx+1, ny+1, nz+1]` for the walking
/// skeleton; the staggering of each component matches the corresponding
/// `E` component (E_x lives at `[..nx, .., ..]`, etc.) but the storage
/// is over-sized for index simplicity.
///
/// Field layout:
/// - `jp_{x,y,z}` — current-step auxiliary field
///   - Drude: J_p^{n+1/2} (polarization current density, A/m²)
///   - Lorentz: P^n (polarization, C/m²)
///   - Debye: unused (Debye uses `p_*` only).
/// - `jp_{x,y,z}_prev` — previous-step auxiliary field
///   - Drude: J_p^{n-1/2}
///   - Lorentz: P^{n-1}
///   - Debye: unused.
/// - `p_{x,y,z}` — Debye polarization P^n (C/m²).
#[derive(Debug, Clone)]
pub struct DispersiveState {
    /// Current-step polarization-current (Drude) or polarization (Lorentz),
    /// x component.
    pub jp_x: Array3<f64>,
    /// Current-step polarization-current (Drude) or polarization (Lorentz),
    /// y component.
    pub jp_y: Array3<f64>,
    /// Current-step polarization-current (Drude) or polarization (Lorentz),
    /// z component.
    pub jp_z: Array3<f64>,
    /// Previous-step polarization-current (Drude) or polarization (Lorentz),
    /// x component.
    pub jp_x_prev: Array3<f64>,
    /// Previous-step polarization-current (Drude) or polarization (Lorentz),
    /// y component.
    pub jp_y_prev: Array3<f64>,
    /// Previous-step polarization-current (Drude) or polarization (Lorentz),
    /// z component.
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
    /// All fields start at zero — the medium is assumed unpolarized at
    /// `t = 0`, consistent with a vacuum-initial-condition simulation.
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

    /// Fused ADE + E-field update.
    ///
    /// Replaces the standard [`crate::update::update_e`] call in the
    /// presence of dispersive materials. The H field must already have
    /// been advanced by [`crate::update::update_h`] before calling this.
    ///
    /// Per-cell logic:
    /// - [`Material::Vacuum`] → standard curl-H update on `ε∞ = 1`.
    /// - [`Material::Drude`] → Taflove §9.4.3 ADE.
    /// - [`Material::Lorentz`] → Taflove §9.5.2 polarization ADE.
    /// - [`Material::Debye`] → Taflove §9.6.2 polarization ADE.
    pub fn update_e_with_dispersion(&mut self, grid: &mut YeeGrid, materials: &MaterialMap) {
        let dt = grid.dt;
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;

        // ---------------- E_x: shape [nx, ny+1, nz+1] ----------------
        // Interior j ∈ [1, ny), k ∈ [1, nz); j == 0, ny and k == 0, nz
        // are PEC faces (or CPML face cells managed separately).
        for i in 0..nx {
            for j in 1..ny {
                for k in 1..nz {
                    let dhz_dy = (grid.hz[(i, j, k)] - grid.hz[(i, j - 1, k)]) / dy;
                    let dhy_dz = (grid.hy[(i, j, k)] - grid.hy[(i, j, k - 1)]) / dz;
                    let curl_h = dhz_dy - dhy_dz;
                    let mat = materials.material_at(i, j, k);
                    let e_old = grid.ex[(i, j, k)];
                    let (e_new, jp_new, jp_prev_new, p_new) = ade_step(
                        mat,
                        e_old,
                        curl_h,
                        self.jp_x[(i, j, k)],
                        self.jp_x_prev[(i, j, k)],
                        self.p_x[(i, j, k)],
                        dt,
                    );
                    grid.ex[(i, j, k)] = e_new;
                    self.jp_x_prev[(i, j, k)] = jp_prev_new;
                    self.jp_x[(i, j, k)] = jp_new;
                    self.p_x[(i, j, k)] = p_new;
                }
            }
        }

        // ---------------- E_y: shape [nx+1, ny, nz+1] ----------------
        for i in 1..nx {
            for j in 0..ny {
                for k in 1..nz {
                    let dhx_dz = (grid.hx[(i, j, k)] - grid.hx[(i, j, k - 1)]) / dz;
                    let dhz_dx = (grid.hz[(i, j, k)] - grid.hz[(i - 1, j, k)]) / dx;
                    let curl_h = dhx_dz - dhz_dx;
                    let mat = materials.material_at(i, j, k);
                    let e_old = grid.ey[(i, j, k)];
                    let (e_new, jp_new, jp_prev_new, p_new) = ade_step(
                        mat,
                        e_old,
                        curl_h,
                        self.jp_y[(i, j, k)],
                        self.jp_y_prev[(i, j, k)],
                        self.p_y[(i, j, k)],
                        dt,
                    );
                    grid.ey[(i, j, k)] = e_new;
                    self.jp_y_prev[(i, j, k)] = jp_prev_new;
                    self.jp_y[(i, j, k)] = jp_new;
                    self.p_y[(i, j, k)] = p_new;
                }
            }
        }

        // ---------------- E_z: shape [nx+1, ny+1, nz] ----------------
        for i in 1..nx {
            for j in 1..ny {
                for k in 0..nz {
                    let dhy_dx = (grid.hy[(i, j, k)] - grid.hy[(i - 1, j, k)]) / dx;
                    let dhx_dy = (grid.hx[(i, j, k)] - grid.hx[(i, j - 1, k)]) / dy;
                    let curl_h = dhy_dx - dhx_dy;
                    let mat = materials.material_at(i, j, k);
                    let e_old = grid.ez[(i, j, k)];
                    let (e_new, jp_new, jp_prev_new, p_new) = ade_step(
                        mat,
                        e_old,
                        curl_h,
                        self.jp_z[(i, j, k)],
                        self.jp_z_prev[(i, j, k)],
                        self.p_z[(i, j, k)],
                        dt,
                    );
                    grid.ez[(i, j, k)] = e_new;
                    self.jp_z_prev[(i, j, k)] = jp_prev_new;
                    self.jp_z[(i, j, k)] = jp_new;
                    self.p_z[(i, j, k)] = p_new;
                }
            }
        }
    }
}

/// Per-cell ADE step: takes the old E sample, the curl-H, and the auxiliary
/// state, returns the new E and updated auxiliary state.
///
/// Returns `(e_new, jp_new, jp_prev_new, p_new)`:
/// - `jp_prev_new` is what gets stored as `jp_*_prev` after this step
///   (= old `jp_*` for Lorentz; old `jp_*` for Drude too — the rolling shift).
/// - `jp_new` is what gets stored as `jp_*` (the freshly-computed value).
/// - `p_new` is the new Debye polarization (untouched for other models).
///
/// All four return values are unused for `Material::Vacuum` (returns the
/// same auxiliaries unchanged).
#[inline]
fn ade_step(
    mat: Material,
    e_old: f64,
    curl_h: f64,
    jp_old: f64,
    jp_prev_old: f64,
    p_old: f64,
    dt: f64,
) -> (f64, f64, f64, f64) {
    match mat {
        Material::Vacuum => {
            let e_new = e_old + (dt / EPS0) * curl_h;
            (e_new, jp_old, jp_prev_old, p_old)
        }
        Material::Drude {
            eps_inf,
            omega_p,
            gamma,
        } => {
            // Taflove §9.4.3. Update J at the new half-step:
            //   J^{n+1/2} = α J^{n-1/2} + β E^n
            // where α and β below; then update E using the J average for
            // second-order accuracy.
            let half = 0.5 * gamma * dt;
            let alpha = (1.0 - half) / (1.0 + half);
            let beta = (EPS0 * omega_p * omega_p * dt) / (1.0 + half);
            let jp_new = alpha * jp_old + beta * e_old;
            let j_avg = 0.5 * (jp_new + jp_old);
            let e_new = e_old + (dt / (eps_inf * EPS0)) * (curl_h - j_avg);
            // Roll: previous slot gets the value that was current going in.
            (e_new, jp_new, jp_old, p_old)
        }
        Material::Lorentz {
            eps_inf,
            delta_eps,
            omega_0,
            delta,
        } => {
            // Taflove §9.5.2 polarization-form ADE:
            //   P^{n+1} = α_L · P^n + β_L · P^{n-1} + γ_L · E^n
            // with
            //   α_L = (2 − ω₀²·dt²) / (1 + δ·dt)
            //   β_L = (δ·dt − 1)   / (1 + δ·dt)
            //   γ_L = (ε₀ · Δε · ω₀² · dt²) / (1 + δ·dt)
            let denom = 1.0 + delta * dt;
            let alpha_l = (2.0 - omega_0 * omega_0 * dt * dt) / denom;
            let beta_l = (delta * dt - 1.0) / denom;
            let gamma_l = (EPS0 * delta_eps * omega_0 * omega_0 * dt * dt) / denom;
            let p_n = jp_old; // we stash P in the jp_* slot for Lorentz.
            let p_nm1 = jp_prev_old;
            let p_np1 = alpha_l * p_n + beta_l * p_nm1 + gamma_l * e_old;
            // E^{n+1} = E^n + (dt/(ε∞ ε₀)) · (curl_H − (P^{n+1} − P^n)/dt)
            //        = E^n + (dt/(ε∞ ε₀)) · curl_H − (1/(ε∞ ε₀))(P^{n+1} − P^n).
            let e_new = e_old + (dt / (eps_inf * EPS0)) * curl_h - (p_np1 - p_n) / (eps_inf * EPS0);
            (e_new, p_np1, p_n, p_old)
        }
        Material::Debye {
            eps_inf,
            delta_eps,
            tau,
        } => {
            // Taflove §9.6.2 first-order ADE on P:
            //   P^{n+1} = P^n + (dt/τ) · (ε₀ · Δε · E^n − P^n)
            // and
            //   E^{n+1} = E^n + (dt/(ε∞ ε₀)) · curl_H − (P^{n+1} − P^n)/(ε∞ ε₀)
            let p_new = p_old + (dt / tau) * (EPS0 * delta_eps * e_old - p_old);
            let e_new =
                e_old + (dt / (eps_inf * EPS0)) * curl_h - (p_new - p_old) / (eps_inf * EPS0);
            (e_new, jp_old, jp_prev_old, p_new)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Vacuum cells in a DispersiveState should reproduce the standard
    /// `update::update_e` output exactly.
    #[test]
    fn vacuum_cells_match_standard_update() {
        use crate::update;

        let n = 8;
        let dx = 1.0e-3;

        let mut grid_a = YeeGrid::vacuum(n, n, n, dx);
        let mut grid_b = grid_a.clone();

        // Seed an arbitrary H field so the curl is non-zero.
        for i in 0..grid_a.hz.shape()[0] {
            for j in 0..grid_a.hz.shape()[1] {
                for k in 0..grid_a.hz.shape()[2] {
                    let v = 1e-3 * ((i + 2 * j + 3 * k) as f64).sin();
                    grid_a.hz[(i, j, k)] = v;
                    grid_b.hz[(i, j, k)] = v;
                }
            }
        }

        // A: standard vacuum update.
        update::update_e(&mut grid_a);

        // B: vacuum MaterialMap + DispersiveState path.
        let materials = MaterialMap::vacuum(n, n, n);
        let mut state = DispersiveState::new(&materials);
        state.update_e_with_dispersion(&mut grid_b, &materials);

        // Compare the centre interior cell of each E array.
        let centre = (n / 2, n / 2, n / 2);
        assert!(
            (grid_a.ex[centre] - grid_b.ex[centre]).abs() < 1e-15,
            "Ex mismatch: vacuum={} dispersive={}",
            grid_a.ex[centre],
            grid_b.ex[centre]
        );
        assert!(
            (grid_a.ey[centre] - grid_b.ey[centre]).abs() < 1e-15,
            "Ey mismatch"
        );
        assert!(
            (grid_a.ez[centre] - grid_b.ez[centre]).abs() < 1e-15,
            "Ez mismatch"
        );
    }

    /// Drude ADE: with γ = 0, ω_p > 0, and zero curl, the J update should
    /// behave like a lossless harmonic oscillator driven by E.
    /// Specifically, with E held constant and ω_p > 0, J should grow.
    #[test]
    fn drude_ade_drives_current_from_e() {
        let n = 4;
        let dx = 1.0e-3;
        let mut grid = YeeGrid::vacuum(n, n, n, dx);
        let mut materials = MaterialMap::vacuum(n, n, n);
        let drude = Material::Drude {
            eps_inf: 1.0,
            omega_p: 2.0 * std::f64::consts::PI * 1e10,
            gamma: 0.0,
        };
        materials.set_box(0, n + 1, 0, n + 1, 0, n + 1, drude);
        let mut state = DispersiveState::new(&materials);

        // Seed E_z with a constant value at an interior cell.
        let probe = (2, 2, 2);
        grid.ez[probe] = 1.0;

        // First step: curl_H = 0 (H is zero), so E should evolve only via
        // the polarization current. J_p starts at 0, so the very first step
        // gives β·E for J_p^{1/2}, and E updates by −dt/(ε∞ε₀) · (β·E/2).
        state.update_e_with_dispersion(&mut grid, &materials);
        assert!(
            state.jp_z[probe] > 0.0,
            "Drude J_p should be positive after one step with positive E, got {}",
            state.jp_z[probe]
        );
        assert!(
            grid.ez[probe] < 1.0,
            "Drude E should have dropped from seeded value, got {}",
            grid.ez[probe]
        );
    }

    /// Debye ADE: with δε > 0 and E held briefly, P should grow towards
    /// `ε₀ · Δε · E` exponentially.
    #[test]
    fn debye_polarization_grows_toward_steady_state() {
        let n = 4;
        let dx = 1.0e-3;
        let mut grid = YeeGrid::vacuum(n, n, n, dx);
        let mut materials = MaterialMap::vacuum(n, n, n);
        let debye = Material::Debye {
            eps_inf: 1.0,
            delta_eps: 10.0,
            tau: 10.0 * grid.dt, // many time steps to relax.
        };
        materials.set_box(0, n + 1, 0, n + 1, 0, n + 1, debye);
        let mut state = DispersiveState::new(&materials);

        let probe = (2, 2, 2);
        grid.ez[probe] = 1.0;
        state.update_e_with_dispersion(&mut grid, &materials);
        assert!(state.p_z[probe] > 0.0, "Debye P should be positive");
        // Step size dt/τ < 1, so P_new < ε₀·Δε·E_old.
        let target = EPS0 * 10.0 * 1.0;
        assert!(
            state.p_z[probe] < target,
            "Debye P should be below steady state target"
        );
    }

    /// Lorentz ADE sanity: with E held constant, P should grow from 0 in
    /// the first step proportional to γ_L · E.
    #[test]
    fn lorentz_polarization_kicks_from_e() {
        let n = 4;
        let dx = 1.0e-3;
        let mut grid = YeeGrid::vacuum(n, n, n, dx);
        let mut materials = MaterialMap::vacuum(n, n, n);
        let lorentz = Material::Lorentz {
            eps_inf: 1.0,
            delta_eps: 1.0,
            omega_0: 2.0 * std::f64::consts::PI * 1e10,
            delta: 1.0e8,
        };
        materials.set_box(0, n + 1, 0, n + 1, 0, n + 1, lorentz);
        let mut state = DispersiveState::new(&materials);

        let probe = (2, 2, 2);
        grid.ez[probe] = 1.0;
        state.update_e_with_dispersion(&mut grid, &materials);
        assert!(
            state.jp_z[probe] > 0.0,
            "Lorentz P should be positive after one step with positive E, got {}",
            state.jp_z[probe]
        );
    }
}
