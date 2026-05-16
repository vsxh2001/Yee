//! Per-cell dispersive material map.
//!
//! Phase 2.fdtd.3 introduces three single-pole dispersive material models in
//! addition to the default lossless vacuum. Conventions follow Taflove &
//! Hagness, *Computational Electrodynamics*, 3rd ed., §9.4–9.6:
//!
//! - **Drude**:    `eps(ω) = eps_inf − ω_p² / (ω² − j γ ω)`  (Taflove eq. 9.7)
//! - **Lorentz**:  `eps(ω) = eps_inf + Δε · ω₀² / (ω₀² − ω² + 2 j δ ω)`
//!   (Taflove eq. 9.8)
//! - **Debye**:    `eps(ω) = eps_inf + Δε / (1 + j ω τ)`  (Taflove eq. 9.6)
//!
//! These are the *complex relative* permittivities; the corresponding ADE
//! update kernels live in [`crate::dispersive`].
//!
//! ## Walking-skeleton scope (Phase 2.fdtd.3)
//!
//! - Single-pole models only; multi-pole Lorentz / dual Debye are future work.
//! - The [`MaterialMap`] allocates one [`Material`] per E-field cell — the
//!   "full-grid" sparse-free layout. A truly sparse representation that only
//!   tracks tagged cells is a Phase 2.fdtd.3.1 optimisation.
//! - Material assignment is performed via axis-aligned [`MaterialMap::set_box`]
//!   regions; richer geometry primitives (sphere, polygon) come later.

use ndarray::Array3;
use num_complex::Complex64;

/// Single-pole dispersive material model.
///
/// The complex relative permittivity ε_r(ω) can be evaluated via
/// [`Material::permittivity`]; the parameter conventions match Taflove &
/// Hagness §9.4–9.6.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Material {
    /// Lossless vacuum: `ε_r = 1`. The default for un-tagged cells.
    Vacuum,

    /// Drude pole.
    ///
    /// `ε_r(ω) = eps_inf − ω_p² / (ω² − j γ ω)`
    ///
    /// (Taflove eq. 9.7). `omega_p` is the plasma frequency (rad/s) and
    /// `gamma` is the collision rate (rad/s).
    Drude {
        /// High-frequency permittivity (dimensionless; ≥ 1 typically).
        eps_inf: f64,
        /// Plasma frequency (rad/s).
        omega_p: f64,
        /// Collision damping rate (rad/s).
        gamma: f64,
    },

    /// Single-pole Lorentz oscillator.
    ///
    /// `ε_r(ω) = eps_inf + Δε · ω₀² / (ω₀² − ω² + 2 j δ ω)`
    ///
    /// (Taflove eq. 9.8). `omega_0` is the resonance frequency (rad/s) and
    /// `delta` is the damping coefficient (rad/s).
    Lorentz {
        /// High-frequency permittivity (dimensionless).
        eps_inf: f64,
        /// Oscillator strength (dimensionless).
        delta_eps: f64,
        /// Resonance frequency (rad/s).
        omega_0: f64,
        /// Damping coefficient (rad/s).
        delta: f64,
    },

    /// Single-pole Debye relaxation.
    ///
    /// `ε_r(ω) = eps_inf + Δε / (1 + j ω τ)`
    ///
    /// (Taflove eq. 9.6). `tau` is the relaxation time (s).
    Debye {
        /// High-frequency permittivity (dimensionless).
        eps_inf: f64,
        /// Static excess permittivity (dimensionless).
        delta_eps: f64,
        /// Relaxation time (s).
        tau: f64,
    },
}

impl Material {
    /// Complex relative permittivity at angular frequency `omega` (rad/s).
    ///
    /// Implements the model formulas above verbatim. Used by both the unit
    /// tests and by the analytical Fresnel reflection check in the
    /// `dispersive` integration test.
    pub fn permittivity(&self, omega: f64) -> Complex64 {
        match *self {
            Material::Vacuum => Complex64::new(1.0, 0.0),
            Material::Drude {
                eps_inf,
                omega_p,
                gamma,
            } => {
                // ε = ε∞ − ω_p² / (ω² − j γ ω)
                let denom = Complex64::new(omega * omega, -gamma * omega);
                Complex64::new(eps_inf, 0.0) - Complex64::new(omega_p * omega_p, 0.0) / denom
            }
            Material::Lorentz {
                eps_inf,
                delta_eps,
                omega_0,
                delta,
            } => {
                // ε = ε∞ + Δε · ω₀² / (ω₀² − ω² + 2 j δ ω)
                let denom = Complex64::new(omega_0 * omega_0 - omega * omega, 2.0 * delta * omega);
                Complex64::new(eps_inf, 0.0)
                    + Complex64::new(delta_eps * omega_0 * omega_0, 0.0) / denom
            }
            Material::Debye {
                eps_inf,
                delta_eps,
                tau,
            } => {
                // ε = ε∞ + Δε / (1 + j ω τ)
                let denom = Complex64::new(1.0, omega * tau);
                Complex64::new(eps_inf, 0.0) + Complex64::new(delta_eps, 0.0) / denom
            }
        }
    }

    /// `true` if this material is *not* [`Material::Vacuum`].
    ///
    /// Hot path: [`crate::dispersive::DispersiveState::update_e_with_dispersion`]
    /// uses this to branch between vacuum and ADE updates per cell.
    #[inline]
    pub fn is_dispersive(&self) -> bool {
        !matches!(self, Material::Vacuum)
    }
}

/// Per-E-field-cell dispersive-material tag.
///
/// Sized to match an "average" E-field cell extent (`[nx+1, ny+1, nz+1]`) so
/// a single index can address every E component. Cells not explicitly set are
/// [`Material::Vacuum`] and fall through to the standard vacuum update.
#[derive(Debug, Clone)]
pub struct MaterialMap {
    /// Number of primary cells along x.
    pub nx: usize,
    /// Number of primary cells along y.
    pub ny: usize,
    /// Number of primary cells along z.
    pub nz: usize,
    /// One [`Material`] per cell. Shape `[nx+1, ny+1, nz+1]`; the union of
    /// every staggered E component shape fits inside.
    pub cells: Array3<Material>,
}

impl MaterialMap {
    /// Build an `nx × ny × nz` map filled with [`Material::Vacuum`].
    pub fn vacuum(nx: usize, ny: usize, nz: usize) -> Self {
        let cells = Array3::from_elem((nx + 1, ny + 1, nz + 1), Material::Vacuum);
        Self { nx, ny, nz, cells }
    }

    /// Tag an inclusive-exclusive axis-aligned box `[i0, i1) × [j0, j1) × [k0, k1)`
    /// with `material`.
    ///
    /// Indices are clamped to the underlying [`Self::cells`] extent; out-of-range
    /// boxes are silently truncated rather than panicking — this matches how
    /// CAD-derived geometry often spills slightly past the simulation domain.
    #[allow(clippy::too_many_arguments)]
    pub fn set_box(
        &mut self,
        i0: usize,
        i1: usize,
        j0: usize,
        j1: usize,
        k0: usize,
        k1: usize,
        material: Material,
    ) {
        let (ni, nj, nk) = self.cells.dim();
        let i1 = i1.min(ni);
        let j1 = j1.min(nj);
        let k1 = k1.min(nk);
        if i0 >= i1 || j0 >= j1 || k0 >= k1 {
            return;
        }
        for i in i0..i1 {
            for j in j0..j1 {
                for k in k0..k1 {
                    self.cells[(i, j, k)] = material;
                }
            }
        }
    }

    /// Read the material tag at cell `(i, j, k)`.
    #[inline]
    pub fn material_at(&self, i: usize, j: usize, k: usize) -> Material {
        self.cells[(i, j, k)]
    }

    /// Number of cells whose material is *not* [`Material::Vacuum`].
    ///
    /// Diagnostic helper — used by tests and by progress logging.
    pub fn dispersive_cell_count(&self) -> usize {
        self.cells.iter().filter(|m| m.is_dispersive()).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn vacuum_permittivity_is_unity() {
        let m = Material::Vacuum;
        let e = m.permittivity(2.0 * PI * 1e9);
        assert!((e.re - 1.0).abs() < 1e-15);
        assert!(e.im.abs() < 1e-15);
    }

    #[test]
    fn drude_static_limit_diverges_lossless() {
        // Lossless Drude (γ = 0) at ω → 0 diverges to −∞ along real axis.
        let m = Material::Drude {
            eps_inf: 1.0,
            omega_p: 2.0 * PI * 1e10,
            gamma: 0.0,
        };
        let e_lo = m.permittivity(2.0 * PI * 1e8);
        let e_hi = m.permittivity(2.0 * PI * 1e12);
        // Below plasma: negative real part.
        assert!(
            e_lo.re < 0.0,
            "Drude eps below ω_p should be negative, got {}",
            e_lo.re
        );
        // Far above plasma: approaches eps_inf from below.
        assert!(e_hi.re > 0.0 && e_hi.re < 1.0);
    }

    #[test]
    fn drude_lossy_has_negative_imag() {
        // Convention ε = ε∞ − ω_p² / (ω² − jγω). Expanding,
        //   1 / (ω² − jγω) = (ω² + jγω) / ((ω²)² + (γω)²)
        // → Im(1/(ω²−jγω)) > 0. Subtracting ω_p²·(...) then gives
        // Im(ε) < 0, the standard "lossy passive medium" sign under the
        // physics-convention exp(−jωt) Fourier transform used in Taflove.
        let m = Material::Drude {
            eps_inf: 1.0,
            omega_p: 2.0 * PI * 1e10,
            gamma: 2.0 * PI * 1e8,
        };
        let e = m.permittivity(2.0 * PI * 1e10);
        assert!(
            e.im < 0.0,
            "lossy Drude should have negative Im(ε), got {}",
            e.im
        );
    }

    #[test]
    fn lorentz_static_value_is_eps_inf_plus_delta_eps() {
        // At ω = 0: ε = ε∞ + Δε · ω₀² / ω₀² = ε∞ + Δε.
        let m = Material::Lorentz {
            eps_inf: 2.25,
            delta_eps: 1.0,
            omega_0: 2.0 * PI * 1e15,
            delta: 1.0e13,
        };
        let e = m.permittivity(0.0);
        assert!((e.re - 3.25).abs() < 1e-10);
        assert!(e.im.abs() < 1e-10);
    }

    #[test]
    fn lorentz_high_freq_limit_is_eps_inf() {
        let m = Material::Lorentz {
            eps_inf: 2.25,
            delta_eps: 1.0,
            omega_0: 2.0 * PI * 1e10,
            delta: 1.0e8,
        };
        // ω ≫ ω₀ ⇒ ε → ε∞.
        let e = m.permittivity(2.0 * PI * 1e14);
        assert!((e.re - 2.25).abs() < 1e-3, "got {}", e.re);
    }

    #[test]
    fn debye_static_and_infinite_limits() {
        let m = Material::Debye {
            eps_inf: 3.0,
            delta_eps: 75.0,
            tau: 8.0e-12,
        };
        // Static: ε∞ + Δε.
        let e0 = m.permittivity(0.0);
        assert!((e0.re - 78.0).abs() < 1e-10);
        assert!(e0.im.abs() < 1e-10);
        // Very high: ε → ε∞, imag → 0.
        let e_hi = m.permittivity(2.0 * PI * 1e16);
        assert!((e_hi.re - 3.0).abs() < 1e-3, "got {}", e_hi.re);
    }

    #[test]
    fn material_map_default_is_vacuum() {
        let map = MaterialMap::vacuum(8, 8, 8);
        assert_eq!(map.material_at(3, 4, 5), Material::Vacuum);
        assert_eq!(map.dispersive_cell_count(), 0);
    }

    #[test]
    fn material_map_set_box_tags_cells() {
        let mut map = MaterialMap::vacuum(8, 8, 8);
        let drude = Material::Drude {
            eps_inf: 1.0,
            omega_p: 1.0,
            gamma: 0.0,
        };
        map.set_box(2, 5, 0, 9, 0, 9, drude);
        assert_eq!(map.material_at(3, 4, 5), drude);
        assert_eq!(map.material_at(1, 4, 5), Material::Vacuum);
        // (2..5) × (0..9) × (0..9) = 3 · 9 · 9 = 243 cells.
        assert_eq!(map.dispersive_cell_count(), 243);
    }

    #[test]
    fn material_map_set_box_clamps_oversized() {
        let mut map = MaterialMap::vacuum(4, 4, 4);
        let drude = Material::Drude {
            eps_inf: 1.0,
            omega_p: 1.0,
            gamma: 0.0,
        };
        // (0..100) clamped to (0..5).
        map.set_box(0, 100, 0, 100, 0, 100, drude);
        assert_eq!(map.dispersive_cell_count(), 5 * 5 * 5);
    }

    #[test]
    fn material_map_set_box_empty_is_noop() {
        let mut map = MaterialMap::vacuum(4, 4, 4);
        let drude = Material::Drude {
            eps_inf: 1.0,
            omega_p: 1.0,
            gamma: 0.0,
        };
        map.set_box(5, 5, 0, 4, 0, 4, drude);
        assert_eq!(map.dispersive_cell_count(), 0);
    }
}
