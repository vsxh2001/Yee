//! Surface-roughness loss models for conductor surfaces.
//!
//! Phase 1.4.0 walking skeleton: three frequency-dependent loss-multiplier
//! models are implemented, each returning a scalar `K(f)` that scales the
//! per-edge conductor loss relative to a perfectly smooth surface.
//!
//! # Limitation (Phase 1.4.0)
//!
//! The walking skeleton applies `K(f)` *uniformly* — either as a uniform
//! per-pair multiplier in [`crate::fill`] or as a post-fill scalar applied
//! to the entire impedance matrix in [`crate::solve`] (this crate uses the
//! latter, simpler approach: see [`crate::solve::apply_roughness`]). The
//! resulting effect on `Z_in` is approximately correct because the matrix is
//! multiplied uniformly, but it does *not* discriminate between edges whose
//! tangential current aligns with the rough surface and edges whose current
//! lies along a smooth co-planar boundary. Phase 1.4.1 will apply roughness
//! per-edge based on edge tangential-current alignment with the conductor
//! surface normal.
//!
//! # References
//!
//! - Hammerstad & Jensen, "Accurate models for microstrip computer-aided
//!   design", *IEEE MTT-S Int. Microwave Symp. Digest* (1980) 407–409.
//! - Groiss et al., "Modelling of conductor surface roughness on multilayer
//!   PCB structures for high frequency applications", *Proc. EMC Europe*
//!   (1996).
//! - Huray, *The Foundations of Signal Integrity*, Wiley (2010) — small-sphere
//!   ("snowball") model.

/// Frequency-dependent conductor surface-roughness loss multiplier.
///
/// Each variant returns a dimensionless `K(f) ≥ 1` from
/// [`RoughnessModel::loss_multiplier`]; `K = 1` means a perfectly smooth
/// surface (no excess loss).
#[derive(Debug, Clone, Copy)]
pub enum RoughnessModel {
    /// Smooth conductor — loss multiplier is identically 1 at every frequency.
    Smooth,
    /// Hammerstad–Jensen: `K_HJ = 1 + (2/π) · arctan(1.4 · (Δ/δ)²)` where
    /// `Δ` is the RMS surface roughness (m) and `δ` is the conductor skin
    /// depth at the operating frequency. Saturates at `K → 2` for very rough
    /// surfaces (Δ ≫ δ).
    HammerstadJensen {
        /// RMS surface roughness Δ, metres.
        delta_rms_m: f64,
    },
    /// Groiss: `K_G = K_HJ · (1 + 2 · (Δ/δ)²)` — same `Δ`/`δ` ratio as
    /// Hammerstad–Jensen but with a multiplicative `(1 + 2·(Δ/δ)²)` factor
    /// that drives a sharper rise than the arctangent saturation alone. Used
    /// where measurements show HJ under-predicts loss at high `Δ/δ`.
    Groiss {
        /// RMS surface roughness Δ, metres.
        delta_rms_m: f64,
    },
    /// Huray small-sphere ("snowball") model, simplified one-tier form:
    ///
    /// ```text
    /// K_Hu = 1 + surface_ratio · (1 + δ/a + (δ/a)²/2) / (1 + δ/a)
    /// ```
    ///
    /// where `a` is the sphere radius (m) and `surface_ratio` lumps the
    /// per-unit-flat-surface sphere coverage `N · 4πa² / A_flat` into a single
    /// dimensionless parameter. Phase 1.4.0 deliberately collapses the full
    /// `(a, N, A_flat)` triple into `(sphere_radius_m, surface_ratio)` so
    /// callers tune two physical knobs; a fully resolved multi-tier sphere
    /// distribution is deferred to Phase 1.4.1+.
    Huray {
        /// Sphere ("snowball") radius `a`, metres.
        sphere_radius_m: f64,
        /// Dimensionless surface ratio `N · 4πa² / A_flat` lumping coverage
        /// density and sphere geometry into a single tunable.
        surface_ratio: f64,
    },
}

impl RoughnessModel {
    /// Compute the loss multiplier `K(f)` at `freq_hz` on a conductor of
    /// conductivity `sigma_s_per_m` (S/m).
    ///
    /// The skin depth used internally is
    /// `δ = sqrt(2 / (ω · μ₀ · σ))` with `ω = 2π · f`. For copper at room
    /// temperature, pass [`SIGMA_COPPER`].
    pub fn loss_multiplier(&self, freq_hz: f64, sigma_s_per_m: f64) -> f64 {
        let omega = std::f64::consts::TAU * freq_hz;
        let skin_depth = (2.0 / (omega * yee_core::units::MU0 * sigma_s_per_m)).sqrt();
        match self {
            Self::Smooth => 1.0,
            Self::HammerstadJensen { delta_rms_m } => {
                let r = delta_rms_m / skin_depth;
                1.0 + (2.0 / std::f64::consts::PI) * (1.4 * r * r).atan()
            }
            Self::Groiss { delta_rms_m } => {
                let r = delta_rms_m / skin_depth;
                let k_hj = 1.0 + (2.0 / std::f64::consts::PI) * (1.4 * r * r).atan();
                k_hj * (1.0 + 2.0 * r * r)
            }
            Self::Huray {
                sphere_radius_m,
                surface_ratio,
            } => {
                let a = *sphere_radius_m;
                let delta = skin_depth;
                let ratio = delta / a;
                let factor = (1.0 + ratio + 0.5 * ratio * ratio) / (1.0 + ratio);
                1.0 + surface_ratio * factor
            }
        }
    }
}

/// Copper conductivity at room temperature, S/m. Default for all Phase 1.4
/// callers that do not pass an explicit `sigma_s_per_m`.
pub const SIGMA_COPPER: f64 = 5.8e7;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smooth_returns_unity() {
        assert_eq!(
            RoughnessModel::Smooth.loss_multiplier(1.0e9, SIGMA_COPPER),
            1.0
        );
    }

    #[test]
    fn hammerstad_jensen_saturates_at_two() {
        // Limit: very rough surface (Δ >> δ) → K → 1 + 2/π · π/2 = 2.
        let very_rough = RoughnessModel::HammerstadJensen {
            delta_rms_m: 1.0e-3,
        };
        let k = very_rough.loss_multiplier(1.0e9, SIGMA_COPPER);
        assert!((k - 2.0).abs() < 0.01, "HJ saturation: K={k}");
    }

    #[test]
    fn skin_depth_monotonic_in_frequency() {
        let r = RoughnessModel::HammerstadJensen {
            delta_rms_m: 1.0e-6,
        };
        let k_low = r.loss_multiplier(1.0e8, SIGMA_COPPER);
        let k_high = r.loss_multiplier(1.0e10, SIGMA_COPPER);
        assert!(k_high > k_low, "K must increase with frequency");
    }

    #[test]
    fn groiss_exceeds_hj_at_high_roughness() {
        let g = RoughnessModel::Groiss {
            delta_rms_m: 5.0e-6,
        };
        let hj = RoughnessModel::HammerstadJensen {
            delta_rms_m: 5.0e-6,
        };
        let k_g = g.loss_multiplier(10.0e9, SIGMA_COPPER);
        let k_hj = hj.loss_multiplier(10.0e9, SIGMA_COPPER);
        assert!(k_g > k_hj);
    }

    #[test]
    fn huray_smooth_limit_when_surface_ratio_zero() {
        let h = RoughnessModel::Huray {
            sphere_radius_m: 1.0e-6,
            surface_ratio: 0.0,
        };
        assert!((h.loss_multiplier(1.0e9, SIGMA_COPPER) - 1.0).abs() < 1e-12);
    }
}
