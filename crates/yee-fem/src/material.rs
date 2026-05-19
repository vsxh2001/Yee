//! Per-tet dispersive material model for the Phase 4.fem.eig.1 FEM
//! eigensolver.
//!
//! Defines a multi-pole [`Material`] (`ε_∞`, real `μ_r`, and a `Vec<MaterialPole>`
//! of Drude / Lorentz / Debye poles) and a [`MaterialDatabase`] that maps
//! [`MaterialTag`]s to per-region materials. The FEM dispersive solver
//! (Phase 4.fem.eig.1 plan steps D4/D5) evaluates `ε(ω)`, `μ(ω)` at each
//! Newton-tracker trial frequency via the database accessors below.
//!
//! ## Relationship to Phase 2.fdtd.3
//!
//! The ADE update kernels in `yee-fdtd::material::Material` model a *single*
//! single-pole Drude / Lorentz / Debye material per cell. The FEM dispersive
//! solver needs a slightly richer surface — multi-pole sums per material
//! (anticipating Phase 4.fem.eig.1.1) and a tag → material lookup keyed by
//! [`MaterialTag`] (mirroring the `yee-mesh` tag plumbing) — so this module
//! defines a peer multi-pole `Material` rather than reusing the FDTD enum
//! directly. ADR-0039 §"Material relocation" anticipates promoting the
//! shared core to `yee-core` once both crates agree on a single surface; v0
//! ships the FEM-side type here in `yee-fem` to avoid widening the lane.
//!
//! ## Sign convention
//!
//! `ε(ω) = ε_∞ + Σ_p [pole contribution]` with the per-pole contributions
//! listed under [`MaterialPole`]. Each contribution carries its own sign:
//! the Drude contribution is `−ω_p² / (ω² + jγω)` (so Re(Drude) < 0 below
//! the plasma frequency), Lorentz is `+ω_p² / (ω_0² − ω² − jγω)`, and
//! Debye is `+(ε_s − ε_∞) / (1 + jωτ)`. The static-limit identity
//! `ε(0) = ε_∞ + (ε_s − ε_∞) = ε_s` for Debye is exact under this
//! convention. The `+jγω` (rather than `−jγω`) damping sign mirrors the
//! engineering `exp(+jωt)` Fourier convention — note this differs from the
//! `yee-fdtd::material::Material` ADE enum's Taflove-style `−jγω` and is
//! therefore *not* a re-export of that type.
//!
//! ## Walking-skeleton scope
//!
//! - `μ_r` is real and frequency-independent in v1 per spec §2 ("magnetic
//!   dispersion μ(ω) is Phase 4.fem.eig.1.2"). [`Material::mu_at`] returns
//!   `Complex64::new(mu_r, 0.0)` regardless of `ω`.
//! - Multi-pole sums are permitted by construction (`poles: Vec<_>`); the
//!   v1 validation gate `fem-eig-002` exercises a single-pole Drude only.
//! - Anisotropy is out of scope; ε is a scalar `Complex64` per material.

use num_complex::Complex64;
use yee_mesh::MaterialTag;

/// A single dispersive pole.
///
/// All fields are real-valued physical parameters. The complex contribution
/// to `ε(ω)` is computed by [`Material::eps_at`] under the convention
/// `ε(ω) = ε_∞ + Σ_p [pole(ω)]` (see module-level docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialPole {
    /// Drude pole.
    ///
    /// Contribution to `ε(ω)`:
    ///
    /// ```text
    ///     − ω_p² / (ω² + j γ ω)
    /// ```
    ///
    /// Below the plasma frequency `ω_p` the real part of this contribution
    /// is negative and dominates `ε_∞`, producing the classic
    /// `Re(ε) < 0` metallic behaviour.
    Drude {
        /// Plasma frequency (rad/s).
        omega_p: f64,
        /// Collision damping rate (rad/s).
        gamma: f64,
    },

    /// Single-pole Lorentz oscillator.
    ///
    /// Contribution to `ε(ω)`:
    ///
    /// ```text
    ///     + ω_p² / (ω_0² − ω² − j γ ω)
    /// ```
    Lorentz {
        /// Resonance frequency (rad/s).
        omega_0: f64,
        /// Oscillator strength (rad/s).
        omega_p: f64,
        /// Damping coefficient (rad/s).
        gamma: f64,
    },

    /// Single-pole Debye relaxation.
    ///
    /// Contribution to `ε(ω)`:
    ///
    /// ```text
    ///     + (ε_s − ε_∞) / (1 + j ω τ)
    /// ```
    ///
    /// Static limit `ω → 0`: contribution → `ε_s − ε_∞`, so
    /// `ε(0) = ε_∞ + (ε_s − ε_∞) = ε_s`.
    Debye {
        /// Static excess permittivity `ε_s − ε_∞` (dimensionless).
        eps_s_minus_eps_inf: f64,
        /// Relaxation time (s).
        tau: f64,
    },
}

impl MaterialPole {
    /// The contribution this pole adds to `ε_∞` at angular frequency
    /// `omega` (rad/s).
    ///
    /// [`Material::eps_at`] accumulates these by computing
    /// `ε(ω) = ε_∞ + Σ_p p.contribution(ω)`. Per-variant sign is encoded
    /// here, not in the outer sum.
    fn contribution(&self, omega: f64) -> Complex64 {
        match *self {
            MaterialPole::Drude { omega_p, gamma } => {
                // − ω_p² / (ω² + j γ ω)
                let denom = Complex64::new(omega * omega, gamma * omega);
                -Complex64::new(omega_p * omega_p, 0.0) / denom
            }
            MaterialPole::Lorentz {
                omega_0,
                omega_p,
                gamma,
            } => {
                // + ω_p² / (ω_0² − ω² − j γ ω)
                let denom = Complex64::new(omega_0 * omega_0 - omega * omega, -gamma * omega);
                Complex64::new(omega_p * omega_p, 0.0) / denom
            }
            MaterialPole::Debye {
                eps_s_minus_eps_inf,
                tau,
            } => {
                // + (ε_s − ε_∞) / (1 + j ω τ)
                let denom = Complex64::new(1.0, omega * tau);
                Complex64::new(eps_s_minus_eps_inf, 0.0) / denom
            }
        }
    }
}

/// A dispersive material: high-frequency permittivity `ε_∞`, real
/// permeability `μ_r`, and a list of dispersive poles.
///
/// The complex permittivity at angular frequency `ω` is
///
/// ```text
///     ε(ω) = ε_∞ − Σ_p [pole_p contribution]
/// ```
///
/// per [`Material::eps_at`]; the per-pole contributions follow the
/// conventions in [`MaterialPole`]. v1 ships real, frequency-independent
/// `μ_r` only — [`Material::mu_at`] is a placeholder returning
/// `Complex64::new(mu_r, 0.0)` (see module-level docs).
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    /// High-frequency relative permittivity (dimensionless, ≥ 1 typically).
    pub eps_inf: f64,
    /// Relative permeability (dimensionless). Real and frequency-independent
    /// in v1 — magnetic dispersion is Phase 4.fem.eig.1.2.
    pub mu_r: f64,
    /// Dispersive poles. Empty for a non-dispersive constant-`ε_∞`
    /// material; multiple entries are summed in [`Material::eps_at`].
    pub poles: Vec<MaterialPole>,
}

impl Default for Material {
    /// Free-space default: `ε_∞ = 1`, `μ_r = 1`, no dispersive poles.
    fn default() -> Self {
        Self {
            eps_inf: 1.0,
            mu_r: 1.0,
            poles: Vec::new(),
        }
    }
}

impl Material {
    /// Complex relative permittivity `ε_r(ω)` at angular frequency `omega`
    /// (rad/s).
    ///
    /// Implements `ε(ω) = ε_∞ + Σ_p [pole_p contribution(ω)]`. Returns
    /// `Complex64::new(eps_inf, 0.0)` when [`Material::poles`] is empty.
    /// Each pole's contribution is signed (see [`MaterialPole`]); the outer
    /// sum is unsigned.
    pub fn eps_at(&self, omega: f64) -> Complex64 {
        let mut eps = Complex64::new(self.eps_inf, 0.0);
        for pole in &self.poles {
            eps += pole.contribution(omega);
        }
        eps
    }

    /// Complex relative permeability `μ_r(ω)` at angular frequency `omega`
    /// (rad/s).
    ///
    /// v1 placeholder: magnetic dispersion is out of scope per spec §2, so
    /// this returns `Complex64::new(self.mu_r, 0.0)` unconditionally. The
    /// `_omega` parameter is reserved for the Phase 4.fem.eig.1.2 magnetic
    /// dispersion lift; the FEM Newton tracker already takes the
    /// `(ε, μ)` pair as `Complex64` so the lift is a pure-arithmetic change
    /// inside this function.
    pub fn mu_at(&self, _omega: f64) -> Complex64 {
        Complex64::new(self.mu_r, 0.0)
    }
}

/// A tag → [`Material`] lookup for the FEM dispersive eigensolver.
///
/// Built up via the [`MaterialDatabase::with_material`] builder; queries
/// resolve a [`MaterialTag`] (the `yee-mesh` tag type) to a `Material` and
/// evaluate its complex `ε(ω)` or `μ(ω)`. Unregistered tags fall back to
/// the [`Material::default`] free-space response (`ε = 1`, `μ = 1`), which
/// matches the FDTD-side convention from `yee-fdtd::material::MaterialMap`.
///
/// ## Storage
///
/// Materials are stored as a `Vec<(MaterialTag, Material)>` rather than a
/// `HashMap` because (a) FEM material counts are O(10) per simulation —
/// linear scan is unconditionally faster than hashing at that scale — and
/// (b) the linear-scan order is reproducible across runs, which simplifies
/// regression tests. The Phase 4.fem.eig.2 production-scale lift may
/// re-evaluate this choice if material counts grow.
#[derive(Debug, Default, Clone)]
pub struct MaterialDatabase {
    materials: Vec<(MaterialTag, Material)>,
}

impl MaterialDatabase {
    /// Construct an empty database.
    ///
    /// Every [`MaterialTag`] lookup against an empty database falls back to
    /// the [`Material::default`] free-space response.
    pub fn new() -> Self {
        Self {
            materials: Vec::new(),
        }
    }

    /// Builder: register `mat` under `tag`, returning `self` so calls can
    /// chain.
    ///
    /// If `tag` is already registered, the new entry is appended rather than
    /// replacing the existing one. The internal lookup returns the *first*
    /// match, so the original entry remains authoritative — keep builder
    /// chains free of duplicates.
    pub fn with_material(mut self, tag: MaterialTag, mat: Material) -> Self {
        self.materials.push((tag, mat));
        self
    }

    /// Lookup helper: returns the first registered [`Material`] whose tag
    /// matches `tag`, or `None` if `tag` is unregistered. Caller-facing
    /// accessors below short-circuit `None` to the free-space response.
    fn lookup(&self, tag: MaterialTag) -> Option<&Material> {
        self.materials
            .iter()
            .find(|(t, _)| *t == tag)
            .map(|(_, m)| m)
    }

    /// Complex relative permittivity for the material tagged `tag` at
    /// angular frequency `omega` (rad/s).
    ///
    /// An unregistered tag returns `Complex64::new(1.0, 0.0)` (free space).
    pub fn eps_at(&self, tag: MaterialTag, omega: f64) -> Complex64 {
        match self.lookup(tag) {
            Some(m) => m.eps_at(omega),
            None => Complex64::new(1.0, 0.0),
        }
    }

    /// Complex relative permeability for the material tagged `tag` at
    /// angular frequency `omega` (rad/s).
    ///
    /// An unregistered tag returns `Complex64::new(1.0, 0.0)` (free space).
    /// v1 placeholder: see [`Material::mu_at`].
    pub fn mu_at(&self, tag: MaterialTag, omega: f64) -> Complex64 {
        match self.lookup(tag) {
            Some(m) => m.mu_at(omega),
            None => Complex64::new(1.0, 0.0),
        }
    }
}
