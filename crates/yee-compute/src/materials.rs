//! Per-cell material maps, interior PEC masks, and boundary selection (E.1).

use crate::spec::{FdtdSpec, len3};

/// Optional per-cell material maps and interior PEC masks.
///
/// Conventions match `yee_fdtd::grid::YeeGrid` exactly: the ε_r / μ_r / σ
/// maps are `[nx+1, ny+1, nz+1]` (oversized so any staggered component can
/// be addressed by its own primary `(i, j, k)`); the PEC masks share their
/// E component's staggered shape and are clamped to zero after each E
/// half-step. All-`None` is the uniform E.0 behaviour.
#[derive(Debug, Clone, Default)]
pub struct Materials {
    /// Per-cell relative permittivity, `[nx+1, ny+1, nz+1]` row-major.
    pub eps_r_cells: Option<Vec<f64>>,
    /// Per-cell relative permeability, `[nx+1, ny+1, nz+1]` row-major.
    pub mu_r_cells: Option<Vec<f64>>,
    /// Per-cell electric conductivity (S/m), `[nx+1, ny+1, nz+1]` row-major.
    /// When `Some`, the E update uses the lossy CA/CB form (Taflove §3.7).
    pub sigma_cells: Option<Vec<f64>>,
    /// Interior PEC mask for `E_x` (shape of `E_x`).
    pub pec_mask_ex: Option<Vec<bool>>,
    /// Interior PEC mask for `E_y` (shape of `E_y`).
    pub pec_mask_ey: Option<Vec<bool>>,
    /// Interior PEC mask for `E_z` (shape of `E_z`).
    pub pec_mask_ez: Option<Vec<bool>>,
}

impl Materials {
    /// Panic unless every present map/mask has the length its shape demands.
    pub(crate) fn validate(&self, spec: &FdtdSpec) {
        let cells = (spec.nx + 1) * (spec.ny + 1) * (spec.nz + 1);
        for (map, name) in [
            (&self.eps_r_cells, "eps_r_cells"),
            (&self.mu_r_cells, "mu_r_cells"),
            (&self.sigma_cells, "sigma_cells"),
        ] {
            if let Some(m) = map {
                assert_eq!(m.len(), cells, "{name} length mismatch");
            }
        }
        for (mask, len, name) in [
            (&self.pec_mask_ex, len3(spec.ex_dims()), "pec_mask_ex"),
            (&self.pec_mask_ey, len3(spec.ey_dims()), "pec_mask_ey"),
            (&self.pec_mask_ez, len3(spec.ez_dims()), "pec_mask_ez"),
        ] {
            if let Some(m) = mask {
                assert_eq!(m.len(), len, "{name} length mismatch");
            }
        }
    }

    /// True when any interior PEC mask is attached. (Consumed by the GPU
    /// backend's arena setup; unused in a CPU-only build.)
    #[cfg_attr(not(feature = "gpu"), allow(dead_code))]
    pub(crate) fn has_mask(&self) -> bool {
        self.pec_mask_ex.is_some() || self.pec_mask_ey.is_some() || self.pec_mask_ez.is_some()
    }
}

/// Outer-boundary treatment for a stepper.
#[derive(Debug, Clone)]
pub enum Boundary {
    /// No boundary phase at all — the raw E.0 kernel semantics (outer
    /// tangential E faces are simply never written). This is what the
    /// bit-exact kernel gate `compute-001` exercises.
    None,
    /// Legacy reflecting PEC box: outer tangential E clamped to zero after
    /// each half-step, mirroring `yee_fdtd::boundary::apply_pec`. Kept as
    /// the reflecting reference for the CPML gate.
    PecBox,
    /// Roden–Gedney 2000 CPML on the enabled axes.
    Cpml(CpmlConfig),
}

/// CPML configuration, mirroring `yee_fdtd::cpml::CpmlParams`.
#[derive(Debug, Clone, Copy)]
pub struct CpmlConfig {
    /// PML thickness in cells on each enabled face. Standard: 10.
    pub npml: usize,
    /// Polynomial grading order. Standard: 3.
    pub m: i32,
    /// Peak conductivity; see [`CpmlConfig::for_spec`] for the standard
    /// `R_0 = 1e-6` recipe.
    pub sigma_max: f64,
    /// Peak coordinate-stretching factor. Standard: 1.0.
    pub kappa_max: f64,
    /// Peak CFS shift parameter. Standard: 0.05.
    pub alpha_max: f64,
    /// Per-axis enable mask `[x, y, z]`; a disabled axis carries no CPML.
    /// Kept in sync with [`CpmlConfig::faces`] by the builders; `faces` is
    /// the source of truth for the CPU stepper.
    pub axes: [bool; 3],
    /// Per-face enable `[[x−, x+], [y−, y+], [z−, z+]]` (A.2, ADR-0192):
    /// a disabled face stays PEC (the outer tangential E is simply never
    /// written) while the opposite face can absorb — e.g. an antenna's
    /// open top over a PEC ground is `[[t, t], [t, t], [f, t]]`. The GPU
    /// backend supports only face-symmetric configs (equal per axis) and
    /// rejects others with `ComputeError::Unsupported`.
    pub faces: [[bool; 2]; 3],
}

impl CpmlConfig {
    /// Standard Roden–Gedney/Taflove parameter set sized to `spec`
    /// (`σ_max = −(m+1)·ln(1e-6) / (2·η₀·npml·dx)`), all axes enabled —
    /// the exact `CpmlParams::for_grid` recipe.
    pub fn for_spec(spec: &FdtdSpec, npml: usize) -> Self {
        use yee_core::units::ETA0;
        let m = 3i32;
        let r0: f64 = 1.0e-6;
        let sigma_max = -(f64::from(m) + 1.0) * r0.ln() / (2.0 * ETA0 * (npml as f64) * spec.dx);
        Self {
            npml,
            m,
            sigma_max,
            kappa_max: 1.0,
            alpha_max: 0.05,
            axes: [true; 3],
            faces: [[true; 2]; 3],
        }
    }

    /// Set the per-axis enable mask (both faces of each axis) and return
    /// `self`.
    #[must_use]
    pub fn with_axes(mut self, axes: [bool; 3]) -> Self {
        self.axes = axes;
        self.faces = [[axes[0]; 2], [axes[1]; 2], [axes[2]; 2]];
        self
    }

    /// Set the per-face enable mask `[[x−, x+], [y−, y+], [z−, z+]]` and
    /// return `self` (A.2). `axes` is kept consistent (an axis counts as
    /// enabled when either of its faces is).
    #[must_use]
    pub fn with_faces(mut self, faces: [[bool; 2]; 3]) -> Self {
        self.faces = faces;
        self.axes = [
            faces[0][0] || faces[0][1],
            faces[1][0] || faces[1][1],
            faces[2][0] || faces[2][1],
        ];
        self
    }

    /// True when every axis has both faces equal. Historical note: until
    /// R.3 this was the only shape the GPU backend could express; both
    /// backends now honor arbitrary per-face masks.
    pub fn faces_are_axis_symmetric(&self) -> bool {
        self.faces.iter().all(|f| f[0] == f[1])
    }
}
