//! # yee-fdtd
//!
//! 3D FDTD on the Yee staggered grid. This crate currently ships the **Phase 2
//! walking skeleton**: a CPU-only, single-threaded, scalar (FP64) Yee solver
//! that demonstrates leapfrog propagation in vacuum on a uniform grid.
//!
//! ## What is included (Phase 2.0 + 2.1)
//!
//! - `YeeGrid` with vacuum constructor, Courant stability limit
//! - Scalar `update_e` / `update_h` kernels (Taflove & Hagness §3)
//! - Gaussian-in-time point source on `E_z`
//! - **CPML absorbing boundary on all six outer faces (Roden & Gedney 2000)**
//!   via [`CpmlState`] / [`CpmlParams`]
//! - Hard PEC fallback in [`boundary::apply_pec`] for cavity-style problems
//! - [`WalkingSkeletonSolver`]: a tiny [`FdtdSolver`] impl that wires it all
//!   together; choose absorbing vs reflecting boundaries via
//!   [`WalkingSkeletonSolver::with_cpml`] / [`WalkingSkeletonSolver::new`]
//!
//! ## What is NOT included
//!
//! - No GPU kernels, no multi-GPU domain decomposition.
//! - No subgridding, no dispersive materials (Drude / Lorentz / Debye).
//! - No conformal (Dey-Mittra) treatment of curved geometry.
//! - No NTFF, no lumped ports, no waveguide ports.
//!
//! These omissions are intentional. The walking skeleton exists so the rest of
//! the workspace (mesh, I/O, CLI, Python bindings) can integrate against a
//! real solver surface while the high-performance kernels are still in
//! development.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod boundary;
pub mod cpml;
pub mod grid;
pub mod ntff;
pub mod sources;
pub mod update;

pub use cpml::{CpmlParams, CpmlState};
pub use grid::YeeGrid;
pub use ntff::{NtffParams, NtffState};

/// FDTD-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Invalid input from the caller (out-of-range size, bad time step, etc.).
    #[error("invalid input: {0}")]
    Invalid(String),

    /// Numerical failure (NaN, divergence, instability).
    #[error("numerical failure: {0}")]
    Numerical(String),
}

/// FDTD-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;

/// Minimal solver-stepping interface.
///
/// Concrete solvers (walking skeleton, CUDA, multi-GPU) implement this. The
/// trait is intentionally tiny — richer features (sources, probes, NTFF) hang
/// off the concrete type.
pub trait FdtdSolver {
    /// Advance the simulation by exactly one Yee time step.
    fn step(&mut self);

    /// Return the simulation time at the start of the next step (seconds).
    fn current_time(&self) -> f64;
}

/// Reference implementation of [`FdtdSolver`] used by the walking skeleton.
///
/// The update order each step is:
///
/// 1. `update_h` — magnetic field leapfrogs forward by `dt`
/// 2. CPML auxiliary update for H (if a [`CpmlState`] is configured), or PEC
///    tangential-E clamp on the outer faces (legacy mode)
/// 3. *(optional)* inject a source via [`Self::step_with_source`]
/// 4. `update_e` — electric field leapfrogs forward by `dt`
/// 5. CPML auxiliary update for E, or PEC clamp again
/// 6. increment step counter / simulation clock
///
/// Either the CPML state or the PEC clamp manages the outer faces — never
/// both. CPML is preferred for open-domain problems; PEC is kept available
/// for cavity-style runs and for the regression test in
/// `tests/fdtd_propagation.rs`.
pub struct WalkingSkeletonSolver {
    grid: YeeGrid,
    step: u64,
    cpml: Option<CpmlState>,
}

impl WalkingSkeletonSolver {
    /// Wrap a [`YeeGrid`] in a fresh solver at `t = 0` with **hard PEC**
    /// outer boundaries (reflecting). For absorbing CPML boundaries see
    /// [`Self::with_cpml`].
    pub fn new(grid: YeeGrid) -> Self {
        Self {
            grid,
            step: 0,
            cpml: None,
        }
    }

    /// Wrap a [`YeeGrid`] in a fresh solver with a CPML absorbing boundary
    /// on all six outer faces, configured via `params`.
    ///
    /// The CPML state is built from `grid.dt`, so callers must not change
    /// `dt` after construction without rebuilding the solver.
    pub fn with_cpml(grid: YeeGrid, params: CpmlParams) -> Self {
        let cpml = CpmlState::new(&grid, params);
        Self {
            grid,
            step: 0,
            cpml: Some(cpml),
        }
    }

    /// Borrow the CPML state, if one was configured via [`Self::with_cpml`].
    pub fn cpml(&self) -> Option<&CpmlState> {
        self.cpml.as_ref()
    }

    /// Immutable view of the underlying grid (e.g. for probing field values).
    pub fn grid(&self) -> &YeeGrid {
        &self.grid
    }

    /// Mutable view of the underlying grid (escape hatch for callers that need
    /// to write into material/source state directly).
    pub fn grid_mut(&mut self) -> &mut YeeGrid {
        &mut self.grid
    }

    /// Current step index (0-based; equals the number of completed steps).
    pub fn step_index(&self) -> u64 {
        self.step
    }

    /// Time step in seconds.
    pub fn dt(&self) -> f64 {
        self.grid.dt
    }

    /// Step the solver while injecting a Gaussian-in-time pulse on `E_z` at
    /// cell `(i, j, k)`.
    ///
    /// The source contribution is added *between* the H and E updates so it
    /// is visible to the next E update through the standard leapfrog timing.
    /// The Gaussian is sampled at the current simulation time (before this
    /// step advances the clock).
    pub fn step_with_source(&mut self, i: usize, j: usize, k: usize, t0: f64, sigma: f64) {
        let t = self.current_time();
        update::update_h(&mut self.grid);
        if let Some(cpml) = self.cpml.as_mut() {
            cpml.update_h(&mut self.grid);
        } else {
            #[allow(deprecated)]
            boundary::apply_pec(&mut self.grid);
        }
        sources::gaussian_pulse_ez(&mut self.grid, i, j, k, t, t0, sigma);
        update::update_e(&mut self.grid);
        if let Some(cpml) = self.cpml.as_mut() {
            cpml.update_e(&mut self.grid);
        } else {
            #[allow(deprecated)]
            boundary::apply_pec(&mut self.grid);
        }
        self.step += 1;
    }

}

impl FdtdSolver for WalkingSkeletonSolver {
    fn step(&mut self) {
        update::update_h(&mut self.grid);
        if let Some(cpml) = self.cpml.as_mut() {
            cpml.update_h(&mut self.grid);
        } else {
            #[allow(deprecated)]
            boundary::apply_pec(&mut self.grid);
        }
        update::update_e(&mut self.grid);
        if let Some(cpml) = self.cpml.as_mut() {
            cpml.update_e(&mut self.grid);
        } else {
            #[allow(deprecated)]
            boundary::apply_pec(&mut self.grid);
        }
        self.step += 1;
    }

    fn current_time(&self) -> f64 {
        self.step as f64 * self.grid.dt
    }
}
