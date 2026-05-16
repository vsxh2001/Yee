//! End-to-end FDTD driver that wires `YeeGrid` + CPML + a `J_z` dipole
//! current source + NTFF into a single public entry point.
//!
//! This is the **Phase 2.fdtd.4** walking-skeleton driver: it composes the
//! pieces already shipped in [`crate::grid`], [`crate::cpml`],
//! [`crate::sources`], and [`crate::ntff`] into a runnable end-to-end pipe
//! that produces a far-field radiation pattern. The driver itself contains
//! no new physics; it exists so callers (CLI, Python bindings, examples)
//! can ask for a radiation pattern without having to assemble the inner
//! state objects by hand.
//!
//! ## Source model
//!
//! The dipole is modelled as a **soft current source on `E_z`** distributed
//! over `dipole_length_cells` adjacent cells along z, centred on
//! `dipole_center_cells`. The time profile is a sinusoid at
//! `source_freq_hz` ramped on with a **Hann (raised-cosine) window** over
//! the first three periods so the source does not ring the grid.
//!
//! ## NTFF sweep
//!
//! After the time loop completes the driver evaluates the accumulated
//! [`NtffState`] at `θ ∈ [0°, 180°]` in 5° steps with `φ = 0` and returns
//! the result normalized so that `max |E_θ| = 1`. For a z-polarized
//! short dipole the expected analytic pattern is `|E_θ| ∝ sin θ`.

use std::f64::consts::TAU;

use crate::cpml::CpmlParams;
use crate::grid::YeeGrid;
use crate::ntff::{NtffParams, NtffState};
use crate::{FdtdSolver, WalkingSkeletonSolver, update};

/// User-facing configuration for [`FdtdDriver`].
///
/// All cell indices are 0-based and refer to the integer-E node lattice of
/// the underlying [`YeeGrid`]. Frequencies are in Hz.
#[derive(Debug, Clone, Copy)]
pub struct FdtdDriverConfig {
    /// Number of FDTD time steps to run.
    pub n_steps: usize,
    /// `(i, j, k)` of the dipole centre cell.
    pub dipole_center_cells: (usize, usize, usize),
    /// Length of the dipole **in cells along z**. The current is injected
    /// on `dipole_length_cells` cells of `E_z` centred on
    /// `dipole_center_cells.2`. A value of `1` reduces to a point source.
    pub dipole_length_cells: usize,
    /// Sinusoid drive frequency (Hz). The NTFF DFT bin tracks the same
    /// frequency.
    pub source_freq_hz: f64,
    /// Distance from the inner edge of the CPML to the NTFF integration
    /// surface, in cells. The actual `NtffParams::box_margin_cells` value
    /// passed downstream is `cpml_thickness_cells + ntff_surface_pad_cells`.
    pub ntff_surface_pad_cells: usize,
    /// CPML thickness on every face, in cells. Typical value 10.
    pub cpml_thickness_cells: usize,
}

/// End-to-end FDTD driver: grid + CPML + dipole source + NTFF in one.
///
/// Use [`FdtdDriver::new`] to construct and [`FdtdDriver::run`] to step
/// the simulation to completion and return the far-field
/// [`RadiationPattern`].
pub struct FdtdDriver {
    solver: WalkingSkeletonSolver,
    ntff: NtffState,
    cfg: FdtdDriverConfig,
}

/// Result of [`FdtdDriver::run`]: a θ-cut of the magnitude of `E_θ` at
/// `φ = 0`.
///
/// `theta_deg[i]` corresponds to `e_theta_phi0[i]`. The vector is
/// normalized so that `max e_theta_phi0 == 1.0`, so callers can compare
/// it directly against any analytic pattern that is also normalized to
/// its maximum (e.g. `sin θ`).
#[derive(Debug, Clone)]
pub struct RadiationPattern {
    /// Polar angles in **degrees**, sampled from 0 to 180 inclusive.
    pub theta_deg: Vec<f64>,
    /// `|E_θ|` at each `theta_deg[i]`, with `φ = 0`, normalized so the
    /// maximum across the cut equals `1.0`.
    pub e_theta_phi0: Vec<f64>,
}

impl FdtdDriver {
    /// Build a driver wrapping `grid` with the configuration in `cfg`.
    ///
    /// The CPML is sized to `cfg.cpml_thickness_cells` on every face and
    /// the NTFF box margin is
    /// `cfg.cpml_thickness_cells + cfg.ntff_surface_pad_cells`.
    ///
    /// # Panics
    ///
    /// Panics if `cfg.source_freq_hz` is non-positive, or if the grid is
    /// too small to host the CPML + NTFF box + source.
    pub fn new(grid: YeeGrid, cfg: FdtdDriverConfig) -> Self {
        assert!(
            cfg.source_freq_hz.is_finite() && cfg.source_freq_hz > 0.0,
            "source_freq_hz must be positive and finite"
        );
        assert!(
            cfg.dipole_length_cells >= 1,
            "dipole_length_cells must be ≥ 1"
        );

        let cpml = CpmlParams::for_grid(&grid, cfg.cpml_thickness_cells);
        let box_margin = cfg.cpml_thickness_cells + cfg.ntff_surface_pad_cells;
        let ntff_params = NtffParams {
            f_probe: cfg.source_freq_hz,
            box_margin_cells: box_margin,
            // Placeholder; the run() method sweeps θ via far_field_at.
            theta_rad: std::f64::consts::FRAC_PI_2,
            phi_rad: 0.0,
        };
        let solver = WalkingSkeletonSolver::with_cpml(grid, cpml);
        let ntff = NtffState::new(solver.grid(), ntff_params);
        Self { solver, ntff, cfg }
    }

    /// Step the FDTD time loop, accumulating NTFF surface currents at
    /// `cfg.source_freq_hz`, then sweep the far field over
    /// `θ ∈ [0°, 180°]` in 5° steps at `φ = 0`.
    ///
    /// Returns the per-θ magnitude of `E_θ` normalized so the maximum
    /// equals `1.0`.
    pub fn run(mut self) -> RadiationPattern {
        let omega = TAU * self.cfg.source_freq_hz;
        let period = 1.0 / self.cfg.source_freq_hz;
        // Hann window over the first 3 periods so the source ramps on
        // smoothly. After the ramp the source is a steady sinusoid.
        let ramp_duration = 3.0 * period;

        let (ci, cj, ck) = self.cfg.dipole_center_cells;
        let len = self.cfg.dipole_length_cells;
        // Build the z-indices covered by the dipole, centred on ck.
        // For odd length we get a symmetric range [ck − (len−1)/2, …,
        // ck + (len−1)/2]; for even length the centre falls between two
        // cells and we bias up.
        let half = len / 2;
        let k_lo = ck.saturating_sub(half);
        let k_hi = (ck + (len - half)).min(self.solver.grid().nz);
        let z_indices: Vec<usize> = (k_lo..k_hi).collect();

        for _ in 0..self.cfg.n_steps {
            let t = self.solver.current_time();
            let window = hann_ramp(t, ramp_duration);
            let drive = window * (omega * t).sin();
            step_with_dipole(&mut self.solver, &mut self.ntff, ci, cj, &z_indices, drive);
        }

        // Sweep θ ∈ [0°, 180°] in 5° steps at φ = 0. NtffState::far_field_at
        // is cheap (it re-projects the already-accumulated DFT-bin
        // currents) so 37 calls cost a negligible fraction of the total
        // run time.
        const STEP_DEG: f64 = 5.0;
        let n = (180.0 / STEP_DEG) as usize + 1;
        let mut theta_deg = Vec::with_capacity(n);
        let mut mags = Vec::with_capacity(n);
        for i in 0..n {
            let deg = (i as f64) * STEP_DEG;
            let theta_rad = deg.to_radians();
            let e = self.ntff.far_field_at(theta_rad, 0.0);
            theta_deg.push(deg);
            mags.push(e.norm());
        }

        // Normalize to max so callers can compare against any analytic
        // pattern (e.g. sin θ) that is also normalized to its maximum.
        let m = mags.iter().cloned().fold(0.0f64, f64::max);
        let e_theta_phi0 = if m > 0.0 {
            mags.iter().map(|v| v / m).collect()
        } else {
            mags
        };

        RadiationPattern {
            theta_deg,
            e_theta_phi0,
        }
    }
}

/// Take one FDTD step on `solver` while injecting `drive` onto every
/// `E_z(i, j, k)` whose `k` is listed in `z_indices`. The dipole occupies
/// a short line of `E_z` cells along z.
///
/// Mirrors the structure of
/// [`crate::WalkingSkeletonSolver::step_with_source_and_ntff`] but injects
/// a custom z-line current rather than a Gaussian point source. The
/// source is added between the H and E updates so the next E update
/// sees it through the standard leapfrog timing.
fn step_with_dipole(
    solver: &mut WalkingSkeletonSolver,
    ntff: &mut NtffState,
    i: usize,
    j: usize,
    z_indices: &[usize],
    drive: f64,
) {
    // 1. H update.
    {
        let (grid, _) = solver.grid_and_cpml_mut();
        update::update_h(grid);
    }
    // 2. CPML correction for H (or PEC fallback if no CPML).
    {
        let (grid, cpml) = solver.grid_and_cpml_mut();
        if let Some(cpml) = cpml {
            cpml.update_h(grid);
        } else {
            #[allow(deprecated)]
            crate::boundary::apply_pec(grid);
        }
    }

    // 3. Inject J_z across the dipole cells (soft current source on E_z,
    // applied between H and E updates so the next E update sees it).
    {
        let (grid, _) = solver.grid_and_cpml_mut();
        for &k in z_indices {
            grid.ez[(i, j, k)] += drive;
        }
    }

    // 4. E update.
    {
        let (grid, _) = solver.grid_and_cpml_mut();
        update::update_e(grid);
    }
    // 5. CPML correction for E.
    {
        let (grid, cpml) = solver.grid_and_cpml_mut();
        if let Some(cpml) = cpml {
            cpml.update_e(grid);
        } else {
            #[allow(deprecated)]
            crate::boundary::apply_pec(grid);
        }
    }

    // 6. Advance clock + sample NTFF at the end-of-step time.
    solver.advance_clock();
    let t_after = solver.current_time();
    ntff.sample(solver.grid(), t_after);
}

/// Hann (raised-cosine) ramp from 0 to 1 over `ramp_duration` seconds.
///
/// Returns `1.0` for `t ≥ ramp_duration` and
/// `0.5 · (1 − cos(π · t / ramp_duration))` during the ramp.
#[inline]
fn hann_ramp(t: f64, ramp_duration: f64) -> f64 {
    if !ramp_duration.is_finite() || ramp_duration <= 0.0 {
        return 1.0;
    }
    if t >= ramp_duration {
        return 1.0;
    }
    if t <= 0.0 {
        return 0.0;
    }
    0.5 * (1.0 - (std::f64::consts::PI * t / ramp_duration).cos())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hann_ramp_endpoints() {
        let d = 1.0;
        assert!((hann_ramp(-0.5, d) - 0.0).abs() < 1e-15);
        assert!((hann_ramp(0.0, d) - 0.0).abs() < 1e-15);
        assert!((hann_ramp(0.5, d) - 0.5).abs() < 1e-12);
        assert!((hann_ramp(1.0, d) - 1.0).abs() < 1e-15);
        assert!((hann_ramp(2.0, d) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn radiation_pattern_sweeps_zero_to_180_in_5_deg() {
        let grid = YeeGrid::vacuum(40, 40, 40, 5.0e-3);
        let cfg = FdtdDriverConfig {
            n_steps: 5,
            dipole_center_cells: (20, 20, 20),
            dipole_length_cells: 3,
            source_freq_hz: 1.0e9,
            ntff_surface_pad_cells: 2,
            cpml_thickness_cells: 8,
        };
        let pat = FdtdDriver::new(grid, cfg).run();
        assert_eq!(pat.theta_deg.len(), 37);
        assert_eq!(pat.e_theta_phi0.len(), 37);
        assert!((pat.theta_deg[0] - 0.0).abs() < 1e-12);
        assert!((pat.theta_deg[36] - 180.0).abs() < 1e-12);
    }
}
