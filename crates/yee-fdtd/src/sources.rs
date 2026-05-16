//! Source helpers for the FDTD walking skeleton.
//!
//! Phase 2.0 shipped a single point-source primitive: a Gaussian-in-time pulse
//! added (soft source) to a chosen cell of `E_z`. Phase 2.fdtd.5 adds a
//! total-field / scattered-field (TF/SF) plane-wave source, see
//! [`PlaneWaveSource`].
//!
//! Hard sources, modal sources, and lumped ports remain Phase 2.1+ work.

use std::f64::consts::TAU;

use yee_core::units::{EPS0, MU0};

use crate::grid::YeeGrid;

/// Cardinal-axis propagation direction for [`PlaneWaveSource`].
///
/// Phase 2.fdtd.5 only implements [`PlaneWaveDirection::PlusX`] (E_z
/// polarized). The other variants are recognized by the constructor but
/// cause [`PlaneWaveSource::correct_h`] and
/// [`PlaneWaveSource::correct_e`] to `unimplemented!()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneWaveDirection {
    /// Propagation along `+x` (E_z polarized in Phase 2.fdtd.5).
    PlusX,
    /// Propagation along `+y` — not implemented in Phase 2.fdtd.5.
    PlusY,
    /// Propagation along `+z` — not implemented in Phase 2.fdtd.5.
    PlusZ,
    /// Propagation along `-x` — not implemented in Phase 2.fdtd.5.
    MinusX,
    /// Propagation along `-y` — not implemented in Phase 2.fdtd.5.
    MinusY,
    /// Propagation along `-z` — not implemented in Phase 2.fdtd.5.
    MinusZ,
}

/// Total-field / scattered-field (TF/SF) plane-wave source.
///
/// Skeleton struct for Phase 2.fdtd.5: holds the TF region bounds, the
/// propagation direction, the source-pulse shape parameters, and an
/// auxiliary 1-D incident-field grid that will propagate the analytical
/// plane wave with the same numerical dispersion the 3D scheme sees on
/// the propagation axis.
///
/// This commit lands only the data layout and a stubbed [`Self::correct_h`]
/// / [`Self::correct_e`] / [`Self::step_incident_h`] /
/// [`Self::step_incident_e`] API surface so callers can wire up
/// [`crate::WalkingSkeletonSolver::step_with_plane_wave`] later in this
/// phase. The actual 1-D kernel and TF/SF coupling-correction math land
/// in follow-up commits.
///
/// # Reference
///
/// Taflove & Hagness, *Computational Electrodynamics* (3rd ed.) §6 and §14.
#[derive(Debug, Clone)]
#[allow(dead_code)] // i0..k1 are pub fields; correction kernels use them in the next commit.
pub struct PlaneWaveSource {
    /// TF region lower x cell index (inclusive).
    pub i0: usize,
    /// TF region upper x cell index (inclusive).
    pub i1: usize,
    /// TF region lower y cell index (inclusive).
    pub j0: usize,
    /// TF region upper y cell index (inclusive).
    pub j1: usize,
    /// TF region lower z cell index (inclusive).
    pub k0: usize,
    /// TF region upper z cell index (inclusive).
    pub k1: usize,
    /// Propagation direction.
    pub direction: PlaneWaveDirection,
    /// Source carrier frequency (Hz).
    pub frequency: f64,
    /// Hanning-window taper length, in time steps.
    pub ramp_steps: usize,

    /// 1-D auxiliary grid: `E_inc` samples. Length = `(i1 - i0) + 2*pad + 1`
    /// for `PlusX`. Index 0 is the source-injection cell; index `pad`
    /// corresponds to the 3D plane `i = i0`.
    inc_e: Vec<f64>,
    /// 1-D auxiliary grid: `H_inc` samples. Length = `inc_e.len() - 1`,
    /// staggered half a cell to the right of each `inc_e` sample.
    inc_h: Vec<f64>,
    /// Number of "lead-in" cells in the 1-D grid before the TF front face.
    pad: usize,
    /// Cell size of the 3D grid along the propagation axis (cached for
    /// incident-grid updates and corrections).
    dx: f64,
    /// Time step of the 3D grid.
    dt: f64,
    /// 1-D incident grid step counter.
    step: usize,
    /// Previous-step value of `inc_e[N - 1]` (far-end cell), used by the
    /// first-order Mur ABC on the 1-D grid.
    mur_prev_end: f64,
    /// Previous-step value of `inc_e[N - 2]` (cell just inside the
    /// far end), used by the first-order Mur ABC.
    mur_prev_inner: f64,
    /// Mur ABC coefficient `(c·dt - dx)/(c·dt + dx)`, cached.
    mur_coeff: f64,
}

impl PlaneWaveSource {
    /// Build a new TF/SF plane-wave source for the given TF region.
    ///
    /// `pad` controls the number of "lead-in" cells in the 1-D auxiliary
    /// grid before the TF front face. A value of `4` is the documented
    /// minimum: the source pulse needs at least a few cells to develop
    /// before its leading edge reaches the TF boundary.
    ///
    /// # Panics
    ///
    /// Panics if any of the region bounds are inverted (`i0 > i1`, etc.)
    /// or if `pad < 1`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        i0: usize,
        i1: usize,
        j0: usize,
        j1: usize,
        k0: usize,
        k1: usize,
        direction: PlaneWaveDirection,
        frequency: f64,
        ramp_steps: usize,
        dx: f64,
        dt: f64,
        pad: usize,
    ) -> Self {
        assert!(i0 <= i1, "PlaneWaveSource: i0 ({i0}) must be ≤ i1 ({i1})");
        assert!(j0 <= j1, "PlaneWaveSource: j0 ({j0}) must be ≤ j1 ({j1})");
        assert!(k0 <= k1, "PlaneWaveSource: k0 ({k0}) must be ≤ k1 ({k1})");
        assert!(pad >= 1, "PlaneWaveSource: pad ({pad}) must be ≥ 1");
        assert!(
            frequency > 0.0 && frequency.is_finite(),
            "PlaneWaveSource: frequency must be positive and finite"
        );
        assert!(
            dx > 0.0 && dx.is_finite(),
            "PlaneWaveSource: dx must be positive and finite"
        );
        assert!(
            dt > 0.0 && dt.is_finite(),
            "PlaneWaveSource: dt must be positive and finite"
        );

        let n_along = match direction {
            PlaneWaveDirection::PlusX | PlaneWaveDirection::MinusX => i1 - i0,
            PlaneWaveDirection::PlusY | PlaneWaveDirection::MinusY => j1 - j0,
            PlaneWaveDirection::PlusZ | PlaneWaveDirection::MinusZ => k1 - k0,
        };
        let inc_n_cells = n_along + 2 * pad + 1;
        let inc_e = vec![0.0; inc_n_cells];
        let inc_h = vec![0.0; inc_n_cells - 1];

        let c0 = yee_core::units::C0;
        let mur_coeff = (c0 * dt - dx) / (c0 * dt + dx);

        Self {
            i0,
            i1,
            j0,
            j1,
            k0,
            k1,
            direction,
            frequency,
            ramp_steps,
            inc_e,
            inc_h,
            pad,
            dx,
            dt,
            step: 0,
            mur_prev_end: 0.0,
            mur_prev_inner: 0.0,
            mur_coeff,
        }
    }

    /// Hanning (raised-cosine) ramp factor for a sinusoidal source.
    ///
    /// Returns `0.5 * (1 - cos(π · n / ramp_steps))` for `n < ramp_steps`
    /// and `1.0` afterwards. Tapering the carrier on with a Hann window
    /// suppresses the broadband click an unramped sinusoid would inject.
    fn ramp(&self) -> f64 {
        if self.ramp_steps == 0 || self.step >= self.ramp_steps {
            1.0
        } else {
            0.5 * (1.0 - (std::f64::consts::PI * self.step as f64 / self.ramp_steps as f64).cos())
        }
    }

    /// Drive value of the source at the current 1-D-grid step: a sinusoid
    /// `sin(2π f n dt)` modulated by the Hann ramp.
    fn source_value(&self) -> f64 {
        let t = self.step as f64 * self.dt;
        self.ramp() * (TAU * self.frequency * t).sin()
    }

    /// Advance `H_inc` by one time step.
    ///
    /// Standard 1-D Yee H-update:
    /// ```text
    /// H_inc_y[m+1/2] += (Δt/(μ₀ Δx)) · (E_inc_z[m+1] - E_inc_z[m])
    /// ```
    /// matches the sign convention of [`crate::update::update_h`] on
    /// the 3D grid (∂H_y/∂t = +(1/μ) ∂E_z/∂x).
    pub fn step_incident_h(&mut self) {
        let coeff = self.dt / (MU0 * self.dx);
        for m in 0..self.inc_h.len() {
            self.inc_h[m] += coeff * (self.inc_e[m + 1] - self.inc_e[m]);
        }
    }

    /// Update `E_inc`, inject the analytic source at the near end, and
    /// apply a first-order Mur ABC at the far end. See
    /// [`Self::step_incident_h`].
    ///
    /// The leapfrog body is:
    /// ```text
    /// E_inc_z[m]   += (Δt/(ε₀ Δx)) · (H_inc_y[m+1/2] - H_inc_y[m-1/2])
    /// E_inc_z[0]   = ramp(n)·sin(2π f n Δt)              (hard source)
    /// E_inc_z[N-1] = E_inc[N-2]^old + κ·(E_inc[N-2]^new - E_inc[N-1]^old)
    /// ```
    /// where κ = (c·Δt - Δx)/(c·Δt + Δx) is the Mur first-order ABC
    /// coefficient.
    pub fn step_incident_e(&mut self) {
        let coeff = self.dt / (EPS0 * self.dx);
        let n = self.inc_e.len();

        let prev_end = self.mur_prev_end;
        let prev_inner = self.mur_prev_inner;

        // Update E_inc[m] for m ∈ [1, n-1) using the freshly-stepped H_inc.
        for m in 1..n - 1 {
            self.inc_e[m] += coeff * (self.inc_h[m] - self.inc_h[m - 1]);
        }
        // Hard source at m=0.
        self.step += 1;
        self.inc_e[0] = self.source_value();

        // First-order Mur ABC at far end for outgoing +x waves.
        let inner_new = self.inc_e[n - 2];
        self.inc_e[n - 1] = prev_inner + self.mur_coeff * (inner_new - prev_end);

        // Save state for next call.
        self.mur_prev_end = self.inc_e[n - 1];
        self.mur_prev_inner = inner_new;
    }

    /// Apply TF/SF corrections to the magnetic field on the box faces.
    /// Stub in this commit; per-face math lands in the "TF/SF corrections"
    /// commit.
    pub fn correct_h(&self, _grid: &mut YeeGrid) {
        match self.direction {
            PlaneWaveDirection::PlusX => { /* stub */ }
            _ => unimplemented!(
                "PlaneWaveDirection::{:?} is not implemented in Phase 2.fdtd.5",
                self.direction
            ),
        }
    }

    /// Apply TF/SF corrections to the electric field on the box faces.
    /// Stub in this commit.
    pub fn correct_e(&self, _grid: &mut YeeGrid) {
        match self.direction {
            PlaneWaveDirection::PlusX => { /* stub */ }
            _ => unimplemented!(
                "PlaneWaveDirection::{:?} is not implemented in Phase 2.fdtd.5",
                self.direction
            ),
        }
    }

    /// Read access to the auxiliary 1-D incident-E grid (mostly for tests).
    pub fn inc_e(&self) -> &[f64] {
        &self.inc_e
    }

    /// Read access to the auxiliary 1-D incident-H grid (mostly for tests).
    pub fn inc_h(&self) -> &[f64] {
        &self.inc_h
    }

    /// Current 1-D step counter.
    pub fn step_count(&self) -> usize {
        self.step
    }
}

// ---- legacy point-source helpers (Phase 2.0) ----

/// Add a Gaussian-time pulse to `E_z(i, j, k)`.
///
/// The injected value is `exp(-((t - t0) / sigma)²)` (a unit-amplitude soft
/// source). The caller controls the time stepping; this function simply
/// *adds* the source contribution to the existing field value.
///
/// # Panics
///
/// Panics if `(i, j, k)` is outside the bounds of `E_z`
/// (shape `[nx+1, ny+1, nz]`).
pub fn gaussian_pulse_ez(
    grid: &mut YeeGrid,
    i: usize,
    j: usize,
    k: usize,
    t: f64,
    t0: f64,
    sigma: f64,
) {
    assert!(
        sigma > 0.0 && sigma.is_finite(),
        "gaussian sigma must be positive and finite"
    );
    let arg = (t - t0) / sigma;
    let amplitude = (-arg * arg).exp();
    grid.ez[(i, j, k)] += amplitude;
}
