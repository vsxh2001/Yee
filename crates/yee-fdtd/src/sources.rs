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
/// Injects a normally-incident plane wave (Phase 2.fdtd.5 only supports
/// `+x` direction with `E_z` polarization and `H_y` carrier) into the
/// total-field region defined by an axis-aligned box of cell indices
/// `[i0..=i1, j0..=j1, k0..=k1]`.
///
/// # Field convention
///
/// - Inside the TF box, stored `E` and `H` are **total** fields.
/// - Outside the TF box, stored `E` and `H` are **scattered** fields.
///
/// Coupling between the regions is implemented as discrete corrections
/// on the `i = i0` and `i = i1` faces, derived from an auxiliary 1-D
/// FDTD incident-field grid that propagates the analytical plane wave
/// with the same numerical dispersion the 3D scheme sees along the
/// propagation axis.
///
/// # Polarization and supported geometry (Phase 2.fdtd.5 / 2.fdtd.5.1)
///
/// For `+x` propagation, `E_z` polarized, the only non-zero incident
/// field components are `E_inc_z(x, t)` and `H_inc_y(x, t)` — incident
/// `H_x`, `H_z`, `E_x`, `E_y` are all identically zero. The discrete
/// Yee stencils that pick up a non-zero incident contribution across
/// the TF/SF boundary are therefore only:
///
/// - `E_z` update at `i = i0` and `i = i1` — uses `H_inc_y` across
///   the `i`-face in `∂H_y/∂x`. **Correction applied here.**
/// - `H_y` update at `i = i0 - 1` and `i = i1` — uses `E_inc_z`
///   across the `i`-face in `∂E_z/∂x`. **Correction applied here.**
///
/// j-face and k-face stencils that *also* cross the TF/SF boundary
/// involve only `H_inc_x` / `H_inc_z` / `E_inc_x` / `E_inc_y`, all of
/// which are zero — so there is no "incident contribution" the
/// correction kernel needs to subtract.
///
/// **However**, the j- and k-face `E_z` discontinuities (TF inside,
/// SF outside) drive the discrete `H_x` / `H_z` updates at those
/// faces and produce a *scattered* field that leaks into the SF
/// region. For slab geometry (`j0 = 0`, `j1 = ny`, `k0 = 0`,
/// `k1 = nz`) those faces sit in CPML and the leakage is absorbed
/// (slab contrast ≈ 2676×, ~68 dB). For a **finite** TF box
/// (smaller than the grid in `y` and / or `z`), those faces are
/// interior and the leakage shows up in the SF region as ~15 dB of
/// residual amplitude (finite-box contrast ≈ 6×). See
/// `tests/plane_wave_finite_box.rs` for the empirical pin.
///
/// **Recommendation:** for high-fidelity TF/SF runs, use slab
/// geometry (the j and k faces in CPML). Finite-box TF/SF for `+x`
/// `E_z` polarization is supported by the existing kernel but with
/// the degraded contrast above; the missing j/k-face *scattered-field*
/// corrections land in Phase 2.fdtd.5.2 / 2.fdtd.5.3 along with
/// oblique-incidence and arbitrary-polarization support.
///
/// # Reference
///
/// Taflove & Hagness, *Computational Electrodynamics* (3rd ed.) §6 and §14.
///
/// # Phase 2.fdtd.5 / 2.fdtd.5.1 limitations
///
/// - Only `PlusX` direction with `E_z` polarization is implemented;
///   other [`PlaneWaveDirection`] variants `unimplemented!()` in the
///   correction kernels.
/// - Only the `i0` / `i1` faces apply corrections. The j- and k-face
///   scattered-field leakage (described above) limits finite-box
///   contrast to ~6× vs the slab's ~2676×. Use slab geometry
///   (`j0 = 0`, `j1 = ny`, `k0 = 0`, `k1 = nz`) for high-fidelity
///   runs; the proper j/k face corrections land in
///   Phase 2.fdtd.5.2 / 2.fdtd.5.3.
/// - The 1-D auxiliary grid uses the same `dx` and `dt` as the 3D grid;
///   for normal incidence this is exact in the limit of the 3D cubic
///   Yee dispersion relation on-axis, but introduces a small mismatch
///   at finite resolution that is well within the `> 10×` TF/SF
///   contrast gate.
/// - The 1-D far-end uses a first-order Mur ABC, sufficient for runs
///   of several hundred steps without spurious 1-D reflections
///   leaking back into the TF region.
#[derive(Debug, Clone)]
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
    /// Call **after** [`crate::update::update_h`] and **after**
    /// [`Self::step_incident_h`] (which advances `H_inc` from the
    /// current `E_inc`).
    ///
    /// # Phase 2.fdtd.5 scope
    ///
    /// Implements `+x` propagation, `E_z` polarized only. Other variants
    /// of [`PlaneWaveDirection`] call `unimplemented!()`.
    pub fn correct_h(&self, grid: &mut YeeGrid) {
        match self.direction {
            PlaneWaveDirection::PlusX => self.correct_h_plus_x(grid),
            _ => unimplemented!(
                "PlaneWaveDirection::{:?} is not implemented in Phase 2.fdtd.5",
                self.direction
            ),
        }
    }

    /// Apply TF/SF corrections to the electric field on the box faces.
    /// Call **after** [`crate::update::update_e`] and **after**
    /// [`Self::step_incident_e`].
    ///
    /// # Phase 2.fdtd.5 scope
    ///
    /// Implements `+x` propagation, `E_z` polarized only. Other variants
    /// of [`PlaneWaveDirection`] call `unimplemented!()`.
    pub fn correct_e(&self, grid: &mut YeeGrid) {
        match self.direction {
            PlaneWaveDirection::PlusX => self.correct_e_plus_x(grid),
            _ => unimplemented!(
                "PlaneWaveDirection::{:?} is not implemented in Phase 2.fdtd.5",
                self.direction
            ),
        }
    }

    /// Map a 3D x-index `i` to a 1-D incident-grid `E_inc` index.
    #[inline]
    fn e_idx(&self, i: usize) -> usize {
        i - self.i0 + self.pad
    }

    /// Map a 3D H_y i-index to a 1-D incident-grid `H_inc` index.
    /// `H_y[i, *, *]` lives at the half-cell `(i + 1/2, *, *)`, so its
    /// 1-D counterpart is `H_inc[i - i0 + pad]`.
    #[inline]
    fn h_idx(&self, i_h: usize) -> usize {
        i_h - self.i0 + self.pad
    }

    // ----------------------------------------------------------------
    // +x propagation, E_z polarization (Phase 2.fdtd.5 / 2.fdtd.5.1)
    //
    // Derivation (Taflove & Hagness §14):
    //
    // For a +x plane wave with E along z and H along y, the incident
    // field has only E_inc_z and H_inc_y non-zero. The only Yee
    // stencils that pick up a non-zero *incident* contribution across
    // the TF/SF boundary are:
    //   - E_z update at i = i0 / i1   (uses H_inc_y across the i-face)
    //   - H_y update at i = i0-1 / i1 (uses E_inc_z across the i-face)
    //
    // The j- and k-face stencils only see zero incident components, so
    // no incident-correction is needed there. They DO, however, see the
    // E_z TF/SF discontinuity and emit a spurious scattered field —
    // see the struct docstring's "Polarization and supported geometry"
    // section and `tests/plane_wave_finite_box.rs` for the empirical
    // measurement (~6× contrast for finite-box vs ~2676× for slab).
    // The j/k scattered-field corrections that recover full contrast
    // are deferred to Phase 2.fdtd.5.2.
    //
    // Front (i = i0):
    //   H_y[i0-1, j, k] is SF (between SF E_z[i0-1] and TF E_z[i0]).
    //   Standard update_h read E_z[i0] (TF) thinking it was SF;
    //   correction: subtract (dt/(μ₀·dx)) · E_inc_z[at i0]  from H_y[i0-1].
    //
    //   E_z[i0, j, k] is TF, but standard update_e read
    //   H_y[i0-1, j, k] (SF) thinking it was TF; correction:
    //   subtract (dt/(ε₀·dx)) · H_inc_y[at i0-1]  from E_z[i0].
    //
    // Back (i = i1):
    //   H_y[i1, j, k] is SF (between TF E_z[i1] and SF E_z[i1+1]).
    //   Standard update_h read E_z[i1] (TF) thinking it was SF;
    //   correction: add (dt/(μ₀·dx)) · E_inc_z[at i1]  to H_y[i1].
    //
    //   E_z[i1, j, k] is TF, but standard update_e read
    //   H_y[i1, j, k] (SF) thinking it was TF; correction:
    //   add (dt/(ε₀·dx)) · H_inc_y[at i1]  to E_z[i1].
    // ----------------------------------------------------------------

    fn correct_h_plus_x(&self, grid: &mut YeeGrid) {
        let coeff = self.dt / (MU0 * self.dx);
        let einc_front = self.inc_e[self.e_idx(self.i0)];
        let einc_back = self.inc_e[self.e_idx(self.i1)];

        // Bounds-check: H_y has shape [nx, ny+1, nz]. Need i0 ≥ 1
        // (so i0-1 is valid) and i1 ≤ nx-1.
        assert!(
            self.i0 >= 1,
            "PlaneWaveSource (+x): i0 must be ≥ 1 (got {})",
            self.i0
        );
        assert!(
            self.i1 < grid.nx,
            "PlaneWaveSource (+x): i1 ({}) must be < grid.nx ({})",
            self.i1,
            grid.nx
        );

        let k_hi = self.k1.min(grid.nz);
        for j in self.j0..=self.j1 {
            for k in self.k0..k_hi {
                grid.hy[(self.i0 - 1, j, k)] -= coeff * einc_front;
            }
        }
        for j in self.j0..=self.j1 {
            for k in self.k0..k_hi {
                grid.hy[(self.i1, j, k)] += coeff * einc_back;
            }
        }
    }

    fn correct_e_plus_x(&self, grid: &mut YeeGrid) {
        let coeff = self.dt / (EPS0 * self.dx);
        assert!(
            self.h_idx(self.i0 - 1) < self.inc_h.len(),
            "PlaneWaveSource (+x): h_idx out of range (logic bug, please report)"
        );
        let hinc_front = self.inc_h[self.h_idx(self.i0 - 1)];
        let hinc_back = self.inc_h[self.h_idx(self.i1)];

        let k_hi = self.k1.min(grid.nz);
        for j in self.j0..=self.j1 {
            for k in self.k0..k_hi {
                grid.ez[(self.i0, j, k)] -= coeff * hinc_front;
            }
        }
        for j in self.j0..=self.j1 {
            for k in self.k0..k_hi {
                grid.ez[(self.i1, j, k)] += coeff * hinc_back;
            }
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
