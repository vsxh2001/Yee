//! Source helpers for the FDTD walking skeleton.
//!
//! Phase 2.0 shipped a single point-source primitive: a Gaussian-in-time pulse
//! added (soft source) to a chosen cell of `E_z`. Phase 2.fdtd.5 adds a
//! total-field / scattered-field (TF/SF) plane-wave source, see
//! [`PlaneWaveSource`].
//!
//! Hard sources, modal sources, and lumped ports remain Phase 2.1+ work.
//!
//! ## Phase 2.fdtd.5.2 design notes — j/k-face SF corrections
//!
//! **Assumption being challenged:** Phase 2.fdtd.5.1 shipped only `i`-face
//! TF/SF corrections, on the reasoning that for a `+x` `E_z`-polarized
//! plane wave the incident `H_inc_x`, `H_inc_z`, `E_inc_x`, `E_inc_y` are
//! all zero, so the j- and k-face stencils carry no spurious **incident**
//! contribution. That argument is correct for the *incident* leg but
//! misses the *scattered* leg: the j-face has a TF-vs-SF discontinuity
//! in `E_z` (it equals `E_inc_z` plus scattered inside the box, just
//! scattered outside), and the standard Yee `H_x` update at `j = j0 - 1`
//! and `j = j1` straddles that discontinuity in its `∂E_z / ∂y` term,
//! mixing TF and SF `E_z` and emitting spurious scattered field into the
//! SF region. Symmetrically, the `E_x` update at `k = k0` and
//! `k = k1 + 1` straddles the z-discontinuity in `H_y` via its
//! `∂H_y / ∂z` term. With those four corrections added, the finite-box
//! configuration's TF/SF contrast jumps from ~6× (Phase 2.fdtd.5.1
//! empirical pin) to >100× (Phase 2.fdtd.5.2 target), and slab
//! geometry — which puts the j/k faces inside CPML — is unaffected.
//!
//! **Which curl stencils need a correction for `+x` `E_z` polarization**
//! (only `E_inc_z(x)` and `H_inc_y(x)` are non-zero):
//!
//! - `H_y` curl has `∂E_z / ∂x` — i-face straddle. (5.1, shipped.)
//! - `H_x` curl has `∂E_z / ∂y` — j-face straddle. (5.2, this commit.)
//! - `E_z` curl has `∂H_y / ∂x` — i-face straddle. (5.1, shipped.)
//! - `E_x` curl has `∂H_y / ∂z` — k-face straddle. (5.2, this commit.)
//!
//! All other E/H components' curls involve only zero-incident pairs
//! (`E_inc_x = E_inc_y = H_inc_x = H_inc_z = 0`), so they need no
//! correction for `+x` `E_z` polarization. Arbitrary-polarization /
//! oblique-incidence support lands in Phase 2.fdtd.5.3+.
//!
//! **Approach:** extend the existing i-face apply pattern to j and k
//! faces. Each face/component pair is one inclusive 2-D loop matching
//! the i-face's `j0..=j1, k0..k_hi` index conventions. The sign of the
//! correction comes from the side of the box and the orientation of the
//! straddled discontinuity, derived in the per-face comment blocks
//! below.
//!
//! **Reference:** Taflove & Hagness, *Computational Electrodynamics*
//! (3rd ed.) §5.10 — 3-D TF/SF for a rectangular Huygens surface.
//!
//! **DoD:** the `tests/plane_wave_finite_box.rs` contrast must rise
//! from ~6× to ≥ 100×, and the slab variant
//! (`tests/plane_wave_propagation.rs`) must not regress below its
//! previous ~2676×.

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
/// # Polarization and supported geometry (Phase 2.fdtd.5 / 2.fdtd.5.1 / 2.fdtd.5.2)
///
/// For `+x` propagation, `E_z` polarized, the only non-zero incident
/// field components are `E_inc_z(x, t)` and `H_inc_y(x, t)` — incident
/// `H_x`, `H_z`, `E_x`, `E_y` are all identically zero. The discrete
/// Yee stencils that pick up a non-zero incident contribution across
/// the TF/SF boundary are therefore exactly four:
///
/// - `E_z` update at `i = i0` and `i = i1` — uses `H_inc_y` across
///   the `i`-face in `∂H_y/∂x`. **Correction applied (5.1).**
/// - `H_y` update at `i = i0 - 1` and `i = i1` — uses `E_inc_z`
///   across the `i`-face in `∂E_z/∂x`. **Correction applied (5.1).**
/// - `H_x` update at `j = j0 - 1` and `j = j1` — uses `E_inc_z`
///   across the `j`-face in `∂E_z/∂y`. **Correction applied (5.2).**
/// - `E_x` update at `k = k0` and `k = k1 + 1` — uses `H_inc_y`
///   across the `k`-face in `∂H_y/∂z`. **Correction applied (5.2).**
///
/// Phase 2.fdtd.5.1 shipped only the first two. With the 5.2 j/k-face
/// additions, finite-box geometry (TF box bounded on all six faces)
/// achieves a contrast ratio well above 100×, comparable to the slab
/// configuration. Slab geometry remains the recommended option when
/// the geometry permits, because the slab j/k faces still sit in CPML
/// and avoid even the discretized-correction round-off.
///
/// # Reference
///
/// Taflove & Hagness, *Computational Electrodynamics* (3rd ed.) §5.10
/// (3-D TF/SF for a rectangular Huygens surface) and §6 / §14.
///
/// # Phase 2.fdtd.5 / 2.fdtd.5.1 / 2.fdtd.5.2 limitations
///
/// - Only `PlusX` direction with `E_z` polarization is implemented;
///   other [`PlaneWaveDirection`] variants `unimplemented!()` in the
///   correction kernels.
/// - All four faces (`i0`, `i1`, `j0/j1`, `k0/k1`) now apply
///   corrections for `+x` `E_z` polarization. Arbitrary polarization
///   and oblique incidence land in Phase 2.fdtd.5.3+.
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
    // +x propagation, E_z polarization (Phase 2.fdtd.5 / 2.fdtd.5.1 / 2.fdtd.5.2)
    //
    // Derivation (Taflove & Hagness §5.10 / §14):
    //
    // For a +x plane wave with E along z and H along y, the incident
    // field has only E_inc_z(x) and H_inc_y(x) non-zero. The four Yee
    // stencils that pick up a non-zero incident contribution across the
    // TF/SF boundary are:
    //
    //   - E_z update at i = i0 / i1   (uses H_inc_y across the i-face,
    //                                  in `∂H_y/∂x` term of E_z curl)
    //   - H_y update at i = i0-1 / i1 (uses E_inc_z across the i-face,
    //                                  in `∂E_z/∂x` term of H_y curl)
    //   - H_x update at j = j0-1 / j1 (uses E_inc_z across the j-face,
    //                                  in `∂E_z/∂y` term of H_x curl)
    //   - E_x update at k = k0 / k1+1 (uses H_inc_y across the k-face,
    //                                  in `∂H_y/∂z` term of E_x curl)
    //
    // The first two are i-face corrections (Phase 2.fdtd.5 / 5.1); the
    // last two are j/k-face corrections (Phase 2.fdtd.5.2).
    //
    // ----- i-face -----
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
    //
    // ----- j-face (Phase 2.fdtd.5.2) -----
    //
    // H_x[i, j, k] update is:
    //     H_x += (dt/μ₀) · (∂E_y/∂z − ∂E_z/∂y)
    //          = (dt/μ₀) · ( (E_y[i,j,k+1]−E_y[i,j,k])/dz
    //                        − (E_z[i,j+1,k]−E_z[i,j,k])/dy )
    // Only the `∂E_z/∂y` term can straddle the j-face TF/SF boundary
    // (E_z is the only incident-bearing field; E_y has no incident).
    //
    // Front (j = j0):
    //   H_x[i, j0-1, k] is SF (between SF E_z[i,j0-1,k] and TF
    //   E_z[i,j0,k]). Standard update_h read E_z[i,j0,k] (TF) as if it
    //   were SF; that put a spurious `−(dt/(μ₀·dy))·E_inc_z(i)` into
    //   H_x[i, j0-1, k] (negative because of the minus sign on
    //   `∂E_z/∂y`). Correction:
    //
    //     H_x[i, j0-1, k]  +=  (dt/(μ₀·dy)) · E_inc_z[at i]
    //
    // Back (j = j1):
    //   H_x[i, j1, k] is SF (between TF E_z[i,j1,k] and SF
    //   E_z[i,j1+1,k]). The TF E_z is now the *subtracted* term in
    //   `(E_z[j1+1] − E_z[j1])/dy`, so the spurious contribution is
    //   `+(dt/(μ₀·dy))·E_inc_z(i)`. Correction:
    //
    //     H_x[i, j1, k]    −=  (dt/(μ₀·dy)) · E_inc_z[at i]
    //
    // Both corrections use `E_inc_z` sampled at the x-index of the H_x
    // cell (`H_x[i,*,*]` lives at integer x = i, and E_z[i,*,*] also
    // lives at integer x = i, so they share `inc_e[e_idx(i)]`).
    //
    // No j-face correction is needed for E_z updates: E_z curl has
    // `∂H_x/∂y`, and H_x has no incident component (H_inc_x = 0).
    //
    // ----- k-face (Phase 2.fdtd.5.2) -----
    //
    // E_x[i, j, k] update is:
    //     E_x += (dt/ε₀) · (∂H_z/∂y − ∂H_y/∂z)
    //          = (dt/ε₀) · ( (H_z[i,j,k]−H_z[i,j-1,k])/dy
    //                        − (H_y[i,j,k]−H_y[i,j,k-1])/dz )
    // Only the `∂H_y/∂z` term can straddle the k-face TF/SF boundary
    // (H_y is the only incident-bearing field; H_z has no incident).
    //
    // Front (k = k0):
    //   E_x[i, j, k0] lives at z = k0 (integer plane), on the boundary
    //   between SF H_y[i,j,k0-1] (z = k0-1/2) and TF H_y[i,j,k0]
    //   (z = k0+1/2). By the same convention used on the i-face — E
    //   nodes on the boundary plane are claimed as TF — E_x[i,j,k0] is
    //   TF. The standard update read H_y[i,j,k0-1] (SF) thinking it was
    //   TF, so it under-read by `−H_inc_y(i+1/2)`; that propagated into
    //   `∂H_y/∂z` as `−(−H_inc_y)/dz = +H_inc_y/dz`, and into the E_x
    //   update with the `−∂H_y/∂z` sign as `−(dt/(ε₀·dz))·H_inc_y(i+1/2)`.
    //   Correction:
    //
    //     E_x[i, j, k0]    +=  (dt/(ε₀·dz)) · H_inc_y[at i+1/2]
    //
    // Back (k = k1 + 1):
    //   E_x[i, j, k1+1] lives at z = k1+1 (integer plane), on the
    //   boundary between TF H_y[i,j,k1] (z = k1+1/2) and SF
    //   H_y[i,j,k1+1] (z = k1+3/2). E_x at this plane is TF by
    //   convention. Standard update read H_y[i,j,k1+1] (SF) as if TF,
    //   under-reading by `−H_inc_y(i+1/2)`; that propagated through
    //   `∂H_y/∂z` and `−∂H_y/∂z` as `+(dt/(ε₀·dz))·H_inc_y(i+1/2)` of
    //   spurious E_x. Correction:
    //
    //     E_x[i, j, k1+1]  −=  (dt/(ε₀·dz)) · H_inc_y[at i+1/2]
    //
    // Both corrections use `H_inc_y` sampled at the x-coord of the
    // E_x cell. E_x[i,*,*] lives at x = i+1/2; H_y[i,*,*] also lives
    // at x = i+1/2; so the 1-D-grid lookup is `inc_h[h_idx(i)]`.
    //
    // No k-face correction is needed for E_z updates (E_z curl has no
    // z-derivative) or for H_x / H_y updates whose z-derivative terms
    // involve E_x or E_y (no incident).
    // ----------------------------------------------------------------

    fn correct_h_plus_x(&self, grid: &mut YeeGrid) {
        self.correct_h_iface_plus_x(grid);
        self.correct_h_jface_plus_x(grid);
    }

    fn correct_e_plus_x(&self, grid: &mut YeeGrid) {
        self.correct_e_iface_plus_x(grid);
        self.correct_e_kface_plus_x(grid);
    }

    fn correct_h_iface_plus_x(&self, grid: &mut YeeGrid) {
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

        // H_y cross-section: all (j, k) where E_z[i0, j, k] is TF, i.e.
        // j ∈ [j0, j1] and k ∈ [k0, k1]. H_y has shape [nx, ny+1, nz],
        // so j up to ny is valid and k up to nz-1 is valid; clamp k1.
        // (Phase 2.fdtd.5.2: switched k from exclusive `k0..k_hi` to
        // inclusive `k0..=k1.min(nz-1)` so the upper-z TF slice gets
        // corrected too; slab geometry — where k1 = nz — is unchanged
        // because min(nz, nz-1) = nz-1 = the original `k_hi - 1`.)
        let k_hi = self.k1.min(grid.nz.saturating_sub(1));
        for j in self.j0..=self.j1.min(grid.ny) {
            for k in self.k0..=k_hi {
                grid.hy[(self.i0 - 1, j, k)] -= coeff * einc_front;
                grid.hy[(self.i1, j, k)] += coeff * einc_back;
            }
        }
    }

    /// Apply the j-face H_x corrections (Phase 2.fdtd.5.2).
    ///
    /// Cancels the spurious `E_inc_z` contribution picked up by the
    /// standard `update_h` `∂E_z/∂y` stencil at `H_x[i, j0-1, k]`
    /// (front face) and `H_x[i, j1, k]` (back face).
    ///
    /// `H_x` has shape `[nx+1, ny, nz]`. The j-face correction is a
    /// no-op when `j0 == 0` (no SF row at `j = j0 - 1` to correct) or
    /// when `j1 >= ny` (no SF row at `j = j1`); both situations
    /// correspond to slab geometry where the j-face sits in CPML.
    fn correct_h_jface_plus_x(&self, grid: &mut YeeGrid) {
        // H_x dy-coefficient. `grid.dy` is the relevant cell size for
        // the `∂E_z/∂y` stencil; the 3D grid is cubic in the walking
        // skeleton (dx = dy = dz), but we use grid.dy explicitly for
        // forward compatibility with non-cubic cells.
        let coeff = self.dt / (MU0 * grid.dy);

        // H_x has shape [nx+1, ny, nz]. The cross-section we correct
        // covers (i, k) ∈ [i0, i1] × [k0, k1]; clamp to valid H_x
        // indices.
        let i_hi = self.i1.min(grid.nx);
        let k_hi = self.k1.min(grid.nz.saturating_sub(1));

        // Front j-face: SF H_x row at j = j0 - 1. Skip when j0 == 0
        // (slab in y — the j-face sits in CPML, no correction needed).
        if self.j0 >= 1 {
            for i in self.i0..=i_hi {
                let einc = self.inc_e[self.e_idx(i)];
                for k in self.k0..=k_hi {
                    grid.hx[(i, self.j0 - 1, k)] += coeff * einc;
                }
            }
        }

        // Back j-face: SF H_x row at j = j1. Skip when j1 >= ny (slab
        // in y — the j-face row at j = j1 is past the H_x j-range).
        if self.j1 < grid.ny {
            for i in self.i0..=i_hi {
                let einc = self.inc_e[self.e_idx(i)];
                for k in self.k0..=k_hi {
                    grid.hx[(i, self.j1, k)] -= coeff * einc;
                }
            }
        }
    }

    fn correct_e_iface_plus_x(&self, grid: &mut YeeGrid) {
        let coeff = self.dt / (EPS0 * self.dx);
        assert!(
            self.h_idx(self.i0 - 1) < self.inc_h.len(),
            "PlaneWaveSource (+x): h_idx out of range (logic bug, please report)"
        );
        let hinc_front = self.inc_h[self.h_idx(self.i0 - 1)];
        let hinc_back = self.inc_h[self.h_idx(self.i1)];

        // E_z cross-section at i = i0 / i1: all (j, k) where this E_z
        // is itself TF. E_z has shape [nx+1, ny+1, nz], so j up to ny
        // and k up to nz-1 are valid. (See `correct_h_iface_plus_x`
        // for the Phase 2.fdtd.5.2 inclusive-k rationale.)
        let k_hi = self.k1.min(grid.nz.saturating_sub(1));
        for j in self.j0..=self.j1.min(grid.ny) {
            for k in self.k0..=k_hi {
                grid.ez[(self.i0, j, k)] -= coeff * hinc_front;
                grid.ez[(self.i1, j, k)] += coeff * hinc_back;
            }
        }
    }

    /// Apply the k-face E_x corrections (Phase 2.fdtd.5.2).
    ///
    /// Cancels the spurious `H_inc_y` contribution picked up by the
    /// standard `update_e` `∂H_y/∂z` stencil at `E_x[i, j, k0]`
    /// (front face) and `E_x[i, j, k1+1]` (back face).
    ///
    /// `E_x` has shape `[nx, ny+1, nz+1]`. The k-face correction is a
    /// no-op when `k0 == 0` (no E_x row above the boundary in z — the
    /// k=0 face is PEC/CPML) or when `k1 + 1 > nz`; both situations
    /// correspond to slab geometry where the k-face sits in CPML.
    fn correct_e_kface_plus_x(&self, grid: &mut YeeGrid) {
        // E_x dz-coefficient.
        let coeff = self.dt / (EPS0 * grid.dz);

        // E_x cross-section to correct: the (i, j) cells where the
        // standard `update_e` `∂H_y/∂z` stencil straddles the k-face.
        // The straddle exists only where H_y on the TF side of the
        // face is itself TF. By the i-face convention, H_y is TF for
        // i ∈ [i0, i1-1] (i1 is the SF back-boundary index for H_y);
        // and TF for j ∈ [j0, j1]. The i-range is therefore one cell
        // **narrower** than the H_x j-face cross-section (because of
        // the half-cell offset between H_y's i-coordinate (x = i+1/2)
        // and E_z's i-coordinate (x = i, integer)). Clamp to valid
        // E_x bounds; `E_x` has shape `[nx, ny+1, nz+1]`.
        let i_hi = self.i1.saturating_sub(1).min(grid.nx.saturating_sub(1));
        let j_hi = self.j1.min(grid.ny);

        // Front k-face: TF E_x slab at k = k0. Skip when k0 == 0
        // because then there is no SF H_y row at k = k0 - 1 (the
        // boundary sits at the grid edge — CPML territory).
        if self.k0 >= 1 && self.i0 <= i_hi {
            for i in self.i0..=i_hi {
                let hinc = self.inc_h[self.h_idx(i)];
                for j in self.j0..=j_hi {
                    grid.ex[(i, j, self.k0)] += coeff * hinc;
                }
            }
        }

        // Back k-face: TF E_x slab at k = k1 + 1. Skip when k1+1 > nz
        // because then there is no E_x row at that k (slab geometry).
        if self.k1 < grid.nz && self.i0 <= i_hi {
            for i in self.i0..=i_hi {
                let hinc = self.inc_h[self.h_idx(i)];
                for j in self.j0..=j_hi {
                    grid.ex[(i, j, self.k1 + 1)] -= coeff * hinc;
                }
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
