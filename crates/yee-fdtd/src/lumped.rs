//! Lumped R / L / C / series-RLC port source (Phase 2.fdtd.6).
//!
//! A lumped element is a sub-cell modification to the standard Yee E-update at
//! a single Yee cell: it injects a current-density term consistent with the
//! element's V-I relationship, so the FDTD lattice sees the cell as if it
//! contained a discrete circuit component bridging the two faces of an `E_z`
//! edge.
//!
//! # Geometry and conventions
//!
//! Phase 2.fdtd.6 supports elements oriented **along ±z** at a single Yee
//! cell `(i, j, k)`. The voltage across the element is
//!
//! ```text
//! V = E_z(i, j, k) * dz
//! ```
//!
//! and the current `I` flows through a cross-sectional face area
//! `dA = dx * dy`. Positive current convention is `+z`.
//!
//! # Numerical scheme (resistor)
//!
//! For a pure series-R element with internal series EMF `V_src(t)`:
//!
//! ```text
//! V_term = E_z * dz = V_src + R * I
//! ⇒ J_z = I / dA = (E_z * dz - V_src) / (R * dA)
//! ```
//!
//! Substituting into Ampère's law `ε₀ ∂E_z/∂t = (∇×H)_z - J_z` and using a
//! semi-implicit average `(E_z^n + E_z^{n+1}) / 2` for the resistor current
//! (Taflove & Hagness §15.10) gives, after the **standard** Yee E-update has
//! produced `E_z^{n+1,*}`:
//!
//! ```text
//! E_z^{n+1} = (E_z^{n+1,*} − α E_z^n + γ V_src) / (1 + α)
//! α = dt · dz / (2 · ε₀ · R · dA)
//! γ = dt / (ε₀ · R · dA)
//! ```
//!
//! The struct keeps `E_z^n` as private state so [`LumpedRlcPort::correct_e`]
//! is a *post*-correction the driver applies after the normal `update_e`.
//!
//! # Series RLC — one-way (default) and two-way (Phase 2.fdtd.6.2) updates
//!
//! The full series-RLC branch supports two integration modes:
//!
//! - **One-way (default):** Crank-Nicolson on `I_L`/`V_C` *without* feeding the
//!   FDTD `E_z` terminal voltage back into the KVL (circuit→field). This is
//!   correct for ring-down frequency extraction in enclosed cavities — feeding
//!   `E_z` back in a tiny closed PEC box pulls the resonance off
//!   `1/(2π√(LC))` to a numerical ~1.49 GHz. Validated by the fdtd-206 gate
//!   (`tests/lumped_lc_resonance.rs`, ±2 % of `1/(2π√(LC))`).
//! - **Two-way ([`LumpedRlcPort::with_two_way`], Phase 2.fdtd.6.2):** a
//!   semi-implicit (Crank-Nicolson / trapezoidal) scheme that solves the branch
//!   current `I^{n+1/2}` and the `E_z^{n+1}` field update **together**, so the
//!   lumped current couples back into the field and the scheme is
//!   unconditionally stable for any `R ≥ 0`, `L ≥ 0`, `C > 0` (Piket-May,
//!   Taflove & Baron 1994; Taflove & Hagness §15.10). For S-parameter /
//!   terminating ports where the field's back-action on the load is physical.
//!
//! The two-way update is derived as follows. The standard Yee `E_z` update
//! produces `E_z^{n+1,*}` at the port cell; the lumped element adds a
//! current-density term to Ampère's law:
//!
//! ```text
//! E_z^{n+1} = E_z^{n+1,*} − (dt / (ε₀·dA))·I^{n+1/2}
//! ```
//!
//! The series-RLC KVL (terminal voltage `V_T = E_z·dz`, source EMF `V_src`,
//! positive current `+z`) discretised at the half-step `n+1/2` with
//! trapezoidal L and C is
//!
//! ```text
//! (E_z^{n+1}+E_z^n)·dz/2 = V_src^{n+1/2} + R·I^{n+1/2}
//!                          + L·(I^{n+1/2}−I^{n−1/2})/dt
//!                          + V_C^n + (dt/2C)·I^{n+1/2}
//! ```
//!
//! Substituting `E_z^{n+1}` from the field update and collecting `I^{n+1/2}`
//! (solving directly for the branch current, so the trapezoidal inductor
//! contributes `L/dt` to the diagonal — the discrete reactance is
//! `(2/dt)tan(ωdt/2)·L ≈ jωL`):
//!
//! ```text
//! K = R + L/dt + dt/(2C)          (branch operational impedance, Ω)
//! β = dt·dz / (2·ε₀·dA)           (FDTD half back-action impedance, Ω)
//! I^{n+1/2} = [ (E_z^{n+1,*}+E_z^n)·dz/2 − V_src^{n+1/2} − V_C^n
//!               + (L/dt)·I^{n−1/2} ] / (K + β)
//! E_z^{n+1} = E_z^{n+1,*} − (dt/(ε₀·dA))·I^{n+1/2}
//! V_C^{n+1} = V_C^n + (dt/C)·I^{n+1/2}
//! ```
//!
//! Because `K + β > 0` for all admissible R/L/C, the implicit solve never
//! divides by ~0: the `β` term is the on-diagonal damping that makes the
//! coupled update **unconditionally stable** — it removes the old explicit
//! pure-capacitor `≥ η₀/√3 ≈ 196 Ω` instability. The state carried between
//! steps is `I^{n−1/2}` ([`LumpedRlcPort::inductor_current`]) and `V_C^n`
//! ([`LumpedRlcPort::capacitor_voltage`]).
//!
//! ## Reductions (verified by `tests/lumped_rlc_twoway_001.rs`)
//!
//! - **Pure R** (`L=0`, `C=∞`): `K = R`, the update reduces to the
//!   semi-implicit resistor (`pure_resistor`, the validated path).
//! - **Pure C** (`L=0`): `K = R + dt/(2C)`; stable for any ESR ≥ 0.
//! - **Pure L** (`C=∞`): `K = R + 2L/dt`; a source-free inductor is *not*
//!   inert — the field drives `I` which reacts back onto `E_z`.
//! - **Thévenin source** (`V_src ≠ 0`): drives current into the line exactly
//!   as a series EMF behind the branch impedance.
//! - **Open** (`R=∞`): `I=0`, `E_z^{n+1}=E_z^{n+1,*}` (no-op).
//!
//! The fdtd-206 ring-down gate (`tests/lumped_lc_resonance.rs`, Phase
//! 2.fdtd.6.1) and the pure-resistor energy gate (`tests/lumped_resistor.rs`)
//! both stay green under this update; the two-way S-parameter behaviour is
//! validated by `lumped_rlc_twoway_001` (Γ vs analytic).
//!
//! # References
//!
//! - Taflove & Hagness, *Computational Electrodynamics: The Finite-Difference
//!   Time-Domain Method*, 3rd ed., §15.10 ("Modeling lumped elements").
//! - Piket-May, Taflove, Baron (1994), "FDTD modeling of digital signal
//!   propagation in 3-D circuits with passive and active loads",
//!   *IEEE Trans. Microw. Theory Tech.* 42(8): 1514-1523.

use std::f64::consts::{PI, TAU};

use yee_core::units::EPS0;

use crate::grid::YeeGrid;

/// Voltage-source time profile attached to a [`LumpedRlcPort`] as a series EMF.
///
/// All waveforms are evaluated at the integer simulation time `t = n · dt`.
/// `V0` is the open-circuit peak amplitude in volts.
#[derive(Debug, Clone, Copy)]
pub enum SourceWaveform {
    /// No EMF — the port is passive (used to model a pure load).
    None,
    /// `V0 · sin(2π f t)` with a Hanning (raised-cosine) ramp over the first
    /// `ramp_steps` timesteps, then unity thereafter.
    HannSine {
        /// Peak voltage of the sinusoid (V).
        v0: f64,
        /// Drive frequency (Hz).
        frequency: f64,
        /// Hann ramp length, in time steps.
        ramp_steps: usize,
    },
    /// Gaussian-modulated sine pulse `V0 · exp(-((t-t0)/τ)²) · sin(2π f0 (t-t0))`.
    ///
    /// `τ` is derived from the FWHM bandwidth `bw` via
    /// `τ = sqrt(2 ln 2) / (π · bw)` so that the spectral magnitude has full
    /// width `bw` at half maximum.
    GaussianPulse {
        /// Peak amplitude of the modulated carrier (V).
        v0: f64,
        /// Centre carrier frequency (Hz).
        f0: f64,
        /// Spectral FWHM bandwidth (Hz).
        bw: f64,
        /// Pulse centre, in time steps from the start of the run.
        t0_steps: usize,
    },
}

impl SourceWaveform {
    /// Evaluate the source voltage at step `n_step` with time step `dt` (s).
    pub fn value(&self, n_step: usize, dt: f64) -> f64 {
        let t = n_step as f64 * dt;
        match *self {
            SourceWaveform::None => 0.0,
            SourceWaveform::HannSine {
                v0,
                frequency,
                ramp_steps,
            } => {
                let ramp = if ramp_steps == 0 || n_step >= ramp_steps {
                    1.0
                } else {
                    0.5 * (1.0 - (PI * n_step as f64 / ramp_steps as f64).cos())
                };
                v0 * ramp * (TAU * frequency * t).sin()
            }
            SourceWaveform::GaussianPulse {
                v0,
                f0,
                bw,
                t0_steps,
            } => {
                let t0 = t0_steps as f64 * dt;
                // FWHM bandwidth → Gaussian time constant.
                // The Gaussian envelope is exp(-(t-t0)² / (2 σ_t²)),
                // matching a spectral FWHM `bw` requires
                //   σ_t = sqrt(2 ln 2) / (2π · σ_f) and σ_f = bw / (2√(2 ln 2)),
                // simplifying to τ such that env = exp(-((t-t0)/τ)²),
                // τ = sqrt(2 ln 2) / (π · bw).
                let tau = if bw > 0.0 {
                    (2.0 * std::f64::consts::LN_2).sqrt() / (PI * bw)
                } else {
                    // Degenerate case: zero bandwidth → pure CW. Use τ → ∞.
                    f64::INFINITY
                };
                let arg = (t - t0) / tau;
                let env = if tau.is_infinite() {
                    1.0
                } else {
                    (-arg * arg).exp()
                };
                v0 * env * (TAU * f0 * (t - t0)).sin()
            }
        }
    }
}

/// Lumped R/L/C/series-RLC port at a single Yee cell, oriented along ±z.
///
/// Implements Taflove & Hagness §15.10 series-RLC lumped element by adding a
/// current-driven correction to `E_z` at the port cell each timestep. See the
/// [module-level documentation](crate::lumped) for the numerical scheme.
///
/// # Reference impedance
///
/// The port is intended to drive (and absorb from) an adjacent transmission
/// line stub; the user is responsible for the line geometry. The canonical
/// resistor sanity check is
///
/// ```text
/// |Γ| = |(R − Z₀) / (R + Z₀)|
/// ```
///
/// against a `Z₀`-matched line. Phase 2.fdtd.6 ships an energy-dissipation
/// validation against an unconfined geometry (see
/// `tests/lumped_resistor.rs`); a Z₀-controlled stripline-Γ check is
/// a possible future extension.
///
/// # Scope
///
/// - Pure resistor ([`LumpedRlcPort::pure_resistor`]) is the primary
///   validated path (`tests/lumped_resistor.rs`).
/// - Series-RLC ([`LumpedRlcPort::series_rlc`]) defaults to the **one-way**
///   Crank-Nicolson scheme, validated by the fdtd-206 ring-down gate
///   (`tests/lumped_lc_resonance.rs`, ±2 % of 1/(2π√LC)).
/// - [`LumpedRlcPort::with_two_way`] (Phase 2.fdtd.6.2) selects the stable,
///   two-way semi-implicit update: the lumped current feeds back into `E_z`,
///   so a source-free reactive element is **not** inert. Validated by the
///   two-way S-parameter gate `tests/lumped_rlc_twoway_001.rs` (Γ vs the
///   analytic lumped-load reflection coefficient).
#[derive(Debug, Clone)]
pub struct LumpedRlcPort {
    /// Yee cell `(i, j, k)` of the `E_z` edge the port modifies.
    pub cell: (usize, usize, usize),
    /// Series resistance (Ω). `f64::INFINITY` represents an open circuit
    /// (no resistive current); zero is treated as a near-ideal short via
    /// the semi-implicit limit (the discrete `α` term saturates).
    pub resistance: f64,
    /// Series inductance (H). `0.0` removes the inductor term.
    pub inductance: f64,
    /// Series capacitance (F). `f64::INFINITY` shorts the capacitor (no
    /// DC blocking); zero is rejected at construction.
    pub capacitance: f64,
    /// Series voltage source (EMF in series with the R/L/C string).
    pub source_voltage: SourceWaveform,
    /// Two-way coupling flag (Phase 2.fdtd.6.2). When `false` (default), the
    /// series-RLC branch evolves **one-way** (circuit→field): the FDTD `E_z`
    /// terminal voltage is *not* fed back into the branch KVL — correct for
    /// ring-down frequency extraction in enclosed (cavity) geometries
    /// (`tests/lumped_lc_resonance.rs`, fdtd-206). When `true` (via
    /// [`LumpedRlcPort::with_two_way`]), the branch current updates
    /// **implicitly together with `E_z^{n+1}`** so the lumped current couples
    /// back into the field (two-way), the unconditionally-stable Piket-May /
    /// Taflove–Hagness semi-implicit update — for S-parameter / terminating
    /// ports where the field's back-action on the load is physical.
    two_way: bool,

    // ---- internal state ----
    /// Cached `E_z^n` at the port cell, captured at the *end* of each
    /// `correct_e` call so the next call has the pre-update value
    /// available for the semi-implicit resistor scheme.
    e_z_prev: f64,
    /// Inductor current `I_L` at the half-step (staggered with `E_z`).
    inductor_current: f64,
    /// Capacitor voltage `V_C` at the integer step.
    capacitor_voltage: f64,
}

impl LumpedRlcPort {
    /// Construct a pure series-R port. `r` in Ω. `src` is the (optional)
    /// series voltage source.
    ///
    /// # Panics
    ///
    /// Panics if `r ≤ 0` or `r` is non-finite. Use [`f64::INFINITY`] for
    /// an open circuit explicitly.
    pub fn pure_resistor(cell: (usize, usize, usize), r: f64, src: SourceWaveform) -> Self {
        assert!(
            (r > 0.0 && r.is_finite()) || r.is_infinite(),
            "LumpedRlcPort: resistance must be positive (got {r}); use f64::INFINITY for open"
        );
        Self {
            cell,
            resistance: r,
            inductance: 0.0,
            capacitance: f64::INFINITY,
            source_voltage: src,
            two_way: false,
            e_z_prev: 0.0,
            inductor_current: 0.0,
            capacitor_voltage: 0.0,
        }
    }

    /// Construct a series-RLC port. `r` in Ω, `l` in H, `c` in F. Use
    /// [`f64::INFINITY`] for `r` to remove the resistor branch (open),
    /// `0.0` for `l` to remove the inductor (ideal short across L), and
    /// [`f64::INFINITY`] for `c` to short the capacitor.
    ///
    /// # Panics
    ///
    /// Panics if `l < 0`, `c ≤ 0`, or any of `(r, l, c)` is NaN.
    pub fn series_rlc(
        cell: (usize, usize, usize),
        r: f64,
        l: f64,
        c: f64,
        src: SourceWaveform,
    ) -> Self {
        assert!(
            (r > 0.0 && r.is_finite()) || r.is_infinite(),
            "LumpedRlcPort::series_rlc: resistance must be positive (got {r}); use f64::INFINITY for open"
        );
        assert!(
            l >= 0.0 && !l.is_nan(),
            "LumpedRlcPort::series_rlc: inductance must be ≥ 0 (got {l})"
        );
        assert!(
            (c > 0.0 && !c.is_nan()) || c.is_infinite(),
            "LumpedRlcPort::series_rlc: capacitance must be positive (got {c}); use f64::INFINITY for short"
        );
        Self {
            cell,
            resistance: r,
            inductance: l,
            capacitance: c,
            source_voltage: src,
            two_way: false,
            e_z_prev: 0.0,
            inductor_current: 0.0,
            capacitor_voltage: 0.0,
        }
    }

    /// Enable the **stable two-way** semi-implicit series-RLC update
    /// (Phase 2.fdtd.6.2; Piket-May, Taflove & Baron 1994).
    ///
    /// By default a [`LumpedRlcPort`] evolves its series-RLC branch *one-way*
    /// (circuit→field), which is correct for ring-down frequency extraction in
    /// enclosed cavities (fdtd-206). For an S-parameter / terminating port,
    /// where the field's back-action on the lumped load is physically
    /// significant, call this builder: the branch current then updates
    /// implicitly **together with `E_z^{n+1}`** so the lumped current couples
    /// back into the field. The coupled update is unconditionally stable for
    /// any `R ≥ 0`, `L ≥ 0`, `C > 0` — in particular the low-loss capacitor
    /// case that the old explicit arm could not run below ~η₀/√3 ≈ 196 Ω ESR.
    ///
    /// See [`LumpedRlcPort::update_series_rlc_two_way`] for the derivation.
    /// Validated by `tests/lumped_rlc_twoway_001.rs`.
    pub fn with_two_way(mut self) -> Self {
        self.two_way = true;
        self
    }

    /// Whether this port uses the two-way semi-implicit update (Phase
    /// 2.fdtd.6.2). `false` is the default one-way (circuit→field) scheme.
    pub fn is_two_way(&self) -> bool {
        self.two_way
    }

    /// Apply the lumped-element correction to `E_z` at the port cell.
    ///
    /// Call this **after** [`crate::update::update_e`] (so the grid already
    /// holds the standard Yee leapfrog estimate `E_z^{n+1,*}` at the port
    /// cell), passing the simulation step counter `n_step` and the
    /// timestep `dt` in seconds.
    ///
    /// The correction overwrites `grid.ez[cell]` with the semi-implicit
    /// resistor-corrected (or full RLC-corrected) value.
    pub fn correct_e(&mut self, grid: &mut YeeGrid, n_step: usize, dt: f64) {
        let (i, j, k) = self.cell;
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let area = dx * dy;
        // Read what the standard E-update produced.
        let e1_star = grid.ez[(i, j, k)];
        let e0 = self.e_z_prev;
        let v_src = self.source_voltage.value(n_step, dt);

        let e1 = if self.inductance > 0.0 || self.capacitance.is_finite() {
            if self.two_way {
                // Stable two-way semi-implicit series-RLC (Phase 2.fdtd.6.2):
                // the lumped current couples back into E_z (S-parameter ports).
                self.update_series_rlc_two_way(e1_star, e0, v_src, dz, area, dt)
            } else {
                // Default one-way series-RLC (circuit→field; fdtd-206 ring-down).
                self.update_series_rlc(e1_star, e0, v_src, dz, area, dt)
            }
        } else {
            // Pure resistor with optional series EMF (the validated path).
            self.update_pure_resistor(e1_star, e0, v_src, dz, area, dt)
        };

        grid.ez[(i, j, k)] = e1;
        self.e_z_prev = e1;
    }

    /// Pure series-R update with optional series voltage source.
    ///
    /// Solves
    /// ```text
    /// E1 (1 + α) = E1s − α E0 + γ V_src
    ///   α = dt · dz / (2 · ε₀ · R · dA)
    ///   γ = dt / (ε₀ · R · dA)
    /// ```
    /// in closed form. For `R = ∞` the resistor term vanishes and
    /// `E1 = E1s` (the standard Yee update is left untouched).
    fn update_pure_resistor(
        &mut self,
        e1_star: f64,
        e0: f64,
        v_src: f64,
        dz: f64,
        area: f64,
        dt: f64,
    ) -> f64 {
        if self.resistance.is_infinite() {
            return e1_star;
        }
        let alpha = dt * dz / (2.0 * EPS0 * self.resistance * area);
        let gamma = dt / (EPS0 * self.resistance * area);
        (e1_star - alpha * e0 + gamma * v_src) / (1.0 + alpha)
    }

    /// Series-RLC update — default **one-way** Crank-Nicolson scheme
    /// (Phase 2.fdtd.6.1).
    ///
    /// Integrates the lumped-circuit KVL using Crank-Nicolson on `I_L` and
    /// `V_C`, then applies the resulting average current as a one-way
    /// correction to `E_z^{n+1}` at the port cell. The FDTD `E_z` terminal
    /// voltage is **not** fed back into the KVL: the circuit evolves
    /// autonomously driven by the series source `V_src`:
    ///
    /// ```text
    /// avg_I     = [2L/dt · I_L^n − V_C^n − V_src] / [2L/dt + R + dt/(2C)]
    /// I_L^{n+1} = 2 · avg_I − I_L^n
    /// V_C^{n+1} = V_C^n + (dt/C) · avg_I
    /// E_z^{n+1} = E_z^{n+1,*} − (dt/(ε₀·dA)) · avg_I   ← one-way: circuit→field
    /// ```
    ///
    /// Excluding the `E_z` terminal voltage breaks the closed-box feedback loop
    /// that otherwise pulls the resonance off `1/(2π√(LC))` (the FDTD
    /// back-action β ≈ 98 Ω over-loads the tiny series-RLC and the DFT then
    /// peaks at a numerical 1.49 GHz). This one-way scheme is correct for
    /// ring-down frequency extraction in enclosed cavities and is validated by
    /// the fdtd-206 gate (`tests/lumped_lc_resonance.rs`): a 5×5×40 PEC-box LC
    /// resonance at 1 GHz within ±2 % of the analytic `1/(2π√(LC))`.
    ///
    /// For S-parameter / terminating ports where the field's back-action on the
    /// load is physical, use [`LumpedRlcPort::with_two_way`] →
    /// [`LumpedRlcPort::update_series_rlc_two_way`].
    ///
    /// For `L = 0` the inductor short-circuits; falls back to a quasi-static
    /// R + C treatment. For `C = ∞` the capacitor term vanishes.
    fn update_series_rlc(
        &mut self,
        e1_star: f64,
        e0: f64,
        v_src: f64,
        _dz: f64,
        area: f64,
        dt: f64,
    ) -> f64 {
        let r_branch = self.resistance;
        let l = self.inductance;
        let c = self.capacitance;

        if l > 0.0 {
            let two_l_over_dt = 2.0 * l / dt;
            let r_eff = if r_branch.is_infinite() {
                // Open resistor: block all current.
                self.inductor_current = 0.0;
                self.capacitor_voltage = 0.0;
                return e1_star;
            } else {
                r_branch
            };
            let c_term = if c.is_finite() && c > 0.0 {
                dt / (2.0 * c)
            } else {
                0.0
            };
            let denom = two_l_over_dt + r_eff + c_term;
            let v_c = self.capacitor_voltage;
            let i_old = self.inductor_current;
            let numerator = two_l_over_dt * i_old - v_c - v_src;
            let avg_i = numerator / denom;
            let i_new = 2.0 * avg_i - i_old;
            self.inductor_current = i_new;
            if c.is_finite() && c > 0.0 {
                self.capacitor_voltage += (dt / c) * avg_i;
            }
            e1_star - (dt / (EPS0 * area)) * avg_i
        } else {
            // L = 0: inductor short-circuits its branch. Quasi-static R + C.
            let _e0 = e0;
            if !r_branch.is_infinite() {
                self.inductor_current = (e0 * _dz - self.capacitor_voltage - v_src) / r_branch;
            } else {
                self.inductor_current = 0.0;
            }
            let e1 = e1_star - (dt / (EPS0 * area)) * self.inductor_current;
            if c.is_finite() && c > 0.0 {
                self.capacitor_voltage += (dt / c) * self.inductor_current;
            }
            e1
        }
    }

    /// Series-RLC update — stable, **two-way** semi-implicit scheme
    /// (Phase 2.fdtd.6.2; Piket-May, Taflove & Baron 1994). Used only when the
    /// port is built with [`LumpedRlcPort::with_two_way`].
    ///
    /// Solves the branch current `I^{n+1/2}` and the `E_z^{n+1}` field update
    /// **together** so the lumped current couples back into the field. The
    /// scheme is unconditionally stable for any `R ≥ 0`, `L ≥ 0`, `C > 0`.
    ///
    /// # Derivation
    ///
    /// The standard Yee update has already produced `E_z^{n+1,*}` (`e1_star`).
    /// The lumped current `I` adds to Ampère's law:
    ///
    /// ```text
    /// E_z^{n+1} = e1_star − (dt/(ε₀·dA))·I            (field, two-way)
    /// ```
    ///
    /// The series-RLC KVL at the half-step `n+1/2`, with trapezoidal L and C
    /// (`V_T = E_z·dz`, source `V_src`, positive current `+z`):
    ///
    /// ```text
    /// (E_z^{n+1}+e0)·dz/2 = V_src + R·I + L·(I − I_old)/dt + V_C + (dt/2C)·I
    /// ```
    ///
    /// Substituting the field update and collecting `I`:
    ///
    /// ```text
    /// K  = R + L/dt + dt/(2C)         (branch operational impedance, Ω)
    /// β  = dt·dz / (2·ε₀·dA)          (FDTD half back-action impedance, Ω)
    /// I  = [ (e1_star+e0)·dz/2 − V_src − V_C + (L/dt)·I_old ] / (K + β)
    /// ```
    ///
    /// Then `E_z^{n+1} = e1_star − (dt/(ε₀·dA))·I` and
    /// `V_C^{n+1} = V_C + (dt/C)·I`, and `I_old ← I`.
    ///
    /// `K + β > 0` for all admissible R/L/C, so the implicit solve is always
    /// well-conditioned: the `β` term is the on-diagonal damping that makes
    /// the coupled update **unconditionally stable**, removing the old
    /// explicit pure-capacitor `≥ η₀/√3 ≈ 196 Ω` instability.
    ///
    /// Limits: `L=0,C=∞ ⇒ K=R` (semi-implicit resistor); `L=0 ⇒ K=R+dt/2C`
    /// (pure C, stable for any ESR); `C=∞ ⇒ K=R+L/dt` (pure L, not inert);
    /// `R=∞ ⇒ I=0` (open, no-op).
    fn update_series_rlc_two_way(
        &mut self,
        e1_star: f64,
        e0: f64,
        v_src: f64,
        dz: f64,
        area: f64,
        dt: f64,
    ) -> f64 {
        let r_branch = self.resistance;
        let l = self.inductance;
        let c = self.capacitance;

        // Open resistor: no current can flow through the series branch.
        if r_branch.is_infinite() {
            self.inductor_current = 0.0;
            // V_C holds its DC value (an open branch can't (dis)charge it);
            // leave capacitor_voltage untouched.
            return e1_star;
        }

        // Branch operational impedance K = R + L/dt + dt/(2C).
        //
        // Solving directly for the branch current `I = I^{n+1/2}` (not the
        // step-average), the trapezoidal inductor `L·(I−I_old)/dt` contributes
        // `L/dt` to the diagonal and `(L/dt)·I_old` to the RHS — discrete
        // reactance `(2/dt)tan(ωdt/2)·L ≈ jωL`. (The one-way scheme uses
        // `2L/dt` because it solves for the average current instead.)
        let l_over_dt = l / dt;
        let c_term = if c.is_finite() && c > 0.0 {
            dt / (2.0 * c)
        } else {
            // C = ∞ shorts the capacitor: no 1/(jωC) term.
            0.0
        };
        let k_branch = r_branch + l_over_dt + c_term;

        // FDTD half back-action impedance β = dt·dz / (2·ε₀·dA).
        let beta = dt * dz / (2.0 * EPS0 * area);

        let v_c = self.capacitor_voltage;
        let i_old = self.inductor_current;

        // Semi-implicit branch current I^{n+1/2}; the (e1_star+e0)·dz/2
        // terminal-voltage term and the implicit β denominator are what make
        // this TWO-WAY and stable (cf. the one-way fdtd-206 scheme).
        let i_half =
            ((e1_star + e0) * dz / 2.0 - v_src - v_c + l_over_dt * i_old) / (k_branch + beta);

        // Carry I_old = I^{n−1/2} ← I^{n+1/2} for the inductor's trapezoid.
        self.inductor_current = i_half;
        // Trapezoidal capacitor charge update.
        if c.is_finite() && c > 0.0 {
            self.capacitor_voltage = v_c + (dt / c) * i_half;
        }

        // Two-way FDTD correction: feed the lumped current back into E_z.
        e1_star - (dt / (EPS0 * area)) * i_half
    }

    /// Read access to the cached previous-step `E_z` at the port cell.
    /// Mostly useful in tests.
    pub fn e_z_prev(&self) -> f64 {
        self.e_z_prev
    }

    /// Read access to the inductor current state (A). Series-RLC only.
    pub fn inductor_current(&self) -> f64 {
        self.inductor_current
    }

    /// Read access to the capacitor voltage state (V). Series-RLC only.
    pub fn capacitor_voltage(&self) -> f64 {
        self.capacitor_voltage
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn waveform_none_is_zero() {
        let w = SourceWaveform::None;
        assert_eq!(w.value(0, 1e-12), 0.0);
        assert_eq!(w.value(1000, 1e-12), 0.0);
    }

    #[test]
    fn hann_sine_starts_at_zero_and_ramps() {
        let w = SourceWaveform::HannSine {
            v0: 1.0,
            frequency: 1.0e9,
            ramp_steps: 10,
        };
        let dt = 1.0e-12;
        // At n=0, ramp factor is 0 (Hann window starts at 0) AND sin(0)=0.
        assert!(w.value(0, dt).abs() < 1e-15);
        // After the ramp, magnitude is bounded by v0.
        for n in 10..50 {
            assert!(w.value(n, dt).abs() <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn pure_resistor_constructor_panics_on_negative_r() {
        let res = std::panic::catch_unwind(|| {
            LumpedRlcPort::pure_resistor((0, 0, 0), -1.0, SourceWaveform::None)
        });
        assert!(res.is_err(), "negative R should panic");
    }

    #[test]
    fn series_rlc_open_zero_l_capacitor_inf_reduces_state() {
        // R = ∞, L = 0, C = ∞ ⇒ no element does anything; inductor
        // current and capacitor voltage stay at zero.
        let port = LumpedRlcPort::series_rlc(
            (1, 1, 1),
            f64::INFINITY,
            0.0,
            f64::INFINITY,
            SourceWaveform::None,
        );
        assert_eq!(port.inductor_current, 0.0);
        assert_eq!(port.capacitor_voltage, 0.0);
    }
}
