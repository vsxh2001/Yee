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
//!
//! # Multi-cell APERTURE lumped port (Phase 2.fdtd.6.9, ADR-0125)
//!
//! The single-cell two-way port above references **one Yee cell** for both its
//! terminal voltage (`V = E_z·dz`, a single edge) and its field back-action
//! (`(dt/(ε₀·dA))·I` with the bare `dA = dx²`). ADR-0124's dx-sweep showed this
//! makes a sharp L‖C resonance impossible: as the grid refines the inductor's
//! two-way back-action collapses as **O(dx²)** while the capacitor freezes at a
//! fixed per-cell short — the realized reactance is dx-dependent and the tank
//! degenerates to a transparent line.
//!
//! The [`LumpedRlcPort::aperture`] constructor (Phase 2.fdtd.6.9) fixes this by
//! referencing the field coupling to the **modal port face** of physical area
//! `A = w·h` (trace width × substrate height), NOT one Yee cell:
//!
//! - **Modal branch voltage** — `V = ∫E_z·dz` over the full substrate height
//!   (all `n_sub` `E_z` edges in a column), averaged over the `w`-direction
//!   columns, instead of one `E_z·dz` edge:
//!   ```text
//!   V_T = (1 / N_col) · Σ_columns ( Σ_height E_z(i,j,k)·dz )
//!   ```
//! - **Aperture-area back-action** — one *aggregate* branch current `I` threads
//!   the whole aperture as a sheet displacement current `J_z = I / A`; every
//!   `E_z` cell in the aperture is corrected with the **physical** area `A`:
//!   ```text
//!   E_z^{n+1}(cell) = E_z^{n+1,*}(cell) − (dt/(ε₀·A)) · I        ∀ cell ∈ aperture
//!   ```
//!   Because `A = w·h` is a fixed physical area (not `dx²`) and `V_T` integrates
//!   the full height (not one edge), the realized `Z_L` is **dx-independent** —
//!   the `O(dx²)` inductor collapse and the dx-frozen capacitor short are both
//!   removed (the root fix, ADR-0125 item 1).
//! - **(y,z) aperture distribution** — the lumped value is the *aggregate*
//!   `R`/`L`/`C` of the element (one branch over the whole face), not the
//!   ad-hoc `C/N`, `N·L` per-cell tiling of ADR-0124. The aperture
//!   normalization (modal `V_T`, sheet `J_z = I/A`) holds the aggregate `Z_L`
//!   fixed independent of how many cells span the face.
//!
//! The two-way coupled solve is the same Piket-May / Taflove–Hagness
//! semi-implicit scheme as the single-cell path, but with `dz → h`
//! (full substrate height for the modal voltage) and `dA → A` (aperture area
//! for the back-action). `K + β > 0` is preserved (unconditional stability),
//! and the **pure-R limit reduces to the validated resistor exactly** per cell
//! — see [`LumpedRlcPort::correct_e_aperture`]. The single-edge
//! `series_rlc` / `pure_resistor` / `with_two_way` path is untouched.

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

/// Modal port-face aperture for the multi-cell aperture lumped port
/// (Phase 2.fdtd.6.9, ADR-0125).
///
/// Describes the `(y, z)` port face the lumped element bridges: the set of
/// `E_z` cells spanning the aperture, the physical aperture area `A = w·h`
/// (trace width × substrate height), and the substrate height `h` the modal
/// branch voltage `V = ∫E_z·dz` integrates over. The element is a **single
/// aggregate branch** (one `R`/`L`/`C`) across the whole face — see the
/// [module-level docs](crate::lumped#multi-cell-aperture-lumped-port-phase-2fdtd69-adr-0125).
#[derive(Debug, Clone)]
pub struct ApertureSpec {
    /// The `E_z` cells `(i, j, k)` spanning the `(y, z)` aperture face. All
    /// share the same x-index `i` (the port plane); the `(j, k)` indices tile
    /// the trace width (`y`) × substrate height (`z`) the mode occupies.
    pub cells: Vec<(usize, usize, usize)>,
    /// Number of `z`-columns (width-direction `y` positions). The modal voltage
    /// averages the per-column height integral over these columns, so an
    /// aperture wider in cells does not multiply the modal `V`.
    pub n_columns: usize,
    /// Physical aperture area `A = w·h` (m²), the modal cross-section the
    /// aggregate branch current threads. The field back-action references this
    /// **physical** area, NOT the per-cell `dx²` — the dx-stability fix.
    pub area: f64,
    /// Physical substrate height `h` (m) the modal branch voltage integrates
    /// over (`V = ∫E_z·dz` over the full height). Used for the KVL terminal
    /// voltage; the per-cell `E_z·dz` sum over a column equals `E_z·h` for a
    /// uniform field.
    pub height: f64,
}

/// Lumped R/L/C/series-RLC port at a single Yee cell, oriented along ±z.
///
/// Implements Taflove & Hagness §15.10 series-RLC lumped element by adding a
/// current-driven correction to `E_z` at the port cell each timestep. See the
/// [module-level documentation](crate::lumped) for the numerical scheme.
///
/// For the **multi-cell aperture** variant (Phase 2.fdtd.6.9, ADR-0125) whose
/// field coupling references the modal port face `A = w·h` rather than a single
/// Yee cell — removing the `O(dx²)` reactance collapse — see
/// [`LumpedRlcPort::aperture`] and [`LumpedRlcPort::correct_e_aperture`].
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

    /// Multi-cell aperture spec (Phase 2.fdtd.6.9, ADR-0125). When `Some`, the
    /// port couples to the modal port face (`V = ∫E_z·dz` over the substrate
    /// height, back-action referenced to the physical area `A = w·h`) instead
    /// of the single `cell`; apply via [`LumpedRlcPort::correct_e_aperture`].
    /// `None` (default) is the single-edge path. An aperture port is always
    /// two-way coupled.
    aperture: Option<ApertureSpec>,

    // ---- internal state ----
    /// Cached `E_z^n` at the port cell (single-edge path) or the modal terminal
    /// voltage `V_T^n` (aperture path), captured at the *end* of each correct
    /// call so the next call has the pre-update value for the semi-implicit
    /// scheme.
    e_z_prev: f64,
    /// Inductor current `I_L` at the half-step (staggered with `E_z`). For the
    /// aperture port this is the aggregate branch current `I^{n+1/2}`.
    inductor_current: f64,
    /// Capacitor voltage `V_C` at the integer step.
    capacitor_voltage: f64,
    /// Aperture path only: the step-centred modal terminal voltage
    /// `V_T = (V_T* + V_T^n)/2` used in the last `correct_e_aperture` call.
    /// Exposed via [`LumpedRlcPort::last_terminal_voltage`] so a bench can read
    /// the port's OWN realized `(V, I)` and form the realized branch impedance
    /// `Z = V/I` directly — without the fragile line de-embed.
    last_terminal_voltage: f64,
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
            aperture: None,
            e_z_prev: 0.0,
            inductor_current: 0.0,
            capacitor_voltage: 0.0,
            last_terminal_voltage: 0.0,
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
            aperture: None,
            e_z_prev: 0.0,
            inductor_current: 0.0,
            capacitor_voltage: 0.0,
            last_terminal_voltage: 0.0,
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

    /// Construct a **multi-cell aperture** series-RLC port (Phase 2.fdtd.6.9,
    /// ADR-0125): one aggregate `R`/`L`/`C` branch bridging the modal port face
    /// described by `aperture` (the `(y, z)` cells, physical area `A = w·h`,
    /// substrate height `h`).
    ///
    /// Unlike [`LumpedRlcPort::series_rlc`] (which references one Yee cell —
    /// `V = E_z·dz`, back-action `(dt/(ε₀·dx²))·I`), this port references the
    /// **modal port face**: the branch voltage is `V = ∫E_z·dz` over the full
    /// substrate height averaged across the width columns, and the field
    /// back-action injects a sheet current `J_z = I/A` referenced to the
    /// **physical** aperture area `A`, NOT the single-cell `dx²`. This removes
    /// the `O(dx²)` inductor collapse (and the dx-frozen capacitor short) that
    /// makes a single-cell port's realized `Z_L` dx-dependent (ADR-0124).
    ///
    /// The port is always two-way coupled. Apply it with
    /// [`LumpedRlcPort::correct_e_aperture`] after the standard `update_e`.
    /// `R`/`L`/`C` are the **aggregate** element values across the whole face
    /// (one branch), not per-cell — the aperture normalization holds the
    /// aggregate `Z_L` fixed regardless of cell count.
    ///
    /// # Panics
    ///
    /// Panics under the same R/L/C validity rules as
    /// [`LumpedRlcPort::series_rlc`], and additionally if the aperture has no
    /// cells, no columns, or a non-finite / non-positive `area` or `height`.
    pub fn aperture(aperture: ApertureSpec, r: f64, l: f64, c: f64, src: SourceWaveform) -> Self {
        assert!(
            (r > 0.0 && r.is_finite()) || r.is_infinite(),
            "LumpedRlcPort::aperture: resistance must be positive (got {r}); use f64::INFINITY for open"
        );
        assert!(
            l >= 0.0 && !l.is_nan(),
            "LumpedRlcPort::aperture: inductance must be ≥ 0 (got {l})"
        );
        assert!(
            (c > 0.0 && !c.is_nan()) || c.is_infinite(),
            "LumpedRlcPort::aperture: capacitance must be positive (got {c}); use f64::INFINITY for short"
        );
        assert!(
            !aperture.cells.is_empty(),
            "LumpedRlcPort::aperture: aperture must span at least one cell"
        );
        assert!(
            aperture.n_columns >= 1,
            "LumpedRlcPort::aperture: aperture must have ≥ 1 column"
        );
        assert!(
            aperture.area.is_finite() && aperture.area > 0.0,
            "LumpedRlcPort::aperture: area must be finite and positive (got {})",
            aperture.area
        );
        assert!(
            aperture.height.is_finite() && aperture.height > 0.0,
            "LumpedRlcPort::aperture: height must be finite and positive (got {})",
            aperture.height
        );
        // The representative `cell` is the first aperture cell (used only as a
        // fallback / for the single-edge `correct_e` if ever mis-called; the
        // aperture path ignores it and iterates `aperture.cells`).
        let cell = aperture.cells[0];
        Self {
            cell,
            resistance: r,
            inductance: l,
            capacitance: c,
            source_voltage: src,
            two_way: true,
            aperture: Some(aperture),
            e_z_prev: 0.0,
            inductor_current: 0.0,
            capacitor_voltage: 0.0,
            last_terminal_voltage: 0.0,
        }
    }

    /// Whether this port is a multi-cell aperture port (Phase 2.fdtd.6.9).
    pub fn is_aperture(&self) -> bool {
        self.aperture.is_some()
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

        let has_l = self.inductance > 0.0;
        let has_c = self.capacitance.is_finite();
        let e1 = if has_l || has_c {
            if self.two_way {
                // Two-way S-parameter ports (Phase 2.fdtd.6.4, ADR-0118): the
                // canonical Taflove–Hagness *per-element* lumped updates. These
                // couple to the field correctly per timestep so a shunt inductor
                // presents `jωL` and a shunt capacitor presents `1/(jωC)` — not
                // the instantaneous-`K` loading that ADR-0117 proved defective.
                match (has_l, has_c) {
                    // Pure capacitor (L = 0): effective-permittivity update.
                    (false, true) => self.update_pure_capacitor_two_way(e1_star, e0, dz, area),
                    // Pure inductor (C = ∞): accumulated branch current.
                    (true, false) => {
                        self.update_pure_inductor_two_way(e1_star, e0, v_src, dz, area, dt)
                    }
                    // Series R-L-C: combined accumulated-current + V_C update.
                    (true, true) => {
                        self.update_series_rlc_two_way(e1_star, e0, v_src, dz, area, dt)
                    }
                    // has_l || has_c guarantees at least one reactive term.
                    (false, false) => unreachable!(),
                }
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

    /// Apply the **multi-cell aperture** lumped correction (Phase 2.fdtd.6.9,
    /// ADR-0125). Call after the standard [`crate::update::update_e`].
    ///
    /// This is the dx-stable counterpart of [`LumpedRlcPort::correct_e`]: the
    /// field coupling references the modal port face (physical area `A = w·h`,
    /// substrate height `h`) instead of one Yee cell.
    ///
    /// # Update equations
    ///
    /// 1. **Modal terminal voltage** — average the per-column height path
    ///    integral `Σ_height E_z^{n+1,*}·dz` over the `N_col` width columns:
    ///    ```text
    ///    V_T* = (1/N_col) · Σ_columns Σ_height E_z^{n+1,*}(i,j,k)·dz
    ///    ```
    ///    For a uniform modal field `V_T* ≈ E_z·h`. We use the post-update
    ///    `E_z^{n+1,*}` (and the cached previous `V_T^n`) symmetrically with the
    ///    single-edge scheme.
    /// 2. **Aggregate branch current** — one branch current `I^{n+1/2}` for the
    ///    whole face, solved by the same semi-implicit Piket-May / Taflove
    ///    scheme as the single-edge two-way path, but with `dz → h`, `dA → A`:
    ///    ```text
    ///    K = R + L/dt + dt/(2C)            (branch impedance, Ω)
    ///    β = dt·h / (2·ε₀·A)              (aperture FDTD half back-action, Ω)
    ///    I = [ (V_T* + V_T^n)/2 − V_src − V_C^n + (L/dt)·I_old ] / (K + β)
    ///    ```
    ///    (Pure-R reduces to `K = R`, β with aperture `A`/`h` — the exact
    ///    aperture resistor; pure-C/pure-L use the same `A`/`h` references.)
    /// 3. **(y,z) aperture back-action** — the branch current threads the whole
    ///    face as a sheet `J_z = I/A`; every aperture cell is corrected with the
    ///    **physical** `A`:
    ///    ```text
    ///    E_z^{n+1}(cell) = E_z^{n+1,*}(cell) − (dt/(ε₀·A))·I    ∀ cell
    ///    ```
    ///
    /// Because `A`, `h`, and `R`/`L`/`C` are all physical (dx-independent), the
    /// realized `Z_L` no longer collapses as `O(dx²)` — the decisive fix
    /// (ADR-0125 item 1). `K + β > 0` keeps it unconditionally stable.
    ///
    /// # Panics
    ///
    /// Panics if the port was not built with [`LumpedRlcPort::aperture`].
    pub fn correct_e_aperture(&mut self, grid: &mut YeeGrid, n_step: usize, dt: f64) {
        let spec = self
            .aperture
            .as_ref()
            .expect("correct_e_aperture called on a non-aperture port")
            .clone();
        let dz = grid.dz;
        let h = spec.height;
        let area = spec.area;
        let n_col = spec.n_columns as f64;
        let v_src = self.source_voltage.value(n_step, dt);

        // (1) Modal terminal voltage from the post-update field: average the
        // per-column height path integral Σ E_z·dz over the width columns.
        let mut v_sum = 0.0;
        for &(i, j, k) in &spec.cells {
            v_sum += grid.ez[(i, j, k)] * dz;
        }
        let v_term_star = v_sum / n_col; // V_T* (modal, post-update)
        let v_term_prev = self.e_z_prev; // cached modal V_T^n (see end)
        // Step-centred modal terminal voltage logged for the realized-impedance
        // probe (a bench reads V/I directly from the port — no line de-embed).
        self.last_terminal_voltage = 0.5 * (v_term_star + v_term_prev);

        // (2) Aggregate branch current, semi-implicit two-way solve with the
        // aperture references (dz→h, dA→A). The pure-R limit reduces exactly to
        // the validated semi-implicit resistor (in the aperture sense).
        let i_branch = self.aperture_branch_current(v_term_star, v_term_prev, v_src, h, area, dt);

        // (3) Distribute the sheet current J_z = I/A back onto every aperture
        // cell with the PHYSICAL area A (not dx²) — the dx-stable back-action.
        let back = (dt / (EPS0 * area)) * i_branch;
        for &(i, j, k) in &spec.cells {
            grid.ez[(i, j, k)] -= back;
        }

        // Cache the modal terminal voltage for the next step's V_T^n. We read
        // it back AFTER the correction so the semi-implicit average uses the
        // realized (corrected) modal voltage, consistent with the single-edge
        // scheme caching the corrected `E_z`.
        let mut v_sum_post = 0.0;
        for &(i, j, k) in &spec.cells {
            v_sum_post += grid.ez[(i, j, k)] * dz;
        }
        self.e_z_prev = v_sum_post / n_col;
    }

    /// Solve the aggregate aperture branch current `I^{n+1/2}` for one step.
    ///
    /// Mirrors the single-edge two-way semi-implicit scheme (Piket-May,
    /// Taflove & Hagness §15.10) but with the modal terminal voltage `V_T`
    /// (= `∫E_z·dz` over the substrate height) and the aperture back-action
    /// impedance `β = dt·h/(2·ε₀·A)`. Carries the inductor current and
    /// capacitor voltage as state. `R = ∞` (open) blocks the branch.
    fn aperture_branch_current(
        &mut self,
        v_term_star: f64,
        v_term_prev: f64,
        v_src: f64,
        h: f64,
        area: f64,
        dt: f64,
    ) -> f64 {
        let r = self.resistance;
        let l = self.inductance;
        let c = self.capacitance;

        // Open resistor: no branch current.
        if r.is_infinite() {
            self.inductor_current = 0.0;
            return 0.0;
        }

        // Branch operational impedance K = R + L/dt + dt/(2C).
        let l_over_dt = if l > 0.0 { l / dt } else { 0.0 };
        let c_term = if c.is_finite() && c > 0.0 {
            dt / (2.0 * c)
        } else {
            0.0
        };
        let k_branch = r + l_over_dt + c_term;
        // Aperture FDTD half back-action impedance β = dt·h/(2·ε₀·A).
        let beta = dt * h / (2.0 * EPS0 * area);

        let v_c = self.capacitor_voltage;
        let i_old = self.inductor_current;
        // Step-centred modal terminal voltage (V_T* + V_T^n)/2 — the same
        // semi-implicit average the single-edge resistor uses, generalised to
        // the modal voltage.
        let v_term_mid = 0.5 * (v_term_star + v_term_prev);
        // I = [ V_T_mid − V_src − V_C + (L/dt)·I_old ] / (K + β).
        let i_branch = (v_term_mid - v_src - v_c + l_over_dt * i_old) / (k_branch + beta);

        self.inductor_current = i_branch;
        if c.is_finite() && c > 0.0 {
            self.capacitor_voltage = v_c + (dt / c) * i_branch;
        }
        i_branch
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

    /// Pure-capacitor two-way update — **canonical Taflove effective
    /// permittivity** (Phase 2.fdtd.6.4, ADR-0118). Used when `two_way` and
    /// `L = 0`, `C` finite.
    ///
    /// A lumped capacitor `C` bridging an `E_z` edge augments that cell's
    /// natural free-space capacitance. The textbook result (Taflove & Hagness
    /// §15.10; Piket-May, Taflove & Baron 1994) is a *local effective
    /// permittivity*
    ///
    /// ```text
    /// ε_eff = ε₀ + C·dz/dA
    /// ```
    ///
    /// so the cell's Ampère update runs at `ε_eff` instead of `ε₀`. The
    /// standard Yee step has already produced `E_z^{n+1,*}` at `ε₀`:
    ///
    /// ```text
    /// E_z^{n+1,*} = E_z^n + (dt/ε₀)·(∇×H)_z
    /// ⇒ (∇×H)_z = (ε₀/dt)·(E_z^{n+1,*} − E_z^n)
    /// ```
    ///
    /// Re-running that same curl term at `ε_eff` gives the canonical capacitor
    /// update purely in terms of the already-computed `E_z^{n+1,*}`:
    ///
    /// ```text
    /// E_z^{n+1} = E_z^n + (dt/ε_eff)·(∇×H)_z
    ///           = E_z^n + (ε₀/ε_eff)·(E_z^{n+1,*} − E_z^n)
    /// ```
    ///
    /// Because `ε_eff ≥ ε₀`, the update can only *raise* the cell capacitance,
    /// so it is unconditionally stable (no CFL penalty). The element presents
    /// `Z_C = 1/(jωC)` to the line: at high frequency the term `ε₀/ε_eff → 0`
    /// freezes the field (near-short), at low frequency `→ 1` (near-open),
    /// exactly the `1/(jωC)` reactance. An optional series ESR is ignored here
    /// (the gate's reactive case drives `R → 0`); a lossy capacitor uses the
    /// series-RLC arm.
    fn update_pure_capacitor_two_way(&mut self, e1_star: f64, e0: f64, dz: f64, area: f64) -> f64 {
        let c = self.capacitance;
        // ε_eff = ε₀ + C·dz/dA. C is finite & > 0 by construction here.
        let eps_eff = EPS0 + c * dz / area;
        let ratio = EPS0 / eps_eff; // ε₀/ε_eff ∈ (0, 1]
        // E_z^{n+1} = E_z^n + (ε₀/ε_eff)·(E_z^{n+1,*} − E_z^n).
        e0 + ratio * (e1_star - e0)
    }

    /// Pure-inductor two-way update — **canonical Taflove accumulated branch
    /// current** (Phase 2.fdtd.6.4, ADR-0118). Used when `two_way` and `L > 0`,
    /// `C = ∞`.
    ///
    /// A lumped inductor `L` bridging an `E_z` edge carries an auxiliary branch
    /// current that *accumulates* the terminal voltage `V = E_z·dz` (the
    /// textbook lumped-L FDTD source, Taflove & Hagness §15.10):
    ///
    /// ```text
    /// I_L^{n+1/2} = I_L^{n−1/2} + (dt·dz/L)·E_z^n           (Faraday: dI/dt = V/L)
    /// E_z^{n+1}   = E_z^{n+1,*} − (dt/(ε₀·dA))·I_L^{n+1/2}  (Ampère: −J back-action)
    /// ```
    ///
    /// The current is integrated **explicitly** from the *present* field `E_z^n`
    /// (`e0`), so — unlike the defective instantaneous-`K` scheme (ADR-0117) —
    /// the inductor is NOT loaded by the huge `L/dt` term in a single step; the
    /// accumulated `I_L` builds over many steps and presents the physical
    /// `Z_L = jωL` to the line. Stable: an inductor adds no CFL constraint.
    /// An optional series ESR is ignored here (the gate drives `R → 0`); a
    /// lossy inductor uses the series-RLC arm. `V_src` enters as a series EMF
    /// (`V = E_z·dz − V_src`).
    fn update_pure_inductor_two_way(
        &mut self,
        e1_star: f64,
        e0: f64,
        v_src: f64,
        dz: f64,
        area: f64,
        dt: f64,
    ) -> f64 {
        let l = self.inductance;
        // Accumulate the terminal voltage onto the branch current:
        //   I_L^{n+1/2} = I_L^{n−1/2} + (dt·dz/L)·E_z^n − (dt/L)·V_src.
        // (V across the inductor = E_z·dz − V_src for a series EMF.)
        let v_term = e0 * dz - v_src;
        self.inductor_current += (dt / l) * v_term;
        // Ampère back-action: feed the branch current into E_z.
        e1_star - (dt / (EPS0 * area)) * self.inductor_current
    }

    /// Series-RLC two-way update — **canonical Taflove combined lumped-RLC `E`
    /// update** (Phase 2.fdtd.6.4, ADR-0118). Used when `two_way` and both
    /// `L > 0` and `C` finite.
    ///
    /// One branch current `I` flows through `R`, `L`, `C` in series, driven by
    /// the terminal voltage `V_T = E_z·dz` (minus any series EMF `V_src`). The
    /// canonical discretisation accumulates the inductor current explicitly
    /// (as in the pure-L arm) while carrying the capacitor voltage `V_C` as a
    /// state and treating `R` semi-implicitly for stability. KVL at step `n`:
    ///
    /// ```text
    /// L·dI/dt = V_T − R·I − V_C − V_src,     dV_C/dt = I/C
    /// ```
    ///
    /// Discretising `L·dI/dt` with the leapfrog increment, `R·I` with the
    /// step-centred average `(I^{n+1/2}+I^{n−1/2})/2`, and `V_C`, `V_T` at the
    /// integer step `n`:
    ///
    /// ```text
    /// (L/dt)(I^{n+1/2} − I^{n−1/2}) = E_z^n·dz − V_src − V_C^n
    ///                                 − (R/2)(I^{n+1/2}+I^{n−1/2})
    /// ⇒ I^{n+1/2} = [ (L/dt − R/2)·I^{n−1/2}
    ///                 + E_z^n·dz − V_src − V_C^n ] / (L/dt + R/2)
    /// ```
    ///
    /// then the Ampère back-action `E_z^{n+1} = E_z^{n+1,*} −
    /// (dt/(ε₀·dA))·I^{n+1/2}` and the trapezoidal capacitor charge
    /// `V_C^{n+1} = V_C^n + (dt/C)·I^{n+1/2}`. The branch presents
    /// `R + jωL + 1/(jωC)` to the line. Reductions: `R→0` recovers the LC
    /// resonator; `R=∞` blocks the branch (open, no-op). `R ≥ 0` keeps the
    /// `L/dt + R/2` denominator strictly positive, so the solve is always
    /// well-conditioned; the explicit-from-`E_z^n` coupling presents the
    /// physical reactance rather than the instantaneous-`K` loading (ADR-0117).
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
            // V_C holds its DC value (an open branch can't (dis)charge it).
            return e1_star;
        }

        let l_over_dt = l / dt;
        // Semi-implicit resistor: R/2 on the diagonal, (L/dt − R/2) on I_old.
        let r_half = 0.5 * r_branch;
        let v_c = self.capacitor_voltage;
        let i_old = self.inductor_current;

        // Branch current accumulated from the present terminal voltage E_z^n·dz
        // (explicit field coupling → physical reactance, not instantaneous-K).
        let v_term = e0 * dz - v_src - v_c;
        let i_half = ((l_over_dt - r_half) * i_old + v_term) / (l_over_dt + r_half);

        self.inductor_current = i_half;
        // Trapezoidal capacitor charge update.
        if c.is_finite() && c > 0.0 {
            self.capacitor_voltage = v_c + (dt / c) * i_half;
        }

        // Ampère back-action: feed the branch current into E_z.
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

    /// Aperture port only: the step-centred modal terminal voltage `V_T`
    /// (volts) from the most recent [`LumpedRlcPort::correct_e_aperture`] call.
    ///
    /// Paired with [`LumpedRlcPort::inductor_current`] (the aggregate branch
    /// current `I`), this lets a bench form the port's OWN realized branch
    /// impedance `Z = V_T(ω)/I(ω)` directly — the physical `R + jωL + 1/(jωC)`
    /// the discrete port presents, independent of the surrounding line `Z₀` and
    /// the fragile shunt de-embed. Used by `tests/aperture_port_001.rs` for the
    /// dx-stability check (ADR-0125).
    pub fn last_terminal_voltage(&self) -> f64 {
        self.last_terminal_voltage
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

    // ---- Aperture port (Phase 2.fdtd.6.9, ADR-0125) ----

    fn one_cell_spec(cell: (usize, usize, usize), dx: f64) -> ApertureSpec {
        // A degenerate "aperture" of one cell with A = dx² and h = dz = dx — so
        // its references collapse onto the single-edge resistor's dA = dx²,
        // dz = dx. Used to prove the exact resistor reduction.
        ApertureSpec {
            cells: vec![cell],
            n_columns: 1,
            area: dx * dx,
            height: dx,
        }
    }

    #[test]
    fn aperture_constructor_validates() {
        let spec = one_cell_spec((1, 1, 1), 1e-3);
        // Empty aperture must panic.
        let mut bad = spec.clone();
        bad.cells.clear();
        assert!(
            std::panic::catch_unwind(|| {
                LumpedRlcPort::aperture(bad, 50.0, 0.0, f64::INFINITY, SourceWaveform::None)
            })
            .is_err(),
            "empty aperture should panic"
        );
        // Non-positive area must panic.
        let mut bad = one_cell_spec((1, 1, 1), 1e-3);
        bad.area = 0.0;
        assert!(
            std::panic::catch_unwind(|| {
                LumpedRlcPort::aperture(bad, 50.0, 0.0, f64::INFINITY, SourceWaveform::None)
            })
            .is_err(),
            "zero area should panic"
        );
        // A valid aperture port reports itself as aperture + two-way.
        let p = LumpedRlcPort::aperture(spec, 50.0, 0.0, f64::INFINITY, SourceWaveform::None);
        assert!(p.is_aperture());
        assert!(p.is_two_way());
    }

    #[test]
    fn aperture_pure_resistor_reduces_to_single_edge_exactly() {
        // A one-cell aperture (A = dx², h = dz = dx) pure-R port must produce
        // the SAME E_z update as the validated single-edge semi-implicit
        // resistor, step for step — the resistor-exact reduction (ADR-0125
        // item 4). Drive both with the same E1*/E0 history.
        let dx = 1.0e-3;
        let dt = 1.5e-12;
        let r = 73.0;
        let cell = (2, 2, 2);

        let grid_dims = (5, 5, 5);
        let mut grid_ap = YeeGrid::vacuum(grid_dims.0, grid_dims.1, grid_dims.2, dx);
        let mut grid_se = YeeGrid::vacuum(grid_dims.0, grid_dims.1, grid_dims.2, dx);

        let mut ap = LumpedRlcPort::aperture(
            one_cell_spec(cell, dx),
            r,
            0.0,
            f64::INFINITY,
            SourceWaveform::None,
        );
        // Single-edge two-way resistor (K = R, the validated semi-implicit
        // limit). Use the public single-edge path.
        let mut se = LumpedRlcPort::series_rlc(cell, r, 0.0, f64::INFINITY, SourceWaveform::None)
            .with_two_way();

        // Feed an identical post-update E1* sequence into both and compare the
        // corrected E_z.
        let e1_star_seq = [1.0, 0.7, -0.4, 0.2, -0.9, 0.55, -0.1];
        for (n, &e1s) in e1_star_seq.iter().enumerate() {
            grid_ap.ez[cell] = e1s;
            grid_se.ez[cell] = e1s;
            ap.correct_e_aperture(&mut grid_ap, n, dt);
            se.correct_e(&mut grid_se, n, dt);
            let a = grid_ap.ez[cell];
            let s = grid_se.ez[cell];
            assert!(
                (a - s).abs() <= 1e-12 * (a.abs().max(s.abs()).max(1.0)),
                "step {n}: aperture pure-R {a} != single-edge two-way resistor {s}"
            );
        }
    }

    #[test]
    fn aperture_open_resistor_is_noop() {
        let dx = 1.0e-3;
        let dt = 1.5e-12;
        let cell = (2, 2, 2);
        let mut grid = YeeGrid::vacuum(5, 5, 5, dx);
        let mut ap = LumpedRlcPort::aperture(
            one_cell_spec(cell, dx),
            f64::INFINITY, // open
            1.0e-9,        // L present but blocked by open R
            f64::INFINITY,
            SourceWaveform::None,
        );
        grid.ez[cell] = 0.42;
        ap.correct_e_aperture(&mut grid, 0, dt);
        assert_eq!(grid.ez[cell], 0.42, "open aperture R must be a no-op");
        assert_eq!(ap.inductor_current(), 0.0);
    }

    #[test]
    fn aperture_reactive_stays_finite_over_many_steps() {
        // Stability: a multi-cell aperture pure-inductor and pure-capacitor must
        // not blow up (K + β > 0 keeps the implicit solve well-conditioned).
        let dx = 1.0e-3;
        let dt = 1.5e-12;
        let port_i = 2;
        let cells: Vec<(usize, usize, usize)> = (1..4)
            .flat_map(|j| (0..3).map(move |k| (port_i, j, k)))
            .collect();
        let n_cols = 3;
        let spec = ApertureSpec {
            cells: cells.clone(),
            n_columns: n_cols,
            area: (3.0 * dx) * (3.0 * dx),
            height: 3.0 * dx,
        };
        for (l, c) in [(2.0e-9, f64::INFINITY), (0.0, 5.0e-13)] {
            let mut grid = YeeGrid::vacuum(5, 5, 5, dx);
            let mut ap = LumpedRlcPort::aperture(spec.clone(), 1e-6, l, c, SourceWaveform::None);
            for n in 0..2000 {
                // Drive a bounded oscillating field at the aperture cells.
                let drive = (n as f64 * 0.3).sin();
                for &cc in &cells {
                    grid.ez[cc] = drive;
                }
                ap.correct_e_aperture(&mut grid, n, dt);
                for &cc in &cells {
                    assert!(
                        grid.ez[cc].is_finite(),
                        "aperture reactive update went non-finite at step {n} (L={l}, C={c})"
                    );
                }
            }
        }
    }
}
