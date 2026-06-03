//! Quasi-TEM microstrip wave-port closures for the open-boundary FEM
//! driven solver (FEM-EM brick B3, ADR-0153).
//!
//! This module is the microstrip analogue of the WR-90 TE_{10}
//! wave-port fixture in `tests/open_boundary_sweep_matrix.rs`: it
//! supplies the `(β(ω), e_t(x))` closure pair that
//! [`crate::open_boundary::PortDefinition::single_mode`] consumes to
//! drive / absorb a microstrip line at a port face.
//!
//! ## Geometry contract
//!
//! The port attaches to a [`crate::microstrip_mesh::layered_microstrip_mesh`]
//! box, whose axes are fixed by that mesh builder:
//!
//! ```text
//!   x ∈ [0, box_w]    substrate / box width
//!   y ∈ [0, line_len] PROPAGATION axis
//!   z ∈ [0, box_h]    substrate-normal: ground plane at z = 0,
//!                      dielectric for z ∈ [0, sub_h], trace top at
//!                      z = sub_h, air above for z ∈ (sub_h, box_h].
//! ```
//!
//! A wave-port face sits on a `y = const` end-cap of the line, so its
//! **transverse plane is the `(x, z)` cross-section** and the
//! propagation direction is `ŷ`. The microstrip quasi-TEM mode's
//! dominant field is the substrate-normal `E_z` in the trace↔ground
//! gap — and crucially `E_z` is *in-plane* on a `(x, z)` port face, so
//! it is representable by the 3-D Whitney-1 edge basis (including the
//! vertical, `z`-running edges) that the FEM scatter projects onto
//! (`scatter_port_face`, modal dot-product). This is exactly the
//! representability that makes the microstrip wave-port well-posed in
//! 3-D FEM where it is **ill-posed for planar MoM** (CLAUDE.md §10:
//! the planar in-plane RWG basis cannot carry the substrate-normal
//! `E_z`); see ADR-0153.
//!
//! ## What this module is (B3) and is not (B4)
//!
//! B3 ships the closure *machinery* + a closure-sanity gate
//! (`tests/microstrip_port_closures.rs`): the analytic
//! Hammerstad-Jensen phase constant β(ω) and an `E_z`-dominant
//! transverse modal shape. The **fidelity** of the analytic modal
//! shape — i.e. whether driving a line with it recovers the correct
//! Hammerstad-Jensen ε_eff end-to-end — is **B4's** job and is *not*
//! asserted here. See [`microstrip_port`] for an honest read on the v1
//! shape's ε_eff prospects.

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_core::units::C0;

use crate::open_boundary::PortDefinition;

/// Quasi-TEM phase constant `β(ω)` of a microstrip line (rad/m).
///
/// ```text
///   β(ω) = (ω / c) · sqrt(ε_eff(w, h, εr))
/// ```
///
/// where `ε_eff` is the Hammerstad-Jensen / Schneider effective
/// permittivity from [`yee_layout::eps_eff`] (validated by the
/// `geo_002_hammerstad` gate in `crates/yee-layout`) and `c` is
/// [`yee_core::units::C0`].
///
/// Unlike the WR-90 TE_{10} guide, the quasi-TEM microstrip mode has
/// **no low-frequency cutoff** — `β` is positive for every `ω > 0`, so
/// (unlike `beta_te10`) there is no below-cutoff clip to `0`.
///
/// # Arguments
///
/// - `w` — trace width, metres.
/// - `h` — substrate height, metres.
/// - `eps_r` — substrate relative permittivity.
/// - `omega` — angular frequency `ω = 2πf`, rad/s.
pub fn beta_microstrip(w: f64, h: f64, eps_r: f64, omega: f64) -> f64 {
    let eps_eff = yee_layout::eps_eff(w, h, eps_r);
    (omega / C0) * eps_eff.sqrt()
}

/// Substrate-normal (`E_z`-dominant) quasi-TEM transverse modal shape,
/// **un-scaled** (the solver L²-normalises — see the
/// [`crate::open_boundary::PortMode::modal_e_t`] doc).
///
/// On a `(x, z)` port face (propagation along `ŷ`) the dominant
/// quasi-TEM field is the substrate-normal `E_z` filling the
/// parallel-plate-like trace↔ground gap. This v1 returns a field whose
/// only non-zero component is `E_z`, with the magnitude profile
///
/// ```text
///   |E_z|(z) = 1                       for 0 ≤ z ≤ sub_h     (gap: ~uniform)
///   |E_z|(z) = exp(−(z − sub_h)/d)     for z >  sub_h        (air: decaying)
/// ```
///
/// with the air-decay length `d = sub_h` (one substrate height — the
/// quasi-TEM fringing field above the trace decays on the order of the
/// substrate thickness). The field is **uniform in `x`** in this v1
/// (no trace-window confinement); see [`modal_e_t_microstrip_windowed`]
/// for an `x`-confined variant and [`microstrip_port`] for the honest
/// B4 fidelity caveat.
///
/// Returns an un-scaled `Vector3` with `E_x = E_y = 0` and `E_z`
/// per the profile above. Below the ground plane (`z < 0`, which does
/// not occur on a valid port face) the field is `0`.
///
/// # Arguments
///
/// - `sub_h` — substrate height (`= h`), metres; sets both the gap top
///   and the air-decay length.
/// - `p` — world-space sample point on the port face.
pub fn modal_e_t_microstrip(sub_h: f64, p: Vector3<f64>) -> Vector3<f64> {
    let ez = ez_profile(sub_h, p.z);
    Vector3::new(0.0, 0.0, ez)
}

/// `x`-confined variant of [`modal_e_t_microstrip`]: the same
/// substrate-normal `E_z` gap/air `z`-profile, multiplied by a smooth
/// raised-cosine taper in `x` centred on the trace and falling to zero
/// outside a `±w` half-window about the trace centre.
///
/// This is a closer analytic stand-in for the true quasi-TEM mode
/// (whose `E_z` is concentrated under the trace and fringes outward),
/// for use by the **B4** end-to-end ε_eff driver where the in-`x`
/// confinement of the source matters. The B3 gate exercises the
/// un-windowed [`modal_e_t_microstrip`]; this variant is provided so a
/// B4 driver can opt into trace-centred confinement without a separate
/// formulation.
///
/// The taper is
///
/// ```text
///   t(x) = ½(1 + cos(π · (x − x_c)/w))   for |x − x_c| ≤ w
///   t(x) = 0                              otherwise
/// ```
///
/// so `t(x_c) = 1` at the trace centre and `t` reaches `0` with zero
/// slope at `x_c ± w`.
///
/// # Arguments
///
/// - `x_c` — trace-centre `x` coordinate (`= box_w / 2` for a centred
///   trace), metres.
/// - `w` — trace width, metres; the taper half-window is `±w`.
/// - `sub_h` — substrate height, metres (gap top + air-decay length).
/// - `p` — world-space sample point on the port face.
pub fn modal_e_t_microstrip_windowed(
    x_c: f64,
    w: f64,
    sub_h: f64,
    p: Vector3<f64>,
) -> Vector3<f64> {
    let dx = (p.x - x_c).abs();
    let taper = if dx <= w {
        0.5 * (1.0 + (PI * (p.x - x_c) / w).cos())
    } else {
        0.0
    };
    let ez = ez_profile(sub_h, p.z) * taper;
    Vector3::new(0.0, 0.0, ez)
}

/// Shared `E_z` magnitude profile in `z` (substrate-normal): unit in
/// the trace↔ground gap `0 ≤ z ≤ sub_h`, exponentially decaying with
/// length `sub_h` in the air above, zero below the ground plane.
#[inline]
fn ez_profile(sub_h: f64, z: f64) -> f64 {
    if z < 0.0 {
        0.0
    } else if z <= sub_h {
        1.0
    } else {
        (-(z - sub_h) / sub_h).exp()
    }
}

/// Construct a quasi-TEM microstrip wave-port [`PortDefinition`].
///
/// Mirrors the WR-90 TE_{10} fixture
/// (`PortDefinition::single_mode(Box::new(beta_te10),
/// Box::new(modal_e_t_te10))`): bundles the analytic
/// [`beta_microstrip`] phase constant and the `E_z`-dominant
/// [`modal_e_t_microstrip`] transverse shape into a single-mode port
/// (incident amplitude `a_inc = 1`).
///
/// # Modal shape (v1) and B4 fidelity prospects
///
/// The v1 modal shape is **substrate-normal `E_z`, uniform in `x`,
/// unit in the trace↔ground gap and exponentially decaying into the
/// air above the trace** (decay length = one substrate height). It is
/// `E_z`-dominant in the gap, which is the qualitatively correct
/// quasi-TEM polarisation, and on a `(x, z)` port face that `E_z` is
/// in-plane and so representable by the 3-D Whitney-1 basis (ADR-0153).
///
/// **Honest read for B4:** this analytic shape is a *coarse* stand-in.
/// Two known idealisations will bound its end-to-end ε_eff accuracy
/// when the B4 driven solve projects the FEM field onto it:
///
/// 1. **No `x`-confinement.** The true quasi-TEM `E_z` is concentrated
///    under the trace and fringes outward over ~one substrate height;
///    this v1 is flat in `x`. The trace-centred
///    [`modal_e_t_microstrip_windowed`] is provided for B4 to opt into
///    `x`-confinement if the flat profile proves inadequate.
/// 2. **Idealised gap/air `z`-profile.** A uniform gap field plus a
///    single-exponential air tail is not the exact conformal-mapping
///    quasi-TEM field; the fringing distribution and the small
///    in-plane `E_x` (transverse-to-trace) component are dropped.
///
/// Either could shift the projected ε_eff away from the
/// Hammerstad-Jensen target. Whether that shift is within B4's
/// tolerance is **B4's gate to decide** — B3 asserts only closure
/// sanity (β value, `E_z`-dominance, finite/non-zero modal
/// self-inner-product). If B4 finds the flat-`x` shape inadequate, the
/// windowed variant (or a refined conformal-mapping field) is the next
/// lever, *not* a change to this β.
///
/// # Arguments
///
/// - `w` — trace width, metres.
/// - `h` — substrate height, metres.
/// - `eps_r` — substrate relative permittivity.
pub fn microstrip_port(w: f64, h: f64, eps_r: f64) -> PortDefinition {
    let beta = move |omega: f64| beta_microstrip(w, h, eps_r, omega);
    let e_t = move |p: Vector3<f64>| modal_e_t_microstrip(h, p);
    PortDefinition::single_mode(Box::new(beta), Box::new(e_t))
}

/// Construct a quasi-TEM microstrip wave-port [`PortDefinition`] with a
/// trace-centred `x`-confinement (raised-cosine taper) on the modal
/// shape.
///
/// Identical β(ω) to [`microstrip_port`], but the modal field uses
/// [`modal_e_t_microstrip_windowed`] so the source `E_z` is confined to
/// a `±w` half-window about the trace centre `x_c`. Intended for the
/// **B4** driven solve where the in-`x` confinement of the modal
/// source may matter for ε_eff fidelity (see [`microstrip_port`]).
///
/// # Arguments
///
/// - `box_w` — substrate / box width, metres; the trace centre is taken
///   as `box_w / 2`.
/// - `w` — trace width, metres.
/// - `h` — substrate height, metres.
/// - `eps_r` — substrate relative permittivity.
pub fn microstrip_port_windowed(box_w: f64, w: f64, h: f64, eps_r: f64) -> PortDefinition {
    let x_c = box_w / 2.0;
    let beta = move |omega: f64| beta_microstrip(w, h, eps_r, omega);
    let e_t = move |p: Vector3<f64>| modal_e_t_microstrip_windowed(x_c, w, h, p);
    PortDefinition::single_mode(Box::new(beta), Box::new(e_t))
}
