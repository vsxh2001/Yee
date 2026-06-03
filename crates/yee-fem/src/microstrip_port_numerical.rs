//! Numerical-eigenmode microstrip wave-port for the open-boundary FEM
//! driven solver (FEM-EM brick N1, ADR-0154).
//!
//! This module is the *high-fidelity* successor to the analytic-shape
//! [`crate::microstrip_port`] (B3): instead of an idealised flat-`E_z`
//! transverse modal shape, it feeds the **true** quasi-TEM transverse
//! field — solved by yee-mom's SHIPPED 2-D cross-section eigensolver
//! ([`yee_mom::ports::NumericalCrossSection`]) — as the FEM port's
//! `modal_e_t`.
//!
//! ## Why a numerical eigenmode (the v1 ~9 % overlap floor)
//!
//! The analytic v1 shape ([`crate::modal_e_t_microstrip`]) is `E_z`-
//! dominant in the trace↔ground gap with a single-exponential air tail,
//! uniform in `x`. It is the *qualitatively* correct quasi-TEM
//! polarisation, but only *partially* overlaps the true eigenmode (the
//! real `E_z` is concentrated under the trace and fringes outward, with
//! a small in-plane `E_x` component the v1 drops). When the driven solve
//! projects the FEM field onto this coarse shape, the modal overlap is
//! low: the B4 straight-line gate (`fem_line_eeff_001`) measured
//! `|S21| ≈ 0.089` (≈ −21 dB/port) — the *phase* is coherent (so the
//! two-length ε_eff is recovered to 0.61 %), but the *amplitude* sits on
//! a modal-overlap floor that caps a clean filter `|S21|`.
//!
//! Feeding the numerical eigenmode lifts that overlap: the de-risk probe
//! that seeded this module (ADR-0154) measured `|S21| ≈ 0.778`,
//! `|S11| ≈ 0.087`, ε_eff to 0.61 % — an ~8.7× amplitude lift over v1
//! with the phase fidelity unchanged.
//!
//! ## Geometry / frame contract
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
//! A wave-port face sits on a `y = const` end-cap, so its transverse
//! plane is the `(x, z)` cross-section and propagation is `ŷ`. The
//! yee-mom `NumericalCrossSection` is a 2-D mesh over `(coord0, coord1) =
//! (x_width, substrate-normal)`, and `e_tangential_at(x, y)` returns the
//! in-cross-section transverse field `[E_x, E_normal]`. The frame map is
//! therefore:
//!
//! ```text
//!   coord0 (x_width)         ↔  FEM x̂
//!   coord1 (substrate-normal) ↔ FEM ẑ
//!   sample e_tangential_at(p.x, p.z) → [ex, e_normal]
//!   emit Vector3::new(ex, 0.0, e_normal)   (E_y = 0: transverse shape)
//! ```
//!
//! This is the validated frame map of the ADR-0154 de-risk probe — the
//! substrate-normal eigenmode component is placed onto the FEM `ẑ`, and
//! the propagation-axis (`ŷ`) component is `0` for a transverse modal
//! shape. On a `(x, z)` port face that substrate-normal field is
//! *in-plane*, so it is representable by the 3-D Whitney-1 edge basis the
//! FEM scatter projects onto — the representability that makes the
//! microstrip wave-port well-posed in 3-D FEM where it is **ill-posed for
//! planar MoM** (CLAUDE.md §10, ADR-0153/0154).
//!
//! ## Cost
//!
//! [`microstrip_port_numerical`] runs the yee-mom 2-D eigensolve **once**
//! per call (sub-second on the internal validated-density mesh) and wraps
//! the resulting mode in an [`Arc`] so the returned port's `modal_e_t`
//! closure takes shared ownership of the heap-allocated eigenmode without
//! a `Copy` bound (the β closure captures only scalars). The β is the
//! analytic Hammerstad-Jensen
//! [`crate::beta_microstrip`] (the eigensolve's own β is *not* used for
//! the absorbing-boundary impedance — only the modal *shape* comes from
//! the eigensolve), so the only changed variable vs the analytic v1 port
//! is the transverse field shape.

use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::Error;

use crate::microstrip_port::beta_microstrip;
use crate::open_boundary::PortDefinition;

/// Internal cross-section subdivisions for the numerical eigensolve.
///
/// Based on the validated `20 × 10` reference density of
/// `yee-mom/tests/eigensolver_microstrip_quasi_tem.rs`, with `ny` bumped
/// to `12` (production density `20 × 12`) so a 1 mm substrate inside a
/// 6 mm box is two cells tall
/// (`dy = 0.5 mm`) and the strip hole is a clean one-cell band — the same
/// `dx ≈ dz ≈ 0.5 mm` aspect the FEM volume cross-section carries
/// (`NX = NZ = 12` over 6 mm). This is the density the ADR-0154 probe
/// validated; it is **not** the coarse `8 × 8` that failed to converge in
/// `mom_002_numerical_waveport.rs`.
const XSEC_NX: usize = 20;
const XSEC_NY: usize = 12;

/// Geometry of a shielded straight-microstrip line, in metres.
///
/// The fields are the physical inputs the numerical-eigenmode port shares
/// with the [`crate::microstrip_mesh::layered_microstrip_mesh`] box the
/// port attaches to: the same `box_w` / `box_h` cross-section, substrate
/// height, trace width, and substrate permittivity. The port's transverse
/// modal field is the quasi-TEM eigenmode of this exact `(box_w × box_h)`
/// FR-4 cross-section, so the eigenmode is the true transverse mode of the
/// cross-section the FEM volume carries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MicrostripPortGeom {
    /// Trace (signal-strip) width along `x`, metres.
    pub trace_w: f64,
    /// Substrate height along `z`, metres: dielectric fills `z ∈ [0, sub_h]`.
    pub sub_h: f64,
    /// Substrate relative permittivity (e.g. 4.4 for FR-4).
    pub eps_r: f64,
    /// Box / substrate width along `x`, metres.
    pub box_w: f64,
    /// Box height along `z`, metres: substrate + air.
    pub box_h: f64,
}

/// Build the `(x, substrate-normal)` shielded-microstrip cross-section
/// [`yee_mesh::TriMesh2D`] for the numerical-eigenmode port, mirroring the
/// strip-as-hole builder validated in
/// `yee-mom/tests/eigensolver_microstrip_quasi_tem.rs` but at this line's
/// physical box: `[0, box_w] × [0, box_h]`, FR-4 (tag `1`) for
/// `coord1 < sub_h`, air (tag `0`) above, with the signal strip a
/// rectangular hole centred in `x` at `coord1 ∈ [sub_h, sub_h + t_strip]`.
///
/// Hole cells are omitted, so their border edges are mesh-boundary = PEC
/// (the signal-strip conductor inside the outer PEC box). The outer box
/// boundary is PEC = the shield + ground plane. A microstrip is a
/// two-conductor line, and a quasi-TEM mode exists only because there are
/// two separated conductors — the strip-as-hole construction supplies the
/// second (signal) conductor.
fn microstrip_cross_section(
    box_w: f64,
    box_h: f64,
    sub_h: f64,
    trace_w: f64,
    t_strip: f64,
    nx: usize,
    ny: usize,
) -> Result<yee_mesh::TriMesh2D, Error> {
    let xs: Vec<f64> = (0..=nx).map(|i| box_w * (i as f64) / (nx as f64)).collect();
    let ys: Vec<f64> = (0..=ny).map(|j| box_h * (j as f64) / (ny as f64)).collect();
    let xc = box_w / 2.0;
    let (sx0, sx1) = (xc - trace_w / 2.0, xc + trace_w / 2.0);
    let (sy0, sy1) = (sub_h, sub_h + t_strip);
    let in_strip = |cx: f64, cy: f64| {
        cx > sx0 - 1e-12 && cx < sx1 + 1e-12 && cy > sy0 - 1e-12 && cy < sy1 + 1e-12
    };

    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for &y in &ys {
        for &x in &xs {
            vertices.push([x, y]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;
    let mut triangles = Vec::new();
    let mut tags = Vec::new();
    for j in 0..ny {
        let yc = 0.5 * (ys[j] + ys[j + 1]);
        for i in 0..nx {
            let xcell = 0.5 * (xs[i] + xs[i + 1]);
            if in_strip(xcell, yc) {
                continue; // hole = signal-strip PEC conductor
            }
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v11 = idx(i + 1, j + 1);
            let v01 = idx(i, j + 1);
            let tag = if yc < sub_h { 1u32 } else { 0u32 };
            triangles.push([v00, v10, v11]);
            tags.push(tag);
            triangles.push([v00, v11, v01]);
            tags.push(tag);
        }
    }
    // `TriMesh2D::new` validates the cross-section (≥1 triangle, in-bounds
    // indices, CCW winding); a degenerate geometry (e.g. `trace_w ≥ box_w`
    // yielding zero cells) surfaces as `yee_mesh::Error`. Map it to the
    // crate-wide `yee_core::Error` and propagate instead of panicking on a
    // library `Result` path.
    yee_mesh::TriMesh2D::new(vertices, triangles, None, Some(tags))
        .map_err(|e| Error::Invalid(format!("microstrip cross-section: {e}")))
}

/// Construct a **numerical-eigenmode** quasi-TEM microstrip wave-port
/// [`PortDefinition`].
///
/// Solves the yee-mom 2-D quasi-TEM cross-section eigenmode of `geom`'s
/// `(box_w × box_h)` FR-4 cross-section **once** at `f_hz`, then returns a
/// single-mode port whose:
///
/// * **β(ω)** is the analytic Hammerstad-Jensen
///   [`crate::beta_microstrip`]`(trace_w, sub_h, eps_r, ω)` — the
///   absorbing-boundary impedance reference, *not* the eigensolve's β; and
/// * **`modal_e_t`** samples the numerical eigenmode's transverse field
///   via the validated frame map (see the module docs):
///   `e_tangential_at(p.x, p.z) → [ex, e_normal]` →
///   `Vector3::new(ex, 0.0, e_normal)`.
///
/// So the only quantity that differs from the analytic [`microstrip_port`]
/// is the transverse modal *shape*, which is what lifts the driven `|S21|`
/// out of the v1 ~0.089 (≈ −21 dB) modal-overlap floor (ADR-0154 probe:
/// `|S21| ≈ 0.778`).
///
/// [`microstrip_port`]: crate::microstrip_port
///
/// The cross-section is built at the internal validated default density
/// (`20 × 12`, one-cell-thick signal strip) — the same density the
/// ADR-0154 probe validated. The eigensolve is sub-second; the returned
/// port shares the solved mode through an [`Arc`].
///
/// # Errors
///
/// Returns a [`yee_core::Error`] (the crate-wide error type) if either the
/// cross-section mesh build fails (a degenerate geometry — e.g.
/// `trace_w ≥ box_w` — yielding an invalid `TriMesh2D`) or the yee-mom
/// quasi-TEM eigensolve fails to surface a propagating mode. The eigensolve
/// already returns [`yee_core::Error`] (propagated unchanged); the mesh
/// builder's `yee_mesh::Error` is mapped to [`yee_core::Error::Invalid`].
pub fn microstrip_port_numerical(
    geom: &MicrostripPortGeom,
    f_hz: f64,
) -> Result<PortDefinition, Error> {
    // One-cell-thick signal strip (matches the probe's
    // `t_strip = box_h / ny`), so the strip hole is a clean one-cell band.
    let t_strip = geom.box_h / (XSEC_NY as f64);
    let mesh = microstrip_cross_section(
        geom.box_w,
        geom.box_h,
        geom.sub_h,
        geom.trace_w,
        t_strip,
        XSEC_NX,
        XSEC_NY,
    )?;

    // Material tags: air tag 0 = (1, 0); FR-4 tag 1 = (eps_r, 0); μ_r = 1.
    let mut eps: HashMap<u32, Complex64> = HashMap::new();
    eps.insert(0u32, Complex64::new(1.0, 0.0)); // air
    eps.insert(1u32, Complex64::new(geom.eps_r, 0.0)); // FR-4
    let mut mu: HashMap<u32, Complex64> = HashMap::new();
    mu.insert(0u32, Complex64::new(1.0, 0.0));
    mu.insert(1u32, Complex64::new(1.0, 0.0));

    let mut mode = yee_mom::ports::NumericalCrossSection::new(mesh, eps, mu).with_quasi_tem();
    // The eigensolve error type is already `yee_core::Error`, so `?`
    // surfaces it unchanged — no error-variant remapping needed.
    mode.solve(f_hz)?;

    // The modal-shape closure takes shared ownership of the solved mode via
    // an `Arc` (cheap to hold, no `Copy` bound); the β closure below captures
    // only scalars and does not touch it.
    let mode = Arc::new(mode);

    let (trace_w, sub_h, eps_r) = (geom.trace_w, geom.sub_h, geom.eps_r);
    let beta = move |omega: f64| beta_microstrip(trace_w, sub_h, eps_r, omega);
    let e_t = move |p: Vector3<f64>| {
        // yee-mom coord0 = x_width, coord1 = substrate-normal. Sample at
        // (x = p.x, substrate-normal = p.z); returns [E_x, E_normal]. Map
        // the substrate-normal component onto FEM ẑ; the propagation-axis
        // (ŷ) component is 0 for a transverse modal shape.
        let et = mode.e_tangential_at(p.x, p.z);
        Vector3::new(et[0], 0.0, et[1])
    };

    Ok(PortDefinition::single_mode(Box::new(beta), Box::new(e_t)))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Line geometry of the `fem_line_eeff_numerical_001` gate: a ~50 Ω FR-4
    // line on 1 mm substrate in a 6 mm × 6 mm box.
    const TRACE_W: f64 = 1.0e-3;
    const SUB_H: f64 = 1.0e-3;
    const EPS_R: f64 = 4.4;
    const BOX_W: f64 = 6.0e-3;
    const BOX_H: f64 = 6.0e-3;

    fn line_geom() -> MicrostripPortGeom {
        MicrostripPortGeom {
            trace_w: TRACE_W,
            sub_h: SUB_H,
            eps_r: EPS_R,
            box_w: BOX_W,
            box_h: BOX_H,
        }
    }

    /// The numerical port's `modal_e_t` must be a finite, non-zero,
    /// `E_z`-dominant field in the trace↔ground gap that decays as the
    /// sample point rises into the air above the trace — the
    /// numerical-eigenmode analogue of
    /// `microstrip_port_closures::modal_e_t_is_ez_dominant_in_gap_and_decays_in_air`.
    ///
    /// FAST: builds the port (a sub-second 2-D eigensolve) and only
    /// *samples* the closure — no FEM volume solve.
    #[test]
    fn numerical_modal_e_t_is_ez_dominant_in_gap_and_decays_in_air() {
        let port = microstrip_port_numerical(&line_geom(), 10.0e9)
            .expect("numerical quasi-TEM cross-section eigensolve must succeed");
        let e_t = &port.modes[0].modal_e_t;

        // A point in the middle of the trace↔ground gap, under the centred
        // trace (x = box_w / 2, z = sub_h / 2).
        let x_trace = BOX_W / 2.0;
        let p_gap = Vector3::new(x_trace, 0.0, SUB_H / 2.0);
        let e_gap = e_t(p_gap);

        assert!(
            e_gap.x.is_finite() && e_gap.y.is_finite() && e_gap.z.is_finite(),
            "in-gap modal field must be finite, got {e_gap:?}"
        );
        let mag_gap = e_gap.norm();
        assert!(
            mag_gap > 0.0,
            "in-gap modal field must be non-zero, got |E| = {mag_gap}"
        );
        // The transverse mode is constructed as Vector3::new(E_x, 0, E_z),
        // so E_y is identically 0; the dominance test is E_z vs the
        // in-plane E_x.
        assert!(
            e_gap.z.abs() > e_gap.x.abs(),
            "in-gap field should be substrate-normal (E_z)-dominant: \
             |E_z| = {} not > |E_x| = {}",
            e_gap.z.abs(),
            e_gap.x.abs()
        );

        // A point well above the trace in air (z = sub_h + 3·sub_h): the
        // fringing field must have decayed below the in-gap magnitude.
        let p_air = Vector3::new(x_trace, 0.0, SUB_H + 3.0 * SUB_H);
        let e_air = e_t(p_air);
        assert!(
            e_air.norm() < mag_gap,
            "field above the trace in air (|E| = {}) should decay below the \
             in-gap value (|E| = {mag_gap})",
            e_air.norm()
        );
    }

    /// The numerical port's β(ω) must equal the analytic Hammerstad-Jensen
    /// phase constant `(ω/c)·sqrt(ε_eff)` derived from `yee_layout::eps_eff`
    /// — the same cross-check the v1
    /// `microstrip_port_closures::beta_matches_hammerstad_jensen` makes
    /// (to `1e-9`). The numerical port deliberately keeps the analytic β
    /// (the absorbing-impedance reference); only the modal *shape* is
    /// numerical.
    #[test]
    fn numerical_port_beta_matches_hammerstad_jensen() {
        use yee_core::units::C0;

        let omega = 2.0 * std::f64::consts::PI * 10.0e9;
        let port = microstrip_port_numerical(&line_geom(), 10.0e9)
            .expect("numerical quasi-TEM cross-section eigensolve must succeed");

        let beta_closure = (port.modes[0].beta_mode)(omega);
        let eps_eff = yee_layout::eps_eff(TRACE_W, SUB_H, EPS_R);
        let beta_expected = (omega / C0) * eps_eff.sqrt();

        assert!(
            (beta_closure - beta_expected).abs() < 1e-9,
            "numerical port β = {beta_closure}, Hammerstad-Jensen expected = \
             {beta_expected}, |diff| = {:e} exceeds 1e-9",
            (beta_closure - beta_expected).abs()
        );
        // Quasi-TEM β exceeds the free-space wavenumber (ε_eff > 1).
        assert!(
            beta_expected > omega / C0,
            "quasi-TEM β = {beta_expected} should exceed k0 = {} (ε_eff > 1)",
            omega / C0
        );
    }
}
