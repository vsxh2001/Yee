//! Automatic meshing + convergence loop (FS.0a, ADR-0204).
//!
//! The market research behind `FULL-SUITE-ROADMAP.md` found manual mesh
//! selection to be the #1 practitioner-cited barrier to open-EM-tool
//! adoption: a novice cannot know the λ/20, substrate-resolution, and
//! feature-resolution rules, and results are sensitive to getting them
//! wrong. This module is that rulebook as code, plus the HFSS-style
//! adaptive-pass loop in FDTD flavour: solve, refine dx, re-solve, stop
//! when the S-curve stops moving. No kernel change — the loop rides the
//! shared [`crate::board`] fixture, so every design flow (gates, studio,
//! Python, WS) gets push-button meshing the same way.

use yee_layout::Layout;

use crate::board::{TwoPortBoardOptions, two_port_board_job};
use crate::{JobEvent, sparams};

/// Speed of light in vacuum, m/s.
const C0: f64 = 299_792_458.0;

/// The meshing rulebook: the largest dx that satisfies every rule.
///
/// - **Wavelength**: `dx ≤ λ_min/20` with `λ_min = c/(f_max·√ε_r)` — the
///   shortest in-dielectric wavelength the drive contains.
/// - **Substrate**: `dx ≤ h/3` — at least three cells across the substrate
///   so the vertical quasi-TEM field is resolved (the S.9 CPML collapse
///   was ultimately a substrate-resolution interaction).
/// - **Feature**: `dx ≤ w_min/2` — at least two cells across the smallest
///   trace width or gap in the layout (the R.4 coupling-floor lesson:
///   under-resolved gaps read wrong couplings, silently).
///
/// The result is clamped to `[1 µm, 1 mm]` — below 1 µm the volumetric
/// FDTD premise itself breaks down for board work (the MMIC caveat in the
/// roadmap), above 1 mm nothing at RF board scale is resolved.
pub fn auto_dx(layout: &Layout, f_max_hz: f64) -> f64 {
    let lambda_min = C0 / (f_max_hz * layout.substrate.eps_r.sqrt());
    let by_wavelength = lambda_min / 20.0;
    let by_substrate = layout.substrate.height_m / 3.0;
    let by_feature = min_feature_m(layout) / 2.0;
    by_wavelength
        .min(by_substrate)
        .min(by_feature)
        .clamp(1e-6, 1e-3)
}

/// The smallest feature the mesh must resolve: the minimum over every
/// trace rectangle's width/height and every inter-trace gap along x/y
/// (axis-aligned bounding-box gap between polygon pairs; the generators
/// in this workspace emit axis-aligned rectangles, so this is exact for
/// them and conservative-ish for arbitrary polygons).
pub fn min_feature_m(layout: &Layout) -> f64 {
    let mut min_f = f64::INFINITY;
    let boxes: Vec<(f64, f64, f64, f64)> = layout
        .traces
        .iter()
        .map(|p| {
            let (mut x0, mut y0, mut x1, mut y1) = (
                f64::INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            );
            for v in &p.verts {
                x0 = x0.min(v.x);
                y0 = y0.min(v.y);
                x1 = x1.max(v.x);
                y1 = y1.max(v.y);
            }
            (x0, y0, x1, y1)
        })
        .collect();
    for &(x0, y0, x1, y1) in &boxes {
        min_f = min_f.min(x1 - x0).min(y1 - y0);
    }
    for (a, &(ax0, ay0, ax1, ay1)) in boxes.iter().enumerate() {
        for &(bx0, by0, bx1, by1) in boxes.iter().skip(a + 1) {
            // Gap along x when the boxes overlap in y, and vice versa.
            let x_gap = (bx0 - ax1).max(ax0 - bx1);
            let y_gap = (by0 - ay1).max(ay0 - by1);
            let y_overlap = ay1.min(by1) - ay0.max(by0);
            let x_overlap = ax1.min(bx1) - ax0.max(bx0);
            if x_gap > 0.0 && y_overlap > 0.0 {
                min_f = min_f.min(x_gap);
            }
            if y_gap > 0.0 && x_overlap > 0.0 {
                min_f = min_f.min(y_gap);
            }
        }
    }
    min_f
}

/// One convergence pass: the dx it ran at and its |S21| curve (dB).
#[derive(Debug, Clone)]
pub struct ConvergencePass {
    /// Cell size of this pass, metres.
    pub dx_m: f64,
    /// Directional |S21| in dB at each requested frequency.
    pub s21_db: Vec<f64>,
}

/// The convergence-loop result.
#[derive(Debug, Clone)]
pub struct Converged {
    /// Every pass, coarsest first; the last is the answer.
    pub passes: Vec<ConvergencePass>,
    /// Max |Δ|S21|| in **linear magnitude** between the final two passes.
    /// Linear, not dB, deliberately: near a deep notch a tiny frequency or
    /// depth shift produces tens of dB of per-bin delta while the linear
    /// change is milliunits — the first gate run measured exactly that
    /// (Δ = 15 dB at a converged 4.900 GHz notch). Commercial adaptive
    /// refinement (HFSS's ΔS) uses the linear metric for the same reason.
    pub final_delta: f64,
    /// Whether `final_delta ≤ tol` within the pass budget. `false` is
    /// reported, never hidden: the caller decides whether an unconverged
    /// answer is usable.
    pub converged: bool,
}

/// Run one two-port measurement (reference + DUT) at the given options;
/// returns the directional |S21| curve.
fn measure(
    layout: &Layout,
    reference: &Layout,
    opts: &TwoPortBoardOptions,
    freqs_hz: &[f64],
) -> Result<Vec<f64>, String> {
    let run = |l: &Layout| -> Result<(Vec<Vec<f64>>, f64, f64), String> {
        let job = two_port_board_job(l, opts)?;
        let (dt, spacing) = (job.dt_s, job.spacing_m);
        let handle = crate::submit(job.spec);
        for event in handle.events() {
            match event {
                JobEvent::Done { result } => return Ok((result.probes, dt, spacing)),
                JobEvent::Error { message } => return Err(message),
                _ => {}
            }
        }
        Err("engine stream ended without a result".into())
    };
    let (ref_p, dt, spacing) = run(reference)?;
    let (dut_p, dt2, _) = run(layout)?;
    if dt != dt2 {
        return Err("passes diverged in dt".into());
    }
    Ok(sparams::directional_transmission_db(
        [&dut_p[3], &dut_p[4], &dut_p[5]],
        [&ref_p[3], &ref_p[4], &ref_p[5]],
        dt,
        spacing,
        freqs_hz,
    ))
}

/// The adaptive-pass loop (FDTD flavour of HFSS's adaptive refinement):
/// starting from `opts.dx_m` (use [`auto_dx`] to seed it), solve the
/// two-port, refine `dx → dx/√2`, and stop when the max per-frequency
/// Δ|S21| stops moving. **Every pass must solve the same physical problem**,
/// so everything the fixture sizes in cells is rescaled to hold its
/// physical size: `n_steps` (constant time window), the CPML margin, the
/// air height under the lid, and the CPML absorber depth. The first loop
/// version scaled none of the last three — at dx₀/2 the lid sat at half
/// height and the absorber was half thickness, and the DUT (which scatters
/// into those boundaries where the reference line doesn't) read a
/// non-physical broadband |S21| up to +10.7 dB. Convergence is judged on
/// the max per-frequency Δ|S21| in
/// **linear magnitude** between consecutive passes is ≤ `tol` (HFSS's
/// ΔS ≈ 0.02 is the commercial reference point; staircased FDTD needs a
/// looser walking-skeleton value — the graded grid of FS.0b tightens it)
/// — or the pass budget runs out, which the result reports honestly.
///
/// Cost note: each pass is 2 solves and the finest pass dominates
/// (cells ×2^1.5 per pass, plus steps ×√2). This is exactly the workload
/// the GPU backend exists for — set `opts.backend` accordingly.
pub fn converge_two_port(
    layout: &Layout,
    reference: &Layout,
    mut opts: TwoPortBoardOptions,
    freqs_hz: &[f64],
    tol: f64,
    max_passes: usize,
) -> Result<Converged, String> {
    assert!(max_passes >= 2, "convergence needs at least two passes");
    assert!(tol > 0.0 && tol.is_finite(), "tol must be positive");
    let lin = |db: f64| 10.0_f64.powf(db / 20.0);
    let base_steps = opts.n_steps as f64 * opts.dx_m;
    // Physical fixture sizes at the starting dx: the loop must vary ONLY the
    // discretization, so the CPML margin and the air height are held constant
    // in metres (their cell counts grow as dx shrinks), not in cells.
    let margin_m = opts.margin_cells as f64 * opts.dx_m;
    let air_above_m = opts.air_above_cells as f64 * opts.dx_m;
    let npml_m = opts.npml as f64 * opts.dx_m;
    let mut passes: Vec<ConvergencePass> = Vec::new();
    let mut final_delta = f64::INFINITY;
    for _ in 0..max_passes {
        // Keep the physical time window constant as dx (and thus dt) shrink.
        opts.n_steps = (base_steps / opts.dx_m).round() as usize;
        opts.margin_cells = (margin_m / opts.dx_m).round() as usize;
        opts.air_above_cells = (air_above_m / opts.dx_m).round() as usize;
        opts.npml = (npml_m / opts.dx_m).round() as usize;
        let s21_db = measure(layout, reference, &opts, freqs_hz)?;
        if let Some(prev) = passes.last() {
            final_delta = s21_db
                .iter()
                .zip(&prev.s21_db)
                .map(|(a, b)| (lin(*a) - lin(*b)).abs())
                .fold(0.0_f64, f64::max);
        }
        passes.push(ConvergencePass {
            dx_m: opts.dx_m,
            s21_db,
        });
        if final_delta <= tol {
            return Ok(Converged {
                passes,
                final_delta,
                converged: true,
            });
        }
        opts.dx_m /= std::f64::consts::SQRT_2;
    }
    Ok(Converged {
        passes,
        final_delta,
        converged: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use yee_layout::{BBox, Point2, Polygon, PortRef, Substrate};

    fn layout_with(traces: Vec<Polygon>, eps_r: f64, h_m: f64) -> Layout {
        let bbox = BBox::from_polygons(&traces);
        Layout {
            substrate: Substrate {
                eps_r,
                height_m: h_m,
                loss_tangent: 0.0,
                metal_thickness_m: 35e-6,
            },
            traces,
            ports: vec![PortRef {
                at: Point2::new(0.0, 0.0),
                width_m: 1e-3,
                ref_impedance_ohm: 50.0,
            }],
            bbox,
        }
    }

    #[test]
    fn each_rule_binds_when_it_is_the_constraint() {
        // Wide line, thick substrate, high f: wavelength rule binds.
        // λ_min at 10 GHz in ε_r 4.4 = 14.3 mm → /20 = 0.715 mm... above
        // the 1 mm... use 30 GHz: λ_min = 4.77 mm → /20 = 0.238 mm.
        let l = layout_with(vec![Polygon::rect(0.0, 0.0, 50e-3, 10e-3)], 4.4, 3e-3);
        let dx = auto_dx(&l, 30.0e9);
        let lam = 299_792_458.0 / (30.0e9 * 4.4_f64.sqrt());
        assert!((dx - lam / 20.0).abs() < 1e-12, "wavelength rule: {dx}");

        // Thin substrate binds: h = 0.3 mm → h/3 = 0.1 mm.
        let l = layout_with(vec![Polygon::rect(0.0, 0.0, 50e-3, 10e-3)], 4.4, 0.3e-3);
        let dx = auto_dx(&l, 5.0e9);
        assert!((dx - 0.1e-3).abs() < 1e-12, "substrate rule: {dx}");

        // Narrow gap binds: two 10 mm-wide lines 0.15 mm apart → 75 µm.
        let l = layout_with(
            vec![
                Polygon::rect(0.0, 0.0, 50e-3, 10e-3),
                Polygon::rect(0.0, 10.15e-3, 50e-3, 10e-3),
            ],
            4.4,
            1.6e-3,
        );
        let dx = auto_dx(&l, 5.0e9);
        assert!((dx - 0.075e-3).abs() < 1e-12, "feature rule: {dx}");
    }

    #[test]
    fn min_feature_finds_widths_and_gaps() {
        let l = layout_with(
            vec![
                Polygon::rect(0.0, 0.0, 20e-3, 1.5e-3),
                Polygon::rect(0.0, 2.1e-3, 20e-3, 1.5e-3), // y-gap 0.6 mm
            ],
            4.4,
            1.6e-3,
        );
        assert!((min_feature_m(&l) - 0.6e-3).abs() < 1e-12);
        // Single wide trace: its own height is the feature.
        let l = layout_with(vec![Polygon::rect(0.0, 0.0, 20e-3, 3e-3)], 4.4, 1.6e-3);
        assert!((min_feature_m(&l) - 3e-3).abs() < 1e-12);
    }

    #[test]
    fn auto_dx_is_clamped() {
        // Absurdly fine demand clamps at 1 µm.
        let l = layout_with(vec![Polygon::rect(0.0, 0.0, 1e-3, 1e-6)], 4.4, 1.6e-3);
        assert_eq!(auto_dx(&l, 5.0e9), 1e-6);
    }
}
