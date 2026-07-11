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
    let boxes = trace_boxes(layout);
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

/// The graded **coarse ceiling** (FS.0b.1, ADR-0210): [`auto_dx`] without
/// the feature rule — `min(λ_min/20, h/3)`, clamped to `[1 µm, 1 mm]`.
///
/// On a graded grid the feature rule moves into the local fine bands
/// ([`auto_spacings`]); keeping it out of the bulk ceiling is the payoff —
/// a single narrow gap no longer drags the whole domain to `feature/2`
/// the way the uniform rulebook must.
pub fn auto_dx_bulk(layout: &Layout, f_max_hz: f64) -> f64 {
    let lambda_min = C0 / (f_max_hz * layout.substrate.eps_r.sqrt());
    (lambda_min / 20.0)
        .min(layout.substrate.height_m / 3.0)
        .clamp(1e-6, 1e-3)
}

/// Axis-aligned bounding box `(x0, y0, x1, y1)` of every trace polygon —
/// exact for the axis-aligned rectangles the workspace generators emit.
fn trace_boxes(layout: &Layout) -> Vec<(f64, f64, f64, f64)> {
    layout
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
        .collect()
}

// ===========================================================================
// Graded mesh rules (FS.0b.1, ADR-0210)
// ===========================================================================

/// Options for [`auto_spacings`]: the fixture geometry (in **metres** —
/// the ADR-0204 loop-hygiene lesson: everything a fixture sizes in cells
/// must be held in metres) plus the grading knobs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradedMeshOptions {
    /// CPML margin added on each x/y side of the layout bbox, metres.
    pub margin_m: f64,
    /// Air above the trace plane, metres.
    pub air_above_m: f64,
    /// CPML absorber depth in cells. The `npml` outermost x/y cells are
    /// kept exactly coarse (the FS.0b.0 uniform-inside-absorbers scope
    /// rule, ADR-0208).
    pub npml: usize,
    /// Maximum cell-to-cell growth ratio of the geometric taper. Must be
    /// in `(1, 1.3]` — 1.3 is the compute-019-certified regime (measured
    /// taper reflection −52.7 dB, ADR-0208).
    pub growth: f64,
    /// Fine half-band around every trace edge, metres. The fringing-field
    /// length scale is the substrate height, the
    /// [`GradedMeshOptions::for_board`] default.
    pub guard_m: f64,
}

impl GradedMeshOptions {
    /// The FS.0a-board-shaped defaults: the `TwoPortBoardOptions::for_band`
    /// margins (34 cells) and absorber depth (10 cells) expressed in
    /// metres at the layout's own [`auto_dx_bulk`], growth 1.3, and a
    /// guard of one substrate height.
    pub fn for_board(layout: &Layout, f_max_hz: f64) -> Self {
        let d0 = auto_dx_bulk(layout, f_max_hz);
        Self {
            margin_m: 34.0 * d0,
            air_above_m: 34.0 * d0,
            npml: 10,
            growth: 1.3,
            guard_m: layout.substrate.height_m,
        }
    }
}

/// The output of [`auto_spacings`]: per-axis primal cell widths plus the
/// origin / z-stack metadata `yee_voxel::voxelize_microstrip_graded`
/// needs, and the two rule outcomes for reports.
#[derive(Debug, Clone, PartialEq)]
pub struct AutoSpacings {
    /// Primal cell widths along x (metres).
    pub dx: Vec<f64>,
    /// Primal cell widths along y (metres).
    pub dy: Vec<f64>,
    /// Primal cell widths along z (metres).
    pub dz: Vec<f64>,
    /// Layout-frame x of grid node 0 (`bbox.min.x − margin_m`, exactly).
    pub x0_m: f64,
    /// Layout-frame y of grid node 0.
    pub y0_m: f64,
    /// Ground-sheet layer index (always 0: the classic floor-ground stack).
    pub k_gnd: usize,
    /// Trace-plane layer index; the substrate fills `k = 0 .. k_top` with
    /// cells of exactly `h / k_top` (the ADR-0108 no-air-gap z-stack).
    pub k_top: usize,
    /// The coarse ceiling the rules produced (= [`auto_dx_bulk`]), metres.
    pub coarse_m: f64,
    /// The fine spacing inside feature bands, metres.
    pub fine_m: f64,
}

impl AutoSpacings {
    /// The `JobSpec`-ready spacing arrays.
    pub fn to_spacings(&self) -> crate::GradedSpacings {
        crate::GradedSpacings {
            dx: self.dx.clone(),
            dy: self.dy.clone(),
            dz: self.dz.clone(),
        }
    }

    /// Total cell count `nx · ny · nz`.
    pub fn cell_count(&self) -> usize {
        self.dx.len() * self.dy.len() * self.dz.len()
    }
}

/// The graded meshing rulebook (FS.0b.1, ADR-0210): per-axis spacing
/// arrays for a layout, refining only where the geometry needs it.
///
/// Rules — the FS.0a [`auto_dx`] rulebook generalized per axis:
///
/// - **Coarse ceiling everywhere:** `coarse = auto_dx_bulk(layout, f_max_hz)`
///   (λ/20-in-dielectric, h/3, clamped) — the feature rule moves into the
///   fine bands, so a narrow gap refines locally instead of dragging the
///   whole domain down.
/// - **Fine bands (x, y):** `fine = min(min_feature/2, coarse/2)` inside
///   `edge ± guard_m` around every trace-AABB edge and across every
///   inter-trace axis gap. The `coarse/2` term is measured, not assumed:
///   a feature-rule-only fine spacing leaves the FS.0a stub board
///   entirely unrefined (min_feature/2 = 1.5 mm > coarse), i.e. the
///   uniform pass-0 mesh whose notch ADR-0204 measured 5.2 % off the
///   converged answer; one halving at the staircase-limited edges is
///   what the uniform convergence trajectory showed sufficient.
/// - **Grading:** geometric ladder with ratio ≤ [`GradedMeshOptions::growth`]
///   between fine and coarse, junction steps included (compute-019-
///   certified regime).
/// - **z:** the substrate is `ceil(h / (coarse/2))` cells of exactly
///   `h / n_sub`; the air above grows geometrically to coarse and stays
///   coarse until `air_above_m` is covered.
/// - **Absorbers:** the `npml` outermost x/y cells are bit-equal coarse
///   (FS.0b.0 scope rule); `JobSpec::dx_m` should stay `coarse`, the
///   nominal spacing the CPML σ_max recipe assumes.
///
/// # Errors
///
/// Returns a message when the options are non-physical (growth outside
/// `(1, 1.3]`, non-positive margins/guard) or when a fine band would run
/// into a CPML absorber (widen `margin_m` — grading inside absorbers is
/// rejected by the kernel, ADR-0208).
pub fn auto_spacings(
    layout: &Layout,
    f_max_hz: f64,
    opts: &GradedMeshOptions,
) -> Result<AutoSpacings, String> {
    if !(opts.growth > 1.0 && opts.growth <= 1.3) {
        return Err(format!(
            "growth ratio {} outside (1, 1.3] (1.3 is the compute-019-certified regime)",
            opts.growth
        ));
    }
    for (name, v) in [
        ("margin_m", opts.margin_m),
        ("air_above_m", opts.air_above_m),
        ("guard_m", opts.guard_m),
    ] {
        if !(v.is_finite() && v > 0.0) {
            return Err(format!("{name} must be positive and finite (got {v})"));
        }
    }
    let coarse = auto_dx_bulk(layout, f_max_hz);
    let fine = (min_feature_m(layout) / 2.0).min(coarse / 2.0).max(1e-6);

    // Fine intervals per axis: every trace-AABB edge ± guard, plus every
    // inter-trace axis gap (the min_feature_m pair idiom).
    let boxes = trace_boxes(layout);
    let g = opts.guard_m;
    let mut ivx: Vec<(f64, f64)> = Vec::new();
    let mut ivy: Vec<(f64, f64)> = Vec::new();
    for &(x0, y0, x1, y1) in &boxes {
        ivx.push((x0 - g, x0 + g));
        ivx.push((x1 - g, x1 + g));
        ivy.push((y0 - g, y0 + g));
        ivy.push((y1 - g, y1 + g));
    }
    for (a, &(ax0, ay0, ax1, ay1)) in boxes.iter().enumerate() {
        for &(bx0, by0, bx1, by1) in boxes.iter().skip(a + 1) {
            if ay1.min(by1) - ay0.max(by0) > 0.0 {
                if bx0 - ax1 > 0.0 {
                    ivx.push((ax1, bx0));
                } else if ax0 - bx1 > 0.0 {
                    ivx.push((bx1, ax0));
                }
            }
            if ax1.min(bx1) - ax0.max(bx0) > 0.0 {
                if by0 - ay1 > 0.0 {
                    ivy.push((ay1, by0));
                } else if ay0 - by1 > 0.0 {
                    ivy.push((by1, ay0));
                }
            }
        }
    }

    let x0_m = layout.bbox.min.x - opts.margin_m;
    let y0_m = layout.bbox.min.y - opts.margin_m;
    let dx = mesh_axis(
        "x",
        x0_m,
        layout.bbox.max.x + opts.margin_m,
        coarse,
        fine,
        opts.npml,
        opts.growth,
        &ivx,
    )?;
    let dy = mesh_axis(
        "y",
        y0_m,
        layout.bbox.max.y + opts.margin_m,
        coarse,
        fine,
        opts.npml,
        opts.growth,
        &ivy,
    )?;
    let (dz, k_top) = mesh_z(
        layout.substrate.height_m,
        coarse,
        opts.growth,
        opts.air_above_m,
    );
    Ok(AutoSpacings {
        dx,
        dy,
        dz,
        x0_m,
        y0_m,
        k_gnd: 0,
        k_top,
        coarse_m: coarse,
        fine_m: fine,
    })
}

/// Ascending geometric ladder strictly between `fine` and `coarse`:
/// `fine·g, fine·g², …` while `< coarse`. Every junction ratio is ≤ `g`
/// by construction (the first would-be size ≥ coarse caps the last step).
fn taper_ladder(fine: f64, coarse: f64, growth: f64) -> Vec<f64> {
    let mut sizes = Vec::new();
    let mut d = fine * growth;
    while d < coarse {
        sizes.push(d);
        d *= growth;
    }
    sizes
}

/// Mesh one axis: march `lo → hi` with `npml` coarse absorber cells at
/// each end, `fine` cells across every (merged) fine interval, geometric
/// tapers between fine and coarse, and coarse fill elsewhere. The fine
/// band may start up to one coarse cell early (the coarse-fill `floor`
/// leftover) — coverage is extended, never clipped; the total length
/// overshoots `hi` by less than one coarse cell (the uniform voxelizer's
/// `ceil` behaviour).
#[allow(clippy::too_many_arguments)]
fn mesh_axis(
    axis: &str,
    lo: f64,
    hi: f64,
    coarse: f64,
    fine: f64,
    npml: usize,
    growth: f64,
    intervals: &[(f64, f64)],
) -> Result<Vec<f64>, String> {
    let ladder = taper_ladder(fine, coarse, growth);
    let ladder_len: f64 = ladder.iter().sum();

    // Merge overlapping / near intervals (a gap that cannot hold both
    // tapers plus slack is bridged with fine cells by merging).
    let mut iv: Vec<(f64, f64)> = intervals.iter().copied().filter(|(a, b)| b > a).collect();
    iv.sort_by(|p, q| p.0.total_cmp(&q.0));
    let mut merged: Vec<(f64, f64)> = Vec::new();
    for (a, b) in iv {
        match merged.last_mut() {
            Some(last) if a <= last.1 + 2.0 * ladder_len + coarse => last.1 = last.1.max(b),
            _ => merged.push((a, b)),
        }
    }

    if merged.is_empty() {
        let n = (((hi - lo) / coarse).ceil() as usize).max(2 * npml + 1);
        return Ok(vec![coarse; n]);
    }

    // Fine bands must clear the absorbers: grading inside CPML layers is
    // rejected by the kernel (FS.0b.0 scope rule, ADR-0208).
    let clear = npml as f64 * coarse + ladder_len + coarse;
    if merged[0].0 < lo + clear || merged.last().unwrap().1 > hi - clear {
        return Err(format!(
            "fine band [{:.4}, {:.4}] mm on {axis} runs into the CPML absorber \
             margin ([{:.4}, {:.4}] mm domain, {npml} coarse layers + taper); \
             widen margin_m",
            merged[0].0 * 1e3,
            merged.last().unwrap().1 * 1e3,
            lo * 1e3,
            hi * 1e3
        ));
    }

    let mut cells: Vec<f64> = Vec::new();
    let mut p = lo;
    let push = |cells: &mut Vec<f64>, p: &mut f64, d: f64, n: usize| {
        for _ in 0..n {
            cells.push(d);
            *p += d;
        }
    };
    push(&mut cells, &mut p, coarse, npml);
    let mut at_fine = false;
    for &(a, b) in &merged {
        // The gap p → a: up-taper (when leaving a fine band), coarse
        // fill, down-taper landing at or just before `a`.
        let up_len = if at_fine { ladder_len } else { 0.0 };
        let n_c = ((a - p - up_len - ladder_len) / coarse).floor().max(0.0) as usize;
        if at_fine {
            for &d in &ladder {
                push(&mut cells, &mut p, d, 1);
            }
        }
        push(&mut cells, &mut p, coarse, n_c);
        for &d in ladder.iter().rev() {
            push(&mut cells, &mut p, d, 1);
        }
        // Fine cells covering [p, b] (p ≤ a plus float slack).
        let n_f = (((b - p) / fine).ceil() as usize).max(1);
        push(&mut cells, &mut p, fine, n_f);
        at_fine = true;
    }
    // Tail: up-taper, then coarse through the far absorber.
    for &d in &ladder {
        push(&mut cells, &mut p, d, 1);
    }
    let n_tail = (((hi - p) / coarse).ceil() as usize).max(npml);
    push(&mut cells, &mut p, coarse, n_tail);
    Ok(cells)
}

/// Mesh the z axis: `n_sub = ceil(h / (coarse/2))` substrate cells of
/// exactly `h / n_sub` (`k_top = n_sub`; the ADR-0108 no-air-gap
/// z-stack), then air growing geometrically to coarse until
/// `air_above_m` is covered (at least one air layer).
fn mesh_z(h: f64, coarse: f64, growth: f64, air_above_m: f64) -> (Vec<f64>, usize) {
    let n_sub = ((h / (coarse / 2.0)).ceil() as usize).max(1);
    let dz_sub = h / n_sub as f64;
    let mut dz = vec![dz_sub; n_sub];
    let mut air = 0.0;
    let mut d = dz_sub;
    while air < air_above_m || dz.len() == n_sub {
        d = (d * growth).min(coarse);
        dz.push(d);
        air += d;
    }
    (dz, n_sub)
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
/// returns the |S21| curve from the **launch-normalized double ratio**
/// `|T_dut| / |T_ref|` with `T = fwd_B/fwd_A` per run
/// ([`sparams::forward_transfer`], the R.2-validated observable).
///
/// Deliberately NOT the single ratio `fwd_B(dut)/fwd_B(ref)`: that assumes
/// both runs launch the same incident wave, and the automesh forensics
/// (ADR-0204) measured that assumption failing at fine dx — the DUT read a
/// clean-fit, non-physical +10.7 dB because its launch differed from the
/// reference's. Normalizing each run by its own plane-A forward wave
/// cancels the launch exactly; the reference division then removes the
/// plane-A→plane-B line factor.
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
    let transfer = |p: &[Vec<f64>]| {
        sparams::forward_transfer(
            [&p[0], &p[1], &p[2]],
            [&p[3], &p[4], &p[5]],
            dt,
            spacing,
            freqs_hz,
        )
    };
    let t_dut = transfer(&dut_p);
    let t_ref = transfer(&ref_p);
    Ok(t_dut
        .iter()
        .zip(&t_ref)
        .map(|(d, r)| 20.0 * (d.0.hypot(d.1) / r.0.hypot(r.1)).log10())
        .collect())
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
    let spacing_m = opts.spacing_cells as f64 * opts.dx_m;
    let mut passes: Vec<ConvergencePass> = Vec::new();
    let mut final_delta = f64::INFINITY;
    for _ in 0..max_passes {
        // Keep the physical time window constant as dx (and thus dt) shrink.
        opts.n_steps = (base_steps / opts.dx_m).round() as usize;
        opts.margin_cells = (margin_m / opts.dx_m).round() as usize;
        opts.air_above_cells = (air_above_m / opts.dx_m).round() as usize;
        opts.npml = (npml_m / opts.dx_m).round() as usize;
        // The probe-triple span is measurement geometry: βd must stay put
        // or the standing-wave fit's conditioning changes pass to pass.
        opts.spacing_cells = (spacing_m / opts.dx_m).round() as usize;
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

    // --- auto_spacings (FS.0b.1, ADR-0210) --------------------------------

    /// The FS.0a stub-notch shape: a long line with an offset open stub.
    fn stub_like_layout() -> Layout {
        layout_with(
            vec![
                Polygon::rect(0.0, 0.0, 100.0e-3, 3.0e-3),
                Polygon::rect(48.5e-3, 3.0e-3, 3.0e-3, 8.0e-3),
            ],
            4.4,
            1.6e-3,
        )
    }

    /// Node positions from a running sum (tolerance-based checks only).
    fn nodes(x0: f64, widths: &[f64]) -> Vec<f64> {
        let mut v = vec![x0];
        for d in widths {
            v.push(v.last().unwrap() + d);
        }
        v
    }

    /// Every cell intersecting `[a, b]` by more than float slack has
    /// width ≤ `fine` (with an ulp-scale tolerance).
    fn band_is_fine(widths: &[f64], node: &[f64], a: f64, b: f64, fine: f64) {
        for (i, d) in widths.iter().enumerate() {
            let overlap = node[i + 1].min(b) - node[i].max(a);
            if overlap > 1e-12 {
                assert!(
                    *d <= fine * (1.0 + 1e-9),
                    "cell {i} ([{:.4}, {:.4}] mm, width {:.4} mm) inside fine band \
                     [{:.4}, {:.4}] mm exceeds fine {:.4} mm",
                    node[i] * 1e3,
                    node[i + 1] * 1e3,
                    d * 1e3,
                    a * 1e3,
                    b * 1e3,
                    fine * 1e3
                );
            }
        }
    }

    #[test]
    fn growth_ratio_is_bounded_everywhere() {
        let l = stub_like_layout();
        let opts = GradedMeshOptions::for_board(&l, 6.0e9);
        let s = auto_spacings(&l, 6.0e9, &opts).unwrap();
        for (axis, arr) in [("dx", &s.dx), ("dy", &s.dy), ("dz", &s.dz)] {
            for (i, w) in arr.windows(2).enumerate() {
                let ratio = (w[1] / w[0]).max(w[0] / w[1]);
                assert!(
                    ratio <= opts.growth * (1.0 + 1e-9),
                    "{axis}[{i}→{}] ratio {ratio} exceeds growth {}",
                    i + 1,
                    opts.growth
                );
            }
        }
    }

    #[test]
    fn fine_bands_cover_every_trace_edge_and_spacings_are_valid() {
        let l = stub_like_layout();
        let opts = GradedMeshOptions::for_board(&l, 6.0e9);
        let s = auto_spacings(&l, 6.0e9, &opts).unwrap();
        // fine = min(min_feature/2, coarse/2) = coarse/2 here (3 mm traces).
        assert!((s.fine_m - s.coarse_m / 2.0).abs() < 1e-15);
        let xn = nodes(s.x0_m, &s.dx);
        let yn = nodes(s.y0_m, &s.dy);
        let g = opts.guard_m;
        for e in [0.0, 100.0e-3, 48.5e-3, 51.5e-3] {
            band_is_fine(&s.dx, &xn, e - g, e + g, s.fine_m);
        }
        for e in [0.0, 3.0e-3, 11.0e-3] {
            band_is_fine(&s.dy, &yn, e - g, e + g, s.fine_m);
        }
        // Every width positive/finite and within [substrate-fine, coarse].
        for arr in [&s.dx, &s.dy, &s.dz] {
            assert!(
                arr.iter()
                    .all(|d| d.is_finite() && *d > 0.0 && *d <= s.coarse_m * (1.0 + 1e-12))
            );
        }
    }

    #[test]
    fn absorber_layers_are_uniform_coarse_and_length_bookkeeping_holds() {
        let l = stub_like_layout();
        let opts = GradedMeshOptions::for_board(&l, 6.0e9);
        let s = auto_spacings(&l, 6.0e9, &opts).unwrap();
        for (axis, arr, lo, hi) in [
            (
                "x",
                &s.dx,
                l.bbox.min.x - opts.margin_m,
                l.bbox.max.x + opts.margin_m,
            ),
            (
                "y",
                &s.dy,
                l.bbox.min.y - opts.margin_m,
                l.bbox.max.y + opts.margin_m,
            ),
        ] {
            // Bit-equal coarse inside both absorbers (FS.0b.0 scope rule).
            assert!(
                arr[..opts.npml].iter().all(|d| *d == s.coarse_m),
                "{axis}-min absorber not uniform coarse"
            );
            assert!(
                arr[arr.len() - opts.npml..]
                    .iter()
                    .all(|d| *d == s.coarse_m),
                "{axis}-max absorber not uniform coarse"
            );
            // Sum spans the domain: covers [lo, hi], overshoot < one coarse
            // cell (the uniform voxelizer's ceil behaviour).
            let span: f64 = arr.iter().sum();
            let want = hi - lo;
            assert!(
                span >= want - 1e-12 && span < want + s.coarse_m,
                "{axis} span {span} vs domain {want}"
            );
        }
        assert_eq!(s.x0_m, l.bbox.min.x - opts.margin_m);
        assert_eq!(s.y0_m, l.bbox.min.y - opts.margin_m);
    }

    #[test]
    fn z_stack_is_exact_substrate_fill_plus_graded_air() {
        let l = stub_like_layout();
        let opts = GradedMeshOptions::for_board(&l, 6.0e9);
        let s = auto_spacings(&l, 6.0e9, &opts).unwrap();
        assert_eq!(s.k_gnd, 0);
        // Substrate: uniform cells summing to h exactly (h / n_sub each).
        let h = l.substrate.height_m;
        assert!(s.k_top >= 6, "expected ≥ 2× the h/3 substrate resolution");
        assert!(s.dz[..s.k_top].iter().all(|d| *d == s.dz[0]));
        let sub: f64 = s.dz[..s.k_top].iter().sum();
        assert!((sub - h).abs() < 1e-12 * h);
        // Air: at least one layer, total ≥ air_above_m, none above coarse.
        let air: f64 = s.dz[s.k_top..].iter().sum();
        assert!(s.dz.len() > s.k_top && air >= opts.air_above_m);
        assert!(s.dz[s.k_top..].iter().all(|d| *d <= s.coarse_m));
    }

    #[test]
    fn single_rect_layout_degenerates_to_near_uniform() {
        // One wide trace, no sub-coarse features: only the four edge
        // bands are fine; everything else is bit-equal coarse.
        let l = layout_with(vec![Polygon::rect(0.0, 0.0, 50e-3, 10e-3)], 4.4, 1.6e-3);
        let opts = GradedMeshOptions::for_board(&l, 6.0e9);
        let s = auto_spacings(&l, 6.0e9, &opts).unwrap();
        let xn = nodes(s.x0_m, &s.dx);
        let g = opts.guard_m;
        // Taper extent past a band: the ladder is shorter than 2 coarse.
        let slack = g + 2.0 * s.coarse_m;
        let edges = [0.0, 50e-3];
        for (i, d) in s.dx.iter().enumerate() {
            let (c0, c1) = (xn[i], xn[i + 1]);
            let near_edge = edges
                .iter()
                .any(|e| c1 > e - slack - s.coarse_m && c0 < e + slack + s.coarse_m);
            if !near_edge {
                assert!(
                    *d == s.coarse_m,
                    "cell {i} away from every edge is not coarse: {d}"
                );
            }
        }
        // Coarse cells dominate the count.
        let coarse_n = s.dx.iter().filter(|d| **d == s.coarse_m).count();
        assert!(
            2 * coarse_n > s.dx.len(),
            "coarse does not dominate: {coarse_n}/{}",
            s.dx.len()
        );
    }

    #[test]
    fn inter_trace_gap_is_fine_throughout() {
        // Two lines 0.4 mm apart in y: the gap binds the fine spacing
        // (min_feature/2 = 0.2 mm < coarse/2 = 0.267 mm) — but only
        // LOCALLY: the coarse ceiling stays h/3, the graded payoff over
        // the uniform rulebook, where feature/2 would cap the whole grid.
        let l = layout_with(
            vec![
                Polygon::rect(0.0, 0.0, 20e-3, 1.5e-3),
                Polygon::rect(0.0, 1.9e-3, 20e-3, 1.5e-3),
            ],
            4.4,
            1.6e-3,
        );
        let opts = GradedMeshOptions::for_board(&l, 6.0e9);
        let s = auto_spacings(&l, 6.0e9, &opts).unwrap();
        assert!((s.fine_m - 0.2e-3).abs() < 1e-15, "fine = {}", s.fine_m);
        assert!(
            (s.coarse_m - 1.6e-3 / 3.0).abs() < 1e-15,
            "coarse must stay the bulk h/3 ceiling, got {}",
            s.coarse_m
        );
        let yn = nodes(s.y0_m, &s.dy);
        band_is_fine(&s.dy, &yn, 1.5e-3, 1.9e-3, s.fine_m);
    }

    #[test]
    fn bad_options_and_absorber_collisions_are_errors() {
        let l = stub_like_layout();
        let mut opts = GradedMeshOptions::for_board(&l, 6.0e9);
        opts.growth = 1.5;
        assert!(
            auto_spacings(&l, 6.0e9, &opts)
                .unwrap_err()
                .contains("growth")
        );
        // A guard so wide the band reaches the absorber must error, not
        // silently grade inside the CPML (the kernel rejects it anyway).
        let mut opts = GradedMeshOptions::for_board(&l, 6.0e9);
        opts.guard_m = opts.margin_m;
        assert!(
            auto_spacings(&l, 6.0e9, &opts)
                .unwrap_err()
                .contains("absorber")
        );
    }
}
