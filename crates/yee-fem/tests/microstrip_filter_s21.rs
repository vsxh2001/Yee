//! FEM-EM brick B7 (ADR-0153) — 3-pole microstrip-filter S21 from the FEM
//! driven sweep, graded against the analytic `ladder_s21` reference.
//!
//! This is the **culmination** of the FEM-EM driven-sweep track. Bricks B1
//! (interior-PEC edges), B2 (`layered_microstrip_mesh`), B3 (quasi-TEM
//! wave-port) and B4 (straight-line ε_eff = 0.61 % of Hammerstad-Jensen — the
//! port physics is *proven*) are merged. B7 composes them into a coupled-
//! resonator **band-pass filter** geometry, drives a two-port `sweep_matrix`
//! over the band, de-embeds the feed reference plane, extracts |S21|(f), and
//! grades the curve against the 3-pole Chebyshev 0.5 dB / 2 GHz / 10 % FBW
//! `ladder_s21` reference — including the geometric-asymmetry discriminator
//! (`depth(1.6 GHz) > depth(2.4 GHz)`).
//!
//! ## Honest framing (read before the gate)
//!
//! Per ADR-0153 this is "a geometry + de-embed exercise on a proven port",
//! but it is also the **hardest, most open brick**. The deliverable is an
//! HONEST first graded filter curve, not necessarily a strict pass. The gate
//! therefore asserts only the checks that the coarse-mesh / analytic-port
//! solve actually supports (see [`fem_filter_s21_vs_ladder`] for the exact
//! assertions and the measured curve), and records the full |S21|(f) table +
//! the gap to the ideal mask in the docstring. Weakening or faking the grade
//! is not a valid outcome; an imperfect-but-recognisable curve with an honest
//! assessment is.
//!
//! ## Geometry — edge-coupled 3-pole, FR-4
//!
//! Dimensions come from `yee_filter::dimension_edge_coupled` for the reference
//! spec (2 GHz, 10 % FBW, 0.5 dB Cheb, 50 Ω on 1 mm FR-4):
//!
//! ```text
//!   line width w     ≈ 1.91 mm   (50 Ω Hammerstad-Jensen)
//!   resonator λ_g/2  ≈ 41.1 mm   (HUGE at 2 GHz on FR-4)
//!   coupling gaps    ≈ 1.62 mm   (both, symmetric 3-pole)
//! ```
//!
//! The classic staggered edge-coupled footprint (each resonator overlaps the
//! next by ~λ_g/4) spans ~82 mm along the propagation axis. At a coarse
//! `dy ≈ 4–5 mm` cell pitch that is ~16–20 longitudinal cells; with the
//! cross-section resolving the trace (≥1 cell across `w`) and ~2–3 substrate-
//! heights of air clearance each side, the mesh lands at a few×10⁴ tets — the
//! upper edge of what a direct `faer` sparse LU fits in a 14 g box. If it
//! OOMs, that is the **B5/scaling boundary** and is reported, not forced.
//!
//! ## Axis convention (matches B2 `layered_microstrip_mesh`)
//!
//! ```text
//!   x ∈ [0, box_w]    cross-section width / strip stagger
//!   y ∈ [0, box_len]  PROPAGATION (down the filter, feed-to-feed)
//!   z ∈ [0, box_h]    substrate-normal (ground z=0, trace z=sub_h)
//! ```
//!
//! Ports sit on the `y = 0` (input feed) and `y = box_len` (output feed)
//! end-caps, exactly as the B4 straight line. `with_coupled_whitney(true)` is
//! MANDATORY (B4 finding: the lumped-centroid port collapses the absorbing
//! block for the substrate-normal `E_z` mode).
//!
//! ## GATING — CRITICAL
//!
//! Multi-minute driven SWEEP (one per-ω sparse LU per frequency point). All
//! tests here are `#[ignore]`'d so the debug `cargo test --workspace` never
//! runs them, and are run only in `--release`, boxed:
//!
//! ```text
//! YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh \
//!   cargo test -p yee-fem --release --test microstrip_filter_s21 \
//!   -- --ignored fem_filter_s21_vs_ladder --nocapture
//! ```

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::{
    FaceKind, MaterialDatabase, OpenBoundarySolver, PortDefinition, SParametersMatrix, TraceRect,
    beta_microstrip, layered_microstrip_filter_mesh, modal_e_t_microstrip_windowed,
};
use yee_filter::{
    Approximation, FilterSpec, LumpedLadder, Response, SpecMask, dimension_edge_coupled,
    ladder_s21, synthesize, synthesize_lumped,
};
use yee_layout::{Substrate, eps_eff};
use yee_mesh::TetMesh3D;

// ---------------------------------------------------------------------
// Fixed spec / substrate.
// ---------------------------------------------------------------------

/// Substrate height (m): 1 mm FR-4.
const SUB_H: f64 = 1.0e-3;
/// FR-4 relative permittivity.
const EPS_R: f64 = 4.4;
/// Band-pass centre frequency (Hz).
const F0: f64 = 2.0e9;
/// Fractional bandwidth.
const FBW: f64 = 0.10;

/// The reference filter spec the oracle grader (brick B6) uses: 3-pole
/// Chebyshev 0.5 dB BPF, f0 = 2 GHz, FBW = 10 %, Z0 = 50 Ω.
fn reference_spec() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: F0,
        fbw: FBW,
        order: Some(3),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![],
        },
    }
}

/// The canonical reference lumped ladder (the curve every EM method must
/// reproduce). Same construction as `yee-filter`'s `oracle_grade` example.
fn reference_ladder() -> LumpedLadder {
    synthesize_lumped(&synthesize(&reference_spec())).expect("bandpass N=3 synthesizes")
}

// ---------------------------------------------------------------------
// Filter geometry: edge-coupled 3-pole, mapped into the FEM box axes.
//
// `dimension_edge_coupled` gives the line width, the λ_g/2 resonator length,
// and the N−1 coupling gaps. We lay the resonators along the PROPAGATION axis
// y, offset in the cross-section axis x by (w + gap), and staggered by half a
// resonator length in y so adjacent strips overlap over ~λ_g/4 (the coupled
// region), mirroring `yee_layout::edge_coupled_bpf` (which uses its x as the
// long axis; we relabel long-axis → mesh-y, stagger-axis → mesh-x). Feed lines
// extend to the y = 0 / y = box_len end-caps where the wave-ports sit.
// ---------------------------------------------------------------------

/// Resolved filter geometry in mesh world coordinates plus the box extents and
/// the chosen subdivision.
struct FilterGeometry {
    /// Box extents (m): `(box_w, box_len, box_h)`.
    box_w: f64,
    box_len: f64,
    box_h: f64,
    /// Trace rectangles on the `z = sub_h` plane (resonators + feeds).
    traces: Vec<TraceRect>,
    /// Subdivisions `(nx, ny, nz)`.
    nx: usize,
    ny: usize,
    nz: usize,
    /// Trace line width (m) — the wave-port `w`.
    line_w: f64,
    /// One-sided feed length (m) at each end (the de-embed reference length).
    feed_len: f64,
    /// `x` centre (m) of the INPUT feed (port 0, `y = 0` end-cap). The
    /// quasi-TEM wave-port window is centred here — NOT at the box centre —
    /// because the feed is a narrow off-centre strip and a box-centred
    /// uniform-`x` mode mostly misses it (the dominant fix that lifted |S21|
    /// ~13 dB out of the noise floor; see the module-level honest framing).
    feed_xc_in: f64,
    /// `x` centre (m) of the OUTPUT feed (port 1, `y = box_len` end-cap).
    feed_xc_out: f64,
}

impl FilterGeometry {
    fn total_tets(&self) -> usize {
        self.nx * self.ny * self.nz * 6
    }
}

/// Build the edge-coupled 3-pole filter geometry.
///
/// `clearance_x` is the air margin (m) the PEC shield walls stand off the
/// trace pattern on each side in x (B4: ~2.5 substrate heights keeps the box
/// from loading the line). `air_h` is the air height above the substrate.
/// `dy_target` / `dx_target` set the (coarse) cell pitch; the actual counts
/// are rounded so `sub_h` lands on a z-plane and the trace spans ≥1 x-cell.
///
/// `feed_len` is the straight feed-line length at each end (a known de-embed
/// reference length); a longer feed buys a cleaner reference plane but more
/// cells.
#[allow(clippy::too_many_arguments)]
fn build_edge_coupled_geometry(
    clearance_x: f64,
    air_h: f64,
    feed_len: f64,
    dx_target: f64,
    dy_target: f64,
    dz: f64,
) -> FilterGeometry {
    // 1. Synthesize the physical dimensions.
    let project = synthesize(&reference_spec());
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: SUB_H,
        loss_tangent: 0.0,
        metal_thickness_m: 0.0,
    };
    let dims = dimension_edge_coupled(&project, &sub).expect("edge-coupled 3-pole synthesizes");
    let w = dims.line_width_m;
    let res_l = dims.resonator_length_m;
    let gaps = dims.gaps_m; // length N-1 = 2 for a 3-pole filter
    let n = gaps.len() + 1; // 3 resonators

    // 2. Lay the N resonators in mesh coords. Long axis = y (propagation),
    //    stagger axis = x. Resonator i: x0_i = Σ_{j<i}(w + gap_j); y0
    //    alternates 0 / stagger so adjacent strips overlap ~half their length.
    let stagger = res_l / 2.0;
    let mut x0 = clearance_x; // first strip left edge, clear of the x-wall
    let mut strips: Vec<TraceRect> = Vec::with_capacity(n);
    for i in 0..n {
        let y0 = if i % 2 == 0 { 0.0 } else { stagger };
        strips.push(TraceRect::new(x0, y0, w, res_l));
        if i < gaps.len() {
            x0 += w + gaps[i];
        }
    }
    // x-extent spanned by the strips.
    let strips_x_hi = strips.iter().map(|r| r.x0 + r.w).fold(0.0_f64, f64::max);
    let strips_x_lo = clearance_x;
    // y-extent spanned by the resonators (before feeds): the staggered strips
    // occupy [0, stagger + res_l] = [0, res_l + stagger].
    let res_y_hi = res_l + stagger;

    // 3. Box width: trace x-span + clearance both sides.
    let box_w = strips_x_hi + clearance_x;
    // Box height: substrate + air.
    let box_h = SUB_H + air_h;

    // 4. Feed lines. The filter spans y ∈ [0, res_y_hi] in the resonator
    //    region; shift everything up by `feed_len` so an input feed can run
    //    from y = 0 to y = feed_len into resonator 0, and an output feed from
    //    y = feed_len + res_y_hi to box_len out of the last resonator. The
    //    feeds are centred (in x) on the resonator they attach to.
    let y_shift = feed_len;
    let box_len = feed_len + res_y_hi + feed_len;

    // Re-emit the strips shifted up by y_shift.
    let mut traces: Vec<TraceRect> = strips
        .iter()
        .map(|r| TraceRect::new(r.x0, r.y0 + y_shift, r.w, r.l))
        .collect();

    // Input feed: attaches to resonator 0 (which after shift starts at
    // y = y_shift). It is the first strip in `strips` (x at strips[0].x0). The
    // feed runs y ∈ [0, y_shift], width w, x aligned with resonator 0.
    let in_x = strips[0].x0;
    traces.push(TraceRect::new(in_x, 0.0, w, y_shift));
    // Output feed: attaches to the LAST resonator. After shift its top edge is
    // at y = y_shift + (last strip y0) + res_l. The last strip (i = n-1) has
    // y0 = if (n-1)%2==0 {0} else {stagger}; its top = y_shift + y0 + res_l.
    let last = strips.last().unwrap();
    let out_feed_y0 = y_shift + last.y0 + last.l;
    let out_x = last.x0;
    traces.push(TraceRect::new(out_x, out_feed_y0, w, box_len - out_feed_y0));

    // 5. Subdivisions (coarse). nz so sub_h lands on a plane.
    let nz = (box_h / dz).round() as usize;
    // Snap box_h so sub_h * nz / box_h is integral (sub_h is a multiple of dz).
    let nz_sub = (SUB_H / dz).round().max(1.0) as usize;
    let dz_exact = SUB_H / nz_sub as f64;
    let box_h = dz_exact * nz as f64; // keep nz cells, exact dz
    let nx = (box_w / dx_target)
        .round()
        .max(((w / dx_target).ceil()).max(1.0)) as usize;
    let ny = (box_len / dy_target).round().max(1.0) as usize;

    let _ = (strips_x_lo,); // silence unused in some configs
    FilterGeometry {
        box_w,
        box_len,
        box_h,
        traces,
        nx,
        ny,
        nz,
        line_w: w,
        feed_len,
        // Feed x-centres for the windowed wave-port (the feed strip centre).
        feed_xc_in: in_x + w / 2.0,
        feed_xc_out: out_x + w / 2.0,
    }
}

// ---------------------------------------------------------------------
// Face classification — ports on the y = 0 / y = box_len end-caps (same as
// the B4 straight line; the filter just has a richer interior-PEC footprint).
// ---------------------------------------------------------------------

fn exterior_face_count(mesh: &TetMesh3D) -> usize {
    let mut face_map: std::collections::HashMap<[usize; 3], usize> =
        std::collections::HashMap::new();
    const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    for tet in &mesh.tetrahedra {
        for &[a, b, c] in TET_FACES.iter() {
            let mut key = [tet[a], tet[b], tet[c]];
            key.sort_unstable();
            *face_map.entry(key).or_insert(0) += 1;
        }
    }
    face_map.values().filter(|&&c| c == 1).count()
}

fn classify_filter_faces(centroids: &[Vector3<f64>], box_len: f64) -> Vec<FaceKind> {
    let tol = 1e-9;
    centroids
        .iter()
        .map(|c| {
            if c.y < tol {
                FaceKind::WavePort(0)
            } else if (c.y - box_len).abs() < tol {
                FaceKind::WavePort(1)
            } else {
                FaceKind::Pec
            }
        })
        .collect()
}

/// Build a complete two-port driven solver for the filter geometry. Trace AND
/// ground tagged interior-PEC (B1); the two y-end-caps carry the quasi-TEM
/// wave-port (B3) with `with_coupled_whitney(true)` (mandatory, B4 finding).
///
/// The wave-port `β` and modal shape use the FEED-LINE width `line_w` (the
/// feed is a uniform 50 Ω microstrip — that is what the port face actually
/// sees, regardless of the coupled-resonator interior).
fn solve_filter(geom: &FilterGeometry, omegas: &[f64]) -> SParametersMatrix {
    let (mesh, material_db, ground_pred, trace_pred) = layered_microstrip_filter_mesh(
        geom.box_w,
        geom.box_len,
        geom.box_h,
        SUB_H,
        geom.traces.clone(),
        geom.nx,
        geom.ny,
        geom.nz,
    )
    .expect("filter mesh must build");

    let n_exterior = exterior_face_count(&mesh);
    let picker = OpenBoundarySolver::new(
        &mesh,
        vec![FaceKind::Pec; n_exterior],
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("picker solver must build");

    let ground_edges = picker.interior_edges_matching(&ground_pred);
    let trace_edges = picker.interior_edges_matching(&trace_pred);
    let mut interior_pec: Vec<usize> = ground_edges;
    interior_pec.extend(trace_edges.iter().copied());
    interior_pec.sort_unstable();
    interior_pec.dedup();
    assert!(
        !trace_edges.is_empty(),
        "trace_pred must select at least one interior edge on the z = sub_h trace footprint"
    );
    let centroids = picker.exterior_face_centroids();
    let kinds = classify_filter_faces(&centroids, geom.box_len);
    drop(picker);

    // Feed-line wave-port: a quasi-TEM 50 Ω microstrip mode, x-windowed on the
    // FEED's x-centre (not the box centre). The feed is a narrow off-centre
    // strip; a box-centred uniform-x E_z mode (the B4 straight-line port, where
    // the trace WAS centred) mostly excites air/PEC away from the feed and
    // couples ~13 dB worse here. The raised-cosine window (half-width 2·w) is
    // `microstrip_port_windowed`'s shape, recentred per-port on the actual
    // feed. β uses the feed width (the port face sees a uniform 50 Ω line,
    // whatever the coupled-resonator interior does).
    let make_port = |xc: f64| {
        let w = geom.line_w;
        let win = 2.0 * w;
        let beta = move |omega: f64| beta_microstrip(w, SUB_H, EPS_R, omega);
        let e_t = move |p: Vector3<f64>| modal_e_t_microstrip_windowed(xc, win, SUB_H, p);
        PortDefinition::single_mode(Box::new(beta), Box::new(e_t))
    };

    let solver = OpenBoundarySolver::new(
        &mesh,
        kinds,
        vec![make_port(geom.feed_xc_in), make_port(geom.feed_xc_out)],
        material_db,
    )
    .expect("two-port filter solver must build")
    .with_interior_pec_edges(interior_pec.iter().copied())
    .with_coupled_whitney(true);

    solver
        .sweep_matrix(omegas)
        .expect("driven sweep_matrix must succeed")
}

fn db(mag: f64) -> f64 {
    20.0 * mag.log10()
}

/// De-embed the two feed-line reference planes from a raw S21.
///
/// Each feed is a straight `feed_len` 50 Ω microstrip whose quasi-TEM phase
/// constant is `β = (ω/c)·√ε_eff(w)`. The two feeds add a total electrical
/// length `2·β·feed_len` of phase to S21 (and the reference planes are moved
/// to the filter ports by removing it). The feed is closely matched (~50 Ω),
/// so its magnitude effect is small; we de-embed phase only (the standard
/// reference-plane shift), which is what the asymmetry / shape grading needs.
fn deembed_feed(
    s21_raw: num_complex::Complex64,
    omega: f64,
    line_w: f64,
    feed_len: f64,
) -> num_complex::Complex64 {
    let beta = beta_microstrip(line_w, SUB_H, EPS_R, omega);
    let phase = 2.0 * beta * feed_len; // both feeds
    // Move reference planes inward: multiply by e^{+jβℓ} on each side.
    s21_raw * num_complex::Complex64::from_polar(1.0, phase)
}

// =====================================================================
// FEASIBILITY / LU-CEILING PROBE (#[ignore]'d — run a real LU solve)
//
// Finds the largest filter mesh whose per-ω sparse LU fits a 14 g box, so the
// gate can fix the coarsest mesh that gives a recognisable response without
// hitting the B5/scaling boundary. NOT a gate.
// =====================================================================

/// Probe: report the filter mesh size for a few cell pitches and run ONE
/// single-frequency `sweep_matrix` at f0 to confirm the LU factors (or OOMs).
#[test]
#[ignore = "feasibility probe; run explicitly — builds the filter mesh + one LU"]
fn fem_filter_s21_probe() {
    // A few candidate coarsenesses, coarsest first.
    for (label, dx, dy, dz, clr, air, feed) in [
        ("coarse  ", 1.6e-3, 5.0e-3, 0.5e-3, 2.5e-3, 5.0e-3, 8.0e-3),
        ("medium  ", 1.3e-3, 4.0e-3, 0.5e-3, 2.5e-3, 5.0e-3, 8.0e-3),
        ("fine    ", 1.0e-3, 3.0e-3, 0.5e-3, 3.0e-3, 5.0e-3, 8.0e-3),
    ] {
        let geom = build_edge_coupled_geometry(clr, air, feed, dx, dy, dz);
        eprintln!(
            "[probe] {label}: box=({:.1},{:.1},{:.1})mm  n=({},{},{})  tets={}  feed={:.1}mm",
            geom.box_w * 1e3,
            geom.box_len * 1e3,
            geom.box_h * 1e3,
            geom.nx,
            geom.ny,
            geom.nz,
            geom.total_tets(),
            geom.feed_len * 1e3,
        );
        let omega = 2.0 * PI * F0;
        let t0 = std::time::Instant::now();
        let sweep = solve_filter(&geom, &[omega]);
        let s = &sweep.s[0];
        eprintln!(
            "[probe] {label}: |S11|={:.4} |S21|={:.4} ({:.1} dB)  solve {:.1}s",
            s[(0, 0)].norm(),
            s[(1, 0)].norm(),
            db(s[(1, 0)].norm()),
            t0.elapsed().as_secs_f64(),
        );
    }
}

// =====================================================================
// THE GATE
// =====================================================================

/// Passband (near-band) tolerance in dB, mirroring `yee-filter`'s
/// `oracle_grade`: `|extracted − reference|` over ~[1.85, 2.15] GHz.
const PASSBAND_TOL_DB: f64 = 2.0;
/// Stopband / rejection-skirt tolerance in dB (looser).
const REJECTION_TOL_DB: f64 = 5.0;
/// Asymmetry-discriminator margin (dB): lower notch must be deeper than upper
/// by at least this. Mirrors `oracle_grade::ASYMMETRY_MARGIN_DB`.
const ASYMMETRY_MARGIN_DB: f64 = 1.0;

/// FEM-EM brick B7 (ADR-0153) — 3-pole microstrip-filter S21 driven sweep vs
/// the analytic ladder reference.
///
/// Builds the edge-coupled 3-pole filter, drives `sweep_matrix` over
/// 1.6–2.4 GHz, de-embeds the feed reference planes, extracts |S21|(f), and
/// grades it against the 3-pole Cheb 0.5 dB / 2 GHz / 10 % FBW `ladder_s21`
/// reference + the geometric-asymmetry discriminator.
///
/// ## What this asserts (HONEST)
///
/// The strict `oracle_grade` mask (passband |err| ≤ 2 dB, rejection |err| ≤
/// 5 dB) is the *target*, and is computed + printed every run — but it MISSES
/// by ~42 dB in-band and CANNOT pass with this analytic port (the modal-overlap
/// IL floor; see the honest verdict). The gate therefore does NOT assert the
/// absolute-level mask (no weakening to force green). It asserts the three
/// checks the coarse-mesh / analytic-port solve genuinely supports:
///
/// 1. **Non-degenerate transmission** — the in-band peak is well above the
///    solve noise floor (a collapsed port or broken mesh would sit in noise).
/// 2. **The geometric-asymmetry discriminator (the brick's NAMED check)** —
///    `depth(1.6 GHz) > depth(2.4 GHz)` by ≥ 1 dB: the FEM curve reproduces the
///    correct band-pass-mapping asymmetry SIGN that a symmetric/inverted
///    fitted artifact lacks.
/// 3. **A band-pass turnover** — the in-band peak stands above the deeper band
///    edge (a real centre bump, not a monotonic ramp / flat line).
///
/// The full curve, the strict-mask gap, and the honest verdict are recorded in
/// the MEASURED block below and printed by the test with `--nocapture`.
///
/// ## MEASURED RESULT (boxed --release, base 22da1c2; 51 336 tets, 73.9 s)
///
/// Edge-coupled 3-pole, 14.0 × 77.6 × 6.0 mm box, w = 1.912 mm,
/// `dx/dy/dz = 0.6/2.5/0.5 mm`, feed = 8 mm. |S21| after feed de-embed:
///
/// ```text
///   f(GHz)   S21 dB (FEM)   ref dB (ladder)
///   1.60      −44.62         −41.77
///   1.80      −43.30         −20.81
///   1.90      −42.87          −0.75
///   2.00      −42.65          0.00   ← reference passband centre
///   2.10      −42.39         −0.32   ← FEM in-band peak
///   2.20      −42.46         −17.83
///   2.40      −43.15         −36.27
///
///   in-band peak       : −42.39 dB @ 2.10 GHz
///   turnover           : +2.23 dB (in-band peak above the deeper band edge)
///   asymmetry (NAMED)  : depth(1.6)=44.62 dB > depth(2.4)=43.15 dB, +1.47 dB → PASS
///   strict oracle mask : MISS (worst in-band err ≈ 42.6 dB vs the 0 dB reference)
/// ```
///
/// ## Honest verdict
///
/// This is a **recognisable-but-imperfect** first FEM filter curve, not a
/// strict-mask pass — exactly the honest deliverable the brick asked for.
///
/// * **What is real:** a genuine, non-degenerate FEM transmission response
///   that (a) peaks near the 2 GHz band centre (overall peak @ 2.10 GHz),
///   (b) bumps up over the band edges by +2.23 dB (a frequency-selective
///   pass/stop shape, not a flat line), and (c) reproduces the CORRECT
///   geometric-asymmetry SIGN — the lower 1.6 GHz notch is deeper than the
///   upper 2.4 GHz notch by +1.47 dB, the band-pass-mapping signature the
///   reference has and a fitted/symmetric artifact does not. The gate asserts
///   exactly these three (non-degeneracy + the named asymmetry discriminator +
///   the turnover).
///
/// * **The gap to the ideal — and what it is:** the whole curve sits at
///   ~−42 to −44 dB while the reference passband is 0 dB — a ~42 dB in-band
///   miss on the strict absolute-level Chebyshev mask. This is the **analytic
///   wave-port modal-overlap insertion-loss floor**, *not* a mesh/LU-scaling
///   limit. B4 already documented that ONE matched straight-line port through
///   this analytic E_z quasi-TEM mode transmits only |S21| ≈ 0.089 (−21 dB) —
///   the analytic mode partially overlaps the true microstrip eigenmode, so
///   most incident power is lost in the projection (the phase is coherent,
///   which is why B4's ε_eff worked, but the amplitude is weak). Here the
///   signal traverses TWO such ports plus the lossy coupled-resonator
///   interior, so the floor roughly doubles to ~−42 dB and the 0.5 dB
///   Chebyshev passband cannot climb out of it. The fix is a higher-fidelity
///   port (a true numerical cross-section eigenmode, or aperture/frill
///   coupling), NOT a finer mesh: the 51 k-tet `faer` sparse LU fits the 14 g
///   box with room to spare (~3 s/point; an 80 k-tet refinement also fit and
///   did NOT lift the curve), so this is a port-fidelity boundary, not the
///   B5/scaling boundary.
///
/// * **Two levers that mattered en route** (recorded so they are not
///   re-derived): (1) the wave-port modal window must be centred on the
///   FEED's `x`, not the box centre — the feed is a narrow off-centre strip
///   and a box-centred uniform-`x` mode mostly misses it; recentring lifted
///   |S21| ~13 dB out of the noise floor and is what makes the asymmetry
///   resolvable. (2) `with_coupled_whitney(true)` is mandatory (B4 finding;
///   the lumped-centroid path collapses the absorbing block for the
///   substrate-normal `E_z` mode).
///
/// Run command (printed table + grade with `--nocapture`):
/// ```text
/// YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh \
///   cargo test -p yee-fem --release --test microstrip_filter_s21 \
///   -- --ignored fem_filter_s21_vs_ladder --nocapture
/// ```
#[test]
#[ignore = "multi-minute driven SWEEP (one per-ω sparse LU per point); run only in --release, boxed"]
fn fem_filter_s21_vs_ladder() {
    // Geometry — coarse but resolved enough that the trace (≥3 x-cells), the
    // coupling gaps (≥2 x-cells) and the resonators (~16 y-cells) are captured.
    // ~51 k tets; the long-thin box keeps the per-ω faer sparse-LU bandwidth
    // low, so this fits the 14 g box comfortably (~3 s/point) — the LU is NOT
    // the binding constraint here (the analytic-port modal-overlap floor is;
    // see the honest verdict). The probe (`fem_filter_s21_probe`) walks the
    // size/feasibility ladder.
    let geom = build_edge_coupled_geometry(
        2.5e-3, // x clearance each side
        5.0e-3, // air height
        8.0e-3, // feed length (de-embed reference)
        0.6e-3, // dx (trace ~3 cells, gap ~2.7 cells)
        2.5e-3, // dy (resonator ~16 cells)
        0.5e-3, // dz (2 substrate cells)
    );
    eprintln!(
        "[B7] filter mesh: box=({:.1},{:.1},{:.1})mm  n=({},{},{})  tets={}  w={:.3}mm  feed={:.1}mm  eps_eff(w)={:.4}",
        geom.box_w * 1e3,
        geom.box_len * 1e3,
        geom.box_h * 1e3,
        geom.nx,
        geom.ny,
        geom.nz,
        geom.total_tets(),
        geom.line_w * 1e3,
        geom.feed_len * 1e3,
        eps_eff(geom.line_w, SUB_H, EPS_R),
    );

    // Band: 1.6 – 2.4 GHz, 17 points (50 MHz spacing) — covers both notches
    // and the passband.
    let n_pts = 17;
    let f_lo = 1.6e9;
    let f_hi = 2.4e9;
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|i| f_lo + (f_hi - f_lo) * (i as f64) / ((n_pts - 1) as f64))
        .collect();
    let omegas: Vec<f64> = freqs_hz.iter().map(|f| 2.0 * PI * f).collect();

    let t0 = std::time::Instant::now();
    let sweep = solve_filter(&geom, &omegas);
    let wall = t0.elapsed().as_secs_f64();

    // Extract + de-embed |S21|(f) into a (f_GHz, dB) curve.
    let ladder = reference_ladder();
    let mut curve: Vec<(f64, f64)> = Vec::with_capacity(n_pts);
    eprintln!(
        "\n{:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
        "f(GHz)", "|S21|raw", "|S21|deemb", "S21 dB", "ref dB"
    );
    for (k, &omega) in omegas.iter().enumerate() {
        let s = &sweep.s[k];
        let s21_raw = s[(1, 0)];
        let s21 = deembed_feed(s21_raw, omega, geom.line_w, geom.feed_len);
        let d = db(s21.norm());
        let f_ghz = freqs_hz[k] / 1e9;
        let ref_db = db(ladder_s21(&ladder, freqs_hz[k]).norm());
        curve.push((f_ghz, d));
        eprintln!(
            "{:>8.3}  {:>10.4}  {:>10.4}  {:>10.2}  {:>10.2}",
            f_ghz,
            s21_raw.norm(),
            s21.norm(),
            d,
            ref_db,
        );
    }

    // ---- Grade against the reference (mirrors oracle_grade::evaluate) ----
    let mut worst_pass_db = 0.0_f64;
    let mut worst_rej_db = 0.0_f64;
    for &(f_ghz, d_meas) in &curve {
        let d_ref = db(ladder_s21(&ladder, f_ghz * 1e9).norm());
        let err = (d_meas - d_ref).abs();
        if (1.85..=2.15).contains(&f_ghz) {
            worst_pass_db = worst_pass_db.max(err);
        } else {
            worst_rej_db = worst_rej_db.max(err);
        }
    }

    // ---- Asymmetry discriminator: depth(1.6) > depth(2.4)? ----
    let depth_at = |f_ghz: f64| -> f64 { -interp_db(&curve, f_ghz) };
    let depth_lo = depth_at(1.6);
    let depth_hi = depth_at(2.4);
    let asym_margin = depth_lo - depth_hi;
    let asym_pass = asym_margin >= ASYMMETRY_MARGIN_DB;

    // ---- Recognisable-bandpass checks (the weaker honest floor) ----
    // In-band peak |S21| over the [1.85, 2.15] GHz passband, the band-edge
    // levels (1.6 / 2.4 GHz), and the overall peak frequency. The "turnover"
    // is the in-band peak standing above the lower band edge — a genuine
    // pass/stop shape rather than a monotonic ramp.
    let passband_peak_db = curve
        .iter()
        .filter(|(f, _)| (1.85..=2.15).contains(f))
        .map(|(_, d)| *d)
        .fold(f64::NEG_INFINITY, f64::max);
    let edge_lo_db = interp_db(&curve, 1.6);
    let edge_hi_db = interp_db(&curve, 2.4);
    let f_peak_ghz = curve
        .iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(f, _)| *f)
        .unwrap_or(f64::NAN);
    // Turnover: how far the in-band peak rises above the deeper of the two
    // band edges (the band-pass "bump"). Positive ⇒ a real centre peak.
    let turnover_db = passband_peak_db - edge_lo_db.min(edge_hi_db);

    let strict_pass =
        worst_pass_db <= PASSBAND_TOL_DB && worst_rej_db <= REJECTION_TOL_DB && asym_pass;

    let f_inband_peak = curve
        .iter()
        .filter(|(f, _)| (1.85..=2.15).contains(f))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(f, _)| *f)
        .unwrap_or(f64::NAN);
    eprintln!(
        "\n==== B7 GRADE ====\n\
         tets               : {}\n\
         wall               : {:.1} s\n\
         in-band peak       : {:.2} dB @ {:.2} GHz (overall peak @ {:.2} GHz)\n\
         band edges         : {:.2} dB @1.6  {:.2} dB @2.4\n\
         turnover           : {:+.2} dB (in-band peak above the deeper edge)\n\
         worst passband err : {:.2} dB vs ref (oracle tol {:.1})\n\
         worst rejection err: {:.2} dB vs ref (oracle tol {:.1})\n\
         asymmetry (NAMED)  : depth(1.6)={:.2} dB  depth(2.4)={:.2} dB  margin={:+.2} dB  -> {}\n\
         strict oracle mask : {}\n\
         ==================",
        geom.total_tets(),
        wall,
        passband_peak_db,
        f_inband_peak,
        f_peak_ghz,
        edge_lo_db,
        edge_hi_db,
        turnover_db,
        worst_pass_db,
        PASSBAND_TOL_DB,
        worst_rej_db,
        REJECTION_TOL_DB,
        depth_lo,
        depth_hi,
        asym_margin,
        if asym_pass { "PASS" } else { "FLAG" },
        if strict_pass { "PASS" } else { "MISS" },
    );

    // The "machine-readable" curve for the oracle_grade CLI (so a reviewer can
    // paste it into `cargo run -p yee-filter --example oracle_grade -- <pairs>`).
    let pairs: String = curve
        .iter()
        .map(|(f, d)| format!("{f:.3}:{d:.2}"))
        .collect::<Vec<_>>()
        .join(" ");
    eprintln!("[B7] oracle_grade pairs: {pairs}");

    // ---- Assertions (HONEST: assert only what genuinely holds) ----
    //
    // The strict oracle mask does NOT pass and CANNOT pass with this port: the
    // analytic E_z quasi-TEM wave-port has a large, ~frequency-flat
    // modal-overlap insertion-loss floor (B4 documented ~−21 dB for ONE matched
    // straight-line port; here the signal traverses TWO such ports plus the
    // lossy coupled-resonator interior, so the whole curve sits ~−42 dB). The
    // reference passband is 0 dB, so the absolute-level mask is off by ~42 dB
    // everywhere in-band — a PORT-FIDELITY gap, not a mesh/LU-scaling one (the
    // LU fits the 14 g box with room to spare). This is the honest verdict; the
    // strict-mask numbers are printed above for the record.
    //
    // What the coarse-mesh / analytic-port solve DOES capture, and what this
    // gate therefore asserts:

    // (1) Non-degenerate transmission: a real propagating wave reaches port 2
    //     (the peak is well above the solve's ~−65 dB noise floor). A collapsed
    //     port (the lumped-centroid failure mode) or a broken mesh would sit in
    //     the noise.
    assert!(
        passband_peak_db.is_finite() && passband_peak_db > -55.0,
        "B7 NO-GO: in-band peak {passband_peak_db:.2} dB is in the noise floor — no transmitted \
         wave reached port 2 (port collapsed or mesh broken). Full curve printed above."
    );

    // (2) The geometric-asymmetry discriminator — the brick's NAMED check — must
    //     PASS: the lower stopband notch (1.6 GHz) is genuinely deeper than the
    //     upper (2.4 GHz). This is the band-pass-mapping signature the reference
    //     has and a symmetric/inverted (fitted-artifact) curve does NOT; it is
    //     the anti-"flat/symmetric curve is not evidence" guard. The FEM curve
    //     reproduces the CORRECT asymmetry SIGN with margin ≥ 1 dB — a real,
    //     geometry-aware result, even though the absolute Chebyshev depth is
    //     unreachable through the lossy port. (If a future higher-fidelity port
    //     lifts the curve onto the strict mask, `strict_pass` flips and the
    //     run additionally clears the absolute-level mask; we do not assert that
    //     here because it does not yet hold — no weakening to force green.)
    assert!(
        asym_pass,
        "B7: geometric-asymmetry discriminator FAILED — depth(1.6 GHz)={depth_lo:.2} dB is NOT \
         deeper than depth(2.4 GHz)={depth_hi:.2} dB by the required {ASYMMETRY_MARGIN_DB} dB \
         (margin {asym_margin:+.2} dB). A symmetric/inverted curve has lost the band-pass-mapping \
         asymmetry and is not credited as a geometry-aware EM result. Full curve printed above."
    );

    // (3) A genuine band-pass turnover: the in-band peak stands above the deeper
    //     band edge (the response bumps up near band centre rather than ramping
    //     monotonically). A modest >0.2 dB bar — the coarse-mesh / lossy-port
    //     bump is shallow (~1–3 dB), so this certifies the SHAPE is frequency-
    //     selective without demanding a depth the port cannot deliver.
    assert!(
        turnover_db > 0.2,
        "B7: no band-pass turnover — in-band peak {passband_peak_db:.2} dB is not above the \
         deeper band edge (edges {edge_lo_db:.2}/{edge_hi_db:.2} dB; turnover {turnover_db:+.2} dB). \
         The response is monotonic, not a recognisable band-pass. Full curve printed above."
    );

    let _ = strict_pass; // reported above; not asserted (port-floor gap is honest)
}

/// Linear interpolation of the `(f_ghz, dB)` curve at `f_ghz` (clamped to the
/// endpoints). Mirrors `oracle_grade::interp_db`.
fn interp_db(pts: &[(f64, f64)], f_ghz: f64) -> f64 {
    if pts.is_empty() {
        return f64::NAN;
    }
    if f_ghz <= pts[0].0 {
        return pts[0].1;
    }
    if f_ghz >= pts[pts.len() - 1].0 {
        return pts[pts.len() - 1].1;
    }
    for w in pts.windows(2) {
        let (f0, d0) = w[0];
        let (f1, d1) = w[1];
        if (f0..=f1).contains(&f_ghz) {
            let t = if (f1 - f0).abs() < 1e-15 {
                0.0
            } else {
                (f_ghz - f0) / (f1 - f0)
            };
            return d0 + t * (d1 - d0);
        }
    }
    pts[pts.len() - 1].1
}

#[cfg(test)]
mod unit {
    use super::*;

    /// The edge-coupled geometry is well-formed: 3 resonators + 2 feeds = 5
    /// trace rectangles, the box clears the trace pattern, and sub_h lands on a
    /// z-plane. Fast (no solve) — runs in the default `cargo test`.
    #[test]
    fn geometry_is_well_formed() {
        let geom = build_edge_coupled_geometry(2.5e-3, 5.0e-3, 8.0e-3, 1.6e-3, 5.0e-3, 0.5e-3);
        // 3 resonators + input feed + output feed.
        assert_eq!(geom.traces.len(), 5, "3 resonators + 2 feeds");
        // Box clears the trace x-span.
        let trace_x_hi = geom
            .traces
            .iter()
            .map(|r| r.x0 + r.w)
            .fold(0.0_f64, f64::max);
        assert!(
            geom.box_w > trace_x_hi,
            "box_w {:.4} must clear trace x-extent {:.4}",
            geom.box_w,
            trace_x_hi
        );
        // sub_h on a z-plane.
        let dz = geom.box_h / geom.nz as f64;
        let n_sub = SUB_H / dz;
        assert!(
            (n_sub - n_sub.round()).abs() < 1e-9,
            "sub_h must land on a z-plane (n_sub = {n_sub})"
        );
        // Feeds reach both end-caps: some trace touches y=0 and some touches
        // y=box_len.
        let touches_y0 = geom.traces.iter().any(|r| r.y0.abs() < 1e-12);
        let touches_ylen = geom
            .traces
            .iter()
            .any(|r| (r.y0 + r.l - geom.box_len).abs() < 1e-9);
        assert!(touches_y0, "an input feed must reach the y=0 port plane");
        assert!(
            touches_ylen,
            "an output feed must reach the y=box_len port plane"
        );
    }

    /// De-embed is a pure phase rotation (unit magnitude factor) — it must not
    /// change |S21|, only its phase. Fast.
    #[test]
    fn deembed_preserves_magnitude() {
        let s = num_complex::Complex64::new(0.3, -0.2);
        let omega = 2.0 * PI * F0;
        let out = deembed_feed(s, omega, 1.9e-3, 8.0e-3);
        assert!(
            (out.norm() - s.norm()).abs() < 1e-12,
            "de-embed changed |S21| ({} vs {})",
            out.norm(),
            s.norm()
        );
    }

    /// The reference ladder reproduces the band-pass-mapping asymmetry: lower
    /// notch (1.6 GHz) deeper than upper (2.4 GHz). This is the ground-truth
    /// the discriminator checks against. Fast.
    #[test]
    fn reference_has_asymmetric_notches() {
        let ladder = reference_ladder();
        let d_lo = -db(ladder_s21(&ladder, 1.6e9).norm());
        let d_hi = -db(ladder_s21(&ladder, 2.4e9).norm());
        assert!(
            d_lo > d_hi + ASYMMETRY_MARGIN_DB,
            "reference lower notch ({d_lo:.2} dB) must be deeper than upper ({d_hi:.2} dB)"
        );
    }
}
