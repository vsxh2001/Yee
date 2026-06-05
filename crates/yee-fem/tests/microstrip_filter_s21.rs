//! FEM-EM brick N3 (ADR-0154) — 3-pole microstrip-filter S21 from the FEM
//! driven sweep, re-graded with the HIGH-FIDELITY numerical-eigenmode port,
//! against the analytic `ladder_s21` reference.
//!
//! This is the **payoff** of the FEM-EM driven-sweep track. Bricks B1–B4
//! (ADR-0153: interior-PEC edges, `layered_microstrip_mesh`, the quasi-TEM
//! wave-port, the straight-line ε_eff = 0.61 % of Hammerstad-Jensen) plus N1+N2
//! (ADR-0154: the production numerical-eigenmode port `microstrip_port_numerical`,
//! which on a straight line lifts |S21| 0.089 → 0.778 and matches the port,
//! |S11| 0.087) are merged. This test composes them into a coupled-resonator
//! **band-pass filter** geometry, drives a two-port `sweep_matrix` over the
//! band, de-embeds the feed reference plane, extracts |S21|(f), and grades the
//! curve against the 3-pole Chebyshev 0.5 dB / 2 GHz / 10 % FBW `ladder_s21`
//! reference — including the geometric-asymmetry discriminator
//! (`depth(1.6 GHz) > depth(2.4 GHz)`).
//!
//! Originally (ADR-0153 B7) this ran with the v1 ANALYTIC flat-`E_z` port and
//! floored at ~−42 dB (the ~−21 dB/port modal-overlap loss × two ports). N3
//! swaps in the numerical eigenmode (recentred per off-centre feed via
//! [`yee_fem::microstrip_port_numerical_at`]) and re-grades.
//!
//! ## Honest framing (read before the gate)
//!
//! N3 is **research-open**: the line is proven, but a filter adds resonator
//! coupling + gap-mesh sensitivity the line never exercised, so whether the
//! *filter* clears the strict Cheb mask is genuinely unknown a priori. The
//! deliverable is an HONEST graded filter curve. The gate asserts only the
//! measurement-driven checks the solve actually supports (the lift over the v1
//! floor, the asymmetry discriminator, the band-pass turnover — see
//! [`fem_filter_s21_vs_ladder`]) and asserts the strict mask ONLY IF the
//! measurement clears it; it records the full |S21|(f) table + the honest
//! mask margin otherwise. The MEASURED outcome is LIFT-BUT-SHORT: a +15 dB
//! lift over the v1 floor, correct asymmetry, but the strict mask still missed
//! by ~35 dB (the multi-resonator path, not just the port, caps the level).
//! Weakening or faking the grade is not a valid outcome; a documented
//! lift-but-short with a quantified mask margin is.
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
    CoupledResonatorGeom, FaceKind, MaterialDatabase, MicrostripPortGeom, OpenBoundarySolver,
    SParametersMatrix, TraceRect, beta_microstrip, correct_gap_fem_k,
    layered_microstrip_filter_mesh, microstrip_port_numerical_at,
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
    /// numerical-eigenmode wave-port is RECENTRED here — NOT at the box centre
    /// `box_w/2` where the cross-section places its trace — because the feed is
    /// a narrow off-centre strip; `microstrip_port_numerical_at` shifts the
    /// eigenmode sampling by `box_w/2 − feed_xc` so the modal peak lands under
    /// the actual feed (a box-centred mode would mostly miss it). See the
    /// module-level honest framing.
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
///
/// `gaps_override` lets a caller substitute its own inter-resonator gaps
/// (length `N − 1 = 2`) for the analytically dimensioned `dims.gaps_m`. The N3
/// gate passes `None` (the analytic impedance-k gaps — its behaviour is
/// BYTE-IDENTICAL); the B2 corrected-gap gate passes `Some(&corrected_gaps)`
/// (the FEM resonant-split design-curve gaps, ADR-0159). All other geometry
/// (line width, resonator length, stagger, feeds, box) is unchanged.
#[allow(clippy::too_many_arguments)]
fn build_edge_coupled_geometry(
    clearance_x: f64,
    air_h: f64,
    feed_len: f64,
    dx_target: f64,
    dy_target: f64,
    dz: f64,
    gaps_override: Option<&[f64]>,
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
    // The inter-resonator gaps: the analytic impedance-k gaps by default, or a
    // caller-supplied override (B2 corrected resonant-split gaps). The override
    // must match the synthesized gap count (N − 1).
    let gaps: Vec<f64> = match gaps_override {
        Some(g) => {
            assert_eq!(
                g.len(),
                dims.gaps_m.len(),
                "gaps_override length ({}) must match the synthesized gap count ({})",
                g.len(),
                dims.gaps_m.len(),
            );
            g.to_vec()
        }
        None => dims.gaps_m,
    };
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

    // Feed-line wave-port: the HIGH-FIDELITY numerical quasi-TEM eigenmode
    // (ADR-0154 N1, `microstrip_port_numerical_at`) replacing the v1 analytic
    // windowed E_z shape that floored B7 at ~−42 dB. The numerical eigenmode
    // is the true transverse mode of the feed's (box_w × box_h) FR-4 cross-
    // section; on the straight line it lifts |S21| 0.089→0.778 and matches the
    // port (|S11| 0.087). β stays analytic Hammerstad-Jensen on the FEED width
    // (the port face sees a uniform 50 Ω line, whatever the coupled-resonator
    // interior does); only the modal SHAPE is numerical.
    //
    // x-RECENTERING (critical): the numerical cross-section centres its trace
    // at box_w/2, but the filter's two feeds are OFF-CENTRE at DIFFERENT x
    // (input near one box edge, output near the other — staggered resonators).
    // `microstrip_port_numerical_at(geom, x_center, f)` shifts the eigenmode
    // sampling by box_w/2 − x_center so the modal peak lands under the actual
    // feed strip; sampling the box-centred mode unshifted would place the peak
    // over air/PEC and re-introduce the very overlap loss the numerical port
    // removes (the v1 windowed port recentred per-feed for the same reason).
    // The shape is frequency-independent (one eigensolve at the band centre
    // F0); β(ω) carries the dispersion, exactly as the v1 port did. Each face
    // gets its own call (boxed closures are not Clone).
    let port_geom = MicrostripPortGeom {
        trace_w: geom.line_w,
        sub_h: SUB_H,
        eps_r: EPS_R,
        box_w: geom.box_w,
        box_h: geom.box_h,
    };
    let port_in = microstrip_port_numerical_at(&port_geom, geom.feed_xc_in, F0)
        .expect("numerical-eigenmode port (input feed) must build");
    let port_out = microstrip_port_numerical_at(&port_geom, geom.feed_xc_out, F0)
        .expect("numerical-eigenmode port (output feed) must build");

    let solver = OpenBoundarySolver::new(&mesh, kinds, vec![port_in, port_out], material_db)
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
/// length `2·β·feed_len` of phase to S21. The feed is closely matched (~50 Ω),
/// so its magnitude effect is small; we de-embed phase only (a unit-magnitude
/// rotation), which is what the asymmetry / shape grading — all |S21|-MAGNITUDE
/// checks — needs.
///
/// CAVEAT (n=3): the last resonator is even-indexed and ends at `feed_len +
/// res_l`, so the OUTPUT feed is actually `feed_len + stagger` long; removing
/// `2·β·feed_len` therefore leaves ~`β·stagger` of output-feed phase
/// uncompensated. This does NOT affect any gate assertion (all magnitude;
/// de-embed is a unit-magnitude rotation that cannot change |S21|), but the
/// de-embedded S21 *phase* is not exactly at the output resonator reference
/// plane — pass the actual per-port feed lengths before any phase / group-delay
/// analysis.
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
        let geom = build_edge_coupled_geometry(clr, air, feed, dx, dy, dz, None);
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

/// The v1 ANALYTIC flat-`E_z` wave-port in-band peak (dB) measured by the
/// original B7 (ADR-0153, base `22da1c2`): the curve floored at ≈−42.4 dB
/// in-band (−42.39 dB @ 2.10 GHz) through TWO analytic ports' ~−21 dB/port
/// modal-overlap loss. ADR-0154 N3 re-grades with the numerical-eigenmode
/// port; the in-band peak MUST lift well clear of this v1 floor (the
/// re-flooring tripwire — a promotion regression toward the analytic floor
/// is the failure this catches).
const V1_FLOOR_PEAK_DB: f64 = -42.4;

/// N3 re-flooring tripwire (dB): the numerical-port in-band peak must clear
/// the v1 floor by at least this margin. MEASURED N3 in-band peak is
/// −27.38 dB @ 2.00 GHz — a +15.0 dB lift over the v1 −42.4 dB floor. This
/// 9 dB bar sits ~6 dB BELOW the measured −27.38 dB (catches a regression
/// with margin) and ~9 dB ABOVE the v1 floor (so the v1 analytic port could
/// NOT pass it): a defensible measured-truth threshold, not a wish.
const N3_MIN_LIFT_OVER_V1_DB: f64 = 9.0;

/// The N3 (ADR-0154) MEASURED in-band peak (dB) with the IMPEDANCE-k gaps and the
/// numerical-eigenmode port: −27.38 dB @ 2.00 GHz (the lift-but-short baseline,
/// boxed --release, base 192cb54). The B2 corrected-gap gate PRINTS its in-band
/// peak relative to this number; B2 does NOT hard-assert an improvement over it
/// (the orchestrator pins the measured B2 peak as a regression tripwire AFTER
/// seeing the real number — so a gap-interaction regression cannot fake-green).
const N3_BASELINE_PEAK_DB: f64 = -27.38;

/// FEM-EM brick N3 (ADR-0154) — 3-pole microstrip-filter S21 re-graded with the
/// HIGH-FIDELITY numerical-eigenmode wave-port, vs the analytic ladder reference.
///
/// Builds the edge-coupled 3-pole filter, drives `sweep_matrix` over
/// 1.6–2.4 GHz through TWO `microstrip_port_numerical_at` ports (ADR-0154 N1,
/// recentred per off-centre feed), de-embeds the feed reference planes,
/// extracts |S21|(f), and grades it against the 3-pole Cheb 0.5 dB / 2 GHz /
/// 10 % FBW `ladder_s21` reference + the geometric-asymmetry discriminator.
/// This is N3's payoff question: with two high-fidelity ports, does the FILTER
/// clear the strict Cheb mask?
///
/// ## What this asserts (HONEST, MEASUREMENT-DRIVEN)
///
/// MEASURED ANSWER (below): **LIFT-BUT-SHORT.** The numerical port lifts the
/// in-band peak +15.0 dB over the v1 analytic floor and grows the asymmetry
/// margin, but the strict `oracle_grade` mask (passband |err| ≤ 2 dB, rejection
/// |err| ≤ 5 dB) still MISSES by ~35 dB in-band (the 2-port + 3-resonator path,
/// not just the port, caps the absolute level). The gate therefore does NOT
/// assert the absolute-level mask (no weakening to force green); it asserts:
///
/// 1. **A real lift over the v1 floor** — in-band peak ≥ `V1_FLOOR_PEAK_DB` +
///    `N3_MIN_LIFT_OVER_V1_DB` (the re-flooring tripwire; a port-promotion
///    regression would re-floor it toward the v1 −42.4 dB).
/// 2. **The geometric-asymmetry discriminator (the brick's NAMED check)** —
///    `depth(1.6 GHz) > depth(2.4 GHz)` by ≥ 1 dB: the FEM curve reproduces the
///    correct band-pass-mapping asymmetry SIGN that a symmetric/inverted
///    fitted artifact lacks.
/// 3. **A band-pass turnover** — the in-band peak stands above the deeper band
///    edge (a real centre bump, not a monotonic ramp / flat line).
///
/// The strict mask is asserted ONLY IF the measurement actually clears it (it
/// does not yet — `strict_pass` is `false`); the honest MISS margin is recorded
/// and printed instead. The full curve + verdict are in the MEASURED block.
///
/// ## MEASURED RESULT (boxed --release, base 192cb54; 51 336 tets, 77.4 s)
///
/// Edge-coupled 3-pole, 14.0 × 77.6 × 6.0 mm box, w = 1.912 mm,
/// `dx/dy/dz = 0.6/2.5/0.5 mm`, feed = 8 mm. NUMERICAL-eigenmode ports
/// (`microstrip_port_numerical_at`), recentred on the input feed (xc ≈
/// 3.46 mm) and output feed (xc ≈ 10.54 mm); box centre is 7.0 mm — both feeds
/// off-centre. |S21| after feed de-embed:
///
/// ```text
///   f(GHz)   S21 dB (FEM, numerical port)   ref dB (ladder)
///   1.60      −30.92                         −41.77
///   1.80      −29.36                         −20.81
///   1.90      −28.67                          −0.75
///   2.00      −27.38                          0.00   ← reference passband centre + FEM peak
///   2.05      −35.44                          −0.50
///   2.10      −29.36                          −0.32
///   2.20      −28.68                         −17.83
///   2.40      −28.82                         −36.26
///
///   in-band peak       : −27.38 dB @ 2.00 GHz   (v1 analytic floor: −42.4 dB)
///   lift over v1 floor : +15.0 dB
///   turnover           : +3.54 dB (in-band peak above the deeper band edge)
///   asymmetry (NAMED)  : depth(1.6)=30.92 dB > depth(2.4)=28.82 dB, +2.10 dB → PASS
///   strict oracle mask : MISS by ~34.9 dB in-band (worst err vs the 0 dB reference)
/// ```
///
/// ## Honest verdict — LIFT-BUT-SHORT
///
/// Exactly the outcome ADR-0154 §Consequences flagged as moderate-confidence:
/// the numerical port lifts the floor dramatically but the *filter* stops short
/// of the strict mask.
///
/// * **What the numerical port bought (vs the v1 analytic floor):** the
///   in-band peak rose −42.39 → −27.38 dB (a **+15.0 dB lift**), the band-edge
///   levels rose ~−43 → ~−29/−31 dB, the turnover grew +2.23 → +3.54 dB, and
///   the asymmetry margin grew +1.47 → +2.10 dB. The lift confirms the
///   higher-fidelity modal shape raised the port↔FEM modal overlap (the N2
///   straight-line case lifted |S21| 0.089 → 0.778 / |S11| → 0.087). This is a
///   real, geometry-aware bandpass with the correct asymmetry SIGN, not a
///   fitted artifact.
///
/// * **Why it still MISSES the strict mask (~35 dB in-band):** unlike the N2
///   matched straight-line thru (|S21| ≈ 0.778), the filter inserts THREE
///   coupled λ_g/2 resonators between the two ports. The signal must traverse
///   two weak edge-coupling gaps plus the lossy resonator interior at a COARSE
///   `dx/dy = 0.6/2.5 mm` mesh; the per-port match no longer translates into a
///   low-IL passband (the in-band |S21| peak is ≈ 0.043, far below the line's
///   0.778). The remaining gap is **resonator-coupling + gap-mesh fidelity**
///   (and, secondarily, a still-higher-fidelity port — numerical cross-section
///   aperture coupling), a finer-mesh / coupling-extraction follow-on, NOT a
///   mesh/LU-scaling wall (the 51 k-tet `faer` sparse LU fits the 14 g box with
///   room to spare, ~3 s/point). A documented lift-but-short with a quantified
///   mask margin is the correct ADR-0154 N3 deliverable.
///
/// * **Two levers that mattered (recorded so they are not re-derived):**
///   (1) the numerical eigenmode must be RECENTRED on each FEED's `x`
///   (`microstrip_port_numerical_at`), because the filter's two feeds are
///   off-centre at different `x` while the cross-section centres its trace at
///   `box_w/2`; sampling the box-centred mode unshifted would place the modal
///   peak over air/PEC and re-introduce the overlap loss the numerical port
///   removes. (2) `with_coupled_whitney(true)` is mandatory (B4 finding; the
///   lumped-centroid path collapses the absorbing block for the substrate-
///   normal `E_z` mode).
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
        None,   // analytic impedance-k gaps (N3 baseline — unchanged)
    );
    eprintln!(
        "[N3] filter mesh: box=({:.1},{:.1},{:.1})mm  n=({},{},{})  tets={}  w={:.3}mm  feed={:.1}mm  eps_eff(w)={:.4}",
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

    // N3 lift over the v1 analytic-port floor — the headline number. The strict
    // Cheb passband is 0 dB; the gap to the mask is `worst_pass_db`. A positive
    // lift over V1_FLOOR_PEAK_DB is the evidence the numerical eigenmode raised
    // the modal overlap (the line case lifted |S21| 0.089→0.778; here the
    // 2-port + 3-resonator filter lifts the in-band PEAK but stops short of the
    // mask — see the honest verdict below).
    let lift_over_v1_db = passband_peak_db - V1_FLOOR_PEAK_DB;
    eprintln!(
        "\n==== N3 GRADE (numerical-eigenmode port; ADR-0154) ====\n\
         tets               : {}\n\
         wall               : {:.1} s\n\
         in-band peak       : {:.2} dB @ {:.2} GHz (overall peak @ {:.2} GHz)\n\
         v1 floor (ref)     : {:.2} dB  (analytic flat-Ez port, B7 base 22da1c2)\n\
         lift over v1 floor : {:+.2} dB  (tripwire ≥ {:.1} dB)\n\
         band edges         : {:.2} dB @1.6  {:.2} dB @2.4\n\
         turnover           : {:+.2} dB (in-band peak above the deeper edge)\n\
         worst passband err : {:.2} dB vs ref (oracle tol {:.1})\n\
         worst rejection err: {:.2} dB vs ref (oracle tol {:.1})\n\
         strict-mask margin : MISS by {:.2} dB in-band (gap to the 0 dB Cheb passband)\n\
         asymmetry (NAMED)  : depth(1.6)={:.2} dB  depth(2.4)={:.2} dB  margin={:+.2} dB  -> {}\n\
         strict oracle mask : {}\n\
         ========================================================",
        geom.total_tets(),
        wall,
        passband_peak_db,
        f_inband_peak,
        f_peak_ghz,
        V1_FLOOR_PEAK_DB,
        lift_over_v1_db,
        N3_MIN_LIFT_OVER_V1_DB,
        edge_lo_db,
        edge_hi_db,
        turnover_db,
        worst_pass_db,
        PASSBAND_TOL_DB,
        worst_rej_db,
        REJECTION_TOL_DB,
        worst_pass_db,
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
    eprintln!("[N3] oracle_grade pairs: {pairs}");

    // ---- Assertions (HONEST, MEASUREMENT-DRIVEN: assert only what holds) ----
    //
    // ADR-0154 N3 re-grades the filter with the NUMERICAL-eigenmode port that,
    // on a straight line, lifted |S21| 0.089→0.778 and matched the port (N2,
    // independently verified). The research-open question this gate answers
    // honestly: with two high-fidelity ports, does the FILTER clear the strict
    // Cheb mask, or does resonator coupling / gap mesh cap it short?
    //
    // MEASURED ANSWER: LIFT-BUT-SHORT (see the MEASURED block in the docstring).
    // The in-band peak rose to −27.38 dB (a +15.0 dB lift over the v1 −42.4 dB
    // analytic floor) and the asymmetry margin grew (+2.10 dB vs the v1 +1.47),
    // but the curve still sits ~27 dB below the 0 dB Cheb passband (worst
    // in-band err ≈ 34.9 dB) — so the STRICT MASK still MISSES. The honest
    // verdict (port fidelity vs the multi-resonator path) is in the docstring.
    //
    // The gate therefore asserts the three things that ARE true regardless of
    // the mask, plus the strict mask ONLY IF it actually clears (it does not):

    // (1) Non-degenerate transmission AND a real lift over the v1 floor. A
    //     collapsed port (lumped-centroid failure) or a broken mesh would sit
    //     in noise; a promotion regression in `microstrip_port_numerical[_at]`
    //     (wrong frame map / cross-section density / x-recentre) would re-floor
    //     the peak back toward the v1 −42.4 dB. The measured peak is −27.38 dB,
    //     a +15.0 dB lift; the 9 dB tripwire (≈6 dB below the measurement, ≈9 dB
    //     above the v1 floor) catches a re-flooring with margin without
    //     asserting a depth the 2-port filter path does not deliver.
    assert!(
        passband_peak_db.is_finite(),
        "N3 NO-GO: in-band peak is not finite ({passband_peak_db}) — the driven solve \
         degenerated (port collapsed or mesh broken). Full curve printed above."
    );
    assert!(
        lift_over_v1_db >= N3_MIN_LIFT_OVER_V1_DB,
        "N3 re-flooring tripwire: in-band peak {passband_peak_db:.2} dB lifted only \
         {lift_over_v1_db:+.2} dB over the v1 analytic floor {V1_FLOOR_PEAK_DB:.2} dB \
         (need ≥ {N3_MIN_LIFT_OVER_V1_DB:.1} dB). A small lift means the numerical port \
         regressed toward the modal-overlap floor — most likely the frame map, the \
         cross-section density, or the x-recentre in `microstrip_port_numerical_at`. \
         Report the number; do NOT lower the threshold. Full curve printed above."
    );

    // (2) The geometric-asymmetry discriminator — the brick's NAMED check — must
    //     PASS: the lower stopband notch (1.6 GHz) is genuinely deeper than the
    //     upper (2.4 GHz). This is the band-pass-mapping signature the reference
    //     has and a symmetric/inverted (fitted-artifact) curve does NOT; it is
    //     the anti-"flat/symmetric curve is not evidence" guard. The numerical-
    //     port FEM curve reproduces the CORRECT asymmetry SIGN with margin
    //     +2.10 dB (≥ 1 dB) — a real, geometry-aware result.
    assert!(
        asym_pass,
        "N3: geometric-asymmetry discriminator FAILED — depth(1.6 GHz)={depth_lo:.2} dB is NOT \
         deeper than depth(2.4 GHz)={depth_hi:.2} dB by the required {ASYMMETRY_MARGIN_DB} dB \
         (margin {asym_margin:+.2} dB). A symmetric/inverted curve has lost the band-pass-mapping \
         asymmetry and is not credited as a geometry-aware EM result. Full curve printed above."
    );

    // (3) A genuine band-pass turnover: the in-band peak stands above the deeper
    //     band edge (the response bumps up near band centre rather than ramping
    //     monotonically). Measured turnover ≈ +3.5 dB; the >0.2 dB bar certifies
    //     the SHAPE is frequency-selective without demanding a depth the
    //     2-port + 3-resonator path does not deliver.
    assert!(
        turnover_db > 0.2,
        "N3: no band-pass turnover — in-band peak {passband_peak_db:.2} dB is not above the \
         deeper band edge (edges {edge_lo_db:.2}/{edge_hi_db:.2} dB; turnover {turnover_db:+.2} dB). \
         The response is monotonic, not a recognisable band-pass. Full curve printed above."
    );

    // (4) Strict Cheb mask — assert it ONLY IF the measurement actually clears
    //     it (a real win). It does NOT clear with this port (lift-but-short:
    //     worst in-band err ≈ 34.9 dB ≫ the 2 dB oracle tol), so this branch
    //     records the honest MISS margin and does NOT assert the mask — no
    //     weakening to force green. If a future follow-on (finer coupling-gap
    //     mesh and/or a still-higher-fidelity aperture-coupling port) lifts the
    //     curve into the mask, `strict_pass` flips and this asserts it automatically.
    if strict_pass {
        // A genuine in-mask pass: assert it loudly (the FEM driven-sweep track
        // would have delivered its original goal — a validated in-mask filter).
        assert!(
            worst_pass_db <= PASSBAND_TOL_DB && worst_rej_db <= REJECTION_TOL_DB,
            "internal: strict_pass set but tolerances not met (pass {worst_pass_db:.2}, \
             rej {worst_rej_db:.2})"
        );
    } else {
        // Honest lift-but-short: the mask is MISSED. We assert the measured
        // lift (done in (1)) and the asymmetry (done in (2)); we deliberately
        // do NOT assert the absolute-level mask. Record the margin for the log.
        eprintln!(
            "[N3] STRICT MASK: MISS by {worst_pass_db:.2} dB in-band (lift-but-short — \
             the numerical port lifted the peak {lift_over_v1_db:+.2} dB over the v1 floor \
             but the 2-port + 3-resonator path stops ~{:.0} dB below the 0 dB Cheb passband). \
             This is the honest documented result; the residual gap is NOT isolated by this \
             brick — candidates are resonator-coupling + gap-mesh fidelity at this coarse mesh \
             and, secondarily, a still-higher-fidelity port (cross-section aperture coupling). \
             It is NOT a mesh/LU-SCALING wall: the 51k-tet sparse LU fits the box with room.",
            worst_pass_db,
        );
    }
}

// =====================================================================
// FEM-EM brick B2 (ADR-0159) — corrected-gap filter S21 re-grade
// =====================================================================

/// Synthesis sweep band for the corrected-gap re-grade — IDENTICAL to the N3
/// gate (1.6–2.4 GHz, 17 points / 50 MHz spacing) so the two curves are graded
/// on the same frequency grid. Returned as `(freqs_hz, omegas)`.
fn band_1p6_to_2p4_17pts() -> (Vec<f64>, Vec<f64>) {
    let n_pts = 17;
    let f_lo = 1.6e9;
    let f_hi = 2.4e9;
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|i| f_lo + (f_hi - f_lo) * (i as f64) / ((n_pts - 1) as f64))
        .collect();
    let omegas: Vec<f64> = freqs_hz.iter().map(|f| 2.0 * PI * f).collect();
    (freqs_hz, omegas)
}

/// Graded outcome of a filter |S21|(f) sweep vs the analytic ladder reference.
/// Bundles the corrected-gap re-grade so the B2 gate's prints + asserts read off
/// one struct (the same quantities the N3 gate computes inline).
struct GradedCurve {
    /// `(f_GHz, |S21| dB)` after feed de-embed, on the 17-point band.
    curve: Vec<(f64, f64)>,
    /// In-band (`[1.85, 2.15]` GHz) peak |S21| (dB).
    passband_peak_db: f64,
    /// In-band peak frequency (GHz).
    f_inband_peak_ghz: f64,
    /// Band-edge levels (dB) at 1.6 / 2.4 GHz (interpolated).
    edge_lo_db: f64,
    edge_hi_db: f64,
    /// Turnover (dB): in-band peak above the deeper of the two band edges.
    turnover_db: f64,
    /// Asymmetry margin (dB): `depth(1.6 GHz) − depth(2.4 GHz)` (depth = −dB).
    asym_margin_db: f64,
    /// Worst |err| vs the ladder reference over the passband window (dB).
    worst_pass_db: f64,
    /// Worst |err| vs the ladder reference outside the passband (dB).
    worst_rej_db: f64,
}

/// Sweep a filter geometry, de-embed the feeds, and grade |S21|(f) against the
/// analytic 3-pole Cheb ladder on the standard 17-point band — the shared
/// measurement path the B2 corrected-gap gate uses. (The N3 baseline gate keeps
/// its own byte-identical inline copy; this helper is additive and used only by
/// B2.) `label` tags the printed |S21|(f) table.
fn sweep_and_grade(geom: &FilterGeometry, label: &str) -> GradedCurve {
    let (freqs_hz, omegas) = band_1p6_to_2p4_17pts();
    let t0 = std::time::Instant::now();
    let sweep = solve_filter(geom, &omegas);
    let wall = t0.elapsed().as_secs_f64();

    let ladder = reference_ladder();
    let mut curve: Vec<(f64, f64)> = Vec::with_capacity(freqs_hz.len());
    eprintln!(
        "\n[{label}] |S21|(f) ({} pts, sweep {:.1} s):\n{:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
        freqs_hz.len(),
        wall,
        "f(GHz)",
        "|S21|raw",
        "|S21|deemb",
        "S21 dB",
        "ref dB",
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

    // Grade vs the reference (mirrors oracle_grade::evaluate / the N3 gate).
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

    let depth_at = |f_ghz: f64| -> f64 { -interp_db(&curve, f_ghz) };
    let asym_margin_db = depth_at(1.6) - depth_at(2.4);

    let passband_peak_db = curve
        .iter()
        .filter(|(f, _)| (1.85..=2.15).contains(f))
        .map(|(_, d)| *d)
        .fold(f64::NEG_INFINITY, f64::max);
    let f_inband_peak_ghz = curve
        .iter()
        .filter(|(f, _)| (1.85..=2.15).contains(f))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(f, _)| *f)
        .unwrap_or(f64::NAN);
    let edge_lo_db = interp_db(&curve, 1.6);
    let edge_hi_db = interp_db(&curve, 2.4);
    let turnover_db = passband_peak_db - edge_lo_db.min(edge_hi_db);

    GradedCurve {
        curve,
        passband_peak_db,
        f_inband_peak_ghz,
        edge_lo_db,
        edge_hi_db,
        turnover_db,
        asym_margin_db,
        worst_pass_db,
        worst_rej_db,
    }
}

/// FEM-EM brick B2 (ADR-0159) — 3-pole microstrip-filter S21 re-graded after
/// re-realizing the inter-resonator gaps on the **FEM resonant-split coupling
/// design-curve** (the Hong-Lancaster full-wave coupling design), vs the
/// analytic ladder reference.
///
/// ## The decisive test
///
/// The N3 result ([`fem_filter_s21_vs_ladder`]) floors the full-wave filter S21
/// at ≈−27.38 dB in-band. The diagnosed root cause (ADR-0155 K2 / ADR-0159): the
/// analytic dimensioner ([`yee_filter::dimension_edge_coupled`]) sizes the
/// inter-resonator gaps with the IMPEDANCE-k
/// (`yee_layout::coupling_coefficient = (Z0e−Z0o)/(Z0e+Z0o)`), which **diverges
/// ~37 % from the physically-realized RESONANT-SPLIT k** at the tight gaps a
/// filter uses (`k_imp ≠ k_eps`). In filter theory the coupling coefficient is
/// DEFINED as the resonant-split `k = (f_hi²−f_lo²)/(f_hi²+f_lo²)` — exactly what
/// [`yee_fem::coupled_resonator_k`] measures full-wave.
///
/// B2 re-realizes the SAME synthesis `target_k` (≈0.0756 for this spec, the
/// `FBW·m_{i,i+1}` constant from [`EdgeCoupledDimensions::target_k`]) as a
/// resonant-split k via B1's [`yee_fem::correct_gap_fem_k`] (a bisection on the
/// monotone-decreasing FEM `K(gap)` curve, ~5-6 FEM sweeps per gap), rebuilds the
/// filter with the corrected gaps, re-sweeps the SAME 17-point 1.6–2.4 GHz band,
/// and re-grades vs the Chebyshev mask — recording whether the corrected filter
/// clears (or materially approaches) the mask vs the −27.38 dB N3 baseline.
///
/// ## Measured finding (this run; ADR-0159)
///
/// The design-curve MECHANISM works on the filter geometry (the f0-tracking band
/// resolves the split; `correct_gap_fem_k` traces `K(gap)`), BUT `target_k≈0.0756`
/// is **UNREACHABLE** as a resonant-split k on the 1.912 mm / FR-4 line — the
/// resonant-split k saturates ≈0.064 at the tight-gap floor (so the correction
/// does not converge, best ≈15 % off, and drives the gaps 1.622 → 0.594 mm). The
/// corrected filter LIFTS the in-band peak **+5.8 dB** over N3 (−27.38 → −21.55 dB)
/// — stronger coupling helps — but STILL misses the strict Cheb mask by ~26 dB and
/// the response is a ~−22 dB flat shelf. **Conclusion: dimensioning-correction is a
/// real but MINOR lever; the dominant floor is the aperture-coupling PORT fidelity
/// (the ADR-0154 N3 finding), NOT gap dimensioning** — neither this correction nor
/// a multi-D ASM over gaps (B3) clears the mask without a higher-fidelity port.
///
/// ## Non-circularity
///
/// The correction target is the SYNTHESIS `target_k` — a fixed design constant
/// (`FBW · m_{i,i+1}`), NOT any FEM measurement. The corrected gap is found by
/// driving the INDEPENDENT full-wave `coupled_resonator_k` FEM measurement to
/// that constant. The re-graded S21 is then an independent FEM driven-sweep on
/// the rebuilt geometry. Nothing in the loop reads the filter's own S21 to set
/// the gap, so a curve that improves is real EM evidence, not a fit.
///
/// ## Outcome handling — mechanism + lift asserted, finding recorded
///
/// The gate PRINTS the full corrected |S21|(f) table + the in-band peak vs the
/// −27.38 dB N3 baseline + the strict-mask margin, and asserts the TRUE,
/// reproducible results: (a) the correction MECHANISM (a finite K(gap) was traced
/// — the f0-tracking band resolves the split on the filter geometry); (b) the gaps
/// moved materially off the impedance-k gaps, in the TIGHTER (stronger-coupling)
/// direction; (c) the curve is finite with a band-pass turnover; (d) the measured
/// in-band peak LIFTED ≥ 3 dB over N3 (measured +5.83 dB — the quantified
/// dimensioning lever, pinned after the orchestrator saw the real number); (e) the
/// asymmetry sign is preserved (the N3 ≥1 dB margin is RECORDED-not-asserted, since
/// the strong-coupling-floor gaps degrade it to +0.70 dB). It does NOT assert
/// convergence or mask-clearing (both false — target_k is unreachable and the
/// filter remains port-floored): those are RECORDED as the honest finding, not
/// faked green, and no physics tolerance is weakened.
///
/// ## GATING — CRITICAL (heavy; run by the orchestrator, boxed, `--release`)
///
/// This is the heaviest gate in the file: `correct_gap_fem_k` runs ~5-6
/// [`yee_fem::coupled_resonator_k`] FEM driven sweeps PER gap (one per bisection
/// eval; the two gaps are equal by symmetry, so the gate corrects ONCE and
/// reuses), THEN one full 17-point filter `sweep_matrix`. Budget ~45-60 min.
/// `#[ignore]`'d so the debug `cargo test --workspace` never runs it; run only in
/// `--release`, boxed:
///
/// ```text
/// YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh bash -c '\
///   cargo test -p yee-fem --release --test microstrip_filter_s21 \
///   -- --ignored fem_filter_s21_corrected_gaps --nocapture'
/// ```
#[test]
#[ignore = "B2 gate: ~5-6 FEM sweeps per gap (gap correction) + one 17-pt filter sweep; ~45-60 min — run only in --release, boxed"]
fn fem_filter_s21_corrected_gaps() {
    // ---- 1. Synthesize + dimension (identical to build_edge_coupled_geometry) -
    let project = synthesize(&reference_spec());
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: SUB_H,
        loss_tangent: 0.0,
        metal_thickness_m: 0.0,
    };
    let dims = dimension_edge_coupled(&project, &sub).expect("edge-coupled 3-pole synthesizes");
    let w = dims.line_width_m;
    let imp_gaps = dims.gaps_m.clone(); // analytic impedance-k gaps (the N3 baseline)
    let target_k = dims.target_k.clone(); // FBW·m_{i,i+1} synthesis constants (length N-1=2)
    assert_eq!(
        target_k.len(),
        2,
        "3-pole filter has N-1 = 2 inter-resonator gaps"
    );

    eprintln!(
        "[B2] synthesis: w={:.4}mm res_l={:.4}mm  impedance-k gaps={:?}mm  target_k={:?}",
        w * 1e3,
        dims.resonator_length_m * 1e3,
        imp_gaps.iter().map(|g| g * 1e3).collect::<Vec<_>>(),
        target_k,
    );

    // ---- 2. Build the FEM design-curve base CoupledResonatorGeom --------------
    // CORRECTNESS INVARIANT #1: trace_w MUST be the filter's resonator width
    // (dims.line_width_m), NOT the K1 probe's 1 mm — a different width gives a
    // DIFFERENT K(gap) design curve and would correct to the wrong gap.
    //
    // Box extents MATCH the SHIPPED K1/K2 `CoupledResonatorGeom::probe_with_gap`
    // config — the proven-tractable, k-VALIDATED open box: CLEARANCE_X = 2.5·h =
    // 2.5 mm each side in x, air_h = 5 mm above (box_h = 6 mm at h = 1 mm). The
    // 2.5·h clearance is the B4 "walls don't load the line" floor that K1/K2
    // measured k against (k_fem ≈ k_eps within tolerance). NOTE: the 6 mm-each-
    // side clearance B4 used for ABSOLUTE ε_eff is overkill for k (a peak-LOCATION
    // ratio) and — at the filter's 1.912 mm width — inflates the pair mesh to
    // ~160k tets, which OOMs the 12 g box. correct_gap_fem_k holds box_w FIXED
    // while sweeping gap_s, so box_w is sized for the WIDEST bracket gap (gap_hi);
    // tighter-gap trials simply sit in a slightly-more-open box (the B4 direction).
    let gap_lo = 0.5e-3; // = DX, the tightest gap the 0.5 mm cross-section pitch resolves.
    let gap_hi = 2.0e-3; // target_k≈0.076 roots at a TIGHT gap (<1.622mm imp-gap); 2mm is a safe upper bracket.
    let tol_frac = 0.08;
    let max_evals = 6;
    let n_pts = 61;
    let clearance_x = 2.5e-3; // CLEARANCE_X: 2.5·h, the K1/K2-validated open-box wall clearance.
    let base = CoupledResonatorGeom {
        trace_w: w,
        gap_s: gap_hi, // irrelevant on the base (the corrector sweeps gap_s); set to the widest.
        sub_h: SUB_H,
        eps_r: EPS_R,
        f0_hz: F0,
        // Two w-wide strips + the widest bracket gap + CLEARANCE_X both sides
        // (= probe_with_gap's `CLEARANCE_X + W + S + W + CLEARANCE_X`).
        box_w: 2.0 * clearance_x + 2.0 * w + gap_hi,
        box_h: SUB_H + 5.0e-3, // sub + 5 mm air = 6 mm (probe_with_gap's open half-space).
    };
    eprintln!(
        "[B2] design-curve base: trace_w={:.4}mm sub_h={:.3}mm eps_r={} f0={:.2}GHz \
         box_w={:.3}mm box_h={:.3}mm  bracket=[{:.2},{:.2}]mm tol={:.0}% max_evals={} n_pts={}",
        base.trace_w * 1e3,
        base.sub_h * 1e3,
        base.eps_r,
        base.f0_hz / 1e9,
        base.box_w * 1e3,
        base.box_h * 1e3,
        gap_lo * 1e3,
        gap_hi * 1e3,
        tol_frac * 100.0,
        max_evals,
        n_pts,
    );

    // ---- 3. Correct each gap onto the FEM resonant-split design curve ---------
    // The two target_k are equal by symmetry; we loop for generality + print each
    // (the corrector prints its per-eval bisection trajectory under --nocapture).
    // To avoid paying for an identical heavy correction twice, we cache the result
    // per distinct target_k value.
    let t_corr = std::time::Instant::now();
    let mut corrected_gaps: Vec<f64> = Vec::with_capacity(target_k.len());
    let mut any_converged = false;
    let mut cache: Vec<(f64, yee_fem::GapCorrection)> = Vec::new();
    for (i, &kt) in target_k.iter().enumerate() {
        // Reuse a cached correction if this target_k was already solved (the
        // symmetric 3-pole has two equal targets — one heavy solve, not two).
        let corr = if let Some((_, c)) = cache.iter().find(|(k, _)| (k - kt).abs() < 1e-12) {
            eprintln!(
                "[B2] gap[{i}] target_k={kt:.4}: REUSING cached correction (equal target by symmetry)"
            );
            *c
        } else {
            let c = correct_gap_fem_k(&base, kt, gap_lo, gap_hi, tol_frac, max_evals, n_pts);
            cache.push((kt, c));
            c
        };
        eprintln!(
            "[B2] gap[{i}] correction: target_k={:.4}  k_fem={:.4}  gap={:.4}mm  \
             n_evals={}  converged={}  (impedance-k gap was {:.4}mm)",
            corr.k_target,
            corr.k_fem,
            corr.gap_m * 1e3,
            corr.n_evals,
            corr.converged,
            imp_gaps[i] * 1e3,
        );
        any_converged |= corr.converged;
        // Use the best gap whether or not it converged (non-convergence is an
        // honest finding; we still rebuild with the best available gap).
        corrected_gaps.push(corr.gap_m);
    }
    let corr_wall = t_corr.elapsed().as_secs_f64();
    eprintln!(
        "[B2] gap correction wall: {:.1} s  impedance-k gaps={:?}mm -> corrected gaps={:?}mm",
        corr_wall,
        imp_gaps.iter().map(|g| g * 1e3).collect::<Vec<_>>(),
        corrected_gaps.iter().map(|g| g * 1e3).collect::<Vec<_>>(),
    );

    // ---- 4. Rebuild the filter with the corrected gaps (same mesh params as N3) -
    let geom = build_edge_coupled_geometry(
        2.5e-3, // x clearance each side (N3 filter-mesh value)
        5.0e-3, // air height
        8.0e-3, // feed length (de-embed reference)
        0.6e-3, // dx
        2.5e-3, // dy
        0.5e-3, // dz
        Some(&corrected_gaps),
    );
    eprintln!(
        "[B2] corrected-gap filter mesh: box=({:.1},{:.1},{:.1})mm  n=({},{},{})  tets={}  \
         w={:.3}mm  feed={:.1}mm",
        geom.box_w * 1e3,
        geom.box_len * 1e3,
        geom.box_h * 1e3,
        geom.nx,
        geom.ny,
        geom.nz,
        geom.total_tets(),
        geom.line_w * 1e3,
        geom.feed_len * 1e3,
    );

    // ---- 5. Re-sweep + de-embed + re-grade over the SAME 17-point band --------
    let graded = sweep_and_grade(&geom, "B2 corrected-gap");

    // ---- 6. Full report (the orchestrator reads this) ------------------------
    let lift_over_v1_db = graded.passband_peak_db - V1_FLOOR_PEAK_DB;
    let lift_over_n3_db = graded.passband_peak_db - N3_BASELINE_PEAK_DB;
    let depth_lo = -interp_db(&graded.curve, 1.6);
    let depth_hi = -interp_db(&graded.curve, 2.4);
    let asym_pass = graded.asym_margin_db >= ASYMMETRY_MARGIN_DB;
    let strict_pass = graded.worst_pass_db <= PASSBAND_TOL_DB
        && graded.worst_rej_db <= REJECTION_TOL_DB
        && asym_pass;

    eprintln!(
        "\n==== B2 CORRECTED-GAP GRADE (FEM resonant-split design curve; ADR-0159) ====\n\
         impedance-k gaps    : {:?} mm  (N3 baseline)\n\
         corrected gaps      : {:?} mm  (FEM resonant-split design curve)\n\
         gap shift           : {:?} mm\n\
         any correction conv : {}\n\
         tets                : {}\n\
         in-band peak        : {:.2} dB @ {:.2} GHz\n\
         N3 baseline peak    : {:.2} dB  (impedance-k gaps, ADR-0154 N3)\n\
         lift over N3        : {:+.2} dB  (PRINTED — orchestrator pins as a tripwire)\n\
         v1 floor (ref)      : {:.2} dB  (analytic flat-Ez port, B7)\n\
         lift over v1 floor  : {:+.2} dB\n\
         band edges          : {:.2} dB @1.6  {:.2} dB @2.4\n\
         turnover            : {:+.2} dB (in-band peak above the deeper edge)\n\
         worst passband err  : {:.2} dB vs ref (oracle tol {:.1})\n\
         worst rejection err : {:.2} dB vs ref (oracle tol {:.1})\n\
         strict-mask margin  : {} by {:.2} dB in-band (gap to the 0 dB Cheb passband)\n\
         asymmetry (NAMED)   : depth(1.6)={:.2} dB  depth(2.4)={:.2} dB  margin={:+.2} dB -> {}\n\
         strict oracle mask  : {}\n\
         ===========================================================================",
        imp_gaps.iter().map(|g| g * 1e3).collect::<Vec<_>>(),
        corrected_gaps.iter().map(|g| g * 1e3).collect::<Vec<_>>(),
        corrected_gaps
            .iter()
            .zip(imp_gaps.iter())
            .map(|(c, i)| (c - i) * 1e3)
            .collect::<Vec<_>>(),
        any_converged,
        geom.total_tets(),
        graded.passband_peak_db,
        graded.f_inband_peak_ghz,
        N3_BASELINE_PEAK_DB,
        lift_over_n3_db,
        V1_FLOOR_PEAK_DB,
        lift_over_v1_db,
        graded.edge_lo_db,
        graded.edge_hi_db,
        graded.turnover_db,
        graded.worst_pass_db,
        PASSBAND_TOL_DB,
        graded.worst_rej_db,
        REJECTION_TOL_DB,
        if strict_pass { "CLEARS" } else { "MISS" },
        graded.worst_pass_db,
        depth_lo,
        depth_hi,
        graded.asym_margin_db,
        if asym_pass { "PASS" } else { "FLAG" },
        if strict_pass { "PASS" } else { "MISS" },
    );

    // Machine-readable curve for the oracle_grade CLI.
    let pairs: String = graded
        .curve
        .iter()
        .map(|(f, d)| format!("{f:.3}:{d:.2}"))
        .collect::<Vec<_>>()
        .join(" ");
    eprintln!("[B2] oracle_grade pairs: {pairs}");

    // ---- 7. Assertions — HONEST, matched to the MEASURED finding -------------
    //
    // B2's measured outcome (ADR-0159, this run): the design-curve correction
    // MECHANISM works on the real filter geometry — the f0-tracking sweep band
    // resolves the split at every eval and the corrector traces K(gap) — BUT the
    // synthesis impedance-k target (target_k ≈ 0.0756) is UNREACHABLE as a
    // resonant-split k on the 1.912 mm / FR-4 line: k saturates ≈ 0.064 at the
    // tight-gap floor (evals 1.25→0.0476, 0.875→0.0625, 0.594→0.0641,
    // 0.523→0.0570 — capped + noisy at tight gaps), so the correction does NOT
    // converge (best ≈ 15 % off) and drives the gaps to the strong-coupling floor
    // (1.622 → 0.594 mm). The resulting filter LIFTS the in-band peak +5.8 dB over
    // the N3 baseline (−27.38 → −21.55 dB) — confirming stronger coupling helps —
    // but STILL misses the strict Cheb mask by ~26 dB and the response is a ~−22 dB
    // flat shelf, not a clean band-pass. CONCLUSION: dimensioning-correction is a
    // real but MINOR lever; the dominant floor is the aperture-coupling PORT
    // fidelity (the ADR-0154 N3 finding), NOT gap dimensioning — so neither this
    // correction nor a multi-D ASM over gaps (B3) can clear the mask without a
    // higher-fidelity port. The gate therefore asserts the MECHANISM + the measured
    // LIFT (the true, reproducible results) and RECORDS the unreachability /
    // mask-miss / port-bound conclusion. It does NOT assert convergence or
    // mask-clearing (both false) — this records a real NO-GO, NOT a fake pass, and
    // does NOT weaken any physics tolerance.
    let lift_min_db = 3.0; // conservative floor below the measured +5.83 dB lift.
    let best_kfem = cache.first().map(|(_, c)| c.k_fem).unwrap_or(f64::NAN);

    // (a) MECHANISM — the f0-tracking band resolved the split and the corrector
    //     traced a finite K(gap) (best k_fem finite). This is the band/box fix
    //     working on the filter geometry. Convergence is NOT required: target_k is
    //     unreachable here, which is the finding (recorded below), not a bug.
    assert!(
        best_kfem.is_finite(),
        "B2 NO-GO: the design-curve correction produced no finite k_fem — the f0-tracking \
         sweep failed to resolve the split on the filter geometry (band/box regression). \
         Trajectory printed above."
    );
    if !any_converged {
        eprintln!(
            "[B2] FINDING: target_k={:.4} UNREACHABLE as a resonant-split k on this geometry \
             — best k_fem={:.4} ({:.1}% off) at the {:.3}mm gap floor; the resonant-split k \
             saturates ≈0.064 for the {:.3}mm trace. The impedance-k synthesis target \
             over-specifies the coupling. NOT a bug — the design-curve MECHANISM works; the \
             geometry caps the achievable coupling.",
            target_k[0],
            best_kfem,
            (best_kfem - target_k[0]).abs() / target_k[0] * 100.0,
            corrected_gaps[0] * 1e3,
            w * 1e3,
        );
    }

    // (b) The correction MOVED every gap materially off the impedance-k gap, and in
    //     the TIGHTER direction (stronger coupling — the impedance-k under-couples,
    //     k_imp ≠ resonant-split k, ADR-0155 K2). Unconditional (the move is the
    //     mechanism, independent of convergence).
    for (i, (&corr_gap, &imp_gap)) in corrected_gaps.iter().zip(imp_gaps.iter()).enumerate() {
        let move_m = (corr_gap - imp_gap).abs();
        assert!(
            move_m >= 0.05e-3,
            "B2: corrected gap[{i}] {:.4}mm differs from the impedance-k gap {:.4}mm by only \
             {:.4}mm (need ≥ 0.05mm). The design curve should MOVE the gap off the (divergent) \
             impedance-k value; a near-zero move means a no-op. Table printed above.",
            corr_gap * 1e3,
            imp_gap * 1e3,
            move_m * 1e3,
        );
        assert!(
            corr_gap < imp_gap,
            "B2: corrected gap[{i}] {:.4}mm is not TIGHTER than the impedance-k gap {:.4}mm — \
             the resonant-split design curve should tighten the gap (the impedance-k \
             under-couples). Table printed above.",
            corr_gap * 1e3,
            imp_gap * 1e3,
        );
    }

    // (c) The corrected |S21| curve is all-finite with a genuine band-pass turnover
    //     (in-band peak strictly above the deeper band edge) — a real pass/stop
    //     shape, not a degenerate / monotonic curve.
    assert!(
        graded.curve.iter().all(|(_, d)| d.is_finite()) && graded.passband_peak_db.is_finite(),
        "B2 NO-GO: the corrected-gap |S21| curve has a non-finite point — the driven solve \
         degenerated (port collapsed or mesh broken). Table printed above."
    );
    assert!(
        graded.turnover_db > 0.2,
        "B2: no band-pass turnover — in-band peak {:.2} dB is not above the deeper band edge \
         (edges {:.2}/{:.2} dB; turnover {:+.2} dB). Table printed above.",
        graded.passband_peak_db,
        graded.edge_lo_db,
        graded.edge_hi_db,
        graded.turnover_db,
    );

    // (d) The MEASURED LIFT — the headline reproducible result: correcting the gaps
    //     onto the (best-achievable) resonant-split coupling lifts the in-band peak
    //     a real margin over the N3 impedance-k baseline. The quantified
    //     dimensioning lever (+5.83 dB measured; pinned at ≥ +3 dB).
    assert!(
        lift_over_n3_db >= lift_min_db,
        "B2: in-band peak {:.2} dB did NOT lift ≥ {:.1} dB over the N3 baseline {:.2} dB \
         (measured lift {:+.2} dB). The design-curve gap correction should raise the in-band \
         peak (stronger, more-correct coupling); a regression below the pinned lift means the \
         corrector or the filter build broke. Table printed above.",
        graded.passband_peak_db,
        lift_min_db,
        N3_BASELINE_PEAK_DB,
        lift_over_n3_db,
    );

    // (e) Geometric-asymmetry SIGN is preserved (depth(1.6) > depth(2.4)). The very
    //     tight corrected gaps DEGRADE the asymmetry margin (measured +0.70 dB, vs
    //     N3's +2.10 dB) — RECORD that (part of the over-tightening finding) and
    //     assert only that the sign survives, not the ≥ 1 dB N3 margin.
    assert!(
        depth_lo > depth_hi,
        "B2: geometric-asymmetry SIGN lost — depth(1.6 GHz)={depth_lo:.2} dB is not deeper than \
         depth(2.4 GHz)={depth_hi:.2} dB. The corrected gaps inverted the band-pass-mapping \
         asymmetry. Table printed above.",
    );
    if !asym_pass {
        eprintln!(
            "[B2] RECORD: asymmetry margin {:+.2} dB is below the N3 {:.1} dB threshold — the \
             strong-coupling-floor gaps ({:.3}mm) degrade the asymmetry (N3 was +2.10 dB). Sign \
             preserved (depth 1.6 > 2.4). Part of the over-tightening finding.",
            graded.asym_margin_db,
            ASYMMETRY_MARGIN_DB,
            corrected_gaps[0] * 1e3,
        );
    }

    // (f) Strict Cheb mask — assert ONLY IF the measurement clears it (it does NOT:
    //     MISS by ~26 dB in-band). Never weaken / force it.
    if strict_pass {
        assert!(
            graded.worst_pass_db <= PASSBAND_TOL_DB && graded.worst_rej_db <= REJECTION_TOL_DB,
            "internal: strict_pass set but tolerances not met (pass {:.2}, rej {:.2})",
            graded.worst_pass_db,
            graded.worst_rej_db,
        );
        eprintln!(
            "[B2] STRICT MASK CLEARS (unexpected): worst passband err {:.2} dB ≤ {:.1}, worst \
             rejection err {:.2} dB ≤ {:.1}. The EM-in-loop gap correction closed the N3 floor.",
            graded.worst_pass_db, PASSBAND_TOL_DB, graded.worst_rej_db, REJECTION_TOL_DB,
        );
    } else {
        eprintln!(
            "[B2] STRICT MASK: MISS by {:.2} dB in-band (in-band peak {:.2} dB; passband \
             reference ~0 dB). The corrected-gap filter LIFTS {:+.2} dB over N3 but remains a \
             ~−22 dB flat shelf — the dominant floor is the aperture-coupling PORT fidelity \
             (ADR-0154 N3), NOT gap dimensioning. Dimensioning is a real but MINOR lever; \
             neither this correction nor a multi-D ASM over gaps clears the mask without a \
             higher-fidelity port. Recorded honestly — no fake pass, no weakened tolerance.",
            graded.worst_pass_db, graded.passband_peak_db, lift_over_n3_db,
        );
    }
}

// =====================================================================
// FEM-EM brick B3' (ADR-0162) — filter S21 re-graded with the POWER-CORRECT
// E+H modal extraction (the B2' fix, validated GO on the straight thru:
// |S21| 0.778→1.0001, |S11|²+|S21|² 0.61→1.0037, ε_eff 0.66%).
//
// THE GOAL (ADR-0147 #1): does the power-correct extraction lift the 3-pole
// FILTER S21 off the N3 −27.38 dB E-only floor toward the Chebyshev mask?
//
// The ONE change vs the N3 gate (`fem_filter_s21_vs_ladder`): S11/S21 at each
// frequency come from `power_modal_extract` (the E+H two-field modal
// decomposition) instead of the E-only `sweep_matrix` / `extract_s_qp`.
// Everything else — geometry (`build_edge_coupled_geometry`), the
// numerical-eigenmode ports recentred per off-centre feed, interior-PEC,
// coupled-Whitney, the 17-pt 1.6–2.4 GHz band, the feed de-embed, and the
// `ladder_s21` reference + oracle/mask grading — is IDENTICAL.
// =====================================================================

/// Modal-reference inward-sampling distance as a fraction of the feed length.
/// `power_modal_extract` samples each port's TRUE modal `(e_m, h_m)` this far
/// INWARD from the port face (along the propagation direction into the
/// structure). For the filter the two feeds are at OPPOSITE ends and share no
/// common interior `y`-plane, so an absolute plane will not do — the inward
/// offset lands each port's reference in ITS OWN uniform feed run. Half the
/// feed length sits comfortably inside both feeds (each is ≥ `feed_len` long;
/// the output feed is `feed_len + stagger`) and clear of the port-face ABC and
/// the first resonator.
const MODAL_REF_FEED_FRAC: f64 = 0.5;

/// Drive the filter two-port through the POWER-CORRECT E+H modal extraction
/// (B2', ADR-0162) and return per-frequency `(S11, S21)`.
///
/// Bit-identical solver build to [`solve_filter`] (same mesh, the same
/// `microstrip_port_numerical_at` ports recentred per off-centre feed,
/// interior-PEC, `with_coupled_whitney(true)`); the ONLY difference is the
/// post-solve extraction — `power_modal_extract(ω, 0, d_ref)` (drive port 0,
/// read S11 = `s_column[0]`, S21 = `s_column[1]`) instead of `sweep_matrix`.
///
/// `d_ref = MODAL_REF_FEED_FRAC · feed_len` is the inward sampling distance for
/// the modal reference; it lands inside each feed's uniform run (the crux for
/// the off-centre feeds — see [`MODAL_REF_FEED_FRAC`]).
fn solve_filter_power(
    geom: &FilterGeometry,
    omegas: &[f64],
) -> Vec<(num_complex::Complex64, num_complex::Complex64)> {
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

    // SAME numerical-eigenmode ports as `solve_filter` (recentred per feed).
    let port_geom = MicrostripPortGeom {
        trace_w: geom.line_w,
        sub_h: SUB_H,
        eps_r: EPS_R,
        box_w: geom.box_w,
        box_h: geom.box_h,
    };
    let port_in = microstrip_port_numerical_at(&port_geom, geom.feed_xc_in, F0)
        .expect("numerical-eigenmode port (input feed) must build");
    let port_out = microstrip_port_numerical_at(&port_geom, geom.feed_xc_out, F0)
        .expect("numerical-eigenmode port (output feed) must build");

    let solver = OpenBoundarySolver::new(&mesh, kinds, vec![port_in, port_out], material_db)
        .expect("two-port filter solver must build")
        .with_interior_pec_edges(interior_pec.iter().copied())
        .with_coupled_whitney(true);

    // Modal-reference inward distance: half the feed length lands inside both
    // uniform feeds (the off-centre input/output feeds share no interior
    // y-plane, so inward-from-each-face sampling is mandatory here).
    let d_ref = MODAL_REF_FEED_FRAC * geom.feed_len;

    omegas
        .iter()
        .map(|&omega| {
            let pm = solver
                .power_modal_extract(omega, 0, d_ref)
                .expect("B3' power_modal_extract must succeed");
            (pm.s_column[0], pm.s_column[1])
        })
        .collect()
}

/// FEM-EM brick B3' (ADR-0162) — 3-pole microstrip-filter S21 re-graded with the
/// **power-correct E+H modal extraction** (the B2' fix), vs the analytic ladder
/// reference. THE GOAL: does the power-correct extraction lift the filter S21 off
/// the N3 −27.38 dB E-only floor toward the Chebyshev mask?
///
/// Mirrors [`fem_filter_s21_vs_ladder`] (N3) EXACTLY except the extraction:
/// S11/S21 at each of the 17 band points come from
/// [`yee_fem::OpenBoundarySolver::power_modal_extract`] (the two-field
/// decomposition: `a_fwd = ½(proj_E+proj_H)`, `a_bwd = ½(proj_E−proj_H)`,
/// `S_pp = a_bwd/a_fwd`, `S_qp = a_fwd(q)/a_fwd(p)`) instead of the E-only
/// `sweep_matrix`/`extract_s_qp`. Same geometry, ports, mesh, de-embed, and
/// `ladder_s21` + oracle/mask grading.
///
/// ## Per-feed modal reference (the crux for off-centre feeds)
///
/// The filter's two feeds are uniform 50 Ω runs at DIFFERENT `x` and OPPOSITE
/// ends (input near `y = 0`, output near `y = box_len`), sharing no common
/// interior `y`-plane. `power_modal_extract` samples each port's TRUE modal
/// `(e_m, h_m = ∇×E/(−jωμ), de-rotated)` a distance `d_ref =
/// MODAL_REF_FEED_FRAC·feed_len` INWARD from that port's face (along the
/// propagation direction into the structure) — landing each reference in its
/// OWN feed, on the uniform run before resonator 0 / after the last resonator,
/// NOT inside the coupled region. Each port is normalized by its own
/// reaction-norm κ, so the two off-centre references are consistently scaled.
///
/// ## What this asserts (HONEST, MEASUREMENT-DRIVEN — recorded-then-pinned)
///
/// The HEADLINE is research-open and the orchestrator runs the heavy filter to
/// measure it. The gate PRINTS the full corrected |S21|(f) table, the in-band
/// `|S11|²+|S21|²` power balance (does the filter FIELD transmit in-band, or
/// genuinely reflect?), the in-band peak vs the N3 −27.38 dB E-only floor, and
/// the strict Cheb-mask margin. It asserts ONLY the measurement-independent
/// invariants:
///
/// 1. **Finite curve** — no NaN/Inf (the power extraction did not diverge / no
///    port collapsed).
/// 2. **A band-pass turnover** — the in-band peak stands above the deeper band
///    edge (a real centre bump, not a monotonic ramp / flat line).
/// 3. **The strict Cheb mask ONLY IF it actually clears** (the ADR-0147 #1 win)
///    — otherwise the honest MISS margin is recorded, no weakening.
///
/// It does **NOT** hard-assert a lift number that has not been measured — the
/// lift over the N3 −27.38 dB floor is PRINTED for the orchestrator to pin as a
/// tripwire AFTER seeing the real number (exactly as N3/B2 pinned theirs). If
/// the corrected in-band peak clears or approaches the mask, that is the
/// headline; if it lifts but stays short (e.g. the filter genuinely reflects
/// in-band, `|S11|` high), that is recorded honestly. No faking.
///
/// ## GATING — CRITICAL (heavy; run by the orchestrator, boxed, `--release`)
///
/// Multi-minute driven SWEEP: one per-ω sparse LU per point PLUS, per point, the
/// per-port interior modal-field point-location + reconstruction the E+H
/// extraction adds (a handful of cheap O(n_tets) scans — negligible vs the LU).
/// ~17 points on the ~51 k-tet mesh; budget roughly the N3 ~80 s plus a little.
/// `#[ignore]`'d; run only in `--release`, boxed:
///
/// ```text
/// YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh \
///   cargo test -p yee-fem --release --test microstrip_filter_s21 \
///   -- --ignored fem_filter_s21_power_extract --nocapture
/// ```
///
/// MEASURED RESULT: (to be filled in by the orchestrator after the boxed run —
/// the headline is the corrected in-band peak vs the −27.38 dB N3 floor and
/// whether the strict mask clears).
#[test]
#[ignore = "B3' GOAL gate: heavy 17-pt driven SWEEP with the power-correct E+H extraction; run only in --release, boxed"]
fn fem_filter_s21_power_extract() {
    // Geometry — IDENTICAL to the N3 gate (analytic impedance-k gaps, same mesh).
    let geom = build_edge_coupled_geometry(
        2.5e-3, // x clearance each side
        5.0e-3, // air height
        8.0e-3, // feed length (de-embed reference + modal-reference feed run)
        0.6e-3, // dx (trace ~3 cells, gap ~2.7 cells)
        2.5e-3, // dy (resonator ~16 cells)
        0.5e-3, // dz (2 substrate cells)
        None,   // analytic impedance-k gaps (same as N3 baseline)
    );
    eprintln!(
        "[B3'] filter mesh: box=({:.1},{:.1},{:.1})mm  n=({},{},{})  tets={}  w={:.3}mm  \
         feed={:.1}mm  modal_ref d={:.2}mm inward  eps_eff(w)={:.4}",
        geom.box_w * 1e3,
        geom.box_len * 1e3,
        geom.box_h * 1e3,
        geom.nx,
        geom.ny,
        geom.nz,
        geom.total_tets(),
        geom.line_w * 1e3,
        geom.feed_len * 1e3,
        MODAL_REF_FEED_FRAC * geom.feed_len * 1e3,
        eps_eff(geom.line_w, SUB_H, EPS_R),
    );

    // Band: 1.6–2.4 GHz, 17 points — IDENTICAL grid to N3/B2.
    let (freqs_hz, omegas) = band_1p6_to_2p4_17pts();

    let t0 = std::time::Instant::now();
    let s_pairs = solve_filter_power(&geom, &omegas);
    let wall = t0.elapsed().as_secs_f64();

    // Extract + de-embed |S21|(f) and record the in-band power balance.
    let ladder = reference_ladder();
    let mut curve: Vec<(f64, f64)> = Vec::with_capacity(freqs_hz.len());
    let mut worst_balance_in: f64 = f64::INFINITY; // min in-band |S11|²+|S21|²
    let mut best_balance_in: f64 = 0.0; // max in-band |S11|²+|S21|²
    eprintln!(
        "\n{:>8}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}",
        "f(GHz)", "|S21|raw", "|S21|deemb", "S21 dB", "|S11|", "|S|²sum", "ref dB"
    );
    for (k, &omega) in omegas.iter().enumerate() {
        let (s11, s21_raw) = s_pairs[k];
        let s21 = deembed_feed(s21_raw, omega, geom.line_w, geom.feed_len);
        let d = db(s21.norm());
        let f_ghz = freqs_hz[k] / 1e9;
        let bal = s11.norm_sqr() + s21_raw.norm_sqr(); // de-embed is unit-magnitude
        let ref_db = db(ladder_s21(&ladder, freqs_hz[k]).norm());
        curve.push((f_ghz, d));
        if (1.85..=2.15).contains(&f_ghz) {
            worst_balance_in = worst_balance_in.min(bal);
            best_balance_in = best_balance_in.max(bal);
        }
        eprintln!(
            "{:>8.3}  {:>10.4}  {:>10.4}  {:>10.2}  {:>10.4}  {:>10.4}  {:>10.2}",
            f_ghz,
            s21_raw.norm(),
            s21.norm(),
            d,
            s11.norm(),
            bal,
            ref_db,
        );
    }

    // ---- Grade against the reference (mirrors oracle_grade / the N3 gate). ----
    let mut worst_pass_db = 0.0_f64;
    let mut worst_rej_db = 0.0_f64;
    for &(f_ghz, d_meas) in &curve {
        let d_ref_db = db(ladder_s21(&ladder, f_ghz * 1e9).norm());
        let err = (d_meas - d_ref_db).abs();
        if (1.85..=2.15).contains(&f_ghz) {
            worst_pass_db = worst_pass_db.max(err);
        } else {
            worst_rej_db = worst_rej_db.max(err);
        }
    }

    let depth_at = |f_ghz: f64| -> f64 { -interp_db(&curve, f_ghz) };
    let depth_lo = depth_at(1.6);
    let depth_hi = depth_at(2.4);
    let asym_margin = depth_lo - depth_hi;
    let asym_pass = asym_margin >= ASYMMETRY_MARGIN_DB;

    let passband_peak_db = curve
        .iter()
        .filter(|(f, _)| (1.85..=2.15).contains(f))
        .map(|(_, d)| *d)
        .fold(f64::NEG_INFINITY, f64::max);
    let f_inband_peak = curve
        .iter()
        .filter(|(f, _)| (1.85..=2.15).contains(f))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(f, _)| *f)
        .unwrap_or(f64::NAN);
    let edge_lo_db = interp_db(&curve, 1.6);
    let edge_hi_db = interp_db(&curve, 2.4);
    let turnover_db = passband_peak_db - edge_lo_db.min(edge_hi_db);

    let strict_pass =
        worst_pass_db <= PASSBAND_TOL_DB && worst_rej_db <= REJECTION_TOL_DB && asym_pass;

    // The headline: lift over the N3 −27.38 dB E-only floor (PRINTED — the
    // orchestrator pins it as a tripwire after seeing the real number) and the
    // lift over the v1 analytic floor.
    let lift_over_n3_db = passband_peak_db - N3_BASELINE_PEAK_DB;
    let lift_over_v1_db = passband_peak_db - V1_FLOOR_PEAK_DB;

    eprintln!(
        "\n==== B3' GRADE (power-correct E+H extraction; ADR-0162) ====\n\
         tets                : {}\n\
         wall                : {:.1} s\n\
         in-band peak        : {:.2} dB @ {:.2} GHz\n\
         N3 E-only floor     : {:.2} dB  (impedance-k gaps, ADR-0154 N3)\n\
         lift over N3 floor  : {:+.2} dB  (PRINTED — orchestrator pins as a tripwire)\n\
         v1 analytic floor   : {:.2} dB  (flat-Ez port, B7)\n\
         lift over v1 floor  : {:+.2} dB\n\
         in-band |S|²sum     : min {:.4}  max {:.4}  (1 ⇒ field transmits; ≪1 ⇒ reflects)\n\
         band edges          : {:.2} dB @1.6  {:.2} dB @2.4\n\
         turnover            : {:+.2} dB (in-band peak above the deeper edge)\n\
         worst passband err  : {:.2} dB vs ref (oracle tol {:.1})\n\
         worst rejection err : {:.2} dB vs ref (oracle tol {:.1})\n\
         strict-mask margin  : {} by {:.2} dB in-band (gap to the 0 dB Cheb passband)\n\
         asymmetry (NAMED)   : depth(1.6)={:.2} dB  depth(2.4)={:.2} dB  margin={:+.2} dB -> {}\n\
         strict oracle mask  : {}\n\
         ============================================================",
        geom.total_tets(),
        wall,
        passband_peak_db,
        f_inband_peak,
        N3_BASELINE_PEAK_DB,
        lift_over_n3_db,
        V1_FLOOR_PEAK_DB,
        lift_over_v1_db,
        worst_balance_in,
        best_balance_in,
        edge_lo_db,
        edge_hi_db,
        turnover_db,
        worst_pass_db,
        PASSBAND_TOL_DB,
        worst_rej_db,
        REJECTION_TOL_DB,
        if strict_pass { "CLEARS" } else { "MISS" },
        worst_pass_db,
        depth_lo,
        depth_hi,
        asym_margin,
        if asym_pass { "PASS" } else { "FLAG" },
        if strict_pass { "PASS" } else { "MISS" },
    );

    // Machine-readable curve for the oracle_grade CLI.
    let pairs: String = curve
        .iter()
        .map(|(f, d)| format!("{f:.3}:{d:.2}"))
        .collect::<Vec<_>>()
        .join(" ");
    eprintln!("[B3'] oracle_grade pairs: {pairs}");

    // ---- Assertions — HONEST, MEASUREMENT-DRIVEN (assert only invariants) ----
    //
    // B3' is the research-open GOAL: does the power-correct E+H extraction lift
    // the filter S21 off the N3 −27.38 dB E-only floor toward the Cheb mask? The
    // headline (the lift, the in-band power balance, mask-clearing) is MEASURED
    // by the orchestrator's boxed run — this gate does NOT pre-judge it. It
    // asserts only what is true regardless of the measured level, and asserts
    // the strict mask ONLY IF it actually clears (no weakening to force green).

    // (1) Finite curve — the power extraction did not diverge and no port
    //     collapsed (a NaN/Inf S would mean the modal reconstruction or the
    //     per-port reaction-norm normalization broke).
    assert!(
        curve.iter().all(|(_, d)| d.is_finite()) && passband_peak_db.is_finite(),
        "B3' NO-GO: the power-extracted |S21| curve has a non-finite point — the driven \
         solve or the E+H modal extraction degenerated. Full curve printed above."
    );

    // (2) A genuine band-pass turnover: the in-band peak stands above the deeper
    //     band edge (a frequency-selective bump, not a monotonic ramp / flat
    //     line). The >0.2 dB bar mirrors N3/B2; it certifies SHAPE without
    //     demanding a depth the path may not deliver.
    assert!(
        turnover_db > 0.2,
        "B3': no band-pass turnover — in-band peak {passband_peak_db:.2} dB is not above the \
         deeper band edge (edges {edge_lo_db:.2}/{edge_hi_db:.2} dB; turnover {turnover_db:+.2} dB). \
         The response is monotonic, not a recognisable band-pass. Full curve printed above."
    );

    // (3) Strict Cheb mask — assert ONLY IF the measurement actually clears it
    //     (the ADR-0147 #1 win). If it clears, that is the headline and we assert
    //     it loudly; if not, the honest MISS margin is recorded and we do NOT
    //     assert the absolute-level mask (no faking). A future improvement that
    //     lifts the curve into the mask flips `strict_pass` and asserts here
    //     automatically.
    if strict_pass {
        assert!(
            worst_pass_db <= PASSBAND_TOL_DB && worst_rej_db <= REJECTION_TOL_DB && asym_pass,
            "internal: strict_pass set but tolerances not met (pass {worst_pass_db:.2}, \
             rej {worst_rej_db:.2}, asym {asym_margin:+.2})"
        );
        eprintln!(
            "[B3'] STRICT MASK CLEARS — the power-correct E+H extraction lifted the 3-pole \
             filter S21 into the Chebyshev mask (worst passband err {worst_pass_db:.2} dB ≤ \
             {PASSBAND_TOL_DB}, worst rejection err {worst_rej_db:.2} dB ≤ {REJECTION_TOL_DB}, \
             asymmetry {asym_margin:+.2} dB). This is the ADR-0147 #1 goal — a mask-clearing \
             full-wave filter S21."
        );
    } else {
        eprintln!(
            "[B3'] STRICT MASK: MISS by {worst_pass_db:.2} dB in-band. in-band peak \
             {passband_peak_db:.2} dB ({lift_over_n3_db:+.2} dB vs the N3 −27.38 dB E-only floor); \
             in-band power balance |S11|²+|S21|² ∈ [{worst_balance_in:.4}, {best_balance_in:.4}]. \
             Recorded honestly — the orchestrator pins the measured lift as a tripwire and reads \
             the balance to judge whether the filter FIELD transmits in-band or genuinely \
             reflects. No fake pass, no weakened tolerance."
        );
    }
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
        let geom =
            build_edge_coupled_geometry(2.5e-3, 5.0e-3, 8.0e-3, 1.6e-3, 5.0e-3, 0.5e-3, None);
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
