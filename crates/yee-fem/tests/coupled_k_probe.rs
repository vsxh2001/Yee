//! DE-RISK PROBE (measurement, NOT a gate, NOT merged) — coupled-microstrip
//! resonator-PAIR coupling coefficient `k` from a frequency-domain FEM S21
//! sweep.
//!
//! ## The one question
//!
//! Does a frequency-domain FEM `sweep_matrix` of a coupled-microstrip
//! resonator PAIR show **two resolvable transmission peaks** whose
//! frequency-split yields a coupling coefficient `k ≈` the analytic reference?
//! If yes, the FEM driven-sweep (ADR-0153/0154) is the wall-free EM back-end
//! for the filter app's F1.1b inter-resonator coupling extraction; the FDTD
//! resonant-split method was abandoned (ADR-0108: no box is simultaneously
//! high-Q and non-confining). The FEM sweep is frequency-domain — no cavity
//! wall, no time-stepping — and `k` depends only on peak *locations* (robust
//! to the absolute |S21| floor that limits FEM transmission, exactly as B4's
//! ε_eff was; see `microstrip_eeff.rs`).
//!
//! ## Two definitions of k (read before reading the verdict)
//!
//! There are two distinct, both-legitimate "coupling coefficients" in play,
//! and they are NOT the same number for full-length parallel coupling:
//!
//! * **`k_split` = (f_hi² − f_lo²)/(f_hi² + f_lo²)** — the *resonator* coupling
//!   from the split of the two resonant peaks (the brief's formula, and what a
//!   filter designer extracts from EM). For two parallel λ_g/2 resonators of
//!   *fixed physical length* the even/odd modes resonate at `f_e =
//!   c/(2L√ε_eff,e)` and `f_o = c/(2L√ε_eff,o)`, so the split is driven by the
//!   even/odd **ε_eff** difference.
//! * **`k_imp` = (Z0e − Z0o)/(Z0e + Z0o)** — the coupled-*line* voltage
//!   coupling ([`yee_layout::coupling_coefficient`], Kirschning-Jansen). This
//!   is the brief's named analytic reference.
//!
//! For full-length parallel coupling `k_imp / k_split` ranges from ≈3.3 (tight
//! gap S/W=0.3) down to ≈1.1 (weak gap S/W=2): they converge only in the weak
//! limit. We therefore pick a **weak gap (S = 2·W = 2 mm)** where the two
//! definitions agree to ≈11 %, so the brief's prescribed comparison
//! (`k_split` formula vs `coupling_coefficient`) is physically meaningful and
//! the coupling is weak enough that the two peaks do not smear together. The
//! probe reports BOTH analytic references for transparency.
//!
//! ## Geometry (matches the shipped FEM line/filter work — FR-4, h = 1 mm)
//!
//! ```text
//!   substrate h   = 1 mm     ε_r = 4.4 (FR-4)
//!   strip W       = 1 mm     (W/h = 1, ~the B4 line)
//!   gap S         = 2 mm     (S/W = 2, weak coupling: k_imp ≈ k_split)
//!   resonator L   = λ_g/2 at f0 = 2.4 GHz using single-line ε_eff ≈ 35.1 mm
//!   box           = 9 × ~49 × 6 mm  (walls ≥ 2.5·h clear in x — B4 box-loading
//!                   finding; 6 mm tall = open-half-space air)
//!   feed coupling = 1 mm end-gap (1 dy cell) tapping each resonator open end
//! ```
//!
//! Axes match `layered_microstrip_filter_mesh` (B2/B7): `x` cross-section,
//! `y` propagation (feed-to-feed), `z` substrate-normal (ground z=0, trace
//! z=sub_h). Ports on the `y=0` / `y=box_len` end-caps.
//!
//! ## Excitation — weakly-coupled feeds
//!
//! TWO ports, each WEAKLY coupled to ONE resonator (weak coupling is essential
//! — over-coupling smears the split). A straight feed line runs from each
//! port plane up to a small **end-gap** before its resonator's open end; the
//! gap turns the feed into a weak capacitive tap. The feeds are off-centre in
//! `x` (each aligned with its own resonator), so each numerical-eigenmode
//! wave-port is RECENTRED on its feed's `x` via
//! [`yee_fem::microstrip_port_numerical_at`] (the filter test's per-feed
//! recentre, ADR-0154 N3). `with_coupled_whitney(true)` is MANDATORY (B4
//! finding: the lumped-centroid port collapses the absorbing block for the
//! substrate-normal `E_z` mode).
//!
//! The 1 mm end-gap is one `dy = 1 mm` cell — resolved by the mesh. (A finer
//! gap at the coarser `dy = 2.5 mm` filter pitch would be sub-cell and
//! unresolved; this probe pays for `dy = 1 mm` to resolve a real coupling
//! gap, landing ≈ 63 k tets in a long-thin box that the per-ω faer sparse LU
//! fits in a 14 g box.)
//!
//! ## GATING — CRITICAL
//!
//! Multi-minute driven SWEEP (one per-ω sparse LU per frequency point). The
//! probe is `#[ignore]`'d so the debug `cargo test --workspace` never runs it;
//! run only in `--release`, boxed:
//!
//! ```text
//! YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh bash -c '\
//!   cargo test -p yee-fem --release --test coupled_k_probe \
//!   -- --ignored fem_coupled_k_probe --nocapture'
//! ```

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::{
    FaceKind, MaterialDatabase, MicrostripPortGeom, OpenBoundarySolver, SParametersMatrix,
    TraceRect, layered_microstrip_filter_mesh, microstrip_port_numerical_at,
};
use yee_layout::{coupled_microstrip, coupling_coefficient, eps_eff};
use yee_mesh::TetMesh3D;

// ---------------------------------------------------------------------
// Fixed substrate / geometry constants.
// ---------------------------------------------------------------------

/// Speed of light (m/s).
const C0: f64 = 299_792_458.0;
/// Substrate height (m): 1 mm FR-4.
const SUB_H: f64 = 1.0e-3;
/// FR-4 relative permittivity.
const EPS_R: f64 = 4.4;
/// Strip width (m).
const W: f64 = 1.0e-3;
/// Edge-to-edge coupling gap between the two resonators (m). S/W = 2 → weak
/// coupling, where `k_imp ≈ k_split` (see module docs).
const S: f64 = 2.0e-3;
/// Resonator centre frequency (Hz). L is set to λ_g/2 here.
const F0: f64 = 2.4e9;

/// PEC-shield clearance each side in x (m). B4: walls ≥ 2.5·h clear or they
/// load the line / pull ε_eff below the open-microstrip HJ value.
const CLEARANCE_X: f64 = 2.5e-3;
/// Air height above the substrate (m): 5 mm → 6 mm box (open-half-space).
const AIR_H: f64 = 5.0e-3;
/// Straight feed-line length at each end (m), before the coupling gap.
const FEED_RUN: f64 = 6.0e-3;
/// Feed-to-resonator end-gap (m): one `dy` cell, weak capacitive tap.
const FEED_GAP: f64 = 1.0e-3;

/// Cross-section pitch (m): dx = 0.5 mm → W = 2 cells, S = 4 cells.
const DX: f64 = 0.5e-3;
/// Propagation pitch (m): dy = 1.0 mm → resolves the 1 mm feed gap (1 cell).
const DY: f64 = 1.0e-3;
/// Substrate-normal pitch (m): dz = 0.5 mm → 2 substrate z-cells.
const DZ: f64 = 0.5e-3;

// ---------------------------------------------------------------------
// Geometry construction.
// ---------------------------------------------------------------------

/// Resolved coupled-pair geometry in mesh world coordinates.
struct CoupledGeometry {
    box_w: f64,
    box_len: f64,
    box_h: f64,
    traces: Vec<TraceRect>,
    nx: usize,
    ny: usize,
    nz: usize,
    /// Resonator length L = λ_g/2 (m).
    res_l: f64,
    /// Strip width = the feed wave-port width (m).
    line_w: f64,
    /// One-sided straight feed length (m) — the de-embed reference length.
    feed_len: f64,
    /// x-centre (m) of the INPUT feed (port 0, y=0 end-cap).
    feed_xc_in: f64,
    /// x-centre (m) of the OUTPUT feed (port 1, y=box_len end-cap).
    feed_xc_out: f64,
}

impl CoupledGeometry {
    fn total_tets(&self) -> usize {
        self.nx * self.ny * self.nz * 6
    }
}

/// Resonator length λ_g/2 at `F0` using the single-line Hammerstad-Jensen
/// ε_eff (the resonator is dominantly a single line; the coupling perturbs the
/// resonance, which is exactly what we measure).
fn resonator_length() -> f64 {
    let eeff = eps_eff(W, SUB_H, EPS_R);
    let lam_g = C0 / F0 / eeff.sqrt();
    lam_g / 2.0
}

/// Build the coupled resonator-pair geometry.
///
/// Two identical parallel λ_g/2 strips along `y`, separated by gap `S` in `x`.
/// Each is open-open (floating). Two straight feeds run from the `y=0` /
/// `y=box_len` port planes up to a `FEED_GAP` end-gap before the near open end
/// of resonator A (input) / resonator B (output) respectively.
///
/// y-layout: `FEED_RUN | FEED_GAP | L (resonators) | FEED_GAP | FEED_RUN`.
fn build_coupled_geometry() -> CoupledGeometry {
    let res_l = resonator_length();

    // x: two strips with a gap S, plus clearance both sides.
    //   resonator A: x ∈ [CLEARANCE_X,           CLEARANCE_X + W]
    //   resonator B: x ∈ [CLEARANCE_X + W + S,    CLEARANCE_X + 2W + S]
    let res_a_x0 = CLEARANCE_X;
    let res_b_x0 = CLEARANCE_X + W + S;
    let pattern_x_hi = res_b_x0 + W;
    let box_w = pattern_x_hi + CLEARANCE_X;

    // y: feeds + gaps + resonators.
    let y_res_lo = FEED_RUN + FEED_GAP;
    let y_res_hi = y_res_lo + res_l;
    let box_len = y_res_hi + FEED_GAP + FEED_RUN;

    // z: substrate + air. Snap so sub_h lands on a z-plane (it does: 1mm/0.5mm).
    let nz = ((SUB_H + AIR_H) / DZ).round() as usize;
    let nz_sub = (SUB_H / DZ).round().max(1.0) as usize;
    let dz_exact = SUB_H / nz_sub as f64;
    let box_h = dz_exact * nz as f64;

    // Traces: two resonators + two feeds.
    let mut traces = Vec::with_capacity(4);
    // Resonator A (coupled to the input feed).
    traces.push(TraceRect::new(res_a_x0, y_res_lo, W, res_l));
    // Resonator B (coupled to the output feed).
    traces.push(TraceRect::new(res_b_x0, y_res_lo, W, res_l));
    // Input feed: y ∈ [0, FEED_RUN], aligned with resonator A. The FEED_GAP
    // separates the feed top (y = FEED_RUN) from resonator A's bottom open end
    // (y = y_res_lo = FEED_RUN + FEED_GAP).
    traces.push(TraceRect::new(res_a_x0, 0.0, W, FEED_RUN));
    // Output feed: y ∈ [box_len − FEED_RUN, box_len], aligned with resonator B.
    // The FEED_GAP separates the feed bottom from resonator B's top open end.
    let out_feed_y0 = box_len - FEED_RUN;
    traces.push(TraceRect::new(res_b_x0, out_feed_y0, W, FEED_RUN));

    let nx = (box_w / DX).round().max((W / DX).ceil().max(1.0)) as usize;
    let ny = (box_len / DY).round().max(1.0) as usize;

    CoupledGeometry {
        box_w,
        box_len,
        box_h,
        traces,
        nx,
        ny,
        nz,
        res_l,
        line_w: W,
        feed_len: FEED_RUN,
        feed_xc_in: res_a_x0 + W / 2.0,
        feed_xc_out: res_b_x0 + W / 2.0,
    }
}

// ---------------------------------------------------------------------
// Face classification — ports on the y = 0 / y = box_len end-caps.
// (Identical to the filter test; the interior just has a different footprint.)
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

fn classify_faces(centroids: &[Vector3<f64>], box_len: f64) -> Vec<FaceKind> {
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

/// Build the two-port driven solver for the coupled pair and run the sweep.
/// Trace (resonators + feeds) AND ground tagged interior-PEC (B1); the two
/// y-end-caps carry the numerical-eigenmode wave-port recentred per feed
/// (ADR-0154 N3) with `with_coupled_whitney(true)` (mandatory, B4).
fn solve_coupled(geom: &CoupledGeometry, omegas: &[f64]) -> SParametersMatrix {
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
    .expect("coupled-pair mesh must build");

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
    let kinds = classify_faces(&centroids, geom.box_len);
    drop(picker);

    // Numerical-eigenmode wave-port, recentred on each off-centre feed (the
    // feeds sit over their resonators, not at box_w/2). β stays analytic-HJ on
    // the feed width; only the modal SHAPE is numerical. Each face needs its
    // own call (boxed closures are not Clone).
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
        .expect("two-port coupled-pair solver must build")
        .with_interior_pec_edges(interior_pec.iter().copied())
        .with_coupled_whitney(true);

    solver
        .sweep_matrix(omegas)
        .expect("driven sweep_matrix must succeed")
}

fn db(mag: f64) -> f64 {
    20.0 * mag.log10()
}

/// Find local maxima in the `(f_ghz, |S21|_linear)` curve that stand above
/// both neighbours, returned as `(f_ghz, mag)` sorted by descending mag.
fn local_maxima(curve: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut peaks = Vec::new();
    for i in 1..curve.len().saturating_sub(1) {
        let (f, m) = curve[i];
        if m > curve[i - 1].1 && m > curve[i + 1].1 {
            peaks.push((f, m));
        }
    }
    peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    peaks
}

// =====================================================================
// THE PROBE (#[ignore]'d — run a real multi-point LU sweep)
// =====================================================================

/// DE-RISK PROBE — coupled-microstrip resonator-pair coupling `k` from a FEM
/// S21 sweep. Builds the pair, drives `sweep_matrix` over a band straddling
/// the even+odd resonances, finds the two transmission peaks, computes
/// `k_split = (f_hi² − f_lo²)/(f_hi² + f_lo²)`, and compares to BOTH analytic
/// references (`k_imp` from `coupling_coefficient`, and an `k_eps` derived from
/// the even/odd ε_eff that `k_split` should track for full-length coupling).
///
/// Prints the full swept |S21|(f) so the shape is auditable, then a verdict
/// block (GO / PARTIAL / NO-GO). This is a MEASUREMENT — it does not assert a
/// "pass"; it asserts only that the sweep ran (finite, non-degenerate) so a
/// broken pipeline surfaces, and prints the decision.
#[test]
#[ignore = "DE-RISK PROBE: multi-minute driven SWEEP (one per-ω sparse LU per point); run only in --release, boxed"]
fn fem_coupled_k_probe() {
    let geom = build_coupled_geometry();

    // ---- Analytic references --------------------------------------------
    let cm = coupled_microstrip(W, S, SUB_H, EPS_R);
    let k_imp = coupling_coefficient(&cm); // (Z0e−Z0o)/(Z0e+Z0o)
    // For two parallel λ_g/2 resonators of fixed length L, the even/odd modes
    // resonate at f_e = c/(2L√ε_eff,e), f_o = c/(2L√ε_eff,o). This is the
    // split k_split *should* reproduce (the physical reference for the FEM
    // measurement), distinct from the coupled-line k_imp.
    let f_even = C0 / (2.0 * geom.res_l * cm.eps_eff_e.sqrt());
    let f_odd = C0 / (2.0 * geom.res_l * cm.eps_eff_o.sqrt());
    let (f_e_lo, f_e_hi) = (f_even.min(f_odd), f_even.max(f_odd));
    let k_eps = (f_e_hi * f_e_hi - f_e_lo * f_e_lo) / (f_e_hi * f_e_hi + f_e_lo * f_e_lo);

    eprintln!(
        "[probe] geometry: W={:.3}mm S={:.3}mm h={:.3}mm eps_r={} L(λg/2)={:.3}mm f0={:.2}GHz\n\
         [probe] box=({:.2},{:.2},{:.2})mm  n=({},{},{})  tets={}  feed={:.1}mm gap={:.1}mm\n\
         [probe] feed xc: in={:.3}mm out={:.3}mm  (box_w/2={:.3}mm)\n\
         [probe] ANALYTIC: Z0e={:.2}Ω Z0o={:.2}Ω  eps_eff_e={:.4} eps_eff_o={:.4}\n\
         [probe]   k_imp  = (Z0e−Z0o)/(Z0e+Z0o)               = {:.4} ({:.2}%)   <- `coupling_coefficient`\n\
         [probe]   k_eps  = even/odd ε_eff split (f_e={:.4},f_o={:.4}GHz) = {:.4} ({:.2}%)   <- what k_split tracks\n\
         [probe]   (k_imp/k_eps ratio = {:.2}; they converge in the weak-gap limit)",
        W * 1e3,
        S * 1e3,
        SUB_H * 1e3,
        EPS_R,
        geom.res_l * 1e3,
        F0 / 1e9,
        geom.box_w * 1e3,
        geom.box_len * 1e3,
        geom.box_h * 1e3,
        geom.nx,
        geom.ny,
        geom.nz,
        geom.total_tets(),
        geom.feed_len * 1e3,
        FEED_GAP * 1e3,
        geom.feed_xc_in * 1e3,
        geom.feed_xc_out * 1e3,
        geom.box_w * 1e3 / 2.0,
        cm.z0e_ohm,
        cm.z0o_ohm,
        cm.eps_eff_e,
        cm.eps_eff_o,
        k_imp,
        k_imp * 100.0,
        f_even / 1e9,
        f_odd / 1e9,
        k_eps,
        k_eps * 100.0,
        k_imp / k_eps,
    );

    // ---- Sweep band: straddle f_even / f_odd with margin --------------------
    // Centre the band on f0 and span generously past both analytic resonances
    // so the two FEM peaks (which may shift from the analytic estimate) are
    // captured with valley + shoulders. ≥40 points across the split.
    let f_lo_hz = 2.10e9;
    let f_hi_hz = 2.70e9;
    let n_pts = 61; // 10 MHz spacing → ~14 points across a ~140 MHz split
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|i| f_lo_hz + (f_hi_hz - f_lo_hz) * (i as f64) / ((n_pts - 1) as f64))
        .collect();
    let omegas: Vec<f64> = freqs_hz.iter().map(|f| 2.0 * PI * f).collect();

    eprintln!(
        "[probe] sweep: {:.3}–{:.3} GHz, {} pts ({:.0} MHz step)",
        f_lo_hz / 1e9,
        f_hi_hz / 1e9,
        n_pts,
        (f_hi_hz - f_lo_hz) / 1e6 / (n_pts as f64 - 1.0),
    );

    let t0 = std::time::Instant::now();
    let sweep = solve_coupled(&geom, &omegas);
    let wall = t0.elapsed().as_secs_f64();

    // ---- Extract |S21|(f). No de-embed needed: k depends only on peak
    //      LOCATIONS, which a unit-magnitude feed-phase rotation cannot move. --
    let mut curve_db: Vec<(f64, f64)> = Vec::with_capacity(n_pts);
    let mut curve_lin: Vec<(f64, f64)> = Vec::with_capacity(n_pts);
    eprintln!("\n{:>8}  {:>12}  {:>10}", "f(GHz)", "|S21| lin", "S21 dB");
    for (i, &_omega) in omegas.iter().enumerate() {
        let s = &sweep.s[i];
        let mag = s[(1, 0)].norm();
        let f_ghz = freqs_hz[i] / 1e9;
        curve_lin.push((f_ghz, mag));
        curve_db.push((f_ghz, db(mag)));
        eprintln!("{:>8.3}  {:>12.6}  {:>10.2}", f_ghz, mag, db(mag));
    }

    // ---- Find the two transmission peaks --------------------------------
    let peaks = local_maxima(&curve_lin);
    eprintln!("\n[probe] local maxima (desc mag): {peaks:?}");

    // Build the verdict from the top-two peaks (if any), requiring a valley
    // meaningfully below both for "resolvable". The measurement outputs
    // (resolvable / k_fem / peak freqs / valley) are mutated per branch; the
    // (verdict, reason) pair is the VALUE of the match so there is no dead
    // initialiser.
    let mut peaks_resolvable = false;
    let mut k_fem = f64::NAN;
    let mut f_lo_peak = f64::NAN;
    let mut f_hi_peak = f64::NAN;
    let mut valley_db = f64::NAN;

    let (verdict, reason): (&str, String) = if peaks.len() >= 2 {
        // Take the two strongest peaks; order by frequency.
        let mut top2 = [peaks[0], peaks[1]];
        top2.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let (fa, ma) = top2[0];
        let (fb, mb) = top2[1];
        f_lo_peak = fa;
        f_hi_peak = fb;

        // Valley = minimum |S21| strictly between the two peaks.
        let valley = curve_lin
            .iter()
            .filter(|(f, _)| *f > fa && *f < fb)
            .map(|(_, m)| *m)
            .fold(f64::INFINITY, f64::min);
        valley_db = db(valley);
        let peak_lo_db = db(ma.min(mb)); // shallower of the two peaks
        // Resolvable ⇒ a clear valley: the dip between peaks sits meaningfully
        // (≥ 1.5 dB) below the SHALLOWER peak.
        let valley_depth_db = peak_lo_db - valley_db;
        peaks_resolvable = valley.is_finite() && valley_depth_db >= 1.5 && fb > fa;

        let f2_lo = fa * 1e9;
        let f2_hi = fb * 1e9;
        k_fem = (f2_hi * f2_hi - f2_lo * f2_lo) / (f2_hi * f2_hi + f2_lo * f2_lo);

        eprintln!(
            "[probe] two strongest peaks: f_lo={fa:.4} GHz ({:.2} dB)  f_hi={fb:.4} GHz ({:.2} dB)\n\
             [probe] valley between = {valley_db:.2} dB  (depth below shallower peak = {valley_depth_db:.2} dB; need ≥1.5)",
            db(ma),
            db(mb),
        );

        if peaks_resolvable {
            // Compare k_fem to BOTH references; the physical match is k_eps
            // (k_split-vs-k_split), with k_imp reported alongside.
            let err_eps = (k_fem - k_eps).abs() / k_eps * 100.0;
            let err_imp = (k_fem - k_imp).abs() / k_imp * 100.0;
            // GO band: ≲30% vs the physical k_eps reference (the like-for-like
            // resonant-split comparison). Report the k_imp error too.
            if err_eps <= 30.0 {
                (
                    "GO",
                    format!(
                        "two clean resolvable peaks; k_fem={k_fem:.4} vs k_eps={k_eps:.4} \
                         ({err_eps:.1}% err, within 30%); vs k_imp={k_imp:.4} ({err_imp:.1}%)"
                    ),
                )
            } else {
                (
                    "PARTIAL",
                    format!(
                        "peaks resolve cleanly but k off: k_fem={k_fem:.4} vs k_eps={k_eps:.4} \
                         ({err_eps:.1}% err > 30%; vs k_imp={k_imp:.4} {err_imp:.1}%). Likely \
                         levers: tune L so f0 centres the band, widen/narrow S, or refine the \
                         coupling-gap mesh"
                    ),
                )
            }
        } else {
            (
                "PARTIAL",
                format!(
                    "two maxima exist but the valley is shallow (depth {valley_depth_db:.2} dB < 1.5) \
                     — peaks not cleanly separated; likely over-coupled (reduce gap) or mesh too \
                     coarse across the split. k_fem (provisional) = {k_fem:.4} vs k_eps {k_eps:.4}"
                ),
            )
        }
    } else if peaks.len() == 1 {
        f_lo_peak = peaks[0].0;
        (
            "NO-GO",
            format!(
                "only ONE transmission peak at {:.4} GHz — the even/odd modes did not split \
                 (coupling too weak to resolve at this band/step, OR the two resonators are not \
                 both excited). Widen the gap-coupling or check both feeds tap.",
                peaks[0].0
            ),
        )
    } else {
        (
            "NO-GO",
            "NO transmission peak in band — the resonances are buried by the |S21| floor \
             or fall outside the swept band. Print the spectrum; re-centre the band on the \
             measured resonance."
                .to_string(),
        )
    };

    // ---- Verdict block --------------------------------------------------
    eprintln!(
        "\n==== COUPLED-K PROBE VERDICT ====\n\
         geometry           : W={:.3}mm S={:.3}mm h={:.3}mm eps_r={} L={:.3}mm  box=({:.1},{:.1},{:.1})mm\n\
         mesh               : tets={}  ({}×{}×{})\n\
         sweep              : {:.3}–{:.3} GHz, {} pts ({:.0} MHz step)\n\
         solve wall         : {:.1} s ({:.2} s/point)\n\
         peaks_resolvable   : {}\n\
         f_lo / f_hi (FEM)  : {:.4} / {:.4} GHz   (valley {:.2} dB)\n\
         k_fem  (split)     : {:.4} ({:.2}%)\n\
         k_eps  (analytic)  : {:.4} ({:.2}%)   <- physical resonant-split ref\n\
         k_imp  (analytic)  : {:.4} ({:.2}%)   <- coupled-line `coupling_coefficient`\n\
         VERDICT            : {}\n\
         reasoning          : {}\n\
         =================================",
        W * 1e3,
        S * 1e3,
        SUB_H * 1e3,
        EPS_R,
        geom.res_l * 1e3,
        geom.box_w * 1e3,
        geom.box_len * 1e3,
        geom.box_h * 1e3,
        geom.total_tets(),
        geom.nx,
        geom.ny,
        geom.nz,
        f_lo_hz / 1e9,
        f_hi_hz / 1e9,
        n_pts,
        (f_hi_hz - f_lo_hz) / 1e6 / (n_pts as f64 - 1.0),
        wall,
        wall / n_pts as f64,
        peaks_resolvable,
        f_lo_peak,
        f_hi_peak,
        valley_db,
        k_fem,
        k_fem * 100.0,
        k_eps,
        k_eps * 100.0,
        k_imp,
        k_imp * 100.0,
        verdict,
        reason,
    );

    // ---- Honest assertion: only that the pipeline ran (NOT a pass). -----
    // A NO-GO is a valid, valuable measurement (it sends F1.1b back to FDTD
    // or the maintainer). We assert ONLY non-degeneracy so a broken solve
    // (collapsed port / NaN sweep) surfaces loudly; we do NOT assert a verdict.
    let any_finite = curve_lin.iter().all(|(_, m)| m.is_finite());
    assert!(
        any_finite,
        "PROBE PIPELINE DEGENERATED: a non-finite |S21| in the sweep — the driven solve broke \
         (port collapsed or mesh invalid). Full spectrum printed above; this is a pipeline \
         failure, not a NO-GO measurement."
    );
}

// =====================================================================
// FAST UNIT CHECKS (no solve — run in the default `cargo test`)
// =====================================================================

#[cfg(test)]
mod unit {
    use super::*;

    /// The coupled-pair geometry is well-formed: 2 resonators + 2 feeds = 4
    /// rectangles, the box clears the trace pattern in x, sub_h lands on a
    /// z-plane, the feed gap is exactly one dy cell, and the feeds reach both
    /// port planes while the resonators stay clear of them.
    #[test]
    fn geometry_is_well_formed() {
        let geom = build_coupled_geometry();
        assert_eq!(geom.traces.len(), 4, "2 resonators + 2 feeds");

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

        // Feed gap is exactly one dy cell (so it is resolved, not sub-cell).
        let dy = geom.box_len / geom.ny as f64;
        assert!(
            (FEED_GAP / dy - 1.0).abs() < 0.05,
            "feed gap ({:.3}mm) should be ~1 dy cell ({:.3}mm) — got {:.2} cells",
            FEED_GAP * 1e3,
            dy * 1e3,
            FEED_GAP / dy,
        );

        // An input feed reaches y=0; an output feed reaches y=box_len.
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

        // The resonators do NOT touch the port planes (they are floating,
        // gap-coupled): no resonator rectangle starts at y=0 or ends at box_len.
        let res_a = &geom.traces[0];
        let res_b = &geom.traces[1];
        assert!(
            res_a.y0 > 1e-9 && (res_a.y0 + res_a.l) < geom.box_len - 1e-9,
            "resonator A must be floating (not touching either port plane)"
        );
        assert!(
            res_b.y0 > 1e-9 && (res_b.y0 + res_b.l) < geom.box_len - 1e-9,
            "resonator B must be floating (not touching either port plane)"
        );

        // The two resonators are gap-separated in x by S.
        let gap_x = res_b.x0 - (res_a.x0 + res_a.w);
        assert!(
            (gap_x - S).abs() < 1e-9,
            "resonator x-gap ({:.3}mm) must equal S ({:.3}mm)",
            gap_x * 1e3,
            S * 1e3
        );
    }

    /// The two analytic k definitions agree to within ~15% at this weak gap
    /// (S/W = 2), which is the premise that makes the k_split-vs-`coupling_
    /// coefficient` comparison meaningful. (Sanity-guards the gap choice; if a
    /// future edit tightens S this fires and flags that the two references
    /// have diverged.)
    #[test]
    fn analytic_k_definitions_agree_at_weak_gap() {
        let cm = coupled_microstrip(W, S, SUB_H, EPS_R);
        let k_imp = coupling_coefficient(&cm);
        let l = resonator_length();
        let f_e = C0 / (2.0 * l * cm.eps_eff_e.sqrt());
        let f_o = C0 / (2.0 * l * cm.eps_eff_o.sqrt());
        let (lo, hi) = (f_e.min(f_o), f_e.max(f_o));
        let k_eps = (hi * hi - lo * lo) / (hi * hi + lo * lo);
        let ratio = k_imp / k_eps;
        assert!(
            (1.0..=1.3).contains(&ratio),
            "at S/W=2 the coupled-line k_imp ({k_imp:.4}) and resonant-split k_eps ({k_eps:.4}) \
             should agree to ~15% (ratio {ratio:.2}); a larger ratio means the gap is too tight \
             for the two definitions to be compared like-for-like"
        );
    }
}
