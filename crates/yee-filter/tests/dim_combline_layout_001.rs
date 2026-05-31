//! dim-combline-layout-001 (Filter Phase F1.2.6): combline board-layout gate.
//!
//! A NON-vacuous geometry gate for [`yee_filter::dimension_combline_layout`]
//! (ADR-0145), the board-layout companion of [`yee_filter::dimension_combline`].
//! It composes an honest **comb** — `N` aligned, short-circuited resonator lines
//! on a common ground spine, capacitively loaded at the open ends, with tapped
//! input/output feeds — and this gate proves the geometry actually consumes the
//! synthesized dimensions rather than emitting an empty or uniform placeholder.
//!
//! Build the demo order-5 combline layout (5-pole 0.1 dB Chebyshev band-pass,
//! FBW = 0.10, θ0 = 45° = π/4 — the same fixture as `dim_combline_001`) and
//! assert:
//!
//! 1. The `N = 5` resonator-line traces are present with dimensions matching
//!    [`yee_filter::dimension_combline`] (each `line_width_m` wide ×
//!    `resonator_length_m` long, tight tolerance).
//! 2. Consecutive resonator-line left edges (sorted) differ by
//!    `line_width_m + gaps_m[i]` — proving the layout uses the REAL solved
//!    per-section gaps, asserted against `dimension_combline`'s own `gaps_m`. The
//!    5-pole Chebyshev coupling matrix is symmetric (M₁₂=M₄₅, M₂₃=M₃₄), so the
//!    pitches mirror about the centre and are NOT all equal; a uniform-gap
//!    placeholder fails both the per-section match and the non-uniformity check.
//! 3. Exactly 2 ports, each `ref_impedance_ohm == spec z0`; bbox width + height
//!    positive + finite; every trace rect has positive extent.
//!
//! Patterned on `dim_stepped_001` / `dim_combline_001` (fixture) and the
//! geometry-gate idiom.

use std::f64::consts::FRAC_PI_4;

use yee_filter::{
    Approximation, FilterSpec, Response, SpecMask, dimension_combline, dimension_combline_layout,
    synthesize,
};
use yee_layout::{Polygon, Substrate};

/// FR-4 representative substrate (matches the combline / hairpin / stepped-Z gates).
fn substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// The demo 5-pole 0.1 dB Chebyshev band-pass at FBW = 0.10 (H&L §5.2.5 combline
/// example — the same fixture `dim_combline_001` uses).
fn spec_5pole_cheb_01db() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.1 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.1,
            return_loss_db: 16.0,
            stopband: vec![(2.4e9, 30.0)],
        },
    }
}

/// Axis-aligned rectangle extent `(x0, y0, w, h)` of a `Polygon::rect` polygon.
fn rect_extent(p: &Polygon) -> (f64, f64, f64, f64) {
    let xs = p.verts.iter().map(|v| v.x);
    let ys = p.verts.iter().map(|v| v.y);
    let x0 = xs.clone().fold(f64::INFINITY, f64::min);
    let x1 = xs.fold(f64::NEG_INFINITY, f64::max);
    let y0 = ys.clone().fold(f64::INFINITY, f64::min);
    let y1 = ys.fold(f64::NEG_INFINITY, f64::max);
    (x0, y0, x1 - x0, y1 - y0)
}

#[test]
fn dim_combline_layout_001_comb_geometry() {
    let spec = spec_5pole_cheb_01db();
    let z0 = spec.z0_ohm;
    let proj = synthesize(&spec);
    let sub = substrate();

    // The dimensions the layout must consume (no recompute in the layout).
    let dims = dimension_combline(&proj, FRAC_PI_4, &sub)
        .expect("N=5 coupled-resonator combline fixture should dimension without error");
    let n = dims.gaps_m.len() + 1;
    assert_eq!(n, 5, "order-5 combline → N = 5 resonators, 4 gaps");

    let w = dims.line_width_m;
    let l = dims.resonator_length_m;

    let layout = dimension_combline_layout(&proj, FRAC_PI_4, &sub)
        .expect("combline layout should compose without error");

    // --- (3, part) substrate carried through, all rects positive extent -----
    assert_eq!(layout.substrate, sub, "layout carries the input substrate");
    for (i, t) in layout.traces.iter().enumerate() {
        let (_, _, tw, th) = rect_extent(t);
        assert!(
            tw.is_finite() && th.is_finite() && tw > 0.0 && th > 0.0,
            "trace {i} must have positive finite extent, got w={tw:.6e} h={th:.6e}"
        );
    }

    // --- (1) N resonator-line traces with dims == dimension_combline --------
    // Identify resonator lines by height == resonator_length_m (the spine is
    // w-tall, cap pads are w×w squares, feeds are w-tall horizontal bars — none
    // share the resonator's l height, so this is an unambiguous discriminator).
    let tol = 1e-9_f64.max(l * 1e-6);
    let mut resonators: Vec<(f64, f64, f64, f64)> = layout
        .traces
        .iter()
        .map(rect_extent)
        .filter(|&(_, _, tw, th)| (th - l).abs() <= tol && (tw - w).abs() <= tol)
        .collect();
    assert_eq!(
        resonators.len(),
        n,
        "expected {n} resonator-line traces (line_width × resonator_length), found {}",
        resonators.len()
    );
    for (i, &(x0, y0, tw, th)) in resonators.iter().enumerate() {
        assert!(
            (tw - w).abs() <= tol,
            "resonator {i} width {tw:.6e} != line_width_m {w:.6e}"
        );
        assert!(
            (th - l).abs() <= tol,
            "resonator {i} height {th:.6e} != resonator_length_m {l:.6e}"
        );
        // Short-circuit end at y = 0 (the comb's grounded end).
        assert!(
            y0.abs() <= tol,
            "resonator {i} short-circuit end y0 = {y0:.6e} should be at y = 0"
        );
        let _ = x0;
    }

    // --- (2) consecutive resonator x-pitch == w + gaps_m[i] (monotone) ------
    // Sort resonator left edges and check each pitch against the SOLVED gap. A
    // uniform-gap placeholder (constant pitch) fails because gaps_m are distinct.
    resonators.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let left_edges: Vec<f64> = resonators.iter().map(|&(x0, ..)| x0).collect();
    let pitch_tol = 1e-9_f64.max(w * 1e-6);
    let mut pitches = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let pitch = left_edges[i + 1] - left_edges[i];
        let expected = w + dims.gaps_m[i];
        pitches.push(pitch);
        assert!(
            (pitch - expected).abs() <= pitch_tol,
            "resonator x-pitch[{i}] = {pitch:.9e} m != line_width_m + gaps_m[{i}] = \
             {expected:.9e} m (w={w:.6e}, gap={:.6e})",
            dims.gaps_m[i]
        );
    }
    // Non-uniform placement: the solved gaps are NOT all equal, so the pitches
    // must differ across sections — a uniform-gap placeholder (constant pitch)
    // fails here. (A 5-pole Chebyshev coupling matrix is symmetric — M₁₂=M₄₅ and
    // M₂₃=M₃₄ — so the gaps mirror, gaps_m = [g0, g1, g1, g0]; they are NOT
    // monotone, but the inner pair differs from the outer pair, which is the
    // real, non-uniform per-section structure we must reproduce.)
    let pitch0 = pitches[0];
    assert!(
        pitches.iter().any(|&p| (p - pitch0).abs() > pitch_tol),
        "x-pitches must not all be equal (a uniform-gap placeholder fails): {pitches:?}"
    );
    // The coupling matrix is symmetric, so the gaps — and hence the pitches —
    // mirror about the centre: pitch[i] == pitch[N-2-i]. (For this 5-pole
    // Chebyshev the inner pair couples more weakly than the outer pair, so
    // pitch[1] > pitch[0]; we assert only the symmetry + the non-uniformity, not
    // a hardcoded inner/outer ordering.)
    for i in 0..pitches.len() {
        let j = pitches.len() - 1 - i;
        assert!(
            (pitches[i] - pitches[j]).abs() <= pitch_tol,
            "symmetric coupling → mirrored pitches: pitch[{i}]={:.9e} != pitch[{j}]={:.9e}",
            pitches[i],
            pitches[j]
        );
    }

    // --- (3) exactly 2 ports at z0, positive + finite bbox ------------------
    assert_eq!(layout.ports.len(), 2, "combline layout has exactly 2 ports");
    for (i, p) in layout.ports.iter().enumerate() {
        assert!(
            (p.ref_impedance_ohm - z0).abs() < 1e-12,
            "port {i} ref impedance {} != spec z0 {z0}",
            p.ref_impedance_ohm
        );
        assert!(
            p.width_m.is_finite() && p.width_m > 0.0,
            "port {i} width {:.6e} must be finite and > 0",
            p.width_m
        );
    }

    let bw = layout.bbox.width();
    let bh = layout.bbox.height();
    assert!(
        bw.is_finite() && bh.is_finite() && bw > 0.0 && bh > 0.0,
        "bbox extent must be positive + finite: w={bw:.6e} h={bh:.6e}"
    );

    // Surface the geometry for the verification log.
    println!(
        "combline layout: N = {n} resonator lines, line_width = {w:.6e} m, \
         resonator_length = {l:.6e} m, bbox = {bw:.6e} × {bh:.6e} m"
    );
    println!("  solved gaps_m   = {:?}", dims.gaps_m);
    println!("  measured pitches= {pitches:?}");
    println!(
        "  expected pitches= {:?}",
        dims.gaps_m.iter().map(|g| w + g).collect::<Vec<_>>()
    );
}
