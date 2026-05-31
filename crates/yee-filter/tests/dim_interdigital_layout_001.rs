//! dim-interdigital-layout-001 (Filter Phase F1.2.8): interdigital board-layout gate.
//!
//! A NON-vacuous geometry gate for [`yee_filter::dimension_interdigital_layout`]
//! (ADR-0149), the board-layout companion of [`yee_filter::dimension_interdigital`].
//! It composes an honest interdigital **comb** — `N` aligned, short-circuited
//! `λ_g/4` resonator lines grounded at ALTERNATING ends between TWO ground rails,
//! with NO loading-cap pads and tapped input/output feeds — and this gate proves
//! the geometry actually consumes the synthesized dimensions AND realizes the
//! three interdigital-distinct features that set it apart from combline:
//!
//!   - TWO ground rails (combline has one common spine),
//!   - NO cap pads (trace count `N + 4`, not combline's `2N + 3`),
//!   - alternating even/odd `y`-origin offset, with no resonator touching both
//!     rails (no accidental short → cavity).
//!
//! A combline-style single-spine / with-pads layout would FAIL parts 2, 3 and 4.
//!
//! Build the demo order-5 interdigital layout (5-pole 0.1 dB Chebyshev band-pass,
//! FBW = 0.10 — the same fixture as `dim_combline_layout_001`, minus the θ0
//! parameter since interdigital is θ = π/2 fixed) on the FR-4 substrate and assert
//! the spec DoD parts 1–6.
//!
//! Patterned on `dim_combline_layout_001` (fixture + geometry-gate idiom).

use yee_filter::{
    Approximation, FilterSpec, Response, SpecMask, dimension_interdigital,
    dimension_interdigital_layout, synthesize,
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

/// The demo 5-pole 0.1 dB Chebyshev band-pass at FBW = 0.10 (same fixture the
/// combline-layout gate uses).
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
fn dim_interdigital_layout_001_comb_geometry() {
    let spec = spec_5pole_cheb_01db();
    let z0 = spec.z0_ohm;
    let proj = synthesize(&spec);
    let sub = substrate();

    // The dimensions the layout must consume (no recompute in the layout).
    let dims = dimension_interdigital(&proj, &sub)
        .expect("N=5 coupled-resonator interdigital fixture should dimension without error");
    let n = dims.gaps_m.len() + 1;
    assert_eq!(n, 5, "order-5 interdigital → N = 5 resonators, 4 gaps");

    let w = dims.line_width_m;
    let l = dims.resonator_length_m;
    // The open-end coupling gap is the neutral fixed default `g_open = w`
    // (ADR-0149); the layout offsets odd resonators up by this amount.
    let g_open = w;

    let layout = dimension_interdigital_layout(&proj, &sub)
        .expect("interdigital layout should compose without error");

    // --- all rects positive extent, substrate carried through ---------------
    assert_eq!(layout.substrate, sub, "layout carries the input substrate");
    for (i, t) in layout.traces.iter().enumerate() {
        let (_, _, tw, th) = rect_extent(t);
        assert!(
            tw.is_finite() && th.is_finite() && tw > 0.0 && th > 0.0,
            "trace {i} must have positive finite extent, got w={tw:.6e} h={th:.6e}"
        );
    }

    // === DoD (3) NO cap pads: trace count == N + 4, NOT combline's 2N + 3 ====
    // N resonator lines + 2 ground rails + 2 feeds. A combline-style layout
    // (1 spine + N w×w cap pads) would have 2N + 3 = 13 traces and FAIL here.
    let expected_traces = n + 4;
    assert_eq!(
        layout.traces.len(),
        expected_traces,
        "interdigital trace count must be N + 4 = {expected_traces} (N lines + 2 rails + 2 \
         feeds, NO cap pads); a combline-style layout (2N+3 = {}) would FAIL here",
        2 * n + 3
    );
    println!(
        "trace-count: {} == N + 4 (= {} lines + 2 rails + 2 feeds, NO cap pads) [combline 2N+3 \
         = {} would FAIL]",
        layout.traces.len(),
        n,
        2 * n + 3
    );

    let tol = 1e-9_f64.max(l * 1e-6);
    let w_tol = 1e-9_f64.max(w * 1e-6);

    // === DoD (2) TWO ground rails (the alternating-ground structure) =========
    // A rail is a horizontal w-tall bar spanning the full comb x-range
    // [0, comb_right]; one below the comb (y < 0) and one above (y > l). Combline
    // has exactly ONE. Identify rails by width >> w (they span all N lines).
    let comb_right = {
        // Reconstruct comb_right from the resonator placement (= last left edge + w).
        let mut x = 0.0_f64;
        for i in 0..n {
            if i == n - 1 {
                break;
            }
            x += w + dims.gaps_m[i];
        }
        x + w
    };
    let rails: Vec<(f64, f64, f64, f64)> = layout
        .traces
        .iter()
        .map(rect_extent)
        .filter(|&(x0, _, tw, th)| {
            (th - w).abs() <= w_tol && x0.abs() <= tol && (tw - comb_right).abs() <= tol
        })
        .collect();
    assert_eq!(
        rails.len(),
        2,
        "interdigital must have exactly TWO ground rails spanning [0, comb_right]; combline has \
         ONE. found {} rail(s)",
        rails.len()
    );
    // One rail below the comb (y < 0), one above (y > l).
    let below: Vec<_> = rails.iter().filter(|&&(_, y0, ..)| y0 < -tol).collect();
    let above: Vec<_> = rails.iter().filter(|&&(_, y0, ..)| y0 > l - tol).collect();
    assert_eq!(
        below.len(),
        1,
        "exactly one ground rail below the comb (y < 0)"
    );
    assert_eq!(
        above.len(),
        1,
        "exactly one ground rail above the comb (y > l)"
    );
    let bottom_rail = below[0];
    let top_rail = above[0];
    // Bottom rail at y ∈ [−w, 0]; top rail at y ∈ [l + g_open, l + g_open + w].
    assert!(
        (bottom_rail.1 - (-w)).abs() <= tol,
        "bottom rail y0 = {:.6e} should be at −w = {:.6e}",
        bottom_rail.1,
        -w
    );
    assert!(
        (top_rail.1 - (l + g_open)).abs() <= tol,
        "top rail y0 = {:.6e} should be at l + g_open = {:.6e}",
        top_rail.1,
        l + g_open
    );
    println!(
        "two-rail: bottom y ∈ [{:.6e}, {:.6e}], top y ∈ [{:.6e}, {:.6e}], both span x ∈ \
         [0, {:.6e}] (combline has ONE spine → would FAIL)",
        bottom_rail.1,
        bottom_rail.1 + bottom_rail.3,
        top_rail.1,
        top_rail.1 + top_rail.3,
        comb_right
    );

    // === DoD (1) N resonator lines, dims from the engine (no recompute drift) =
    // Resonator lines have height == resonator_length_m (= λ_g/4) and width == w.
    // The rails are w-tall and span comb_right wide, the feeds are w-tall — none
    // share the resonator's l height, so this is an unambiguous discriminator.
    let mut resonators: Vec<(f64, f64, f64, f64)> = layout
        .traces
        .iter()
        .map(rect_extent)
        .filter(|&(_, _, tw, th)| (th - l).abs() <= tol && (tw - w).abs() <= w_tol)
        .collect();
    assert_eq!(
        resonators.len(),
        n,
        "expected {n} resonator-line traces (line_width × resonator_length), found {}",
        resonators.len()
    );
    for (i, &(_, _, tw, th)) in resonators.iter().enumerate() {
        assert!(
            (tw - w).abs() <= w_tol,
            "resonator {i} width {tw:.6e} != line_width_m {w:.6e}"
        );
        assert!(
            (th - l).abs() <= tol,
            "resonator {i} height {th:.6e} != resonator_length_m {l:.6e} (λ_g/4)"
        );
    }

    // === DoD (4) Alternating connectivity (even/odd y-origin offset) =========
    // Sort by left edge so index = position in the comb. Even-index resonators
    // sit at y0 = 0 (grounded bottom, open top gapped g_open below the top rail);
    // odd-index resonators sit at y0 = g_open (grounded top, open bottom gapped
    // g_open above the bottom rail). NO resonator may touch BOTH rails. A combline
    // layout (all y0 = 0) would FAIL the odd-index check.
    resonators.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (i, &(_, y0, _, th)) in resonators.iter().enumerate() {
        let expected_y0 = if i % 2 == 0 { 0.0 } else { g_open };
        assert!(
            (y0 - expected_y0).abs() <= tol,
            "resonator {i} y-origin {y0:.6e} != expected {expected_y0:.6e} (even → 0 grounded \
             bottom, odd → g_open grounded top); a combline all-y0=0 layout FAILS here"
        );
        // No resonator touches BOTH rails (no accidental short → cavity).
        let bottom_edge = y0;
        let top_edge = y0 + th;
        let touches_bottom = bottom_edge.abs() <= tol; // shares bottom rail's y=0 edge
        let touches_top = (top_edge - (l + g_open)).abs() <= tol; // shares top rail's y=l+g_open edge
        assert!(
            !(touches_bottom && touches_top),
            "resonator {i} touches BOTH rails (bottom_edge={bottom_edge:.6e}, \
             top_edge={top_edge:.6e}) — that is an accidental short → cavity; the offset is wrong"
        );
        // It must touch exactly one rail (the grounded end).
        assert!(
            touches_bottom ^ touches_top,
            "resonator {i} must be grounded at exactly one rail (bottom={touches_bottom}, \
             top={touches_top})"
        );
    }
    println!(
        "alternating-offset: y-origins = {:?} (even→0, odd→g_open={:.6e}); no resonator touches \
         both rails (no short)",
        resonators.iter().map(|&(_, y0, ..)| y0).collect::<Vec<_>>(),
        g_open
    );

    // === DoD (5) Solved per-section pitch + symmetry =========================
    // Centre-to-centre pitch i→i+1 == w + gaps_m[i] (the real solved gaps), and
    // symmetric (gaps_m palindrome) for the symmetric Chebyshev coupling. A
    // uniform-gap placeholder (constant pitch) fails the non-uniformity check.
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
    let pitch0 = pitches[0];
    assert!(
        pitches.iter().any(|&p| (p - pitch0).abs() > pitch_tol),
        "x-pitches must not all be equal (a uniform-gap placeholder fails): {pitches:?}"
    );
    // Symmetric coupling → mirrored pitches: pitch[i] == pitch[N-2-i] (gaps_m a
    // palindrome), exactly as the combline gate asserts.
    for i in 0..pitches.len() {
        let j = pitches.len() - 1 - i;
        assert!(
            (pitches[i] - pitches[j]).abs() <= pitch_tol,
            "symmetric coupling → mirrored pitches: pitch[{i}]={:.9e} != pitch[{j}]={:.9e}",
            pitches[i],
            pitches[j]
        );
    }

    // === DoD (6) Two Z0-referenced ports =====================================
    assert_eq!(
        layout.ports.len(),
        2,
        "interdigital layout has exactly 2 ports"
    );
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
        "interdigital layout: N = {n} resonator lines (λ_g/4), line_width = {w:.6e} m, \
         resonator_length = {l:.6e} m, bbox = {bw:.6e} × {bh:.6e} m"
    );
    println!("  solved gaps_m   = {:?}", dims.gaps_m);
    println!("  measured pitches= {pitches:?}");
}
