//! lumped-pcb-001 (Filter Phase F2.2): the LC ladder places onto a valid board.
//!
//! Synthesize the committed Chebyshev 0.5 dB N=5 BPF (f0 = 2 GHz, FBW = 0.10,
//! Z0 = 50 Ω — the `lumped_001` fixture), realize it as a lumped LC ladder via
//! [`synthesize_lumped`], place it on a board with [`lumped_board`] (0603
//! footprints on FR-4), and assert the placement is geometrically valid:
//!
//! 1. exactly `2·N = 10` placements (an L + a C footprint per resonator), all
//!    ref-des unique;
//! 2. the layout carries at least `2` copper pads per placement (each footprint
//!    is two pads) plus the signal line + ground rail;
//! 3. **no pad overlaps** — every pair of pad rectangles is axis-aligned-disjoint
//!    (the placement spacing is valid);
//! 4. the bounding box is finite, positive-area, and contains every pad;
//! 5. series footprints sit on the signal line; shunt footprints drop toward the
//!    ground rail at `y = 0` (asserted by branch).
//!
//! Pure-geometry, deterministic, NO FDTD. This is the published-benchmark gate
//! for the lumped-LC → PCB placement; do NOT weaken it.

use yee_filter::{
    Approximation, BranchKind, FilterSpec, Footprint, Response, SpecMask, lumped_board, synthesize,
    synthesize_lumped,
};
use yee_layout::{Point2, Polygon, Substrate};

/// Chebyshev 0.5 dB N=5 bandpass spec (clone of the lumped-001 fixture).
fn fixture() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 40.0)],
        },
    }
}

/// FR-4 substrate (εr 4.4, h 1.6 mm) — the project's reference board.
fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// Axis-aligned rectangle `(min, max)` in metres.
#[derive(Clone, Copy, Debug)]
struct Rect {
    min: Point2,
    max: Point2,
}

/// Extract the axis-aligned bounding rect of a polygon (our pads/traces are all
/// axis-aligned rects, so this is exact).
fn rect_of(poly: &Polygon) -> Rect {
    let mut min = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for v in &poly.verts {
        min.x = min.x.min(v.x);
        min.y = min.y.min(v.y);
        max.x = max.x.max(v.x);
        max.y = max.y.max(v.y);
    }
    Rect { min, max }
}

/// Do two axis-aligned rects overlap with positive area? Edge-touching
/// (shared boundary, zero overlap area) is **not** an overlap.
fn rects_overlap(a: &Rect, b: &Rect) -> bool {
    let x_overlap = (a.min.x.max(b.min.x)) < (a.max.x.min(b.max.x)) - 1e-15;
    let y_overlap = (a.min.y.max(b.min.y)) < (a.max.y.min(b.max.y)) - 1e-15;
    x_overlap && y_overlap
}

#[test]
fn lumped_pcb_001() {
    let spec = fixture();
    let proj = synthesize(&spec);
    let n = proj.prototype.order();
    assert_eq!(n, 5, "fixture is order N=5");

    let ladder = synthesize_lumped(&proj).expect("N=5 bandpass fixture should synthesize");
    let board = lumped_board(&ladder, &fr4(), Footprint::Smd0603);

    // --- (1) 2·N placements, unique ref-des -------------------------------
    assert_eq!(
        board.placements.len(),
        2 * n,
        "expected an L + a C footprint per resonator (2·N = {})",
        2 * n
    );
    let mut ref_des: Vec<&str> = board
        .placements
        .iter()
        .map(|p| p.ref_des.as_str())
        .collect();
    let unique_count = {
        let mut v = ref_des.clone();
        v.sort_unstable();
        v.dedup();
        v.len()
    };
    assert_eq!(unique_count, ref_des.len(), "ref-des must be unique");
    // Ladder order: L1, C1, L2, C2, …
    let expected: Vec<String> = (1..=n)
        .flat_map(|k| [format!("L{k}"), format!("C{k}")])
        .collect();
    ref_des.clear();
    assert_eq!(
        board
            .placements
            .iter()
            .map(|p| p.ref_des.clone())
            .collect::<Vec<_>>(),
        expected,
        "ref-des should be L1,C1,L2,C2,… in ladder order"
    );

    // --- (2) pad count -----------------------------------------------------
    // Each footprint is two pads; the signal line + ground rail add more rects.
    assert!(
        board.layout.traces.len() >= 2 * board.placements.len(),
        "expected >= 2 pads per placement ({} traces, {} placements)",
        board.layout.traces.len(),
        board.placements.len()
    );

    // --- (3) no pad overlap ------------------------------------------------
    // Take every trace's axis-aligned rect (pads AND line/rail are all rects)
    // and assert every distinct pair is disjoint in positive area. The signal
    // line abuts the series pads and the ground-rail edge, but only touches
    // them (shared boundary, zero overlap area), so the strict-overlap test
    // passes. This is the placement-validity gate.
    let rects: Vec<Rect> = board.layout.traces.iter().map(rect_of).collect();
    for i in 0..rects.len() {
        for j in (i + 1)..rects.len() {
            assert!(
                !rects_overlap(&rects[i], &rects[j]),
                "trace rects {i} and {j} overlap: {:?} vs {:?}",
                rects[i],
                rects[j]
            );
        }
    }

    // --- (4) finite, positive-area bbox containing all pads ---------------
    let bb = board.layout.bbox;
    assert!(
        bb.min.x.is_finite()
            && bb.min.y.is_finite()
            && bb.max.x.is_finite()
            && bb.max.y.is_finite(),
        "bbox must be finite: {bb:?}"
    );
    assert!(
        bb.width() > 0.0 && bb.height() > 0.0,
        "bbox must have positive area: {bb:?}"
    );
    for (i, r) in rects.iter().enumerate() {
        assert!(
            r.min.x >= bb.min.x - 1e-12
                && r.max.x <= bb.max.x + 1e-12
                && r.min.y >= bb.min.y - 1e-12
                && r.max.y <= bb.max.y + 1e-12,
            "trace {i} {r:?} not contained in bbox {bb:?}"
        );
    }

    // --- (5) series on the line, shunt toward ground ---------------------
    let series_y: Vec<f64> = board
        .placements
        .iter()
        .filter(|p| p.kind == BranchKind::Series)
        .map(|p| p.center_m.1)
        .collect();
    let shunt_y: Vec<f64> = board
        .placements
        .iter()
        .filter(|p| p.kind == BranchKind::Shunt)
        .map(|p| p.center_m.1)
        .collect();
    // Shunt-first N=5 → resonators Shunt,Series,Shunt,Series,Shunt.
    assert_eq!(
        series_y.len(),
        4,
        "N=5 shunt-first → 2 series resonators × 2 footprints"
    );
    assert_eq!(
        shunt_y.len(),
        6,
        "N=5 shunt-first → 3 shunt resonators × 2 footprints"
    );
    // Every series footprint sits strictly above every shunt footprint, and
    // every shunt footprint sits strictly between the ground rail (y≈0) and the
    // signal line.
    for &sy in &series_y {
        assert!(sy > 0.0, "series footprint y must be positive: {sy}");
        for &hy in &shunt_y {
            assert!(
                hy < sy,
                "shunt footprint {hy} must drop below series footprint {sy}"
            );
            assert!(hy > 0.0, "shunt footprint y must be above the rail: {hy}");
        }
    }
}
