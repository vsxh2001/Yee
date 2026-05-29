//! `gerber-003` — board-outline (Edge.Cuts) structural-validity gate.
//!
//! For a known small `Layout`, assert the emitted outline Gerber: begins with
//! the `%FSLAX46Y46*%` then `%MOMM*%` header; defines exactly one aperture
//! (`%ADD…`); has exactly one `D02*` pen-up move and at least four `D01*`
//! pen-down draws (a closed rectangle: three sides + the explicit close);
//! contains NO `G36*`/`G37*` (the outline is a stroked contour, not a region
//! fill); and ends with `M02*`.
//!
//! Like `gerber-001`, this is the I/O-structural analogue of a physics gate —
//! Gerber is an interchange format, so its gate is structural validity.

use yee_export::{OutlineOptions, layout_to_gerber_outline};
use yee_layout::{BBox, Layout, Polygon, Substrate};

/// A small, deterministic single-rectangle layout with a known bbox.
fn sample_layout() -> Layout {
    let traces = vec![Polygon::rect(2.0e-3, 1.0e-3, 10.0e-3, 5.0e-3)];
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    Layout {
        substrate,
        bbox: BBox::from_polygons(&traces),
        traces,
        ports: vec![],
    }
}

#[test]
fn gerber_003_outline_structure() {
    let layout = sample_layout();
    let gerber = layout_to_gerber_outline(&layout, &OutlineOptions::default());

    // Header order: %FSLAX46Y46*% then %MOMM*% at the very top.
    assert!(
        gerber.starts_with("%FSLAX46Y46*%\n%MOMM*%\n"),
        "header must begin with the format then units statements, got:\n{}",
        &gerber[..gerber.len().min(80)]
    );

    // Exactly one aperture definition.
    let n_add = gerber.matches("%ADD").count();
    assert_eq!(
        n_add, 1,
        "outline must define exactly one aperture (%ADD…), found {n_add}"
    );

    // Exactly one D02* pen-up move (to the first corner).
    let n_d02 = gerber.matches("D02*").count();
    assert_eq!(n_d02, 1, "expected exactly one D02* move, found {n_d02}");

    // At least four D01* draws (three sides + explicit close).
    let n_d01 = gerber.matches("D01*").count();
    assert!(
        n_d01 >= 4,
        "expected at least four D01* draws (closed rectangle), found {n_d01}"
    );

    // The outline is stroked, NOT region-filled: no G36*/G37*.
    assert_eq!(
        gerber.matches("G36*").count(),
        0,
        "outline must not contain any G36* (it is stroked, not filled)"
    );
    assert_eq!(
        gerber.matches("G37*").count(),
        0,
        "outline must not contain any G37* (it is stroked, not filled)"
    );

    // Ends with the end-of-file code.
    assert!(
        gerber.trim_end().ends_with("M02*"),
        "file must end with M02*, got tail:\n{}",
        &gerber[gerber.len().saturating_sub(40)..]
    );
}
