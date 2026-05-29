//! `gerber-001` — RS-274X structural-validity gate.
//!
//! For a known small `Layout`, assert the emitted Gerber: begins with the
//! `%FSLAX46Y46*%` then `%MOMM*%` header; contains exactly one `G36*` and one
//! `G37*` per polygon (so the region count equals the layout polygon count);
//! defines at least one aperture (`%ADD…`); and ends with `M02*`.
//!
//! This is the I/O analogue of the Touchstone round-trip gate — Gerber is an
//! interchange format, so its gate is structural validity, not a physics
//! benchmark.

use yee_export::{GerberOptions, layout_to_gerber};
use yee_layout::{EdgeCoupledParams, EdgeCoupledSection, Substrate, edge_coupled_bpf};

/// A small, deterministic 2-section edge-coupled BPF on FR-4.
fn sample_layout() -> yee_layout::Layout {
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    let params = EdgeCoupledParams {
        substrate,
        sections: vec![
            EdgeCoupledSection {
                length_m: 10.0e-3,
                width_m: 1.0e-3,
                gap_m: 0.3e-3,
            },
            EdgeCoupledSection {
                length_m: 10.0e-3,
                width_m: 1.0e-3,
                gap_m: 0.3e-3,
            },
        ],
        feed_width_m: 3.0e-3,
        feed_length_m: 5.0e-3,
    };
    edge_coupled_bpf(&params)
}

#[test]
fn gerber_001_structure() {
    let layout = sample_layout();
    let n_polys = layout.traces.len();
    assert!(n_polys >= 1, "sample layout should have at least one trace");

    let gerber = layout_to_gerber(&layout, &GerberOptions::default());

    // Header order: %FSLAX46Y46*% then %MOMM*% at the very top.
    assert!(
        gerber.starts_with("%FSLAX46Y46*%\n%MOMM*%\n"),
        "header must begin with the format then units statements, got:\n{}",
        &gerber[..gerber.len().min(80)]
    );

    // At least one aperture definition.
    assert!(
        gerber.contains("%ADD"),
        "output must define at least one aperture (%ADD…)"
    );

    // Exactly one G36*/G37* region per polygon.
    let n_g36 = gerber.matches("G36*").count();
    let n_g37 = gerber.matches("G37*").count();
    assert_eq!(
        n_g36, n_polys,
        "expected one G36* per polygon ({n_polys}), found {n_g36}"
    );
    assert_eq!(
        n_g37, n_polys,
        "expected one G37* per polygon ({n_polys}), found {n_g37}"
    );

    // Ends with the end-of-file code.
    assert!(
        gerber.trim_end().ends_with("M02*"),
        "file must end with M02*, got tail:\n{}",
        &gerber[gerber.len().saturating_sub(40)..]
    );
}
