//! Gate `studio-design-e2e-001` (R.5, ADR-0198): the studio's scripted
//! design flow, driven headlessly at the command layer — the same
//! `design_filter_impl` the `design_filter` Tauri command wraps — with the
//! export artifacts **byte-checked** against the underlying validated
//! libraries:
//!
//! 1. the `.s2p` artifact round-trips through `yee_io::touchstone` (read →
//!    re-render → byte-identical), so the studio always emits a valid,
//!    passive Touchstone file;
//! 2. the Gerber artifacts are byte-identical to `yee_export`'s output for
//!    the same synthesized layout — the studio adds no bytes of its own;
//! 3. the design response is the validated coupling-matrix response
//!    (passband ≈ 0 dB at f0, stopband rejected, S11 dips in-band).

use yee_studio_app::design::{FilterDesignRequest, design_filter_impl};

fn request() -> FilterDesignRequest {
    // The R.4 gate's stack: h = 0.8 mm FR-4, where the qe→tap is realizable.
    serde_json::from_value(serde_json::json!({
        "f0_hz": 5.0e9,
        "fbw": 0.22,
        "order": 3,
        "eps_r": 4.4,
        "height_m": 0.8e-3,
    }))
    .expect("request deserializes")
}

#[test]
fn scripted_design_flow_produces_bytechecked_artifacts() {
    let resp = design_filter_impl(&request()).expect("design flow failed");

    // ---- 1. Touchstone artifact: valid + byte-stable round-trip ----
    let dir = std::env::temp_dir().join("studio_design_e2e");
    std::fs::create_dir_all(&dir).unwrap();
    let s2p_path = dir.join("design.s2p");
    std::fs::write(&s2p_path, &resp.s2p).unwrap();
    let back = yee_io::touchstone::read(&s2p_path).expect("studio .s2p failed to read back");
    assert_eq!(back.n_ports, 2);
    let re_rendered = yee_io::touchstone::to_string(&back).expect("re-render failed");
    assert_eq!(
        resp.s2p, re_rendered,
        "studio .s2p is not byte-stable through the touchstone reader"
    );

    // ---- 2. Gerber artifacts: byte-identical to yee-export's own output ----
    let spec = yee_filter::FilterSpec {
        response: yee_filter::Response::Bandpass,
        approximation: yee_filter::Approximation::Butterworth,
        f0_hz: 5.0e9,
        fbw: 0.22,
        order: Some(3),
        z0_ohm: 50.0,
        mask: yee_filter::SpecMask {
            passband_ripple_db: 3.0,
            return_loss_db: 10.0,
            stopband: vec![],
        },
    };
    let substrate = yee_layout::Substrate {
        eps_r: 4.4,
        height_m: 0.8e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let project = yee_filter::synthesize(&spec);
    let dims = yee_filter::dimension_hairpin_with_fold(&project, &substrate, 2.0).unwrap();
    let layout = yee_layout::hairpin_bpf_sections(&yee_layout::HairpinSectionParams {
        substrate,
        arm_length_m: dims.arm_length_m,
        line_width_m: dims.line_width_m,
        fold_spacing_m: dims.fold_spacing_m,
        gaps_m: dims.gaps_m.clone(),
        tap_offset_m: dims.tap_offset_m,
        feed_width_m: dims.line_width_m,
        feed_length_m: dims.arm_length_m,
    });
    let copper = yee_export::layout_to_gerber(&layout, &yee_export::GerberOptions::default());
    let outline =
        yee_export::layout_to_gerber_outline(&layout, &yee_export::OutlineOptions::default());
    assert_eq!(
        resp.gerber_copper, copper,
        "studio copper Gerber diverged from yee-export's output"
    );
    assert_eq!(
        resp.gerber_outline, outline,
        "studio outline Gerber diverged from yee-export's output"
    );

    // ---- 3. The design response is a real band-pass design ----
    let f0 = 5.0e9;
    let at = |f: f64| -> usize {
        resp.freqs_hz
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| ((**a) - f).abs().partial_cmp(&((**b) - f).abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap()
    };
    let i0 = at(f0);
    let i_stop = at(f0 * (1.0 + 1.8 * 0.22));
    assert!(
        resp.s21_db[i0] > -1.0,
        "design |S21|(f0) = {} dB — passband missing",
        resp.s21_db[i0]
    );
    assert!(
        resp.s21_db[i_stop] < -15.0,
        "design |S21| at the stopband edge = {} dB — no rejection",
        resp.s21_db[i_stop]
    );
    assert!(
        resp.s11_db[i0] < -6.0,
        "design |S11|(f0) = {} dB — no in-band match",
        resp.s11_db[i0]
    );
    // Dimensions surfaced for display are the synthesized ones.
    assert_eq!(resp.gaps_m.len(), 2);
    assert!((resp.tap_offset_m - dims.tap_offset_m).abs() < 1e-15);
}

#[test]
fn unrealizable_spec_reports_a_designer_grade_error() {
    // A narrow FBW on the thick 1.6 mm board with a wide fold: qe = 1/FBW
    // = 10 puts the tap beyond the fold-shortened arm — the flow must
    // surface the dims' TapNotRealizable message (with the realizable qe
    // range), not panic. (The R.6 corner correction lengthens arms, so the
    // pre-R.6 trigger — fbw 0.10 at the default fold — became realizable;
    // fold_widths 3.5 re-creates the wall.)
    let mut req = request();
    req.height_m = 1.6e-3;
    req.fbw = 0.10;
    req.fold_widths = 3.5;
    let err = design_filter_impl(&req).expect_err("narrow-FBW wide-fold stack should be rejected");
    assert!(
        err.contains("realizable range"),
        "unexpected error text: {err}"
    );
}
