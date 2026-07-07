//! Gate `studio-antenna-e2e-001` (R.5c, ADR-0203): the antenna design
//! flow at the command layer — instant closed-form design with
//! byte-checked Gerber artifacts (the R.5 philosophy), plus the
//! designer-grade error path. The full-wave verify's physics is gated in
//! `yee-engine` (A.0–A.3); a reduced-fidelity pipe run keeps the studio
//! side honest.

use yee_studio_app::antenna::{AntennaDesignRequest, design_antenna_impl, verify_antenna_impl};

fn request() -> AntennaDesignRequest {
    serde_json::from_value(serde_json::json!({
        "f0_hz": 2.45e9,
        "eps_r": 4.4,
        "height_m": 1.6e-3,
    }))
    .expect("request deserializes")
}

#[test]
fn design_produces_balanis_dims_and_bytechecked_gerbers() {
    let resp = design_antenna_impl(&request()).expect("design failed");

    // Dims are the Balanis closed forms.
    let dims = yee_layout::patch_antenna_dims(2.45e9, 4.4, 1.6e-3);
    assert_eq!(resp.width_m, dims.width_m);
    assert_eq!(resp.length_m, dims.length_m);
    // Default inset is the A.3-measured 0.25·L, not the closed-form seed.
    assert!((resp.inset_m - 0.25 * dims.length_m).abs() < 1e-15);

    // Gerbers byte-identical to yee-export's own output for the same layout.
    let substrate = yee_layout::Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let layout =
        yee_layout::inset_fed_patch_with_depth(2.45e9, &substrate, 50.0, 0.25 * dims.length_m);
    assert_eq!(
        resp.gerber_copper,
        yee_export::layout_to_gerber(&layout, &yee_export::GerberOptions::default())
    );
    assert_eq!(
        resp.gerber_outline,
        yee_export::layout_to_gerber_outline(&layout, &yee_export::OutlineOptions::default())
    );
}

#[test]
fn unphysical_spec_reports_an_error() {
    let mut req = request();
    req.inset_frac = 0.6; // beyond the patch centre
    let err = design_antenna_impl(&req).expect_err("inset beyond L/2 must be rejected");
    assert!(err.contains("inset_frac"), "{err}");
}

#[test]
fn verify_pipe_streams_progress_and_returns_a_finite_curve() {
    // Reduced fidelity: coarse grid, short run — a pipe exercise (the
    // physics is gated by engine-antenna-001..004).
    let req = serde_json::from_value(serde_json::json!({
        "design": { "f0_hz": 2.45e9, "eps_r": 4.4, "height_m": 1.6e-3 },
        "dx_m": 0.9e-3,
        "n_steps": 700,
    }))
    .expect("verify request deserializes");
    let mut events = 0usize;
    let resp = verify_antenna_impl(&req, &mut |p| {
        assert_eq!(p.phase, "antenna");
        events += 1;
    })
    .expect("verify pipe failed");
    assert!(events > 0, "no progress streamed");
    assert_eq!(resp.freqs_hz.len(), resp.s11_db.len());
    assert!(resp.s11_db.iter().all(|v| v.is_finite()));
    assert!(
        resp.dip_db <= 0.0,
        "|S11| cannot exceed 0 dB: {}",
        resp.dip_db
    );
}
