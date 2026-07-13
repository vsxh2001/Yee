//! Gate `match-em-001` (FS.6.2b, ADR-0219): **matching synthesized from a
//! measured Γ improves the measured S11** — the FS.6 roadmap gate, closing
//! the loop measure → synthesize → regenerate layout → re-measure.
//!
//! Scenario: an edge-fed 2.45 GHz patch (deliberately mismatched — the
//! A.0 topology) on a long feed. Run 1 measures the complex reflection
//! Γ and the line's β at the reference plane A (a probe triple partway
//! down the feed, `sparams::fit_standing_wave` conventions). The
//! single-stub synthesis (`yee_layout::single_stub_match`, the
//! stub-match-001-certified construction) places a shunt open stub at
//! `x_A − d`; run 2 re-measures with the stub in place. The improvement
//! is judged at the port-side triple P, which sits generator-side of
//! every possible stub position by construction — |Γ| is plane-invariant
//! on the lossless feed, so P's before/after compare cleanly.
//!
//! `#[ignore]`'d (2 release FDTD solves, ~2.6 M cells each):
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_stub_match -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, Substrate, edge_fed_patch, single_stub_match};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const F0_HZ: f64 = 2.45e9;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const BW_HZ: f64 = 2.0e9;
const N_STEPS: usize = 9000;
/// Long feed: the synthesized stub lands up to λ_g/2 ≈ 33.6 mm behind
/// the reference plane A, which itself sits clear of the port triple.
const FEED_LEN_M: f64 = 58.0e-3;
/// Reference plane A: distance from the port along the feed.
const A0_OFFSET_M: f64 = 46.0e-3;
/// Port-side triple P: clear of the aperture's evanescent zone.
const P0_OFFSET_M: f64 = 3.0e-3;
/// Probe-triple spacing, cells.
const SP_CELLS: usize = 12;

fn fr4() -> Substrate {
    Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// The A.0 edge-fed patch with the feed stretched to `FEED_LEN_M`.
fn long_feed_patch() -> Layout {
    let mut dut = edge_fed_patch(F0_HZ, &fr4(), Z0_OHM);
    let w = dut.ports[0].width_m;
    dut.traces[0] = Polygon::rect(-FEED_LEN_M, -w / 2.0, FEED_LEN_M, w);
    dut.ports[0].at = Point2::new(-FEED_LEN_M, 0.0);
    dut.bbox = BBox::from_polygons(&dut.traces);
    dut
}

/// Add the synthesized shunt open stub (feed-width, jutting +y) at
/// layout-frame `x_stub`, overlapping the feed by one cell.
fn with_stub(mut dut: Layout, x_stub: f64, l_open: f64) -> Layout {
    let w = dut.ports[0].width_m;
    dut.traces.push(Polygon::rect(
        x_stub - w / 2.0,
        w / 2.0 - DX_M,
        w,
        l_open + DX_M,
    ));
    dut.bbox = BBox::from_polygons(&dut.traces);
    dut
}

struct Measured {
    /// Γ (re, im) at plane A and P, and the fitted β at f0, plus |Γ(f0)|
    /// magnitudes at both planes.
    gamma_a: (f64, f64),
    gamma_p: (f64, f64),
    beta_rad_m: f64,
    /// Layout-frame x of probe A0 (the Γ_a reference plane).
    x_a0: f64,
}

/// Voxelize, run, and wave-split one layout (probes: P triple then A
/// triple, both ordered along +x).
fn measure(layout: &Layout) -> Measured {
    let model = voxelize_microstrip(
        layout,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let (i_port, j_strip, k_top) = model.port_cells[0];
    let k_probe = k_top.saturating_sub(1).max(1);
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let port_x = layout.ports[0].at.x;
    let i_p0 = i_for(port_x + P0_OFFSET_M);
    let i_a0 = i_for(port_x + A0_OFFSET_M);
    assert!(i_a0 > i_p0 + 2 * SP_CELLS, "triples overlap");

    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < w_feed / 2.0 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);

    let materials = MaterialsSpec {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ex: model
            .grid
            .pec_mask_ex
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ey: model
            .grid
            .pec_mask_ey
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        ..MaterialsSpec::default()
    };
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * BW_HZ)) / dt).ceil() as usize;
    let mk_probe = |i: usize| ProbeSpec {
        component: "ez".into(),
        cell: (i, j_strip, k_probe),
    };
    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: DX_M,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![AperturePortSpec {
            i: i_port,
            j_lo,
            j_hi,
            k_lo: 0,
            k_top,
            resistance_ohm: Z0_OHM,
            v0: 1.0,
            f0_hz: F0_HZ,
            bw_hz: BW_HZ,
            t0_steps,
            record: false,
        }],
        probes: vec![
            mk_probe(i_p0),
            mk_probe(i_p0 + SP_CELLS),
            mk_probe(i_p0 + 2 * SP_CELLS),
            mk_probe(i_a0),
            mk_probe(i_a0 + SP_CELLS),
            mk_probe(i_a0 + 2 * SP_CELLS),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        spacings: None,
        backend: BackendChoice::Cpu,
    };
    let handle = yee_engine::submit(spec);
    let p = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result.probes),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");

    let spacing = SP_CELLS as f64 * dx;
    let split_at = |k0: usize| {
        let v: Vec<(f64, f64)> = (k0..k0 + 3)
            .map(|k| sparams::single_bin_dft(&p[k], dt, F0_HZ))
            .collect();
        sparams::fit_standing_wave(v[0], v[1], v[2], spacing)
    };
    let sp_p = split_at(0);
    let sp_a = split_at(3);
    let cdiv = |a: (f64, f64), b: (f64, f64)| {
        let n = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / n, (a.1 * b.0 - a.0 * b.1) / n)
    };
    Measured {
        gamma_a: cdiv(sp_a.bwd, sp_a.fwd),
        gamma_p: cdiv(sp_p.bwd, sp_p.fwd),
        beta_rad_m: sp_a.beta_rad_m,
        x_a0: x0 + i_a0 as f64 * dx,
    }
}

#[test]
#[ignore = "slow: 2 release FDTD solves (~2.6 M cells each); match-em-001 gate (FS.6.2b) — run with --release --ignored"]
fn synthesized_stub_match_improves_measured_s11() {
    let mag = |c: (f64, f64)| f64::hypot(c.0, c.1);
    let db = |x: f64| 20.0 * x.log10();

    // Run 1: the mismatched patch — measure Γ and β at plane A.
    let dut = long_feed_patch();
    eprintln!("match-em-001: unmatched run (feed {} mm)", FEED_LEN_M * 1e3);
    let m0 = measure(&dut);
    let g_a = mag(m0.gamma_a);
    let g_p0 = mag(m0.gamma_p);
    eprintln!(
        "  unmatched: |Γ_A| = {g_a:.4} ({:.2} dB), |Γ_P| = {g_p0:.4} ({:.2} dB), \
         β = {:.2} rad/m (HJ-expected ~{:.2})",
        db(g_a),
        db(g_p0),
        m0.beta_rad_m,
        2.0 * PI * F0_HZ * yee_layout::eps_eff(dut.ports[0].width_m, H_M, EPS_R).sqrt()
            / 299_792_458.0
    );
    assert!(
        g_a >= 0.35,
        "match-em-001: the edge-fed patch is not meaningfully mismatched \
         (|Γ| = {g_a:.3}) — matching it proves nothing"
    );

    // Synthesize the stub from the MEASURED Γ and β.
    let m = single_stub_match(m0.gamma_a, m0.beta_rad_m);
    let x_stub = m0.x_a0 - m.d_m;
    eprintln!(
        "  synthesis: d = {:.3} mm, l_open = {:.3} mm (b = {:+.3}) → stub at x = {:.3} mm",
        m.d_m * 1e3,
        m.l_open_m * 1e3,
        m.b,
        x_stub * 1e3
    );
    // The stub must land on the feed, generator-side of triple A and
    // load-side of triple P (both by fixture construction).
    let port_x = dut.ports[0].at.x;
    assert!(
        x_stub > port_x + P0_OFFSET_M + (2 * SP_CELLS) as f64 * DX_M + 2.0e-3,
        "stub {x_stub:.4} m collides with the port triple"
    );
    assert!(x_stub < m0.x_a0 - 1.0e-3, "stub collides with plane A");

    // Run 2: the matched layout.
    eprintln!("match-em-001: matched run");
    let m1 = measure(&with_stub(dut, x_stub, m.l_open_m));
    let g_p1 = mag(m1.gamma_p);
    eprintln!(
        "  matched:   |Γ_P| = {g_p1:.4} ({:.2} dB) — improvement {:.2} dB",
        db(g_p1),
        db(g_p0) - db(g_p1)
    );

    // The roadmap gate: the synthesized match improves the measured S11.
    // (Pinned from the first green run; the stub also radiates and
    // couples, so the TL-ideal null is not expected.)
    assert!(
        db(g_p1) <= db(g_p0) - 6.0,
        "match-em-001 FAILED: matched |S11| {:.2} dB vs unmatched {:.2} dB \
         (need ≥ 6 dB improvement)",
        db(g_p1),
        db(g_p0)
    );
}
