//! Gate `engine-import-twin-001` (FS.3.2c, ADR-0229): the Gerber import
//! chain — export → bytes → import → [`Layout`] — reproduces a natively
//! built board closely enough that the SAME full-wave measurement lands
//! on the same physics. Round-trip byte-identity (`gerber-rt-001`) only
//! proves the file layer; this gate proves the chain end to end:
//! Gerber file → import → outline → `Layout` → voxelize → measurement.
//!
//! **Native twin**: the S.6 `sparams_stub_notch` scenario (FR-4
//! microstrip line + a Hammerstad-corrected λ/4 open stub — a textbook
//! bandstop notch) rebuilt here as the trace geometry. It is exported to
//! Gerber bytes with [`yee_export::layout_to_gerber`], reimported with
//! [`yee_export::gerber_to_layout`] — the exact helper the studio's
//! `import_gerber` command already uses
//! (`studio/src-tauri/src/import.rs`; Gerber carries no stackup/ports, so
//! the caller supplies both, same as the studio's import dialog) — and
//! measured through the R.5b [`yee_engine::board::two_port_board_job`]
//! builder, run **independently** from each `Layout` (no sharing of the
//! native grid).
//!
//! **Structural assert** (the `gerber-rt-001` tolerance, 0.5 nm — half
//! the 4.6 fixed-point quantum — not raw `==`: a stub length derived from
//! λ_g/4 is generically irrational and does not itself sit on the
//! nanometre grid, so the reimported vertex is the *nearest representable
//! nanometre*, not the pre-export float): every trace polygon reproduces
//! vertex-for-vertex.
//!
//! **Measured assert**: `two_port_board_job`'s voxelizer rounds every
//! coordinate to the nearest cell of the (0.3 mm) grid — five to six
//! orders of magnitude coarser than the nanometre quantization the import
//! adds — so the two Layouts voxelize to the *same* grid, and the FDTD
//! pipeline is otherwise deterministic (no RNG, no reduction-order
//! sensitivity: per-cell explicit updates). Bit-identical |S21| curves
//! are therefore the expected, not merely hoped-for, outcome, and that is
//! what is measured (see the assert below for the pinned deltas).
//!
//! `#[ignore]`'d (4 release FDTD solves — the `engine-miter-001` budget):
//!
//! ```bash
//! cargo test -p yee-engine --release --test import_twin -- --ignored --nocapture
//! ```

use yee_engine::board::{TwoPortBoardOptions, reference_through_line, two_port_board_job};
use yee_engine::{BackendChoice, JobEvent, JobSpec, sparams};
use yee_export::{GerberOptions, gerber_to_layout, layout_to_gerber};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.3e-3;

/// Hammerstad microstrip open-end length correction ΔL (mirrors
/// `sparams_stub_notch.rs`).
fn open_end_delta_l(w_m: f64, h_m: f64, e_eff: f64) -> f64 {
    let u = w_m / h_m;
    0.412 * h_m * ((e_eff + 0.3) * (u + 0.264)) / ((e_eff - 0.258) * (u + 0.8))
}

/// Build the native S.6 stub-notch board: a 3λ_g FR-4 microstrip feed
/// line with a mid-line open stub sized λ_g/4 − ΔL (Hammerstad
/// open-end-corrected), so closed forms predict the notch at 5 GHz.
/// Identical trace geometry to `sparams_stub_notch::stub_job(true)`.
fn native_stub_layout() -> Layout {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - open_end_delta_l(W_M, H_M, e_eff);

    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    let bbox = BBox::from_polygons(&[line.clone(), stub.clone()]);
    Layout {
        substrate: Substrate {
            eps_r: EPS_R,
            height_m: H_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces: vec![line, stub],
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: 50.0,
            },
            PortRef {
                at: Point2::new(l_m - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: 50.0,
            },
        ],
        bbox,
    }
}

fn board_opts() -> TwoPortBoardOptions {
    TwoPortBoardOptions {
        dx_m: DX_M,
        margin_cells: 34,
        air_above_cells: 34,
        n_steps: 9000,
        f0_hz: F0_HZ,
        bw_hz: 0.8 * F0_HZ,
        z0_ohm: 50.0,
        spacing_cells: 12,
        npml: 10,
        backend: BackendChoice::Cpu,
    }
}

fn run(spec: JobSpec) -> Vec<Vec<f64>> {
    let handle = yee_engine::submit(spec);
    handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result.probes),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event")
}

/// Launch-normalized double-ratio |S21| (ADR-0204) for one board layout,
/// via two runs (DUT + its through-line reference) both derived from
/// `board` alone.
fn s21_lin(board: &Layout, opts: &TwoPortBoardOptions, freqs: &[f64]) -> Vec<f64> {
    let reference = reference_through_line(board);
    let dut_job = two_port_board_job(board, opts).expect("dut job build failed");
    let ref_job = two_port_board_job(&reference, opts).expect("reference job build failed");
    assert_eq!(
        (dut_job.spec.nx, dut_job.spec.ny, dut_job.spec.nz),
        (ref_job.spec.nx, ref_job.spec.ny, ref_job.spec.nz),
        "DUT and reference must share a grid"
    );
    let (dt, spacing) = (dut_job.dt_s, dut_job.spacing_m);
    let dut_p = run(dut_job.spec);
    let ref_p = run(ref_job.spec);
    let transfer = |p: &[Vec<f64>]| {
        sparams::forward_transfer(
            [&p[0], &p[1], &p[2]],
            [&p[3], &p[4], &p[5]],
            dt,
            spacing,
            freqs,
        )
    };
    let t_dut = transfer(&dut_p);
    let t_ref = transfer(&ref_p);
    t_dut
        .iter()
        .zip(&t_ref)
        .map(|(d, r)| d.0.hypot(d.1) / r.0.hypot(r.1))
        .collect()
}

#[test]
#[ignore = "slow: 4 release FDTD solves; engine-import-twin-001 gate (FS.3.2c) — run with --release --ignored"]
fn imported_board_measures_the_same_notch_as_its_native_twin() {
    let native = native_stub_layout();

    // --- Structural twin: export -> Gerber bytes -> import -> Layout ---
    let gerber = layout_to_gerber(&native, &GerberOptions::default());
    let twin = gerber_to_layout(&gerber, native.substrate, native.ports.clone())
        .expect("re-import of our own writer output must succeed");

    assert_eq!(
        twin.traces.len(),
        native.traces.len(),
        "engine-import-twin-001 FAILED: polygon count changed on import"
    );
    for (t, n) in twin.traces.iter().zip(&native.traces) {
        assert_eq!(
            t.verts.len(),
            n.verts.len(),
            "engine-import-twin-001 FAILED: vertex count changed on import"
        );
        for (a, b) in t.verts.iter().zip(&n.verts) {
            // gerber-rt-001 tolerance: half the 4.6 fixed-point quantum
            // (1 nm), not raw `==` — the pre-export float is generically
            // not itself nanometre-aligned (see module docs).
            assert!(
                (a.x - b.x).abs() < 0.5e-9 && (a.y - b.y).abs() < 0.5e-9,
                "engine-import-twin-001 FAILED: vertex ({}, {}) vs native ({}, {})",
                a.x,
                a.y,
                b.x,
                b.y
            );
        }
    }

    // --- Measured twin: same measurement pipeline, run independently ---
    let opts = board_opts();
    // The ADR-0216 criterion band, mirroring the other stub-notch gates.
    let freqs: Vec<f64> = (0..=64).map(|n| 3.0e9 + n as f64 * 50.0e6).collect();

    eprintln!("engine-import-twin-001: native twin");
    let s21_native = s21_lin(&native, &opts, &freqs);
    eprintln!("engine-import-twin-001: imported twin");
    let s21_twin = s21_lin(&twin, &opts, &freqs);

    let db = |x: f64| 20.0 * x.log10();
    let min_of = |v: &[f64]| {
        v.iter()
            .enumerate()
            .map(|(n, x)| (n, *x))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .expect("empty")
    };
    let (n_notch, notch_native) = min_of(&s21_native);
    let notch_twin = s21_twin[n_notch];
    let f_notch = freqs[n_notch];
    let max_abs_delta = s21_native
        .iter()
        .zip(&s21_twin)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);

    eprintln!(
        "engine-import-twin-001: notch at {:.3} GHz | native {:.2} dB, imported twin {:.2} dB \
         | max |Δ|S21|| across the band = {:.3e} (linear)",
        f_notch / 1e9,
        db(notch_native),
        db(notch_twin),
        max_abs_delta,
    );

    // (a) The twin must still be a genuine notch (regression floor,
    // shared with the other stub-notch gates: ≥ 8 dB deep).
    assert!(
        db(notch_twin) <= -8.0,
        "engine-import-twin-001 FAILED: imported twin's notch only {:.1} dB deep (need ≤ -8 dB)",
        db(notch_twin)
    );
    // (b) The twins measure BIT-IDENTICAL |S21| curves. This is the
    // expected outcome, not merely hoped for: the voxelizer rounds every
    // coordinate to the nearest 0.3 mm cell, five to six orders of
    // magnitude coarser than the nanometre quantization import adds, so
    // native and twin voxelize to the identical grid and the (RNG-free,
    // per-cell-explicit) FDTD pipeline reproduces identical fields.
    // Measured 2026-07-24: max |Δ|S21|| = 0.0 (exact) over the full band.
    // A nonzero delta here would mean the import changed something the
    // voxelizer is actually sensitive to (vertex order, polygon winding,
    // a coordinate landing on a cell boundary) and would need root-cause
    // before this assert is loosened.
    assert_eq!(
        max_abs_delta, 0.0,
        "engine-import-twin-001 FAILED: native vs imported-twin |S21| differ by up to {max_abs_delta:.3e} \
         (linear) — vertex-exact import was expected to voxelize bit-identically; root-cause before \
         weakening this assert (vertex order? winding? a coordinate landing on a cell boundary?)"
    );
}
