//! Board-level two-port measurement plans (R.5b, ADR-0199).
//!
//! Promotes the voxelize → [`JobSpec`] pattern that every board gate
//! (S.8–S.12, R.0–R.4) carried as copy-pasted test code into a library
//! API: given a two-port [`yee_layout::Layout`] (feeds along ±x, ports at
//! the feed ends — what every `yee-layout`/`yee-filter` generator emits),
//! build the S.9/S.10-certified measurement fixture — CPML-xy walls,
//! aperture ports on the feed cross-sections, and two 3-probe triples for
//! the S.12 directional standing-wave observables. The studio's verify
//! command and future gates share this builder, so the fixture cannot
//! drift between them.

use serde::{Deserialize, Serialize};
use yee_layout::{Layout, Polygon};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

use crate::{AperturePortSpec, BackendChoice, BoundarySpec, JobSpec, MaterialsSpec, ProbeSpec};

/// Options for [`two_port_board_job`]. Defaults mirror the R.4 gate
/// scenario except for the drive centre/bandwidth, which have no sensible
/// universal default and must be set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TwoPortBoardOptions {
    /// Uniform cell size, metres.
    pub dx_m: f64,
    /// CPML margin cells on each x/y side.
    pub margin_cells: usize,
    /// Air cells above the substrate.
    pub air_above_cells: usize,
    /// Time steps to run.
    pub n_steps: usize,
    /// Drive centre frequency, Hz.
    pub f0_hz: f64,
    /// Drive −3 dB bandwidth, Hz.
    pub bw_hz: f64,
    /// Port/system impedance, ohms.
    pub z0_ohm: f64,
    /// Probe-triple spacing in cells (choose βd well inside (0, π)).
    pub spacing_cells: usize,
    /// CPML absorber depth in cells (inside the margin). Callers that vary
    /// `dx_m` across runs (the automesh convergence loop) must scale this
    /// with 1/dx so the absorber keeps its physical thickness — a
    /// cells-thin CPML at fine dx reflects long wavelengths.
    pub npml: usize,
    /// Backend to run on.
    pub backend: BackendChoice,
}

impl TwoPortBoardOptions {
    /// The R.4-gate-shaped defaults for a given drive band.
    pub fn for_band(f0_hz: f64, bw_hz: f64) -> Self {
        Self {
            dx_m: 0.2e-3,
            margin_cells: 34,
            air_above_cells: 34,
            n_steps: 13000,
            f0_hz,
            bw_hz,
            z0_ohm: 50.0,
            spacing_cells: 12,
            npml: 10,
            backend: BackendChoice::Cpu,
        }
    }
}

/// A ready-to-submit measurement job plus the constants needed to
/// post-process its probes with [`crate::sparams`].
#[derive(Debug, Clone)]
pub struct TwoPortBoardJob {
    /// The job. Probes 0–2 are triple A (input feed), 3–5 triple B
    /// (output feed), both ordered along +x.
    pub spec: JobSpec,
    /// The grid's dt (also set on the spec), seconds.
    pub dt_s: f64,
    /// Probe-triple spacing, metres.
    pub spacing_m: f64,
}

/// Build the straight-`Z₀`-through-line reference for a two-port layout:
/// the same substrate, ports, and bbox (→ the identical voxel grid), with
/// one straight trace at the port height spanning port to port.
pub fn reference_through_line(dut: &Layout) -> Layout {
    let p0 = dut.ports[0].at;
    let p1 = dut.ports[1].at;
    let w = dut.ports[0].width_m;
    Layout {
        substrate: dut.substrate,
        traces: vec![Polygon::rect(p0.x, p0.y - w / 2.0, p1.x - p0.x, w)],
        ports: dut.ports.clone(),
        bbox: dut.bbox,
    }
}

/// Voxelize a two-port board layout and express one measurement run as a
/// [`JobSpec`] with the S.9/S.10/S.12 fixture (CPML-xy + PEC ground/lid,
/// aperture ports, two directional probe triples).
///
/// # Errors
///
/// Returns a message when the layout has fewer than two ports, when the
/// feed band rasterizes to zero cells, or when the probe triples do not
/// fit on the feeds (feeds shorter than `2.4 mm + 2·spacing`).
pub fn two_port_board_job(
    layout: &Layout,
    opts: &TwoPortBoardOptions,
) -> Result<TwoPortBoardJob, String> {
    if layout.ports.len() < 2 {
        return Err("two_port_board_job needs a two-port layout".into());
    }
    let model = voxelize_microstrip(
        layout,
        &VoxelOptions {
            dx_m: opts.dx_m,
            xy_margin_cells: opts.margin_cells,
            air_above_cells: opts.air_above_cells,
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let k_top = model.port_cells[0].2;
    let load_cell = model.port_cells[1];
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - opts.margin_cells as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;

    // Aperture / probe j band: the feed width centred on the port height.
    let tap_y = layout.ports[0].at.y;
    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - opts.margin_cells as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx - tap_y).abs() < w_feed / 2.0 };
    let j_lo = (0..ny)
        .find(|&j| in_band(j))
        .ok_or("feed band rasterized to zero cells")?;
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);
    if j_hi <= j_lo {
        return Err("aperture band empty".into());
    }
    let j_strip = (j_lo + j_hi) / 2;

    // Probe triples on the feeds, ordered along +x.
    let spacing_m = opts.spacing_cells as f64 * dx;
    let clearance = 2.4e-3;
    let i_a0 = i_for(layout.ports[0].at.x + clearance);
    let i_b0 = i_for(layout.ports[1].at.x - clearance - 2.0 * spacing_m);
    let feed_span = (layout.ports[1].at.x - layout.ports[0].at.x).abs();
    if feed_span < 2.0 * (clearance + 2.0 * spacing_m) {
        return Err(format!(
            "feeds too short for the probe triples: port span {feed_span:.4} m \
             vs required {:.4} m",
            2.0 * (clearance + 2.0 * spacing_m)
        ));
    }

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

    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt()
        / (std::f64::consts::PI * opts.bw_hz))
        / dt)
        .ceil() as usize;

    let mk_probe = |i: usize| ProbeSpec {
        component: "ez".into(),
        cell: (i, j_strip, k_probe),
    };
    let mk_port = |i: usize, v0: f64| AperturePortSpec {
        i,
        j_lo,
        j_hi,
        k_lo: 0,
        k_top,
        resistance_ohm: opts.z0_ohm,
        v0,
        f0_hz: opts.f0_hz,
        bw_hz: opts.bw_hz,
        t0_steps,
        record: false,
    };
    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: opts.dx_m,
        n_steps: opts.n_steps,
        // Side-wall CPML, PEC ground/lid (S.9) — the board-level boundary.
        boundary: BoundarySpec::Cpml {
            npml: opts.npml,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![
            mk_port(model.port_cells[0].0, 1.0),
            mk_port(load_cell.0, 0.0),
        ],
        probes: vec![
            mk_probe(i_a0),
            mk_probe(i_a0 + opts.spacing_cells),
            mk_probe(i_a0 + 2 * opts.spacing_cells),
            mk_probe(i_b0),
            mk_probe(i_b0 + opts.spacing_cells),
            mk_probe(i_b0 + 2 * opts.spacing_cells),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        spacings: None,
        backend: opts.backend,
    };
    Ok(TwoPortBoardJob {
        spec,
        dt_s: dt,
        spacing_m,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use yee_layout::{BBox, Point2, PortRef, Substrate};

    fn line_layout(len_m: f64) -> Layout {
        let w = 1.5e-3;
        let traces = vec![Polygon::rect(0.0, -w / 2.0, len_m, w)];
        let bbox = BBox::from_polygons(&traces);
        Layout {
            substrate: Substrate {
                eps_r: 4.4,
                height_m: 0.8e-3,
                loss_tangent: 0.0,
                metal_thickness_m: 35e-6,
            },
            traces,
            ports: vec![
                PortRef {
                    at: Point2::new(0.0, 0.0),
                    width_m: w,
                    ref_impedance_ohm: 50.0,
                },
                PortRef {
                    at: Point2::new(len_m, 0.0),
                    width_m: w,
                    ref_impedance_ohm: 50.0,
                },
            ],
            bbox,
        }
    }

    #[test]
    fn plan_carries_the_certified_fixture() {
        let layout = line_layout(30.0e-3);
        let opts = TwoPortBoardOptions::for_band(5.0e9, 4.0e9);
        let job = two_port_board_job(&layout, &opts).unwrap();
        // CPML-xy + PEC z (the S.9 board boundary).
        match &job.spec.boundary {
            BoundarySpec::Cpml { npml, axes, faces } => {
                assert_eq!(*npml, 10);
                assert_eq!(*axes, [true, true, false]);
                assert!(faces.is_none());
            }
            other => panic!("wrong boundary: {other:?}"),
        }
        // Two aperture ports (drive + matched load), six probes in two
        // ordered triples.
        assert_eq!(job.spec.aperture_ports.len(), 2);
        assert!(job.spec.aperture_ports[0].v0 > 0.0);
        assert_eq!(job.spec.aperture_ports[1].v0, 0.0);
        assert_eq!(job.spec.probes.len(), 6);
        let i_of = |n: usize| job.spec.probes[n].cell.0;
        assert_eq!(i_of(1) - i_of(0), opts.spacing_cells);
        assert_eq!(i_of(2) - i_of(1), opts.spacing_cells);
        assert_eq!(i_of(4) - i_of(3), opts.spacing_cells);
        assert!(i_of(3) > i_of(2), "triple B must sit beyond triple A");
        assert!((job.spacing_m - opts.spacing_cells as f64 * opts.dx_m).abs() < 1e-15);
        // The reference shares grid dims with the DUT by construction.
        let reference = reference_through_line(&layout);
        let ref_job = two_port_board_job(&reference, &opts).unwrap();
        assert_eq!(
            (job.spec.nx, job.spec.ny, job.spec.nz),
            (ref_job.spec.nx, ref_job.spec.ny, ref_job.spec.nz)
        );
        assert_eq!(job.dt_s, ref_job.dt_s);
    }

    #[test]
    fn short_feeds_are_rejected() {
        let layout = line_layout(5.0e-3);
        let err =
            two_port_board_job(&layout, &TwoPortBoardOptions::for_band(5.0e9, 4.0e9)).unwrap_err();
        assert!(err.contains("too short"), "{err}");
    }
}
