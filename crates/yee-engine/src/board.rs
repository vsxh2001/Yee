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
        thin_wires: vec![],
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

// ---------------------------------------------------------------------------
// FS.0b.2a (ADR-0210 addendum): the graded two-port fixture — the
// engine-graded-001 gate's hand-rolled setup promoted to a library API so
// every graded measurement (filter/antenna verify, the FS.0b.2 converge
// integration) shares one certified builder.

use crate::automesh::{AutoSpacings, GradedMeshOptions, auto_spacings};
use yee_voxel::{GradedMicrostripModel, GradedVoxelGrid, voxelize_microstrip_graded};

/// Options for [`two_port_board_jobs_graded`].
#[derive(Debug, Clone)]
pub struct GradedBoardOptions {
    /// The graded-mesh rulebook knobs (margins, npml, growth, guard).
    pub mesh: GradedMeshOptions,
    /// Drive centre frequency, Hz.
    pub f0_hz: f64,
    /// Drive −3 dB bandwidth, Hz.
    pub bw_hz: f64,
    /// Port/system impedance, ohms.
    pub z0_ohm: f64,
    /// Probe-triple spacing in **coarse** cells (probes sit on
    /// bit-equal-coarse runs; see [`GradedTwoPortBoardJob`]).
    pub spacing_cells: usize,
    /// Time steps; `None` applies the engine-graded-001 rule
    /// `9000 · 0.3 mm / fine` (the FS.0a physical window at the fine
    /// spacing).
    pub n_steps: Option<usize>,
    /// Backend to run on.
    pub backend: BackendChoice,
}

impl GradedBoardOptions {
    /// Fixture defaults for a layout and drive band (mirrors
    /// [`TwoPortBoardOptions::for_band`] plus the graded rulebook's
    /// [`GradedMeshOptions::for_board`]).
    pub fn for_board(layout: &Layout, f_max_hz: f64, f0_hz: f64, bw_hz: f64) -> Self {
        Self {
            mesh: GradedMeshOptions::for_board(layout, f_max_hz),
            f0_hz,
            bw_hz,
            z0_ohm: 50.0,
            spacing_cells: 12,
            n_steps: None,
            backend: BackendChoice::Cpu,
        }
    }
}

/// A graded measurement job plus its post-processing constants. Probes
/// 0–2 are triple A (input feed), 3–5 triple B (output feed), both on
/// stretches of **bit-equal coarse** cells (`fit_standing_wave` requires
/// equal spacing), ordered along +x — the same layout as the uniform
/// [`TwoPortBoardJob`], so [`crate::sparams`] post-processing is
/// identical.
#[derive(Debug, Clone)]
pub struct GradedTwoPortBoardJob {
    /// The ready-to-submit job (graded `spacings` attached).
    pub spec: JobSpec,
    /// The graded Courant dt (also set on the spec), seconds.
    pub dt_s: f64,
    /// Probe-triple spacing, metres (`spacing_cells · coarse`).
    pub spacing_m: f64,
    /// Total cell count of the graded grid.
    pub cells: usize,
}

/// First index `i` such that a full `2·sp`-cell probe span starting at
/// node `i` lies on bit-equal-`coarse` cells with `nodes[i] ≥ x_min` and
/// `nodes[i + 2·sp] ≤ x_max`.
fn coarse_run_from_left(
    nodes: &[f64],
    widths: &[f64],
    coarse: f64,
    sp: usize,
    x_min: f64,
    x_max: f64,
) -> Option<usize> {
    (0..widths.len().saturating_sub(2 * sp)).find(|&i| {
        nodes[i] >= x_min
            && nodes[i + 2 * sp] <= x_max
            && widths[i..i + 2 * sp].iter().all(|d| *d == coarse)
    })
}

/// Last such index (searching from the right).
fn coarse_run_from_right(
    nodes: &[f64],
    widths: &[f64],
    coarse: f64,
    sp: usize,
    x_max: f64,
) -> Option<usize> {
    (0..widths.len().saturating_sub(2 * sp))
        .rev()
        .find(|&i| nodes[i + 2 * sp] <= x_max && widths[i..i + 2 * sp].iter().all(|d| *d == coarse))
}

/// Build the **(DUT, reference)** graded measurement pair: one
/// [`auto_spacings`] grid derived from the DUT, both layouts voxelized on
/// it. Returning the pair from one call encodes the ADR-0204
/// same-physical-problem lesson in the API shape — a caller cannot
/// accidentally measure the reference on its own (different) grid.
///
/// # Errors
///
/// Propagates [`auto_spacings`] failures and reports fixture-geometry
/// failures (feed band empty, no coarse probe stretch, triples overlap).
pub fn two_port_board_jobs_graded(
    dut: &Layout,
    f_max_hz: f64,
    opts: &GradedBoardOptions,
) -> Result<(GradedTwoPortBoardJob, GradedTwoPortBoardJob), String> {
    if dut.ports.len() < 2 {
        return Err("two_port_board_jobs_graded needs a two-port layout".into());
    }
    let spac = auto_spacings(dut, f_max_hz, &opts.mesh)?;
    let grid = GradedVoxelGrid {
        dx_m: spac.dx.clone(),
        dy_m: spac.dy.clone(),
        dz_m: spac.dz.clone(),
        x0_m: spac.x0_m,
        y0_m: spac.y0_m,
        k_gnd: spac.k_gnd,
        k_top: spac.k_top,
    };
    let reference = reference_through_line(dut);
    let dut_model = voxelize_microstrip_graded(dut, &grid);
    let ref_model = voxelize_microstrip_graded(&reference, &grid);
    if dut_model.dims != ref_model.dims {
        return Err("DUT and reference voxelized to different dims on one grid".into());
    }

    let fixture = graded_fixture(dut, &spac, &dut_model, opts)?;
    let dut_job = graded_job(&dut_model, &spac, &fixture, opts);
    let ref_job = graded_job(&ref_model, &spac, &fixture, opts);
    Ok((dut_job, ref_job))
}

/// The shared fixture geometry (probe/port placement, time base) —
/// derived once from the DUT grid and applied to both runs.
struct GradedFixture {
    j_lo: usize,
    j_hi: usize,
    j_strip: usize,
    k_probe: usize,
    i_a0: usize,
    i_b0: usize,
    dt: f64,
    n_steps: usize,
    t0_steps: usize,
    spacing_m: f64,
}

fn graded_fixture(
    dut: &Layout,
    spac: &AutoSpacings,
    model: &GradedMicrostripModel,
    opts: &GradedBoardOptions,
) -> Result<GradedFixture, String> {
    let (nx, ny, _) = model.dims;

    // Aperture / probe j band: the feed width centred on the port height,
    // from the true graded cell centres.
    let tap_y = dut.ports[0].at.y;
    let w_feed = dut.ports[0].width_m;
    let yc = |j: usize| (model.y_nodes_m[j] + model.y_nodes_m[j + 1]) / 2.0;
    let in_band = |j: usize| (yc(j) - tap_y).abs() < w_feed / 2.0;
    let j_lo = (0..ny)
        .find(|&j| in_band(j))
        .ok_or("feed band rasterized to zero cells")?;
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);
    let j_strip = (j_lo + j_hi) / 2;
    let k_probe = spac.k_top.saturating_sub(1).max(1);

    // Probe triples on bit-equal-coarse stretches, clear of the ports.
    let sp = opts.spacing_cells;
    let spacing_m = sp as f64 * spac.coarse_m;
    let clearance = 2.4e-3;
    let i_a0 = coarse_run_from_left(
        &model.x_nodes_m,
        &spac.dx,
        spac.coarse_m,
        sp,
        dut.ports[0].at.x + clearance,
        dut.bbox.max.x,
    )
    .ok_or("no uniform-coarse stretch for probe triple A")?;
    let i_b0 = coarse_run_from_right(
        &model.x_nodes_m,
        &spac.dx,
        spac.coarse_m,
        sp,
        dut.ports[1].at.x - clearance,
    )
    .ok_or("no uniform-coarse stretch for probe triple B")?;
    if i_b0 <= i_a0 + 2 * sp {
        return Err("probe triples overlap: feeds too short for the graded fixture".into());
    }
    let _ = nx;

    // Time base: 0.9× the graded Courant limit (every axis minimum is the
    // fine spacing), physical window at the fine spacing.
    let min_d = |a: &[f64]| a.iter().copied().fold(f64::INFINITY, f64::min);
    let (mx, my, mz) = (min_d(&spac.dx), min_d(&spac.dy), min_d(&spac.dz));
    const C0: f64 = 299_792_458.0;
    let dt = 0.9 / (C0 * (1.0 / (mx * mx) + 1.0 / (my * my) + 1.0 / (mz * mz)).sqrt());
    let n_steps = opts
        .n_steps
        .unwrap_or_else(|| (9000.0 * 0.3e-3 / mx).round() as usize);
    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt()
        / (std::f64::consts::PI * opts.bw_hz))
        / dt)
        .ceil() as usize;

    Ok(GradedFixture {
        j_lo,
        j_hi,
        j_strip,
        k_probe,
        i_a0,
        i_b0,
        dt,
        n_steps,
        t0_steps,
        spacing_m,
    })
}

fn graded_job(
    model: &GradedMicrostripModel,
    spac: &AutoSpacings,
    fx: &GradedFixture,
    opts: &GradedBoardOptions,
) -> GradedTwoPortBoardJob {
    let (nx, ny, nz) = model.dims;
    let materials = MaterialsSpec {
        eps_r_cells: Some(model.eps_r_cells.as_slice().unwrap().to_vec()),
        pec_mask_ex: Some(model.pec_mask_ex.as_slice().unwrap().to_vec()),
        pec_mask_ey: Some(model.pec_mask_ey.as_slice().unwrap().to_vec()),
        ..MaterialsSpec::default()
    };
    let mk_probe = |i: usize| ProbeSpec {
        component: "ez".into(),
        cell: (i, fx.j_strip, fx.k_probe),
    };
    let mk_port = |i: usize, v0: f64| AperturePortSpec {
        i,
        j_lo: fx.j_lo,
        j_hi: fx.j_hi,
        k_lo: 0,
        k_top: spac.k_top,
        resistance_ohm: opts.z0_ohm,
        v0,
        f0_hz: opts.f0_hz,
        bw_hz: opts.bw_hz,
        t0_steps: fx.t0_steps,
        record: false,
    };
    let sp = opts.spacing_cells;
    let spec = JobSpec {
        nx,
        ny,
        nz,
        // The nominal spacing feeding the CPML sigma_max recipe: the
        // absorbing layers are exactly coarse (ADR-0208 scope rule).
        dx_m: spac.coarse_m,
        n_steps: fx.n_steps,
        boundary: BoundarySpec::Cpml {
            npml: opts.mesh.npml,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![
            mk_port(model.port_cells[0].0, 1.0),
            mk_port(model.port_cells[1].0, 0.0),
        ],
        thin_wires: vec![],
        probes: vec![
            mk_probe(fx.i_a0),
            mk_probe(fx.i_a0 + sp),
            mk_probe(fx.i_a0 + 2 * sp),
            mk_probe(fx.i_b0),
            mk_probe(fx.i_b0 + sp),
            mk_probe(fx.i_b0 + 2 * sp),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(fx.dt),
        spacings: Some(spac.to_spacings()),
        backend: opts.backend,
    };
    GradedTwoPortBoardJob {
        spec,
        dt_s: fx.dt,
        spacing_m: fx.spacing_m,
        cells: nx * ny * nz,
    }
}

#[cfg(test)]
mod graded_tests {
    use super::*;
    use yee_layout::{BBox, Point2, PortRef, Substrate};

    /// The engine-graded-001 stub scenario, structurally.
    fn stub_layout() -> Layout {
        let sub = Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        };
        let w = 3.0e-3;
        let l = 66.0e-3;
        let line = Polygon::rect(0.0, 0.0, l, w);
        let stub = Polygon::rect(l / 2.0 - w / 2.0, w, w, 8.0e-3);
        let traces = vec![line, stub];
        let bbox = BBox::from_polygons(&traces);
        Layout {
            substrate: sub,
            traces,
            ports: vec![
                PortRef {
                    at: Point2::new(0.5e-3, w / 2.0),
                    width_m: w,
                    ref_impedance_ohm: 50.0,
                },
                PortRef {
                    at: Point2::new(l - 0.5e-3, w / 2.0),
                    width_m: w,
                    ref_impedance_ohm: 50.0,
                },
            ],
            bbox,
        }
    }

    #[test]
    fn graded_fixture_builds_the_shared_grid_pair() {
        let layout = stub_layout();
        let opts = GradedBoardOptions::for_board(&layout, 6.0e9, 5.0e9, 4.0e9);
        let (dut, reference) =
            two_port_board_jobs_graded(&layout, 6.0e9, &opts).expect("fixture failed");

        // Same grid for both runs (the ADR-0204 lesson, structurally).
        assert_eq!(
            (dut.spec.nx, dut.spec.ny, dut.spec.nz),
            (reference.spec.nx, reference.spec.ny, reference.spec.nz)
        );
        assert_eq!(dut.dt_s, reference.dt_s);
        assert_eq!(dut.spec.spacings, reference.spec.spacings);
        let spac = dut
            .spec
            .spacings
            .as_ref()
            .expect("graded spacings attached");

        // dt = 0.9x the graded Courant limit of the FINE spacing.
        let min_d = |a: &[f64]| a.iter().copied().fold(f64::INFINITY, f64::min);
        let (mx, my, mz) = (min_d(&spac.dx), min_d(&spac.dy), min_d(&spac.dz));
        let c0 = 299_792_458.0;
        let dt_courant = 0.9 / (c0 * (1.0 / (mx * mx) + 1.0 / (my * my) + 1.0 / (mz * mz)).sqrt());
        assert_eq!(dut.dt_s, dt_courant);

        // Probe triples sit on bit-equal-coarse runs: consecutive probe
        // i-gaps are exactly spacing_cells, and every crossed cell width
        // equals the coarse spacing (== the spec's nominal dx_m).
        let coarse = dut.spec.dx_m;
        for t in [0usize, 3] {
            let i0 = dut.spec.probes[t].cell.0;
            let i2 = dut.spec.probes[t + 2].cell.0;
            assert_eq!(i2 - i0, 2 * opts.spacing_cells, "triple stride");
            assert!(
                spac.dx[i0..i2].iter().all(|d| *d == coarse),
                "probe triple {t} not on bit-equal coarse cells"
            );
        }
        // Physical probe spacing is spacing_cells x coarse.
        assert_eq!(dut.spacing_m, opts.spacing_cells as f64 * coarse);

        // The graded grid must be meaningfully smaller than the uniform
        // grid at the same fine spacing would be (sanity, not the gate).
        assert_eq!(dut.cells, dut.spec.nx * dut.spec.ny * dut.spec.nz);
        let uniform_at_fine = ((spac.dx.iter().sum::<f64>() / mx).round()
            * (spac.dy.iter().sum::<f64>() / my).round()
            * (spac.dz.iter().sum::<f64>() / mz).round()) as usize;
        assert!(
            dut.cells < uniform_at_fine / 2,
            "graded {} vs uniform-at-fine {uniform_at_fine}: payoff lost",
            dut.cells
        );
    }
}
