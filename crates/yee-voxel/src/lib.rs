//! # yee-voxel
//!
//! Native bridge from a planar microstrip [`yee_layout::Layout`] to a
//! material-assigned [`yee_fdtd::YeeGrid`] (Filter Phase F1.1a, ADR-0091).
//!
//! [`voxelize_microstrip`] rasterizes the layout's top-metal polygons onto a
//! cubic Yee grid: a PEC ground plane at the bottom cell layer, a dielectric
//! substrate slab of `ε_r = layout.substrate.eps_r`, a one-cell-thick PEC
//! top-metal layer where a trace polygon covers the cell centre
//! (point-in-polygon ray-cast), and air above. The result is a `YeeGrid` with
//! per-cell `ε_r` and tangential `Ex`+`Ey` PEC masks already attached (a
//! horizontal PEC sheet zeroes the in-plane field, not the normal `Ez`), ready
//! for the F1.1b k/Q_e extraction step.
//!
//! This crate does **no** EM time-stepping — building the grid assigns
//! materials only, so it runs in milliseconds.
//!
//! ## WASM-safety boundary (ADR-0089)
//!
//! `yee-layout` deliberately has no `yee-fdtd` dependency so it stays
//! WASM-safe. `yee-voxel` is the separate **native** crate that bridges the
//! two; it depends on both.
//!
//! ## Z-stack (`z` increasing upward)
//!
//! The ground plane sits on the plane `z = 0` and the trace on `z = k_top·dx`,
//! so the dielectric fills the `n_sub` `E_z` edges spanning the ground-to-trace
//! gap (`k = 0 .. n_sub`) with **no air series gap** at the ground.
//!
//! - `k = 0` (plane `z = 0`): ground plane — tangential-E PEC over the whole
//!   x-y extent.
//! - `E_z` edges `k = 0 .. n_sub` (`z ∈ [0, n_sub·dx]`): substrate,
//!   `n_sub = round(height_m / dx)` (≥ 1), `ε_r = layout.substrate.eps_r`.
//! - `k = k_top = n_sub` (plane `z = n_sub·dx ≈ h`): top-metal layer —
//!   tangential-E PEC where a trace polygon covers the cell centre.
//! - `E_z` edges `k = k_top .. k_top + air_above_cells`: air (`ε_r = 1`).
//! - `nz = k_top + 1 + air_above_cells`.
//!
//! ## Example
//!
//! ```
//! use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
//! use yee_voxel::{voxelize_microstrip, VoxelOptions};
//!
//! let substrate = Substrate {
//!     eps_r: 4.4,
//!     height_m: 1.6e-3,
//!     loss_tangent: 0.0,
//!     metal_thickness_m: 35e-6,
//! };
//! let trace = Polygon::rect(0.0, 0.0, 3.0e-3, 20.0e-3);
//! let traces = vec![trace];
//! let bbox = BBox::from_polygons(&traces);
//! let layout = Layout {
//!     substrate,
//!     traces,
//!     ports: vec![PortRef {
//!         at: Point2::new(1.5e-3, 0.0),
//!         width_m: 3.0e-3,
//!         ref_impedance_ohm: 50.0,
//!     }],
//!     bbox,
//! };
//! let opts = VoxelOptions { dx_m: 0.5e-3, xy_margin_cells: 4, air_above_cells: 8 };
//! let model = voxelize_microstrip(&layout, &opts);
//! assert_eq!(model.port_cells.len(), 1);
//! ```

use ndarray::Array3;
use yee_fdtd::{LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver, YeeGrid};
use yee_layout::{Layout, Point2, Polygon};

/// Lumped-LC FDTD EM simulation of a synthesized filter board (Filter Phase
/// F2.3, ADR-0115): place each ladder L/C as a [`yee_fdtd::LumpedRlcPort`] on
/// the voxelized board, drive/sense two ports, and extract `|S21|(f)`.
pub mod lumped_sim;
pub use lumped_sim::{LumpedSimConfig, SERIES_ESR_OHM, simulate_lumped_board};

/// Voxelization parameters for [`voxelize_microstrip`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoxelOptions {
    /// Isotropic cell size `dx = dy = dz`, metres.
    pub dx_m: f64,
    /// Air margin (in cells) added around the layout bounding box in `x` and `y`.
    pub xy_margin_cells: usize,
    /// Number of air layers placed above the top-metal layer.
    pub air_above_cells: usize,
}

/// A material-assigned FDTD model produced from a planar microstrip layout.
#[derive(Debug)]
pub struct MicrostripModel {
    /// The Yee grid with per-cell `ε_r` and tangential `Ex`+`Ey` PEC masks attached.
    pub grid: YeeGrid,
    /// Grid cell dimensions `(nx, ny, nz)`.
    pub dims: (usize, usize, usize),
    /// Cell size used to build the grid, metres (echoes [`VoxelOptions::dx_m`]).
    pub dx_m: f64,
    /// Each layout port mapped to its `(i, j, k_top)` grid cell.
    pub port_cells: Vec<(usize, usize, usize)>,
}

/// Voxelize a planar microstrip [`Layout`] into a material-assigned
/// [`YeeGrid`].
///
/// Builds a cubic-cell grid sized from `layout.bbox` plus
/// [`VoxelOptions::xy_margin_cells`] of air in `x`/`y`, with the z-stack
/// described in the [crate] docs. Returns the assigned grid together with its
/// dimensions, the cell size, and the per-port cell indices.
///
/// # Panics
///
/// Panics if [`VoxelOptions::dx_m`] is not positive and finite, or if the
/// resulting grid would have a zero dimension (an empty / degenerate bounding
/// box with no margin).
pub fn voxelize_microstrip(layout: &Layout, opts: &VoxelOptions) -> MicrostripModel {
    let dx = opts.dx_m;
    assert!(
        dx.is_finite() && dx > 0.0,
        "VoxelOptions::dx_m must be positive and finite"
    );

    // --- X-Y extent: bbox padded by `margin` cells of air on every side. ---
    let margin = opts.xy_margin_cells as f64 * dx;
    let x0 = layout.bbox.min.x - margin;
    let x1 = layout.bbox.max.x + margin;
    let y0 = layout.bbox.min.y - margin;
    let y1 = layout.bbox.max.y + margin;

    let nx = ((x1 - x0) / dx).ceil() as usize;
    let ny = ((y1 - y0) / dx).ceil() as usize;
    assert!(
        nx > 0 && ny > 0,
        "voxelize_microstrip: degenerate x-y extent (nx={nx}, ny={ny}); \
         increase xy_margin_cells or check the layout bbox"
    );

    // --- Z-stack. ---
    //
    // The ground plane (tangential-E PEC) sits on the plane `z = 0` (the `k = 0`
    // staggered node). The trace (tangential-E PEC) sits on the plane
    // `z = k_top·dx`. The dielectric must fill the *entire* ground-to-trace gap,
    // which is `n_sub` cubic cells thick, so the trace is at `k_top = n_sub`
    // (giving a ground-to-trace spacing of exactly `n_sub·dx ≈ h`) and the
    // dielectric fills the `E_z` edges `k = 0 .. n_sub`.
    //
    // (Earlier the trace sat at `k_top = 1 + n_sub` with the dielectric filling
    // only `k = 1 ..= n_sub`, which left a one-cell *air* gap between the ground
    // and the substrate — a series air capacitance that drove the FDTD-measured
    // ε_eff ~20 % too low: the propagation gate fdtd-line-eeff-001 measured
    // ε_eff ≈ 2.5 vs analytic 3.33 until this gap was closed, after which it
    // measured 3.31, ADR-0108.)
    let n_sub = ((layout.substrate.height_m / dx).round() as usize).max(1);
    let k_top = n_sub; // ground (k=0) + n_sub dielectric cells -> top-metal layer
    let nz = k_top + 1 + opts.air_above_cells;

    let eps_r_sub = layout.substrate.eps_r;

    // --- Material arrays at YeeGrid's exact required shapes. ---
    // `with_eps_r_cells` requires `(nx+1, ny+1, nz+1)`.
    // A horizontal PEC sheet (the ground plane and the metal traces) zeroes the
    // TANGENTIAL field — `Ex` and `Ey` — on its plane, NOT the normal `Ez`. So
    // mask the two in-plane components at their staggered node positions:
    //   `with_pec_mask_ex` requires `(nx, ny+1, nz+1)`; `Ex` node at ((i+0.5)dx, j·dx).
    //   `with_pec_mask_ey` requires `(nx+1, ny, nz+1)`; `Ey` node at (i·dx, (j+0.5)dx).
    let mut eps = Array3::<f64>::from_elem((nx + 1, ny + 1, nz + 1), 1.0);
    let mut pec_ex = Array3::<bool>::from_elem((nx, ny + 1, nz + 1), false);
    let mut pec_ey = Array3::<bool>::from_elem((nx + 1, ny, nz + 1), false);

    let in_trace = |x: f64, y: f64| {
        layout
            .traces
            .iter()
            .any(|p| point_in_polygon(Point2 { x, y }, p))
    };

    // Substrate dielectric on the `E_z` edges spanning the ground-to-trace gap,
    // `k = 0 .. n_sub` (i.e. `0 ..= k_top − 1`); air `ε_r = 1.0` elsewhere is the
    // array default. Filling from `k = 0` (the cell directly above the ground
    // plane) leaves NO air series gap at the ground — see the z-stack note above.
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..n_sub {
                eps[(i, j, k)] = eps_r_sub;
            }
        }
    }

    // Tangential `Ex` PEC: ground plane (k = 0, whole layer) + traces
    // (k = k_top, where the `Ex` node ((i+0.5)dx, j·dx) lies under a trace).
    for i in 0..nx {
        for j in 0..=ny {
            pec_ex[(i, j, 0)] = true;
            let x = x0 + (i as f64 + 0.5) * dx;
            let y = y0 + j as f64 * dx;
            if in_trace(x, y) {
                pec_ex[(i, j, k_top)] = true;
            }
        }
    }

    // Tangential `Ey` PEC: ground plane (k = 0) + traces (k = k_top), with the
    // `Ey` node at (i·dx, (j+0.5)dx).
    for i in 0..=nx {
        for j in 0..ny {
            pec_ey[(i, j, 0)] = true;
            let x = x0 + i as f64 * dx;
            let y = y0 + (j as f64 + 0.5) * dx;
            if in_trace(x, y) {
                pec_ey[(i, j, k_top)] = true;
            }
        }
    }

    // --- Map layout ports to grid cells at the top-metal layer. ---
    let port_cells = layout
        .ports
        .iter()
        .map(|port| {
            let i = (((port.at.x - x0) / dx).floor() as isize).clamp(0, nx as isize - 1) as usize;
            let j = (((port.at.y - y0) / dx).floor() as isize).clamp(0, ny as isize - 1) as usize;
            (i, j, k_top)
        })
        .collect();

    let grid = YeeGrid::vacuum(nx, ny, nz, dx)
        .with_eps_r_cells(eps)
        .with_pec_mask_ex(pec_ex)
        .with_pec_mask_ey(pec_ey);

    MicrostripModel {
        grid,
        dims: (nx, ny, nz),
        dx_m: dx,
        port_cells,
    }
}

/// Test whether point `p` lies inside polygon `poly` via the standard
/// even-odd ray-cast (crossing-number) rule.
///
/// A horizontal ray is cast to `+x`; the point is inside when it crosses an odd
/// number of polygon edges. The polygon is treated as implicitly closed
/// (last vertex → first). Robust for the axis-aligned rectangular traces this
/// crate targets; points exactly on an edge are reported consistently (no
/// special handling) which is acceptable for cell-centre sampling.
fn point_in_polygon(p: Point2, poly: &Polygon) -> bool {
    let verts = &poly.verts;
    let n = verts.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = verts[i];
        let vj = verts[j];
        // Does edge (vj -> vi) straddle the horizontal line y = p.y, and is the
        // edge's x at y = p.y to the right of p.x?
        if (vi.y > p.y) != (vj.y > p.y) {
            let x_cross = (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y) + vi.x;
            if p.x < x_cross {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

// ===========================================================================
// FDTD propagation-based microstrip ε_eff driver (Filter Phase F1.1b.1,
// ADR-0108)
// ===========================================================================
//
// `run_line_eeff` is the SHIPPED full-wave gate for F1.1b.1. It supersedes the
// resonant-split method above, which PR #1 (7 CI/local iterations) proved
// unworkable for a microstrip on an *open* domain:
//
//   * a small hard-PEC box CONFINES the fringing / air-gap fields that set the
//     even/odd ε_eff split → the split collapses (k ≈ 0.02);
//   * a large hard-PEC box becomes a resonant CAVITY whose box modes swamp the
//     microstrip resonances → argmax picks box modes (k ≈ 0.01);
//   * open CPML walls remove both pathologies but KILL the resonator Q (the
//     λ/2 line radiates into the absorber) → zero detectable peaks.
//
// There is no box that is simultaneously high-Q and non-confining/
// non-resonant, so a *resonant* split is the wrong observable. The robust,
// textbook FDTD coupled-line characterization is a PROPAGATION measurement:
// drive a long, NON-resonant line, terminate both propagation ends in matched
// (CPML) loads, and read the phase velocity of the traveling wave directly off
// two probe planes a known distance apart. No Q, no cavity modes, no
// peak-picking.
//
// `run_line_eeff` measures the single-line ε_eff this way (the validated
// walking-skeleton gate `fdtd-line-eeff-001`). The same machinery generalizes
// to the even/odd coupled split via `run_coupled_line_eeff` (drive two coupled
// strips in-phase → even, anti-phase → odd; gate `fdtd-line-eeff-coupled-001`).

/// Configuration for the propagation-based ε_eff drivers [`run_line_eeff`] and
/// [`run_coupled_line_eeff`].
///
/// The defaults target a single straight FR-4 microstrip line driven at 5 GHz
/// (a short guided wavelength → fewer cells per λ_g → a cheap grid) with a
/// modulated-Gaussian lumped port and a long line terminated in stable hard-PEC
/// walls (the DFT is time-gated to the forward pulse, before the far-wall
/// reflection returns — see [`Self::gate_steps`]). Two probe planes a known
/// distance apart sample `E_z`; the single-bin DFT phase advance between them at
/// the drive centre frequency gives the phase velocity → `ε_eff`. The defaults
/// are tuned to converge the `fdtd-line-eeff-001` gate within ≤ 15 % on a
/// tractable grid (validated locally in the bounded dev container, ADR-0108).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineRunConfig {
    /// Isotropic Yee cell size `dx = dy = dz`, metres.
    pub dx_m: f64,
    /// Air margin (in cells) around the layout bounding box in `x`/`y`
    /// (forwarded to [`VoxelOptions::xy_margin_cells`]). Must keep the strip a
    /// few cells clear of the CPML region.
    pub xy_margin_cells: usize,
    /// Air layers above the top metal (forwarded to
    /// [`VoxelOptions::air_above_cells`]).
    pub air_above_cells: usize,
    /// Drive centre frequency `f0` (Hz). The single-bin DFT phase is read at
    /// this frequency; the guided wavelength `λ_g = c / (f0·√ε_eff)` sets the
    /// probe spacing in wavelengths.
    pub f0_hz: f64,
    /// Fractional FWHM bandwidth of the modulated-Gaussian drive, as a fraction
    /// of `f0` (`bw = freq_span · f0`). Broadband enough to launch a clean
    /// traveling pulse but narrow enough that `f0` dominates the spectrum at the
    /// probes.
    pub freq_span: f64,
    /// Total number of FDTD time steps. Must be long enough for the launched
    /// pulse to fully transit both probe planes and clear the far CPML, so the
    /// single-bin DFT at each probe integrates the whole passage.
    pub n_steps: usize,
    /// Series resistance (Ω) of the drive port. A moderate value (≈ the line
    /// impedance) launches the wave without a strong reflection at the feed.
    pub port_resistance_ohm: f64,
    /// Peak drive voltage (V) of the Gaussian-modulated pulse.
    pub drive_v0: f64,
    /// Number of time steps over which the single-bin DFT at each probe is
    /// **integrated**, counted from `t = 0`. This *time-gates* the measurement
    /// to the FORWARD traveling pulse, before the reflection off the far line
    /// end returns to the probes.
    ///
    /// The line is terminated in hard-PEC outer walls (stable — see below), so
    /// the far end reflects. A long enough line (several λ_g of clearance past
    /// the downstream probe) keeps that reflection out of the gate window: the
    /// forward pulse fully transits both probes, and the DFT integrates only it.
    /// `None` integrates the whole `n_steps` record (use only when the line is
    /// long enough that no reflection ever returns within `n_steps`).
    ///
    /// **Why PEC walls + a time gate, not CPML:** the determined PR #1 method
    /// was open CPML terminations, but `CpmlParams::for_grid` applies CPML on
    /// all six faces and is **late-time unstable** for a microstrip whose PEC
    /// ground plane and high-ε substrate run *into* the boundary region
    /// (container measurement, ADR-0108: fields diverge ~1e13 by ~2.5 k steps).
    /// A PEC box is unconditionally stable; gating the DFT to the forward
    /// passage recovers a clean, reflection-free traveling-wave phase — the
    /// standard "long line, short look" FDTD line characterization.
    pub gate_steps: Option<usize>,
    /// Absolute `x` position (metres, in layout coordinates) of probe plane A —
    /// the *upstream* probe. Placed a fraction of a wavelength past the feed so
    /// the launch transient has settled into a clean traveling wave.
    pub probe_a_x_m: f64,
    /// Absolute `x` position (metres, in layout coordinates) of probe plane B —
    /// the *downstream* probe. `probe_b_x_m > probe_a_x_m`; the phase advance is
    /// read between A and B, so the spacing `probe_b_x_m − probe_a_x_m` should
    /// be ≈ λ_g/4 … λ_g/2 — large enough to resolve cleanly, small enough that
    /// the true advance is `< 2π` (unambiguous; no phase wrap).
    pub probe_b_x_m: f64,
}

impl Default for LineRunConfig {
    /// Walking-skeleton defaults (single FR-4 line at 5 GHz, cheap grid). See
    /// the per-field docs; container-validated to pass `fdtd-line-eeff-001`
    /// within ≤ 15 % (ADR-0108).
    fn default() -> Self {
        Self {
            // 0.4 mm cells: ~4 substrate cells through a 1.6 mm FR-4 slab and
            // ~80 cells per guided wavelength at 5 GHz (λ_g ≈ 33 mm). Fine
            // enough to resolve the phase advance between the two probe planes.
            dx_m: 0.4e-3,
            xy_margin_cells: 14,
            air_above_cells: 16,
            // 5 GHz: short λ_g → fewer cells per wavelength → a cheaper grid for
            // a given probe spacing (and well within FR-4's quasi-static band).
            f0_hz: 5.0e9,
            // 80 % FWHM: a SHORT launched pulse (≈ 810 steps) so it fully
            // transits both probes before the far-wall reflection returns —
            // wide enough to fit inside the time gate, narrow enough that f0
            // still dominates the probe spectra.
            freq_span: 0.8,
            // The forward pulse reaches the downstream probe within a few
            // hundred steps and is ~800 steps long; 2.5 k steps capture its full
            // passage well before the gate cutoff.
            n_steps: 2_500,
            // ≈ line impedance: launch the wave with a weak feed reflection.
            port_resistance_ohm: 50.0,
            drive_v0: 1.0,
            // Time-gate the DFT to the forward passage (the test sets this to
            // just before the far-PEC reflection returns to the downstream
            // probe — see `gate_steps` docs). `None` here; the test always sets
            // it explicitly to match its line length.
            gate_steps: None,
            // Probe planes in absolute layout-x metres. The default geometry is
            // a single ~6 λ_g line (λ_g ≈ 33 mm at 5 GHz); place the probes near
            // the middle, a third of a wavelength apart → Δx ≈ λ_g/3, a ~120°
            // advance that is unambiguous (< 2π) and well-resolved. The test
            // overrides these to match its line length.
            probe_a_x_m: 82.0e-3,
            probe_b_x_m: 93.0e-3,
        }
    }
}

/// Result of a propagation-based ε_eff measurement ([`run_line_eeff`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineRunResult {
    /// FDTD-extracted effective relative permittivity
    /// `ε_eff = (c / v_p)²`, from the measured phase velocity.
    pub eps_eff: f64,
    /// Measured phase velocity `v_p = ω·Δx / Δφ` (m/s).
    pub v_p: f64,
    /// Unwrapped phase advance `Δφ = φ_A − φ_B` between the two probe planes at
    /// the drive centre frequency (radians; positive for a `+x`-traveling
    /// wave).
    pub delta_phi: f64,
    /// Probe-plane separation `Δx` along the line (metres).
    pub delta_x: f64,
}

/// Run a full-wave FDTD solve on a single straight microstrip line and extract
/// its effective relative permittivity `ε_eff` from the measured **phase
/// velocity** of a traveling wave.
///
/// This is the SHIPPED walking-skeleton gate for Filter Phase F1.1b.1
/// (`fdtd-line-eeff-001`, ADR-0108): the first full-wave EM solve in the
/// filter-design pipeline. It composes shipped primitives —
/// [`voxelize_microstrip`] (layout→grid) and `yee-fdtd`'s [`LumpedRlcPort`]
/// drive + [`WalkingSkeletonSolver`] time-stepping.
///
/// # Method
///
/// A *non-resonant* propagation measurement (the textbook FDTD line
/// characterization; replaces the unworkable resonant-split method — see the
/// module comment and PR #1):
///
/// 1. Voxelize the line with hard-PEC outer walls (unconditionally stable; CPML
///    is late-time unstable for a microstrip whose PEC ground / high-ε substrate
///    run into the boundary — ADR-0108). The line is made several λ_g long so
///    the far-wall reflection returns *after* the forward pulse has cleared the
///    probes.
/// 2. Drive `E_z` at one end with a modulated-Gaussian lumped port, launching a
///    `+x`-traveling pulse.
/// 3. Record `E_z` in the substrate under the strip centre at two probe planes A
///    and B a known distance `Δx` apart along the line.
/// 4. Time-gated single-bin DFT at `f0` at each probe (gated to the forward
///    passage via [`LineRunConfig::gate_steps`], before the reflection returns)
///    → complex phasors; the phase advance `Δφ = φ_A − φ_B` (positive,
///    downstream lags) gives the phase velocity `v_p = ω·Δx / Δφ` and hence
///    `ε_eff = (c / v_p)²`.
///
/// # Panics
///
/// Panics if the layout has no ports (a drive port is required), if the probe
/// positions do not satisfy `probe_b_x_m > probe_a_x_m`, or if the resulting
/// probe planes collapse to the same grid column.
pub fn run_line_eeff(layout: &Layout, cfg: &LineRunConfig) -> LineRunResult {
    assert!(
        !layout.ports.is_empty(),
        "run_line_eeff: need ≥ 1 port (a drive port)"
    );
    assert!(
        cfg.probe_b_x_m > cfg.probe_a_x_m,
        "run_line_eeff: require probe_b_x_m ({}) > probe_a_x_m ({})",
        cfg.probe_b_x_m,
        cfg.probe_a_x_m
    );

    // --- 1. Voxelize: layout -> material-assigned grid + port cells. --------
    let opts = VoxelOptions {
        dx_m: cfg.dx_m,
        xy_margin_cells: cfg.xy_margin_cells,
        air_above_cells: cfg.air_above_cells,
    };
    let model = voxelize_microstrip(layout, &opts);
    let (nx, _ny, _nz) = model.dims;
    let drive_cell = model.port_cells[0];
    let (_i_drive, j_strip, k_top) = drive_cell;
    let dt = model.grid.dt;
    let dx = model.dx_m;

    // Probe `E_z` in the SUBSTRATE, under the strip — that is where the
    // quasi-TEM mode's dominant vertical field lives (between the trace at
    // `k_top` and the ground at `k = 0`). The `E_z` node at `(i, j, k)` spans
    // `z ∈ [k·dx, (k+1)·dx]`, so the substrate cell directly beneath the trace
    // is `k_top − 1`. (The node at `k_top` itself sits in the air just above
    // the metal, where the field is far weaker.)
    let k_probe = k_top.saturating_sub(1).max(1);

    // Map the two probe planes from absolute layout-x to grid columns. The
    // voxelizer's grid origin is `x0 = bbox.min.x − xy_margin_cells·dx` (see
    // `voxelize_microstrip`), so the column for layout-x `xp` is
    // `round((xp − x0)/dx)`. Sample `E_z` in the substrate under the strip
    // centre column `j_strip` at `k_probe`. The phase advance is read between
    // these two interior columns; the probe spacing `Δx = probe_b_x_m −
    // probe_a_x_m` is chosen (by the caller) < λ_g so the advance is
    // unambiguous.
    let x0 = layout.bbox.min.x - cfg.xy_margin_cells as f64 * dx;
    let i_for = |xp: f64| -> usize {
        (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize
    };
    let i_a = i_for(cfg.probe_a_x_m);
    let i_b = i_for(cfg.probe_b_x_m);
    assert!(
        i_b > i_a,
        "run_line_eeff: probe planes collapsed to the same column (i_a = i_b = {i_a}); \
         widen the probe spacing or refine dx"
    );
    let probe_a = (i_a, j_strip, k_probe);
    let probe_b = (i_b, j_strip, k_probe);
    let delta_x = (i_b - i_a) as f64 * dx;

    // --- 2. Hard-PEC outer walls (unconditionally stable). The CPML path was
    //        found late-time unstable for a microstrip whose PEC ground / high-ε
    //        substrate run into the boundary (ADR-0108); a PEC box plus a
    //        time-gated DFT (`gate_steps`) on the FORWARD pulse recovers a
    //        clean, reflection-free traveling-wave phase instead. --------------
    let mut solver = WalkingSkeletonSolver::new(model.grid);

    // --- 3. Drive port: modulated-Gaussian launch at the −x end. ------------
    let bw = cfg.freq_span * cfg.f0_hz;
    // Centre the pulse ~3.5 time-constants in so its t = 0 tail is negligible.
    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw))
        / dt)
        .ceil() as usize;
    let wave = SourceWaveform::GaussianPulse {
        v0: cfg.drive_v0,
        f0: cfg.f0_hz,
        bw,
        t0_steps,
    };
    let mut port = LumpedRlcPort::pure_resistor(drive_cell, cfg.port_resistance_ohm, wave);

    // --- 4. Time-step; integrate the single-bin DFT at each probe over the
    //        time gate (the forward passage, before the far-wall reflection
    //        returns). The step body matches the canonical `step_with_source`
    //        order: H + boundary, E + boundary, then the lumped-port correction
    //        (after the E boundary, per `LumpedRlcPort::correct_e`'s call site),
    //        then advance the clock. -------------------------------------------
    let mut acc = [0.0_f64; 4]; // [reA, imA, reB, imB] single-bin DFT
    let omega = 2.0 * std::f64::consts::PI * cfg.f0_hz;
    let gate = cfg.gate_steps.unwrap_or(cfg.n_steps).min(cfg.n_steps);
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        solver.apply_cpml_e();
        port.correct_e(solver.grid_mut(), n, dt);

        solver.advance_clock();

        if n < gate {
            let grid = solver.grid();
            let ez_a = grid.ez[probe_a];
            let ez_b = grid.ez[probe_b];
            // Single-bin DFT accumulation at f0. The gate confines the record to
            // the forward traveling pulse (finite support, fully contained), so
            // the rectangular window introduces no high-Q sidelobe comb.
            let phase = omega * n as f64 * dt;
            let (c, s) = (phase.cos(), phase.sin());
            acc[0] += ez_a * c;
            acc[1] -= ez_a * s;
            acc[2] += ez_b * c;
            acc[3] -= ez_b * s;
        }
    }

    // --- 5. Phase advance A → B → phase velocity → ε_eff. -------------------
    let phi_a = acc[1].atan2(acc[0]);
    let phi_b = acc[3].atan2(acc[2]);
    // A +x-traveling wave lags downstream, so φ decreases from A to B; the
    // phase advance Δφ = φ_A − φ_B is positive. Reduce into (0, 2π) so a wrap
    // across the atan2 branch cut does not flip the sign (the probe spacing is
    // chosen < λ_g so the true advance is < 2π).
    let mut delta_phi = phi_a - phi_b;
    while delta_phi <= 0.0 {
        delta_phi += 2.0 * std::f64::consts::PI;
    }
    while delta_phi > 2.0 * std::f64::consts::PI {
        delta_phi -= 2.0 * std::f64::consts::PI;
    }

    let v_p = omega * delta_x / delta_phi;
    let eps_eff = (C0_M_S / v_p).powi(2);

    let mag_a = (acc[0] * acc[0] + acc[1] * acc[1]).sqrt();
    let mag_b = (acc[2] * acc[2] + acc[3] * acc[3]).sqrt();
    eprintln!(
        "[fdtd-line-eeff DIAG] probes i_a={i_a} i_b={i_b} (j={j_strip}, k={k_probe}), \
         Δx={:.3} mm | φ_A={:.4} (|A|={:.3e}) φ_B={:.4} (|B|={:.3e}) | \
         Δφ={:.4} rad, v_p={:.4e} m/s, ε_eff={:.4}",
        delta_x * 1e3,
        phi_a,
        mag_a,
        phi_b,
        mag_b,
        delta_phi,
        v_p,
        eps_eff,
    );

    LineRunResult {
        eps_eff,
        v_p,
        delta_phi,
        delta_x,
    }
}

/// Speed of light in vacuum (m/s), for the `ε_eff = (c / v_p)²` conversion.
const C0_M_S: f64 = 299_792_458.0;

/// Result of a coupled-line even/odd ε_eff measurement
/// ([`run_coupled_line_eeff`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupledLineResult {
    /// Even-mode (in-phase drive) effective permittivity `ε_eff,e`.
    pub eps_eff_e: f64,
    /// Odd-mode (anti-phase drive) effective permittivity `ε_eff,o`.
    pub eps_eff_o: f64,
    /// The even/odd ε_eff-split `k = (ε_eff,e − ε_eff,o)/(ε_eff,e + ε_eff,o)`.
    /// (The even mode concentrates more field in the substrate, so normally
    /// `ε_eff,e ≥ ε_eff,o` and `k ≥ 0`.)
    pub k_split: f64,
}

/// Run two full-wave FDTD propagation solves on a *coupled* microstrip pair and
/// extract the even- and odd-mode effective permittivities `ε_eff,e` /
/// `ε_eff,o` from their phase velocities (Filter Phase F1.1b.1 coupled
/// follow-on, gate `fdtd-line-eeff-coupled-001`, ADR-0108).
///
/// # Method
///
/// The same reflection-free propagation measurement as [`run_line_eeff`],
/// applied to two parallel edge-coupled strips, run twice:
///
/// 1. **Even mode** — drive *both* strips in phase (`+v0`, `+v0`). The symmetric
///    excitation launches only the even supermode; the parity-matched probe
///    `E_z(strip1) + E_z(strip2)` reconstructs it. Its phase velocity →
///    `ε_eff,e`.
/// 2. **Odd mode** — drive the strips anti-phase (`+v0`, `−v0`). The
///    antisymmetric excitation launches only the odd supermode; the
///    parity-matched probe `E_z(strip1) − E_z(strip2)` reconstructs it (a SUM
///    would cancel). Its phase velocity → `ε_eff,o`.
///
/// Each run is the stable PEC-box + time-gated-DFT propagation measurement of
/// [`run_line_eeff`]; the split `k = (ε_eff,e − ε_eff,o)/(ε_eff,e + ε_eff,o)`
/// is the physically-correct even/odd ε_eff-split (NOT the impedance coupling —
/// see PR #1 root cause, ADR-0108).
///
/// # Panics
///
/// Panics if the layout has fewer than two ports (one drive per strip), or for
/// the same probe-geometry preconditions as [`run_line_eeff`].
pub fn run_coupled_line_eeff(
    layout: &Layout,
    cfg: &LineRunConfig,
    probe_b_offset_m: f64,
) -> CoupledLineResult {
    assert!(
        layout.ports.len() >= 2,
        "run_coupled_line_eeff: need ≥ 2 ports (one drive per strip); got {}",
        layout.ports.len()
    );
    let eps_eff_e = coupled_mode_eeff(layout, cfg, probe_b_offset_m, false);
    let eps_eff_o = coupled_mode_eeff(layout, cfg, probe_b_offset_m, true);
    // Order so the higher-ε (even) mode is first; the split is then ≥ 0.
    let (eps_eff_e, eps_eff_o) = if eps_eff_e >= eps_eff_o {
        (eps_eff_e, eps_eff_o)
    } else {
        (eps_eff_o, eps_eff_e)
    };
    let k_split = (eps_eff_e - eps_eff_o) / (eps_eff_e + eps_eff_o);
    CoupledLineResult {
        eps_eff_e,
        eps_eff_o,
        k_split,
    }
}

/// Drive a coupled microstrip pair into one supermode (even if `anti_phase ==
/// false`, odd if `true`) and return that mode's propagation ε_eff via the
/// time-gated phase-velocity measurement.
///
/// `probe_b_offset_m` is the downstream-probe offset from probe A (≈ λ_g/3); the
/// upstream probe A sits at `cfg.probe_a_x_m`. Both strips are probed in the
/// substrate (`k_top − 1`); the parity-matched combination
/// `E_z(strip1) ± E_z(strip2)` isolates the supermode.
fn coupled_mode_eeff(
    layout: &Layout,
    cfg: &LineRunConfig,
    probe_b_offset_m: f64,
    anti_phase: bool,
) -> f64 {
    let opts = VoxelOptions {
        dx_m: cfg.dx_m,
        xy_margin_cells: cfg.xy_margin_cells,
        air_above_cells: cfg.air_above_cells,
    };
    let model = voxelize_microstrip(layout, &opts);
    let (nx, _ny, _nz) = model.dims;
    let cell0 = model.port_cells[0];
    let cell1 = model.port_cells[1];
    let (_i0, j0, k_top) = cell0;
    let (_i1, j1, _k1) = cell1;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - cfg.xy_margin_cells as f64 * dx;
    let i_for = |xp: f64| -> usize {
        (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize
    };
    let i_a = i_for(cfg.probe_a_x_m);
    let i_b = i_for(cfg.probe_a_x_m + probe_b_offset_m);
    assert!(
        i_b > i_a,
        "run_coupled_line_eeff: probe planes collapsed (i_a = i_b = {i_a})"
    );
    let delta_x = (i_b - i_a) as f64 * dx;

    // Stable PEC outer walls; the lateral PEC walls are far from the strips.
    let mut solver = WalkingSkeletonSolver::new(model.grid);

    let bw = cfg.freq_span * cfg.f0_hz;
    let t0_steps = ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw))
        / dt)
        .ceil() as usize;
    let make_wave = |v0: f64| SourceWaveform::GaussianPulse {
        v0,
        f0: cfg.f0_hz,
        bw,
        t0_steps,
    };
    let v1 = if anti_phase {
        -cfg.drive_v0
    } else {
        cfg.drive_v0
    };
    let mut port0 =
        LumpedRlcPort::pure_resistor(cell0, cfg.port_resistance_ohm, make_wave(cfg.drive_v0));
    let mut port1 = LumpedRlcPort::pure_resistor(cell1, cfg.port_resistance_ohm, make_wave(v1));

    // Probe both strips at each plane; the parity-matched combination
    // (strip0 + strip1 for even, strip0 − strip1 for odd) reconstructs the
    // supermode, mirroring the drive parity (a mismatched combination would
    // cancel the mode of interest).
    let probe0_a = (i_a, j0, k_probe);
    let probe1_a = (i_a, j1, k_probe);
    let probe0_b = (i_b, j0, k_probe);
    let probe1_b = (i_b, j1, k_probe);

    let mut acc = [0.0_f64; 4]; // [reA, imA, reB, imB]
    let omega = 2.0 * std::f64::consts::PI * cfg.f0_hz;
    let gate = cfg.gate_steps.unwrap_or(cfg.n_steps).min(cfg.n_steps);
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        solver.apply_cpml_e();
        port0.correct_e(solver.grid_mut(), n, dt);
        port1.correct_e(solver.grid_mut(), n, dt);

        solver.advance_clock();

        if n < gate {
            let grid = solver.grid();
            let comb = |p0: (usize, usize, usize), p1: (usize, usize, usize)| {
                if anti_phase {
                    grid.ez[p0] - grid.ez[p1]
                } else {
                    grid.ez[p0] + grid.ez[p1]
                }
            };
            let ez_a = comb(probe0_a, probe1_a);
            let ez_b = comb(probe0_b, probe1_b);
            let phase = omega * n as f64 * dt;
            let (c, s) = (phase.cos(), phase.sin());
            acc[0] += ez_a * c;
            acc[1] -= ez_a * s;
            acc[2] += ez_b * c;
            acc[3] -= ez_b * s;
        }
    }

    let phi_a = acc[1].atan2(acc[0]);
    let phi_b = acc[3].atan2(acc[2]);
    let mut delta_phi = phi_a - phi_b;
    while delta_phi <= 0.0 {
        delta_phi += 2.0 * std::f64::consts::PI;
    }
    while delta_phi > 2.0 * std::f64::consts::PI {
        delta_phi -= 2.0 * std::f64::consts::PI;
    }
    let v_p = omega * delta_x / delta_phi;
    let eps_eff = (C0_M_S / v_p).powi(2);

    let mode = if anti_phase { "odd " } else { "even" };
    eprintln!(
        "[fdtd-line-eeff-coupled DIAG] {mode}: Δx={:.3} mm, Δφ={:.4} rad, \
         v_p={:.4e} m/s, ε_eff={:.4}",
        delta_x * 1e3,
        delta_phi,
        v_p,
        eps_eff,
    );
    eps_eff
}
