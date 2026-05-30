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
//! - `k = 0`: ground plane — PEC over the whole x-y extent.
//! - `k = 1 ..= n_sub`: substrate, `n_sub = round(height_m / dx)` (≥ 1),
//!   `ε_r = layout.substrate.eps_r`.
//! - `k = k_top = 1 + n_sub`: top-metal layer — PEC where a trace polygon
//!   covers the cell centre; `ε_r = 1` elsewhere.
//! - `k = k_top + 1 ..= k_top + air_above_cells`: air (`ε_r = 1`).
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
use yee_filter::extract_coupling;
use yee_layout::{Layout, Point2, Polygon};

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
    let n_sub = ((layout.substrate.height_m / dx).round() as usize).max(1);
    let k_top = 1 + n_sub; // ground (k=0) + n_sub substrate layers -> top-metal layer
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

    // Substrate dielectric per cell, k = 1 ..= n_sub (air ε_r = 1.0 elsewhere is
    // the array default).
    for i in 0..nx {
        for j in 0..ny {
            for k in 1..=n_sub {
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
// FDTD coupled-resonator driver (Filter Phase F1.1b.1, ADR-0108)
// ===========================================================================
//
// `run_coupled_pair` is the first full-wave EM solve in the filter-design
// pipeline. It voxelizes a two-resonator microstrip [`Layout`], drives one
// resonator with a weakly-coupled lumped port carrying a broadband
// Gaussian-modulated pulse, probes the second resonator's vertical `E_z`,
// single-bin-DFTs the probe time series to a magnitude spectrum, and inverts
// the two split resonance peaks to a coupling coefficient `k` via the shipped
// [`yee_filter::extract_coupling`] (`k = (f_odd² − f_even²)/(f_odd² + f_even²)`).
//
// The validation gate (`tests/fdtd_coupling_001.rs`, `#[ignore]`'d) cross-checks
// the FDTD `k` against the analytic Kirschning-Jansen coupled-line reference
// (`yee_layout::coupled_microstrip` / `coupling_coefficient`). The FDTD run is
// multi-minute, so the gate runs only in a dedicated CI `--release` job — never
// in the default `cargo test`. See ADR-0108.

/// Configuration for the [`run_coupled_pair`] FDTD coupled-resonator driver.
///
/// The defaults target a *walking-skeleton* run: a coarse cubic grid and a
/// modest step count chosen for a tractable CI wall-time rather than tight
/// accuracy (the gate tolerance is loose, ≤ 15 %). The physics knobs
/// (resolution, step count, drive bandwidth) are deliberately conservative and
/// documented; they have **not** been FDTD-tuned on the dev box (memory- and
/// CPU-constrained — see ADR-0108), so CI is the first place the chosen values
/// are exercised end to end.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupledRunConfig {
    /// Isotropic Yee cell size `dx = dy = dz`, metres. Smaller is more accurate
    /// but quadratically (×3 in 3-D) more expensive; the default is a coarse
    /// resolution sufficient to resolve the split resonances on a tractable
    /// grid.
    pub dx_m: f64,
    /// Air margin (in cells) around the layout bounding box in `x`/`y`
    /// (forwarded to [`VoxelOptions::xy_margin_cells`]).
    pub xy_margin_cells: usize,
    /// Air layers above the top metal (forwarded to
    /// [`VoxelOptions::air_above_cells`]).
    pub air_above_cells: usize,
    /// Synchronous resonator centre frequency `f0` (Hz). The DFT magnitude
    /// spectrum is scanned over `[f0·(1 − span), f0·(1 + span)]`; the two split
    /// resonances of a coupled pair straddle `f0`, so the scan window must
    /// bracket both.
    pub f0_hz: f64,
    /// Fractional half-width of the DFT scan window around `f0`
    /// (`f_lo = f0·(1 − span)`, `f_hi = f0·(1 + span)`). Wide enough to capture
    /// both split peaks for the coupling levels this driver targets.
    pub freq_span: f64,
    /// Number of linearly-spaced candidate frequencies in the single-bin DFT
    /// scan. Sets the frequency resolution `Δf = (f_hi − f_lo)/(n_freq_bins−1)`
    /// and hence the precision with which the two split peaks are located.
    pub n_freq_bins: usize,
    /// Total number of FDTD time steps. The simulated time `n_steps · dt` sets
    /// the DFT frequency resolution floor `1/(n_steps · dt)`; the run must be
    /// long enough for the two split resonances to separate cleanly.
    pub n_steps: usize,
    /// Series resistance (Ω) of both lumped ports. A *large* value makes the
    /// ports weakly coupled probes (high loaded Q, sharp split peaks) rather
    /// than matched loads that would damp the resonances away.
    pub port_resistance_ohm: f64,
    /// Peak drive voltage (V) of the Gaussian-modulated pulse fed into the
    /// driven port's series EMF.
    pub drive_v0: f64,
}

impl Default for CoupledRunConfig {
    /// Walking-skeleton defaults (coarse grid, modest step count). See the
    /// per-field docs for the rationale; these are CI-validated, not
    /// dev-box-tuned (ADR-0108).
    fn default() -> Self {
        Self {
            // ~0.2 mm cells. For a ~1.6 mm FR-4 substrate this is ~8 cells
            // through the dielectric (n_sub) and ~5 across a 1 mm coupling gap.
            // The first CI run at 0.4 mm (4 substrate / 2.5 gap cells) measured
            // εeff ≈ 2.5 (vs analytic 3.56) and a 2.1 % split (vs ~19 %) → k
            // 0.021 vs 0.173 (88 % err): the coarse grid under-resolved both the
            // dielectric loading and the gap fields. 0.2 mm roughly doubles the
            // resolution of both; CI re-stresses whether it lands within 15 %.
            dx_m: 0.2e-3,
            xy_margin_cells: 6,
            air_above_cells: 8,
            // 2.4 GHz synchronous centre (the project's canonical ISM band; the
            // gate builds its resonators around this).
            f0_hz: 2.4e9,
            // ±35 % brackets both split peaks for k up to a few tenths.
            freq_span: 0.35,
            n_freq_bins: 600,
            // Halving dx halves the CFL-stable dt, so the step count is doubled
            // (vs the original 40k at 0.4 mm) to preserve the simulated time
            // window — the DFT resolution floor 1/(n_steps·dt) — that lets the
            // two split peaks resolve. CI stresses whether it is sufficient.
            n_steps: 80_000,
            // Weak coupling: ~10× a 50 Ω system impedance keeps loaded Q high.
            port_resistance_ohm: 500.0,
            drive_v0: 1.0,
        }
    }
}

/// Result of a [`run_coupled_pair`] FDTD coupled-resonator solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupledRunResult {
    /// Extracted inter-resonator coupling coefficient
    /// `k = (f_odd² − f_even²)/(f_odd² + f_even²)` (from
    /// [`yee_filter::extract_coupling`]).
    pub k: f64,
    /// Even-mode split-resonance frequency (Hz) — the **lower** of the two
    /// peaks (the even mode has the higher effective permittivity, hence the
    /// lower frequency).
    pub f_even: f64,
    /// Odd-mode split-resonance frequency (Hz) — the **upper** of the two peaks.
    pub f_odd: f64,
}

/// Run a full-wave FDTD solve on a coupled microstrip-resonator [`Layout`] and
/// extract its inter-resonator coupling coefficient `k`.
///
/// This is the integration driver for Filter Phase F1.1b.1 (ADR-0108): it
/// composes primitives that all already ship —
/// [`voxelize_microstrip`] (the layout→grid bridge), `yee-fdtd`'s
/// [`LumpedRlcPort`] drive/probe + [`WalkingSkeletonSolver`] time-stepping, and
/// [`yee_filter::extract_coupling`] (the split-peak `k` inversion).
///
/// # Method
///
/// 1. **Voxelize** the two-resonator layout into a material-assigned grid
///    (PEC ground + PEC traces + substrate slab). The default
///    [`WalkingSkeletonSolver::new`] gives hard-PEC outer walls, i.e. a closed
///    resonant box — appropriate for a weak-coupling resonance measurement.
/// 2. **Drive** the first layout port with a [`LumpedRlcPort`] carrying a
///    broadband [`SourceWaveform::GaussianPulse`] centred on
///    [`CoupledRunConfig::f0_hz`]; the second port is a passive (open-EMF)
///    probe of the same large resistance. Both ports act on the vertical `E_z`
///    edge between trace and ground at their cell — the natural microstrip
///    feed orientation. A *large* series resistance keeps the loaded Q high so
///    the two split resonances stay sharp.
/// 3. **Time-step** with a custom step body (`update_h_only` → `apply_cpml_h`
///    → port `correct_e` → `update_e_only` → `apply_cpml_e` → `advance_clock`),
///    recording the probe-port `E_z` each step.
/// 4. **DFT** the probe time series with a single-bin (Goertzel-style)
///    frequency scan over `[f0·(1 − span), f0·(1 + span)]` — the same idiom as
///    `yee-fdtd`'s `cavity_resonance.rs` / `ntff.rs` accumulators; no FFT
///    dependency.
/// 5. **Extract** the two split peaks and invert them to `k` via
///    [`yee_filter::extract_coupling`].
///
/// # Performance
///
/// The FDTD run is multi-minute even on the coarse default grid; do **not**
/// call this on a constrained machine. Its validation gate
/// (`fdtd-coupling-001`) is `#[ignore]`'d and runs only in a dedicated CI
/// `--release` job (ADR-0108).
///
/// # Panics
///
/// Panics if the layout has fewer than two ports (a coupled-pair drive needs a
/// drive port and a probe port), or if [`extract_coupling`] cannot find two
/// distinct split peaks in the scanned spectrum (the run did not resolve the
/// resonance split — usually too coarse a grid, too few steps, or a scan window
/// that does not bracket both peaks).
pub fn run_coupled_pair(layout: &Layout, cfg: &CoupledRunConfig) -> CoupledRunResult {
    assert!(
        layout.ports.len() >= 2,
        "run_coupled_pair: need ≥ 2 ports (a drive port and a probe port); got {}",
        layout.ports.len()
    );

    // --- 1. Voxelize: layout -> material-assigned grid + port cells. --------
    let opts = VoxelOptions {
        dx_m: cfg.dx_m,
        xy_margin_cells: cfg.xy_margin_cells,
        air_above_cells: cfg.air_above_cells,
    };
    let model = voxelize_microstrip(layout, &opts);
    let drive_cell = model.port_cells[0];
    let probe_cell = model.port_cells[1];

    let dt = model.grid.dt;
    let mut solver = WalkingSkeletonSolver::new(model.grid);

    // --- 2. Ports: a weakly-coupled Gaussian drive + a passive probe. -------
    //
    // The drive carries a broadband Gaussian-modulated pulse centred on f0 with
    // a spectral FWHM covering the whole scan window, so it excites both split
    // resonances. The probe is the same large-R port with no EMF — its
    // `correct_e` still applies the resistive back-reaction (a weak load), and
    // we read the probe-cell `E_z` directly from the grid each step.
    let bw = 2.0 * cfg.freq_span * cfg.f0_hz;
    let drive_wave = SourceWaveform::GaussianPulse {
        v0: cfg.drive_v0,
        f0: cfg.f0_hz,
        bw,
        // Centre the pulse a few characteristic times in so its t=0 tail is
        // negligible. The Gaussian time constant is τ = √(2 ln2)/(π·bw); place
        // t0 at ~3.5 τ expressed in steps.
        t0_steps: ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw))
            / dt)
            .ceil() as usize,
    };
    let mut drive_port =
        LumpedRlcPort::pure_resistor(drive_cell, cfg.port_resistance_ohm, drive_wave);
    let mut probe_port =
        LumpedRlcPort::pure_resistor(probe_cell, cfg.port_resistance_ohm, SourceWaveform::None);

    // --- 3. Time-step with a custom body; record the probe E_z. -------------
    //
    // The body mirrors the cavity_resonance.rs custom step: H half-step + PEC
    // outer-wall clamp (the no-CPML fall-through of `apply_cpml_h`), then the
    // lumped-port corrections between H and E, then the E half-step + clamp,
    // then advance the clock. The port `correct_e` runs after `update_e_only`
    // (it overwrites the standard Yee E_z estimate at the port cell), matching
    // its documented call site.
    let mut probe_series: Vec<f64> = Vec::with_capacity(cfg.n_steps);
    for n in 0..cfg.n_steps {
        solver.update_h_only();
        solver.apply_cpml_h();

        solver.update_e_only();
        drive_port.correct_e(solver.grid_mut(), n, dt);
        probe_port.correct_e(solver.grid_mut(), n, dt);
        solver.apply_cpml_e();

        solver.advance_clock();

        let ez = solver.grid().ez[probe_cell];
        probe_series.push(ez);
    }

    // --- 4. Single-bin DFT scan over the split-resonance window. ------------
    let f_lo = cfg.f0_hz * (1.0 - cfg.freq_span);
    let f_hi = cfg.f0_hz * (1.0 + cfg.freq_span);
    let n_bins = cfg.n_freq_bins.max(5);
    let df = (f_hi - f_lo) / (n_bins - 1) as f64;

    let mut freqs = Vec::with_capacity(n_bins);
    let mut mag = Vec::with_capacity(n_bins);
    for bin in 0..n_bins {
        let f = f_lo + bin as f64 * df;
        let omega = 2.0 * std::f64::consts::PI * f;
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;
        for (m, &x) in probe_series.iter().enumerate() {
            let phase = omega * m as f64 * dt;
            re += x * phase.cos();
            im -= x * phase.sin();
        }
        freqs.push(f);
        mag.push((re * re + im * im).sqrt());
    }

    // --- 5. Invert the two split peaks to k. --------------------------------
    let extraction = extract_coupling(&freqs, &mag).expect(
        "run_coupled_pair: extract_coupling found < 2 split peaks in the scanned spectrum; \
         the FDTD run did not resolve the resonance split (grid too coarse, too few steps, \
         or scan window does not bracket both peaks)",
    );

    CoupledRunResult {
        k: extraction.k,
        f_even: extraction.f_lo_hz,
        f_odd: extraction.f_hi_hz,
    }
}
