//! Transport-agnostic simulation job API over the `yee-compute` engine
//! (phase S.0, ADR-0179).
//!
//! One serde protocol, many transports: the Tauri studio calls
//! [`submit`] in-process; `yee-server` (S.1) will forward the same
//! [`JobSpec`] / [`JobEvent`] types over WebSocket; the CLI can print the
//! event stream. S.0 covers driven FDTD jobs (any boundary, soft sources +
//! resistive ports, probes) with progress streaming and cooperative
//! cancellation; S.5 adds per-cell materials + interior PEC masks and an
//! explicit dt ([`MaterialsSpec`], [`JobSpec::dt_s`]), so a voxelized
//! layout runs over the same protocol. Dispersive-material specs are a
//! later slice.
//!
//! ```
//! use yee_engine::{BackendChoice, BoundarySpec, JobEvent, JobSpec, SourceSpec, ProbeSpec};
//!
//! let spec = JobSpec {
//!     nx: 12, ny: 12, nz: 12, dx_m: 1e-3, n_steps: 40,
//!     boundary: BoundarySpec::Pec,
//!     sources: vec![SourceSpec::GaussianEz { cell: (6, 6, 6), t0_steps: 8.0, sigma_steps: 3.0 }],
//!     ports: vec![], aperture_ports: vec![], thin_wires: vec![],
//!     probes: vec![ProbeSpec { component: "ez".into(), cell: (8, 6, 6) }],
//!     slice: None, ntff: None, materials: None, dt_s: None, spacings: None,
//!     backend: BackendChoice::Cpu,
//! };
//! let handle = yee_engine::submit(spec);
//! let mut done = false;
//! for event in handle.events() {
//!     if let JobEvent::Done { result } = event {
//!         assert_eq!(result.probes[0].len(), 40);
//!         done = true;
//!     }
//! }
//! assert!(done);
//! ```

pub mod automesh;
pub mod board;
pub mod farfield;
pub mod sparams;

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::thread;

use serde::{Deserialize, Serialize};

use yee_compute::{
    AperturePort, Boundary, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials,
    Probe, ResistivePort, SoftSource, ThinWire, Waveform,
};

/// Which backend executes the job.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BackendChoice {
    /// Rayon FP64 CPU backend.
    Cpu,
    /// wgpu FP32 GPU backend (errors if no adapter).
    Gpu,
    /// GPU with silent CPU fallback.
    Auto,
}

/// Outer-boundary selection.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum BoundarySpec {
    /// Raw kernels, no boundary phase.
    None,
    /// Reflecting PEC box.
    Pec,
    /// Roden–Gedney CPML with `npml` layers on the enabled axes.
    Cpml {
        /// PML thickness in cells.
        npml: usize,
        /// Per-axis enable `[x, y, z]` (S.9, ADR-0186); a disabled axis
        /// keeps its faces PEC. Defaults to all-on (the pre-S.9 wire
        /// format, with no `axes` key, still deserializes). Board-level
        /// verify wants `[true, true, false]`: absorbing side walls, PEC
        /// ground/lid — a thin substrate stack otherwise sits *inside*
        /// the z-min absorber and the line mode is eaten (ADR-0185/0186).
        #[serde(default = "all_axes")]
        axes: [bool; 3],
        /// Optional per-face enable `[[x−, x+], [y−, y+], [z−, z+]]`
        /// (A.2, ADR-0192) — overrides `axes` when present. A disabled
        /// face stays PEC; an antenna's open top over a PEC ground is
        /// `[[true, true], [true, true], [false, true]]`. Both backends
        /// honor arbitrary face masks (GPU since R.3, ADR-0196).
        #[serde(default)]
        faces: Option<[[bool; 2]; 3]>,
    },
}

/// The pre-S.9 default: CPML on every axis.
fn all_axes() -> [bool; 3] {
    [true; 3]
}

/// Source description (times in steps so the spec is dt-agnostic).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SourceSpec {
    /// Soft Gaussian on `E_z`: `exp(−((t−t0)/σ)²)`.
    GaussianEz {
        /// Source cell.
        cell: (usize, usize, usize),
        /// Pulse centre, in steps.
        t0_steps: f64,
        /// Pulse width, in steps.
        sigma_steps: f64,
    },
    /// Soft **modulated** Gaussian on `E_z` (FS.2b): zero-DC drive for
    /// absolute far-field references — a baseband Gaussian's near-DC
    /// content rings in a CPML box far past any practical window and
    /// leaks into single-bin far-field DFTs (measured: ±40 % scatter in
    /// the engine-scale-001 Hertzian ratios before this variant).
    GaussianPulseEz {
        /// Source cell.
        cell: (usize, usize, usize),
        /// Peak amplitude (field units per step).
        v0: f64,
        /// Carrier (Hz).
        f0_hz: f64,
        /// FWHM bandwidth (Hz).
        bw_hz: f64,
        /// Pulse centre, in steps.
        t0_steps: usize,
    },
}

/// Per-axis nonuniform **primal cell widths** (FS.0b.0, ADR-0208) — the
/// serde mirror of `yee_compute::GradedSpacings`. Lengths must equal
/// `nx`/`ny`/`nz` (metres); every width must be positive and finite.
/// `None` on [`JobSpec::spacings`] is the uniform `dx_m` grid, and the
/// pre-FS.0b wire format (no `spacings` key) still deserializes.
///
/// Walking-skeleton scope: **CPU-only** (`backend: "gpu"` is rejected,
/// `"auto"` falls back to CPU); spacing must be uniform inside the CPML
/// layers of every absorbing face; `ntff` is not supported on a graded
/// grid. `dx_m` stays the nominal spacing feeding the CPML σ_max recipe —
/// keep the absorbing-layer spacing equal to it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GradedSpacings {
    /// Primal cell widths along x (length `nx`, metres).
    pub dx: Vec<f64>,
    /// Primal cell widths along y (length `ny`, metres).
    pub dy: Vec<f64>,
    /// Primal cell widths along z (length `nz`, metres).
    pub dz: Vec<f64>,
}

impl GradedSpacings {
    /// Translate to the compute-layer type.
    fn to_compute(&self) -> yee_compute::GradedSpacings {
        yee_compute::GradedSpacings {
            dx: self.dx.clone(),
            dy: self.dy.clone(),
            dz: self.dz.clone(),
        }
    }
}

/// Resistive drive port description.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PortSpec {
    /// Port `E_z` cell.
    pub cell: (usize, usize, usize),
    /// Resistance (Ω).
    pub resistance_ohm: f64,
    /// Peak EMF (V).
    pub v0: f64,
    /// Carrier (Hz).
    pub f0_hz: f64,
    /// FWHM bandwidth (Hz).
    pub bw_hz: f64,
    /// Pulse centre, in steps.
    pub t0_steps: usize,
}

/// Per-cell material maps and interior PEC masks (S.5): the serde mirror
/// of [`yee_compute::Materials`], so a voxelized layout (e.g.
/// `yee_voxel::voxelize_microstrip`) can travel over the job protocol.
///
/// Conventions match `yee_fdtd::grid::YeeGrid`: the ε_r / μ_r / σ maps are
/// `[nx+1, ny+1, nz+1]` row-major; each PEC mask has its E component's
/// staggered shape. Lengths are validated at submission — a mismatch is a
/// [`JobEvent::Error`], never a panic.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MaterialsSpec {
    /// Per-cell relative permittivity, `[nx+1, ny+1, nz+1]` row-major.
    #[serde(default)]
    pub eps_r_cells: Option<Vec<f64>>,
    /// Per-cell relative permeability, `[nx+1, ny+1, nz+1]` row-major.
    #[serde(default)]
    pub mu_r_cells: Option<Vec<f64>>,
    /// Per-cell electric conductivity (S/m), `[nx+1, ny+1, nz+1]` row-major.
    #[serde(default)]
    pub sigma_cells: Option<Vec<f64>>,
    /// Interior PEC mask for `E_x` (shape of `E_x`).
    #[serde(default)]
    pub pec_mask_ex: Option<Vec<bool>>,
    /// Interior PEC mask for `E_y` (shape of `E_y`).
    #[serde(default)]
    pub pec_mask_ey: Option<Vec<bool>>,
    /// Interior PEC mask for `E_z` (shape of `E_z`).
    #[serde(default)]
    pub pec_mask_ez: Option<Vec<bool>>,
    /// Resistive-sheet surface resistance `R_s` (Ω/square) on the masked
    /// `E_x`/`E_y` trace edges (R.0b, ADR-0202): planar conductor loss at
    /// the design frequency, `R_s = √(π f_ref μ₀/σ)`
    /// (`yee_voxel::surface_resistance_ohm`). `None`/`0` → lossless PEC.
    /// CPU-only for now (the GPU backend rejects it; `Auto` falls back).
    #[serde(default)]
    pub sheet_r_ohm: Option<f64>,
}

impl MaterialsSpec {
    /// Check every present map/mask against the shape `fdtd_spec` demands,
    /// then convert. Errors (instead of panicking like
    /// `yee_compute::Materials::validate`) because job specs arrive from
    /// untrusted transports.
    fn into_materials(self, fdtd_spec: &FdtdSpec) -> Result<Materials, String> {
        let cells = (fdtd_spec.nx + 1) * (fdtd_spec.ny + 1) * (fdtd_spec.nz + 1);
        for (map, name) in [
            (&self.eps_r_cells, "eps_r_cells"),
            (&self.mu_r_cells, "mu_r_cells"),
            (&self.sigma_cells, "sigma_cells"),
        ] {
            if let Some(m) = map
                && m.len() != cells
            {
                return Err(format!(
                    "materials.{name} has {} entries, expected {cells} \
                     ([nx+1, ny+1, nz+1] row-major)",
                    m.len()
                ));
            }
        }
        let dims_len = |d: (usize, usize, usize)| d.0 * d.1 * d.2;
        for (mask, len, name) in [
            (
                &self.pec_mask_ex,
                dims_len(fdtd_spec.ex_dims()),
                "pec_mask_ex",
            ),
            (
                &self.pec_mask_ey,
                dims_len(fdtd_spec.ey_dims()),
                "pec_mask_ey",
            ),
            (
                &self.pec_mask_ez,
                dims_len(fdtd_spec.ez_dims()),
                "pec_mask_ez",
            ),
        ] {
            if let Some(m) = mask
                && m.len() != len
            {
                return Err(format!(
                    "materials.{name} has {} entries, expected {len} \
                     (the component's staggered shape)",
                    m.len()
                ));
            }
        }
        if let Some(r) = self.sheet_r_ohm
            && !(r.is_finite() && r >= 0.0)
        {
            return Err(format!(
                "materials.sheet_r_ohm must be finite and non-negative, got {r}"
            ));
        }
        Ok(Materials {
            eps_r_cells: self.eps_r_cells,
            mu_r_cells: self.mu_r_cells,
            sigma_cells: self.sigma_cells,
            pec_mask_ex: self.pec_mask_ex,
            pec_mask_ey: self.pec_mask_ey,
            pec_mask_ez: self.pec_mask_ez,
            sheet_r_ohm: self.sheet_r_ohm,
        })
    }
}

/// Multi-cell **aperture port** description (S.10, ADR-0187): one
/// aggregate resistive branch bridging the modal port face — the `(y, z)`
/// band `j in [j_lo, j_hi)` x substrate column `k in [0, k_top)` at plane
/// `i`. The engine derives the modal geometry from the grid
/// (`h = k_top*dx`, `A = (j_hi - j_lo)*dx*h`) and runs the validated
/// `yee_fdtd::LumpedRlcPort::aperture` scheme (ported bit-exact, gate
/// compute-014). This is the dx-stable port a single-cell [`PortSpec`]
/// cannot approximate on a multi-cell substrate; a naive series stack of
/// single-cell ports was measured worse and rejected (ADR-0187). Runs on
/// both backends (GPU since R.3, ADR-0196; gate compute-015).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AperturePortSpec {
    /// Port-plane x-index shared by every aperture cell.
    pub i: usize,
    /// Width band start (inclusive), `E_z` y-index.
    pub j_lo: usize,
    /// Width band end (exclusive).
    pub j_hi: usize,
    /// Ground level in cells (FS.1a.1b, ADR-0205): the drive/measure
    /// column spans `E_z` edges `k = k_lo .. k_top`. `0` (the serde
    /// default — full wire back-compat) is the classic ground-on-the-
    /// domain-floor stack; a lifted stack (`voxelize_microstrip_open`)
    /// passes its `k_gnd`.
    #[serde(default)]
    pub k_lo: usize,
    /// Trace level in cells: the column's top `E_z` edge (exclusive
    /// bound is the trace plane; ground sheet at `k = k_lo`).
    pub k_top: usize,
    /// Aggregate branch resistance (Ohm).
    pub resistance_ohm: f64,
    /// Peak EMF (V); `0.0` = passive matched load.
    pub v0: f64,
    /// Carrier (Hz).
    pub f0_hz: f64,
    /// FWHM bandwidth (Hz).
    pub bw_hz: f64,
    /// Pulse centre, in steps.
    pub t0_steps: usize,
    /// Record the per-step `(v_terminal, i_branch)` pair for this port
    /// (FS.2a, ADR-0207) — the accepted-power observables, returned in
    /// [`JobResult::port_records`]. CPU-only; the GPU backend rejects
    /// recording ports with `Unsupported` (Auto falls back).
    #[serde(default)]
    pub record: bool,
}

/// Z-directed thin-wire subcell (FS.1c, ADR-0228): the **Holland & Simpson
/// (1981) in-cell-inductance** model for a wire much thinner than a cell —
/// see [`yee_compute::ThinWire`]'s doc for the full citation and the
/// discrete update it implements. A wire spans `E_z` cells `k ∈ [k_lo,
/// k_hi)` at the fixed `(i, j)` column; `feed_k`, if set, is a literal gap
/// in the wire (no metal, no in-cell inductance) left for a normal
/// [`AperturePortSpec`] (single-cell, `j_lo = j`, `j_hi = j + 1`, `k_lo =
/// feed_k`, `k_top = feed_k + 1`) to drive — the mom-001 delta-gap idiom.
/// CPU backend only: the GPU backend rejects a job whose drive carries any
/// thin wires with a named `Unsupported` error (`Auto` falls back to CPU).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThinWireSpec {
    /// Wire's fixed `x` grid index (the `E_z` column it threads).
    pub i: usize,
    /// Wire's fixed `y` grid index (the `E_z` column it threads).
    pub j: usize,
    /// First `E_z` cell index (`k`) the wire occupies (inclusive).
    pub k_lo: usize,
    /// One past the last `E_z` cell index the wire occupies (exclusive).
    pub k_hi: usize,
    /// Physical wire radius (m); must be well under half the local
    /// transverse cell size (`dx`/`dy`, assumed square).
    pub radius_m: f64,
    /// Optional delta-gap feed cell in `[k_lo, k_hi)`.
    #[serde(default)]
    pub feed_k: Option<usize>,
}

/// Far-field request (A.2, ADR-0192): accumulate a surface DFT at `f_hz`
/// during the run (the validated `yee_fdtd::NtffState` transform) and
/// return `|E|` for each requested `(θ, φ)` direction (radians, θ from
/// +z). The integration box sits `margin_cells` inside every face except
/// optionally the bottom: `k_min` overrides the bottom face for grounded
/// antennas (e.g. `k_min = 1` hugs the ground plane; documented
/// equivalence-surface approximation — the bottom face crosses the
/// substrate). **CPU-only**: `backend: "gpu"` fails with an error event,
/// `"auto"` falls back to CPU.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NtffSpec {
    /// DFT frequency (Hz).
    pub f_hz: f64,
    /// Box margin from every outer face, in cells (keep ≥ npml + a few).
    pub margin_cells: usize,
    /// Optional bottom-face override (`k` index), for grounded antennas.
    #[serde(default)]
    pub k_min: Option<usize>,
    /// Far-field directions `(θ, φ)` in radians.
    pub directions: Vec<(f64, f64)>,
}

/// Field-slice request: one z-plane of a component, returned with the
/// result for visualization (S.3).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SliceSpec {
    /// `"ex"`, `"ey"`, or `"ez"`.
    pub component: String,
    /// z-index (`k`) of the plane, in the component's staggered indexing.
    pub k: usize,
}

/// A field slice payload (row-major `[ni, nj]`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldSlice {
    /// First dimension (i extent).
    pub ni: usize,
    /// Second dimension (j extent).
    pub nj: usize,
    /// Row-major values.
    pub data: Vec<f64>,
}

/// Probe description.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProbeSpec {
    /// `"ex"`, `"ey"`, or `"ez"`.
    pub component: String,
    /// Probed cell.
    pub cell: (usize, usize, usize),
}

/// A simulation job: uniform-vacuum driven FDTD (S.0 scope).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JobSpec {
    /// Cells along x.
    pub nx: usize,
    /// Cells along y.
    pub ny: usize,
    /// Cells along z.
    pub nz: usize,
    /// Cubic cell size (m).
    pub dx_m: f64,
    /// Steps to run.
    pub n_steps: usize,
    /// Outer boundary.
    pub boundary: BoundarySpec,
    /// Soft sources.
    pub sources: Vec<SourceSpec>,
    /// Resistive ports.
    pub ports: Vec<PortSpec>,
    /// Multi-cell aperture ports (S.10; both backends since R.3).
    #[serde(default)]
    pub aperture_ports: Vec<AperturePortSpec>,
    /// Z-directed thin-wire subcells (FS.1c; CPU backend only — the GPU
    /// backend rejects a drive that carries any, `Auto` falls back).
    #[serde(default)]
    pub thin_wires: Vec<ThinWireSpec>,
    /// Probes (recorded every step).
    pub probes: Vec<ProbeSpec>,
    /// Optional final-field slice to return for visualization.
    #[serde(default)]
    pub slice: Option<SliceSpec>,
    /// Optional far-field accumulation (A.2; CPU-only).
    #[serde(default)]
    pub ntff: Option<NtffSpec>,
    /// Optional per-cell materials + interior PEC masks (S.5). `None` is
    /// the S.0 uniform-vacuum behaviour.
    #[serde(default)]
    pub materials: Option<MaterialsSpec>,
    /// Optional explicit time step (s), overriding the default 0.9×Courant
    /// step — e.g. the dt a voxelizer computed for its grid. Must be
    /// positive and at most the Courant limit for `dx_m`.
    #[serde(default)]
    pub dt_s: Option<f64>,
    /// Optional per-axis graded primal spacings (FS.0b.0/0b.2 — see
    /// [`GradedSpacings`]; runs on both backends since ADR-0214, with
    /// NTFF-on-graded and z-taper-straddling aperture ports rejected).
    /// `None` is the uniform `dx_m` grid. When set,
    /// the default time step comes from the **minimum** spacing per axis
    /// (0.9× the graded Courant limit), and `dt_s` is validated against
    /// that limit.
    #[serde(default)]
    pub spacings: Option<GradedSpacings>,
    /// Backend selection.
    pub backend: BackendChoice,
}

/// Result of a finished job.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JobResult {
    /// Backend that actually ran (`"cpu"` / `"gpu"`).
    pub backend: String,
    /// Time step (s) — pairs with per-step probe samples.
    pub dt_s: f64,
    /// One series per probe, one sample per completed step.
    pub probes: Vec<Vec<f64>>,
    /// The requested final-field slice, if any.
    #[serde(default)]
    pub slice: Option<FieldSlice>,
    /// `|E_far|` per requested [`NtffSpec::directions`] entry, if any.
    #[serde(default)]
    pub far_field: Option<Vec<f64>>,
    /// Per-step `(v_src, v_terminal, i_branch)` per aperture port,
    /// present when any port set [`AperturePortSpec::record`] (FS.2a;
    /// empty inner vec for non-recording ports).
    #[serde(default)]
    pub port_records: Option<Vec<Vec<(f64, f64, f64)>>>,
    /// Steps completed (equals the request unless cancelled).
    pub steps_done: usize,
}

/// Progress / completion events streamed while a job runs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "event", rename_all = "snake_case")]
pub enum JobEvent {
    /// Periodic progress (~every 2 % of the run).
    Progress {
        /// Steps completed so far.
        step: usize,
        /// Total steps requested.
        total: usize,
    },
    /// The job finished (or was cancelled after a whole number of chunks).
    Done {
        /// The result payload.
        result: JobResult,
    },
    /// The job failed.
    Error {
        /// Human-readable failure.
        message: String,
    },
}

/// Handle to a submitted job: consume [`JobHandle::events`], call
/// [`JobHandle::cancel`] to stop after the current chunk.
pub struct JobHandle {
    events: Receiver<JobEvent>,
    cancel: Arc<AtomicBool>,
}

impl JobHandle {
    /// The event stream. Iterate to completion; the final event is always
    /// `Done` or `Error`.
    pub fn events(&self) -> impl Iterator<Item = JobEvent> + '_ {
        self.events.iter()
    }

    /// Request cooperative cancellation; the job emits `Done` with the
    /// steps completed so far.
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    /// A detachable canceller for use when the handle itself moves to a
    /// worker (e.g. yee-server cancels on client disconnect).
    pub fn canceller(&self) -> JobCanceller {
        JobCanceller {
            cancel: Arc::clone(&self.cancel),
        }
    }
}

/// Cloneable cancellation handle detached from a [`JobHandle`].
#[derive(Clone)]
pub struct JobCanceller {
    cancel: Arc<AtomicBool>,
}

impl JobCanceller {
    /// Request cooperative cancellation.
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

/// Validate a spec and spawn it on a worker thread.
pub fn submit(spec: JobSpec) -> JobHandle {
    let (tx, rx) = channel();
    let cancel = Arc::new(AtomicBool::new(false));
    let flag = Arc::clone(&cancel);
    thread::spawn(move || run_job(spec, &tx, &flag));
    JobHandle { events: rx, cancel }
}

fn build_drive(
    spec: &JobSpec,
    fdtd_spec: &FdtdSpec,
    dt: f64,
    graded: Option<&yee_compute::GradedSpacings>,
) -> Result<Drive, String> {
    let mut drive = Drive::default();
    for s in &spec.sources {
        match *s {
            SourceSpec::GaussianEz {
                cell,
                t0_steps,
                sigma_steps,
            } => drive.soft_sources.push(SoftSource {
                component: EComponent::Ez,
                cell,
                waveform: Waveform::Gaussian {
                    t0: t0_steps * dt,
                    sigma: sigma_steps * dt,
                },
            }),
            SourceSpec::GaussianPulseEz {
                cell,
                v0,
                f0_hz,
                bw_hz,
                t0_steps,
            } => drive.soft_sources.push(SoftSource {
                component: EComponent::Ez,
                cell,
                waveform: Waveform::GaussianPulse {
                    v0,
                    f0: f0_hz,
                    bw: bw_hz,
                    t0_steps,
                },
            }),
        }
    }
    for p in &spec.ports {
        drive.ports.push(ResistivePort {
            cell: p.cell,
            resistance: p.resistance_ohm,
            waveform: Waveform::GaussianPulse {
                v0: p.v0,
                f0: p.f0_hz,
                bw: p.bw_hz,
                t0_steps: p.t0_steps,
            },
        });
    }
    let ez_dims = fdtd_spec.ez_dims();
    for (n, a) in spec.aperture_ports.iter().enumerate() {
        if a.j_hi <= a.j_lo || a.k_top == 0 {
            return Err(format!(
                "aperture_ports[{n}]: empty aperture (j band [{}, {}), k_top {})",
                a.j_lo, a.j_hi, a.k_top
            ));
        }
        if a.i >= ez_dims.0 || a.j_hi > ez_dims.1 || a.k_top > ez_dims.2 {
            return Err(format!(
                "aperture_ports[{n}]: outside the E_z grid {ez_dims:?}"
            ));
        }
        if a.k_lo >= a.k_top {
            return Err(format!(
                "aperture_ports[{n}]: k_lo {} must sit below k_top {}",
                a.k_lo, a.k_top
            ));
        }
        if !(a.resistance_ohm > 0.0 && a.resistance_ohm.is_finite()) {
            return Err(format!(
                "aperture_ports[{n}]: resistance must be positive and finite (got {})",
                a.resistance_ohm
            ));
        }
        let cells: Vec<(usize, usize, usize)> = (a.j_lo..a.j_hi)
            .flat_map(|j| (a.k_lo..a.k_top).map(move |k| (a.i, j, k)))
            .collect();
        let n_columns = a.j_hi - a.j_lo;
        // Physical aperture geometry (FS.0b.0 × FS.1a.1b): the column
        // spans k_lo..k_top — spacing sums on a graded grid, scalar
        // products on the uniform one.
        let (height, area) = match graded {
            None => {
                let height = (a.k_top - a.k_lo) as f64 * fdtd_spec.dz;
                (height, n_columns as f64 * fdtd_spec.dy * height)
            }
            Some(g) => {
                if a.j_hi > g.dy.len() {
                    return Err(format!(
                        "aperture_ports[{n}]: j band [{}, {}) exceeds the {} primal \
                         y cells on a graded grid",
                        a.j_lo,
                        a.j_hi,
                        g.dy.len()
                    ));
                }
                let height: f64 = g.dz[a.k_lo..a.k_top].iter().sum();
                let width: f64 = g.dy[a.j_lo..a.j_hi].iter().sum();
                (height, width * height)
            }
        };
        drive.aperture_ports.push(AperturePort {
            cells,
            n_columns,
            area,
            height,
            resistance: a.resistance_ohm,
            waveform: Waveform::GaussianPulse {
                v0: a.v0,
                f0: a.f0_hz,
                bw: a.bw_hz,
                t0_steps: a.t0_steps,
            },
            record: a.record,
        });
    }
    for w in &spec.thin_wires {
        drive.thin_wires.push(ThinWire {
            i: w.i,
            j: w.j,
            k_lo: w.k_lo,
            k_hi: w.k_hi,
            radius_m: w.radius_m,
            feed_k: w.feed_k,
        });
    }
    for p in &spec.probes {
        let component = match p.component.as_str() {
            "ex" => EComponent::Ex,
            "ey" => EComponent::Ey,
            "ez" => EComponent::Ez,
            other => return Err(format!("unknown probe component {other:?}")),
        };
        drive.probes.push(Probe {
            component,
            cell: p.cell,
        });
    }
    Ok(drive)
}

/// Extract the requested z-plane from the final fields.
fn extract_slice(
    spec: &JobSpec,
    fdtd_spec: &FdtdSpec,
    fields: &yee_compute::Fields,
) -> Result<Option<FieldSlice>, String> {
    let Some(s) = &spec.slice else {
        return Ok(None);
    };
    let (dims, data): ((usize, usize, usize), &[f64]) = match s.component.as_str() {
        "ex" => (fdtd_spec.ex_dims(), &fields.ex),
        "ey" => (fdtd_spec.ey_dims(), &fields.ey),
        "ez" => (fdtd_spec.ez_dims(), &fields.ez),
        other => return Err(format!("unknown slice component {other:?}")),
    };
    if s.k >= dims.2 {
        return Err(format!("slice k = {} out of range (< {})", s.k, dims.2));
    }
    let mut out = Vec::with_capacity(dims.0 * dims.1);
    for i in 0..dims.0 {
        for j in 0..dims.1 {
            out.push(data[(i * dims.1 + j) * dims.2 + s.k]);
        }
    }
    Ok(Some(FieldSlice {
        ni: dims.0,
        nj: dims.1,
        data: out,
    }))
}

/// Chunked execution with progress events (~2 % granularity, min 1 step).
fn run_job(mut spec: JobSpec, tx: &Sender<JobEvent>, cancel: &AtomicBool) {
    // Specs arrive over untrusted transports (yee-server): every malformed
    // field must come back as an Error event, not a worker-thread panic
    // (a panic closes the channel with no terminal event).
    if spec.nx == 0 || spec.ny == 0 || spec.nz == 0 || !(spec.dx_m.is_finite() && spec.dx_m > 0.0) {
        let _ = tx.send(JobEvent::Error {
            message: format!(
                "invalid grid: {}x{}x{} cells at dx = {} m",
                spec.nx, spec.ny, spec.nz, spec.dx_m
            ),
        });
        return;
    }
    let mut fdtd_spec = FdtdSpec::vacuum(spec.nx, spec.ny, spec.nz, spec.dx_m);
    // Graded spacings (FS.0b.0): validate, then retime the default dt from
    // the minimum spacing per axis (same 0.9× Courant recipe as the scalar
    // path — bit-identical for constant arrays).
    let graded = spec.spacings.as_ref().map(GradedSpacings::to_compute);
    if let Some(g) = &graded {
        if let Err(message) = g.validate(&fdtd_spec) {
            let _ = tx.send(JobEvent::Error { message });
            return;
        }
        fdtd_spec.dt = 0.9 * g.courant_limit();
        if spec.ntff.is_some() {
            let _ = tx.send(JobEvent::Error {
                message: "ntff far-field on a graded grid is not supported (FS.0b.0)".into(),
            });
            return;
        }
    }
    if let Some(dt) = spec.dt_s {
        let limit = graded
            .as_ref()
            .map_or_else(|| fdtd_spec.courant_limit(), |g| g.courant_limit());
        if !(dt.is_finite() && dt > 0.0 && dt <= limit) {
            let _ = tx.send(JobEvent::Error {
                message: format!("dt_s = {dt} s is outside (0, {limit:.6e}] (Courant limit)"),
            });
            return;
        }
        fdtd_spec.dt = dt;
    }
    let materials = match spec.materials.take() {
        Some(m) => match m.into_materials(&fdtd_spec) {
            Ok(m) => m,
            Err(message) => {
                let _ = tx.send(JobEvent::Error { message });
                return;
            }
        },
        None => Materials::default(),
    };
    // NTFF (A.2): validate the box before running; CPU-only.
    let ntff_bounds = match &spec.ntff {
        None => None,
        Some(n) => {
            let m = n.margin_cells;
            let k0 = n.k_min.unwrap_or(m);
            let bounds = (
                m,
                spec.nx.saturating_sub(m),
                m,
                spec.ny.saturating_sub(m),
                k0,
                spec.nz.saturating_sub(m),
            );
            if bounds.0 >= bounds.1 || bounds.2 >= bounds.3 || bounds.4 >= bounds.5 {
                let _ = tx.send(JobEvent::Error {
                    message: format!(
                        "ntff box {bounds:?} has non-positive extent on grid \
                         {}x{}x{}",
                        spec.nx, spec.ny, spec.nz
                    ),
                });
                return;
            }
            if !(n.f_hz.is_finite() && n.f_hz > 0.0) || n.directions.is_empty() {
                let _ = tx.send(JobEvent::Error {
                    message: "ntff needs a positive f_hz and at least one direction".into(),
                });
                return;
            }
            Some(bounds)
        }
    };
    let boundary = match spec.boundary {
        BoundarySpec::None => Boundary::None,
        BoundarySpec::Pec => Boundary::PecBox,
        BoundarySpec::Cpml { npml, axes, faces } => {
            let mut config = CpmlConfig::for_spec(&fdtd_spec, npml).with_axes(axes);
            if let Some(f) = faces {
                config = config.with_faces(f);
            }
            // FS.0b.0 scope: spacing must be uniform inside the CPML layers
            // of every absorbing face.
            if let Some(g) = &graded
                && let Err(message) = g.validate_cpml_layers(npml, config.faces)
            {
                let _ = tx.send(JobEvent::Error { message });
                return;
            }
            Boundary::Cpml(config)
        }
    };
    let drive = match build_drive(&spec, &fdtd_spec, fdtd_spec.dt, graded.as_ref()) {
        Ok(d) => d,
        Err(message) => {
            let _ = tx.send(JobEvent::Error { message });
            return;
        }
    };

    let chunk = (spec.n_steps / 50).max(1);
    let total = spec.n_steps;
    let dt_s = fdtd_spec.dt;

    // GPU path (or auto → GPU when an adapter exists).
    #[cfg(feature = "gpu")]
    if matches!(spec.backend, BackendChoice::Gpu) && spec.ntff.is_some() {
        // Auto falls through to the CPU path; only an explicit GPU request
        // fails.
        let _ = tx.send(JobEvent::Error {
            message: "ntff far-field is CPU-only (A.2, ADR-0192); use backend \"cpu\" or \"auto\""
                .into(),
        });
        return;
    }
    // Graded grids run on the GPU as of FS.0b.2 (ADR-0214): the kernels
    // multiply by per-cell inverse spacings; set_spacings attaches them
    // after construction (aperture-port z-taper and NTFF-on-graded are
    // rejected inside set_spacings with Unsupported — Auto falls back).
    #[cfg(feature = "gpu")]
    if matches!(spec.backend, BackendChoice::Gpu | BackendChoice::Auto) && spec.ntff.is_none() {
        let gpu_engine = yee_compute::GpuFdtd::with_drive(
            fdtd_spec,
            Fields::zero(&fdtd_spec),
            materials.clone(),
            boundary.clone(),
            drive.clone(),
            total,
        )
        .and_then(|mut engine| {
            if let Some(g) = &graded {
                engine.set_spacings(g)?;
            }
            Ok(engine)
        });
        match gpu_engine {
            Ok(mut engine) => {
                let mut done = 0usize;
                while done < total && !cancel.load(Ordering::Relaxed) {
                    let n = chunk.min(total - done);
                    if let Err(e) = engine.step_n(n) {
                        let _ = tx.send(JobEvent::Error {
                            message: e.to_string(),
                        });
                        return;
                    }
                    done += n;
                    let _ = tx.send(JobEvent::Progress { step: done, total });
                }
                let slice = match engine
                    .read_fields()
                    .map_err(|e| e.to_string())
                    .and_then(|f| extract_slice(&spec, &fdtd_spec, &f))
                {
                    Ok(s) => s,
                    Err(message) => {
                        let _ = tx.send(JobEvent::Error { message });
                        return;
                    }
                };
                match engine.read_probes() {
                    Ok(probes) => {
                        let _ = tx.send(JobEvent::Done {
                            result: JobResult {
                                backend: "gpu".into(),
                                dt_s,
                                probes,
                                slice,
                                far_field: None,
                                port_records: None,
                                steps_done: done,
                            },
                        });
                    }
                    Err(e) => {
                        let _ = tx.send(JobEvent::Error {
                            message: e.to_string(),
                        });
                    }
                }
                return;
            }
            Err(
                yee_compute::ComputeError::NoAdapter | yee_compute::ComputeError::Unsupported(_),
            ) if matches!(spec.backend, BackendChoice::Auto) => {} // fall through to CPU
            Err(e) => {
                let _ = tx.send(JobEvent::Error {
                    message: e.to_string(),
                });
                return;
            }
        }
    }
    #[cfg(not(feature = "gpu"))]
    if matches!(spec.backend, BackendChoice::Gpu) {
        let _ = tx.send(JobEvent::Error {
            message: "yee-engine was built without the `gpu` feature".into(),
        });
        return;
    }

    // CPU path.
    let mut engine = CpuFdtd::with_drive(
        fdtd_spec,
        Fields::zero(&fdtd_spec),
        materials,
        boundary,
        drive,
    );
    // Attach the graded spacings (already validated above, so the panicking
    // compute-layer asserts cannot fire from a malformed job spec).
    if let Some(g) = &graded {
        engine.set_spacings(g);
    }
    // NTFF accumulator (A.2): the validated yee-fdtd transform sampling
    // the engine's fields through a host-side grid adapter each step
    // (the compute-010 pattern). The scratch grid mirrors the run's dt.
    let mut ntff = ntff_bounds.map(|bounds| {
        let n = spec.ntff.as_ref().expect("bounds imply spec");
        let mut scratch = yee_fdtd::YeeGrid::vacuum(spec.nx, spec.ny, spec.nz, spec.dx_m);
        scratch.dt = fdtd_spec.dt;
        let state = yee_fdtd::NtffState::with_bounds(
            &scratch,
            yee_fdtd::NtffParams {
                f_probe: n.f_hz,
                box_margin_cells: n.margin_cells,
                theta_rad: 0.0,
                phi_rad: 0.0,
            },
            bounds,
        );
        (scratch, state)
    });
    let mut done = 0usize;
    while done < total && !cancel.load(Ordering::Relaxed) {
        let n = chunk.min(total - done);
        match &mut ntff {
            None => engine.step_n(n),
            Some((scratch, state)) => {
                for _ in 0..n {
                    engine.step_n(1);
                    let f = engine.fields();
                    scratch.ex.as_slice_mut().unwrap().copy_from_slice(&f.ex);
                    scratch.ey.as_slice_mut().unwrap().copy_from_slice(&f.ey);
                    scratch.ez.as_slice_mut().unwrap().copy_from_slice(&f.ez);
                    scratch.hx.as_slice_mut().unwrap().copy_from_slice(&f.hx);
                    scratch.hy.as_slice_mut().unwrap().copy_from_slice(&f.hy);
                    scratch.hz.as_slice_mut().unwrap().copy_from_slice(&f.hz);
                    state.sample(scratch, engine.current_time());
                }
            }
        }
        done += n;
        let _ = tx.send(JobEvent::Progress { step: done, total });
    }
    let slice = match extract_slice(&spec, &fdtd_spec, engine.fields()) {
        Ok(s) => s,
        Err(message) => {
            let _ = tx.send(JobEvent::Error { message });
            return;
        }
    };
    let far_field = ntff.map(|(_, state)| {
        spec.ntff
            .as_ref()
            .expect("state implies spec")
            .directions
            .iter()
            .map(|&(theta, phi)| state.far_field_at(theta, phi).norm())
            .collect()
    });
    let _ = tx.send(JobEvent::Done {
        result: JobResult {
            backend: "cpu".into(),
            dt_s,
            probes: engine.probe_series().to_vec(),
            slice,
            far_field,
            port_records: if engine.aperture_records().iter().any(|r| !r.is_empty()) {
                Some(engine.aperture_records().to_vec())
            } else {
                None
            },
            steps_done: done,
        },
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cavity_spec(backend: BackendChoice) -> JobSpec {
        JobSpec {
            nx: 12,
            ny: 12,
            nz: 12,
            dx_m: 1e-3,
            n_steps: 60,
            boundary: BoundarySpec::Pec,
            sources: vec![SourceSpec::GaussianEz {
                cell: (6, 6, 6),
                t0_steps: 8.0,
                sigma_steps: 3.0,
            }],
            ports: vec![],
            aperture_ports: vec![],
            thin_wires: vec![],
            probes: vec![ProbeSpec {
                component: "ez".into(),
                cell: (8, 6, 6),
            }],
            slice: None,
            ntff: None,
            materials: None,
            dt_s: None,
            spacings: None,
            backend,
        }
    }

    /// A small heterogeneous scenario (S.5): ε_r block + interior PEC mask
    /// + explicit dt on an 8³ grid.
    fn heterogeneous_spec() -> JobSpec {
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.nx = 8;
        spec.ny = 8;
        spec.nz = 8;
        spec.n_steps = 25;
        spec.sources = vec![SourceSpec::GaussianEz {
            cell: (4, 4, 4),
            t0_steps: 6.0,
            sigma_steps: 2.0,
        }];
        spec.probes = vec![ProbeSpec {
            component: "ez".into(),
            cell: (6, 4, 4),
        }];
        let cells = 9 * 9 * 9;
        let mut eps = vec![1.0; cells];
        for (n, e) in eps.iter_mut().enumerate() {
            if n % 9 < 4 {
                *e = 4.4; // lower-z substrate-ish block (k = innermost index)
            }
        }
        let ez_len = 9 * 9 * 8;
        let mut mask = vec![false; ez_len];
        mask[(2 * 9 + 2) * 8 + 5] = true; // one interior PEC E_z edge
        spec.materials = Some(MaterialsSpec {
            eps_r_cells: Some(eps),
            pec_mask_ez: Some(mask),
            ..MaterialsSpec::default()
        });
        let vacuum = FdtdSpec::vacuum(8, 8, 8, spec.dx_m);
        spec.dt_s = Some(0.8 * vacuum.courant_limit());
        spec
    }

    fn run_to_done(spec: JobSpec) -> JobResult {
        let handle = submit(spec);
        handle
            .events()
            .find_map(|e| match e {
                JobEvent::Done { result } => Some(result),
                JobEvent::Error { message } => panic!("job failed: {message}"),
                _ => None,
            })
            .expect("no Done event")
    }

    fn expect_error(spec: JobSpec) -> String {
        let handle = submit(spec);
        handle
            .events()
            .find_map(|e| match e {
                JobEvent::Error { message } => Some(message),
                JobEvent::Done { .. } => panic!("job unexpectedly succeeded"),
                _ => None,
            })
            .expect("no terminal event")
    }

    #[test]
    fn spec_round_trips_through_json() {
        let spec = cavity_spec(BackendChoice::Auto);
        let json = serde_json::to_string(&spec).unwrap();
        let back: JobSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, back);
    }

    #[test]
    fn materials_spec_round_trips_and_defaults() {
        // Full round trip with materials + dt attached.
        let spec = heterogeneous_spec();
        let json = serde_json::to_string(&spec).unwrap();
        let back: JobSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, back);
        // Pre-S.5 specs (no materials/dt_s keys) still deserialize.
        let legacy = serde_json::to_string(&cavity_spec(BackendChoice::Cpu))
            .unwrap()
            .replace(",\"materials\":null", "")
            .replace(",\"dt_s\":null", "");
        let back: JobSpec = serde_json::from_str(&legacy).unwrap();
        assert!(back.materials.is_none() && back.dt_s.is_none());
    }

    #[test]
    fn cpml_axes_default_and_round_trip() {
        // Pre-S.9 wire format (no `axes` key) still deserializes: all-on.
        let legacy: BoundarySpec = serde_json::from_str(r#"{"kind":"cpml","npml":10}"#).unwrap();
        assert_eq!(
            legacy,
            BoundarySpec::Cpml {
                npml: 10,
                axes: [true; 3],
                faces: None
            }
        );
        // Pre-A.2 wire format (no `faces` key) still parses: None.
        let legacy2: BoundarySpec =
            serde_json::from_str(r#"{"kind":"cpml","npml":10,"axes":[true,true,false]}"#).unwrap();
        assert_eq!(
            legacy2,
            BoundarySpec::Cpml {
                npml: 10,
                axes: [true, true, false],
                faces: None
            }
        );
        // Explicit per-axis enables round-trip.
        let side_walls = BoundarySpec::Cpml {
            npml: 8,
            axes: [true, true, false],
            faces: Some([[true, true], [true, true], [false, true]]),
        };
        let json = serde_json::to_string(&side_walls).unwrap();
        let back: BoundarySpec = serde_json::from_str(&json).unwrap();
        assert_eq!(side_walls, back);
    }

    #[test]
    fn heterogeneous_job_is_bit_exact_vs_direct_engine() {
        let spec = heterogeneous_spec();
        let result = run_to_done(spec.clone());
        assert_eq!(result.dt_s, spec.dt_s.unwrap());

        // The same scenario straight on yee-compute, bypassing the protocol.
        let mut fdtd_spec = FdtdSpec::vacuum(spec.nx, spec.ny, spec.nz, spec.dx_m);
        fdtd_spec.dt = spec.dt_s.unwrap();
        let materials = spec
            .materials
            .clone()
            .unwrap()
            .into_materials(&fdtd_spec)
            .unwrap();
        let drive = build_drive(&spec, &fdtd_spec, fdtd_spec.dt, None).unwrap();
        let mut direct = CpuFdtd::with_drive(
            fdtd_spec,
            Fields::zero(&fdtd_spec),
            materials,
            Boundary::PecBox,
            drive,
        );
        direct.step_n(spec.n_steps);

        let reference = &direct.probe_series()[0];
        assert_eq!(result.probes[0].len(), reference.len());
        assert!(reference.iter().any(|v| *v != 0.0), "probe stayed silent");
        for (a, b) in result.probes[0].iter().zip(reference) {
            assert!(a == b, "protocol path diverged from direct engine");
        }
    }

    #[test]
    fn malformed_specs_error_instead_of_panicking() {
        // Wrong-length eps map.
        let mut spec = heterogeneous_spec();
        spec.materials.as_mut().unwrap().eps_r_cells = Some(vec![1.0; 7]);
        assert!(expect_error(spec).contains("eps_r_cells"));
        // Wrong-length PEC mask.
        let mut spec = heterogeneous_spec();
        spec.materials.as_mut().unwrap().pec_mask_ez = Some(vec![false; 3]);
        assert!(expect_error(spec).contains("pec_mask_ez"));
        // Courant-violating dt.
        let mut spec = heterogeneous_spec();
        spec.dt_s = Some(1.0);
        assert!(expect_error(spec).contains("Courant"));
        // Zero-cell grid.
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.nx = 0;
        assert!(expect_error(spec).contains("invalid grid"));
    }

    /// Constant spacings arrays for a spec's grid (uniform `dx_m`).
    fn constant_spacings(spec: &JobSpec) -> GradedSpacings {
        GradedSpacings {
            dx: vec![spec.dx_m; spec.nx],
            dy: vec![spec.dx_m; spec.ny],
            dz: vec![spec.dx_m; spec.nz],
        }
    }

    #[test]
    fn spacings_round_trip_and_default() {
        // Full round trip with graded spacings attached.
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.spacings = Some(constant_spacings(&spec));
        let json = serde_json::to_string(&spec).unwrap();
        let back: JobSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, back);
        // Pre-FS.0b specs (no spacings key) still deserialize to None.
        let legacy = serde_json::to_string(&cavity_spec(BackendChoice::Cpu))
            .unwrap()
            .replace(",\"spacings\":null", "");
        assert!(!legacy.contains("spacings"), "replace missed the key");
        let back: JobSpec = serde_json::from_str(&legacy).unwrap();
        assert!(back.spacings.is_none());
    }

    /// FS.0b.0 gate at the protocol level: constant spacings arrays equal
    /// to `dx_m` are bit-identical to the scalar path — including the
    /// default dt (0.9× Courant from the minimum spacing).
    #[test]
    fn constant_spacings_match_scalar_path_bit_exactly() {
        let scalar = run_to_done(cavity_spec(BackendChoice::Cpu));
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.spacings = Some(constant_spacings(&spec));
        let graded = run_to_done(spec);
        assert!(
            scalar.dt_s == graded.dt_s,
            "graded default dt diverged: {} vs {}",
            scalar.dt_s,
            graded.dt_s
        );
        assert_eq!(graded.backend, "cpu");
        assert!(scalar.probes[0].iter().any(|v| *v != 0.0));
        for (n, (a, b)) in scalar.probes[0].iter().zip(&graded.probes[0]).enumerate() {
            assert!(a == b, "step {n}: scalar {a:e} != graded {b:e}");
        }
    }

    #[test]
    fn malformed_spacings_error_instead_of_panicking() {
        // Wrong array length.
        let mut spec = cavity_spec(BackendChoice::Cpu);
        let mut g = constant_spacings(&spec);
        g.dy.pop();
        spec.spacings = Some(g);
        assert!(expect_error(spec).contains("spacings.dy"));
        // Non-positive width.
        let mut spec = cavity_spec(BackendChoice::Cpu);
        let mut g = constant_spacings(&spec);
        g.dz[3] = -1.0e-3;
        spec.spacings = Some(g);
        assert!(expect_error(spec).contains("spacings.dz[3]"));
        // Grading inside an absorbing CPML layer (FS.0b.0 scope rule).
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.boundary = BoundarySpec::Cpml {
            npml: 4,
            axes: [true; 3],
            faces: None,
        };
        let mut g = constant_spacings(&spec);
        g.dx[1] = 0.5 * spec.dx_m;
        spec.spacings = Some(g);
        assert!(expect_error(spec).contains("CPML layer"));
        // dt_s above the graded (minimum-spacing) Courant limit, even
        // though it satisfies the scalar dx_m limit.
        let mut spec = cavity_spec(BackendChoice::Cpu);
        let mut g = constant_spacings(&spec);
        for d in &mut g.dx {
            *d *= 0.25; // fine axis tightens the limit
        }
        spec.spacings = Some(g);
        let scalar_limit = FdtdSpec::vacuum(spec.nx, spec.ny, spec.nz, spec.dx_m).courant_limit();
        spec.dt_s = Some(0.9 * scalar_limit);
        assert!(expect_error(spec).contains("Courant"));
        // NTFF on a graded grid is FS.0b.0-unsupported.
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.spacings = Some(constant_spacings(&spec));
        spec.ntff = Some(NtffSpec {
            f_hz: 10.0e9,
            margin_cells: 3,
            k_min: None,
            directions: vec![(0.0, 0.0)],
        });
        assert!(expect_error(spec).contains("graded"));
    }

    /// Graded grids are CPU-only (FS.0b.0): an explicit GPU request fails
    /// with the Unsupported message before any adapter query; auto falls
    /// back to the CPU backend.
    #[cfg(feature = "gpu")]
    #[test]
    fn graded_runs_on_gpu_or_auto_falls_back() {
        // FS.0b.2 (ADR-0214): graded grids run on the GPU backend. Auto
        // must always complete — on "gpu" where an adapter exists, on
        // "cpu" otherwise (the compute-020 uniform-fill scenario).
        let mut spec = cavity_spec(BackendChoice::Auto);
        spec.spacings = Some(constant_spacings(&spec));
        let result = run_to_done(spec);
        assert!(
            result.backend == "gpu" || result.backend == "cpu",
            "unexpected backend {}",
            result.backend
        );
        assert!(result.probes[0].iter().any(|v| *v != 0.0));
    }

    #[test]
    fn cpu_job_streams_progress_and_completes() {
        let handle = submit(cavity_spec(BackendChoice::Cpu));
        let mut progress = 0usize;
        let mut result = None;
        for event in handle.events() {
            match event {
                JobEvent::Progress { step, total } => {
                    assert!(step <= total);
                    progress += 1;
                }
                JobEvent::Done { result: r } => result = Some(r),
                JobEvent::Error { message } => panic!("job failed: {message}"),
            }
        }
        let result = result.expect("no Done event");
        assert!(progress >= 10, "too few progress events: {progress}");
        assert_eq!(result.steps_done, 60);
        assert_eq!(result.probes.len(), 1);
        assert_eq!(result.probes[0].len(), 60);
        assert!(result.probes[0].iter().any(|v| *v != 0.0));
        assert_eq!(result.backend, "cpu");
    }

    #[test]
    fn slice_is_returned_when_requested() {
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.slice = Some(SliceSpec {
            component: "ez".into(),
            k: 6,
        });
        let handle = submit(spec);
        let result = handle
            .events()
            .find_map(|e| match e {
                JobEvent::Done { result } => Some(result),
                JobEvent::Error { message } => panic!("job failed: {message}"),
                _ => None,
            })
            .expect("no Done event");
        let slice = result.slice.expect("no slice returned");
        assert_eq!((slice.ni, slice.nj), (13, 13)); // ez dims: [nx+1, ny+1, nz]
        assert_eq!(slice.data.len(), 13 * 13);
        assert!(slice.data.iter().any(|v| *v != 0.0));
    }

    #[test]
    fn cancellation_stops_early() {
        let mut spec = cavity_spec(BackendChoice::Cpu);
        spec.n_steps = 100_000; // long enough that cancel lands mid-run
        let handle = submit(spec);
        let mut result = None;
        for event in handle.events() {
            match event {
                JobEvent::Progress { step, .. } if step > 0 => handle.cancel(),
                JobEvent::Done { result: r } => result = Some(r),
                JobEvent::Error { message } => panic!("job failed: {message}"),
                _ => {}
            }
        }
        let result = result.expect("no Done event");
        assert!(
            result.steps_done < 100_000,
            "cancel had no effect ({} steps)",
            result.steps_done
        );
    }

    #[test]
    fn auto_backend_always_completes() {
        let handle = submit(cavity_spec(BackendChoice::Auto));
        let done = handle
            .events()
            .find_map(|e| match e {
                JobEvent::Done { result } => Some(result),
                JobEvent::Error { message } => panic!("job failed: {message}"),
                _ => None,
            })
            .expect("no Done event");
        assert!(done.backend == "cpu" || done.backend == "gpu");
        assert!(done.probes[0].iter().any(|v| *v != 0.0));
    }
}
