//! Transport-agnostic simulation job API over the `yee-compute` engine
//! (phase S.0, ADR-0179).
//!
//! One serde protocol, many transports: the Tauri studio calls
//! [`submit`] in-process; `yee-server` (S.1) will forward the same
//! [`JobSpec`] / [`JobEvent`] types over WebSocket; the CLI can print the
//! event stream. The S.0 walking skeleton covers driven vacuum FDTD jobs
//! (any boundary, soft sources + resistive ports, probes) with progress
//! streaming and cooperative cancellation; materials/dispersion plumb
//! through in a later slice.
//!
//! ```
//! use yee_engine::{BackendChoice, BoundarySpec, JobEvent, JobSpec, SourceSpec, ProbeSpec};
//!
//! let spec = JobSpec {
//!     nx: 12, ny: 12, nz: 12, dx_m: 1e-3, n_steps: 40,
//!     boundary: BoundarySpec::Pec,
//!     sources: vec![SourceSpec::GaussianEz { cell: (6, 6, 6), t0_steps: 8.0, sigma_steps: 3.0 }],
//!     ports: vec![],
//!     probes: vec![ProbeSpec { component: "ez".into(), cell: (8, 6, 6) }],
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

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::thread;

use serde::{Deserialize, Serialize};

use yee_compute::{
    Boundary, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials, Probe,
    ResistivePort, SoftSource, Waveform,
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
    /// Roden–Gedney CPML with `npml` layers on every face.
    Cpml {
        /// PML thickness in cells.
        npml: usize,
    },
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
    /// Probes (recorded every step).
    pub probes: Vec<ProbeSpec>,
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
}

/// Validate a spec and spawn it on a worker thread.
pub fn submit(spec: JobSpec) -> JobHandle {
    let (tx, rx) = channel();
    let cancel = Arc::new(AtomicBool::new(false));
    let flag = Arc::clone(&cancel);
    thread::spawn(move || run_job(spec, &tx, &flag));
    JobHandle { events: rx, cancel }
}

fn build_drive(spec: &JobSpec, dt: f64) -> Result<Drive, String> {
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

/// Chunked execution with progress events (~2 % granularity, min 1 step).
fn run_job(spec: JobSpec, tx: &Sender<JobEvent>, cancel: &AtomicBool) {
    let fdtd_spec = FdtdSpec::vacuum(spec.nx, spec.ny, spec.nz, spec.dx_m);
    let boundary = match spec.boundary {
        BoundarySpec::None => Boundary::None,
        BoundarySpec::Pec => Boundary::PecBox,
        BoundarySpec::Cpml { npml } => Boundary::Cpml(CpmlConfig::for_spec(&fdtd_spec, npml)),
    };
    let drive = match build_drive(&spec, fdtd_spec.dt) {
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
    if matches!(spec.backend, BackendChoice::Gpu | BackendChoice::Auto) {
        match yee_compute::GpuFdtd::with_drive(
            fdtd_spec,
            Fields::zero(&fdtd_spec),
            Materials::default(),
            boundary.clone(),
            drive.clone(),
            total,
        ) {
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
                match engine.read_probes() {
                    Ok(probes) => {
                        let _ = tx.send(JobEvent::Done {
                            result: JobResult {
                                backend: "gpu".into(),
                                dt_s,
                                probes,
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
            Err(yee_compute::ComputeError::NoAdapter)
                if matches!(spec.backend, BackendChoice::Auto) => {} // fall through to CPU
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
        Materials::default(),
        boundary,
        drive,
    );
    let mut done = 0usize;
    while done < total && !cancel.load(Ordering::Relaxed) {
        let n = chunk.min(total - done);
        engine.step_n(n);
        done += n;
        let _ = tx.send(JobEvent::Progress { step: done, total });
    }
    let _ = tx.send(JobEvent::Done {
        result: JobResult {
            backend: "cpu".into(),
            dt_s,
            probes: engine.probe_series().to_vec(),
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
            probes: vec![ProbeSpec {
                component: "ez".into(),
                cell: (8, 6, 6),
            }],
            backend,
        }
    }

    #[test]
    fn spec_round_trips_through_json() {
        let spec = cavity_spec(BackendChoice::Auto);
        let json = serde_json::to_string(&spec).unwrap();
        let back: JobSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, back);
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
