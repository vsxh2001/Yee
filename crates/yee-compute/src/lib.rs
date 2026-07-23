//! Portable GPU/CPU execution layer for Yee grid solvers (ADR-0175, phase E.0).
//!
//! `yee-compute` owns *execution* — flat field buffers, parallel kernels,
//! device selection, readback — while `yee-fdtd` remains the scalar FP64
//! *reference implementation* every backend here is gated against:
//!
//! - [`CpuFdtd`]: rayon-parallel FP64 kernels, **bit-exact** against
//!   `yee_fdtd::update::{update_h, update_e}` for the uniform lossless arm
//!   (gate `compute-001`).
//! - [`GpuFdtd`] (feature `gpu`, default-on): wgpu compute (Vulkan / Metal /
//!   DX12) running FP32 WGSL kernels, gated against the CPU backend within
//!   FP32 accumulation tolerance (gate `compute-002`; self-skips when no
//!   adapter is present, runs for real on the GPU nightly).
//!
//! E.1 (ADR-0176) added Roden–Gedney CPML absorbing boundaries, per-cell
//! ε_r / μ_r / σ maps (lossy CA/CB update), interior PEC masks, the legacy
//! PEC box, and a Gaussian soft source on the CPU backend — configured via
//! [`Materials`] + [`Boundary`] through the `with_config` constructors, and
//! gated by `compute-003` (bit-exact), `compute-004` (CPML ≥ 30 dB), and
//! `compute-005` (GPU parity). Sources/ports as engine primitives are E.2;
//! dispersive ADE + NTFF are E.5 — see `ENGINE-STUDIO-ROADMAP.md`. Specs:
//! `docs/superpowers/specs/2026-07-05-gpu-engine-web-studio-design.md` and
//! `docs/superpowers/specs/2026-07-06-e1-cpml-materials-design.md`.
//!
//! FS.0b.0 (ADR-0208) adds per-axis nonuniform primal spacings to the CPU
//! backend ([`GradedSpacings`] + [`CpuFdtd::set_spacings`]): H updates
//! divide by the primal cell width, E updates by the dual spacing, gated
//! bit-exact on uniform arrays (`compute-018`) and by a measured graded
//! interface-reflection floor (`compute-019`). FS.0b.2 (ADR-0214) brings
//! the same capability to the GPU backend (`GpuFdtd::set_spacings`): the
//! WGSL kernels multiply by per-cell **inverse** primal/dual spacings from
//! a packed storage buffer whose uniform fill is bit-equal to the retired
//! scalar `inv_dx/dy/dz` uniforms, gated bit-for-bit on uniform arrays
//! against the GPU's own scalar path (`compute-020`) and cross-backend on
//! the compute-019 taper scenario (`compute-021`). NTFF-DFT + graded and
//! z-taper-straddling aperture ports are rejected with
//! [`ComputeError::Unsupported`].
//!
//! # Example
//!
//! ```
//! use yee_compute::{Fields, FdtdEngine, FdtdSpec};
//!
//! let spec = FdtdSpec::vacuum(16, 16, 16, 1e-3);
//! let init = Fields::with_gaussian_ez(&spec, (8, 8, 8), 2.0);
//! let mut engine = FdtdEngine::new_cpu(spec, init);
//! engine.step_n(10).unwrap();
//! let fields = engine.read_fields().unwrap();
//! assert!(fields.hx.iter().any(|v| *v != 0.0));
//! ```

mod cpml;
mod cpu;
mod dispersive;
mod drive;
mod engine;
mod error;
mod fields;
mod materials;
mod spec;

#[cfg(feature = "gpu")]
mod gpu;

pub use cpu::CpuFdtd;
pub use dispersive::{DispersiveMap, DispersiveMaterial};
pub use drive::{
    AperturePort, Drive, EComponent, HComponent, HProbe, Probe, ResistivePort, SoftSource, Waveform,
};
pub use engine::FdtdEngine;
pub use error::ComputeError;
pub use fields::Fields;
pub use materials::{Boundary, CpmlConfig, Materials};
pub use spec::{FdtdSpec, GradedSpacings};

#[cfg(feature = "gpu")]
pub use gpu::GpuFdtd;
