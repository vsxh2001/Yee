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
//! The E.0 walking-skeleton scope is a uniform lossless vacuum inside a PEC
//! box: no CPML, per-cell materials, sources, or NTFF yet — those are phases
//! E.1/E.2 in `ENGINE-STUDIO-ROADMAP.md`. The spec lives at
//! `docs/superpowers/specs/2026-07-05-gpu-engine-web-studio-design.md`.
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

mod cpu;
mod engine;
mod error;
mod fields;
mod spec;

#[cfg(feature = "gpu")]
mod gpu;

pub use cpu::CpuFdtd;
pub use engine::FdtdEngine;
pub use error::ComputeError;
pub use fields::Fields;
pub use spec::FdtdSpec;

#[cfg(feature = "gpu")]
pub use gpu::GpuFdtd;
