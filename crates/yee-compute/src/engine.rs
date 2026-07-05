//! Backend-dispatching engine front.

use crate::cpu::CpuFdtd;
use crate::error::ComputeError;
use crate::fields::Fields;
use crate::spec::FdtdSpec;

#[cfg(feature = "gpu")]
use crate::gpu::GpuFdtd;

/// An FDTD stepper on a runtime-selected backend.
///
/// Construct explicitly with [`FdtdEngine::new_cpu`] / [`FdtdEngine::new_gpu`],
/// or let [`FdtdEngine::new_auto`] try the GPU and fall back to the CPU.
#[derive(Debug)]
pub enum FdtdEngine {
    /// Rayon-parallel FP64 CPU backend.
    Cpu(CpuFdtd),
    /// wgpu FP32 GPU backend (feature `gpu`).
    #[cfg(feature = "gpu")]
    Gpu(GpuFdtd),
}

impl FdtdEngine {
    /// Build on the CPU backend (always available, FP64).
    pub fn new_cpu(spec: FdtdSpec, fields: Fields) -> Self {
        Self::Cpu(CpuFdtd::new(spec, fields))
    }

    /// Build on the GPU backend (FP32). Fails with
    /// [`ComputeError::NoAdapter`] when no GPU is present, or
    /// [`ComputeError::GpuNotEnabled`] when the crate was built without the
    /// `gpu` feature.
    pub fn new_gpu(spec: FdtdSpec, fields: Fields) -> Result<Self, ComputeError> {
        #[cfg(feature = "gpu")]
        {
            GpuFdtd::new(spec, fields).map(Self::Gpu)
        }
        #[cfg(not(feature = "gpu"))]
        {
            let _ = (spec, fields);
            Err(ComputeError::GpuNotEnabled)
        }
    }

    /// Try the GPU first, silently falling back to the CPU backend.
    pub fn new_auto(spec: FdtdSpec, fields: Fields) -> Self {
        #[cfg(feature = "gpu")]
        {
            match GpuFdtd::new(spec, fields.clone()) {
                Ok(gpu) => Self::Gpu(gpu),
                Err(_) => Self::new_cpu(spec, fields),
            }
        }
        #[cfg(not(feature = "gpu"))]
        {
            Self::new_cpu(spec, fields)
        }
    }

    /// The problem description this engine was built from.
    pub fn spec(&self) -> &FdtdSpec {
        match self {
            Self::Cpu(cpu) => cpu.spec(),
            #[cfg(feature = "gpu")]
            Self::Gpu(gpu) => gpu.spec(),
        }
    }

    /// Short backend identifier (`"cpu"` / `"gpu"`), for diagnostics.
    pub fn backend_name(&self) -> &'static str {
        match self {
            Self::Cpu(_) => "cpu",
            #[cfg(feature = "gpu")]
            Self::Gpu(_) => "gpu",
        }
    }

    /// Advance the state by `n` leapfrog steps.
    pub fn step_n(&mut self, n: usize) -> Result<(), ComputeError> {
        match self {
            Self::Cpu(cpu) => {
                cpu.step_n(n);
                Ok(())
            }
            #[cfg(feature = "gpu")]
            Self::Gpu(gpu) => gpu.step_n(n),
        }
    }

    /// Fetch the current field state (a copy; the GPU backend reads back
    /// through a staging buffer and widens FP32 → FP64).
    pub fn read_fields(&mut self) -> Result<Fields, ComputeError> {
        match self {
            Self::Cpu(cpu) => Ok(cpu.fields().clone()),
            #[cfg(feature = "gpu")]
            Self::Gpu(gpu) => gpu.read_fields(),
        }
    }
}
