//! # yee-cuda
//!
//! Thin safe layer over [`cudarc`] for the rest of the Yee workspace. Everything that
//! touches the CUDA driver, NVRTC, cuBLAS, cuSOLVER, cuSPARSE, cuFFT, or NCCL goes
//! through this crate. The goal is one swap point if `cudarc` ever breaks under us.
//!
//! Phase 0 scope:
//! - Device enumeration (`Device::list`)
//! - Context/stream RAII handles
//! - NVRTC kernel compilation helpers
//! - A "hello world" stencil kernel (lives in `kernels/`)
//!
//! Phase 1 scope: cuSOLVER `Zgetrf`/`Zgetrs` wrappers, cuBLAS GEMM helpers.
//! Phase 2 scope: NCCL boundary-exchange helpers for FDTD domain decomposition.
//!
//! Build with `--features cuda` to compile against a real CUDA toolkit. Without the
//! feature the crate still builds (stub APIs return [`Error::NotEnabled`]) so CI on
//! CPU-only hosts is green.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// CUDA-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Built without the `cuda` feature.
    #[error("yee-cuda built without the `cuda` feature; rebuild with --features cuda")]
    NotEnabled,
    /// Underlying cudarc / driver error.
    #[error("cuda error: {0}")]
    Driver(String),
}

/// CUDA-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;

/// Stub device handle. Phase 0 returns [`Error::NotEnabled`] without the feature.
#[derive(Debug, Clone)]
pub struct Device {
    /// Ordinal of the device as returned by the CUDA driver.
    pub ordinal: usize,
    /// Human-readable name.
    pub name: String,
}

impl Device {
    /// Enumerate visible CUDA devices. Phase 0 stub.
    pub fn list() -> Result<Vec<Self>> {
        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "cuda")]
        {
            // TODO(phase-0): use `cudarc::driver::CudaContext::count` + `name`.
            Err(Error::Driver(
                "device enumeration not yet implemented".into(),
            ))
        }
    }
}
