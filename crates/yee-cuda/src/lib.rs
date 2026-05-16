//! # yee-cuda
//!
//! Thin safe layer over `cudarc` for the rest of the Yee workspace. Everything that
//! touches the CUDA driver, NVRTC, cuBLAS, cuSOLVER, cuSPARSE, cuFFT, or NCCL goes
//! through this crate. The goal is one swap point if `cudarc` ever breaks under us.
//!
//! Phase 0 scope:
//! - Device enumeration ([`Device::list`])
//! - NVRTC kernel compilation helper ([`nvrtc::compile`])
//! - Internal [`backend::Backend`] trait that is the intended swap point for
//!   Phase 1
//! - A "hello world" kernel (lives in `kernels/`)
//!
//! Phase 1 scope: cuSOLVER `Zgetrf`/`Zgetrs` wrappers, cuBLAS GEMM helpers.
//! Phase 2 scope: NCCL boundary-exchange helpers for FDTD domain decomposition.
//!
//! Build with `--features cuda` to compile against a real CUDA toolkit. Without the
//! feature the crate still builds (stub APIs return [`Error::NotEnabled`]) so CI on
//! CPU-only hosts is green.

// Phase 1.5 wires cuSOLVER's dense LU into `cusolver.rs`. The cudarc
// `result::*` + `sys::*` functions for cuSOLVER are `unsafe fn`, so the
// cusolver module locally opts back into unsafe. The rest of the crate still
// denies it.
#![deny(unsafe_code)]
#![warn(missing_docs)]

pub mod backend;
pub mod cusolver;
pub mod nvrtc;

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

/// Static information about a single CUDA device.
///
/// Populated by [`Device::list`] when the crate is built with `--features cuda`
/// and a CUDA driver is reachable.
#[derive(Debug, Clone)]
pub struct Device {
    /// Ordinal of the device as returned by the CUDA driver.
    pub ordinal: usize,
    /// Human-readable device name (e.g. `"NVIDIA RTX 4090"`).
    pub name: String,
    /// CUDA compute capability as `(major, minor)`.
    pub compute_cap: (u8, u8),
    /// Total global memory in bytes.
    pub mem_total_bytes: u64,
}

impl Device {
    /// Enumerate visible CUDA devices.
    ///
    /// # Errors
    ///
    /// - [`Error::NotEnabled`] if the crate was built without `--features cuda`.
    /// - [`Error::Driver`] if the underlying CUDA driver call fails for any
    ///   reason other than "no devices visible" (which yields `Ok(vec![])`).
    ///
    /// # Behaviour matrix
    ///
    /// | feature `cuda` | devices visible | return                     |
    /// |----------------|-----------------|----------------------------|
    /// | off            | n/a             | `Err(Error::NotEnabled)`   |
    /// | on             | 0               | `Ok(vec![])`               |
    /// | on             | n               | `Ok(Vec<Device>)` len `n`  |
    pub fn list() -> Result<Vec<Self>> {
        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "cuda")]
        {
            // Route through the internal backend trait so Phase 1 can swap the
            // binding behind a single seam.
            use crate::backend::{Backend, CudarcBackend};

            let count = CudarcBackend::device_count()?;
            let mut out = Vec::with_capacity(count);
            for i in 0..count {
                let (name, compute_cap, mem_total_bytes) = CudarcBackend::device_props(i)?;
                out.push(Device {
                    ordinal: i,
                    name,
                    compute_cap,
                    mem_total_bytes,
                });
            }
            Ok(out)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn device_struct_constructs_and_reads() {
        let d = Device {
            ordinal: 0,
            name: "TestGPU".to_string(),
            compute_cap: (8, 9),
            mem_total_bytes: 24 * 1024 * 1024 * 1024,
        };
        assert_eq!(d.ordinal, 0);
        assert_eq!(d.name, "TestGPU");
        assert_eq!(d.compute_cap, (8, 9));
        assert_eq!(d.mem_total_bytes, 24 * 1024 * 1024 * 1024);
    }

    #[cfg(not(feature = "cuda"))]
    #[test]
    fn list_without_feature_is_not_enabled() {
        match Device::list() {
            Err(Error::NotEnabled) => {}
            other => panic!("expected Err(NotEnabled), got {other:?}"),
        }
    }

    #[cfg(not(feature = "cuda"))]
    #[test]
    fn nvrtc_compile_without_feature_is_not_enabled() {
        match nvrtc::compile("extern \"C\" __global__ void k(){}", "k") {
            Err(Error::NotEnabled) => {}
            other => panic!("expected Err(NotEnabled), got {other:?}"),
        }
    }

    /// On a CUDA-enabled build, `Device::list()` must succeed and the empty-Vec
    /// case (no GPU on a toolkit-installed host) is treated as success, not an
    /// error. We do not assert on the length because CI hosts vary.
    #[cfg(feature = "cuda")]
    #[test]
    fn list_with_feature_is_ok() {
        let res = Device::list();
        assert!(res.is_ok(), "Device::list returned {res:?}");
    }
}
