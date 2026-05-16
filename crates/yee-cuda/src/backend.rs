//! Internal backend trait: the single swap point if cudarc breaks under us.
//!
//! Phase 0 keeps the rest of the crate calling cudarc directly through small
//! helper modules ([`crate::nvrtc`] and [`Device`](crate::Device)). This trait
//! exists so Phase 1 can route those call sites through a single seam without
//! churning the public API.
//
// TODO(phase-1): route `Device::list` and `nvrtc::compile` through this trait
// so swapping the underlying binding only touches `CudarcBackend`.

use crate::Result;

/// Operations the CUDA layer needs from its underlying binding.
///
/// All methods are associated functions because Phase 0 has no per-instance
/// state. Phase 1 may convert these to `&self` methods once a context handle
/// becomes part of the backend.
pub trait Backend {
    /// Number of CUDA devices visible to the process.
    fn device_count() -> Result<usize>;

    /// Static properties of device `i`: name, compute capability `(major, minor)`,
    /// and total global memory in bytes.
    fn device_props(i: usize) -> Result<(String, (u8, u8), u64)>;

    /// Compile a CUDA-C source string to PTX bytes.
    fn nvrtc_compile(src: &str, name: &str) -> Result<Vec<u8>>;
}

/// Cudarc-backed [`Backend`] implementation.
///
/// Available only when the crate is built with `--features cuda`.
#[cfg(feature = "cuda")]
pub struct CudarcBackend;

#[cfg(feature = "cuda")]
impl Backend for CudarcBackend {
    fn device_count() -> Result<usize> {
        use crate::Error;
        use cudarc::driver::CudaContext;
        let n = CudaContext::device_count().map_err(|e| Error::Driver(format!("{e}")))?;
        Ok(n.max(0) as usize)
    }

    fn device_props(i: usize) -> Result<(String, (u8, u8), u64)> {
        use crate::Error;
        use cudarc::driver::CudaContext;
        let ctx = CudaContext::new(i).map_err(|e| Error::Driver(format!("{e}")))?;
        let name = ctx.name().map_err(|e| Error::Driver(format!("{e}")))?;
        let (major, minor) = ctx
            .compute_capability()
            .map_err(|e| Error::Driver(format!("{e}")))?;
        let mem = ctx.total_mem().map_err(|e| Error::Driver(format!("{e}")))?;
        Ok((name, (major as u8, minor as u8), mem as u64))
    }

    fn nvrtc_compile(src: &str, name: &str) -> Result<Vec<u8>> {
        crate::nvrtc::compile(src, name)
    }
}
