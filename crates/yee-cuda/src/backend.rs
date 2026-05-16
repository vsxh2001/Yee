//! Internal backend trait: the single swap point if cudarc breaks under us.
//!
//! Phase 0 routes [`Device::list`](crate::Device::list) and [`nvrtc::compile`](crate::nvrtc::compile)
//! through this trait so swapping the underlying binding only touches
//! `CudarcBackend` (defined below under `#[cfg(feature = "cuda")]`).
//
// Phase-1 deferral: the trait uses associated functions, which means it is
// NOT object-safe — callers cannot hold a `Box<dyn Backend>`. That is fine for
// Phase 0 because we statically pick `CudarcBackend` at compile time. When
// Phase 1 needs runtime backend selection (e.g. a mock backend for tests, or
// an alternative binding), convert these to `&self` methods that take a
// shared context handle.

/// Maximum number of CUDA devices we are willing to enumerate.
///
/// The CUDA driver theoretically returns an `i32`, so an out-of-band value
/// (negative or absurdly large) almost certainly indicates a corrupted runtime
/// or memory smash. 256 is well above any realistic multi-GPU host (DGX-class
/// boxes top out at 8–16) and well below `i32::MAX`, so it makes a reasonable
/// guard rail. Gated to the cuda-feature build because that is the only path
/// that consults it.
#[cfg(feature = "cuda")]
const MAX_DEVICES: i32 = 256;

use crate::Result;

/// Operations the CUDA layer needs from its underlying binding.
///
/// All methods are associated functions because Phase 0 has no per-instance
/// state. Phase 1 may convert these to `&self` methods once a context handle
/// becomes part of the backend (see the module-level note about object
/// safety).
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
        let raw = CudaContext::device_count().map_err(|e| Error::Driver(format!("{e}")))?;
        if !(0..=MAX_DEVICES).contains(&raw) {
            return Err(Error::Driver(format!(
                "device_count returned implausible value: {raw} (expected 0..={MAX_DEVICES})"
            )));
        }
        // `raw` is in `0..=MAX_DEVICES`, so the cast is lossless on any
        // platform where `usize` is at least 16 bits.
        Ok(raw as usize)
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

        // Guarded conversions: surface a `Driver` error rather than silently
        // truncating if the underlying driver ever returns a value that does
        // not fit. All realistic compute caps fit in u8 (Hopper is 9.0) and
        // 64-bit hosts widen `usize` to `u64` losslessly, but defending here
        // keeps the pattern uniform across the backend surface.
        let major = u8::try_from(major).map_err(|_| {
            Error::Driver(format!("compute capability major {major} overflows u8"))
        })?;
        let minor = u8::try_from(minor).map_err(|_| {
            Error::Driver(format!("compute capability minor {minor} overflows u8"))
        })?;
        let mem = u64::try_from(mem)
            .map_err(|_| Error::Driver(format!("total_mem {mem} overflows u64")))?;

        Ok((name, (major, minor), mem))
    }

    fn nvrtc_compile(src: &str, name: &str) -> Result<Vec<u8>> {
        crate::nvrtc::compile(src, name)
    }
}
