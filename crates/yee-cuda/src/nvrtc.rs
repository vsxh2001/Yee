//! NVRTC kernel-compilation helpers.
//!
//! Phase 0 surface area is a single function, [`compile`], that turns a CUDA-C
//! source string into PTX bytes. Without the `cuda` feature it returns
//! [`Error::NotEnabled`] so CPU-only CI hosts can still build the crate.

use crate::{Error, Result};

/// Compile a CUDA-C source string to PTX.
///
/// `name` is the logical program name reported in nvrtc diagnostics (it does
/// **not** need to match any kernel name inside `src`).
///
/// # Errors
///
/// - [`Error::NotEnabled`] if the crate was built without `--features cuda`.
/// - [`Error::Driver`] if nvrtc rejects the source or the runtime cannot be
///   reached (e.g. the CUDA toolkit is absent on the host).
pub fn compile(src: &str, name: &str) -> Result<Vec<u8>> {
    #[cfg(not(feature = "cuda"))]
    {
        let _ = (src, name);
        Err(Error::NotEnabled)
    }
    #[cfg(feature = "cuda")]
    {
        use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};

        let opts = CompileOptions {
            name: Some(name.to_string()),
            ..Default::default()
        };
        // Use `Display` to match the formatter convention in `backend.rs`.
        // cudarc's `CompileError` Display impl currently delegates to Debug,
        // so this loses no diagnostic detail (it just keeps the surface
        // uniform across the crate).
        let ptx = compile_ptx_with_opts(src, opts)
            .map_err(|e| Error::Driver(format!("nvrtc compile failed: {e}")))?;
        // `Ptx::as_bytes` returns `Some` for image-kind Ptx produced by the
        // compile path. We copy into an owned `Vec<u8>` so the caller does not
        // depend on cudarc types in its signature.
        let bytes = ptx
            .as_bytes()
            .ok_or_else(|| Error::Driver("nvrtc returned non-image Ptx".into()))?
            .to_vec();
        Ok(bytes)
    }
}
