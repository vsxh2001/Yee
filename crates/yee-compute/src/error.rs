//! Error type shared by every `yee-compute` backend.

/// Errors surfaced by compute backends.
///
/// The CPU backend is infallible in the E.0 scope; every variant here comes
/// from GPU device discovery, creation, or readback, plus the feature-gate
/// stub used when the crate is built with `--no-default-features`.
#[derive(Debug, thiserror::Error)]
pub enum ComputeError {
    /// No wgpu adapter (GPU) is available on this machine.
    #[error("no compatible GPU adapter found")]
    NoAdapter,

    /// The adapter was found but a device could not be created.
    #[error("GPU device request failed: {0}")]
    Device(String),

    /// Mapping a staging buffer back to the host failed.
    #[error("GPU readback failed: {0}")]
    Readback(String),

    /// The crate was built without the `gpu` feature.
    #[error("yee-compute was built without the `gpu` feature")]
    GpuNotEnabled,
}
