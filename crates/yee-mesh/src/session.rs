//! Safe `Session` skeleton wrapping a Gmsh API session.
//!
//! Without the `gmsh` feature this module is a typed placeholder: every method
//! returns [`Error::NotEnabled`]. With the feature on, the bodies call into the
//! in-tree bindgen-generated FFI (see `build.rs`).
//!
//! The crate-level `#![forbid(unsafe_code)]` is feature-gated off for the
//! `gmsh` build because the FFI calls in this file are inherently `unsafe`.
//! All `unsafe` is localized to the [`ffi`] submodule and to the bodies of
//! `Session` methods immediately around the raw calls; data shaping and error
//! mapping happens in safe code.

use std::path::Path;

use crate::{Error, Result, TriMesh};

/// Raw FFI bindings to `gmshc.h`, generated at build time by `bindgen` and
/// included from `$OUT_DIR/bindings.rs`. When `$GMSH_SDK_ROOT` is unset the
/// build script writes an empty stub there — feature-gated builds on hosts
/// without the SDK will compile this module but fail to link any FFI call.
#[cfg(feature = "gmsh")]
#[allow(unsafe_code, non_upper_case_globals, non_camel_case_types, non_snake_case)]
#[allow(dead_code)]
pub(crate) mod ffi {
    pub use std::os::raw::{c_int, c_void};

    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

/// Owns a single Gmsh API session. Dropping the session releases Gmsh's
/// global state via `gmshFinalize`.
///
/// In the no-feature build the struct is an opaque zero-sized handle so that
/// callers can write type-correct code against the API surface even when the
/// crate is built without `gmsh`.
#[derive(Debug)]
pub struct Session {
    // Marker used to track whether `gmshFinalize` must be called on drop.
    // In the feature-on build this is set to `true` after a successful
    // `gmshInitialize`. In the no-feature build it is unused but kept so the
    // struct shape is identical across feature configurations.
    initialized: bool,
}

impl Session {
    /// Open a new Gmsh session.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`] and
    /// performs no work. With the feature on, calls `gmshInitialize`.
    pub fn new() -> Result<Self> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            let mut ierr: ffi::c_int = 0;
            // SAFETY: `gmshInitialize` accepts a null `argv` when `argc == 0`
            // per the upstream C API. `ierr` is a valid mutable pointer to a
            // stack-allocated int.
            unsafe {
                ffi::gmshInitialize(
                    0,
                    std::ptr::null_mut(),
                    /* readConfigFiles = */ 0,
                    /* run = */ 0,
                    &mut ierr,
                );
            }
            if ierr != 0 {
                return Err(Error::Gmsh(ierr as i32));
            }
            Ok(Self { initialized: true })
        }
    }

    /// Import a STEP file into the current session's OCC kernel.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`] and
    /// performs no work. With the `gmsh` feature, this panics in Phase 0
    /// (FFI wiring deferred to Phase 1).
    pub fn import_step(&mut self, _path: &Path) -> Result<()> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 1.mesh.0 follow-up commit fills this in.
            todo!("Session::import_step: bindgen FFI wiring lands in follow-up commit")
        }
    }

    /// Mesh the loaded geometry up to the requested dimension.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`] and
    /// performs no work. With the `gmsh` feature, this panics in Phase 0
    /// (FFI wiring deferred to Phase 1).
    pub fn mesh(&mut self, _dim: i32) -> Result<()> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 1.mesh.0 follow-up commit fills this in.
            todo!("Session::mesh: bindgen FFI wiring lands in follow-up commit")
        }
    }

    /// Extract surface triangles from the current mesh.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`] and
    /// performs no work. With the `gmsh` feature, this panics in Phase 0
    /// (FFI wiring deferred to Phase 1).
    pub fn tris(&self) -> Result<TriMesh> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 1.mesh.0 follow-up commit fills this in.
            todo!("Session::tris: bindgen FFI wiring lands in follow-up commit")
        }
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        #[cfg(not(feature = "gmsh"))]
        {
            let _ = self.initialized;
        }
        #[cfg(feature = "gmsh")]
        {
            if !self.initialized {
                return;
            }
            let mut ierr: ffi::c_int = 0;
            // SAFETY: `ierr` is a valid mutable pointer; `gmshFinalize` is
            // safe to call after a successful `gmshInitialize`.
            unsafe {
                ffi::gmshFinalize(&mut ierr);
            }
            if ierr != 0 {
                // Drop must not panic; log instead.
                tracing::error!(ierr = ierr as i32, "gmshFinalize returned non-zero");
            }
        }
    }
}
