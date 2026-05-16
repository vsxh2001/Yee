//! Safe `Session` skeleton wrapping a Gmsh API session.
//!
//! Without the `gmsh` feature this module is a typed placeholder: every method
//! returns [`Error::NotEnabled`]. With the feature on, the bodies will be
//! filled in during Phase 1 by calling into the in-tree bindgen-generated FFI;
//! for Phase 0 they remain `todo!()`. This keeps downstream crates compiling
//! against the type surface on hosts without a Gmsh SDK.

use std::path::Path;

// `Error` is only referenced by the no-feature stubs (which return
// `Error::NotEnabled`). When the `gmsh` feature is on, every method body
// is `todo!()` and never names `Error`, so importing it would be unused.
#[cfg(not(feature = "gmsh"))]
use crate::Error;
use crate::{Result, TriMesh};

/// Owns a single Gmsh API session. Dropping the session releases Gmsh's
/// global state.
///
/// The struct deliberately holds no observable state in the no-feature build:
/// it is a zero-sized handle so that callers can write type-correct code
/// against the API surface even when the crate is built without `gmsh`.
#[derive(Debug)]
pub struct Session {
    // Private field keeps the struct opaque and forces use of `Session::new`.
    _private: (),
}

impl Session {
    /// Open a new Gmsh session.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`] and
    /// performs no work. It never panics.
    pub fn new() -> Result<Self> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 0: real FFI wiring deferred to Phase 1.
            todo!("Session::new: bindgen FFI wiring lands in Phase 1")
        }
    }

    /// Import a STEP file into the current session's OCC kernel.
    pub fn import_step(&mut self, _path: &Path) -> Result<()> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 0: real FFI wiring deferred to Phase 1.
            todo!("Session::import_step: bindgen FFI wiring lands in Phase 1")
        }
    }

    /// Mesh the loaded geometry up to the requested dimension.
    pub fn mesh(&mut self, _dim: i32) -> Result<()> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 0: real FFI wiring deferred to Phase 1.
            todo!("Session::mesh: bindgen FFI wiring lands in Phase 1")
        }
    }

    /// Extract surface triangles from the current mesh.
    pub fn tris(&self) -> Result<TriMesh> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // Phase 0: real FFI wiring deferred to Phase 1.
            todo!("Session::tris: bindgen FFI wiring lands in Phase 1")
        }
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        // Without the feature there is no native session to tear down.
        // With the feature, Phase 1 will call `gmshFinalize` here.
    }
}
