//! # yee-mesh
//!
//! Meshing for Yee. The `gmsh` feature enables in-tree bindgen-generated FFI against
//! the upstream `gmshc.h` (Gmsh 4.15+). Without the feature this crate exposes only
//! the safe data structures (triangle mesh, hex grid) so downstream crates compile
//! on hosts without a Gmsh SDK.
//!
//! The pre-existing `rgmsh` crate is unmaintained since 2019 and targets Gmsh 4.4.1;
//! we generate fresh bindings (see `build.rs` in Phase 0).

// Phase 0 forbade `unsafe_code` crate-wide. Phase 1.mesh.0 must call into
// the Gmsh C API via `bindgen`, which is inherently `unsafe`. We narrow the
// forbid to the no-feature build (where the crate is pure data structures)
// and rely on localized `#[allow(unsafe_code)]` inside the FFI submodule.
#![cfg_attr(not(feature = "gmsh"), forbid(unsafe_code))]
#![warn(missing_docs)]

use nalgebra::Vector3;

mod session;

pub use session::Session;

/// Meshing errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Built without the `gmsh` feature.
    #[error("yee-mesh built without the `gmsh` feature; rebuild with --features gmsh")]
    NotEnabled,
    /// Gmsh returned a non-zero status from its C API.
    #[error("gmsh error code {0}")]
    Gmsh(i32),
    /// Geometry validation failed before meshing.
    #[error("invalid geometry: {0}")]
    Invalid(String),
}

/// Mesh-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;

/// A planar triangle mesh — the primary input to `yee-mom`.
#[derive(Debug, Default, Clone)]
pub struct TriMesh {
    /// Vertices in world coordinates.
    pub vertices: Vec<Vector3<f64>>,
    /// Triangle indices into `vertices`, 3 per face.
    pub triangles: Vec<[u32; 3]>,
    /// Per-triangle physical-group tag (used to map ports, materials).
    pub tags: Vec<u32>,
}

impl TriMesh {
    /// Build a `TriMesh` after validating its invariants.
    ///
    /// Currently the only invariant enforced is that there is exactly one
    /// tag per triangle. Returns [`Error::Invalid`] with a descriptive
    /// message when `triangles.len() != tags.len()`.
    ///
    /// Phase 1 will additionally validate that each triangle index is
    /// `< vertices.len()`.
    pub fn new(
        vertices: Vec<Vector3<f64>>,
        triangles: Vec<[u32; 3]>,
        tags: Vec<u32>,
    ) -> Result<Self> {
        if triangles.len() != tags.len() {
            return Err(Error::Invalid(format!(
                "triangles ({}) and tags ({}) must have equal length",
                triangles.len(),
                tags.len()
            )));
        }
        Ok(Self {
            vertices,
            triangles,
            tags,
        })
    }

    /// Number of triangles.
    pub fn n_tris(&self) -> usize {
        self.triangles.len()
    }
}
