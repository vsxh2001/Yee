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

    /// Import a STEP file into the current session's OCC kernel and
    /// synchronize the OCC representation with the Gmsh model.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`].
    pub fn import_step(&mut self, _path: &Path) -> Result<()> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // The C API wants a NUL-terminated path. `Path::to_string_lossy`
            // can elide invalid UTF-8 — acceptable here because Gmsh itself
            // expects ASCII/UTF-8 paths.
            let path_str = _path.to_string_lossy().into_owned();
            let c_path = std::ffi::CString::new(path_str.as_bytes())
                .map_err(|_| Error::Invalid("STEP path contains interior NUL byte".to_string()))?;
            // Empty `format` lets Gmsh infer from the file extension.
            let c_fmt = std::ffi::CString::new("").expect("empty CString cannot fail");

            let mut out_dim_tags: *mut ffi::c_int = std::ptr::null_mut();
            let mut out_dim_tags_n: usize = 0;
            let mut ierr: ffi::c_int = 0;
            // SAFETY: all out-pointers point to local stack slots; the input
            // strings outlive the call (they're owned `CString`s held by this
            // function). Gmsh allocates `out_dim_tags` and we free it via
            // `gmshFree` below.
            unsafe {
                ffi::gmshModelOccImportShapes(
                    c_path.as_ptr(),
                    &mut out_dim_tags,
                    &mut out_dim_tags_n,
                    c_fmt.as_ptr(),
                    &mut ierr,
                );
            }
            // Free even on failure if Gmsh allocated.
            let _out_dim_tags_guard = GmshBuffer::from_raw(out_dim_tags as *mut ffi::c_void);
            if ierr != 0 {
                return Err(Error::Gmsh(ierr as i32));
            }
            tracing::debug!(
                imported_dim_tags = out_dim_tags_n / 2,
                "gmshModelOccImportShapes succeeded"
            );

            // Sync OCC → Gmsh model so subsequent mesh ops see the geometry.
            self.synchronize()?;
            Ok(())
        }
    }

    /// Synchronize the OCC kernel representation with the Gmsh model.
    /// Required after OCC geometry edits (import, primitives, booleans)
    /// before meshing.
    ///
    /// Without the `gmsh` feature this returns [`Error::NotEnabled`].
    pub fn synchronize(&mut self) -> Result<()> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            let mut ierr: ffi::c_int = 0;
            // SAFETY: `ierr` is a valid mutable pointer.
            unsafe {
                ffi::gmshModelOccSynchronize(&mut ierr);
            }
            if ierr != 0 {
                return Err(Error::Gmsh(ierr as i32));
            }
            Ok(())
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
            let mut ierr: ffi::c_int = 0;
            // SAFETY: scalar arg by value; `ierr` is a valid mutable pointer.
            unsafe {
                ffi::gmshModelMeshGenerate(_dim as ffi::c_int, &mut ierr);
            }
            if ierr != 0 {
                return Err(Error::Gmsh(ierr as i32));
            }
            Ok(())
        }
    }

    /// Extract surface triangles from the current mesh.
    ///
    /// Builds a vertex array from all mesh nodes (re-indexed densely from
    /// Gmsh's sparse `nodeTag` space) and a triangle index array from all
    /// 3-node triangle elements (Gmsh element type 2). Per-triangle tags are
    /// zeroed; physical-group propagation lands in Phase 1.mesh.1.
    pub fn tris(&self) -> Result<TriMesh> {
        #[cfg(not(feature = "gmsh"))]
        {
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "gmsh")]
        {
            // --- Nodes ---
            let mut node_tags_ptr: *mut usize = std::ptr::null_mut();
            let mut node_tags_n: usize = 0;
            let mut coord_ptr: *mut f64 = std::ptr::null_mut();
            let mut coord_n: usize = 0;
            let mut param_ptr: *mut f64 = std::ptr::null_mut();
            let mut param_n: usize = 0;
            let mut ierr: ffi::c_int = 0;
            // SAFETY: all out-pointers point to local stack slots; Gmsh
            // allocates the buffers and we free each below.
            unsafe {
                ffi::gmshModelMeshGetNodes(
                    &mut node_tags_ptr,
                    &mut node_tags_n,
                    &mut coord_ptr,
                    &mut coord_n,
                    &mut param_ptr,
                    &mut param_n,
                    /* dim = */ -1,
                    /* tag = */ -1,
                    /* includeBoundary = */ 0,
                    /* returnParametricCoord = */ 0,
                    &mut ierr,
                );
            }
            // RAII-style free: guards run on every exit path.
            let _node_tags_guard = GmshBuffer::from_raw(node_tags_ptr as *mut ffi::c_void);
            let _coord_guard = GmshBuffer::from_raw(coord_ptr as *mut ffi::c_void);
            let _param_guard = GmshBuffer::from_raw(param_ptr as *mut ffi::c_void);
            if ierr != 0 {
                return Err(Error::Gmsh(ierr as i32));
            }
            if coord_n != node_tags_n * 3 {
                return Err(Error::Invalid(format!(
                    "gmshModelMeshGetNodes: coord_n ({coord_n}) != 3 * node_tags_n ({node_tags_n})"
                )));
            }

            // SAFETY: `node_tags_ptr` / `coord_ptr` are valid for
            // `node_tags_n` / `coord_n` elements respectively per the call.
            let node_tags_slice: &[usize] = if node_tags_n == 0 {
                &[]
            } else {
                unsafe { std::slice::from_raw_parts(node_tags_ptr, node_tags_n) }
            };
            let coord_slice: &[f64] = if coord_n == 0 {
                &[]
            } else {
                unsafe { std::slice::from_raw_parts(coord_ptr, coord_n) }
            };

            // Build dense vertex array + sparse `nodeTag` → dense index map.
            // Gmsh tags are 1-indexed but not necessarily contiguous.
            let mut vertices: Vec<nalgebra::Vector3<f64>> = Vec::with_capacity(node_tags_n);
            let mut tag_to_index: std::collections::HashMap<usize, u32> =
                std::collections::HashMap::with_capacity(node_tags_n);
            for (i, &tag) in node_tags_slice.iter().enumerate() {
                let off = i * 3;
                vertices.push(nalgebra::Vector3::new(
                    coord_slice[off],
                    coord_slice[off + 1],
                    coord_slice[off + 2],
                ));
                tag_to_index.insert(tag, i as u32);
            }

            // --- Triangle elements (Gmsh element type 2 = 3-node triangle) ---
            const GMSH_ELEMENT_TYPE_TRIANGLE: ffi::c_int = 2;
            let mut elem_tags_ptr: *mut usize = std::ptr::null_mut();
            let mut elem_tags_n: usize = 0;
            let mut elem_node_tags_ptr: *mut usize = std::ptr::null_mut();
            let mut elem_node_tags_n: usize = 0;
            let mut ierr2: ffi::c_int = 0;
            // SAFETY: out-pointers point to local slots; Gmsh allocates each
            // buffer and we free via the guards below.
            unsafe {
                ffi::gmshModelMeshGetElementsByType(
                    GMSH_ELEMENT_TYPE_TRIANGLE,
                    &mut elem_tags_ptr,
                    &mut elem_tags_n,
                    &mut elem_node_tags_ptr,
                    &mut elem_node_tags_n,
                    /* tag = */ -1,
                    /* task = */ 0,
                    /* numTasks = */ 1,
                    &mut ierr2,
                );
            }
            let _elem_tags_guard = GmshBuffer::from_raw(elem_tags_ptr as *mut ffi::c_void);
            let _elem_node_tags_guard =
                GmshBuffer::from_raw(elem_node_tags_ptr as *mut ffi::c_void);
            if ierr2 != 0 {
                return Err(Error::Gmsh(ierr2 as i32));
            }
            if elem_node_tags_n != elem_tags_n * 3 {
                return Err(Error::Invalid(format!(
                    "gmshModelMeshGetElementsByType(2): elem_node_tags_n ({elem_node_tags_n}) \
                     != 3 * elem_tags_n ({elem_tags_n})"
                )));
            }

            // SAFETY: `elem_node_tags_ptr` is valid for `elem_node_tags_n`.
            let elem_node_tags_slice: &[usize] = if elem_node_tags_n == 0 {
                &[]
            } else {
                unsafe { std::slice::from_raw_parts(elem_node_tags_ptr, elem_node_tags_n) }
            };

            let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(elem_tags_n);
            for tri in elem_node_tags_slice.chunks_exact(3) {
                let mut idx = [0_u32; 3];
                for (slot, &t) in idx.iter_mut().zip(tri.iter()) {
                    *slot = *tag_to_index.get(&t).ok_or_else(|| {
                        Error::Invalid(format!(
                            "triangle references unknown node tag {t} not returned by \
                             gmshModelMeshGetNodes"
                        ))
                    })?;
                }
                triangles.push(idx);
            }
            let tags = vec![0_u32; triangles.len()];

            TriMesh::new(vertices, triangles, tags)
        }
    }
}

/// RAII guard that frees a Gmsh-allocated buffer via `gmshFree` on drop.
/// Holding null is a no-op so call sites can construct one unconditionally.
#[cfg(feature = "gmsh")]
struct GmshBuffer {
    ptr: *mut ffi::c_void,
}

#[cfg(feature = "gmsh")]
impl GmshBuffer {
    fn from_raw(ptr: *mut ffi::c_void) -> Self {
        Self { ptr }
    }
}

#[cfg(feature = "gmsh")]
impl Drop for GmshBuffer {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            // SAFETY: pointer was allocated by Gmsh and not yet freed; the
            // documented release function is `gmshFree`.
            unsafe { ffi::gmshFree(self.ptr) };
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
