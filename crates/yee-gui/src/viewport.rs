// Allow temporarily-unused items: the wgpu paint callback and tab wiring that
// consume `Vertex`, `TriMeshCpu`, and `thin_cylinder` land in follow-up
// commits inside the same phase 1.gui.1 series.
#![allow(dead_code)]

//! CPU-side mesh data + cylinder fixture for the Phase 1.gui.1 3D viewport.
//!
//! This first slice of the viewport module sets up the data plumbing the wgpu
//! paint callback will consume in the next commit:
//!
//! - [`Vertex`] — `Pod`/`Zeroable` interleaved position + normal layout matching
//!   the WGSL shader. Includes [`Vertex::layout`] returning the
//!   `wgpu::VertexBufferLayout` so the render pipeline can be built from one
//!   source of truth.
//! - [`TriMeshCpu`] — vertex/index buffers ready for GPU upload, plus the
//!   axis-aligned bounding-box max dimension used to seed the orbit camera.
//! - [`thin_cylinder`] — hand-coded triangulation of a cylinder's lateral
//!   surface, mirroring the `yee-mom` test fixture so this crate does not take
//!   a dev-fixtures dependency on the mom crate.
//!
//! The wgpu render resources, paint callback, camera state, and shader are
//! introduced in follow-up commits.

use bytemuck::{Pod, Zeroable};
use egui_wgpu::wgpu;
use glam::Vec3;

// ----------------------------------------------------------------------------
// Vertex layout
// ----------------------------------------------------------------------------

/// Vertex layout for the triangle-mesh pipeline: interleaved position + normal.
///
/// Per-vertex normals are computed by averaging adjacent face normals at mesh
/// construction time, so a single index buffer supports both flat- and
/// smooth-shaded fragments. Marked `Pod` + `Zeroable` so the whole vertex
/// slice can be cast straight into a `wgpu::Buffer` via `bytemuck`.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Pod, Zeroable)]
pub struct Vertex {
    /// Object-space position in metres.
    pub position: [f32; 3],
    /// Object-space normal (unit length, or `+Z` for degenerate fans).
    pub normal: [f32; 3],
}

impl Vertex {
    /// `wgpu` vertex-buffer layout matching this struct.
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 12,
                    shader_location: 1,
                },
            ],
        }
    }
}

// ----------------------------------------------------------------------------
// CPU mesh
// ----------------------------------------------------------------------------

/// CPU-side triangle mesh ready for GPU upload.
#[derive(Debug, Clone)]
pub struct TriMeshCpu {
    /// Interleaved position/normal vertices.
    pub vertices: Vec<Vertex>,
    /// Triangle indices into [`Self::vertices`] (CCW front faces).
    pub indices: Vec<u32>,
    /// Largest axis-aligned bounding-box extent, used to seed the camera
    /// distance so the mesh fits the default view.
    pub bbox_max_dim: f32,
}

impl TriMeshCpu {
    /// Construct a mesh from raw triangle data, computing per-vertex normals
    /// as the area-weighted average of incident triangle normals.
    fn from_triangles(positions: &[Vec3], triangles: &[[u32; 3]]) -> Self {
        let mut normals = vec![Vec3::ZERO; positions.len()];
        for tri in triangles {
            let p0 = positions[tri[0] as usize];
            let p1 = positions[tri[1] as usize];
            let p2 = positions[tri[2] as usize];
            // Area-weighted normal: |e1 × e2| ∝ 2·area, so just accumulate the
            // raw cross product and normalise at the end.
            let n = (p1 - p0).cross(p2 - p0);
            normals[tri[0] as usize] += n;
            normals[tri[1] as usize] += n;
            normals[tri[2] as usize] += n;
        }
        let vertices: Vec<Vertex> = positions
            .iter()
            .zip(normals.iter())
            .map(|(p, n)| {
                let nn = if n.length_squared() > 0.0 {
                    n.normalize()
                } else {
                    Vec3::Z
                };
                Vertex {
                    position: [p.x, p.y, p.z],
                    normal: [nn.x, nn.y, nn.z],
                }
            })
            .collect();

        let indices: Vec<u32> = triangles.iter().flat_map(|t| t.iter().copied()).collect();

        let (mut min, mut max) = (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY));
        for p in positions {
            min = min.min(*p);
            max = max.max(*p);
        }
        let extent = max - min;
        let bbox_max_dim = extent.x.max(extent.y).max(extent.z);

        Self {
            vertices,
            indices,
            bbox_max_dim,
        }
    }
}

/// Triangulate the lateral surface of a cylinder (no end caps).
///
/// The cylinder's axis is along `z`, centred at the origin. `length_m` is the
/// total length; `radius_m` is the cylinder radius. `n_axial` is the number of
/// axial segments (rings of triangles between adjacent z-cuts); `n_around` is
/// the number of segments around the circumference.
///
/// Two triangles are produced per `(axial × around)` cell, so the total
/// triangle count is `2 · n_axial · n_around`. The vertex count is
/// `(n_axial + 1) · n_around`.
///
/// This is a direct mirror of the `yee-mom` test fixture's triangulation
/// strategy, intentionally duplicated rather than imported so the GUI crate
/// does not take a dev-fixtures dependency.
///
/// # Panics
///
/// Panics if `n_axial < 2` or `n_around < 3`.
pub fn thin_cylinder(length_m: f32, radius_m: f32, n_axial: usize, n_around: usize) -> TriMeshCpu {
    assert!(n_axial >= 2, "n_axial must be >= 2");
    assert!(n_around >= 3, "n_around must be >= 3");

    let mut positions: Vec<Vec3> = Vec::with_capacity((n_axial + 1) * n_around);
    let dz = length_m / (n_axial as f32);
    let z0 = -length_m / 2.0;
    let dtheta = std::f32::consts::TAU / (n_around as f32);

    for i in 0..=n_axial {
        let z = z0 + (i as f32) * dz;
        for j in 0..n_around {
            let theta = (j as f32) * dtheta;
            positions.push(Vec3::new(radius_m * theta.cos(), radius_m * theta.sin(), z));
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_axial * n_around);
    for i in 0..n_axial {
        for j in 0..n_around {
            let j_next = (j + 1) % n_around;
            let a = (i * n_around + j) as u32;
            let b = (i * n_around + j_next) as u32;
            let c = ((i + 1) * n_around + j_next) as u32;
            let d = ((i + 1) * n_around + j) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
        }
    }

    TriMeshCpu::from_triangles(&positions, &triangles)
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thin_cylinder_triangle_count_matches_formula() {
        let n_axial = 12;
        let n_around = 16;
        let mesh = thin_cylinder(1.0, 0.01, n_axial, n_around);
        // 2 triangles per (axial, around) cell, 3 indices per triangle.
        assert_eq!(mesh.indices.len(), 2 * n_axial * n_around * 3);
        assert_eq!(mesh.vertices.len(), (n_axial + 1) * n_around);
    }

    #[test]
    fn test_thin_cylinder_normals_unit_length() {
        let mesh = thin_cylinder(2.0, 0.05, 8, 24);
        for (i, v) in mesh.vertices.iter().enumerate() {
            let n = Vec3::from(v.normal);
            let mag = n.length();
            assert!(
                (mag - 1.0).abs() < 1e-6,
                "vertex {i} normal magnitude = {mag}"
            );
        }
    }

    #[test]
    fn test_bbox_max_dim_correct() {
        // r = 1, length = 2: bbox is x ∈ [-1, 1], y ∈ [-1, 1], z ∈ [-1, 1].
        // Max extent on any axis = 2.
        let mesh = thin_cylinder(2.0, 1.0, 4, 32);
        assert!(
            (mesh.bbox_max_dim - 2.0).abs() < 1e-5,
            "bbox_max_dim = {}",
            mesh.bbox_max_dim
        );
    }
}
