//! Wgpu-rendered 3D triangle-mesh viewport for the Phase 1.gui.1 walking-skeleton shell.
//!
//! This module provides:
//!
//! - A CPU-side mesh representation ([`TriMeshCpu`]) with positions + per-vertex
//!   normals, suitable for direct upload into a `wgpu::Buffer`.
//! - A hand-coded thin-cylinder mesher ([`thin_cylinder`]) used as the default
//!   walking-skeleton geometry — same triangulation strategy as the mom-001
//!   dipole fixture, copied here (per the Track-D dispatch) so this crate does
//!   not take a dev-fixture dependency on `yee-mom`.
//! - An orbit camera state + wireframe toggle ([`ViewportState`]).
//! - An [`egui_wgpu::CallbackTrait`] implementation ([`MeshCallback`]) that
//!   builds GPU resources lazily on the first `prepare` call and draws an
//!   indexed triangle list with flat-shaded diffuse lighting.
//!
//! ## Scope
//!
//! Phase 1.gui.1 is intentionally a walking skeleton: one fixed mesh, no
//! picking, no current-density colouring, no scene management. The viewport's
//! sole job is to prove the egui → wgpu paint-callback bridge end-to-end.
//!
//! See `crate::app` for how this module is wired into the dock layout.

use bytemuck::{Pod, Zeroable};
use egui_wgpu::wgpu;
use glam::{Mat4, Vec3};

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
// Viewport / camera state
// ----------------------------------------------------------------------------

/// CPU-side state for the Mesh 3D tab: orbit-camera angles + the active mesh.
///
/// Defaults (set by [`ViewportState::new`]): yaw=30°, pitch=-20°,
/// distance = `mesh.bbox_max_dim · 3` (clamped to ≥ 0.1 m so a degenerate mesh
/// can't drop the camera onto the origin). Wireframe overlay is off by default.
#[derive(Debug, Clone)]
pub struct ViewportState {
    /// Camera yaw in degrees (around world +Z).
    pub camera_yaw_deg: f32,
    /// Camera pitch in degrees (positive = look down, negative = look up).
    pub camera_pitch_deg: f32,
    /// Camera distance from the origin in metres.
    pub camera_dist: f32,
    /// If `true`, the fragment shader draws an additional wireframe overlay.
    pub wireframe: bool,
    /// Active mesh.
    pub mesh: TriMeshCpu,
}

impl ViewportState {
    /// Build a fresh viewport state around the given mesh, with the
    /// default orbit-camera orientation (yaw=30°, pitch=-20°, dist=3·bbox).
    pub fn new(mesh: TriMeshCpu) -> Self {
        let dist = (mesh.bbox_max_dim * 3.0).max(0.1);
        Self {
            camera_yaw_deg: 30.0,
            camera_pitch_deg: -20.0,
            camera_dist: dist,
            wireframe: false,
            mesh,
        }
    }

    /// Compute the world-space camera position from yaw/pitch/dist.
    pub fn camera_position(&self) -> Vec3 {
        let yaw = self.camera_yaw_deg.to_radians();
        let pitch = self.camera_pitch_deg.to_radians();
        let cp = pitch.cos();
        Vec3::new(
            self.camera_dist * cp * yaw.cos(),
            self.camera_dist * cp * yaw.sin(),
            self.camera_dist * pitch.sin(),
        )
    }

    /// Build the combined view-projection matrix for the current camera.
    pub fn view_proj(&self, aspect: f32) -> Mat4 {
        let eye = self.camera_position();
        let view = Mat4::look_at_rh(eye, Vec3::ZERO, Vec3::Z);
        let proj = Mat4::perspective_rh(
            45f32.to_radians(),
            aspect.max(0.01),
            (self.camera_dist * 0.01).max(1e-4),
            (self.camera_dist * 10.0).max(1.0),
        );
        proj * view
    }
}

// ----------------------------------------------------------------------------
// GPU resources + paint callback
// ----------------------------------------------------------------------------

/// Uniform block passed to both vertex and fragment stages: combined MVP,
/// world-space camera position (for the directional light), and a wireframe
/// toggle.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct Uniforms {
    mvp: [[f32; 4]; 4],
    camera_pos: [f32; 4], // xyz + 1 padding lane
    wireframe: [f32; 4],  // x = 0/1, yzw padding
}

/// GPU resources for the triangle-mesh viewport. Allocated lazily on the first
/// paint and reused thereafter; the uniform buffer is rewritten each frame
/// because it's tiny, and the vertex/index buffers are static for the lifetime
/// of the (currently fixed) mesh.
pub struct MeshRenderResources {
    /// Vertex buffer (positions + normals).
    pub vertex_buf: wgpu::Buffer,
    /// Index buffer (`u32`).
    pub index_buf: wgpu::Buffer,
    /// Uniform buffer (one [`Uniforms`] block).
    pub uniform_buf: wgpu::Buffer,
    /// Bind group binding the uniform buffer.
    pub bind_group: wgpu::BindGroup,
    /// Render pipeline (triangle list, depth-less, alpha-blended).
    pub pipeline: wgpu::RenderPipeline,
    /// Number of indices in [`Self::index_buf`].
    pub n_indices: u32,
}

impl MeshRenderResources {
    fn new(device: &wgpu::Device, target_format: wgpu::TextureFormat, mesh: &TriMeshCpu) -> Self {
        use wgpu::util::DeviceExt;

        let vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("yee-gui.viewport.vertices"),
            contents: bytemuck::cast_slice(&mesh.vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });

        let index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("yee-gui.viewport.indices"),
            contents: bytemuck::cast_slice(&mesh.indices),
            usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
        });

        let uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-gui.viewport.uniforms"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("yee-gui.viewport.bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("yee-gui.viewport.bg"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buf.as_entire_binding(),
            }],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("yee-gui.viewport.shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("yee-gui.viewport.pl"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            // wgpu 29 renamed `push_constant_ranges` to `immediate_size`
            // (and changed semantics: bytes of immediate data, not ranges).
            // No immediate data is used by this pipeline.
            immediate_size: 0,
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("yee-gui.viewport.pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                buffers: &[Vertex::layout()],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                // No back-face culling: the cylinder mesh is open at both ends
                // and the camera will see inside surfaces while orbiting.
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            // wgpu 29 renamed `multiview: Option<NonZeroU32>` to
            // `multiview_mask: Option<NonZeroU32>`; same semantics, new name.
            multiview_mask: None,
            cache: None,
        });

        Self {
            vertex_buf,
            index_buf,
            uniform_buf,
            bind_group,
            pipeline,
            n_indices: mesh.indices.len() as u32,
        }
    }
}

/// Egui paint callback that uploads per-frame uniforms and draws the mesh.
///
/// Lifetime: the callback is constructed every frame inside `TabViewer::ui`,
/// stamped with the current viewport state snapshot (camera + wireframe
/// toggle) and the pixel size of the paint rect. The CPU mesh is cloned each
/// frame; the GPU resources are cached inside `egui_wgpu::CallbackResources`
/// across frames so the actual upload happens once at startup.
pub struct MeshCallback {
    /// Snapshot of the mesh — cloned each frame because the callback is owned
    /// by the egui paint list and lives until egui hands it back to wgpu.
    pub mesh: TriMeshCpu,
    /// MVP for this frame.
    pub mvp: Mat4,
    /// World-space camera position (for the directional light).
    pub camera_pos: Vec3,
    /// `true` if the wireframe overlay should be drawn this frame.
    pub wireframe: bool,
}

impl egui_wgpu::CallbackTrait for MeshCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        _screen_descriptor: &egui_wgpu::ScreenDescriptor,
        _egui_encoder: &mut wgpu::CommandEncoder,
        resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        // Lazily build GPU resources the first time we paint.
        if resources.get::<MeshRenderResources>().is_none() {
            // We need the surface format to build the pipeline; egui_wgpu
            // doesn't pass it to `prepare`, so we use the eframe desktop
            // default of `Bgra8UnormSrgb`. If eframe ever switches the
            // surface format (e.g. for HDR), the pipeline rebuild will
            // surface as a wgpu validation error at the first draw, which is
            // the right place to handle it.
            let target_format = wgpu::TextureFormat::Bgra8UnormSrgb;
            let res = MeshRenderResources::new(device, target_format, &self.mesh);
            resources.insert(res);
        }

        if let Some(res) = resources.get::<MeshRenderResources>() {
            let uniforms = Uniforms {
                mvp: self.mvp.to_cols_array_2d(),
                camera_pos: [self.camera_pos.x, self.camera_pos.y, self.camera_pos.z, 1.0],
                wireframe: [if self.wireframe { 1.0 } else { 0.0 }, 0.0, 0.0, 0.0],
            };
            queue.write_buffer(&res.uniform_buf, 0, bytemuck::bytes_of(&uniforms));
        }
        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(res) = resources.get::<MeshRenderResources>() else {
            return;
        };
        render_pass.set_pipeline(&res.pipeline);
        render_pass.set_bind_group(0, &res.bind_group, &[]);
        render_pass.set_vertex_buffer(0, res.vertex_buf.slice(..));
        render_pass.set_index_buffer(res.index_buf.slice(..), wgpu::IndexFormat::Uint32);
        render_pass.draw_indexed(0..res.n_indices, 0, 0..1);
    }
}

// ----------------------------------------------------------------------------
// Shader
// ----------------------------------------------------------------------------

/// WGSL shader for the triangle-mesh pipeline.
///
/// The vertex stage transforms `position` by the MVP and passes a world-space
/// normal + position to the fragment stage. The fragment stage computes
/// Lambertian diffuse against a directional light aligned with the camera
/// (so the surface is always lit from the viewer's perspective), with an
/// ambient term of 0.2. When the `wireframe` uniform is non-zero, silhouette
/// fragments are darkened — a coarse stand-in for true wireframe rendering,
/// which would need a second draw call with `PolygonMode::Line` (left out of
/// Phase 1.gui.1 because wgpu requires the `POLYGON_MODE_LINE` feature flag
/// which isn't available on every backend).
const SHADER: &str = r#"
struct Uniforms {
    mvp: mat4x4<f32>,
    camera_pos: vec4<f32>,
    wireframe: vec4<f32>,
};

@group(0) @binding(0) var<uniform> u: Uniforms;

struct VsIn {
    @location(0) position: vec3<f32>,
    @location(1) normal:   vec3<f32>,
};

struct VsOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_pos: vec3<f32>,
    @location(1) normal:    vec3<f32>,
};

@vertex
fn vs_main(in: VsIn) -> VsOut {
    var out: VsOut;
    out.clip_pos = u.mvp * vec4<f32>(in.position, 1.0);
    out.world_pos = in.position;
    out.normal = in.normal;
    return out;
}

@fragment
fn fs_main(in: VsOut) -> @location(0) vec4<f32> {
    // Directional light from the camera toward the origin.
    let light_dir = normalize(u.camera_pos.xyz);
    let n = normalize(in.normal);
    let ambient = 0.2;
    let diffuse = max(dot(n, light_dir), 0.0);
    let intensity = clamp(ambient + (1.0 - ambient) * diffuse, 0.0, 1.0);
    let base = vec3<f32>(0.72, 0.78, 0.92);
    var color = base * intensity;
    if (u.wireframe.x > 0.5) {
        // Coarse wireframe overlay: emphasise silhouettes by darkening
        // grazing-angle fragments so facet edges read more clearly.
        let edge = 1.0 - abs(dot(n, light_dir));
        color = mix(color, vec3<f32>(0.05, 0.05, 0.08), edge);
    }
    return vec4<f32>(color, 1.0);
}
"#;

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

    #[test]
    fn test_camera_position_default_orientation() {
        let mesh = thin_cylinder(1.0, 0.01, 4, 8);
        let s = ViewportState::new(mesh);
        let p = s.camera_position();
        // Sanity: distance from origin matches camera_dist.
        let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        assert!(
            (r - s.camera_dist).abs() < 1e-5,
            "camera distance drift: {r} vs {}",
            s.camera_dist
        );
    }
}
