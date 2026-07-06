//! wgpu compute backend (FP32 WGSL kernels, portable: Vulkan / Metal / DX12).
//!
//! E.1 layout: **arena buffers**. Materials, masks, 12 CPML ψ arrays, and
//! profiles as individual bindings would blow WebGPU's default limit of 8
//! storage buffers per stage, so everything is packed into five arenas
//! (fields, coefficients, ψ, profiles, masks) plus one uniform — see the
//! header of `shaders/fdtd.wgsl` for the exact packing order, which the
//! shader re-derives from the grid dims.
//!
//! Materials are not branched on in the shader: the host materializes four
//! per-cell coefficient maps in f64 — `ca`/`cb` (lossy Taflove CA/CB form,
//! with `ca = 1` for lossless cells so `e = 1·e + cb·curl` reproduces the
//! plain add exactly), `ce_cpml` (`Δt/(ε₀ε_r)`, the coefficient the
//! reference CPML pass uses regardless of σ), and `ch` (`Δt/(μ₀μ_r)`) —
//! then narrows them once to f32.
//!
//! A full step is six fused bulk+CPML dispatches (3 H, then 3 E) plus three
//! mask clamps when masks are attached, all in one compute pass: WebGPU
//! orders storage-buffer writes between dispatches, so no explicit barriers
//! are needed. `step_n` submits in chunks to bound command-buffer size.
//! The FP64 host state is narrowed on upload and widened on readback; gates
//! `compute-002`/`compute-005` bound the accumulated FP32 error against the
//! CPU backend. The PEC box is enforced host-side by zeroing the outer
//! tangential E faces at upload — no kernel ever writes those faces, so the
//! clamp holds for the whole run.

use std::sync::mpsc;

use yee_core::units::{EPS0, MU0};

use crate::cpml::make_profiles;
use crate::cpu::apply_pec_box;
use crate::error::ComputeError;
use crate::fields::Fields;
use crate::materials::{Boundary, Materials};
use crate::spec::{FdtdSpec, len3};

/// Steps encoded per queue submission. Keeps a single command buffer to a
/// few hundred dispatches so watchdog-limited platforms don't kill long runs.
const STEPS_PER_SUBMIT: usize = 64;

const WORKGROUP: u32 = 4;

/// Uniform block mirrored by `struct Params` in `shaders/fdtd.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    nx: u32,
    ny: u32,
    nz: u32,
    npml: u32,
    axes_mask: u32,
    has_cpml: u32,
    has_mask: u32,
    _pad0: u32,
    inv_dx: f32,
    inv_dy: f32,
    inv_dz: f32,
    _pad1: f32,
}

/// GPU FDTD stepper (per-cell materials, CPML or PEC box, interior masks).
#[derive(Debug)]
pub struct GpuFdtd {
    spec: FdtdSpec,
    device: wgpu::Device,
    queue: wgpu::Queue,
    bind_group: wgpu::BindGroup,
    /// Update pipelines in dispatch order: hx, hy, hz, ex, ey, ez.
    update_pipelines: [wgpu::ComputePipeline; 6],
    /// Mask clamp pipelines (ex, ey, ez); dispatched only when `has_mask`.
    clamp_pipelines: [wgpu::ComputePipeline; 3],
    has_mask: bool,
    field_arena: wgpu::Buffer,
    adapter_name: String,
}

impl GpuFdtd {
    /// Build a uniform-vacuum stepper with no boundary phase (raw E.0
    /// semantics), on the first available adapter.
    pub fn new(spec: FdtdSpec, fields: Fields) -> Result<Self, ComputeError> {
        Self::with_config(spec, fields, Materials::default(), Boundary::None)
    }

    /// Build a stepper with per-cell materials / masks and an outer-boundary
    /// treatment (E.1), on the first available adapter (high-performance
    /// preference).
    ///
    /// Returns [`ComputeError::NoAdapter`] when the machine has no
    /// compatible GPU — callers (and the parity gates) treat that as
    /// "skip", not as failure.
    ///
    /// # Panics
    ///
    /// Panics if any field, material, or mask buffer length disagrees with
    /// the spec's staggered shapes (same contract as `CpuFdtd`).
    pub fn with_config(
        spec: FdtdSpec,
        mut fields: Fields,
        materials: Materials,
        boundary: Boundary,
    ) -> Result<Self, ComputeError> {
        let field_lens = [
            len3(spec.ex_dims()),
            len3(spec.ey_dims()),
            len3(spec.ez_dims()),
            len3(spec.hx_dims()),
            len3(spec.hy_dims()),
            len3(spec.hz_dims()),
        ];
        assert_eq!(fields.ex.len(), field_lens[0], "ex length mismatch");
        assert_eq!(fields.ey.len(), field_lens[1], "ey length mismatch");
        assert_eq!(fields.ez.len(), field_lens[2], "ez length mismatch");
        assert_eq!(fields.hx.len(), field_lens[3], "hx length mismatch");
        assert_eq!(fields.hy.len(), field_lens[4], "hy length mismatch");
        assert_eq!(fields.hz.len(), field_lens[5], "hz length mismatch");
        materials.validate(&spec);

        // The PEC box is a host-side invariant: zero the outer tangential E
        // faces once; no kernel ever writes them afterwards.
        let (cpml_config, pec_box) = match boundary {
            Boundary::None => (None, false),
            Boundary::PecBox => (None, true),
            Boundary::Cpml(config) => (Some(config), false),
        };
        if pec_box {
            apply_pec_box(&spec, &mut fields);
        }

        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .map_err(|_| ComputeError::NoAdapter)?;
        let adapter_name = adapter.get_info().name.clone();

        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))
                .map_err(|e| ComputeError::Device(e.to_string()))?;

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("yee-compute fdtd"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/fdtd.wgsl").into()),
        });

        // Bindings: 0 = uniform Params; 1 = fields (rw); 2 = coeffs (r);
        // 3 = psi (rw); 4 = profiles (r); 5 = masks (r).
        let storage = |binding: u32, read_only: bool| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("yee-compute fdtd"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                storage(1, false),
                storage(2, true),
                storage(3, false),
                storage(4, true),
                storage(5, true),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("yee-compute fdtd"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });
        let make_pipeline = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let update_pipelines = [
            make_pipeline("update_hx"),
            make_pipeline("update_hy"),
            make_pipeline("update_hz"),
            make_pipeline("update_ex"),
            make_pipeline("update_ey"),
            make_pipeline("update_ez"),
        ];
        let clamp_pipelines = [
            make_pipeline("clamp_ex"),
            make_pipeline("clamp_ey"),
            make_pipeline("clamp_ez"),
        ];

        let make_storage_buffer = |label: &str, data: &[f32], writable_by_shader: bool| {
            let mut usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
            if writable_by_shader {
                usage |= wgpu::BufferUsages::COPY_SRC;
            }
            let buffer = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: std::mem::size_of_val(data) as u64,
                usage,
                mapped_at_creation: false,
            });
            queue.write_buffer(&buffer, 0, bytemuck::cast_slice(data));
            buffer
        };

        // --- field arena (f64 → f32, packing order ex|ey|ez|hx|hy|hz) ---
        let narrowed: Vec<f32> = [
            &fields.ex, &fields.ey, &fields.ez, &fields.hx, &fields.hy, &fields.hz,
        ]
        .into_iter()
        .flatten()
        .map(|v| *v as f32)
        .collect();
        let field_arena = make_storage_buffer("yee-compute fields", &narrowed, true);

        // --- coefficient arena: ca | cb | ce_cpml | ch, all [nx+1,ny+1,nz+1],
        // materialized in f64 with the exact reference formulas ---
        let n_cells = (spec.nx + 1) * (spec.ny + 1) * (spec.nz + 1);
        let mut coeffs = vec![0.0f32; 4 * n_cells];
        let (ca_arr, rest) = coeffs.split_at_mut(n_cells);
        let (cb_arr, rest) = rest.split_at_mut(n_cells);
        let (ce_arr, ch_arr) = rest.split_at_mut(n_cells);
        let eps_cells = materials.eps_r_cells.as_deref();
        let sig_cells = materials.sigma_cells.as_deref();
        let mu_cells = materials.mu_r_cells.as_deref();
        for cell in 0..n_cells {
            let eps_r = eps_cells.map_or(spec.eps_r, |e| e[cell]);
            let mu_r = mu_cells.map_or(spec.mu_r, |m| m[cell]);
            let (ca, cb) = match sig_cells {
                None => (1.0, spec.dt / (EPS0 * eps_r)),
                Some(sig) => {
                    let denom = 2.0 * EPS0 * eps_r + sig[cell] * spec.dt;
                    (
                        (2.0 * EPS0 * eps_r - sig[cell] * spec.dt) / denom,
                        2.0 * spec.dt / denom,
                    )
                }
            };
            ca_arr[cell] = ca as f32;
            cb_arr[cell] = cb as f32;
            ce_arr[cell] = (spec.dt / (EPS0 * eps_r)) as f32;
            ch_arr[cell] = (spec.dt / (MU0 * mu_r)) as f32;
        }
        let coeff_arena = make_storage_buffer("yee-compute coeffs", &coeffs, false);

        // --- ψ arena + profile buffer (dummies when CPML is off) ---
        let (psi_data, profile_data, npml, axes_mask) = match cpml_config {
            None => (vec![0.0f32; 1], vec![0.0f32; 1], 0u32, 0u32),
            Some(config) => {
                let psi_len = 2 * (field_lens[0] + field_lens[1] + field_lens[2])
                    + 2 * (field_lens[3] + field_lens[4] + field_lens[5]);
                let ((b, c, kappa), (b_h, c_h, kappa_h)) = make_profiles(&config, spec.dt);
                let profiles: Vec<f32> = [b, c, kappa, b_h, c_h, kappa_h]
                    .into_iter()
                    .flatten()
                    .map(|v| v as f32)
                    .collect();
                let axes_mask = config
                    .axes
                    .iter()
                    .enumerate()
                    .fold(0u32, |m, (i, &on)| if on { m | (1 << i) } else { m });
                (
                    vec![0.0f32; psi_len],
                    profiles,
                    config.npml as u32,
                    axes_mask,
                )
            }
        };
        let psi_arena = make_storage_buffer("yee-compute psi", &psi_data, false);
        let profile_buffer = make_storage_buffer("yee-compute profiles", &profile_data, false);

        // --- mask arena: ex | ey | ez as u32 (dummy when no masks) ---
        let has_mask = materials.has_mask();
        let mask_data: Vec<u32> = if has_mask {
            let mut data = Vec::with_capacity(field_lens[0] + field_lens[1] + field_lens[2]);
            for (mask, len) in [
                (&materials.pec_mask_ex, field_lens[0]),
                (&materials.pec_mask_ey, field_lens[1]),
                (&materials.pec_mask_ez, field_lens[2]),
            ] {
                match mask {
                    None => data.extend(std::iter::repeat_n(0u32, len)),
                    Some(m) => data.extend(m.iter().map(|&b| u32::from(b))),
                }
            }
            data
        } else {
            vec![0u32; 1]
        };
        let mask_arena = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-compute masks"),
            size: std::mem::size_of_val(mask_data.as_slice()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&mask_arena, 0, bytemuck::cast_slice(&mask_data));

        // --- params ---
        let params = Params {
            nx: spec.nx as u32,
            ny: spec.ny as u32,
            nz: spec.nz as u32,
            npml,
            axes_mask,
            has_cpml: u32::from(cpml_config.is_some()),
            has_mask: u32::from(has_mask),
            _pad0: 0,
            inv_dx: (1.0 / spec.dx) as f32,
            inv_dy: (1.0 / spec.dy) as f32,
            inv_dz: (1.0 / spec.dz) as f32,
            _pad1: 0.0,
        };
        let params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-compute params"),
            size: std::mem::size_of::<Params>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&params_buffer, 0, bytemuck::bytes_of(&params));

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("yee-compute fdtd"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: field_arena.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: coeff_arena.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: psi_arena.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: profile_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: mask_arena.as_entire_binding(),
                },
            ],
        });

        Ok(Self {
            spec,
            device,
            queue,
            bind_group,
            update_pipelines,
            clamp_pipelines,
            has_mask,
            field_arena,
            adapter_name,
        })
    }

    /// The problem description this stepper was built from.
    pub fn spec(&self) -> &FdtdSpec {
        &self.spec
    }

    /// Adapter name (diagnostics; e.g. printed by the parity gates).
    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    /// Dispatch extents per update pipeline (hx, hy, hz, ex, ey, ez), and
    /// per clamp pipeline (the E extents, last three).
    fn dispatch_extents(&self) -> [(usize, usize, usize); 6] {
        [
            self.spec.hx_dims(),
            self.spec.hy_dims(),
            self.spec.hz_dims(),
            self.spec.ex_dims(),
            self.spec.ey_dims(),
            self.spec.ez_dims(),
        ]
    }

    /// Advance the state by `n` leapfrog steps (3 H, 3 E, then the mask
    /// clamps when masks are attached), submitting in [`STEPS_PER_SUBMIT`]
    /// chunks.
    pub fn step_n(&mut self, n: usize) -> Result<(), ComputeError> {
        let extents = self.dispatch_extents();
        let groups = |len: usize| (len as u32).div_ceil(WORKGROUP);
        let mut remaining = n;
        while remaining > 0 {
            let chunk = remaining.min(STEPS_PER_SUBMIT);
            remaining -= chunk;
            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("yee-compute step"),
                });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("yee-compute step"),
                    timestamp_writes: None,
                });
                pass.set_bind_group(0, &self.bind_group, &[]);
                for _ in 0..chunk {
                    for (pipeline, extent) in self.update_pipelines.iter().zip(extents) {
                        pass.set_pipeline(pipeline);
                        pass.dispatch_workgroups(
                            groups(extent.0),
                            groups(extent.1),
                            groups(extent.2),
                        );
                    }
                    if self.has_mask {
                        for (pipeline, extent) in self.clamp_pipelines.iter().zip(&extents[3..6]) {
                            pass.set_pipeline(pipeline);
                            pass.dispatch_workgroups(
                                groups(extent.0),
                                groups(extent.1),
                                groups(extent.2),
                            );
                        }
                    }
                }
            }
            self.queue.submit(Some(encoder.finish()));
        }
        Ok(())
    }

    /// Copy all six components back to the host, widened to FP64.
    pub fn read_fields(&mut self) -> Result<Fields, ComputeError> {
        let widened = self.read_field_arena()?;
        let lens = [
            len3(self.spec.ex_dims()),
            len3(self.spec.ey_dims()),
            len3(self.spec.ez_dims()),
            len3(self.spec.hx_dims()),
            len3(self.spec.hy_dims()),
            len3(self.spec.hz_dims()),
        ];
        let mut iter = widened.into_iter();
        let mut take = |len: usize| iter.by_ref().take(len).collect::<Vec<f64>>();
        Ok(Fields {
            ex: take(lens[0]),
            ey: take(lens[1]),
            ez: take(lens[2]),
            hx: take(lens[3]),
            hy: take(lens[4]),
            hz: take(lens[5]),
        })
    }

    fn read_field_arena(&self) -> Result<Vec<f64>, ComputeError> {
        let size = self.field_arena.size();
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-compute staging"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("yee-compute readback"),
            });
        encoder.copy_buffer_to_buffer(&self.field_arena, 0, &staging, 0, size);
        self.queue.submit(Some(encoder.finish()));

        let (tx, rx) = mpsc::channel();
        staging
            .slice(..)
            .map_async(wgpu::MapMode::Read, move |result| {
                let _ = tx.send(result);
            });
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .map_err(|e| ComputeError::Readback(e.to_string()))?;
        rx.recv()
            .map_err(|e| ComputeError::Readback(e.to_string()))?
            .map_err(|e| ComputeError::Readback(e.to_string()))?;

        let widened = {
            let view = staging.slice(..).get_mapped_range();
            bytemuck::cast_slice::<u8, f32>(&view)
                .iter()
                .map(|v| f64::from(*v))
                .collect()
        };
        staging.unmap();
        Ok(widened)
    }
}
