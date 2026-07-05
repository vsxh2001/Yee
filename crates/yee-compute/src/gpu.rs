//! wgpu compute backend (FP32 WGSL kernels, portable: Vulkan / Metal / DX12).
//!
//! Field storage is FP32 (industry-standard for FDTD; consumer-GPU FP64 runs
//! at 1/32–1/64 rate). The FP64 host state is narrowed on upload and widened
//! on readback; gate `compute-002` bounds the accumulated FP32 error against
//! the CPU backend. One WGSL module carries six entry points — one per
//! staggered component — and a full step is six dispatches (3 H, then 3 E)
//! in a single compute pass: WebGPU orders storage-buffer writes between
//! dispatches (each is its own usage scope), so no explicit barriers are
//! needed. `step_n` submits in chunks to bound command-buffer size and avoid
//! device timeouts on long runs.

use std::sync::mpsc;

use yee_core::units::{EPS0, MU0};

use crate::error::ComputeError;
use crate::fields::Fields;
use crate::spec::{FdtdSpec, len3};

/// Steps encoded per queue submission. Keeps a single command buffer to a few
/// hundred dispatches so watchdog-limited platforms don't kill long runs.
const STEPS_PER_SUBMIT: usize = 64;

const WORKGROUP: u32 = 4;

/// Uniform block mirrored by `struct Params` in `shaders/fdtd.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    nx: u32,
    ny: u32,
    nz: u32,
    _pad0: u32,
    ch: f32,
    ce: f32,
    inv_dx: f32,
    inv_dy: f32,
    inv_dz: f32,
    _pad1: f32,
    _pad2: f32,
    _pad3: f32,
}

/// GPU FDTD stepper (uniform lossless vacuum, PEC box).
#[derive(Debug)]
pub struct GpuFdtd {
    spec: FdtdSpec,
    device: wgpu::Device,
    queue: wgpu::Queue,
    bind_group: wgpu::BindGroup,
    /// Pipelines in dispatch order: hx, hy, hz, ex, ey, ez.
    pipelines: [wgpu::ComputePipeline; 6],
    /// Field buffers in the same order as the WGSL bindings 1–6:
    /// ex, ey, ez, hx, hy, hz.
    field_buffers: [wgpu::Buffer; 6],
    adapter_name: String,
}

impl GpuFdtd {
    /// Build a stepper on the first available adapter (high-performance
    /// preference), uploading `fields` as the initial state.
    ///
    /// Returns [`ComputeError::NoAdapter`] when the machine has no compatible
    /// GPU — callers (and the `compute-002` gate) treat that as "skip", not
    /// as failure.
    ///
    /// # Panics
    ///
    /// Panics if any buffer length disagrees with the spec's staggered
    /// shapes (same contract as `CpuFdtd::new`).
    pub fn new(spec: FdtdSpec, fields: Fields) -> Result<Self, ComputeError> {
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

        // Bindings: 0 = uniform Params, 1..=6 = ex, ey, ez, hx, hy, hz.
        let layout_entries: Vec<wgpu::BindGroupLayoutEntry> = (0..7)
            .map(|binding| wgpu::BindGroupLayoutEntry {
                binding,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: if binding == 0 {
                        wgpu::BufferBindingType::Uniform
                    } else {
                        wgpu::BufferBindingType::Storage { read_only: false }
                    },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            })
            .collect();
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("yee-compute fdtd"),
            entries: &layout_entries,
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
        let pipelines = [
            make_pipeline("update_hx"),
            make_pipeline("update_hy"),
            make_pipeline("update_hz"),
            make_pipeline("update_ex"),
            make_pipeline("update_ey"),
            make_pipeline("update_ez"),
        ];

        // Coefficients in f64 (identical formulas to the CPU backend), then
        // narrowed once for the shader.
        let params = Params {
            nx: spec.nx as u32,
            ny: spec.ny as u32,
            nz: spec.nz as u32,
            _pad0: 0,
            ch: (spec.dt / (MU0 * spec.mu_r)) as f32,
            ce: (spec.dt / (EPS0 * spec.eps_r)) as f32,
            inv_dx: (1.0 / spec.dx) as f32,
            inv_dy: (1.0 / spec.dy) as f32,
            inv_dz: (1.0 / spec.dz) as f32,
            _pad1: 0.0,
            _pad2: 0.0,
            _pad3: 0.0,
        };
        let params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-compute params"),
            size: std::mem::size_of::<Params>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&params_buffer, 0, bytemuck::bytes_of(&params));

        let host_buffers: [(&[f64], usize, &str); 6] = [
            (&fields.ex, len3(spec.ex_dims()), "ex"),
            (&fields.ey, len3(spec.ey_dims()), "ey"),
            (&fields.ez, len3(spec.ez_dims()), "ez"),
            (&fields.hx, len3(spec.hx_dims()), "hx"),
            (&fields.hy, len3(spec.hy_dims()), "hy"),
            (&fields.hz, len3(spec.hz_dims()), "hz"),
        ];
        let field_buffers = host_buffers.map(|(host, expected_len, name)| {
            assert_eq!(host.len(), expected_len, "{name} length mismatch");
            let narrowed: Vec<f32> = host.iter().map(|v| *v as f32).collect();
            let buffer = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(name),
                size: (narrowed.len() * std::mem::size_of::<f32>()) as u64,
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST
                    | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });
            queue.write_buffer(&buffer, 0, bytemuck::cast_slice(&narrowed));
            buffer
        });

        let bind_entries: Vec<wgpu::BindGroupEntry> = std::iter::once(&params_buffer)
            .chain(field_buffers.iter())
            .enumerate()
            .map(|(binding, buffer)| wgpu::BindGroupEntry {
                binding: binding as u32,
                resource: buffer.as_entire_binding(),
            })
            .collect();
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("yee-compute fdtd"),
            layout: &bind_group_layout,
            entries: &bind_entries,
        });

        Ok(Self {
            spec,
            device,
            queue,
            bind_group,
            pipelines,
            field_buffers,
            adapter_name,
        })
    }

    /// The problem description this stepper was built from.
    pub fn spec(&self) -> &FdtdSpec {
        &self.spec
    }

    /// Adapter name (diagnostics; e.g. printed by the parity gate).
    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    /// Dispatch extents per pipeline, in `pipelines` order.
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

    /// Advance the state by `n` leapfrog steps (3 H dispatches, then 3 E,
    /// per step), submitting in [`STEPS_PER_SUBMIT`] chunks.
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
                    for (pipeline, extent) in self.pipelines.iter().zip(extents) {
                        pass.set_pipeline(pipeline);
                        pass.dispatch_workgroups(
                            groups(extent.0),
                            groups(extent.1),
                            groups(extent.2),
                        );
                    }
                }
            }
            self.queue.submit(Some(encoder.finish()));
        }
        Ok(())
    }

    /// Copy all six components back to the host, widened to FP64.
    pub fn read_fields(&mut self) -> Result<Fields, ComputeError> {
        let mut out: [Vec<f64>; 6] = Default::default();
        for (slot, buffer) in out.iter_mut().zip(self.field_buffers.iter()) {
            *slot = self.read_buffer(buffer)?;
        }
        let [ex, ey, ez, hx, hy, hz] = out;
        Ok(Fields {
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        })
    }

    fn read_buffer(&self, buffer: &wgpu::Buffer) -> Result<Vec<f64>, ComputeError> {
        let size = buffer.size();
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
        encoder.copy_buffer_to_buffer(buffer, 0, &staging, 0, size);
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
