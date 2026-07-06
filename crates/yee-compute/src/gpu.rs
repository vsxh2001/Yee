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
use crate::dispersive::{DispersiveMap, ade_coeffs};
use crate::drive::{Drive, EComponent};
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
    has_dispersion: u32,
    inv_dx: f32,
    inv_dy: f32,
    inv_dz: f32,
    /// ω·Δt for the on-GPU NTFF DFT accumulator (E.5b); 0.0 disables it.
    dft_omega_dt: f32,
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
    /// Drive pipelines: inject_soft, apply_ports, record_probes, bump_step.
    drive_pipelines: [wgpu::ComputePipeline; 4],
    has_mask: bool,
    dft_enabled: bool,
    /// accumulate_dft pipeline (E.5b).
    dft_pipeline: wgpu::ComputePipeline,
    field_arena: wgpu::Buffer,
    psi_arena: wgpu::Buffer,
    drive_data: wgpu::Buffer,
    drive_counts: (usize, usize, usize), // (n_soft, n_ports, n_probes)
    max_steps: usize,
    steps_taken: usize,
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
        fields: Fields,
        materials: Materials,
        boundary: Boundary,
    ) -> Result<Self, ComputeError> {
        Self::with_drive(spec, fields, materials, boundary, Drive::default(), 0)
    }

    /// Build a driven stepper (E.2). All drive amplitudes for up to
    /// `max_steps` steps are precomputed host-side (in f64, then narrowed)
    /// and uploaded once, so stepping stays fully chunked with zero per-step
    /// host round-trips; an on-GPU step counter indexes the tables and the
    /// probe output region. `step_n` beyond `max_steps` panics when a drive
    /// is attached.
    ///
    /// # Panics
    ///
    /// Panics on shape mismatches, out-of-bounds drive cells, or a non-empty
    /// drive with `max_steps == 0`.
    pub fn with_drive(
        spec: FdtdSpec,
        fields: Fields,
        materials: Materials,
        boundary: Boundary,
        drive: Drive,
        max_steps: usize,
    ) -> Result<Self, ComputeError> {
        Self::build(
            spec, fields, materials, boundary, drive, max_steps, None, None,
        )
    }

    /// Build a stepper with the on-GPU NTFF DFT accumulator enabled (E.5b):
    /// after every completed step the six fields are accumulated into a
    /// running full-field DFT phasor at `f_probe_hz` (phase from the on-GPU
    /// step counter, so stepping stays fully chunked). Read the phasor pair
    /// back once with [`GpuFdtd::read_dft_fields`] and feed it to the
    /// reference `NtffState` via two synthetic samples (its accumulation is
    /// linear). Requires `max_steps > 0` (the step counter lives in the
    /// drive buffer).
    ///
    /// # Panics
    ///
    /// Panics if `f_probe_hz` is non-positive, `max_steps == 0`, or on any
    /// shape mismatch.
    #[allow(clippy::too_many_arguments)]
    pub fn with_ntff_dft(
        spec: FdtdSpec,
        fields: Fields,
        materials: Materials,
        boundary: Boundary,
        drive: Drive,
        max_steps: usize,
        f_probe_hz: f64,
    ) -> Result<Self, ComputeError> {
        assert!(
            f_probe_hz > 0.0 && f_probe_hz.is_finite(),
            "f_probe_hz must be positive"
        );
        assert!(max_steps > 0, "with_ntff_dft needs max_steps > 0");
        Self::build(
            spec,
            fields,
            materials,
            boundary,
            drive,
            max_steps,
            None,
            Some(f_probe_hz),
        )
    }

    /// Build a dispersive stepper (E.5c): the E half-step runs the unified
    /// ADE update (coefficients precomputed host-side in f64 into the coeff
    /// arena; aux state appended to the psi arena). Exclusive with per-cell
    /// eps/sigma maps and with CPML (the reference dispersive path carries
    /// neither).
    ///
    /// # Panics
    ///
    /// Panics on the exclusivity violations above or any shape mismatch.
    pub fn with_dispersive(
        spec: FdtdSpec,
        fields: Fields,
        materials: Materials,
        boundary: Boundary,
        drive: Drive,
        max_steps: usize,
        map: &DispersiveMap,
    ) -> Result<Self, ComputeError> {
        assert!(
            materials.eps_r_cells.is_none() && materials.sigma_cells.is_none(),
            "dispersive map is exclusive with per-cell eps/sigma (reference semantics)"
        );
        assert!(
            !matches!(boundary, Boundary::Cpml(_)),
            "GPU dispersion is exclusive with CPML in E.5c (use PecBox/None)"
        );
        map.validate(&spec);
        Self::build(
            spec,
            fields,
            materials,
            boundary,
            drive,
            max_steps,
            Some(map),
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        spec: FdtdSpec,
        mut fields: Fields,
        materials: Materials,
        boundary: Boundary,
        drive: Drive,
        max_steps: usize,
        dispersive: Option<&DispersiveMap>,
        dft_f_probe: Option<f64>,
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
        drive.validate(&spec);
        assert!(
            drive.is_empty() || max_steps > 0,
            "GpuFdtd::with_drive: a non-empty drive needs max_steps > 0"
        );

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
                storage(6, true),
                storage(7, false),
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
        let drive_pipelines = [
            make_pipeline("inject_soft"),
            make_pipeline("apply_ports"),
            make_pipeline("record_probes"),
            make_pipeline("bump_step"),
        ];
        let dft_pipeline = make_pipeline("accumulate_dft");

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
        if let Some(map) = dispersive {
            // Append the six unified-ADE maps: ce | c0 | c1 | c2 | d_new | d_old.
            coeffs.reserve(6 * n_cells);
            let mut push_map = |f: &dyn Fn(crate::dispersive::AdeCoeffs) -> f64| {
                for cell in 0..n_cells {
                    coeffs.push(f(ade_coeffs(map.cells[cell], spec.dt)) as f32);
                }
            };
            push_map(&|c| c.ce);
            push_map(&|c| c.c0);
            push_map(&|c| c.c1);
            push_map(&|c| c.c2);
            push_map(&|c| c.q);
            push_map(&|c| c.s);
        }
        let coeff_arena = make_storage_buffer("yee-compute coeffs", &coeffs, false);

        // --- ψ arena + profile buffer (dummies when CPML is off) ---
        let (mut psi_data, profile_data, npml, axes_mask) = match cpml_config {
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
        if dispersive.is_some() {
            // Six aux1/aux2 maps appended after the CPML block (or dummy).
            psi_data.extend(std::iter::repeat_n(0.0f32, 6 * n_cells));
        }
        if dft_f_probe.is_some() {
            // Full-field DFT phasor pair (re | im) after everything else.
            let total: usize = field_lens.iter().sum();
            psi_data.extend(std::iter::repeat_n(0.0f32, 2 * total));
        }
        let psi_arena = make_storage_buffer("yee-compute psi", &psi_data, true);
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

        // --- drive buffers: index table + state/amplitude/probe data.
        // Layout mirrored by the drv_* helpers in fdtd.wgsl. ---
        let (n_soft, n_ports, n_probes) = (
            drive.soft_sources.len(),
            drive.ports.len(),
            drive.probes.len(),
        );
        let mut drv_idx: Vec<u32> = vec![
            n_soft as u32,
            n_ports as u32,
            n_probes as u32,
            max_steps as u32,
        ];
        for s in &drive.soft_sources {
            drv_idx
                .push((s.component.arena_offset(&spec) + s.component.flat(&spec, s.cell)) as u32);
        }
        for port in &drive.ports {
            drv_idx.push(
                (EComponent::Ez.arena_offset(&spec) + EComponent::Ez.flat(&spec, port.cell)) as u32,
            );
        }
        for probe in &drive.probes {
            drv_idx.push(
                (probe.component.arena_offset(&spec) + probe.component.flat(&spec, probe.cell))
                    as u32,
            );
        }
        let drive_index = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-compute drive index"),
            size: std::mem::size_of_val(drv_idx.as_slice()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&drive_index, 0, bytemuck::cast_slice(&drv_idx));

        // drv_data: [counter, port state ×n, alpha ×n, gamma ×n,
        //            amps (max_steps × n_soft), vsrc (max_steps × n_ports),
        //            probes (max_steps × n_probes)] — all f64-precomputed.
        let mut drv_data: Vec<f32> =
            Vec::with_capacity(1 + 3 * n_ports + max_steps * (n_soft + n_ports + n_probes));
        drv_data.push(0.0); // step counter
        drv_data.extend(std::iter::repeat_n(0.0, n_ports)); // e_z_prev
        let area = spec.dx * spec.dy;
        for port in &drive.ports {
            let alpha = spec.dt * spec.dz / (2.0 * EPS0 * port.resistance * area);
            drv_data.push(alpha as f32);
        }
        for port in &drive.ports {
            let gamma = spec.dt / (EPS0 * port.resistance * area);
            drv_data.push(gamma as f32);
        }
        for n in 0..max_steps {
            for s in &drive.soft_sources {
                drv_data.push(s.waveform.value(n, spec.dt) as f32);
            }
        }
        for n in 0..max_steps {
            for port in &drive.ports {
                drv_data.push(port.waveform.value(n, spec.dt) as f32);
            }
        }
        drv_data.extend(std::iter::repeat_n(0.0, max_steps * n_probes));
        let drive_data = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("yee-compute drive data"),
            size: std::mem::size_of_val(drv_data.as_slice()) as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        queue.write_buffer(&drive_data, 0, bytemuck::cast_slice(&drv_data));

        // --- params ---
        let params = Params {
            nx: spec.nx as u32,
            ny: spec.ny as u32,
            nz: spec.nz as u32,
            npml,
            axes_mask,
            has_cpml: u32::from(cpml_config.is_some()),
            has_mask: u32::from(has_mask),
            has_dispersion: u32::from(dispersive.is_some()),
            inv_dx: (1.0 / spec.dx) as f32,
            inv_dy: (1.0 / spec.dy) as f32,
            inv_dz: (1.0 / spec.dz) as f32,
            dft_omega_dt: dft_f_probe.map_or(0.0, |f| (std::f64::consts::TAU * f * spec.dt) as f32),
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
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: drive_index.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: drive_data.as_entire_binding(),
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
            drive_pipelines,
            has_mask,
            dft_enabled: dft_f_probe.is_some(),
            dft_pipeline,
            field_arena,
            psi_arena,
            drive_data,
            drive_counts: (n_soft, n_ports, n_probes),
            max_steps,
            steps_taken: 0,
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
        let (n_soft, n_ports, n_probes) = self.drive_counts;
        let has_drive = n_soft + n_ports + n_probes > 0 || self.dft_enabled;
        if has_drive {
            assert!(
                self.steps_taken + n <= self.max_steps,
                "GpuFdtd: driven run exceeds max_steps ({} + {n} > {})",
                self.steps_taken,
                self.max_steps
            );
        }
        let extents = self.dispatch_extents();
        let groups = |len: usize| (len as u32).div_ceil(WORKGROUP);
        let lane_groups = |count: usize| (count as u32).div_ceil(64).max(1);
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
                    // H half-step (3 components).
                    for (pipeline, extent) in self.update_pipelines[..3].iter().zip(extents) {
                        pass.set_pipeline(pipeline);
                        pass.dispatch_workgroups(
                            groups(extent.0),
                            groups(extent.1),
                            groups(extent.2),
                        );
                    }
                    // Soft sources between the half-steps (reference order).
                    if n_soft > 0 {
                        pass.set_pipeline(&self.drive_pipelines[0]);
                        pass.dispatch_workgroups(lane_groups(n_soft), 1, 1);
                    }
                    // E half-step (3 components), CPML fused.
                    for (pipeline, extent) in self.update_pipelines[3..].iter().zip(&extents[3..6])
                    {
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
                    // Ports, probes, then the counter bump (reference order).
                    if n_ports > 0 {
                        pass.set_pipeline(&self.drive_pipelines[1]);
                        pass.dispatch_workgroups(lane_groups(n_ports), 1, 1);
                    }
                    if n_probes > 0 {
                        pass.set_pipeline(&self.drive_pipelines[2]);
                        pass.dispatch_workgroups(lane_groups(n_probes), 1, 1);
                    }
                    if self.dft_enabled {
                        let total: usize = [
                            len3(self.spec.ex_dims()),
                            len3(self.spec.ey_dims()),
                            len3(self.spec.ez_dims()),
                            len3(self.spec.hx_dims()),
                            len3(self.spec.hy_dims()),
                            len3(self.spec.hz_dims()),
                        ]
                        .iter()
                        .sum();
                        pass.set_pipeline(&self.dft_pipeline);
                        pass.dispatch_workgroups(lane_groups(total), 1, 1);
                    }
                    if has_drive {
                        pass.set_pipeline(&self.drive_pipelines[3]);
                        pass.dispatch_workgroups(1, 1, 1);
                    }
                }
            }
            self.queue.submit(Some(encoder.finish()));
        }
        self.steps_taken += n;
        Ok(())
    }

    /// Read back the recorded probe series (one `Vec` per [`Drive::probes`]
    /// entry, `steps_taken` samples each), widened to FP64.
    pub fn read_probes(&mut self) -> Result<Vec<Vec<f64>>, ComputeError> {
        let (n_soft, n_ports, n_probes) = self.drive_counts;
        if n_probes == 0 {
            return Ok(Vec::new());
        }
        let widened = self.read_f32_buffer(&self.drive_data)?;
        let probe_base = 1 + 3 * n_ports + self.max_steps * (n_soft + n_ports);
        let mut series = vec![Vec::with_capacity(self.steps_taken); n_probes];
        for step in 0..self.steps_taken {
            for (q, out) in series.iter_mut().enumerate() {
                out.push(widened[probe_base + step * n_probes + q]);
            }
        }
        Ok(series)
    }

    /// Read back the on-GPU NTFF DFT phasor pair `(Ê_re, Ê_im)` — the
    /// running full-field DFT at the `with_ntff_dft` probe frequency,
    /// accumulated as `Σ F·cos(ωt)` / `−Σ F·sin(ωt)` over completed steps
    /// (no Δt factor; the reference `NtffState` supplies it when the pair
    /// is fed back through two synthetic samples).
    ///
    /// # Panics
    ///
    /// Panics unless the stepper was built with [`GpuFdtd::with_ntff_dft`].
    pub fn read_dft_fields(&mut self) -> Result<(Fields, Fields), ComputeError> {
        assert!(self.dft_enabled, "read_dft_fields needs with_ntff_dft");
        let widened = self.read_f32_buffer(&self.psi_arena)?;
        let total: usize = [
            len3(self.spec.ex_dims()),
            len3(self.spec.ey_dims()),
            len3(self.spec.ez_dims()),
            len3(self.spec.hx_dims()),
            len3(self.spec.hy_dims()),
            len3(self.spec.hz_dims()),
        ]
        .iter()
        .sum();
        // DFT region sits at the arena tail (after CPML/dispersion blocks).
        let base = widened.len() - 2 * total;
        let re = self.split_fields(&widened[base..base + total]);
        let im = self.split_fields(&widened[base + total..]);
        Ok((re, im))
    }

    /// Split a packed ex|ey|ez|hx|hy|hz slice into a [`Fields`].
    fn split_fields(&self, packed: &[f64]) -> Fields {
        let lens = [
            len3(self.spec.ex_dims()),
            len3(self.spec.ey_dims()),
            len3(self.spec.ez_dims()),
            len3(self.spec.hx_dims()),
            len3(self.spec.hy_dims()),
            len3(self.spec.hz_dims()),
        ];
        let mut iter = packed.iter().copied();
        let mut take = |len: usize| iter.by_ref().take(len).collect::<Vec<f64>>();
        Fields {
            ex: take(lens[0]),
            ey: take(lens[1]),
            ez: take(lens[2]),
            hx: take(lens[3]),
            hy: take(lens[4]),
            hz: take(lens[5]),
        }
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
        self.read_f32_buffer(&self.field_arena)
    }

    /// Copy a whole f32 storage buffer back to the host, widened to FP64.
    fn read_f32_buffer(&self, source: &wgpu::Buffer) -> Result<Vec<f64>, ComputeError> {
        let size = source.size();
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
        encoder.copy_buffer_to_buffer(source, 0, &staging, 0, size);
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
