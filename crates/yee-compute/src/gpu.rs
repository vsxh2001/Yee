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
//! A full step is two fused bulk+CPML dispatches (`update_h` then
//! `update_e`, each computing all three staggered components per cell —
//! FS.7.1, ADR-0224) plus three mask clamps when masks are attached, all in
//! one compute pass: WebGPU orders storage-buffer writes between dispatches,
//! so no explicit barriers are needed. `step_n` submits in chunks to bound
//! command-buffer size.
//! The FP64 host state is narrowed on upload and widened on readback; gates
//! `compute-002`/`compute-005` bound the accumulated FP32 error against the
//! CPU backend. The PEC box is enforced host-side by zeroing the outer
//! tangential E faces at upload — no kernel ever writes those faces, so the
//! clamp holds for the whole run.

use std::sync::mpsc;
use std::sync::{Mutex, MutexGuard};

use yee_core::units::{EPS0, MU0};

use crate::cpml::make_profiles;
use crate::cpu::apply_pec_box;
use crate::dispersive::{DispersiveMap, ade_coeffs};
use crate::drive::{Drive, EComponent};
use crate::error::ComputeError;
use crate::fields::Fields;
use crate::materials::{Boundary, Materials};
use crate::spec::{FdtdSpec, GradedSpacings, SpacingArrays, len3};

/// Steps encoded per queue submission. Keeps a single command buffer to a
/// few hundred dispatches so watchdog-limited platforms don't kill long runs.
//
// FS.7.0 task 3 (ADR-0223) root-caused the 128^3->224^3 throughput decline
// (2408 -> 1841 -> 1401 -> 934 Mcells/s, i.e. -23.5%/-23.9%/-33.3% per +32
// grid step; RTX 5060 Ti) against this const and found it is NOT the cause:
//   (a) chunking: swept STEPS_PER_SUBMIT in {8,16,32,64,128,256} at every
//       grid size in the bench sweep. 128^3/192^3 Mcells/s stayed flat to
//       within run-to-run noise (~0.3%) across the whole 32x range (e.g.
//       192^3: 1400.9/1401.2/1401.1/1401.0/1400.5/1399.9) -- submission /
//       encoder overhead is not the bottleneck (a single compute pass
//       already covers a whole chunk, so this also rules out per-pass
//       bind-group overhead, hypothesis (b)).
//   (d) power/thermal: `nvidia-smi dmon` during a run showed clocks
//       reaching ~2600-2700 MHz (near the ~3090 MHz boost cap), temps
//       34-45 C, power draw peaking ~107 W -- no throttling.
//   (c) the >128 MiB single-binding path (field arena crosses 128 MiB
//       between 160^3=93.8 MiB and 192^3=162.0 MiB; coeff arena crosses it
//       between 192^3=109.7 MiB and 224^3=173.8 MiB, see `78cd12f`): the
//       decline is smooth across both crossings, not a step -- the 128->160
//       Mcells/s ratio (-23.5%, no crossing) and the 160->192 ratio (-23.9%,
//       crosses the field-arena threshold) are statistically the same drop.
//       No extra driver-side penalty was measured at the crossing.
// Net: the decline is consistent with the per-step working set (fields +
// coeffs, ~10 MiB at 64^3 growing to ~431 MiB at 224^3) outstripping
// on-chip cache reuse as the grid grows -- a memory-hierarchy/roofline
// effect, not a bug in this crate, and not fixable by a chunk-size change.
// STEPS_PER_SUBMIT stays 64 (no measured value beats it); see ADR-0223 for
// the full sweep and analysis.
const STEPS_PER_SUBMIT: usize = 64;

// Workgroup shape for the volume kernels in `shaders/fdtd.wgsl` (post-FS.7.1
// fusion: `update_h`, `update_e`, `clamp_ex/ey/ez` — 5, was 9 pre-fusion;
// `@workgroup_size(WORKGROUP_X, WORKGROUP_Y, WORKGROUP_Z)` there — keep the
// declarations in lockstep).
//
// FS.7.0 (ADR-0223) measured (4,4,4) as the best shape UNDER THE OLD
// gid.x->i mapping (flat-x lost 4-5x there because gid.x drove the
// slowest-varying array axis). FS.7.1 (ADR-0224) remapped gid<->linearization
// (gid.x -> k, gid.y -> j, gid.z -> i, adjacent threads -> adjacent
// k-fastest memory) and re-measured from scratch, RTX 5060 Ti, 128^3/192^3
// Mcells/s (median of 3x200-step reps, post-remap):
//   (4,4,4)   [control]: 3166 / 2998  (was the ADR-0223 winner)
//   (64,1,1)            : 3461 / 3143  (+9.3% / +4.8%)
//   (32,2,2)  [this]    : 3449 / 3144  (+8.9% / +4.9%)
//   (16,4,4)            : 3470 / 3152  (+9.6% / +5.1%)
// All three flat-ish post-remap shapes are within noise of each other at
// 128^3/192^3/224^3 (now that gid.x tracks the fast axis, the previous
// 4-5x flat-x penalty is gone), but they diverge sharply at the small
// grids where the per-dispatch overhead matters more: (32,2,2) clearly
// wins 64^3/96^3 (7569/9379 Mcells/s, reproduced) vs (64,1,1) (6943/8642)
// and (16,4,4) (6723/8098). (32,2,2) is the only shape that is at or near
// the top across all 6 grid sizes, so it is kept as the hardcoded shape
// (no runtime knob per YAGNI). See ADR-0224 for the full 6-grid before/
// after table.
const WORKGROUP_X: u32 = 32;
const WORKGROUP_Y: u32 = 2;
const WORKGROUP_Z: u32 = 2;

/// Uniform block mirrored by `struct Params` in `shaders/fdtd.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    nx: u32,
    ny: u32,
    nz: u32,
    npml: u32,
    /// Per-face CPML enable (R.3): bit `2·axis + side`, side 0 = min face,
    /// side 1 = max face.
    faces_mask: u32,
    has_cpml: u32,
    has_mask: u32,
    has_dispersion: u32,
    /// ω·Δt for the on-GPU NTFF DFT accumulator (E.5b); 0.0 disables it.
    dft_omega_dt: f32,
}

/// Serializes wgpu `Instance`/`Device` creation against teardown, process-wide.
///
/// Root cause (teardown-report.md): two `GpuFdtd`s built on independent
/// `wgpu::Instance`s whose Vulkan devices get destroyed on different threads
/// at the same time corrupt process-global state inside the NVIDIA Vulkan
/// ICD (`libnvidia-glcore.so`) — observed as a SIGSEGV/SIGABRT
/// (`free(): invalid size` / `double free or corruption`) inside
/// `vkDestroyDevice`, ~80% of multi-threaded `cargo test -p yee-engine
/// --release` runs. This is a driver bug, not a wgpu or yee-compute logic
/// error, so the fix is a documented mitigation: never let two
/// create-or-destroy sequences run concurrently. Only the brief
/// creation/teardown window is locked — `step_n` and friends run fully
/// concurrent across `GpuFdtd`s, so this doesn't serialize actual compute.
static GPU_LIFECYCLE_LOCK: Mutex<()> = Mutex::new(());

fn gpu_lifecycle_guard() -> MutexGuard<'static, ()> {
    GPU_LIFECYCLE_LOCK
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// The actual GPU resource bundle — split out from [`GpuFdtd`] so `Drop` can
/// take the [`gpu_lifecycle_guard`] before tearing any of it down; see
/// `GPU_LIFECYCLE_LOCK`.
#[derive(Debug)]
pub struct GpuResources {
    spec: FdtdSpec,
    device: wgpu::Device,
    queue: wgpu::Queue,
    bind_group: wgpu::BindGroup,
    /// Update pipelines in dispatch order: H (fused hx/hy/hz), E (fused
    /// ex/ey/ez) — FS.7.1, ADR-0224.
    update_pipelines: [wgpu::ComputePipeline; 2],
    /// Mask clamp pipelines (ex, ey, ez); dispatched only when `has_mask`.
    clamp_pipelines: [wgpu::ComputePipeline; 3],
    /// Drive pipelines: inject_soft, apply_ports, record_probes,
    /// bump_step, apply_aperture_ports, record_h_probes (FS.4.2a).
    drive_pipelines: [wgpu::ComputePipeline; 6],
    has_mask: bool,
    dft_enabled: bool,
    /// accumulate_dft pipeline (E.5b).
    dft_pipeline: wgpu::ComputePipeline,
    field_arena: wgpu::Buffer,
    psi_arena: wgpu::Buffer,
    drive_data: wgpu::Buffer,
    /// Packed inverse primal/dual spacings (FS.0b.2) — uniform fill at
    /// build, refreshed in place by [`GpuFdtd::set_spacings`].
    spacing_buffer: wgpu::Buffer,
    /// `(npml, faces)` when CPML is active, for the graded scope check.
    cpml_layers: Option<(usize, [[bool; 2]; 3])>,
    has_dispersion: bool,
    /// Host copy of the drive: `set_spacings` recomputes the spacing-derived
    /// port constants (resistive α/γ, aperture `vcoef`) from it.
    drive: Drive,
    // (n_soft, n_ports, n_probes, n_aperture, n_hprobes)
    drive_counts: (usize, usize, usize, usize, usize),
    max_steps: usize,
    steps_taken: usize,
    adapter_name: String,
}

/// GPU FDTD stepper (per-cell materials, CPML or PEC box, interior masks).
///
/// A thin wrapper around `Option<GpuResources>`: the `Option` lets `Drop`
/// pull the resources out and drop them explicitly while holding
/// `GPU_LIFECYCLE_LOCK`, which a plain field-by-field auto-drop can't do
/// (the lock would have to be a struct field, held for the stepper's whole
/// life, not just teardown). All other methods reach fields/methods on
/// `GpuResources` through `Deref`/`DerefMut`.
#[derive(Debug)]
pub struct GpuFdtd(Option<GpuResources>);

impl std::ops::Deref for GpuFdtd {
    type Target = GpuResources;
    fn deref(&self) -> &GpuResources {
        self.0.as_ref().expect("GpuFdtd used after being dropped")
    }
}

impl std::ops::DerefMut for GpuFdtd {
    fn deref_mut(&mut self) -> &mut GpuResources {
        self.0.as_mut().expect("GpuFdtd used after being dropped")
    }
}

impl Drop for GpuFdtd {
    fn drop(&mut self) {
        let _guard = gpu_lifecycle_guard();
        self.0.take(); // vkDestroyDevice etc. run here, still under the lock
    }
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
        if materials.sheet_r_ohm.is_some_and(|r| r > 0.0) {
            return Err(ComputeError::Unsupported(
                "resistive-sheet conductor loss (R.0b) is not on the GPU yet",
            ));
        }
        if drive.aperture_ports.iter().any(|p| p.record) {
            return Err(ComputeError::Unsupported(
                "aperture-port (v, i) recording (FS.2a) is not on the GPU yet",
            ));
        }
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

        // Instance/adapter/device stand-up races against any other
        // `GpuFdtd`'s teardown in the NVIDIA Vulkan ICD — see
        // `GPU_LIFECYCLE_LOCK`. Held only across this bring-up block, not
        // the whole stepper lifetime.
        let gpu_lifecycle_guard = gpu_lifecycle_guard();
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .map_err(|_| ComputeError::NoAdapter)?;
        let adapter_name = adapter.get_info().name.clone();

        // Lift only the buffer-size caps to what the adapter actually
        // supports: wgpu's *default* limits cap storage-buffer bindings at
        // 128 MiB, which rejects grids whose fields arena exceeds that
        // (≈192³ and up) on hardware that could hold them easily. All other
        // limits stay at the WebGPU defaults so the browser-compute seam is
        // unaffected (a browser adapter reports its own, already-clamped
        // limits).
        let adapter_limits = adapter.limits();
        let required_limits = wgpu::Limits {
            max_storage_buffer_binding_size: adapter_limits.max_storage_buffer_binding_size,
            max_buffer_size: adapter_limits.max_buffer_size,
            ..wgpu::Limits::default()
        };
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            required_limits,
            ..Default::default()
        }))
        .map_err(|e| ComputeError::Device(e.to_string()))?;
        drop(gpu_lifecycle_guard);

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("yee-compute fdtd"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/fdtd.wgsl").into()),
        });

        // Bindings: 0 = uniform Params; 1 = fields (rw); 2 = coeffs (r);
        // 3 = psi (rw); 4 = profiles (r); 5 = masks (r); 6 = drive index (r);
        // 7 = drive data (rw); 8 = inverse spacings (r; FS.0b.2 — the 8th
        // storage buffer, exactly the WebGPU default per-stage limit).
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
                storage(8, true),
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
        let update_pipelines = [make_pipeline("update_h"), make_pipeline("update_e")];
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
            make_pipeline("apply_aperture_ports"),
            make_pipeline("record_h_probes"),
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
        let (mut psi_data, profile_data, npml, faces_mask) = match cpml_config {
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
                // Per-face enable bits (R.3): bit 2·axis for the min face,
                // bit 2·axis + 1 for the max face.
                let faces_mask = config.faces.iter().enumerate().fold(0u32, |m, (a, sides)| {
                    let lo = if sides[0] { 1u32 << (2 * a) } else { 0 };
                    let hi = if sides[1] { 1u32 << (2 * a + 1) } else { 0 };
                    m | lo | hi
                });
                (
                    vec![0.0f32; psi_len],
                    profiles,
                    config.npml as u32,
                    faces_mask,
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
        let cpml_layers = cpml_config.map(|c| (c.npml, c.faces));

        // --- inverse-spacing buffer (FS.0b.2): uniform fill at build; the
        // entries are bit-equal to the retired scalar `Params.inv_*` values,
        // so the uniform path is unchanged bit-for-bit (compute-020).
        // `set_spacings` refreshes the contents in place (lengths are
        // functions of the spec, so the bind group never changes). ---
        let spacing_buffer = make_storage_buffer(
            "yee-compute spacings",
            &SpacingArrays::uniform(&spec).inverse_f32(),
            false,
        );

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
        let (n_soft, n_ports, n_probes, n_aperture, n_hprobes) = (
            drive.soft_sources.len(),
            drive.ports.len(),
            drive.probes.len(),
            drive.aperture_ports.len(),
            drive.h_probes.len(),
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
        // Aperture-port cell table (R.3), append-only so the accessors above
        // keep their offsets: [n_ap, n_cells ×n_ap, cells_start ×n_ap,
        // flat E_z field-offsets ...] (starts relative to the cells base).
        drv_idx.push(n_aperture as u32);
        for port in &drive.aperture_ports {
            drv_idx.push(port.cells.len() as u32);
        }
        let mut start = 0u32;
        for port in &drive.aperture_ports {
            drv_idx.push(start);
            start += port.cells.len() as u32;
        }
        for port in &drive.aperture_ports {
            for &cell in &port.cells {
                drv_idx.push(
                    (EComponent::Ez.arena_offset(&spec) + EComponent::Ez.flat(&spec, cell)) as u32,
                );
            }
        }
        // H-probe field-arena offsets (FS.4.2a), appended after the
        // aperture-port block so every accessor above keeps its offset (same
        // append-only convention as the aperture block itself): [n_hp,
        // h-probe field-arena offsets ×n_hp]. `HComponent::arena_offset`
        // already includes the full E-block prefix, so this indexes the same
        // packed `fields` arena the E probes do.
        drv_idx.push(n_hprobes as u32);
        for probe in &drive.h_probes {
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
        //            probes (max_steps × n_probes),
        //            then the R.3 aperture block: v_prev ×n_ap, vcoef ×n_ap,
        //            g ×n_ap, back ×n_ap, ap_vsrc (max_steps × n_ap),
        //            then the FS.4.2a H-probe block: h-probes
        //            (max_steps × n_hp)] — all f64-precomputed.
        let mut drv_data: Vec<f32> = Vec::with_capacity(
            1 + 3 * n_ports
                + max_steps * (n_soft + n_ports + n_probes + n_aperture + n_hprobes)
                + 4 * n_aperture,
        );
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
        // Aperture-port constants (R.3), mirroring the CPU update exactly:
        // vcoef = dz/n_col (modal V from the cell sum), g = 1/(R + β) with
        // β = dt·h/(2ε₀A) (0 for an open port — the CPU's infinite-R arm),
        // back = dt/(ε₀A) (sheet back-action per unit branch current).
        drv_data.extend(std::iter::repeat_n(0.0, n_aperture)); // v_prev
        for port in &drive.aperture_ports {
            drv_data.push((spec.dz / port.n_columns as f64) as f32);
        }
        for port in &drive.aperture_ports {
            let beta = spec.dt * port.height / (2.0 * EPS0 * port.area);
            let g = if port.resistance.is_finite() {
                1.0 / (port.resistance + beta)
            } else {
                0.0
            };
            drv_data.push(g as f32);
        }
        for port in &drive.aperture_ports {
            drv_data.push((spec.dt / (EPS0 * port.area)) as f32);
        }
        for n in 0..max_steps {
            for port in &drive.aperture_ports {
                drv_data.push(port.waveform.value(n, spec.dt) as f32);
            }
        }
        // H-probe output region (FS.4.2a): one sample per H-probe per step,
        // written by `record_h_probes` (same "step-indexed slot" layout as
        // the E-probe region above).
        drv_data.extend(std::iter::repeat_n(0.0, max_steps * n_hprobes));
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
            faces_mask,
            has_cpml: u32::from(cpml_config.is_some()),
            has_mask: u32::from(has_mask),
            has_dispersion: u32::from(dispersive.is_some()),
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
                wgpu::BindGroupEntry {
                    binding: 8,
                    resource: spacing_buffer.as_entire_binding(),
                },
            ],
        });

        Ok(Self(Some(GpuResources {
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
            spacing_buffer,
            cpml_layers,
            has_dispersion: dispersive.is_some(),
            drive,
            drive_counts: (n_soft, n_ports, n_probes, n_aperture, n_hprobes),
            max_steps,
            steps_taken: 0,
            adapter_name,
        })))
    }

    /// The problem description this stepper was built from.
    pub fn spec(&self) -> &FdtdSpec {
        &self.spec
    }

    /// Adapter name (diagnostics; e.g. printed by the parity gates).
    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    /// Attach per-axis nonuniform primal spacings (FS.0b.2, ADR-0214),
    /// mirroring [`crate::CpuFdtd::set_spacings`]: the H kernels use the
    /// primal cell width at the H sample, the E kernels the dual spacing.
    /// Refreshes the inverse-spacing buffer and the spacing-derived port
    /// constants (resistive α/γ from the local dual-transverse area ×
    /// primal dz, aperture `vcoef` from the local primal dz) in place.
    /// Constant arrays equal to the scalar `spec.dx/dy/dz` are bit-exact
    /// against the GPU's own scalar path (gate `compute-020`). Call before
    /// stepping.
    ///
    /// # Errors
    ///
    /// Returns [`ComputeError::Unsupported`] for the GPU-capability gaps:
    /// the on-GPU NTFF DFT accumulator (`with_ntff_dft` — downstream NTFF
    /// surface integration assumes a uniform grid), or an aperture port
    /// whose cells straddle non-uniform z-spacing (the per-port modal
    /// coefficient is one scalar; FS.0b.1 `auto_spacings` substrates are
    /// uniform in z, so engine apertures are unaffected).
    ///
    /// # Panics
    ///
    /// Panics on invalid spacings (wrong lengths, non-positive widths),
    /// grading inside an absorbing CPML layer, a dispersive map (uniform
    /// only, FS.0b.0), `dt` above the graded Courant limit, or when called
    /// after stepping — the same contract as the CPU backend. Untrusted
    /// callers should pre-flight with [`GradedSpacings::validate`] /
    /// [`GradedSpacings::validate_cpml_layers`].
    pub fn set_spacings(&mut self, graded: &GradedSpacings) -> Result<(), ComputeError> {
        if let Err(e) = graded.validate(&self.spec) {
            panic!("set_spacings: {e}");
        }
        if let Some((npml, faces)) = self.cpml_layers
            && let Err(e) = graded.validate_cpml_layers(npml, faces)
        {
            panic!("set_spacings: {e}");
        }
        assert!(
            !self.has_dispersion,
            "set_spacings: dispersive ADE materials are uniform-grid only (FS.0b.0)"
        );
        assert!(
            self.spec.dt <= graded.courant_limit(),
            "set_spacings: dt = {} s exceeds the graded Courant limit {} s \
             (use the minimum spacing per axis)",
            self.spec.dt,
            graded.courant_limit()
        );
        assert!(
            self.steps_taken == 0,
            "set_spacings: call before stepping (the GPU state has already \
             advanced {} steps on the previous spacings)",
            self.steps_taken
        );
        if self.dft_enabled {
            return Err(ComputeError::Unsupported(
                "the on-GPU NTFF DFT accumulator on a graded grid (NTFF surface \
                 integration assumes a uniform grid; FS.0b.2 keeps NTFF+graded CPU-rejected)",
            ));
        }
        let sp = SpacingArrays::graded(graded);

        // Aperture `vcoef = dz/n_columns` is one scalar per port: exact only
        // when every cell of the port shares one primal dz. Reject
        // z-taper-straddling ports BEFORE any buffer is touched, so a failed
        // call leaves the stepper on its previous (uniform) spacings.
        let mut vcoefs = Vec::with_capacity(self.drive.aperture_ports.len());
        for port in &self.drive.aperture_ports {
            let dz0 = sp.z.primal[port.cells[0].2];
            if port.cells.iter().any(|c| sp.z.primal[c.2] != dz0) {
                return Err(ComputeError::Unsupported(
                    "an aperture port spanning non-uniform z-spacing on the GPU \
                     (the per-port modal coefficient is one scalar); keep the \
                     aperture inside a uniform-dz region",
                ));
            }
            vcoefs.push((dz0 / port.n_columns as f64) as f32);
        }

        self.queue.write_buffer(
            &self.spacing_buffer,
            0,
            bytemuck::cast_slice(&sp.inverse_f32()),
        );

        // Resistive-port α/γ with the local cell sizes at the port cell
        // (dual transverse spacings × primal dz) — the exact CPU formulas,
        // in f64, at the α|γ region of `drv_data` (indices [1+n, 1+3n)).
        let n_ports = self.drive.ports.len();
        if n_ports > 0 {
            let mut alpha_gamma = Vec::with_capacity(2 * n_ports);
            for port in &self.drive.ports {
                let (ci, cj, ck) = port.cell;
                let area = sp.x.dual[ci] * sp.y.dual[cj];
                let dz_c = sp.z.primal[ck];
                alpha_gamma
                    .push((self.spec.dt * dz_c / (2.0 * EPS0 * port.resistance * area)) as f32);
            }
            for port in &self.drive.ports {
                let (ci, cj, _) = port.cell;
                let area = sp.x.dual[ci] * sp.y.dual[cj];
                alpha_gamma.push((self.spec.dt / (EPS0 * port.resistance * area)) as f32);
            }
            self.queue.write_buffer(
                &self.drive_data,
                ((1 + n_ports) * std::mem::size_of::<f32>()) as u64,
                bytemuck::cast_slice(&alpha_gamma),
            );
        }

        // Aperture `vcoef` region: dd_ap_base + n_aperture (see the drv_data
        // layout comment in `build` and the dd_ap_* helpers in fdtd.wgsl).
        if !vcoefs.is_empty() {
            let (n_soft, n_ports, n_probes, n_aperture, _n_hprobes) = self.drive_counts;
            let base =
                1 + 3 * n_ports + self.max_steps * (n_soft + n_ports + n_probes) + n_aperture;
            self.queue.write_buffer(
                &self.drive_data,
                (base * std::mem::size_of::<f32>()) as u64,
                bytemuck::cast_slice(&vcoefs),
            );
        }
        Ok(())
    }

    /// Dispatch extent shared by both fused update pipelines (`update_h`,
    /// `update_e`): the cell-centered union `[nx+1, ny+1, nz+1]` of the
    /// three staggered per-component shapes each fuses over (FS.7.1,
    /// ADR-0224). The tuple is `(dim_i, dim_j, dim_k)` — the shader's own
    /// axis order, not the gid order (`gid.x` covers `k`, `gid.z` covers
    /// `i`, see the `dispatch_workgroups` calls below).
    fn fused_extent(&self) -> (usize, usize, usize) {
        self.spec.cell_dims()
    }

    /// Dispatch extents per clamp pipeline (ex, ey, ez) — clamps stay
    /// per-component (each touches only its own field + mask arena, no
    /// fusion benefit). Same `(dim_i, dim_j, dim_k)` axis-order note as
    /// [`GpuFdtd::fused_extent`].
    fn clamp_extents(&self) -> [(usize, usize, usize); 3] {
        [
            self.spec.ex_dims(),
            self.spec.ey_dims(),
            self.spec.ez_dims(),
        ]
    }

    /// Advance the state by `n` leapfrog steps (3 H, 3 E, then the mask
    /// clamps when masks are attached), submitting in [`STEPS_PER_SUBMIT`]
    /// chunks.
    pub fn step_n(&mut self, n: usize) -> Result<(), ComputeError> {
        let (n_soft, n_ports, n_probes, n_aperture, n_hprobes) = self.drive_counts;
        let has_drive =
            n_soft + n_ports + n_probes + n_aperture + n_hprobes > 0 || self.dft_enabled;
        if has_drive {
            assert!(
                self.steps_taken + n <= self.max_steps,
                "GpuFdtd: driven run exceeds max_steps ({} + {n} > {})",
                self.steps_taken,
                self.max_steps
            );
        }
        let fused_extent = self.fused_extent();
        let clamp_extents = self.clamp_extents();
        let groups_x = |len: usize| (len as u32).div_ceil(WORKGROUP_X);
        let groups_y = |len: usize| (len as u32).div_ceil(WORKGROUP_Y);
        let groups_z = |len: usize| (len as u32).div_ceil(WORKGROUP_Z);
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
                    // H half-step (hx/hy/hz fused into one dispatch, FS.7.1).
                    pass.set_pipeline(&self.update_pipelines[0]);
                    pass.dispatch_workgroups(
                        groups_x(fused_extent.2),
                        groups_y(fused_extent.1),
                        groups_z(fused_extent.0),
                    );
                    // Soft sources between the half-steps (reference order).
                    if n_soft > 0 {
                        pass.set_pipeline(&self.drive_pipelines[0]);
                        pass.dispatch_workgroups(lane_groups(n_soft), 1, 1);
                    }
                    // E half-step (ex/ey/ez fused into one dispatch, CPML
                    // fused in too, FS.7.1).
                    pass.set_pipeline(&self.update_pipelines[1]);
                    pass.dispatch_workgroups(
                        groups_x(fused_extent.2),
                        groups_y(fused_extent.1),
                        groups_z(fused_extent.0),
                    );
                    if self.has_mask {
                        for (pipeline, extent) in self.clamp_pipelines.iter().zip(&clamp_extents) {
                            pass.set_pipeline(pipeline);
                            pass.dispatch_workgroups(
                                groups_x(extent.2),
                                groups_y(extent.1),
                                groups_z(extent.0),
                            );
                        }
                    }
                    // Ports, probes, then the counter bump (reference order).
                    if n_ports > 0 {
                        pass.set_pipeline(&self.drive_pipelines[1]);
                        pass.dispatch_workgroups(lane_groups(n_ports), 1, 1);
                    }
                    // Aperture ports after the single-cell ports (R.3),
                    // matching the CPU/reference ordering.
                    if n_aperture > 0 {
                        pass.set_pipeline(&self.drive_pipelines[4]);
                        pass.dispatch_workgroups(lane_groups(n_aperture), 1, 1);
                    }
                    if n_probes > 0 {
                        pass.set_pipeline(&self.drive_pipelines[2]);
                        pass.dispatch_workgroups(lane_groups(n_probes), 1, 1);
                    }
                    // H probes (FS.4.2a), recorded alongside the E probes —
                    // same reference order, same step-indexed slot scheme.
                    if n_hprobes > 0 {
                        pass.set_pipeline(&self.drive_pipelines[5]);
                        pass.dispatch_workgroups(lane_groups(n_hprobes), 1, 1);
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
        let (n_soft, n_ports, n_probes, _n_aperture, _n_hprobes) = self.drive_counts;
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

    /// Read back the recorded H-probe series (FS.4.2a; one `Vec` per
    /// [`Drive::h_probes`] entry, `steps_taken` samples each), widened to
    /// FP64. Mirrors [`GpuFdtd::read_probes`], offset past the
    /// aperture-port block (`dd_ap_base` + 4·n_ap + max_steps·n_ap in
    /// `fdtd.wgsl`'s `dd_hp_base`).
    pub fn read_h_probes(&mut self) -> Result<Vec<Vec<f64>>, ComputeError> {
        let (n_soft, n_ports, n_probes, n_aperture, n_hprobes) = self.drive_counts;
        if n_hprobes == 0 {
            return Ok(Vec::new());
        }
        let widened = self.read_f32_buffer(&self.drive_data)?;
        let hp_base = 1
            + 3 * n_ports
            + self.max_steps * (n_soft + n_ports + n_probes)
            + 4 * n_aperture
            + self.max_steps * n_aperture;
        let mut series = vec![Vec::with_capacity(self.steps_taken); n_hprobes];
        for step in 0..self.steps_taken {
            for (q, out) in series.iter_mut().enumerate() {
                out.push(widened[hp_base + step * n_hprobes + q]);
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

    /// Block until the device has finished all work submitted so far.
    ///
    /// `step_n` only *submits* command buffers — wgpu queues are
    /// asynchronous, so it returns as soon as the driver has accepted the
    /// work, not when the GPU has finished executing it. Naively timing a
    /// bare `step_n` call therefore measures submission overhead, not
    /// device time (confirmed: it reads a multi-hundred-x "speedup" that
    /// does not scale with grid size). `sync()` is the benchmark/sequencing
    /// seam: submit an empty command buffer (so there is always something
    /// to wait on even with zero prior submissions this call) and poll the
    /// device to completion, reusing the same wait idiom `read_fields`'s
    /// readback already blocks on. `read_fields`/`read_probes` do not need
    /// an extra `sync()` call — their buffer-map already waits for the
    /// device to be idle.
    pub fn sync(&self) -> Result<(), ComputeError> {
        let encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("yee-compute sync"),
            });
        self.queue.submit(Some(encoder.finish()));
        self.wait_idle()
    }

    /// Poll the device until all submitted work has retired. Shared by
    /// [`GpuFdtd::sync`] and the buffer-readback path.
    fn wait_idle(&self) -> Result<(), ComputeError> {
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .map_err(|e| ComputeError::Readback(e.to_string()))?;
        Ok(())
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
        self.wait_idle()?;
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
