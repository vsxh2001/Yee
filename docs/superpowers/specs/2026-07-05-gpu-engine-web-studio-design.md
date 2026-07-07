# GPU/CPU compute engine + web studio — design (E.0 walking skeleton + track architecture)

**Date:** 2026-07-05
**ADR:** ADR-0175
**Plan:** `docs/superpowers/plans/2026-07-05-gpu-engine-web-studio.md`
**Roadmap:** `ENGINE-STUDIO-ROADMAP.md` (phases E.* and S.*)

## 1. Goal

Re-center Yee on two products:

- **Part 1 — the engine.** A fast Rust simulation engine that uses GPU *and* CPU. Portable GPU
  compute via wgpu/WGSL (Vulkan/Metal/DX12, later WebGPU), multi-threaded CPU via rayon, with the
  existing cudarc/cuSOLVER dense-LU lane kept as-is. FDTD stepping is the beachhead workload.
- **Part 2 — the studio.** A modern web-technology UI (Tauri 2 shell + React/TypeScript frontend;
  small bundles, Rust host process) that drives the engine in the background through a
  transport-agnostic job protocol, usable both in-process (desktop) and over WebSocket
  (`yee-server`, browser).

This spec fixes the architecture for both and fully specifies the **E.0 walking skeleton**
(the only part built in the first increment).

## 2. Current state (verified 2026-07-05)

- `yee-fdtd`: scalar FP64 triple loops over `ndarray::Array3<f64>`; no rayon/SIMD/GPU; the `cuda`
  feature is wired in Cargo.toml but referenced by zero source lines.
- `yee-cuda`: cuSOLVER Zgetrf/Zgetrs + NVRTC behind an associated-function `Backend` trait; no CPU
  impl; it is a cudarc-churn insulation seam, not a dispatch layer.
- No wgpu compute and no `.wgsl` files exist anywhere; `yee-gui/src/viewport.rs` is render-only.
- No service layer exists (no tokio/axum/websocket/tauri anywhere in the workspace);
  `yee cli run` is a Phase-0 stub; ADR-0089's `yee-server` was never built.
- The studio is `yee-studio-web` (Dioxus→WASM, GitHub Pages), EM-solver-free by design; its
  Verify stage is the ideal coupling-matrix response.

## 3. Target architecture

```
            ┌────────────────────────────── studio (Part 2) ──────────────────────────────┐
            │  Tauri 2 shell (Rust host)              browser deployment                  │
            │  React + TS + Vite frontend             (same frontend bundle)              │
            └───────────────┬─────────────────────────────────┬───────────────────────────┘
                            │ in-process calls + events       │ WebSocket JSON + binary frames
                    ┌───────▼──────────┐              ┌───────▼──────────┐
                    │   yee-engine     │◄─────────────│    yee-server    │  (axum)
                    │  job API (serde) │   wraps      └──────────────────┘
                    │ submit/progress/ │
                    │ cancel/results   │
                    └───────┬──────────┘
            ┌───────────────┼──────────────────────────────┐
     ┌──────▼─────┐  ┌──────▼──────┐                ┌──────▼─────┐
     │ yee-compute│  │ yee-fdtd    │ (reference     │ yee-mom /  │
     │ CPU: rayon │  │ scalar      │  kernels +     │ yee-fem /  │
     │ GPU: wgpu  │  │ CPML/NTFF/  │  physics not   │ yee-cuda   │ (cuSOLVER LU lane, unchanged)
     │ (FP32 WGSL)│  │ dispersive) │  yet ported)   └────────────┘
     └────────────┘  └─────────────┘
```

- **`yee-compute`** owns *execution*: field buffers, kernels, device selection, readback. It does
  **not** own physics policy; `yee-fdtd` remains the reference implementation and the home of the
  not-yet-ported features (CPML, dispersive, NTFF, sources, subgrid).
- **`yee-engine`** owns *orchestration*: a `Job` = scene + solver config + outputs requested;
  progress streaming (step counts, field probes, partial sweeps); cancellation; results (fields,
  S-params, patterns) as serde types. One protocol, two transports (in-process, WebSocket).
- **Studio** is a pure client of the `yee-engine` protocol. Desktop = Tauri commands/events;
  browser = the same JSON over WebSocket against `yee-server`.

### Precision policy

- CPU backend: FP64 (matches every existing gate).
- GPU backend: FP32 storage and arithmetic (industry-standard for FDTD; consumer-GPU FP64 is
  1/32–1/64 rate). Gates that need FP64 stay on CPU until phase E.3 defines a mixed/FP64-capable
  GPU policy (shader-f64 where supported, or compensated summation where it matters).

### Device selection

Runtime enumeration: try wgpu adapter (HighPerformance preference) → fall back to CPU. The
selection is a library call, never a compile-time choice; `gpu` is a default-on cargo feature so
`--no-default-features` yields the CPU-only crate (WASM-safe seam for later browser compute).

## 4. E.0 walking skeleton — precise contract

**In scope:** uniform lossless vacuum (`eps_r`/`mu_r` scalars, σ = 0), PEC outer box (the same
implicit boundary as `yee-fdtd`'s skipped outer tangential E faces), Gaussian-ball `E_z` initial
condition, N-step leapfrog, full six-field readback. **Out of scope:** CPML, per-cell materials,
sources, lossy CA/CB, NTFF, subgrid — these are E.1/E.2.

### 4.1 Types (`crates/yee-compute`)

```rust
/// Uniform-vacuum FDTD problem description (E.0 scope).
pub struct FdtdSpec { nx, ny, nz: usize, dx, dy, dz, dt: f64, eps_r, mu_r: f64 }

/// Six flat field buffers (row-major [i][j][k], ndarray-default order),
/// staggered shapes identical to yee-fdtd's YeeGrid.
pub struct Fields { ex, ey, ez, hx, hy, hz: Vec<f64> }

pub enum FdtdEngine { Cpu(CpuFdtd), Gpu(GpuFdtd) }  // #[cfg(feature="gpu")] on the Gpu arm
impl FdtdEngine {
    fn new_cpu(spec, fields) -> Self;
    fn new_gpu(spec, fields) -> Result<Self, ComputeError>;  // Err when no adapter
    fn step_n(&mut self, n: usize) -> Result<(), ComputeError>;
    fn read_fields(&mut self) -> Result<Fields, ComputeError>;
    fn backend_name(&self) -> &'static str;
}
```

Flat indexing matches `ndarray`'s default row-major layout exactly:
`idx(i, j, k) = (i * dim_j + j) * dim_k + k` per staggered component shape. This makes CPU↔GPU
parity testable index-for-index and `yee-fdtd` interop a `as_slice()` copy.

### 4.2 CPU backend (`CpuFdtd`)

- The six update kernels are line-for-line ports of `yee-fdtd::update::{update_h, update_e}`
  (uniform lossless arm only), on flat `Vec<f64>`, parallelized with rayon by slabbing the
  outermost `i` index (`par_chunks_mut` on the target component). Each cell's update is
  independent within a half-step, so parallelization must be **bit-exact** vs the scalar
  reference — same per-cell expression, same operation order.
- E updates skip outer tangential faces exactly as the reference does (PEC box).

### 4.3 GPU backend (`GpuFdtd`)

- wgpu 29 (workspace-pinned), `pollster` to block on adapter/device futures.
- Six FP32 storage buffers (one per component), one uniform buffer for dims + precomputed
  coefficients (`dt/(μ₀μ_r·d)` and `dt/(ε₀ε_r·d)` per axis, computed in f64 host-side then cast).
- One WGSL module, **six entry points** (`update_hx` … `update_ez`), each dispatched over its own
  staggered extent with in-shader bounds checks; workgroup size 4×4×4. E entry points implement
  the interior-only ranges (PEC faces untouched, matching the reference).
- A step = 6 dispatches (3 H then 3 E). WebGPU guarantees writes from one dispatch are visible to
  subsequent dispatches (separate usage scopes), so all six ride one compute pass. `step_n`
  submits in chunks (≤ 64 steps per command buffer) to bound encoder size and avoid device
  timeouts.
- Readback: copy to a mapped staging buffer, `device.poll(wait)`, widen f32 → f64 into `Fields`.
- Upload: `Fields` (f64) narrowed to f32 at construction.

### 4.4 Validation gates

- **`compute-001` (CPU vs reference, bit-exact)** — `tests/cpu_reference_parity.rs`. Build a
  `YeeGrid::vacuum(24, 20, 22, 1e-3)` (deliberately non-cubic dims to catch index swaps), inject
  a Gaussian ball into `ez` on both sides, run 25 steps of reference `update_h`/`update_e` vs
  `CpuFdtd::step_n(25)`, assert **max |Δ| == 0.0** on all six components. Also assert the field
  actually propagated (energy moved beyond the initial ball) so an all-zeros bug can't pass.
- **`compute-002` (GPU vs CPU, FP32 tolerance)** — `tests/gpu_cpu_parity.rs`. Same scenario,
  100 steps. Per component: L2 error < 1e-4 and L∞ < 1e-3, both normalized by the component's
  **field-family norm** (E family or H family, from the CPU FP64 result). Family normalization,
  not per-component: a pure-E_z pulse excites H_z only at second order (norm ~10³ below H_x), so
  a per-component relative test amplifies plain FP32 round-off into spurious failures — measured
  on llvmpipe, all components sit at ~4e-7 family-relative. E and H families stay separate
  because they differ by η₀ in units. If `FdtdEngine::new_gpu` reports no adapter, print a
  SKIPPED notice and return — green on adapterless CI, real on the GPU nightly runner (and on
  Mesa llvmpipe when present).
- Both gates are fast (< a few seconds) and run non-`#[ignore]`'d in the default workspace test.

### 4.5 CI

- Default `ci.yml` picks the crate up automatically (workspace test); `compute-002` self-skips.
- `gpu-nightly.yml` gains a step: `cargo test -p yee-compute --release` on the self-hosted GPU
  runner, where `compute-002` exercises a real adapter.

## 5. Later phases (specified here, built later — see ENGINE-STUDIO-ROADMAP.md)

- **E.1** CPML + per-cell ε/σ/PEC masks on both backends (gate: existing `cpml_reflection`
  scenario reproduced on `yee-compute` within tolerance).
- **E.2** Sources/lumped ports + a `FdtdDriver`-equivalent so `fdtd-line-eeff-001` can run on the
  engine (gate: ε_eff parity with the `yee-fdtd` driver).
- **E.3** Precision policy phase: FP64-on-GPU where `SHADER_F64` exists; error-budget doc.
- **E.4** Performance: `yee-bench` gains `compute_step` benches; target ≥ 20× scalar-CPU
  throughput on a mid-range discrete GPU at 128³, ≥ (cores·0.6)× for rayon CPU.
- **S.0** `yee-engine` job API (serde protocol, in-process executor, progress/cancel).
- **S.1** `yee-server` (axum WS) exposing S.0; CLI `yee serve`.
- **S.2** Tauri 2 + React studio shell speaking the S.0 protocol in-process (own ADR at start).
- **S.3** Field/3D visualization (three.js viewport, plotly S-params) fed by engine streams.
- **S.4** Filter-studio parity audit vs `yee-studio-web`; Dioxus retirement decision (own ADR).

## 6. Risks

- **wgpu API churn** — mitigated: workspace already pins wgpu 29 (egui-wgpu 0.34 hard-requires
  it); `yee-compute` uses the same pin, bumping remains a workspace-wide decision.
- **FP32 accumulation drift over long runs** — bounded by gate tolerances now; E.3 owns the
  long-run policy. Skeleton runs are O(100) steps.
- **No GPU in default CI** — parity gate self-skips; the GPU nightly runner is the real gate,
  same posture as the CUDA lane.
- **Two studios during transition** — Dioxus studio is feature-frozen but deployed until S.4
  parity; the freeze is explicit in ADR-0175.
