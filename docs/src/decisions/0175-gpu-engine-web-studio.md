# ADR-0175: GPU/CPU compute engine (`yee-compute`) + engine-service → web studio direction

**Status:** Accepted
**Date:** 2026-07-05
**Related:** ADR-0006 (cudarc pin + `yee-cuda::Backend` swap seam), ADR-0089 (filter-app architecture:
light flow in WASM, heavy EM on a native `yee-server` — the server half was never built), ADR-0110 /
ADR-0130 (Dioxus studio, eframe view retired), ADR-0108 (F1.1b.1 — first full-wave FDTD gate in the
filter loop), `ENGINE-STUDIO-ROADMAP.md` (the phase plan this ADR anchors).

---

## Context

The project is re-centering on two products:

1. **A fast Rust simulation engine** that exploits both GPU and CPU. Today the FDTD solver
   (`yee-fdtd`) is a deliberate walking skeleton: scalar FP64 triple loops over
   `ndarray::Array3<f64>`, single-threaded, no SIMD, no rayon, and a `cuda` feature flag that no
   source file actually uses. The only GPU code in the repo is `yee-cuda`'s cuSOLVER dense LU
   (the MoM matrix path) behind a narrow, associated-function `Backend` trait with no CPU impl —
   an insulation seam against cudarc churn, not a compute-dispatch abstraction. There is no wgpu
   compute usage and no `.wgsl` file anywhere; the only wgpu in-tree is the render-only viewport
   in `yee-gui`.
2. **A web-technology studio UI** that drives that engine in the background. The current studio
   (`yee-studio-web`, Dioxus → WASM, deployed to GitHub Pages) deliberately excludes every EM
   solver from its dep graph (ADR-0089), so its Verify stage shows the ideal coupling-matrix
   response, never a full-wave solve. The "native `yee-server` the web client calls" that
   ADR-0089 assumed was specified but never built: no crate in the workspace depends on tokio,
   axum, websockets, or any RPC layer. `yee-cli run` is still a Phase-0 stub.

The gap is therefore symmetric: the engine has no parallel/GPU execution path, and the studio has
no engine to call.

### Why wgpu compute (and not CUDA-first) for the engine

- **Portability is the product.** "Uses GPU and CPU" must mean *any* GPU a designer owns — Metal
  on macOS, Vulkan on Linux, DX12 on Windows, and eventually WebGPU in the browser. wgpu covers
  all of these from one WGSL kernel source. CUDA covers exactly one vendor and requires a
  toolchain install (which is why the `cuda` feature has been dormant since Phase 1.5).
- **FDTD is the right first workload.** The Yee update is embarrassingly parallel (every cell's
  H update reads only E, and vice versa), bandwidth-bound, and already the hot path of the only
  full-wave gate in the filter loop (`fdtd-line-eeff-001`). It is also the workload where FP32 is
  standard industry practice (commercial FDTD solvers run FP32 by default; consumer GPUs execute
  FP64 at 1/32–1/64 rate), so a portable FP32 GPU path is both fast and conventional.
- **wgpu is already in the workspace** (wgpu 29 via egui-wgpu), so no new major dependency tree.

### Why the cuSOLVER path stays

Dense complex LU (MoM) is a different animal: FP64-sensitive, compute-bound, and served by a
vendor library we already wrap. `yee-cuda` keeps that lane; `yee-compute` does not replace it.

### Why Tauri 2 + web frontend (and what happens to Dioxus)

The user-facing direction is a web UI with a small footprint — "React with Electron, or something
more modern to reduce bundle size". That is Tauri 2: system webview (single-digit-MB installers vs
Electron's ~100 MB+), and — decisive for us — **the host process is Rust**, so the engine links
in-process with zero serialization for field data, and the same engine-service protocol can be
re-exposed over WebSocket for a pure-browser deployment. The frontend is React + TypeScript +
Vite (mainstream ecosystem for the plotting/3D/component depth a studio needs: three.js, plotly,
etc.), which is a deliberate break from the all-Rust Dioxus lineage (ADR-0110/0130) **for the
studio app only**. The deployed Dioxus studio stays live and untouched until the new studio
reaches feature parity (see Consequences); this ADR supersedes ADR-0089's UI-architecture half
*prospectively*, not retroactively.

## Decision

1. **New crate `crates/yee-compute`** — the portable GPU/CPU execution layer for grid solvers.
   - A per-solver engine interface (first: FDTD stepping) with two backends behind one API:
     `CpuFdtd` (rayon-parallel, FP64) and `GpuFdtd` (wgpu compute, WGSL, FP32), selected at
     runtime by device enumeration with CPU fallback.
   - Flat, GPU-layout-identical field buffers (row-major, matching `ndarray`'s default order) so
     CPU↔GPU parity is testable index-for-index.
   - Feature `gpu` (wgpu) is **default-on** (no external toolchain needed — CLAUDE.md's
     default-off rule applies to toolchain-gated features only); `--no-default-features` builds
     the CPU-only crate (the future WASM/browser-compute seam).
   - `yee-cuda` remains the cuSOLVER LU lane and cudarc swap seam; `yee-compute` does not absorb
     it in the walking skeleton.
2. **Validation gates for the skeleton** (per §4 house rules — no solver feature without a gate):
   - `compute-001`: `CpuFdtd` matches `yee-fdtd`'s scalar `update_h`/`update_e` **bit-for-bit**
     (uniform lossless vacuum, PEC box, Gaussian-ball initial condition, N steps). The rayon
     parallelization is over independent cells, so exact agreement is required, not approximate.
   - `compute-002`: `GpuFdtd` (FP32) matches `CpuFdtd` (FP64) within FP32 accumulation tolerance
     on the same scenario. The test **self-skips** (with a printed notice) when no wgpu adapter
     exists, and runs for real on the GPU nightly runner.
3. **Engine-service + studio track** (specified now, built in later phases; see
   `ENGINE-STUDIO-ROADMAP.md`): `yee-engine` (transport-agnostic job API: submit scene+solver
   config, stream progress, cancel, fetch results — serde types shared by all transports), then
   `yee-server` (axum WebSocket/JSON exposure of `yee-engine`), then the Tauri 2 + React studio
   consuming the same protocol in-process. New UI stack decisions land as their own ADRs when
   those phases start; this ADR fixes the direction.

## Consequences

- FDTD gains its first parallel execution path; `yee-fdtd`'s scalar kernels become the *reference
  implementation* that `yee-compute` is gated against — they are not deleted or modified.
- FP32-on-GPU / FP64-on-CPU is now an explicit, tested policy rather than an aspiration. Gates
  that need FP64 keep running on the CPU path until an FP64-capable GPU policy phase (E.3).
- The Dioxus studio (`yee-studio-web`) is **feature-frozen but deployed** during the transition;
  no new stages land there once the Tauri studio starts (S.2). Removal requires a follow-up ADR
  after parity.
- The `cuda` feature on `yee-fdtd` stays dormant; if a CUDA FDTD path ever lands it goes through
  `yee-compute`'s backend seam, not ad-hoc.
- **Scope of the walking skeleton (E.0, this ADR's shipped increment):** uniform lossless vacuum
  + PEC box only. CPML, per-cell materials, dispersive ADE, sources/ports, NTFF are E.1/E.2
  phases — the skeleton proves the dispatch/parity/readback plumbing, nothing else.

## Outcome (E.0 — SHIPPED, this branch)

`crates/yee-compute` landed with `FdtdSpec`/`FdtdEngine` (runtime backend selection via
`FdtdEngine::new_cpu` / `FdtdEngine::new_gpu`), `CpuFdtd` (rayon over z-slabs, FP64), `GpuFdtd`
(wgpu 29 compute, six WGSL entry points — one per staggered field component — FP32 storage
buffers, chunked command submission, staging-buffer readback), and both gates:
`tests/cpu_reference_parity.rs` (`compute-001`, bit-exact vs `yee-fdtd`) and
`tests/gpu_cpu_parity.rs` (`compute-002`, self-skipping). GPU nightly runs `compute-002` on real
hardware. Spec: `docs/superpowers/specs/2026-07-05-gpu-engine-web-studio-design.md`; plan:
`docs/superpowers/plans/2026-07-05-gpu-engine-web-studio.md`; phase ledger:
`ENGINE-STUDIO-ROADMAP.md`.

## References

- Taflove & Hagness §3.6 (Yee update), `crates/yee-fdtd/src/update.rs` (reference kernels)
- wgpu 29 / WGSL compute; WebGPU implicit inter-dispatch synchronization guarantees
- ADR-0089 (the never-built `yee-server` this direction finally delivers)
