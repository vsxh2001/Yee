# Engine + Studio Roadmap (GPU/CPU compute engine → web studio)

Direction set by **ADR-0175** (2026-07-05). This is the third top-level roadmap, alongside
`ROADMAP.md` (core EM solvers, Phases 0–4) and `FILTER-DESIGN-ROADMAP.md` (filter application).
It tracks the two-part re-centering of the project:

- **Part 1 — the engine (E.\*):** a fast Rust simulation engine that uses GPU *and* CPU.
  Portable wgpu/WGSL compute + rayon CPU in the new `crates/yee-compute`, with `yee-fdtd`'s
  scalar kernels kept as the validated reference and `yee-cuda`'s cuSOLVER LU lane unchanged.
- **Part 2 — the studio (S.\*):** an engine-service protocol (`yee-engine` → `yee-server`) and a
  modern web-technology studio (Tauri 2 shell + React/TypeScript frontend) that drives the
  engine in the background, in-process on desktop and over WebSocket in the browser.

Spec: `docs/superpowers/specs/2026-07-05-gpu-engine-web-studio-design.md`
Plan: `docs/superpowers/plans/2026-07-05-gpu-engine-web-studio.md`

Conventions match the other roadmaps: every phase ships behind a machine-checkable validation
gate; walking-skeleton first; phases get ADRs when they make a decision worth recording.

---

## Part 1 — Engine track (E.*)

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| **E.0** | `yee-compute` walking skeleton: `FdtdSpec`/`Fields`/`FdtdEngine`, rayon FP64 `CpuFdtd`, wgpu/WGSL FP32 `GpuFdtd`, uniform lossless vacuum + PEC box | `compute-001` (CPU **bit-exact** vs `yee-fdtd` scalar reference, 25 steps, non-cubic grid); `compute-002` (GPU vs CPU, rel-L2 < 1e-4 / L∞ < 1e-3, 100 steps; self-skips without adapter, real on GPU nightly) | **SHIPPED** (ADR-0175, this branch) |
| **E.1** | CPML + per-cell ε_r/μ_r/σ + interior PEC masks + legacy PEC box on both backends; GPU arena-buffer layout (5 storage bindings — inside WebGPU browser limits) | `compute-003` (CPU **bit-exact** vs reference, heterogeneous + CPML + masks, both boundary modes); `compute-004` (CPML reflection: **69.3 dB** measured vs ≥ 30 dB target); `compute-005` (GPU vs CPU on the full E.1 scenario: ~2e-7 E / ~3e-6 H family-rel on llvmpipe; CPML holds 210× less ‖H‖ than PEC) | **SHIPPED** (ADR-0176) |
| **E.2** | Sources, lumped ports, driver-equivalent step loop; first real workload end-to-end | `fdtd-line-eeff-001` ε_eff via `yee-compute` within the existing ±15% HJ gate; CPU↔GPU ε_eff agreement < 0.5% | queued |
| **E.3** | Precision policy: FP64-on-GPU where `SHADER_F64` exists, error-budget doc, long-run drift bounds | drift gate over ≥ 10⁴ steps against CPU FP64 | queued |
| **E.4** | Performance: `yee-bench` `compute_step` benches, workgroup/occupancy tuning, CPU SIMD pass | ≥ 20× scalar-CPU throughput on mid-range dGPU at 128³; rayon CPU ≥ 0.6·cores× scaling | queued |
| **E.5** | Dispersive ADE + NTFF on the engine (full `yee-fdtd` feature parity) | existing dispersive/NTFF gates reproduced via `yee-compute` | queued |

Non-goals for E.*: replacing `yee-cuda`'s cuSOLVER LU lane (stays as-is); MoM/FEM assembly on
wgpu (revisit after E.4 with data).

## Part 2 — Studio track (S.*)

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| **S.0** | `yee-engine` crate: transport-agnostic job protocol (submit/progress/cancel/results, serde), threaded in-process executor | serde round-trip + end-to-end in-process FDTD job test | queued |
| **S.1** | `yee-server` (axum): WebSocket exposure of S.0; `yee serve` CLI | end-to-end WS job test in CI | queued |
| **S.2** | Tauri 2 + React/TS/Vite studio shell speaking the S.0 protocol in-process (kickoff ADR required) | app builds in CI; job submit→progress→result panel works against a stub job | queued |
| **S.3** | Visualization: 3D viewport (three.js) + S-param/Smith plots fed by engine streams | golden-image or DOM-level smoke gates | queued |
| **S.4** | Filter-studio parity audit vs `yee-studio-web`; Dioxus retirement decision (own ADR) | parity checklist green | queued |

Standing decision during S.*: **`yee-studio-web` (Dioxus) is feature-frozen but stays deployed**
until S.4 concludes (ADR-0175). `yee-gui` (egui EM-analysis shell) is unaffected by this track.

---

*Last updated: 2026-07-06 — E.1 shipped (ADR-0176): CPML + per-cell materials + PEC masks on
both backends; GPU moved to arena buffers (WebGPU-limit-safe); gates compute-003/004/005 green
(69.3 dB reflection reduction; GPU parity ~2e-7 on llvmpipe). E.0 (ADR-0175) shipped the day
before. Next: E.2 (sources/ports + `fdtd-line-eeff-001` on the engine) or S.0 (`yee-engine`
job API).*
