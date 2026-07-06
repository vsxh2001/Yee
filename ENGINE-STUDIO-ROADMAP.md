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
| **E.2** | Drive layer: `SoftSource`/`ResistivePort`/`Probe`/`Drive` on both backends (GPU: whole-run f64-precomputed tables + on-GPU step counter → zero per-step host round-trips) | `compute-007` driven step **bit-exact** vs reference; `compute-006` cavity TE₁₀₁ vs **analytic Pozar**: CPU −0.063 %, GPU −0.063 %, CPU↔GPU 0.0000 %; `compute-008` line-eeff on the engine vs **Hammerstad–Jensen**: 0.132 % (≤ 15 % gate), 88.6 s release | **SHIPPED** (ADR-0177) |
| **E.3** | Precision policy: FP32-GPU/FP64-CPU characterized (WGSL has no f64 — SHADER_F64 unreachable without SPIR-V passthrough; noted) | `compute-009` drift over 10⁴ energy-conserving steps: 3e-6…2e-5 family-rel (√N random-walk), 100× inside the 1e-3 gate | **SHIPPED** (ADR-0177) |
| **E.4** | Performance: `yee-bench` `compute_step` (scalar vs rayon CPU vs GPU) landed; container numbers recorded | 4-core container: rayon scales 2.2× internally but nets **0.78×** vs scalar (flat-buffer kernel ~2.8× slower single-thread — bounds-checked idx arithmetic). Row-sliced kernels landed (ADR-0179): single-thread −27 %, 4-thread ≈ scalar (bandwidth-bound container); bit-exact gates unchanged. Real-hardware numbers via the GPU nightly bench; the 20×-dGPU target remains to be certified there | **CLOSED** (ADR-0179; hardware numbers pending nightly) |
| **E.5a** | Far-field on the engine: engine steps, reference `NtffState` consumes fields via host adapter | `compute-010` vs **analytic sin θ**: broadside/endfire 327.9 dB (≥ 20 dB gate) | **SHIPPED** (ADR-0177) |
| **E.5b** | On-GPU full-field DFT phasor accumulation (`accumulate_dft` kernel, psi-arena tail, on-GPU step counter — zero per-step readback); reference `NtffState` projects via two synthetic samples | `compute-013`: GPU-resident dipole — **315.4 dB** analytic null, broadside matches the CPU path to **2.9e-7** | **SHIPPED** (ADR-0179) |
| **E.5c** | Dispersive ADE (Drude/Lorentz/Debye) on both backends: verbatim CPU port; unified-ADE GPU form folded into the coeff/psi arenas | `compute-011` **bit-exact** vs `yee_fdtd::dispersive` (four-arm scenario); `compute-012` differential GPU gate (ADE ≤ 20× standard-pair error, measured ≤ 6×; drift-class backstop) | **SHIPPED** (ADR-0179) |

Non-goals for E.*: replacing `yee-cuda`'s cuSOLVER LU lane (stays as-is); MoM/FEM assembly on
wgpu (revisit after E.4 with data).

## Part 2 — Studio track (S.*)

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| **S.0** | `yee-engine` crate: serde `JobSpec`/`JobEvent`/`JobResult` protocol + threaded chunked executor with progress streaming, cooperative cancel, cpu/gpu/auto backend selection | 4 unit tests + doctest: serde round-trip, progress stream, cancellation, auto-backend | **SHIPPED** (ADR-0179) |
| **S.1** | `yee-server` (axum): WebSocket exposure of S.0; `yee serve` CLI | end-to-end WS job test in CI | queued |
| **S.2** | Tauri 2 + React/TS/Vite studio shell (`studio/`, outside the root workspace) speaking S.0 in-process: `run_job` command + `job://progress` events + probe SVG plot. Frontend 47.9 kB gzipped | walking skeleton verified in-container: `cargo check` (webkit2gtk) + `npm run build` green; interactive run + CI wiring are the S.2 follow-on | **SKELETON SHIPPED** (ADR-0179) |
| **S.3** | Visualization: 3D viewport (three.js) + S-param/Smith plots fed by engine streams | golden-image or DOM-level smoke gates | queued |
| **S.4** | Filter-studio parity audit vs `yee-studio-web`; Dioxus retirement decision (own ADR) | parity checklist green | queued |

Standing decision during S.*: **`yee-studio-web` (Dioxus) is feature-frozen but stays deployed**
until S.4 concludes (ADR-0175). `yee-gui` (egui EM-analysis shell) is unaffected by this track.

---

*Last updated: 2026-07-06 (later) — engine track COMPLETE through E.5 (ADR-0179): E.4 closed
(row-sliced kernels), E.5b shipped (on-GPU NTFF accumulation, 315.4 dB / 2.9e-7 cross-backend),
E.5c shipped (dispersive ADE, bit-exact CPU + differential GPU gate). Python bindings
`yee.compute` shipped (ADR-0178). Studio track underway: S.0 `yee-engine` job API SHIPPED,
S.2 Tauri 2 + React skeleton SHIPPED (47.9 kB gzipped frontend; cargo check + vite build green
in-container). Next: S.1 `yee-server` (axum WS), S.3 visualization, S.4 Dioxus parity audit.
Earlier: E.2/E.3/E.5a (ADR-0177), E.1 (ADR-0176), E.0 (ADR-0175).*
