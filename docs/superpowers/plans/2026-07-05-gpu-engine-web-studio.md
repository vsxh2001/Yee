# GPU/CPU compute engine + web studio — implementation plan

**Spec:** `docs/superpowers/specs/2026-07-05-gpu-engine-web-studio-design.md` (ADR-0175)

## Increment 1 — E.0 walking skeleton (this branch)

1. **Docs first**: ADR-0175 + spec + this plan + `ENGINE-STUDIO-ROADMAP.md`; register the ADR in
   `docs/src/SUMMARY.md`; CLAUDE.md gets `yee-compute` in the workspace map and a pointer to the
   new roadmap.
2. **Crate scaffold**: `crates/yee-compute` (workspace member + workspace dep), `gpu` feature
   default-on pulling `wgpu` + `pollster` + `bytemuck`; always-on deps `rayon`, `thiserror`,
   `yee-core` (EPS0/MU0). `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`.
   Dev-deps: `yee-fdtd`, `ndarray` (reference parity test only).
3. **`spec.rs` / `fields.rs`**: `FdtdSpec` (+ staggered shape helpers + Courant assert),
   `Fields` (flat row-major buffers + `idx` helpers + Gaussian-ball constructor used by tests).
4. **`cpu.rs`**: `CpuFdtd` — rayon slab-parallel ports of the uniform lossless `update_h`/
   `update_e` arms; PEC-box interior ranges identical to the reference.
5. **`gpu.rs` + `shaders/fdtd.wgsl`**: `GpuFdtd` per spec §4.3 (six entry points, uniform params,
   chunked submission, staging readback). Graceful `Err(ComputeError::NoAdapter)` path.
6. **`engine.rs`**: `FdtdEngine` enum front (new_cpu / new_gpu / step_n / read_fields /
   backend_name).
7. **Gates**: `tests/cpu_reference_parity.rs` (compute-001, bit-exact) and
   `tests/gpu_cpu_parity.rs` (compute-002, tolerance + self-skip). Both non-ignored, fast.
8. **CI**: add `cargo test -p yee-compute --release` step to `gpu-nightly.yml`.
9. **Verify**: `cargo fmt --check --all`; `cargo clippy -p yee-compute --all-targets -- -D
   warnings`; `cargo test -p yee-compute`; `cargo check -p yee-compute --no-default-features`;
   `cargo check --workspace`.

**DoD:** all of step 9 green locally; `compute-001` bit-exact; `compute-002` green-or-skipped;
docs registered; pushed to the feature branch.

## Increment 2 — E.1 CPML + materials on the engine

- Port CPML (split-field ψ arrays as additional buffers) and per-cell ε/σ/PEC-mask arms to both
  backends. Gate: reproduce `cpml_reflection`'s ≥30 dB criterion via `yee-compute`.
- Pattern file: `crates/yee-fdtd/src/cpml.rs` + `crates/yee-fdtd/tests/cpml_reflection.rs`.

## Increment 3 — E.2 driver parity

- Sources (soft point / plane wave), lumped ports, step-loop driver equivalent to `FdtdDriver`.
- Gate: `fdtd-line-eeff-001` scenario ε_eff via `yee-compute` within the existing ±15% HJ gate,
  plus CPU-vs-GPU agreement on the extracted ε_eff to < 0.5%.

## Increment 4 — S.0/S.1 engine service

- `crates/yee-engine`: `Job`, `JobEvent` (Progress/Partial/Done/Error), `EngineHandle`
  (submit/cancel/subscribe), threaded executor over `yee-compute` + existing solvers; serde
  round-trip tests.
- `crates/yee-server`: axum, `GET /healthz`, `WS /v1/jobs` (submit → event stream); `yee serve`
  CLI subcommand. Gate: end-to-end WS test running a small FDTD job.

## Increment 5 — S.2+ Tauri studio

- `studio/` (Tauri 2 + React + TS + Vite) with its own ADR at kickoff: project scaffold, engine
  job panel (submit/progress/cancel), S-param plot from a finished job. npm toolchain note goes
  into CLAUDE.md §7 when this lands.
- Dioxus `yee-studio-web` is feature-frozen from this point (ADR-0175); parity audit + retirement
  decision is S.4 with its own ADR.

## Lanes (for multi-track dispatch)

- Engine lane: `crates/yee-compute/**` (+ later `crates/yee-engine/**`, `crates/yee-server/**`)
- Reference lane (read-only from engine lane): `crates/yee-fdtd/**`
- CI lane: `.github/workflows/**`
- Docs lane: `docs/**`, `ENGINE-STUDIO-ROADMAP.md`, `CLAUDE.md`
- Studio lane (later): `studio/**`

Escape hatch (standard): blocked > 15 min → surface and stop.
