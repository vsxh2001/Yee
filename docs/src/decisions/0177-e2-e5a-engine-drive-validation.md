# ADR-0177: E.2–E.5a — engine drive layer + strong-reference validation sweep

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0175 (E.0), ADR-0176 (E.1), `ENGINE-STUDIO-ROADMAP.md`.

---

## Context

After E.1, `yee-compute` could step real materials in open domains but had no way to *drive* a
problem (sources/ports) or *observe* it (probes), so no end-to-end workload could run on the
engine. This increment adds the drive layer and then walks the engine through the strongest
known-result references available to the project, per phase plan E.2/E.3/E.4/E.5.

## Decision

**E.2 — drive layer.** New public surface: `Waveform` (`Gaussian`,
`GaussianPulse` — verbatim ports of `sources::gaussian_pulse_ez` and
`SourceWaveform::GaussianPulse`), `SoftSource` (any E component, injected between the H and E
half-steps), `ResistivePort` (the validated pure-resistor semi-implicit lumped update from
`LumpedRlcPort`, applied after the E boundary phase), `Probe`, and `Drive`, consumed by
`with_drive` constructors. On the **GPU** the entire drive is precomputed host-side in f64 for
the whole run and uploaded once (amplitude/EMF tables + per-port α, γ); an **on-GPU step
counter** (bumped by a 1-thread dispatch at the end of each step) indexes the tables and the
probe output region, so driven stepping stays fully chunked with **zero per-step host
round-trips** — 30 000-step driven runs submit in 64-step chunks exactly like undriven ones.
Two new storage bindings (drive index + drive data) bring the layout to 7 storage + 1 uniform,
still inside WebGPU's 8-per-stage default limit.

**E.3 — precision policy.** FP32-on-GPU / FP64-on-CPU stands. WGSL has no `f64` (naga
implements the WebGPU spec), so "FP64 on GPU via SHADER_F64" is **not** a WGSL-reachable path —
revisit only if a SPIR-V passthrough path is ever justified. Instead the policy is
*characterized*: gate `compute-009` bounds long-run drift (closed PEC box, energy-conserving,
so round-off cannot hide in decay).

**E.4 — benchmarks.** `yee-bench` gains `compute_step` (scalar reference vs rayon CPU vs GPU
incl. readback at 64³). Numbers from this 4-core container are recorded below and are **mixed**
— reported as measured, with the optimization follow-up left open.

**E.5a — far-field on the engine.** The engine owns the stepping; `yee_fdtd::NtffState`
remains the reference near-to-far transform, consuming the engine's fields through a host-side
grid adapter each step. A first-class (GPU-accumulated) NTFF port and the dispersive-ADE port
are explicitly **queued, not shipped** (roadmap E.5b/E.5c): no engine consumer needs them yet —
the filter pipeline's EM-in-loop is non-dispersive FR-4 — and E.1's recipe (bit-exact CPU port
+ arena-fused GPU kernels) is the established path when they're pulled.

## Gates and measured results (all green; strong references in bold)

| Gate | Reference | Result |
|---|---|---|
| `compute-007` (fast, default suite) | `WalkingSkeletonSolver` + `LumpedRlcPort` driven step | **bit-exact** (fields + probe series, max Δ == 0.0) |
| `compute-006` cavity TE₁₀₁, `#[ignore]` | **Analytic Pozar §6.3** `f₁₀₁ = (c/2)√(a⁻²+d⁻²)` | CPU **−0.063 %**, GPU (llvmpipe) **−0.063 %**, CPU↔GPU peak Δ **0.0000 %** (30 000 driven GPU steps, chunked) |
| `compute-008` line-eeff, `#[ignore]` | **Hammerstad–Jensen / Pozar ε_eff** (`yee_layout::eps_eff`), the fdtd-line-eeff-001 scenario | ε_eff 3.3293 vs 3.3249 → **0.132 %** (≤ 15 % gate), 726×78×40 grid, 3 367 steps, **88.6 s** release on 4 cores (the scalar original needs a memory-boxed multi-minute run) |
| `compute-009` FP32 drift, `#[ignore]` | FP64 CPU backend over 10⁴ steps | family-rel L2 3.2e-6 … 1.8e-5 — matches the √N random-walk prediction, 100× inside the 1e-3 gate |
| `compute-010` NTFF dipole, `#[ignore]` | **Analytic sin θ pattern** (endfire null) | broadside/endfire **327.9 dB** (≥ 20 dB gate) |

**E.4 numbers (this 4-core container, 64³ vacuum step):** scalar reference 3.87 ms; rayon CPU
11.1 ms @ 1 thread / 4.97 ms @ 4 threads (2.2× internal scaling, **0.78× vs scalar**); GPU on
llvmpipe (software Vulkan) 770 ms per 64 steps + readback — not a throughput datapoint.
**Finding:** the flat-buffer kernels are ~2.8× slower single-threaded than the `ndarray`
reference (bounds-checked `idx3` arithmetic on every neighbor access is the prime suspect).
Optimization (iterator/row-hoisted inner loops, keeping bit-exactness) is the open E.4
follow-up; the roadmap row stays **partial** until re-measured on real hardware.

## Consequences

- The engine now runs driven, probed, open-domain workloads end-to-end — including the filter
  pipeline's full-wave gate — with every layer certified against either the reference
  implementation (bit-exact) or a published/analytic result.
- CI: new `compute-engine-gates` release job runs the `#[ignore]`'d gates on hosted runners
  (GPU-dependent ones self-skip); the GPU nightly now runs `--include-ignored` so
  compute-002/005/006/009 bite on real hardware.
- `GpuFdtd::with_drive` requires `max_steps` up front (tables sized for the whole run); driven
  runs beyond it panic rather than silently reading garbage.
- E.5b (GPU NTFF), E.5c (dispersive ADE), and the E.4 kernel optimization are the queued
  engine work; the track is otherwise ready for S.0 (`yee-engine` job API).
