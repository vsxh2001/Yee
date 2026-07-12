# ADR-0214: FS.0b.2-GPU — graded (nonuniform) grid on the wgpu/WGSL backend

**Date:** 2026-07-12. **Status:** accepted. **Phase:** FS.0b.2-GPU
(`FULL-SUITE-ROADMAP.md` FS.0 row). **Spec:**
`docs/superpowers/specs/2026-07-12-fs0b2-gpu-graded-design.md`. **Plan:**
`docs/superpowers/plans/2026-07-12-fs0b2-gpu-graded.md`. **Predecessors:**
ADR-0208 (graded CPU kernel), ADR-0210 (graded rules + voxelizer).

## Context

FS.0b.0/0b.1 delivered graded spacings on the CPU backend end to end, but
the heavy board solves that motivated grading run fastest on the GPU
backend, which still multiplied every curl difference by the three scalar
`Params.inv_dx/dy/dz` uniforms. ADR-0208 anticipated this step: "FS.0b.2+
must thread the primal/dual arrays through the WGSL bindings and move the
Unsupported rejection into `gpu.rs`."

## Decision

- **One packed inverse-spacing storage buffer** (binding 8):
  `inv_xp[nx] | inv_yp[ny] | inv_zp[nz] | inv_xd[nx+1] | inv_yd[ny+1] |
  inv_zd[nz+1]`, offsets re-derived in the shader from `nx/ny/nz` (the
  arena idiom). Every update kernel multiplies by the per-cell inverse;
  the divisor **indexing** is the FS.0b.0 CPU kernel's, verbatim
  (H updates → inverse primal at the H sample, E updates → inverse dual).
  The fused CPML ψ corrections reuse the same curl variables, so they
  inherit the graded divisors exactly as `cpml.rs` does — zero extra CPML
  code (the FS.0b.0 scope rule keeps spacing uniform inside absorbers
  anyway).
- **Inverses (multiply), not widths (divide) — deliberately the opposite
  of ADR-0208's CPU choice.** The CPU divides by spacing values because it
  must stay bit-exact against `yee-fdtd`'s dividing FP64 reference. The
  GPU kernel has always multiplied by inverse scalars produced as
  `(1.0 / spec.dx) as f32`; filling the buffer with the identical f64
  expression makes the uniform fill **bit-equal to the old scalars**, so
  the scalar path is unchanged bit-for-bit (every pre-existing FP32 gate
  unaffected by construction) and gate `compute-020` can assert exact
  equality rather than an epsilon. Uploading widths and dividing would
  have perturbed every uniform GPU result by an ulp for no benefit.
- **`Params.inv_dx/inv_dy/inv_dz` removed** (host + WGSL structs) — dead
  scalars alongside a live buffer invite divergence.
- **`GpuFdtd::set_spacings(&GradedSpacings) -> Result<(), ComputeError>`**
  mirrors `CpuFdtd::set_spacings`: panics on invalid input (validate,
  `validate_cpml_layers` against the build-time npml/faces, dispersive
  exclusion, graded Courant check — the CPU contract) plus a
  before-stepping assert (the GPU state is already device-resident);
  returns `ComputeError::Unsupported` for backend-capability gaps (below).
  On success it refreshes, in place via `queue.write_buffer` (buffer
  lengths are functions of the spec, so the bind group never changes):
  the inverse-spacing buffer, the resistive-port α/γ constants (local
  dual-transverse area × primal dz — the exact CPU formulas, recomputed in
  f64 at the stored port cells), and the aperture-port `vcoef` (local
  primal dz / n_columns). The uniform-fill recomputation reproduces the
  build-time constants bit-exactly (same f64 expressions, same order).
- **Scope rejections** (recorded per the brief):
  - **NTFF-DFT + graded → `Unsupported`.** The on-GPU DFT accumulator
    itself is spacing-independent, but the downstream `NtffState` surface
    integration assumes a uniform grid — the same reason FS.0b.0 rejected
    NTFF+graded at the engine. Rejected at the `GpuFdtd` source now.
  - **Aperture port straddling non-uniform z-spacing → `Unsupported`.**
    The aperture kernel folds the modal `∫E_z dz` into one per-port scalar
    `vcoef`; that factoring is exact only when every cell of the port
    shares one primal dz. The CPU's per-cell `dzp[k]` weighting is not
    reproduced; instead `set_spacings` recomputes `vcoef` from the port's
    single local dz and rejects z-taper-straddling ports. This covers the
    real consumer: FS.0b.1 `auto_spacings` substrates are exactly uniform
    in z (`n_sub` cells of `h/n_sub`), and engine apertures live in the
    substrate.
  - Dispersive + graded stays mutually excluded (FS.0b.0 rule, both
    backends); graded-inside-CPML stays rejected by
    `validate_cpml_layers`; resistive-sheet loss (R.0b) remains off the
    GPU entirely (pre-existing).
- **Binding budget:** the spacing buffer is the **8th storage buffer per
  stage — exactly the WebGPU default limit**
  (`maxStorageBuffersPerShaderStage = 8`). The next GPU feature that needs
  a buffer must pack into an existing arena (psi/coeffs/profiles), not add
  a binding.
- **CI:** no `gpu-nightly.yml` change. Its existing
  `cargo test -p yee-compute --release -- --include-ignored --nocapture`
  step runs with default features (`gpu` is default-on) and picks up both
  new gates.

## Measured results (gates)

All numbers below were measured on **llvmpipe (LLVM 20.1.2)** — a software
adapter; the GPU nightly runner re-certifies on real hardware.

- **compute-020** (`gpu_graded_parity.rs::gpu_graded_uniform_parity`,
  fast, non-ignored): the compute-018 drive scenario (CPML npml 4, soft
  source + resistive port + aperture port + 2 probes, 20×16×12, 150
  steps). Uniform-filled `GradedSpacings` vs the GPU's own scalar run:
  **bit-for-bit — 0 differing elements across all six components and
  every probe sample** (the documented "bit-for-bit" arm of the gate, not
  the epsilon fallback). Vs the FP64 CPU backend: max family-rel
  L2 = 4.71e-5 (hx), max family-rel L∞ = 2.75e-4 (hx) — inside the
  compute-002 tolerances (1e-4 / 1e-3).
- **compute-021** (`gpu_graded_parity.rs::gpu_graded_taper_parity`,
  `#[ignore]`, release; picked up by the nightly's `--include-ignored`):
  the compute-019 taper scenario (0.5 → 0.25 → 0.5 mm, ratio 1.122/cell,
  CPML all faces, graded-Courant dt, 560 steps, ~260k cells).
  **Measured 2026-07-12:** CPU↔GPU probe-series parity rel L∞ = 4.709e-6,
  rel L2 = 5.853e-6 (normalized by the CPU trace peak / L2) — **pinned at
  1e-4** (~20× headroom). GPU-side grading reflection via the
  uniform-reference difference method: incident peak 2.721e-3, reflection
  6.322e-6 → **−52.68 dB** — identical to the CPU's ADR-0208 figure at
  the printed precision, under the same pinned −48 dB floor.
  Early-window isolation Δ = **0.000e0 — exactly zero**: before the
  wavefront reaches the grading, the graded and uniform-reference GPU
  runs execute bit-identical FP32 operations (the inverse-spacing values
  are bit-equal upstream), so the difference method isolates the grading
  reflection on the GPU exactly as it does on the CPU. Runtime ~31 s on
  llvmpipe (both gates, release).

## Consequences

- The GPU backend accepts graded jobs; the engine-level rejection
  (`yee-engine/src/lib.rs`, "graded grid (FS.0b) is not on the GPU yet")
  is now obsolete and must be lifted by the dispatcher — out of this
  track's lane. When lifted, graded+GPU engine jobs get resistive and
  aperture ports for free; NTFF+graded remains rejected at both layers.
- `yee-fdtd` (the scalar reference) remains uniform-only; graded GPU runs
  are gated against the graded CPU backend, which is itself gated against
  the uniform reference plus the measured reflection floor (ADR-0208).
- The storage-buffer budget is exhausted; see the binding note above.
