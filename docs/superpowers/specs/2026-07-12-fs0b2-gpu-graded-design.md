# FS.0b.2-GPU — graded (nonuniform) grid on the wgpu/WGSL GPU backend

**Phase:** FS.0b.2-GPU (full-suite track, `FULL-SUITE-ROADMAP.md` FS.0 row).
**Plan:** `docs/superpowers/plans/2026-07-12-fs0b2-gpu-graded.md`. **ADR:**
`docs/src/decisions/0214-fs0b2-gpu-graded.md`. **Predecessors:** ADR-0208
(graded CPU kernel, FS.0b.0), ADR-0210 (graded rules + voxelizer, FS.0b.1).

## Why

FS.0b.0/0b.1 delivered graded spacings end to end on the **CPU** backend:
`GradedSpacings` on the job protocol, `CpuFdtd::set_spacings`, primal/dual
divisors in every kernel, mesh rules, and the graded voxelizer. The payoff —
refining only the staircase-limited feature bands — is real, but the heavy
board solves that motivated it (ADR-0204: the next uniform pass costs ~2.4 h
at 19 M cells) run fastest on the GPU backend, which still assumes uniform
spacing (`Params.inv_dx/dy/dz` scalars) and is guarded by an engine-level
rejection (`yee-engine/src/lib.rs`: "graded grid (FS.0b) is not on the GPU
yet"). FS.0b.2-GPU threads the per-axis spacing arrays through the WGSL
kernels so a graded job can run on the GPU, gated by uniform-fill parity and
by the compute-019 taper scenario cross-backend.

## Design

### 1. Spacing buffer (binding 8)

One new read-only storage buffer holds the six **inverse** spacing arrays,
packed in a fixed order the shader re-derives from `nx/ny/nz` (the arena
idiom):

```
inv_sp: inv_xp[nx] | inv_yp[ny] | inv_zp[nz]
      | inv_xd[nx+1] | inv_yd[ny+1] | inv_zd[nz+1]
```

`*p` = primal cell widths, `*d` = dual spacings, both from the FS.0b.0
`SpacingArrays` (study `spec.rs`: dual = mean of adjacent primals in the
interior, the single adjacent primal at domain edges). Inverses are computed
host-side in **f64** (`1.0 / d`) and narrowed once to f32 — exactly how the
old `Params.inv_dx = (1.0 / spec.dx) as f32` scalar was produced.

**Why inverses (multiply) and not widths (divide), unlike the CPU kernel?**
ADR-0208 chose literal division on the CPU because the CPU must stay
bit-exact against `yee-fdtd`'s dividing reference. The GPU kernel has always
*multiplied by precomputed inverse* scalars; uploading inverses produced by
the identical f64 expression makes the uniform fill **bit-equal to the old
scalars**, so every existing FP32 gate (compute-002/005/…) and the new
uniform-fill gate are bit-for-bit unchanged by construction — the strongest
available "no regression" statement on a tolerance-gated backend. Uploading
widths and dividing would perturb every uniform result by an ulp and demote
the uniform-parity gate to an epsilon test.

The buffer is created at build time with the uniform fill (so the bind group
never changes); `set_spacings` refreshes its contents via
`queue.write_buffer` (all lengths are functions of the spec, so no rebuild).
`Params.inv_dx/inv_dy/inv_dz` are **removed** (host struct + WGSL struct) —
dead scalars invite divergence.

Binding budget: this is the 8th storage buffer per stage — exactly the
WebGPU default limit (`maxStorageBuffersPerShaderStage = 8`). The next
buffer must be packed into an existing arena; noted in the ADR.

### 2. WGSL kernels — divisor indexing mirrors `cpu.rs` exactly

Six accessors (`inv_xp(i)`, `inv_yp(j)`, `inv_zp(k)`, `inv_xd(i)`,
`inv_yd(j)`, `inv_zd(k)`) replace the three `p.inv_*` uniforms at every use.
The index mapping is the FS.0b.0 CPU kernel's, verbatim:

| kernel      | curl term  | divisor (CPU)  | WGSL factor  |
|-------------|-----------|----------------|--------------|
| `update_hx` | ∂E_y/∂z   | `dzp[k]`       | `inv_zp(k)`  |
| `update_hx` | ∂E_z/∂y   | `dyp[j]`       | `inv_yp(j)`  |
| `update_hy` | ∂E_z/∂x   | `dxp[i]`       | `inv_xp(i)`  |
| `update_hy` | ∂E_x/∂z   | `dzp[k]`       | `inv_zp(k)`  |
| `update_hz` | ∂E_x/∂y   | `dyp[j]`       | `inv_yp(j)`  |
| `update_hz` | ∂E_y/∂x   | `dxp[i]`       | `inv_xp(i)`  |
| `update_ex` | ∂H_z/∂y   | `dyd[j]`       | `inv_yd(j)`  |
| `update_ex` | ∂H_y/∂z   | `dzd[k]`       | `inv_zd(k)`  |
| `update_ey` | ∂H_x/∂z   | `dzd[k]`       | `inv_zd(k)`  |
| `update_ey` | ∂H_z/∂x   | `dxd[i]`       | `inv_xd(i)`  |
| `update_ez` | ∂H_y/∂x   | `dxd[i]`       | `inv_xd(i)`  |
| `update_ez` | ∂H_x/∂y   | `dyd[j]`       | `inv_yd(j)`  |

H updates read **primal** (curl-E spans one primary cell), E updates read
**dual** (curl-H spans the distance between H samples). The fused CPML ψ
corrections reuse the same curl variables, so they inherit the graded
divisors exactly as `cpml.rs` does on the CPU — no separate CPML work
(FS.0b.0's scope rule keeps spacing uniform *inside* absorbing layers
anyway, enforced by `validate_cpml_layers`).

### 3. `GpuFdtd::set_spacings(&GradedSpacings) -> Result<(), ComputeError>`

Mirrors `CpuFdtd::set_spacings` with the same validation ladder, plus
backend-capability rejections in the R.3 `Unsupported` idiom:

- **Panics** (mirroring the CPU API's contract for invalid input):
  `GradedSpacings::validate` failure, `validate_cpml_layers` failure
  (npml/faces captured at build), `dt > graded.courant_limit()`, dispersive
  map attached (the CPU asserts the same), or called after stepping
  (`steps_taken > 0` — GPU buffers are already history-laden).
- **`Err(ComputeError::Unsupported)`** for GPU-capability gaps:
  - **NTFF DFT** (`with_ntff_dft`): the on-GPU accumulator itself is
    spacing-independent, but the downstream `NtffState` surface integration
    assumes a uniform grid (the same reason the engine rejects NTFF+graded),
    so graded+DFT is rejected at the source per the brief.
  - **Aperture port spanning non-uniform dz:** the aperture kernel folds the
    modal `∫E_z dz` into one per-port scalar `vcoef = dz/n_columns`; that
    factoring is exact only when every cell of the port shares one primal
    dz. `set_spacings` recomputes `vcoef` from the port's (single) local dz
    and rejects ports whose cells straddle a z-taper. This covers the real
    consumer: FS.0b.1's `auto_spacings` makes the substrate exactly uniform
    in z (`n_sub` cells of `h/n_sub`), and engine apertures live in the
    substrate.
- **Refreshes:** the spacing buffer (inverses of `SpacingArrays::graded`),
  per-port resistive `alpha`/`gamma` (local dual-x·dual-y area × primal dz —
  the CPU's exact formulas, computed in f64 at the stored port cells), and
  aperture `vcoef` — all via `queue.write_buffer` at the known `drv_data`
  offsets. The stored `Drive` copy provides cells/resistances.

### 4. Scope (exactly FS.0b.0's)

Unchanged from the CPU backend: graded-inside-CPML rejected by
`validate_cpml_layers`; dispersive+graded mutually excluded at attach time;
NTFF+graded rejected (here at `GpuFdtd`, since the DFT accumulator lives in
this backend). Soft sources, probes, resistive ports, aperture ports
(uniform-dz), per-cell materials, PEC masks, per-face CPML all work graded.
The engine-level rejection at `yee-engine/src/lib.rs` (~line 830) is lifted
by the dispatcher **after** this track merges — out of this lane.

### 5. Gates

Both in `crates/yee-compute/tests/gpu_graded_parity.rs`,
`#![cfg(feature = "gpu")]`, self-skipping on `ComputeError::NoAdapter`
(the `gpu_cpu_parity.rs` idiom — prints SKIPPED, returns green; the GPU
nightly is the certifying environment).

- **compute-020 `gpu_graded_uniform_parity`** (fast, non-ignored): the
  compute-018 drive scenario (asymmetric dims, CPML, soft source +
  resistive port + aperture port + probes) run three ways — CPU FP64 with
  uniform-filled `GradedSpacings`, GPU with the same spacings, GPU scalar
  (no `set_spacings`). Asserts (a) GPU-graded vs GPU-scalar **bit-for-bit**
  (exact equality on all six widened components and every probe sample —
  achievable because the uniform inverse fill is bit-equal to the old
  scalars; this is the off-by-one tripwire), and (b) GPU-graded vs CPU
  within the compute-002 tolerances (family-rel L2 < 1e-4, L∞ < 1e-3).
- **compute-021 `gpu_graded_taper_parity`** (`#[ignore]`, release; picked
  up by the nightly's `--include-ignored`): the compute-019 taper scenario
  (0.5→0.25→0.5 mm geometric taper, CPML all faces, graded-Courant dt) run
  graded on both backends; the full 560-sample probe time series must match
  within measured-then-pinned FP32 tolerances (normalized by the CPU peak),
  and the GPU trace must stay finite (graded-GPU stability). If the local
  adapter (llvmpipe) permits, also measure the GPU-vs-GPU-uniform
  reflection level against the −48 dB compute-019 floor and record it in
  the ADR.

### 6. CI

`gpu-nightly.yml` already runs
`cargo test -p yee-compute --release -- --include-ignored --nocapture` with
default features (`gpu` is default-on), which picks up both gates — **no
workflow change needed** unless measurement shows otherwise.

## Out of scope

Graded NTFF on GPU (rejected, above); dispersive+graded (both backends,
FS.0b.0 rule); aperture ports straddling a z-taper; graded resistive-sheet
loss (R.0b is not on the GPU at all); lifting the engine-level rejection
(dispatcher, post-merge); `yee-fdtd` reference changes (stays uniform-only).
