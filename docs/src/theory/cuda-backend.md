# CUDA Backend — Theory of Operation

This page is the theory-of-operation reference for Yee's CUDA layer,
implemented in the `yee-cuda` crate. Same audience as the planar-MoM
and FDTD pages (an engineer reading source code with a textbook open),
same conventions (plain-text math, inline citations, source-file
references in inline code).

Unlike the physics chapters, the subject here is half numerics and half
software architecture: a trait that exists for migration insurance, a
dense LU path chosen because nothing iterative would converge, a
feature gate that has to stay green on a laptop without a CUDA toolkit,
and a CI workflow that only certifies the GPU path when a self-hosted
runner is available.

## 1. Overview

Yee uses CUDA in one place that ships today and several places on the
roadmap:

- **Shipping (Phase 1.5):** dense complex-double LU for the planar-MoM
  impedance matrix (`cusolverDnZgetrf` + `cusolverDnZgetrs`).
- **Planned (Phase 2.fdtd):** FDTD `E`/`H` and CPML update kernels via
  NVRTC.
- **Planned (Phase 4):** NCCL halo exchange for multi-GPU FDTD.

All three paths sit behind a single trait, `yee_cuda::backend::Backend`,
defined in `crates/yee-cuda/src/backend.rs`. No other crate in the
workspace imports `cudarc` directly — the `Backend` trait is the only
public surface. The Cargo `cuda` feature on `yee-cuda` toggles between
the real cudarc-backed implementation and a stub that returns
`Error::NotEnabled` from every entry point, so the workspace builds
green on CI hosts without a CUDA toolkit (the default case for CPU-only
laptops and the public `ubuntu-latest` GitHub runner).

Why this much ceremony for one solver call? `cudarc` is pre-alpha and
treats `0.x → 0.y` as breaking, so the underlying binding is the
volatile part of the stack. Wrapping it in a Yee-owned trait localises
that volatility to one file. ADR-0006
(`docs/src/decisions/0006-cudarc-prealpha-pin.md`) records the
trade-off in full.

## 2. The `Backend` trait abstraction

Three observable facts about `cudarc` set the design:

1. **It is the only currently shippable Rust CUDA binding.** `cust`
   was unmaintained from 2022 through early 2025 and as of 2026 still
   has no `crates.io` release; `rustacuda` has been abandoned since
   2021. `cudarc` ships, runs in production at Candle (Hugging Face)
   and `mistral.rs`, and covers driver + cuBLAS + cuSOLVER + cuRAND +
   cuFFT + NCCL + NVRTC.
2. **It self-describes as pre-alpha and exercises that label.** The
   `0.13 → 0.14` bump reorganised device-pointer types and renamed
   half a dozen public functions; `0.16 → 0.17` changed how
   `LaunchAsync` is parameterised. By Cargo convention every
   `0.x → 0.y` step *is* breaking, and the maintainer acts
   accordingly.
3. **Yee's GPU surface is narrow.** Phase 1.5 needs two cuSOLVER
   calls, device enumeration, NVRTC compilation, and a handful of
   stream/memcpy primitives — small enough to abstract without
   inventing a parallel CUDA framework.

The trait, paraphrased from `crates/yee-cuda/src/backend.rs`:

```text
trait Backend {
    fn device_count() -> Result<usize>;
    fn device_props(i: usize) -> Result<(String, (u8, u8), u64)>;
    fn nvrtc_compile(src: &str, name: &str) -> Result<Vec<u8>>;
    fn cusolver_zgetrf(a: &[Complex64], n: usize)
        -> Result<(Vec<Complex64>, Vec<i32>)>;
    fn cusolver_zgetrs(lu: &[Complex64], pivots: &[i32],
                       b: &[Complex64], n: usize, nrhs: usize)
        -> Result<Vec<Complex64>>;
}
```

The concrete implementation, `CudarcBackend`, is gated by
`#[cfg(feature = "cuda")]` and is the only type in the workspace that
imports `cudarc::*`. The Phase 0 trait uses associated functions
(no `&self`) because no per-instance state is required yet; Phase 1
upgrades to `&self` methods on a context handle when runtime backend
selection (e.g. a mock backend for non-GPU CI smoke tests) becomes
useful.

If `cudarc 0.20` lands tomorrow with a breaking change to
`CudaDevice`, the PR to absorb it touches `backend.rs` and
`cusolver.rs` and nothing else; nothing in `yee-mom`, `yee-fdtd`, or
future GPU consumers changes. The same property holds for a migration
to the rebooted `cust`, a vendored fork, or an entirely different
runtime.

## 3. cuSOLVER LU path: `Zgetrf` → `Zgetrs`

The planar-MoM impedance matrix `Z` is dense, complex-symmetric (not
Hermitian — `Z = Z^T` in the symmetric Galerkin formulation, not
`Z = Z^H`), and on the validation gates we care about (mom-001
half-wave dipole, 4224 RWG unknowns at the published nightly mesh) it
is small enough that direct factorisation dominates iterative solvers
in both wall-clock and reliability.

### 3.1 Why direct, not iterative

MoM impedance matrices are a textbook bad case for Krylov-subspace
methods. The matrix is **dense** (each iteration costs `~2n²` flops,
the same as one column of a direct factorisation), **complex-symmetric
but not Hermitian** (ruling out CG and MINRES; leaving GMRES, BiCGStab,
QMR), **ill-conditioned at resonance** (the very regime where the
solver is most useful, inheriting the cavity-mode null space of the
underlying EFIE), and **lacks a cheap preconditioner** — incomplete
LU and algebraic multigrid are sparse-matrix tools. The MoM
literature's response (MLFMM, H-matrices, AIM, FIPWA) collapses the
dense matrix-vector product to `O(n log n)` and only then does GMRES
become competitive. That family is a Phase 4+ research programme.

For `n ≲ 10 000` on a single GPU, `Zgetrf` finishes in a few seconds
with `O(n³)` flops, the factorisation is reusable across right-hand
sides (frequency or port sweeps hit `Zgetrs` only), and no convergence
question survives. Direct wins on every axis until either `n` outgrows
GPU memory (~`60 000` on a 40 GB device storing `16 n²` bytes of LU
plus workspace) or a fast-multipole compression scheme arrives.
Neither has happened yet.

### 3.2 What `Zgetrf` and `Zgetrs` actually do

`Zgetrf` is the complex-double Gaussian elimination with partial
pivoting that LAPACK has shipped since 1992 and cuSOLVER has shipped
on the GPU since CUDA 7. It factors a column-major `n × n` matrix
`A` into

```text
P A = L U
```

where `L` is unit lower-triangular, `U` is upper-triangular, and `P`
is the row permutation produced by partial pivoting (returned as a
1-indexed LAPACK pivot vector). The cost is `(2/3) n³` complex
multiply-adds plus `O(n²)` pivot swaps. cuSOLVER's implementation is
a right-looking blocked variant that batches the trailing-matrix
update into `zgemm` calls, so the asymptotic cost is dominated by
cuBLAS-3 BLAS rather than the panel factorisation.

Pivoting matters: without it even a perfectly conditioned system can
produce an exponentially-growing `U` and floating-point garbage
(Trefethen and Bau §21 work the classic Wilkinson example). The
partial-pivoting growth factor is bounded by `2^(n-1)` in the worst
case but is empirically `O(n)` — Higham §9 gives the error model and
empirical-growth literature.

`Zgetrs` consumes `(L, U, P)` and one or more right-hand sides `B`
and emits `X = A⁻¹ B` in three triangular solves: apply `P` to `B`,
forward-solve `L Y = P B`, back-solve `U X = Y`. Cost is `O(n²)` per
right-hand side. Reusing one `Zgetrf` across many `Zgetrs` is the
whole point: a 401-point frequency sweep pays one `O(n³)`
factorisation and 401 `O(n²)` solves.

`crates/yee-cuda/src/cusolver.rs` wraps both calls behind a safe
`DenseLuComplex` type that owns the device-side LU factors, pivot
vector, cuSOLVER dense handle, and stream so `solve` can be called
repeatedly without rebuilding any of them. The module locally opts
back into `unsafe_code` because cuSOLVER's `sys::*` functions are
`unsafe fn`; each `unsafe` block carries a local `SAFETY` comment.
The rest of the crate still `#![deny(unsafe_code)]`.

## 4. The `cuda` feature gate

The `cuda` Cargo feature on `yee-cuda` is **off by default**. With
the feature off, the crate still compiles, exports the same public
types (`Device`, `Error`, `DenseLuComplex`, `Backend`), and every
constructor or method that would have hit the CUDA driver returns
`Error::NotEnabled`. The behaviour is summarised on `Device::list`'s
doc-comment table:

```text
feature `cuda`   devices visible   return
off              n/a               Err(Error::NotEnabled)
on               0                 Ok(vec![])
on               n                 Ok(Vec<Device>) len n
```

Why not auto-detect at runtime? `cudarc` is an optional Cargo
dependency precisely to avoid pulling the CUDA driver crate (and its
compile times) into every CI build; without the feature off-by-default
the workspace would not build on the public `ubuntu-latest` GitHub
runner. Runtime detection would itself require linking against the
CUDA driver library unconditionally — the failure mode we are trying
to avoid. And `Error::NotEnabled` is a precise, grep-able signal that
"this code path was not exercised on this host", strictly more useful
than the silent empty-list fallback an auto-detector would produce.

The pattern matches `gmsh` on `yee-mesh`: features needing an
external SDK or driver default off and return a `NotEnabled` stub on
the no-feature path.

## 5. Hardware-gated CI

A green `cargo test --workspace` on a CPU-only host certifies the
`Error::NotEnabled` paths and nothing else. The actual GPU path —
matrix uploads, cuSOLVER calls, LU residuals — is exercised only by
tests marked `#[ignore]` that require `--include-ignored` and the
`cuda` feature:

```text
cargo test -p yee-cuda --features cuda --release -- --include-ignored
```

These tests live in `crates/yee-cuda/src/cusolver.rs` and include a
2 × 2 identity smoke test plus a 64 × 64 random Hermitian
positive-definite residual check (`||A x − b|| / ||b|| < 1e-10`).
They are not run on the public CI.

The certification path is `.github/workflows/gpu-nightly.yml`. The
workflow is gated by the repository variable `YEE_GPU_RUNNER_ENABLED`:

```text
if: vars.YEE_GPU_RUNNER_ENABLED == 'true'
runs-on: ${{ vars.YEE_GPU_RUNS_ON || 'self-hosted-gpu-placeholder' }}
```

If the variable is unset (the default for any fork), the workflow
no-ops and a fork without GPU hardware sees no red nightly runs. To
enable: register a self-hosted GitHub Actions runner with a CUDA
12.4+ GPU and the `gpu` label, then set
`YEE_GPU_RUNNER_ENABLED = "true"` in *Settings → Secrets and
variables → Actions → Variables*. The workflow runs `nvcc --version`,
builds with `--features cuda --release`, and runs the ignored tests
nightly at 04:00 UTC (manually dispatch-able from the Actions tab).

A self-hosted runner is the only path to actually certifying the GPU
code: GitHub's hosted runners do not expose GPUs, the CUDA driver
needs kernel-mode access incompatible with shared infrastructure, and
PTXAS-only emulators do not validate cuSOLVER's blocked-LU kernels.
The result is a two-tier model: CPU CI proves the `NotEnabled` paths
and that the public API compiles; the nightly GPU job proves the
numerics. Treat the green-tick on a CPU-only CI run accordingly.

## 6. Numerical precision: why `complex<f64>`, not `complex<f32>`

cuSOLVER ships both single-precision (`Cgetrf`) and double-precision
(`Zgetrf`) complex-LU paths. We use double throughout
(`num_complex::Complex64`, equivalently
`cuDoubleComplex = struct { double x, y; }`). The reason is the
validation gate, not raw precision worship.

mom-001 demands `Z ≈ 87 + j41 Ω` to within ±5% on `Re(Z)` and ±10% on
`Im(Z)` against NEC-4. The matrix has 4224 unknowns at the production
mesh. Standard backwards-error analysis for partial-pivoting LU
(Higham §9) bounds the relative residual by

```text
||A x − b|| / (||A|| ||x||)  ≤  ρ_n · γ_n · u
```

with `u` unit round-off, `γ_n ≈ n u / (1 − n u)`, and `ρ_n` the
growth factor. For `n = 4224` and IEEE `binary32` (`u ≈ 6 × 10⁻⁸`),
`n u ≈ 2.5 × 10⁻⁴` *before* the growth factor; empirical
partial-pivoting growth `ρ_n ≈ O(n)` pushes the residual into the
`10⁻³` range — already at or beyond the ±5% Re-tolerance once port
de-embedding amplifies it. The `Im(Z)` budget is tighter still
because the imaginary part is a small difference of larger numbers
near resonance (capacitive minus inductive storage).

For IEEE `binary64` (`u ≈ 1.1 × 10⁻¹⁶`), `n u ≈ 4.6 × 10⁻¹³`, four
to five orders of magnitude below the validation tolerance. The
residual budget vanishes. mom-001 passes at f64 with the published
NEC-4 figure of merit; no evidence it would pass at f32 on any mesh
dense enough to resolve the gap-region near-field.

The throughput trade is real but small at our scale. On H100/A100
f64 is 1/2 of f32 in the cuBLAS-3 kernels that dominate `Zgetrf`; on
consumer RTX 40-series it is 1/32 or worse. Even at the consumer
ratio, a 4224 × 4224 `Zgetrf` finishes under a minute and is dwarfed
by matrix-fill. The f32 wins would be invisible against fill time
and the backwards-error losses catastrophic against the gate. Worth
re-evaluating only when the fill step itself moves onto the GPU.

## 7. Limitations and known gaps

- **Single-GPU only.** Both cuSOLVER calls target device 0. Multi-GPU
  LU (MAGMA, ScaLAPACK-on-GPU, cuSOLVER MultiGPU) is post-Phase 4.
  The single-GPU memory ceiling (~`n ≲ 60 000` on a 40 GB device) is
  the Phase 1.x upper bound.
- **No iterative path.** No CG, GMRES, BiCGStab, QMR, or
  fast-multipole compression. Problems above the LU ceiling get
  either a faster fill path, MLFMM, or FDTD instead.
- **Round-tripped host memory.** The trait-level entry points take
  and return host buffers; `DenseLuComplex` keeps LU on-device
  between solves but the trait surface does not. Phase 2 introduces
  a device-buffer type so `yee-mom` can hand `Z` directly to cuSOLVER
  without the host round-trip.
- **GPU path is feature-gated, not auto-detected.** A laptop with a
  CUDA-capable GPU but `cuda` off will *not* use the GPU.

## 8. References

- Trefethen, L. N., and Bau, D. *Numerical Linear Algebra.* SIAM,
  1997. (§20 Gaussian elimination, §21 pivoting, §23 stability of
  LU.)
- Higham, N. J. *Accuracy and Stability of Numerical Algorithms.*
  2nd ed. SIAM, 2002. (§9 LU factorisation and linear equations;
  growth factor; backward error bounds.)
- Anderson, E., et al. *LAPACK Users' Guide.* 3rd ed. SIAM, 1999.
  (`Zgetrf` / `Zgetrs` algorithmic specification.)
- NVIDIA. *cuSOLVER Library Documentation.*
  <https://docs.nvidia.com/cuda/cusolver/>. (`cusolverDnZgetrf`,
  `cusolverDnZgetrs`, workspace query convention, info return code
  semantics.)
- Coreylowman. `cudarc` crate.
  <https://github.com/coreylowman/cudarc>. (Pre-alpha status,
  release history.)
- ADR-0006 (`docs/src/decisions/0006-cudarc-prealpha-pin.md`) —
  `cudarc` pinning rationale and `Backend` trait migration
  guarantee.
- ADR-0005 (`docs/src/decisions/0005-nec4-vs-balanis-mom-001.md`) —
  mom-001 NEC-4 reference value and tolerance budget, which sets the
  f64 floor used in §6.
