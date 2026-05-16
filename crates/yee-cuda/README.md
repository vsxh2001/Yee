# yee-cuda

> Single CUDA host-binding layer for the whole Yee workspace.

All Yee CUDA access goes through this crate. Direct `cudarc` imports outside `yee-cuda` are forbidden — that gives us **one swap point** if cudarc breaks (see [TECH_STACK.md § Known watch-outs](../../TECH_STACK.md#known-watch-outs)).

## Why a wrapper layer

- `cudarc` self-describes as "pre-alpha" and has shipped breaking minor releases (0.13 → 0.14).
- Some libraries (cuSPARSE) only have `result`/`sys` modules in cudarc; we wrap them once, here.
- Future-proofs against a swap to a different backend (e.g. a future `cust` rebirth) without touching solver code.

## Scope

### Phase 0
- `Device::list()` — enumerate visible CUDA devices via the driver API
- Context + stream RAII handles
- NVRTC compile + module load helpers (CUDA-C → PTX → loaded module)
- "Hello world" stencil kernel in `kernels/hello.cu`
- CI builds for CUDA 12.4 **and** 13.0

### Phase 1
- cuSOLVER dense LU: `Zgetrf`, `Zgetrs` wrappers for MoM
- cuBLAS GEMM/GEMV wrappers for MoM matrix-fill aggregation
- Iterative GMRES kernel scaffolding (n ≥ 50k path)
- Multi-GPU dense LU via cuSOLVERMg behind a feature flag

### Phase 2
- NCCL boundary exchange for FDTD domain decomposition
- Shared-memory tiled E/H update kernels in `kernels/fdtd/`

### Phase 3
- cuDNN / cuBLASLt mixed-precision GEMM for surrogate training/inference

## Feature flags

| Flag | Effect |
|------|--------|
| (none) | CPU-only build; APIs return `Error::NotEnabled`. CI safe. |
| `cuda` | Link against CUDA Toolkit; real device enumeration and kernel compile. |

The `dynamic-loading` mode of cudarc means **no CUDA libs needed at build time**, only at runtime — keeps CI lean.

## Validation

See [`validation/README.md`](validation/README.md). Hardware-gated: skipped on hosts without an NVIDIA driver.

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
