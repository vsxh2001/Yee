# yee-cuda — Roadmap

## Phase 0 (months 0–6)
- [ ] Build matrix: CUDA 12.4 + 13.0, Linux + Windows
- [ ] `Device::list()` returning name, ordinal, compute capability, mem total
- [ ] `Context`, `Stream`, `Module` RAII handles
- [ ] NVRTC helpers: `compile_kernel(src, name) -> Module`
- [ ] `kernels/hello.cu` — trivial stencil that proves the pipeline
- [ ] Smoke test in `validation/`: enumerate device + run hello kernel + verify output
- [ ] Internal `backend` trait so cudarc can be swapped (defense against pre-alpha churn)

## Phase 1 (months 6–18)
- [ ] cuSOLVER `Zgetrf` + `Zgetrs` safe wrappers (complex double dense LU)
- [ ] cuBLAS `Zgemm` + `Zgemv` wrappers
- [ ] Iterative solver kernels: SpMV with cuSPARSE, GMRES preconditioner
- [ ] Multi-GPU dense LU via cuSOLVERMg (behind `multi-gpu` feature)
- [ ] cuFFT helpers for near-to-far-field projection prep
- [ ] Performance harness: GEMM GFLOPS vs published H100/A100 numbers

## Phase 2 (months 18–30)
- [ ] FDTD update kernels in `kernels/fdtd/`: E and H steps with shared-memory tiling
- [ ] NCCL boundary exchange wrappers
- [ ] Mixed precision (FP32 fields, FP64 accumulators) helpers
- [ ] Domain-decomp grid utilities

## Phase 3 (months 30–42)
- [ ] cuBLASLt mixed-precision Tensor-Core GEMM for surrogate inference
- [ ] Optional cuDNN integration when surrogate models need it

## Validation gates per phase
- Phase 0: `validation/hello-stencil` passes on a real GPU; CI gates only run on a self-hosted runner with NVIDIA driver.
- Phase 1: `cusolverDnZgetrf` reproduces a known 4×4 → 256×256 LU within numerical tolerance; GEMM GFLOPS within 20% of published peak.
- Phase 2: FDTD update kernel benchmark beats openEMS on the same problem size on the same GPU class (NVIDIA RTX 4090 reference).

## Watch-outs
- `cudarc` pre-alpha — pin exact minor version (`=0.19.x`) and route all access through this crate's API.
- CUDA 13.0 dropped some legacy targets — verify compute capability ≥ 7.0 (Volta).
