# yee-cuda — Validation

Hardware-gated. CI skips on hosts without an NVIDIA driver; full run on self-hosted GPU runners.

## Cases

| ID | Description | Reference | Tolerance | Phase |
|----|-------------|-----------|-----------|-------|
| `cuda-001` | Device enumeration returns ≥1 device on GPU host | `nvidia-smi -L` count | exact | 0 |
| `cuda-002` | `hello.cu` stencil: `out[i] = in[i] + 1.0f` | analytical | bit-exact | 0 |
| `cuda-003` | NVRTC compile of intentionally-broken kernel | returns `Error::Driver` with parse-error message | — | 0 |
| `cuda-004` | cuSOLVER `Zgetrf` on a 256×256 random Hermitian | host `faer` LU | ‖A−LU‖ / ‖A‖ ≤ 1e-12 | 1 |
| `cuda-005` | cuBLAS `Zgemm` GFLOPS at n=4096 | published GPU peak | within 30% | 1 |
| `cuda-006` | FDTD update-kernel throughput | openEMS on same hardware | beats it | 2 |
| `cuda-007` | NCCL ring all-reduce across 2 GPUs | analytical | bit-exact | 2 |

## Running

```bash
# Phase 0
cargo test -p yee-cuda --features cuda

# Phase 1+ (requires GPU)
cargo test -p yee-cuda --features cuda --release -- --include-ignored
```

## CI

- PR builds: `--no-default-features` only (must build clean on every host).
- Nightly self-hosted GPU runner: `--features cuda` end-to-end including benches.
