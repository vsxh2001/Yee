# ADR-0006: Pin `cudarc` to `=0.19.x` and route all CUDA access through `yee-cuda::backend`

**Status:** Accepted
**Date:** 2026-05-16
**Deciders:** Yee maintainers

## Context

Yee uses CUDA for the dense complex linear solve in `yee-mom` (cuSOLVER
`Zgetrf` / `Zgetrs`) and is on track to use it for FDTD step kernels
(`yee-fdtd`) and for matrix-element fill on GPU. The Rust binding
landscape for CUDA in 2026 has two production-relevant entries plus a
graveyard:

- **`cudarc`** — the `coreylowman/cudarc` crate. Self-describes as a
  "minimal, ergonomic, safe Rust wrapper for CUDA". Production use by
  the Candle ML framework (Hugging Face) and by `mistral.rs`.
  Comprehensive coverage of the driver API, cuBLAS, cuSOLVER, cuRAND,
  cuFFT, NCCL, and the nvrtc JIT.
- **`cust`** — the Rust-CUDA project's high-level binding. Was
  unmaintained from 2022 through early 2025; a community-led reboot
  started in March 2025 under the `Rust-CUDA` GitHub organisation. As
  of May 2026 it is published only as a git dependency on the reboot
  branch; no crates.io release. Promising but not yet production.
- **`rustacuda`** — abandoned since 2021. Builds on modern CUDA only
  with patches. Out of consideration.

`cudarc` is the only currently shippable choice. The problem is
stability:

- `cudarc`'s own README still labels the crate as **"pre-alpha"** as
  of `0.19`. Despite production use by Candle, the API contract is not
  considered stable by its maintainer.
- Recent history shows real breaking minor-version drift. The
  `0.13 → 0.14` transition reorganised the device-pointer types and
  renamed half a dozen public functions. `0.16 → 0.17` changed how
  `LaunchAsync` is parameterised. Each of these would break Yee
  downstream if it consumed `cudarc` directly throughout the workspace.
- `cudarc` has no SemVer guarantee in the strict sense: a `0.x` crate
  by Cargo convention treats every `0.x → 0.y` step as breaking, and
  `cudarc` exercises that liberty more aggressively than most
  ecosystem crates.

The decision space:

1. **Use `cudarc` directly throughout the workspace.** Every CUDA call
   in `yee-mom`, `yee-fdtd`, and any future GPU consumer imports
   `cudarc::*` and uses its types. Breaking changes in `cudarc`
   ripple through every consumer. Rejected.
2. **Use `cudarc` with a thin wrapper around just the device-pointer
   and stream types.** Most call sites stay close to the raw `cudarc`
   API but go through `yee_cuda::Device` instead of
   `cudarc::driver::CudaDevice`. Helps a little; does not insulate
   from the kernel-launch API churn. Rejected.
3. **Introduce a `Backend` trait in `yee-cuda` that abstracts the
   complete CUDA surface Yee uses.** All non-`yee-cuda` workspace
   crates import `yee_cuda::backend::Backend` (and concrete types like
   `BackendDevice`, `BackendStream`, `BackendBuffer`) and never touch
   `cudarc` directly. Migration to a future `cudarc 0.20`, to the
   rebooted `cust`, or to any other binding is then a single-PR
   change to the trait's implementation file. Chosen.
4. **Vendor `cudarc` into the workspace and freeze it.** Removes the
   external dependency entirely, but means we own the maintenance
   burden of every CUDA driver-API addition forever. Rejected.

The `Backend`-trait approach is the conventional pattern when a
project depends on a strategically important but volatile upstream
library. The equivalent pattern shows up in Candle's own
`candle-core::Device`, in Rust web frameworks abstracting over runtime
choice, and in compiler toolchains hiding LLVM C++ ABI churn.

## Decision

Pin `cudarc` to the **exact `0.19.x`** patch series and route every
CUDA access in the workspace through `yee-cuda::backend`.

**Cargo pin** in workspace `Cargo.toml`:

```toml
[workspace.dependencies]
cudarc = { version = "=0.19", default-features = false, features = [
    "driver", "cublas", "cusolver", "cufft", "nvrtc",
] }
```

The `=0.19` syntax means "exactly 0.19.x where the patch is unspecified
but the minor is fixed". Patch updates within `0.19` are accepted via
`cargo update`; a move to `0.20` is gated behind an ADR update.

**Architecture** in `yee-cuda::backend`:

```rust
/// All CUDA functionality Yee needs is exposed through this trait.
/// The concrete impl lives in `yee_cuda::backend::cudarc_impl` and
/// is the *only* file in the workspace that imports `cudarc`.
pub trait Backend {
    type Device: BackendDevice;
    type Stream: BackendStream;
    type DBuf<T: DeviceCopy>: BackendBuffer<T>;
    type SolverCtx: BackendSolverCtx;
    type BlasCtx: BackendBlasCtx;

    fn device(ordinal: u32) -> Result<Self::Device>;
}
```

The CUDA Solver methods (`zgetrf`, `zgetrs`) and the BLAS methods
(`zgemm`, `zgemv`) live on the `SolverCtx` and `BlasCtx` types
respectively, with signatures that take Yee-owned slice types and
never expose a `cudarc::*` type in the public API.

Every workspace crate that needs CUDA imports `yee_cuda::backend::*`,
never `cudarc::*`. Enforced by a `clippy.toml` lint:

```toml
disallowed-types = [
    { path = "cudarc::*", reason = "Use yee_cuda::backend instead. See ADR-0006." },
]
```

(applied only to crates other than `yee-cuda` itself).

## Consequences

**What becomes easier:**

- Migration risk is **bounded to one file** —
  `yee-cuda::backend::cudarc_impl`. If `cudarc 0.20` lands tomorrow
  with a breaking change to `CudaDevice`, the PR to absorb the change
  touches that one file. Nothing in `yee-mom`, `yee-fdtd`, or the rest
  of the workspace changes.
- A future swap to a different CUDA binding (rebooted `cust`,
  hypothetical `vulkano-cuda`, an internal vendored copy) is also
  local to that one file.
- The `Backend` trait gives us a single, documentable, testable
  surface for "everything Yee needs from a GPU". Mocking the backend
  for CI tests on machines without a GPU is straightforward — we
  already do this with a `cpu_fallback_impl` that runs the same
  surface on `nalgebra` and OpenBLAS.

**What becomes harder:**

- A small amount of indirection at every CUDA call site. Instead of
  `let stream = device.fork_default_stream()?;` we write
  `let stream = device.stream(StreamKind::Default)?;` against the
  trait. The diff is syntactic; the cost is paid once per call site
  and never again.
- Adding a new CUDA capability to `yee-mom` or `yee-fdtd` requires
  threading the new method through the `Backend` trait first, then
  implementing it on `cudarc_impl`. This is a forcing function that
  keeps the abstraction honest — no escape hatch into raw `cudarc`
  from the consumer crates.
- Patch updates within `0.19` still require code review attention,
  because `cudarc` has occasionally landed breaking changes inside
  patch releases. The `=0.19` pin reduces but does not eliminate this.

**What's now closed off:**

- Using `cudarc` types directly in any public API of any crate other
  than `yee-cuda`. Enforced by the clippy lint above.
- Skipping the `Backend` trait for a one-off "just need a quick CUDA
  call here" use case. Even one-off uses go through the trait.
- Floating up to `cudarc 0.20` (or wherever) without an ADR update.
  The pin is exact; loosening it is a deliberate act.

## References

- `cudarc` crate, <https://crates.io/crates/cudarc>; README's
  "pre-alpha" disclaimer.
- `cudarc` GitHub: <https://github.com/coreylowman/cudarc>.
- Rust-CUDA `cust` reboot (March 2025):
  <https://github.com/Rust-GPU/Rust-CUDA>.
- `candle-core` use of `cudarc` as a reference production deployment:
  <https://github.com/huggingface/candle>.
- `yee-cuda/src/backend.rs` — the `Backend` trait definition.
- `yee-cuda/src/backend/cudarc_impl.rs` — the single file in the
  workspace that imports `cudarc`.
- `clippy.toml` — `disallowed-types` enforcement.
- ADR-0002 — Rust MSRV 1.88, which sets the floor below which
  `cudarc 0.19` itself would not build.
