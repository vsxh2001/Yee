# FS.0b.2-GPU implementation plan — graded grid on the wgpu backend

**Spec:** `docs/superpowers/specs/2026-07-12-fs0b2-gpu-graded-design.md`.

## Steps

1. **`yee-compute/src/spec.rs`** — a small `pub(crate)` helper on
   `SpacingArrays` that emits the packed **inverse** f32 vector
   (`inv_xp|inv_yp|inv_zp|inv_xd|inv_yd|inv_zd`, inverses computed in f64,
   narrowed once), shared by the GPU build and refresh paths.
2. **`yee-compute/src/shaders/fdtd.wgsl`** — binding 8
   (`var<storage, read> inv_sp`), the six index accessors, and the divisor
   swap in all six update kernels per the spec's mapping table (H = primal,
   E = dual; the fused CPML corrections reuse the same curl variables and
   need no separate edit). Remove `inv_dx/inv_dy/inv_dz` from `Params`.
3. **`yee-compute/src/gpu.rs`** — mirror the `Params` change; create the
   spacing buffer at build with the uniform fill; extend the bind-group
   layout/group with binding 8; stash `npml`/`faces`, the dispersive flag,
   and the `Drive` copy on the struct; implement
   `set_spacings(&GradedSpacings) -> Result<(), ComputeError>`
   (validation panics mirroring `CpuFdtd`, `Unsupported` for NTFF-DFT and
   z-taper-straddling aperture ports, `write_buffer` refreshes for the
   spacing buffer + resistive-port `alpha`/`gamma` + aperture `vcoef`).
4. **Gate `compute-020`** —
   `yee-compute/tests/gpu_graded_parity.rs::gpu_graded_uniform_parity`
   (fast): compute-018 drive scenario, three runs; GPU-graded ≡ GPU-scalar
   bit-for-bit; GPU-graded vs CPU-FP64 within compute-002 tolerances.
   Self-skipping (`NoAdapter` → SKIPPED + green).
5. **Gate `compute-021`** — same file,
   `gpu_graded_taper_parity` (`#[ignore]`, release): compute-019 taper
   scenario graded on both backends; probe-series parity within
   measured-then-pinned FP32 tolerances; finiteness; if runtime allows,
   the GPU-side reflection level vs the −48 dB floor. Measure on llvmpipe
   first, then pin.
6. **Docs** — ADR-0214 (decisions, measured numbers, scope rejections,
   binding-budget note) + `SUMMARY.md` line + `lib.rs` module-doc update
   (the "CPU-only" FS.0b.0 paragraph). Check `gpu-nightly.yml`: its
   existing `--include-ignored` yee-compute step already covers both gates
   — add nothing unless that turns out false.

## Verification

```sh
cargo test -p yee-compute --features gpu \
  && cargo test -p yee-compute --release --features gpu \
       --test gpu_graded_parity -- --include-ignored --nocapture \
  && cargo clippy -p yee-compute --all-targets --features gpu -- -D warnings \
  && cargo fmt --check --all
```
