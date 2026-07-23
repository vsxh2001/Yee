# FS.7.0 — GPU performance walking skeleton: reproducible bench + first optimization passes

**Date:** 2026-07-23 · **Track:** FS.7 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-compute/**` (+ docs row update)

## Context

FS.7 ("Performance leadership") was queued pending GPU hardware. This machine has an
RTX 5060 Ti 16 GB (448 GB/s GDDR7, Vulkan). A scratch benchmark (2026-07-23, session
T4) measured the wgpu backend at 2.70 Gcells/s @ 64³, 2.48 @ 128³, 1.43 @ 192³ —
against the roadmap's published bar of **3405 Mcells/s (gprMax on Pascal)**. Roofline:
FDTD is memory-bandwidth-bound; at ~130–160 B/cell·step the backend already runs at
~75–85 % of this card's 448 GB/s. Beating the bar needs either less traffic per
cell·step or the honest conclusion that this card's bus can't carry it (the bar was
set on a 732 GB/s HBM2 part) — **a documented NO-GO is an acceptable outcome; a faked
number is not.**

Two measured defects feed this increment:
- `GpuFdtd::step_n` is submit-only (async). Naive timing reads submission overhead
  (~900× "speedup" artifact). There is no public device-wait; benchmarking required a
  `read_fields()` sync hack.
- 192³ throughput drops to 1.43 Gcells/s (vs 2.48 at 128³) — unexplained; suspects are
  `STEPS_PER_SUBMIT` chunking, per-pass bind-group overhead, or the >128 MiB-binding
  path (limits fix `78cd12f` unlocked the size, dip remains).

## Deliverables

1. **`GpuFdtd::sync()`** — public blocking device-wait (poll until queue idle).
   Removes the read-back-to-sync hack; documented as the benchmark/sequencing seam.
2. **`crates/yee-compute/examples/bench.rs`** — the reproducible benchmark the README
   claim will cite: vacuum grid sweep {64³, 96³, 128³, 160³, 192³, 224³}, warm-up,
   sync-correct timing (submit + `sync()`), ≥3 reps w/ median, prints a markdown table
   (grid, Mcells·step/s CPU & GPU, speedup, est. GB/s) and `--json` for CI artifacts.
   Includes a 64³ CPU-vs-GPU rel-L2 sanity gate (≤1e-3, gpu_cpu_parity idiom) so the
   timed physics is verified physics.
3. **Workgroup-shape tuning** — current kernels use `@workgroup_size(4,4,4)`;
   x-contiguous memory favors flat-x shapes. Measure ≥3 shapes, keep the winner.
   Bit-exactness constraint: one thread = one cell, per-cell arithmetic unchanged ⇒
   compute-018/020 uniform bit-exact gates must stay green (verify, don't assume).
4. **192³ dip root-cause** — instrument, isolate (chunk size × grid size matrix),
   fix if software, document if hardware (e.g. cache-tiling cliff).
5. **Docs** — ADR-0223 (numbers, lessons, honest bar verdict), FS.7 row update.

## Non-goals

FP64 GPU, CUDA backend, kernel fusion beyond what the dip fix needs, multi-GPU,
CI perf-regression gates on hosted runners (no GPU there; numbers are artifacts from
the GPU nightly only). E/H pass fusion and a uniform-coeff specialized fast-path are
FS.7.1 candidates if the bar remains unmet.

## Validation gates

- `compute-002/018/020` and all existing yee-compute gates stay green on real hardware.
- New `bench_sanity` (inside bench example): 64³ rel L2 ≤ 1e-3 CPU vs GPU.
- Perf numbers are **reported, not gated** in per-PR CI (machine-dependent); the bench
  binary itself must build in the no-GPU CI (self-skip on NoAdapter, existing idiom).
