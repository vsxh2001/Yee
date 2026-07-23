# ADR-0223: FS.7.0 — GPU performance walking skeleton: reproducible bench + honest bar verdict

**Date:** 2026-07-23 · **Status:** accepted · **Track:** FS.7 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-23-fs7-0-gpu-bench-design.md`

## Context

FS.7 ("Performance leadership") was queued pending GPU hardware. This machine
has an RTX 5060 Ti 16 GB (448 GB/s GDDR7, Vulkan). A prior scratch benchmark
(session T4, `.superpowers/sdd/task-4-report.md`) found two defects that
blocked honest measurement: `GpuFdtd::step_n` is submit-only (no device-wait,
so naive timing reads ~900× submission-overhead noise, not compute), and 192³
throughput dropped to 1.4 Gcells/s vs 2.4 Gcells/s at 128³ for no documented
reason. Both are cleared up below.

## Decision

1. **`GpuFdtd::sync()`** (`crates/yee-compute/src/gpu.rs`) — public blocking
   device-wait, extracted from the existing readback poll into a shared
   `wait_idle()` used by both `sync()` and `read_f32_buffer` (no behavior
   change to readback). This is the benchmark/sequencing seam; timing
   methodology is warm-up 10 steps + `sync()`, then time `step_n(200)` +
   `sync()`, 3 reps, median — readback excluded from the stepping figure.
2. **`crates/yee-compute/examples/bench.rs`** — the reproducible bench:
   vacuum grid sweep {64,96,128,160,192,224}³, CPU (`CpuFdtd`) vs GPU
   (`GpuFdtd`), Mcells·step/s, GB/s (from a documented 144 B/cell/step
   traffic-model assumption), readback ms, `--json`, self-skip on
   `NoAdapter`, per-size `catch_unwind` guard for alloc/validation failures,
   and a 64³ CPU-vs-GPU per-E-component rel L2 ≤ 1e-3 sanity gate.
3. **Workgroup-shape tuning: kept `(4,4,4)` — a documented negative
   result.** `(64,1,1)`, `(32,2,2)`, `(8,8,4)` were measured against the
   baseline at 128³/192³; all three lost (−80.5 %, −78.3 %, −13.8 % at
   128³ respectively). Root cause (measured, not guessed): `fdtd.wgsl`
   linearizes `(i,j,k)` with `k` fastest-varying, but every kernel entry
   point maps `gid.x → i` (the *slowest* axis) directly. Widening the
   workgroup along `gid.x` therefore widens the per-thread memory stride
   instead of improving coalescing. A shape that actually exploited the
   contiguous `k` axis would need remapping `gid.x → k` in every kernel's
   index arithmetic — out of scope for a workgroup-shape task, flagged as
   the real lever for a future increment. No runtime knob was added
   (YAGNI); the winning shape is the pre-existing hardcoded value.
4. **192³ (and onward) throughput dip: root-caused as a memory roofline
   effect, not fixed** (fixing it means restructuring the arena buffers,
   explicitly out of scope this increment). All four plan-listed
   hypotheses were tested by direct measurement on the RTX 5060 Ti and
   refuted:
   - **`STEPS_PER_SUBMIT` chunking** — swept `{8,16,32,64,128,256}`, flat
     to within ~0.3–0.5 % noise across a 32× range at both 128³ and 192³.
   - **Per-pass bind-group/encoder overhead** — same evidence as above (one
     compute pass covers a whole chunk; an 8× increase in pass-boundary
     frequency changed nothing measurable).
   - **>128 MiB single-binding driver slow-path** — the field arena crosses
     128 MiB between 160³ (93.8 MiB) and 192³ (162.0 MiB); if that boundary
     added a distinct slow path, 160³→192³ should show an *extra* drop
     beyond trend. It doesn't: 128³→160³ is −23.5 %, 160³→192³ is −23.9 %
     — the same trend, no discontinuity at the crossing.
   - **Power/thermal** — `nvidia-smi dmon` during a full sweep: core clock
     2600–2700 MHz (near the 3090 MHz boost cap, not throttled), 34–45 °C,
     peak ~107 W. No throttling.

   What the decline actually is: GPU Mcells/s falls monotonically and
   increasingly steeply from 96³ onward (+4.9 %, then −15.8 %, −23.5 %,
   −23.9 %, −33.3 % step-to-step), tracking the per-step working set
   (field + coefficient arenas) growing from ~10 MiB at 64³ to ~431 MiB at
   224³ — the working set outgrowing on-chip cache/L2 reuse as the grid
   grows. This is a general memory-hierarchy/roofline effect for
   bandwidth-bound stencil codes, not a defect in this crate's dispatch
   code. All investigation numbers are recorded in a source comment above
   `STEPS_PER_SUBMIT` in `gpu.rs` so a future FS.7.1 doesn't redo this work.

## Measured results (definitive table, `2f8bc4f`, real adapter)

```
$ cargo run -p yee-compute --release --example bench
adapter: NVIDIA GeForce RTX 5060 Ti
| grid     |      CPU Mc/s |      GPU Mc/s |  speedup |     GB/s |  readback ms |
|----------|---------------|---------------|----------|----------|--------------|
| 64^3     |         941.9 |        2728.8 |    2.90x |      393 |         5.86 |
| 96^3     |         687.2 |        2864.1 |    4.17x |      412 |        20.16 |
| 128^3    |         336.0 |        2411.2 |    7.18x |      347 |        41.54 |
| 160^3    |         296.8 |        1841.3 |    6.20x |      265 |        85.74 |
| 192^3    |         296.6 |        1401.8 |    4.73x |      202 |       148.04 |
| 224^3    |         294.4 |         934.2 |    3.17x |      135 |       231.27 |
sanity (64^3 per-E-component rel L2 vs CPU): PASS (worst rel L2 = 9.745e-7 <= 1e-3)
```

(Task 1/2/3 sessions independently reproduced this table to within ~0.3–0.5 %
run-to-run noise; this run re-confirms it at the final HEAD, no disagreement
to reconcile.)

## Verdict vs the 3405 Mcells/s bar (gprMax on Pascal, 732 GB/s HBM2)

- **Peak measured throughput: 2864 Mcells/s @ 96³ — 84.1 % of the 3405
  Mcells/s bar. Bar NOT MET on this hardware.**
- **Roofline context: at that same 96³ sweet spot the backend already
  moves ~412 GB/s of this card's 448 GB/s bus (91.9 % utilization)**, under
  the bench's documented 144 B/cell/step working-set-floor traffic model.
  The bar's reference card has 732 GB/s (1.63× this card's bus) — a
  bandwidth-naive rescale of 3405 Mcells/s to 448 GB/s would predict
  ~2082 Mcells/s, which this backend already *exceeds* (2864), meaning the
  per-GB/s efficiency here is not the gap; the gap is that FDTD as currently
  dispatched is memory-bandwidth-bound and this card simply has less bus
  than the reference part.
- **What FS.7.1 would need to close the gap on this hardware**: less
  traffic per cell·step, not more parallelism or a different dispatch
  shape (workgroup tuning and chunking were both measured flat/negative
  above) — specifically the two candidates the spec named as FS.7.1
  material: **E/H pass fusion** (avoid round-tripping the H fields through
  VRAM between the two half-steps) and a **uniform-coefficient fast path**
  (skip the per-cell coefficient-map reads/that arena entirely when the
  whole grid is vacuum/uniform, cutting the ~144 B/cell/step model's
  coefficient-read term). A bigger-bus card (a 732+ GB/s part) would also
  close the gap without any code change, per the roofline math above — but
  that's not a software lever.
- **Honest conclusion: NO-GO on this hardware as measured, not a fake
  number.** The backend is bandwidth-bound and already runs the vacuum FDTD
  kernels at ~92 % of this card's bus at the best-fit grid size; beating a
  bar set on a card with 1.63× more bandwidth needs either less traffic per
  cell·step (FS.7.1) or better hardware, not further dispatch tuning (both
  workgroup shape and submit-chunking were measured and found flat/negative
  in this increment).

## Non-goals / queued

FP64 GPU, CUDA backend, multi-GPU, CI perf-regression gates on hosted
runners (no GPU there — numbers are GPU-nightly artifacts only). FS.7.1
candidates if the bar is to be chased further: E/H pass fusion, a
uniform-coefficient specialized fast path, and — orthogonal to any code
change — a card with more memory bandwidth. The `gid.x → k` kernel
index-remap flagged in the workgroup-tuning investigation (a materially
larger change than shape tuning) is also queued there.
