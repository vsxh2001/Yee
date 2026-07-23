# ADR-0224: FS.7.1 — GPU kernel optimization: index remap + pass fusion (bar MET)

**Date:** 2026-07-23 · **Status:** accepted · **Track:** FS.7 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-23-fs7-1-kernel-optimization-design.md`
**Plan:** `docs/superpowers/plans/2026-07-23-fs7-1-kernel-optimization.md`
**Predecessor:** FS.7.0 (ADR-0223, merge `ec03634`) — bench shipped; bar (3405 Mcells/s) NOT MET:
peak 2864 Mcells/s @ 96³ (84.1 % of bar), ~91.9 % of the RTX 5060 Ti's 448 GB/s bus at that grid
size. Two levers queued: `gid.x → k` index remap, and E/H pass fusion.

## Context

ADR-0223 root-caused the FS.7.0 shortfall to two things: (1) every kernel entry point mapped
`gid.x → i`, the *slowest*-varying array axis in the k-fastest linearization, so adjacent threads
in a workgroup touched memory `ny·nz` elements apart — systematically uncoalesced; and (2) a step
dispatches 6 separate volume kernels (hx, hy, hz, ex, ey, ez), each re-streaming overlapping field
data from VRAM, with the working set outgrowing on-chip cache from 128³ onward (monotonic decline
2411→934 Mcells/s, 96³→224³). Both were named as the FS.7.1 candidates and are exactly what this
track executed, plus a pre-existing no-default-features clippy hygiene fix (Task 0, not
perf-relevant, see `fs71-task-1-report.md` for confirmation it was already landed on this branch).

## Decision

### 1. `gid.x → k` index remap (Task 1, commit `d4370c9`)

All 9 volume kernels (`update_hx/hy/hz`, `update_ex/ey/ez`, `clamp_ex/ey/ez`) changed from
`let i = gid.x; let j = gid.y; let k = gid.z;` to `let k = gid.x; let j = gid.y; let i = gid.z;`
— pure re-indexing, no per-cell arithmetic touched. `gpu.rs` dispatch math updated to match
(`groups_x` now covers `nz`, `groups_z` covers `nx`). Verified bit-exact immediately, before any
measurement, per the hard constraint.

**Measured (median of 3×200-step reps, `(4,4,4)` control shape, before vs after remap only):**

| grid  | pre-remap (ADR-0223) | post-remap, `(4,4,4)` | Δ |
|-------|----------------------|------------------------|---|
| 64³   | 2728.8 | 3574.1 | +31.0 % |
| 96³   | 2864.1 | 3976.9 | +38.9 % |
| 128³  | 2411.2 | 3165.9 | +31.3 % |
| 160³  | 1841.3 | 3004.6 | +63.2 % |
| 192³  | 1401.8 | 2997.7 | +113.9 % |
| 224³  | 934.2  | 2984.1 | +219.4 % |

The remap alone erased essentially all of the 128³→224³ decline ADR-0223 root-caused to
cache/roofline — throughput is now flat ~2980–3170 Mc/s across that whole range instead of
falling 2410→934.

### 2. Post-remap workgroup shape re-tune (Task 1, commit `ea2888c`)

FS.7.0's shape table (which found `(4,4,4)` beating all flat-x candidates) is invalidated by the
remap — re-measured `{(4,4,4)` control, `(64,1,1)`, `(32,2,2)`, `(16,4,4)}` from scratch across
all 6 grids, 2 runs each for reproducibility:

| grid  | `(4,4,4)` control | `(64,1,1)` | `(32,2,2)` | `(16,4,4)` |
|-------|--------------------|------------|------------|------------|
| 64³   | 3574.1 / 3549.8 | 6943.4 | 7566.5 / 7569.1 | 6722.8 |
| 96³   | 3976.9 / 3975.3 | 8641.7 | 9408.0 / 9379.4 | 8097.5 |
| 128³  | 3165.9 / 3167.5 | 3460.6 | 3449.1 / 3448.3 | 3469.5 |
| 160³  | 3004.6 / 3003.3 | 3199.9 | 3200.4 / 3201.9 | 3197.4 |
| 192³  | 2997.7 / 2931.1 | 3142.5 | 3143.6 / 3143.2 | 3152.1 |
| 224³  | 2984.1 / 2984.1 | 3153.8 | 3148.9 / 3148.8 | 3140.0 |

All three flat-ish shapes are within noise of each other at 128³–224³ (the previous 4–5× flat-x
penalty is gone now that `gid.x` tracks the fast axis) but diverge sharply at 64³/96³, where
`(32,2,2)` clearly wins and reproduces. `(32,2,2)` kept as the hardcoded shape (no runtime knob,
per YAGNI) — it is the only shape at-or-near the top across the whole sweep. Comment block above
`WORKGROUP_X/Y/Z` in `gpu.rs` rewritten with these numbers.

### 3. E/H pass fusion 6→2 (Task 2, commit `1051c5f`)

Fused the 3 H-update kernels into one `update_h` entry point (each thread computes hx, hy, hz for
its cell) and the 3 E-update kernels into one `update_e`, dispatched over the union extent
`cell_dims() = (nx+1, ny+1, nz+1)` (added to `FdtdSpec`). Per-component bodies kept as separate
WGSL functions (`do_update_hx/hy/hz`, `do_update_ex/ey/ez`) called unconditionally from the fused
entry point — necessary because the dispersive branch's `return` would otherwise abort sibling
components, not just itself; every float expression inside each function is byte-identical to the
pre-fusion body (no reassociation). `update_pipelines: [ComputePipeline; 6] → [2]`; dispatch order
(H, then E) unchanged; clamp/drive/DFT lanes untouched. **No revert needed — positive at every
grid size measured:**

| grid  | before (remap+tune, `(32,2,2)`) | after (fused) | Δ |
|-------|----------------------------------|----------------|---|
| 64³   | 7575.4 | 11672.2 | +54.1 % |
| 96³   | 9384.5 | 12784.7 | +36.2 % |
| 128³  | 3444.9 | 4444.2  | +29.0 % |
| 160³  | 3200.6 | 4519.0  | +41.2 % |
| 192³  | 3142.1 | 4686.0  | +49.1 % |
| 224³  | 3149.3 | 4686.9  | +48.8 % |

Largest in the 128³–224³ DRAM-bound regime the spec predicted (+29–49 %), but also strongly
positive at the small cache-resident grids (+36–54 %) — fusion cut dispatch/traffic overhead
everywhere, not just where the working set exceeds cache.

### Bit-exactness (unmodified gates, every commit)

`compute-018`/`graded_uniform_bitexact` (bit-exact PEC/CPML), `compute-020`
(`gpu_graded_parity`, "bit-for-bit PASS"), `compute-021` (graded-taper reflection
`-52.68 dB`, unchanged across all three commits — expected, since neither lever changes
per-cell arithmetic), `compute-002` (`gpu_cpu_parity`, FP32 tolerance) — all green after
every commit in this track, confirmed again at final HEAD below. Full
`cargo test -p yee-compute --release -- --include-ignored` green throughout (22 test-result
lines, 0 failed), exercising CPML, heterogeneous materials, per-face CPML, dispersive ADE
(the branch with the `return`-inside-fusion restructuring), aperture ports, NTFF, and the
driven microstrip line. No gate file was ever modified — only `gpu.rs`, `shaders/fdtd.wgsl`,
and `spec.rs` changed across Tasks 0–2.

## Definitive measured results (final HEAD `1051c5f`, real adapter)

Ambient GPU state before the definitive sweep: idle (`14 MiB used / 16311 MiB, 0 % util, 32-39
°C, 180 MHz`). Two earlier attempts at this final sweep (not used as evidence) caught transient
contention from an unrelated background process (`nvidia-smi --query-compute-apps` showed a
`python` process at up to 4.3 GiB / 67 % util mid-sweep, `readback ms` spiking to 915 ms at 160³
and GPU Mc/s dropping to ~2670–2780 at 128³–192³ in that run) — a real, documented example of
measurement noise from GPU sharing, not a regression; discarded and re-run once the GPU returned
to idle. The table below is 3 reps run back-to-back on an idle GPU (medians used where they
diverge):

```
$ cargo run -p yee-compute --release --example bench
adapter: NVIDIA GeForce RTX 5060 Ti
| grid     |      CPU Mc/s |      GPU Mc/s |  speedup |     GB/s |  readback ms |
|----------|---------------|---------------|----------|----------|--------------|
| 64^3     |         886.6 |       11672.7 |   13.17x |     1681 |         6.85 |
| 96^3     |         501.5 |       12776.2 |   25.47x |     1840 |        21.81 |
| 128^3    |         251.0 |        4442.0 |   17.70x |      640 |        45.95 |
| 160^3    |         236.4 |        4518.4 |   19.11x |      651 |        91.22 |
| 192^3    |         239.4 |        4685.9 |   19.58x |      675 |       156.96 |
| 224^3    |         237.8 |        4682.5 |   19.69x |      674 |       246.49 |
sanity (64^3 per-E-component rel L2 vs CPU): PASS (worst rel L2 = 9.745e-7 <= 1e-3)
```

| grid  | run1 (above) | run2  | run3  | median | vs FS.7.0 baseline | vs 3405 bar |
|-------|--------------|-------|-------|--------|---------------------|-------------|
| 64³   | 11672.7 | 11535.4 | 11536.6 | 11536.6 | ×4.23 | 338.8 % |
| 96³   | 12776.2 | 12665.0 | 12598.9 | 12665.0 | ×4.42 | **372.0 %** |
| 128³  | 4442.0  | 4442.3  | 4444.0  | 4442.3  | ×1.84 | 130.5 % |
| 160³  | 4518.4  | 4521.2  | 4519.0  | 4519.0  | ×2.45 | 132.7 % |
| 192³  | 4685.9  | 4190.2  | 4685.1  | 4685.1  | ×3.34 | 137.6 % |
| 224³  | 4682.5  | 4687.2  | 4685.0  | 4685.0  | ×5.02 | 137.6 % |

(run2's 192³ = 4190.2 is a mild single-run dip, not reproduced in run1/run3; still well clear of
the bar and not treated as a regression.) All three runs' sanity gate (`64³` per-E-component
rel L2 vs CPU) `PASS`, worst rel L2 ≈ 9.7e-7, on every run. These numbers match Task 2's
post-fusion table to within run-to-run noise (as expected — Task 3 made no functional change),
confirming reproducibility at final HEAD.

## Verdict vs the 3405 Mcells/s bar (gprMax on Pascal, 732 GB/s HBM2)

- **Peak measured throughput: 12665 Mcells/s @ 96³ — 372.0 % of the 3405 Mcells/s bar. Bar MET,
  by a wide margin, reversing FS.7.0's NO-GO.**
- **The bar is now met at every measured grid size, not just the peak**: the worst point in the
  sweep (128³, the smallest grid past the cache-resident regime) still clears the bar at 130.5 %.
  FS.7.0's decline from 96³→224³ (2864→934 Mc/s, falling below the bar past ~100³) is gone; the
  post-remap+fusion curve is flat-to-rising from 128³ onward (4442→4685 Mc/s).
- **Estimated GB/s vs the 448 GB/s nameplate bus**: the bench's traffic model
  (`BYTES_PER_CELL_STEP = 144.0`, a fixed per-cell-step working-set estimate, unchanged by either
  lever) now reads **above the 448 GB/s nameplate at every grid size** — 1681–1840 GB/s at
  64³/96³ (as already flagged by Task 1: those grids' whole working set is cache-resident across
  the 200-step timed window) and, new in this track, **640–675 GB/s at 128³–224³ too** (vs
  452–500 GB/s at those same sizes after Task 1's remap alone, and vs 135–347 GB/s pre-remap in
  ADR-0223). The fixed 144 B/cell/step assumption was calibrated against the pre-fusion,
  6-dispatch traffic pattern; fusion's whole point is that co-located H/E components reuse data
  in registers/L1/L2 across what used to be 3 separate re-streams from DRAM, so the constant now
  systematically **overstates** true external bus traffic at every grid size, not just the
  small ones. This is an honest reading, not a claim the card exceeds its own spec: the model's
  "GB/s" column stopped being a good proxy for external bandwidth once the kernel stopped
  matching the streaming assumption it was calibrated against.
- **Honest conclusion: GO.** Both queued FS.7.0 levers (index remap, pass fusion) closed the gap
  entirely on this hardware — no bigger-bus card was needed. The remap fixed uncoalesced access
  (the mechanical cause of the workgroup-shape flatness FS.7.0 found); fusion cut redundant DRAM
  round-trips for the H/E field data (the mechanical cause of the 128³+ decline). Neither lever
  needed a revert; both measured positive at every one of the 6 grid sizes.

## What FS.7.2 would need

The bar is met with margin everywhere measured, so FS.7.2 is not required to publish a passing
number. If pursued anyway (diminishing-returns territory), the remaining internal shape to chase
is the **96³→128³ cliff**: throughput drops from ~12665 to ~4442 Mc/s (a 2.85× step) exactly
where the per-step working set stops fitting in on-chip cache — a genuine roofline effect, not a
dispatch defect (same diagnosis pattern as ADR-0223, now one tier larger). Closing that step
without new hardware needs **explicit shared-memory (workgroup-local) tiling**: stage each
workgroup's halo region in `workgroup`-address-space memory once, reuse it across the H and E
components computed by threads in that tile, and cut the per-cell DRAM traffic below what the
current per-thread-independent-fetch fused kernel achieves. This is exactly the item the FS.7.1
spec queued conditionally ("only if remap+fusion both land and the bar is within ~10 %, else
FS.7.2") — both landed and the bar is not just within 10 %, it is cleared outright, so tiling is
now optional upside rather than a gap-closer. The uniform-coefficient fast path named in FS.7.0's
non-goals remains explicitly out of scope (benchmark-gaming risk on a vacuum-only workload,
spec's own words) unless a real non-benchmark workload justifies it. A bigger-bus GPU remains a
non-software lever that would trivially extend headroom further but is not needed to clear 3405.

## Non-goals / queued

Unchanged from FS.7.0/FS.7.1's spec: FP64 GPU, CUDA backend, multi-GPU, CI perf-regression gates
on hosted runners (no GPU there — numbers stay GPU-nightly artifacts). Shared-memory tiling and
the uniform-coefficient fast path are now purely optional upside (see above), not gap-closers.
