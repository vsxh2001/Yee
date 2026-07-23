# FS.7.1 — GPU kernel optimization: index remap + pass fusion

**Date:** 2026-07-23 · **Track:** FS.7 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-compute/**` (+ docs)
**Predecessor:** FS.7.0 (ADR-0223, merge `ec03634`) — bench shipped; bar (3405 Mcells/s) NOT MET:
peak 2864 @ 96³, ~91.9 % estimated bus utilization, decline at ≥128³ root-caused to cache/roofline.

## The two queued levers (from ADR-0223)

1. **gid↔linearization index remap.** The shader's arrays are k-fastest (C order over nz), but
   `gid.x` maps to `i` — the *slowest* axis. Threads within a workgroup enumerate x fastest, so
   adjacent threads today touch memory `ny·nz` elements apart: systematically uncoalesced. This is
   the measured reason flat-x shapes lost 4–5×. Remap so `gid.x → k` (adjacent threads → adjacent
   memory), swap the dispatch-dimension math accordingly, then re-tune the workgroup shape — after
   the remap, flat-x shapes {(64,1,1), (32,2,2), (16,4,4)} vs (4,4,4) must be re-measured from
   scratch (FS.7.0's shape table is invalidated by the remap; do not reuse its conclusion).
2. **Pass fusion 6→2.** A step currently runs 6 volume dispatches (hx,hy,hz then ex,ey,ez). At
   ≥128³ the working set exceeds on-chip cache, so each dispatch re-streams overlapping data from
   DRAM (e.g. every H component is read by two separate E kernels at different offsets). Fuse the
   3 H updates into one kernel (each thread updates hx,hy,hz of its cell) and the 3 E updates into
   one — per-component arithmetic and FP32 operation order UNCHANGED, only co-located in one
   thread. Component independence makes this safe: H components read only E; E components read
   only H + their own array. CPML ψ updates ride along unchanged (per-component ψ arrays).

Expected effect is largest exactly where FS.7.0 measured the decline (128³–224³ DRAM-bound
regime). Either lever may also move the 96³ peak past the 3405 bar — or not; honest verdict again.

## Deliverables

0. Hygiene (pre-existing, surfaced in FS.7.0): `cargo clippy -p yee-compute --all-targets
   --no-default-features -- -D warnings` currently fails with 4 dead-code errors
   (`dispersive.rs::{AdeCoeffs, ade_coeffs}`, `drive.rs::{arena_offset, is_empty}` — consumed
   only from gpu-feature code). Fix by cfg-gating (`#[cfg_attr(not(feature = "gpu"), allow(dead_code))]`
   or `#[cfg(feature = "gpu")]` as appropriate — smallest correct gate wins).
1. Index remap + post-remap shape re-tune (keep the measured winner hardcoded).
2. H-fusion + E-fusion (6→2 volume dispatches; mask clamps / drive / DFT lanes untouched).
3. Re-run the FS.7.0 bench; ADR-0224 with full before/after tables + honest bar verdict;
   FS.7 roadmap row update.

## Hard constraints

- **Bit-exactness is the contract.** compute-018/020 (uniform bit-exact) and
  `graded_uniform_bitexact` must pass UNMODIFIED after both levers. Fusion and remap change
  which thread computes what, never what is computed per cell (same FP32 ops, same order).
  If a lever cannot preserve bit-exactness, that lever stops and the reason is documented —
  the gates are not negotiable.
- Full `cargo test -p yee-compute --release -- --include-ignored` green on the real adapter
  after each lever (all GPU output must show `adapter 'NVIDIA GeForce RTX 5060 Ti'`).
- Each lever measured independently (bench before/after per commit); a lever that measures
  flat or negative is reverted and recorded — negative results are deliverables.
- Perf numbers honest; the bar may still be unreachable on a 448 GB/s bus.

## Non-goals

Uniform-coefficient specialized fast path (vacuum-only benchmark gaming risk — revisit only
with a real-workload justification); FP64; CUDA; multi-GPU; shared-memory tiling (only if
remap+fusion both land and the bar is within ~10 %, else FS.7.2).
