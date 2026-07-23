# FS.7.1 implementation plan — index remap + pass fusion

Spec: `docs/superpowers/specs/2026-07-23-fs7-1-kernel-optimization-design.md`
Branch: `feature/fs7.1-kernel-opt` (from `main` @ `ec03634`+)
Lane: `crates/yee-compute/**`, plus Task 3's doc files (`docs/src/decisions/0224-*.md`,
`docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`).

## Global constraints (bind every task)

- Real GPU present (RTX 5060 Ti); GPU evidence must print
  `adapter 'NVIDIA GeForce RTX 5060 Ti'` — a SKIPPED/NoAdapter run is a task failure.
- **Bit-exact gates pass UNMODIFIED after every commit**: `cargo test -p yee-compute --release
  --test graded_uniform_bitexact --test gpu_graded_parity -- --include-ignored` (compute-018,
  compute-020 "bit-for-bit PASS", compute-021). Plus `gpu_cpu_parity` (compute-002).
  Never weaken any assertion/tolerance anywhere.
- Measurement instrument: `cargo run -p yee-compute --release --example bench` (FS.7.0;
  sync-correct, 3-rep medians). Record full tables in the task report; before/after per lever.
- A lever measuring flat/negative within noise → revert the functional change, keep the
  measured record (comment + report), like FS.7.0 Task 2 did for workgroup shapes.
- `cargo clippy -p yee-compute --all-targets -- -D warnings`, same with
  `--no-default-features` (green after Task 0), `cargo fmt --check --all`, before each commit.
- Commit style `yee-compute: <subject>` ≤72 chars; body why. End body with exactly:
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FFo41x449XDGJ7Xyds4L7M

## Task 0 — no-default-features clippy hygiene

- Fix the 4 pre-existing dead-code errors under
  `cargo clippy -p yee-compute --all-targets --no-default-features -- -D warnings`:
  `src/dispersive.rs` (`AdeCoeffs`, `ade_coeffs`), `src/drive.rs` (`arena_offset`, `is_empty`).
  They are consumed only from `gpu`-feature code. Smallest correct gate:
  `#[cfg_attr(not(feature = "gpu"), allow(dead_code))]` on the items (keeps them compiled +
  type-checked in the no-gpu build) — unless inspection shows a cleaner existing idiom in the
  crate; match house style.
- Verify: both clippy invocations green; `cargo test -p yee-compute --no-default-features`
  and default `cargo test -p yee-compute` green. One commit.

## Task 1 — gid↔linearization remap + shape re-tune

- `crates/yee-compute/src/shaders/fdtd.wgsl`: all 9 volume kernels derive (i,j,k) from
  `global_invocation_id` with gid.x→i today. Remap so **gid.x→k, gid.y→j, gid.z→i**
  (adjacent threads → adjacent k = adjacent memory in the k-fastest linearization). Update
  `gpu.rs` dispatch math (`groups_x/y/z` closures: x now covers nz, z covers nx) and the
  bounds checks. Per-cell arithmetic must not change — this is a pure re-indexing.
- Immediately verify bit-exactness (the constraint suite above) BEFORE measuring — a remap
  bug shows up there first.
- Re-tune workgroup shape post-remap: measure {(4,4,4) as control, (64,1,1), (32,2,2),
  (16,4,4)} at 128³ and 192³ with the bench, 3-rep medians. FS.7.0's shape table is
  invalidated by the remap — do not reuse its conclusions. Keep the winner hardcoded;
  update the WORKGROUP_X/Y/Z comment block (numbers + why, replacing the stale table's
  conclusions with a pointer to both ADRs).
- Full sweep (all 6 grids) before/after in the report. One commit (or two: remap, then tune).

## Task 2 — pass fusion 6→2

- Fuse hx/hy/hz update kernels into one `update_h` kernel (each thread computes all three
  components for its cell, including their CPML ψ contributions), and ex/ey/ez into one
  `update_e`. **Preserve per-component FP32 arithmetic and operation order exactly** — copy
  the existing component bodies verbatim into the fused kernel (shared index math may be
  factored; float expressions may not be reassociated). Component safety: H reads only E;
  E reads only H + its own component's array; ψ arrays are per-component.
- `gpu.rs`: `update_pipelines` [6] → [2]; dispatch order H then E unchanged; mask clamp /
  drive / DFT pipelines untouched. Delete the six old kernels (no dead shader code).
- Verify bit-exactness suite FIRST, then full `cargo test -p yee-compute --release --
  --include-ignored`, then measure (full 6-grid sweep before/after).
- Watch register pressure: if the fused kernel measures SLOWER (occupancy loss can beat
  traffic savings), revert per the global revert rule and record — that is a valid outcome.
- One commit.

## Task 3 — bench, ADR-0224, roadmap row

- Definitive full-sweep bench run at final HEAD (3 reps; note ambient GPU state via
  nvidia-smi before running).
- `docs/src/decisions/0224-fs71-kernel-opt.md` (follow 0223's structure): per-lever
  before/after tables, the post-remap shape table, fused-kernel outcome incl. any revert,
  final peak vs the 3405 Mcells/s bar with the honest verdict, estimated GB/s vs 448,
  what FS.7.2 would need (shared-memory tiling, or hardware). Add the SUMMARY.md line
  after ADR-0223. Update the FS.7 row in `FULL-SUITE-ROADMAP.md`.
- Lane exception: exactly these doc files.
