# FS.0b.2a — graded two-port board fixture (design)

**Date:** 2026-07-12 · **Track:** FS.0 (`FULL-SUITE-ROADMAP.md`) · **ADR:** 0210 addendum

## Problem

The FS.0b.1 gate (`engine_graded_notch.rs`) hand-rolls ~120 lines of
fixture: grid derivation, shared-grid voxelization, coarse-run probe
placement, Courant dt, JobSpec assembly. Every future graded measurement
(filter verify, antenna verify, the FS.0b.2 converge integration) would
copy it. Extract it into `yee_engine::board`.

## Design

- `GradedBoardOptions { mesh: GradedMeshOptions, f0_hz, bw_hz, z0_ohm,
  spacing_cells, n_steps: Option<usize>, backend }` — `n_steps: None`
  applies the notch-gate rule (`9000 · 0.3 mm / fine`, the FS.0a physical
  window at the fine spacing).
- **`two_port_board_jobs_graded(dut, f_max_hz, &opts) ->
  Result<(GradedTwoPortBoardJob, GradedTwoPortBoardJob), String>`** —
  returns the **(DUT, reference)** pair in one call, both voxelized on the
  DUT-derived grid: the ADR-0204 same-physical-problem lesson lives in the
  API shape, not in caller discipline.
- `GradedTwoPortBoardJob { spec, dt_s, spacing_m, cells }`; probes 0–2 /
  3–5 are the A/B triples on **bit-equal-coarse runs** (fit_standing_wave
  needs equal spacing), same layout as the uniform fixture so
  `sparams::forward_transfer` post-processing is identical.
- The gate test refactors onto the fixture — the existing
  `engine-graded-001` CI release job then exercises it end-to-end with
  unchanged asserts (the physics numbers must NOT move: same grid, same
  probes, same dt ⇒ bit-identical solves).

## Gates

- Instant structural test in `board.rs` tests: fixture on the stub layout
  → DUT/ref dims equal, spacings attached, dt = 0.9× graded Courant,
  probe triples on bit-equal coarse cells, ports at the voxelizer's port
  cells.
- `engine-graded-001` (existing release CI job) re-verified through the
  fixture — same measured numbers (4.900 GHz @ −37.2 dB, ratio 0.190).

## Lane

`crates/yee-engine/**`, `docs/**`. (The GPU half of FS.0b.2 runs as a
parallel worktree track in `crates/yee-compute/**` — ADR-0214.)
