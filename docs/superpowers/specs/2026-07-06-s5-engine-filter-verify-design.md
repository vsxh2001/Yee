# S.5 — Engine-powered filter verify (walking skeleton): materials on the job protocol

**Date:** 2026-07-06
**Phase:** S.5 (ENGINE-STUDIO-ROADMAP), the convergence path defined by ADR-0181.
**Plan:** `docs/superpowers/plans/2026-07-06-s5-engine-filter-verify.md`

## Problem

ADR-0181 closed S.4 with a defined convergence path between the two studios: the filter
flow submits **full-wave jobs over the S.0 protocol** instead of porting UI. That path is
blocked on one gap: `yee_engine::JobSpec` (S.0) only describes **uniform-vacuum** jobs.
A real layout — the output of `yee_voxel::voxelize_microstrip` — needs per-cell ε_r maps
and interior PEC masks (trace + ground), plus the voxelizer's own `dt`. `yee-compute`
has supported all of that since E.1; the job protocol simply doesn't carry it.

## Design

### Protocol extension (additive, backward-compatible)

`JobSpec` gains two `#[serde(default)]` optional fields, so every existing client
(Tauri studio, `yee-server` WS clients, stored JSON specs) keeps working unchanged:

- `materials: Option<MaterialsSpec>` — a serde mirror of `yee_compute::Materials`:
  `eps_r_cells` / `mu_r_cells` / `sigma_cells` (`[nx+1, ny+1, nz+1]` row-major, the
  YeeGrid convention) and `pec_mask_ex/ey/ez` (each E component's staggered shape).
- `dt_s: Option<f64>` — explicit time-step override (the voxelizer computes its own
  Courant dt; probe series are only meaningful against the dt that actually ran —
  `JobResult.dt_s` already reports it).

**Validation is the engine's job, not the caller's:** wrong-length maps/masks and
non-positive or Courant-violating dt produce a `JobEvent::Error`, never a panic —
`yee_compute::Materials::validate` panics, which is fine for a library API but not for
a spec that arrives over a WebSocket.

Both backends receive the materials unchanged (`CpuFdtd::with_drive` /
`GpuFdtd::with_drive` already take `Materials`); no `yee-compute` changes.

### Voxel → spec conversion stays at the call site

The walking skeleton does **not** add a `Layout → JobSpec` public bridge crate/API.
The gate test converts `VoxelModel` → `MaterialsSpec` inline (≈10 lines), exactly as
`compute-008` converts it to `yee_compute::Materials`. A public bridge (its natural
home: a `yee-voxel` feature or a thin `yee-verify` crate, since `yee-voxel` must not
depend on `yee-engine` by default and vice versa) is deferred until the studio or the
filter CLI actually consumes it — walking-skeleton first.

## Validation gates

- **engine-verify-001** (`crates/yee-engine/tests/verify_line_eeff.rs`, `#[ignore]`,
  release CI): the exact `compute-008` / `fdtd-line-eeff-001` scenario — dimensioned
  FR-4 microstrip line (W = 3 mm, h = 1.6 mm, ε_r = 4.4, L ≈ 6 λ_g at 5 GHz),
  voxelized, 50 Ω resistive-port drive, hard-PEC box, time-gated two-probe
  phase-velocity → ε_eff — but expressed **as a `JobSpec` and run through
  `submit()`/`JobEvent`/`JobResult`**. Assert ≤ 15 % vs the published
  Hammerstad–Jensen closed form (`yee_layout::eps_eff`), the original gate's band.
  This certifies the full protocol chain a filter-verify client would use.
- **Fast, non-ignored** (in-crate): (a) serde round-trip of a spec with materials +
  dt override; (b) **parity**: a small heterogeneous job (ε_r block + PEC mask + dt
  override) through `submit()` is bit-identical to a direct `CpuFdtd::with_drive` run
  of the same scenario; (c) error paths: wrong-length `eps_r_cells`, wrong-length
  mask, and `dt_s ≤ 0` each yield `JobEvent::Error` (no panic, no hang).

## Non-goals

S-parameter extraction / spec-mask comparison on the engine (that is the F1.3 verify
gate's job, and it can now be built on this protocol); a `Layout` field in `JobSpec`
(the protocol stays domain-agnostic); dispersive-material specs over the wire
(engine supports ADE, protocol exposure deferred until a client needs it); studio UI
for material jobs.
