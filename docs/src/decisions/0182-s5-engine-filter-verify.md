# ADR-0182: S.5 engine-powered verify — materials on the job protocol

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0181 (defined this convergence path), ADR-0179/0180 (S.0 protocol +
`yee-server`), ADR-0176 (E.1 materials on `yee-compute`), ADR-0108 (the
`fdtd-line-eeff-001` scenario and its Hammerstad–Jensen reference).
**Spec:** `docs/superpowers/specs/2026-07-06-s5-engine-filter-verify-design.md`

## Context

ADR-0181 closed S.4 with a convergence path instead of a UI port: the filter flow submits
**full-wave jobs over the S.0 protocol**. The blocker was that `yee_engine::JobSpec` only
described uniform-vacuum jobs, while a voxelized layout (`yee_voxel::voxelize_microstrip`)
needs per-cell ε_r maps, interior PEC masks (trace + ground), and the voxelizer's own dt —
all supported by `yee-compute` since E.1 but not carried by the protocol.

## Decision

1. **`JobSpec` gains two `#[serde(default)]` fields** (additive — every existing client,
   stored spec, and WS peer keeps working):
   - `materials: Option<MaterialsSpec>` — a serde mirror of `yee_compute::Materials`
     (ε_r/μ_r/σ cell maps in the `[nx+1, ny+1, nz+1]` YeeGrid convention; per-component
     staggered PEC masks), handed unchanged to both the CPU and GPU constructors.
   - `dt_s: Option<f64>` — explicit time-step override, bounded by the Courant limit.
2. **Malformed specs produce `JobEvent::Error`, never a panic.** `yee_compute::
   Materials::validate` panics (fine for a library, fatal for a server: a worker-thread
   panic closes the event channel with no terminal event). The engine now pre-validates
   map/mask lengths, dt bounds, and grid dimensions and reports failures as events.
3. **Voxel → spec conversion stays at the call site** (≈10 lines in the gate test). A
   public `Layout → JobSpec` bridge is deferred until the studio or filter CLI consumes
   it, to avoid coupling `yee-voxel` and `yee-engine` ahead of need.

## Gates

- **engine-verify-001** (`yee-engine/tests/verify_line_eeff.rs`, `#[ignore]`, release CI):
  the `compute-008` / `fdtd-line-eeff-001` scenario — dimensioned FR-4 microstrip line,
  voxelized, 50 Ω resistive-port drive, time-gated two-probe phase velocity → ε_eff —
  expressed **as a `JobSpec`** and run through `submit()`/`JobEvent`/`JobResult`.
  Measured **ε_eff err 0.132 % vs the published Hammerstad–Jensen closed form**
  (≤ 15 % gate) — identical to compute-008, as it must be: the protocol adds no physics.
- **Fast (non-ignored)**: serde round-trip incl. materials/dt + legacy-spec (missing
  keys) deserialization; a heterogeneous job through `submit()` **bit-identical** to a
  direct `CpuFdtd::with_drive` run; four error paths (bad map length, bad mask length,
  Courant-violating dt, zero-cell grid) each yielding an `Error` event.
- CI: the `compute-engine-gates` job now also runs
  `cargo test -p yee-engine --release -- --include-ignored`.

## Consequences

Any S.0/S.1 client — Tauri studio, `yee-server` WebSocket peers, Python — can now run
full-wave jobs on real voxelized geometry. The F1.3 filter verify gate (S-parameters vs
spec mask) can be built directly on this protocol; S-parameter extraction, dispersive
material specs over the wire, and a public layout→spec bridge are the follow-ons.
