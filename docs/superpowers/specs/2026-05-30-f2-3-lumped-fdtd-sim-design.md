# Filter Phase F2.3 — lumped-LC FDTD EM simulation — Design Spec

**ADR:** ADR-0115 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

The lumped-LC goal names **EM simulation**. F2.3 = full-wave **FDTD** of the
lumped-LC board: voxelize the F2.2 board, place each L/C as a `LumpedRlcPort`
element on the grid, drive/sense two ports, and extract **S21(f)** — then
**cross-validate** the FDTD response against the analytic circuit `ladder_s21`
(F2.0). This is the lumped analogue of the distributed `fdtd-line-eeff` gate.
Heavy (multi-minute FDTD) → validated in the bounded container, CI-gated.

## Modeling (the key subtlety)

`yee_fdtd::LumpedRlcPort::series_rlc(cell, r, l, c, src)` is a **series** R-L-C at
one `E_z` cell (`l=0` → pure capacitor; `c=∞` → pure inductor). Map the ladder:

- **Series-branch resonator** (series L–C, in the signal path): ONE
  `series_rlc(cell, r≈0⁺, L, C)` at the in-line gap cell of that footprint.
- **Shunt-branch resonator** (parallel L–C, line→ground): TWO elements at the
  shunt cell — `series_rlc(cell, r, L, c=∞)` (pure inductor) **in parallel with**
  `series_rlc(cell, r, l=0, C)` (pure capacitor). Both correct `E_z` at the same
  cell ⇒ they act as a parallel L‖C admittance (correct shunt-resonator topology).
- Use a small finite `r` (e.g. component ESR, or a tiny 1e-3 Ω) since `series_rlc`
  requires `r>0`.

Drive a `LumpedRlcPort` source at the input line end, terminate the output with a
matched `Z0` resistive port; voxelize via `yee_voxel::voxelize_microstrip` on the
F2.2 `lumped_board` Layout; time-step, single-bin DFT the input/output → S21
(reuse the F1.1b.1 driver patterns: CPML or PEC box + time-gating per what worked
there — note the propagation/εeff lesson + the voxelizer z-stack fix).

## Changes (`crates/yee-voxel/**` + `.github/workflows/ci.yml`)

- `crates/yee-voxel/Cargo.toml`: add `yee-filter = { workspace = true }` (for
  `LumpedLadder` + `lumped_board`). yee-filter is WASM-safe pure-math; no cycle
  (it does not dep yee-voxel/yee-fdtd).
- `crates/yee-voxel/src/lib.rs` (or a `lumped_sim.rs` module): `pub fn
  simulate_lumped_board(ladder: &LumpedLadder, substrate: &Substrate, cfg:
  &LumpedSimConfig) -> Vec<(f64, f64)>` returning `(freq_hz, |S21|)` over a sweep
  (or a richer `S21Point`). Builds the board (`lumped_board`), voxelizes, places
  the per-branch lumped elements at the mapped cells, drives/senses, runs FDTD,
  DFTs → S21. `LumpedSimConfig` (grid, n_steps, freq sweep, seed-free).

## DoD (machine-checkable)

1. `cargo fmt --check --all` exit 0; `cargo clippy -p yee-voxel --all-targets -- -D warnings` exit 0 (LOCAL, light).
2. The `#[ignore]`'d gate `fdtd_lumped_001` compiles and is structurally sound.
3. **CI release job** runs `fdtd_lumped_001` and is **GREEN on the branch before
   merge** (CLAUDE.md §4): the FDTD `|S21|` matches the analytic `ladder_s21`
   within a **loose** tolerance (target: passband-vs-stopband behaviour correct —
   |S21|≈0 dB in-band within a few dB, and ≥ ~20 dB rejection at the stopband
   point; an exact-match tol is not expected from a coarse FDTD lumped model).
   The cross-validation (FDTD ≈ circuit) is the published-benchmark.
4. yee-filter / yee-fdtd / yee-layout unchanged (consume their public API);
   WASM-safe crates untouched.

## Validation workflow (use the container — DON'T run FDTD on the host)

Build + iterate the gate LOCALLY in the bounded container
(`YEE_BOX_DIR=worktrees/<lane> scripts/yee-box.sh cargo test -p yee-voxel
--release -- --ignored fdtd_lumped_001 --nocapture`) until green; then push the
branch + a CI release job and confirm GREEN before merge. Never weaken the gate to
a no-op; if the coarse-grid FDTD can't match even a loose tol, surface the
measured vs analytic S21 (the F1.1b.1 precedent: such mismatches are often a
voxelizer/boundary issue, not the synthesis).

## Out of scope

Parasitic SRF/ESR sweeps (F2.1b parts); tolerance-over-FDTD (F2.4 is circuit-level);
KiCad-native footprints (F2.2b); the UI EM panel; multi-port (S11/S22) beyond S21
for the skeleton.

## Why now

It is the goal's named "EM simulation" component, the last lumped engine brick,
and the `LumpedRlcPort` + voxelizer + the F1.1b.1 FDTD-driver patterns all exist.
Heavy but bounded; the container makes it iterable.
