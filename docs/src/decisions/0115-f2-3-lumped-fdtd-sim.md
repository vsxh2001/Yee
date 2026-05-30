# ADR-0115: Filter Phase F2.3 — lumped-LC FDTD EM simulation

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0111 (F2.0 LC ladder + `ladder_s21`), ADR-0114 (F2.2 lumped board),
ADR-0108 (F1.1b.1 FDTD driver patterns + voxelizer z-stack), ADR-0017
(`LumpedRlcPort`), the lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

The lumped-LC goal names **EM simulation**. The circuit `ladder_s21` (F2.0) and
the board (F2.2) exist, and `yee-fdtd` has `LumpedRlcPort::series_rlc` + `yee-voxel`
voxelizes microstrip + the F1.1b.1 FDTD-driver patterns are proven — but nothing
full-wave-simulates the lumped board.

## Decision

Add `yee_voxel::simulate_lumped_board(&LumpedLadder, &Substrate, &LumpedSimConfig)
-> Vec<(f64,f64)>` (freq, |S21|): voxelize the F2.2 `lumped_board`, place each
ladder element as `LumpedRlcPort`(s) on the grid — **series branch** = one
`series_rlc(L,C)` in-line; **shunt branch** = pure-L (`c=∞`) ‖ pure-C (`l=0`) at
the shunt cell (parallel L‖C topology) — drive/sense two ports, time-step, DFT →
S21. yee-voxel gains a `yee-filter` dep (WASM-safe, no cycle).

Gate `fdtd_lumped_001` (`#[ignore]`'d, CI `--release` release job): the FDTD
`|S21|` matches the analytic `ladder_s21` within a **loose** tolerance
(in-band ≈ 0 dB within a few dB, ≥ ~20 dB stopband rejection) — cross-validating
the lumped EM sim against the circuit model. Iterated in the bounded container;
GREEN in CI on the branch **before merge** (CLAUDE.md §4); gate never weakened to
a no-op.

## Consequences

**Ships:** full-wave FDTD of the lumped-LC filter board — the goal's "EM
simulation" component; the lumped track's EM-verified response. With F2.0/F2.1/
F2.2/F2.4, the lumped engine is then complete (synth → parts/BOM → PCB → EM →
tolerance); only the polished-UI lumped track remains (rides on the Dioxus studio).

**Gate:** `fdtd_lumped_001` GREEN in the `fdtd-lumped-gate` CI job before merge.

**Not in scope:** SRF/ESR parasitic sweeps (F2.1b); FDTD-based tolerance;
KiCad-native footprints (F2.2b); the UI EM panel; multi-port beyond S21.

---

## References
- `docs/superpowers/specs/2026-05-30-f2-3-lumped-fdtd-sim-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-3-lumped-fdtd-sim.md`.
