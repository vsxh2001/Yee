# F1.1b.1 — FDTD coupled-resonator driver — Plan

**Spec:** `2026-05-30-f1-1b-1-fdtd-coupled-resonator-driver-design.md` · **ADR:** ADR-0108

> **Status: queued for a dedicated implementation tick.** This plan locks the
> design + the CI-routed validation strategy. The heavy implementation (a
> multi-minute FDTD gate that must pass in CI before merge) is intentionally NOT
> run in the long current session — dispatch it from a fresh loop-tick.

## Lane (implementation tick)
`crates/yee-voxel/**` (driver + the `#[ignore]`'d gate) AND
`.github/workflows/ci.yml` (the release FDTD job). Do NOT edit `yee-layout` /
`yee-filter` / `yee-fdtd` internals — consume their public API. Adding
`yee-filter` as a `yee-voxel` dep is allowed (it is native-only already). Out of
lane → finding.

## Base
New worktree off `main` at dispatch time (re-fetch first — cloud-race lesson
[[feedback-cloud-routine-races-local-loop]]). Branch `feature/f1-1b-1-fdtd-coupled-driver`.

## Pattern / context files (READ first)
- `crates/yee-voxel/src/lib.rs` — `voxelize_microstrip`; mirror its style; add the
  driver beside it.
- `yee-fdtd`: the `LumpedRlcPort` series-RLC drive + `inductor_current` /
  `capacitor_voltage` readout; the single-bin DFT at `cavity_resonance.rs:273`;
  the decay-fit Q at `cavity_q.rs:140`; per-cell ε_r `with_eps_r_cells` + PEC
  masks. (grep these exact symbols.)
- `crates/yee-filter/.../extract*` — `extract_coupling` (the `(f2²−f1²)/(f2²+f1²)`
  formula) + `extract_q_ringdown`.
- `crates/yee-layout/.../coupled_microstrip` + `coupling_coefficient` — the
  analytic k reference the gate compares against.
- the `mom-001` release gate wiring in `.github/workflows/ci.yml` + its
  `#[ignore]` test — mirror BOTH the `#[ignore]` idiom and the release-job idiom.

## Steps (implementation tick)
1. `yee-voxel`: `run_coupled_pair(&Layout, &CoupledRunConfig) -> CoupledRunResult`
   (voxelize → LumpedRlcPorts → fdtd run → DFT/FFT → split freqs → k). Document
   all public items.
2. `crates/yee-voxel/tests/fdtd_coupling_001.rs` — `#[ignore]`'d; builds a known
   coupled pair, runs the driver, asserts FDTD k vs analytic k within ≤ 15 %.
3. `.github/workflows/ci.yml` — a job (release) that runs
   `cargo test -p yee-voxel --release -- --ignored fdtd_coupling_001` on a GitHub
   runner. Budget the wall-time (coarse grid; cap timesteps).

## Verify (implementation tick)
- LOCAL (light): `cargo fmt --check --all`; `cargo clippy -p yee-voxel
  --all-targets --jobs 2 -- -D warnings`. Do NOT run the FDTD gate locally (multi-
  minute + OOM-prone box).
- CI (the real gate): push the branch; confirm the new release FDTD job is GREEN
  on the branch; ONLY THEN merge (CLAUDE.md §4 — never merge a never-passed solver
  gate; never weaken the gate to triviality).

## Escape hatch
Blocked > 30 min on the physics (FDTD k won't land within a sane loose tolerance;
the split resonances are not cleanly separable on a tractable grid; CI wall-time
explodes) → STOP and surface the measured k, the grid/step budget, and the
separation. Do NOT merge an unvalidated driver; do NOT weaken the gate to a
no-op; do NOT run the FDTD locally to brute-force it.

## Done when
DoD 1–5 of the spec pass; the CI FDTD gate is GREEN on the branch before merge;
`git diff --stat base..HEAD` = `crates/yee-voxel/**` + `.github/workflows/ci.yml`
(+ the committed docs); WASM-safe crates untouched.
