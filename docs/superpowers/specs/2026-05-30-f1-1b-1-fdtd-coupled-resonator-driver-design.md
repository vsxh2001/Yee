# F1.1b.1 — FDTD coupled-resonator driver (full filter simulation) — Design Spec

**ADR:** ADR-0108 · **Date:** 2026-05-30 · **Status:** Accepted (design; implementation queued)

## Goal
Deliver the **full-wave simulation** component of the product goal: an end-to-end
FDTD driver that takes a pair of coupled microstrip half-wave resonators,
voxelizes them, excites and runs the FDTD, and extracts the inter-resonator
coupling coefficient `k` — validated against the shipped analytic coupled-line
reference (F1.1b.gate, ADR-0094). This is the first real EM solve in the
filter-design pipeline (everything before it is closed-form synthesis/geometry).

## Why this is the (2) "full filter simulation" brick
The pipeline today is closed-form: synth → coupling matrix → dimensional
synthesis → layout → manufacturing files. None of it solves Maxwell. F1.1b.1 adds
the FDTD verification loop: it confirms that the *dimensioned geometry* actually
realizes the *target coupling*, by full-wave simulation. It composes primitives
that ALL already exist (per the Explore recon recorded in
[[project-filter-design-final-goal]] and ROADMAP §9):
- `yee_voxel::voxelize_microstrip(&Layout) → YeeGrid` (F1.1a, ADR-0091).
- `yee-fdtd` per-cell ε_r (`with_eps_r_cells`), PEC masks, `LumpedRlcPort`
  series-RLC drive + `inductor_current`/`capacitor_voltage` readout.
- the single-bin DFT pattern at `cavity_resonance.rs:273`; the decay-fit Q at
  `cavity_q.rs:140`.
- `yee_filter::extract_coupling` (k from the two split peaks
  `(f2²−f1²)/(f2²+f1²)`) + `extract_q_ringdown` (F1.1b.0, ADR-0093).
- the analytic reference: `yee_layout::coupled_microstrip` k (F1.1b.gate, ADR-0094).

## Crate / placement
A new driver function in a crate that may depend on BOTH `yee-fdtd` and
`yee-voxel` and `yee-filter` — i.e. **`yee-voxel`** (already deps yee-layout +
yee-fdtd; add yee-filter) OR a new `yee-fdtd-driver` crate. **Decision:** put it
in `yee-voxel` (it already bridges layout→grid; adding the run+extract keeps the
EM-driver story in one native crate). It must NOT touch `yee-layout`/`yee-filter`
WASM-safety (those stay pure; the driver is native-only, like yee-voxel).

## Proposed API (refine in implementation)
```rust
// yee-voxel
pub struct CoupledRunConfig { /* grid resolution, f0, span, n_steps, ports */ }
pub struct CoupledRunResult { pub k: f64, pub f_even: f64, pub f_odd: f64 /* , Qe? */ }
pub fn run_coupled_pair(layout: &Layout, cfg: &CoupledRunConfig) -> CoupledRunResult;
```
Driver steps: voxelize the two-resonator `Layout` → install two `LumpedRlcPort`s
(weak coupling probes) → time-step yee-fdtd → single-bin DFT at a sweep of
frequencies (or FFT the port readout) → locate the two split resonances
`f_even`/`f_odd` → `k = (f_odd² − f_even²)/(f_odd² + f_even²)`.

## Validation gate (the hard part — CI-routed)
- **`fdtd-coupling-001`**: build a coupled-microstrip pair whose ANALYTIC k (via
  `coupled_microstrip` / `coupling_coefficient`) is known; run the FDTD driver;
  assert the FDTD-extracted k matches the analytic k within a **loose** tolerance
  (target ≤ 15 % for the walking skeleton — FDTD on a coarse grid is approximate;
  tighten later). This is the published-benchmark cross-check required by
  CLAUDE.md §4.
- **The FDTD run is multi-minute** ⇒ the gate is `#[ignore]`'d for the default
  `cargo test` (like `mom-001`/`fem_eig_003`), and runs in a **dedicated CI job in
  `--release`** (mirror the structure of the `mom-001` release gate / the
  `gpu-nightly` gating pattern). It must run on a GitHub runner, NEVER on the
  local memory-constrained box.

## CLAUDE.md §4 compliance (CRITICAL — do not skip)
A solver feature MUST NOT merge with a never-passed gate. Therefore the
implementation tick MUST: push the feature branch → let the new CI release job
run the FDTD gate → confirm it is GREEN on the branch → only THEN merge. The gate
being `#[ignore]`'d locally is fine; the gate being un-run anywhere is NOT. If the
CI FDTD job cannot be made to pass within budget, the increment STOPS and surfaces
(do not merge an unvalidated solver, do not weaken the gate to triviality).

## DoD (implementation tick)
1. `cargo clippy -p yee-voxel --all-targets -- -D warnings` exit 0 (local, light).
2. `cargo fmt --check --all` exit 0.
3. The `#[ignore]`'d `fdtd-coupling-001` gate compiles and is structurally sound
   (asserts FDTD k vs analytic k within tolerance).
4. A CI job runs the gate in `--release` and is **GREEN on the feature branch
   before merge** (the real validation).
5. `yee-layout` / `yee-filter` WASM-safety unaffected (no native dep added to
   their default path).

## Out of scope
Qe / external-coupling extraction in the same gate (a follow-on F1.1b.2);
surrogate-BO EM-in-loop (F1.2.1); the full Swanson-hairpin end-to-end filter
gate (later); GPU FDTD. Walking-skeleton: ONE coupled pair, ONE loose k gate.

## Why this is queued, not done-in-session
The increment requires a push→CI-FDTD→merge cycle (the gate cannot run locally
under "use less cpu" + the OOM-constrained box), and the FDTD driver is the
heaviest single increment in the filter track. It is split out so the
implementation runs in a fresh session/loop-tick with the CI cycle managed
properly, rather than risking a long-session collapse.
