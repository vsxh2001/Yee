# ADR-0108: F1.1b.1 — FDTD coupled-resonator driver (full filter simulation)

**Status:** Accepted (design; implementation queued for a dedicated tick)
**Date:** 2026-05-30
**Related:** ADR-0091 (F1.1a yee-voxel), ADR-0093 (F1.1b.0 extract), ADR-0094
(F1.1b.gate analytic coupled-line k), `FILTER-DESIGN-ROADMAP.md`,
[[project-filter-design-final-goal]]

---

## Context

The product goal has three components: a working web UI, **full filter
simulation**, and KiCad export. Web UI (App.1.2b, ADR-0107) and KiCad export
(ADR-0100/0103/0105/0106) are substantially shipped. The remaining gap is full
simulation: the pipeline is entirely closed-form (synthesis → dimensional
synthesis → layout → manufacturing files) and never solves Maxwell. F1.1b.1 adds
the first full-wave FDTD verification: simulate a coupled-resonator pair and
extract its coupling coefficient `k`, cross-checked against the shipped analytic
coupled-line reference (ADR-0094).

All the primitives exist (voxelizer ADR-0091; FDTD ports / DFT / decay-Q;
`extract_coupling` ADR-0093; analytic k ADR-0094) — F1.1b.1 is the *integration
driver* that wires them into one end-to-end EM solve.

## Decision

Add `run_coupled_pair(&Layout, &CoupledRunConfig) -> CoupledRunResult` to
`yee-voxel` (the native crate that already bridges `Layout → YeeGrid` and deps
`yee-fdtd`; it gains a `yee-filter` dep for `extract_coupling`). It voxelizes a
two-resonator layout, drives it with `LumpedRlcPort`s, runs the FDTD, finds the
two split resonances via the single-bin DFT, and computes
`k = (f_odd² − f_even²)/(f_odd² + f_even²)`.

**Validation (CI-routed):** a `fdtd-coupling-001` gate asserts FDTD-extracted k
matches analytic k within a loose tolerance (≤ 15 % for the skeleton). Because the
FDTD run is multi-minute and the local dev box is memory-constrained (and "use
less CPU" bars local heavy runs), the gate is `#[ignore]`'d for default
`cargo test` and runs in a **dedicated CI `--release` job** on a GitHub runner —
the same "route heavy compute to CI" strategy used for the App.1.2b `trunk build`
(ADR-0107) and consistent with the `mom-001` release-gate pattern.

**CLAUDE.md §4 compliance:** the feature MUST NOT merge until its gate is GREEN
somewhere. The implementation tick pushes the branch, lets the CI FDTD job run,
confirms green on the branch, and only then merges. A never-run `#[ignore]`'d gate
is not acceptable for merge, and the gate must not be weakened to triviality.

**Queued, not done in-session:** this is the heaviest single filter-track
increment and requires a push→CI-FDTD→merge cycle. It is split into its own tick
so the CI cycle is managed properly in a fresh context, rather than risking a
long-session collapse (an explicit standing constraint).

## Consequences

**Ships (next tick):** the first full-wave EM solve in the filter pipeline — FDTD
verification that a dimensioned coupled pair realizes its target coupling, gated
in CI against the analytic reference. This is the "full filter simulation" goal
component's walking skeleton.

**Gate:** `fdtd-coupling-001` (FDTD k vs analytic k, ≤ 15 %), `#[ignore]`'d
locally, GREEN in a CI release job before merge.

**Not in scope:** Qe extraction (F1.1b.2), surrogate-BO EM-in-loop (F1.2.1), the
full Swanson-hairpin end-to-end filter gate, GPU FDTD.

---

## References
- ADR-0091 / ADR-0093 / ADR-0094 (the primitives this integrates).
- `docs/superpowers/specs/2026-05-30-f1-1b-1-fdtd-coupled-resonator-driver-design.md`;
  `docs/superpowers/plans/2026-05-30-f1-1b-1-fdtd-coupled-resonator-driver.md`.
