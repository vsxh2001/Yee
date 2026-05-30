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

---

## Update (2026-05-30) — what actually SHIPPED (resonant → propagation)

The originally-planned **resonant coupled-resonator split** method above proved
unworkable, and the investigation uncovered a deeper bug. What merged to `main`
(`afd1eff`) is a **propagation-based** ε_eff gate. Full blow-by-blow on PR #1.

**Root cause found (the real bug behind the whole saga):** `yee-voxel::
voxelize_microstrip` left a **one-cell air gap between the ground plane and the
dielectric** (trace PEC at `k = 1 + n_sub`, dielectric only at `k = 1..=n_sub`).
That series air capacitance dragged the FDTD microstrip ε_eff ~20 % low (≈2.5 vs
the analytic 3.33), silently corrupting every resonance/coupling measurement.
**Fixed:** dielectric fills `k = 0..n_sub`, trace PEC at `k_top = n_sub`, so
ground→trace is `n_sub·dx = h` of pure dielectric. `voxel_001` pins the corrected
z-stack.

**Why resonant-split was abandoned:** no box is simultaneously high-Q and
non-confining. A small closed PEC box confines the fringing/air-gap fields that
set the even/odd ε_eff difference (split too small); a large PEC box becomes a
resonant cavity (box modes swamp the spectrum); an open CPML box collapses the
resonator Q (no detectable peaks). Also `CpmlParams::for_grid` (all-six-face) is
late-time unstable for a microstrip whose PEC ground / high-ε substrate run into
the boundary.

**What shipped (the gate):** a driven, **time-gated propagation** measurement of
the microstrip phase velocity → ε_eff. `run_line_eeff` drives one end of a long
straight line, records `Ez` at two probe planes Δx apart (Δx < λ_g, so no 2π
phase ambiguity), extracts the phase advance via a single-bin DFT gated to the
forward pulse before the far-wall reflection returns: `v_p = ω·Δx/Δφ`,
`ε_eff = (c/v_p)²`. `run_coupled_line_eeff` does the same in-phase (even) /
anti-phase (odd) for the coupled even/odd ε_eff. Gates (both `#[ignore]`'d, ≤15 %,
CI `--release`): **`fdtd-line-eeff-001`** (single-line ε_eff vs Hammerstad-Jensen:
measured 3.329 vs 3.325, **0.13 %**) and **`fdtd-line-eeff-coupled-001`** (even
4.74 %, odd 1.25 %). This is the first validated full-wave EM solve in the filter
pipeline — the "full filter simulation" product-goal component.

**Enabler:** the bounded Docker dev container (`scripts/yee-box.sh`, main
`6bcd026`) let the multi-minute FDTD gate be iterated **locally** (host-safe),
which is how the voxelizer bug was found. The earlier "never run FDTD locally"
constraint is superseded by the box.

**Tooling note:** the spec/plan filenames still carry the original
`-coupled-resonator-driver` name; their *content* is the resonant approach — read
this Update + PR #1 for the shipped propagation method.
