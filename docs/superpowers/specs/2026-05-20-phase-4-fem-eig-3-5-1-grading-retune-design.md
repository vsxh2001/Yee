# Phase 4.fem.eig.3.5.1 — CFS-PML grading-parameter retune

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3.5.1 (Complex Frequency Shifted PML grading-parameter
retune; ablation across the H1/H2/H3 hypothesis tree to retire the
fem-eig-003 strict absorption-floor band and the fem-eig-006 magnitude
gate that the OOOOOOOOO P5 measurement leaves ~30 dB above target).
**Depends on:** Phase 4.fem.eig.3.5 (CFS-PML implementation P1-P7
shipped under OOOOOOOOO; `AbcOrder::CfsPml(PmlConfig)` variant,
`extend_mesh_with_pml`, anisotropic per-tet `ε_tensor(ω)` assembly,
default `(σ_max, α_max, κ_max=5, m=3, thickness=6)` grading).
**Blocks:** retirement of the `#[ignore]` on
`fem_eig_003_strict_absorption_floor_gate`,
`fem_eig_003_strict_passive_bound_continuum_limit`, and
`fem_eig_006_magnitude_bounded` (one ignored gate per fixture became
two on fem-eig-003 with the OOOOOOOOO P5 measurement; total three).

## 1. Goal

Identify the CFS-PML grading-parameter configuration that retires the
fem-eig-003 strict `[-60, -40] dB` absorption band and the fem-eig-006
`|S_{11}(30 GHz)| < 0.1` magnitude gate, then lock in the winning
parameters as the new `PmlConfig::default()` and the new
`PmlConfig::resolved(freq_hz, mesh_meta)` formulae. Both strict gates
flip from `#[ignore]`'d to CI-default without weakening tolerances.

The OOOOOOOOO P5 escape-hatch measurement leaves fem-eig-003 at
`|S_{11}| ∈ [0.281, 0.423]` ≡ `[-11.0, -7.48] dB` across 8-12 GHz on
the existing `(24, 12, 36)` cavity + 6-cell PML shell (~72 k extended
tets, runtime ~140 s in `--release`). This is ~10 dB better in dB than
the Phase 4.fem.eig.3 2nd-order Engquist-Majda baseline
(`[-2.22e-2, -2.86e-5] dB`; ~35 dB above the v3-spec `[-45, -35] dB`
window) but **~30 dB above** the v3.5 spec §6 `[-60, -40] dB` target
window. fem-eig-006 `|S_{11}(30 GHz)| = 0.926` — the CFS `α_α > 0`
causality canary passes (the value is finite, not `NaN`/`Inf`) but the
magnitude gate misses by ~9.3×. The ~30 dB miss is structural, not a
mesh-density artefact: NNNNNNNNN already showed mesh refinement
saturates at `~ 2× dB/level`, which would require ~14 further mesh
doublings to close on its own — clearly the wrong knob.

This phase performs a focused ablation across three hypotheses against
two production fixtures and ships the winning parameter set, plus the
per-axis `h_α` heuristic upgrade.

## 2. Background

OOOOOOOOO shipped Phase 4.fem.eig.3.5 P1-P7 against ADR-0043 with the
literature-recommended defaults: `kappa_max = 5` (Roden-Gedney 2000
Table I; calibrated for FDTD waveguide-discontinuity benchmarks),
polynomial grading order `m = 3`, `thickness_cells = 6`, and the
analytic
`σ_max = (m + 1) / (150 π h_cell √ε_r)`,
`α_max = 2 π f_centre ε_0` calibrations from Roden-Gedney 2000 §III/IV.
These defaults retired the OOOOOOOOO P4 PML-end-to-end smoke
(`pml_assembly_matches_scalar_on_zero_thickness`, `pml_assembly_finite_at_dc`)
and produced finite, sub-unity `|S_{11}|` on both production fixtures —
so the implementation is correct. **The defaults are sub-optimal for
the frequency-domain FEM regime**, three candidate root causes:

**H1 — single-`h_cell` heuristic mis-predicts mesh spacing on
aspect-ratio cells.** `PmlConfig::resolved(freq_hz, h_cell)` computes
the band-centre `σ_max` from a *single* characteristic cell length
`h_cell ≈ mean_tet_edge_length`. For a high-aspect-ratio cavity like
the WR-90 stub (broad wall 22.86 mm = 24 cells of 0.952 mm; narrow wall
10.16 mm = 12 cells of 0.847 mm; axial 30 mm = 36 cells of 0.833 mm)
the three axes' `h_α` differ by ~14 %. Roden-Gedney 2000 §III §IV
prescribe a *per-axis* `σ_α_max = (m + 1) / (150 π h_α √ε_r)` because
the optimal PML conductivity depends on the wavelength sampled by the
*α-axis* discretisation in the α-PML. The single-`h_cell` resolver
predicts a `σ_max` that is correct on average but wrong on every axis
individually, with the error compounding through `Λ(ω)` to push the
absorption floor up by a measurable margin. fem-eig-006 (100 : 10 : 1
extreme aspect) amplifies the same effect by ~10×.

**H2 — `κ_max = 5` is FDTD-calibrated.** Roden-Gedney 2000 Table I's
`κ_max = 5` benchmark is the FDTD time-domain regime where the
real-coordinate stretch trades against the courant stability bound. In
frequency-domain FEM there is no stability bound and `κ` only affects
the real-axis decay of propagating modes versus reflections off the
PEC truncation surface. Berenger 2002 §V parameter sweeps suggest
`κ_max ∈ [1.5, 3]` for FD regime — `κ_max = 5` *over-stretches* the
PML, pushing physical-wavelength modes toward the PEC outer surface
where small numerical artefacts amplify into reflection.

**H3 — polynomial grading order `m = 3` over-ramps on a 6-cell shell.**
The Roden-Gedney 2000 polynomial grading `σ(d) = σ_max (d/D)^m`
concentrates the absorption near the outer truncation surface. With
`m = 3` and a thin 6-cell shell, the first 3 cells contribute
`(0.5/6)^3 + (1.5/6)^3 + (2.5/6)^3 ≈ 0.094` of the total absorption
budget — 90.6 % of the absorption happens in the outer 3 cells.
Berenger 2002 §IV parameter sweep shows `m = 2` with a slightly thicker
8-cell shell distributes absorption more uniformly and produces
~5-10 dB better reflection floors on equivalent-cell-budget meshes.

The three hypotheses are **independent knobs** with potentially
overlapping effect: any of H1, H2, H3 alone may be sufficient; some
combination almost certainly is. This spec sets up an ablation grid
that decomposes the search.

Reference: Berenger, J.-P., "Numerical reflection from FDTD-PMLs: a
comparison of the split PML with the unsplit and CFS PMLs," *IEEE
Transactions on Antennas and Propagation* 50(3) (March 2002),
pp. 258-265. DOI 10.1109/8.999615. The canonical CFS-PML parameter
sweep study; figures 4-7 sweep `(σ_max, κ_max, m, thickness)` on a 2D
canonical waveguide benchmark and identify the
"`κ_max ∈ [1.5, 3], m ∈ {2, 3}`" basin for off-normal incidence — the
exact regime fem-eig-003 + fem-eig-006 occupy.

## 3. Hypothesis tree

The decision logic is a three-leaf depth-first walk: cheapest
hypothesis first, ship on first retire, fall through to deeper levels
only if the prior leaf misses.

### H1 — per-axis `h_α` back-inference

Replace `PmlConfig::resolved(freq_hz, h_cell)` with
`PmlConfig::resolved(freq_hz, &PmlMeshMeta)` where

```rust
pub struct PmlMeshMeta {
    /// Axis-aligned bounding box extents (m), one per axis.
    pub extents: [f64; 3],
    /// Per-axis cell count of the original cavity mesh.
    pub cell_counts: [usize; 3],
}
```

and compute per-axis `h_α = extents[α] / cell_counts[α]`. The new
resolver derives per-axis `σ_α_max`, `α_α_max` for use by the per-axis
`s_α(ω)` evaluator (which already exists in
`open_boundary.rs::pml_stretching_lambda` — its `s_for(d_α)` closure
becomes parameterised by `α`).

**Decision criterion:** if H1 alone (with the OOOOOOOOO defaults
`κ_max = 5, m = 3, thickness = 6`) moves the fem-eig-003 worst-case
`|S_{11}|` below `-25 dB` *and* moves fem-eig-006 below `0.5`, **ship
H1 standalone**. The `-25 dB` bar is the half-way point between the
OOOOOOOOO measurement (`-7.5 dB`) and the spec gate (`-40 dB`); if H1
clears half the gap on its own it dominates the simpler-knob
hypotheses.

### H2 — `κ_max` sweep with per-axis `h_α` (H1+H2)

If H1 alone misses, sweep `κ_max ∈ {1, 1.5, 2, 3, 5, 7}` with H1
enabled (per-axis `h_α`) and `m = 3, thickness = 6` fixed at the
OOOOOOOOO defaults. The smallest `κ_max` that retires the fem-eig-003
`[-60, -40] dB` band on a swept basis is the candidate default.

**Decision criterion:** if a `κ_max` value in `{1.5, 2, 3}` clears the
fem-eig-003 swept band *and* dominates `κ_max ∈ {5, 7}` on the same
fixture (i.e. monotone-or-better dB at smaller `κ_max`), **ship the
new `κ_max` default**. Ship the smallest value that achieves the gate
(parsimony — smaller `κ_max` minimises PEC-truncation-surface
sensitivity).

### H3 — polynomial order + thickness sweep (H1+H2+H3)

If H1+H2 still misses, sweep `m ∈ {2, 3, 4}` with `thickness_cells ∈
{6, 8, 10}` covariate, with the H1+H2 winners fixed. The pair
`(m, thickness)` that retires both fem-eig-003 and fem-eig-006 strict
gates *with the smallest extended-mesh tet count* wins. The
extended-tet penalty per added shell layer on fem-eig-003 is `~ 1 728`
tets / layer (24 × 12 face × 6 Kuhn); on fem-eig-006 (the 100 : 10 : 1
fixture, with PML on one short face only) is `~ 600` tets / layer.

**Decision criterion:** the `(m, thickness)` pair that **simultaneously
retires both fem-eig-003 and fem-eig-006 strict bands** ships. If
multiple pairs retire, pick the smallest `(m × thickness)` product
(parsimony).

If no `(m, thickness)` pair retires both fixtures, see §7 risks: ship
per-fixture overrides via the `pml_config` builder kwarg and
**document the constraint** rather than weaken either gate.

## 4. Ablation grid

Maximum size:

- H1 standalone: **1 configuration** (per-axis `h_α`, defaults
  otherwise).
- H2 sweep: **6 configurations** (`κ_max ∈ {1, 1.5, 2, 3, 5, 7}`).
- H3 sweep: **3 × 3 = 9 configurations** (`m × thickness`).

Two fixtures × (1 + 6 + 9) = **2 × 16 = 32 runs**. Loose upper bound
"2 × 4 × 6 = 48" in the brief allows for retries / one-off
investigations.

**Stopping rule:** the sweep tool runs each configuration on
fem-eig-003 first (cheaper at ~70 k tets) and only re-runs fem-eig-006
if fem-eig-003 retires. The first `(κ_max, m, thickness)` that
retires **both** fixtures with `|S_{11}| < -40 dB` worst-case on
fem-eig-003 and `|S_{11}| < 0.1` on fem-eig-006 ends the sweep.
fem-eig-003 wall-time is ~140 s `--release`; fem-eig-006 is faster
(~30 s; single-frequency, smaller mesh). Worst-case sweep wall-time
~ (32 × 140 s) ≈ 75 min in `--release`.

## 5. Public API

No new types. The only signature changes are:

1. `PmlConfig::resolved(freq_hz: f64, h_cell: f64) -> Self` becomes
   `PmlConfig::resolved(freq_hz: f64, mesh_meta: &PmlMeshMeta) -> Self`.
   The single-`h_cell` overload is removed (it was an internal
   helper).
2. `OpenBoundarySolver::with_cfs_pml(self, config: PmlConfig) -> Self`
   internally constructs a `PmlMeshMeta` from the input mesh at
   builder time. Callers' public signature unchanged.
3. `PmlConfig::default()` returns the winning `(κ_max, m,
   thickness_cells)` from the §4 ablation (TBD per §3 decision tree).
4. `pml_stretching_lambda` (internal) takes per-axis
   `(σ_α_max, α_α_max, h_α)` triples instead of a single
   `(σ_max, α_max, h_cell)`. The per-axis grading polynomial is
   evaluated with the α-axis cell size `h_α`.

`PmlConfig` field semantics unchanged: `sigma_max`, `alpha_max`,
`kappa_max`, `m`, `thickness_cells`. The mesh-derived per-axis
parameters live in a new private `ResolvedPmlConfig` carrier internal
to `open_boundary.rs`.

Python binding: `yee.fem.solve_open_cavity`'s `pml_config` kwarg
semantics unchanged. The Python user does not need to pass mesh
extents; the Rust builder derives them from the input mesh.

## 6. Validation

Un-ignore exactly three strict gates **after** the §4 ablation
identifies a winning parameter set that retires both fixtures:

1. `fem_eig_003_strict_absorption_floor_gate` in
   `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
   (currently `#[ignore]`'d with OOOOOOOOO P5 measurement docstring
   `|S_{11}| ∈ [0.281, 0.423]`).
2. `fem_eig_003_strict_passive_bound_continuum_limit` in the same
   file (passive bound `|S_{11}| < 1` strictly; currently
   `#[ignore]`'d but already comfortably satisfied at the OOOOOOOOO
   defaults — the un-ignore is gated on the same retune retiring the
   absorption-floor gate).
3. `fem_eig_006_magnitude_bounded` in
   `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
   (`|S_{11}(30 GHz)| < 0.1`; OOOOOOOOO P5 measured `0.926`).

The un-ignore happens *only* after a sweep CSV verifies both fixtures
pass under the new defaults. If the §3 decision tree exhausts without
retiring both fixtures (§7 risk a), the un-ignore is **not**
performed; instead the per-fixture override path lands as a documented
escape hatch and the strict gates stay `#[ignore]`'d, with a queue to
Phase 4.fem.eig.3.5.2 for `α_α(d)` grading.

## 7. Risks

- **(a) 2D-ill-posed defaults.** Worst case, the fem-eig-003 (24 ×
  12 × 36 cuboid, moderately-anisotropic) and fem-eig-006 (100 × 10 × 1
  extreme-anisotropic) fixtures require *different* `(κ_max, m,
  thickness)` to retire, and no single default set retires both
  simultaneously. Mitigation: the
  `OpenBoundarySolver::with_cfs_pml(cfg)` builder already accepts a
  user-supplied `PmlConfig`; ship the ablation result as the **default**
  best-of-both-fixtures compromise and document per-fixture override
  values in `docs/src/tutorials/07-fem-open-cavity.md`. Both strict
  gates use the override values explicitly, so the un-ignore proceeds.

- **(b) `α_α(d)` grading is the other candidate cause.** Berenger 2002
  §VI shows that for canonical evanescent-mode problems the CFS
  frequency-shift parameter `α_α` is also profitably **graded** as
  `α(d) = α_max · (1 - d/D)^m_α` rather than constant. The
  OOOOOOOOO defaults use a constant `α_max`. This ablation is *not*
  in v3.5.1 scope; deferred to v3.5.2 if H1+H2+H3 retune still
  misses by > 5 dB.

- **(c) Sweep wall-time dominates CI budget.** 32 ablation runs at
  ~140 s `--release` = 75 min total; this is **not** a CI gate but a
  one-off design exploration. The sweep tool lives in
  `tools/cfs_pml_grading_sweep.rs` as a yee-validation example
  binary, runs locally (or in a manual CI workflow), and emits a CSV
  that the R3 step analyses. Only the *winning* defaults land in CI
  as the (single) un-ignored gate.

- **(d) Per-axis `h_α` resolver alters the OOOOOOOOO measurement
  baseline.** The OOOOOOOOO measurement was taken with the
  single-`h_cell` resolver; switching to per-axis `h_α` shifts the
  baseline before any other knob moves. The R2 sweep tool re-measures
  the H1-standalone baseline as configuration `(H1=on, κ_max=5, m=3,
  thickness=6)` and records it in the CSV; the dB-improvement
  attribution between H1 / H2 / H3 leaves this baseline as the
  reference.

- **(e) Complex-symmetric stiffness preservation.** Per-axis `Λ(ω)`
  remains diagonal in the global frame (the per-axis
  `(σ_α, κ_α, α_α)` only sharpen the axis-specific entries of `Λ`;
  cross-axis entries stay zero). The complex-LDLᵀ path in
  `faer::sparse::FaerLuSolver<Complex64>` is preserved bit-for-bit;
  ADR-0043 decision (4) is unaffected.

## 8. Lane

Spec file:
`docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune-design.md`

Implementation lane (declared here for the R1-R5 plan):

- `crates/yee-fem/src/open_boundary.rs` — `PmlConfig::resolved`
  signature change, per-axis `s_α(ω)` evaluator, `PmlMeshMeta`
  derivation in `with_cfs_pml`.
- `crates/yee-fem/src/lib.rs` — re-export `PmlMeshMeta`.
- `crates/yee-fem/tests/pml_open_boundary_assembly.rs` — extend
  `pml_assembly_matches_scalar_on_zero_thickness` to also cover the
  new per-axis path with `thickness = 0` (no-op equivalence).
- `tools/cfs_pml_grading_sweep.rs` *(create)* — yee-validation example
  binary that runs the §4 ablation grid and emits a CSV.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` —
  un-ignore both strict gates (after R3 analysis).
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs` —
  un-ignore the magnitude gate (after R3 analysis).
- `docs/src/tutorials/07-fem-open-cavity.md` — note on the new
  defaults and per-fixture override pattern.
- `ROADMAP.md` — Phase 4.fem.eig.3.5.1 entry from planned to
  shipped.

Out of lane: `yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`, `yee-cuda`,
`yee-plotters`, `yee-fdtd`, `yee-surrogate`.

## 9. References

- Berenger, J.-P., "Numerical reflection from FDTD-PMLs: a comparison
  of the split PML with the unsplit and CFS PMLs," *IEEE Transactions
  on Antennas and Propagation* 50(3) (March 2002), pp. 258-265.
  DOI 10.1109/8.999615. The canonical CFS-PML parameter sweep study.
- Roden, J. A. and Gedney, S. D., "Convolutional PML (CPML): An
  efficient FDTD implementation of the CFS-PML for arbitrary media",
  *IEEE Microwave and Wireless Components Letters* 10(5) (May 2000),
  pp. 27-29. The Table-I / §III/IV defaults the OOOOOOOOO baseline
  inherits.
- Kuzuoglu, M. and Mittra, R., "Frequency dependence of the
  constitutive parameters of causal perfectly matched anisotropic
  absorbers", *IEEE MWCL* 6(12) (1996), pp. 447-449. The CFS
  modification (`α_α > 0`) v3.5 implements.
- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
  — Phase 4.fem.eig.3.5 parent spec (CFS-PML implementation).
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-cfs-pml.md`
  — parent plan (P1-P7 ladder OOOOOOOOO shipped).
- `docs/src/decisions/0043-phase-4-fem-eig-3-5-cfs-pml-scope.md` —
  parent ADR (§risks queues this spec under "PML grading parameter
  sensitivity").
- `docs/src/decisions/0044-phase-4-fem-eig-3-5-1-grading-retune.md`
  — this spec's scope ADR.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  §"OOOOOOOOO P5 status" — the `|S_{11}| ∈ [0.281, 0.423]`
  measurement that motivates this retune.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
  §"OOOOOOOOO P5 status" — the `|S_{11}(30 GHz)| = 0.926`
  measurement.
- CLAUDE.md §3, §4, §10.
