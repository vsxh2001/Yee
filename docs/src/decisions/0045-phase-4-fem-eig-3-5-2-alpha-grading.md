# ADR-0045: Phase 4.fem.eig.3.5.2 scope — `alpha_alpha(d)` polynomial grading + extended H3 thickness ablation, retires fem-eig-003 + fem-eig-006 strict gates

## Status

Accepted — 2026-05-20 (spec + plan; implementation queued under the
Phase 4.fem.eig.3.5.2 S1-S5 ladder).

## Context

Phase 4.fem.eig.3.5.1 (ADR-0044) shipped under Track QQQQQQQQQ
against base SHA `56386f1`. The R1-R5 ladder landed the per-axis
`h_α` resolver, the `PmlMeshMeta` carrier, the
`crates/yee-validation/examples/cfs_pml_grading_sweep.rs` example
binary, and the partial H1/H2/H3 ablation grid. The R3 comment
block at `crates/yee-fem/src/open_boundary.rs:319-346` records the
escape-hatch outcome:

- R1 per-axis `h_α` alone moved fem-eig-003 worst-case from
  `-7.48 dB` (OOOOOOOOO baseline) to `-21.74 dB` — ~14 dB
  improvement.
- H2 `kappa_max` sweep across `{1, 1.5, 2, 3, 5}` clustered
  within ~1 dB at `~ -22 dB` worst-case. `kappa_max` is not the
  binding constraint at this mesh resolution (consistent with
  Berenger 2002 §V).
- H3 most-aggressive probe `(kappa=2, m=4, thickness=10)`
  reached `[-58.13, -35.45] dB`. Band minimum already inside the
  `[-60, -40] dB` target window; worst-case still ~5 dB short of
  the `-40 dB` retire threshold.

QQQQQQQQQ R3 elected the escape-hatch path per ADR-0044
decision 5: ship the per-axis resolver as the only behavioural
change, leave `PmlConfig::default()` at the OOOOOOOOO values
`(kappa_max=5, m=3, thickness=6)`, leave the three strict
`#[ignore]`'s in place on the production gates. QQQQQQQQQ R4
updated the gate docstrings with the new H3-probe baseline and
queued v3.5.2 per ADR-0044 §Consequences "`α_α(d)` polynomial
grading is the only remaining CFS-PML knob not exhausted by
v3.5.1".

Two candidate knobs are now in scope, both queued by ADR-0044
§Consequences + spec §7 (b):

1. **`alpha_alpha(d)` polynomial grading.** Berenger 2002 §VI
   shows the constant-`alpha_max` formulation produces a
   discontinuity at the cavity-PML inner boundary that lifts the
   worst-case reflection floor by ~5-10 dB. Grading
   `alpha_alpha(d) = alpha_max · (1 - d/D)^n` from `alpha_max`
   at the inner boundary to `0` at the outer truncation smooths
   the discontinuity; the §VI canonical sweep reports
   ~10-20 dB at the inner boundary and ~5-10 dB worst-case.
2. **Extended H3 thickness.** The v3.5.1 grid capped at
   `thickness_cells = 10`. Linear extrapolation from the v3.5.1
   measurements (~7 dB / 2-cell layer) suggests
   `thickness = 16` could land near `-46 dB` worst-case —
   inside the strict band but only by ~6 dB.

The 5 dB shortfall after v3.5.1 exactly matches the magnitude
that Berenger 2002 §VI attributes to the constant-`alpha_max`
artefact, making `alpha_alpha(d)` grading the most likely
binding constraint. Extended thickness is the parallel candidate
that targets absorption-budget magnitude rather than
inner-boundary smoothness; the two knobs are independent, and
QQQQQQQQQ probe data is insufficient to decide a priori.

## Decision

Phase 4.fem.eig.3.5.2 runs an 18-configuration H4 ablation grid
across the two new knobs simultaneously against both production
fixtures, picks the first quadruple `(kappa_max, m, thickness_cells,
alpha_grading_order)` that retires both strict gates, and ships
it as the new `PmlConfig::default()`.

Five load-bearing decisions:

1. **Combined ablation, single phase.** `alpha_alpha(d)` grading
   and extended thickness are independent knobs but the joint
   18-row grid (`kappa=2 × m ∈ {3, 4} × thickness ∈ {12, 14, 16}
   × alpha_grading_order ∈ {0, 1, 2}`) costs only ~42 min
   `--release` (worst case) and gives a clean decision in one
   sweep instead of two sequential phases. The QQQQQQQQQ H3
   probe finding identified `kappa=2` as the local optimum, so
   the v3.5.2 grid holds `kappa_max` fixed at 2.
2. **First retire wins; ship the simplest quadruple.** Per spec
   §6 stopping rule, the sweep emits `WINNER` on the first row
   retiring both strict gates and exits. If multiple rows would
   retire on a continued sweep, the first-retire choice
   parsimoniously minimises `m × thickness × alpha_grading_order`
   among retiring configurations.
3. **`alpha_grading_order: usize` is a new `PmlConfig` field,
   default `0`.** The `0` default recovers v3.5.1 constant-`alpha_max`
   behaviour bit-for-bit, so callers explicitly passing a full
   `PmlConfig` see no behavioural change. The
   `PmlConfig::default()` quadruple is the only place the new
   field's non-zero value lands in CI.
4. **API addition only; no public-type changes.** No new types,
   no signature changes on `with_cfs_pml`, no Python kwarg
   restructure. The Python `yee.fem.solve_open_cavity`'s
   `pml_config` kwarg gains the optional `alpha_grading_order`
   field. The Rust `PmlConfig` struct grows by one field; all
   existing call sites continue to work via the
   `..Default::default()` struct-update pattern.
5. **Strict gates un-ignore only after retire.** If the §6
   decision tree exhausts without a `WINNER`, the
   `alpha_grading_order` API addition ships as a default-`0`
   no-op, the three strict gate docstrings refresh with the
   v3.5.2 measurement, and Phase 4.fem.eig.3.5.3 (rotated PML
   or fem-eig-006-specific tuning) is queued. We do not weaken
   any gate tolerance.

CPU-only, single-threaded, FP64 complex. No GPU. No new
dependencies. Same `faer::sparse::FaerLuSolver<Complex64>`
surface; the depth-graded `alpha_alpha(d)` lives inside the
per-quadrature-point `s_for(axis, d_alpha)` closure of
`pml_stretching_lambda` — `Λ(ω)` stays diagonal in the global
frame so complex-LDLᵀ is preserved (ADR-0043 decision (4) and
ADR-0044 decision both carried).

## Consequences

- **`PmlConfig` grows by one field**: `alpha_grading_order:
  usize`, default `0`. The default value preserves v3.5.1
  behaviour bit-for-bit. Existing call sites that set every
  field explicitly must add `alpha_grading_order: 0` (or use
  `..Default::default()`); the in-tree call sites are
  `crates/yee-validation/examples/cfs_pml_grading_sweep.rs` and
  the unit tests in `crates/yee-fem/tests/pml_open_boundary_assembly.rs`,
  both updated as part of the S1/S2 ladder.
- **`PmlConfig::default()` ships with new `(kappa_max, m,
  thickness_cells, alpha_grading_order)` quadruple** selected by
  the S3 H4 ablation analysis (or stays at the QQQQQQQQQ
  escape-hatch values on §6 exhaustion). The choice is annotated
  with the winning sweep CSV row number and the post-v3.5.2
  measurement in a comment block immediately above the existing
  v3.5.1 comment block at `crates/yee-fem/src/open_boundary.rs:319-346`;
  the v3.5.1 block stays in place as historical record.
- **fem-eig-003 + fem-eig-006 strict gates flip from `#[ignore]`'d
  to CI-default** on the common-case retire path, without
  weakening any tolerance. Three gates total:
  `fem_eig_003_strict_absorption_floor_gate`,
  `fem_eig_003_strict_passive_bound_continuum_limit`,
  `fem_eig_006_magnitude_bounded`.
- **`crates/yee-validation/examples/cfs_pml_grading_sweep.rs`
  grows by one grid stage** (H4: 18 rows). The pre-v3.5.2 CSV
  rows continue to record `alpha_grading_order = 0` for
  backward-compatibility with v3.5.1 sweep snapshots.
- **PML cell-count growth.** `thickness = 16` extended meshes
  grow from ~72 k to ~195 k tets on fem-eig-003 (~2.7× v3.5.1).
  Per-sweep wall-time grows from ~140 s to ~380 s `--release`.
  The full 18-row design exploration sweep still fits inside
  ~2 h; CI-default ungated test runs use only the winning
  configuration (~5 h for the full 50-point fem-eig-003 sweep).
  If CI wall-time becomes a constraint, downsampling to a
  10-point sweep recovers ~75 min and the worst-case strict-band
  criterion remains robust.
- **`alpha_alpha(d)` causality-canary risk** (§7 a in the spec).
  Ramping `alpha_alpha(D) → 0` at the outer truncation removes
  the Kuzuoglu-Mittra 1996 causality canary from the outermost
  cells. Mitigation: the existing `denom.norm_sqr() <= f64::MIN_POSITIVE`
  guard in `pml_stretching_lambda::s_for` covers the
  `alpha_alpha(D) = 0 ∧ ω → 0` joint limit — the closure
  returns `Complex64::new(kappa_d, 0.0)` instead of NaN. The S1
  unit tests assert this behaviour.

## References

- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading-design.md`
  — design spec.
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading.md`
  — S1-S5 implementation plan.
- ADR-0044 — Phase 4.fem.eig.3.5.1 scope (this ADR's parent;
  §Consequences and spec §7 (b) queue this ADR).
- ADR-0043 — Phase 4.fem.eig.3.5 scope (grandparent; CFS-PML
  implementation).
- ADR-0042 — Phase 4.fem.eig.3 scope (great-grandparent).
- Berenger, J.-P., "Numerical reflection from FDTD-PMLs: a
  comparison of the split PML with the unsplit and CFS PMLs,"
  *IEEE Transactions on Antennas and Propagation* 50(3) (March
  2002), pp. 258-265, DOI 10.1109/8.999615 — §VI
  `alpha_alpha(d) = alpha_max · (1 - d/D)^n` canonical sweep
  with `n ∈ {1, 2, 3}`.
- Roden, J. A. and Gedney, S. D., "Convolutional PML (CPML)",
  *IEEE MWCL* 10(5) (May 2000) — `sigma_max` / `kappa_max`
  Table-I defaults inherited via v3.5 / v3.5.1.
- Kuzuoglu-Mittra 1996 *IEEE MWCL* 6(12) — CFS causality
  modification `alpha_alpha > 0` motivating the §7 (a)
  causality-canary risk.
- `crates/yee-fem/src/open_boundary.rs:319-346` — QQQQQQQQQ R3
  comment block recording the `[-58.13, -35.45] dB` H3
  most-aggressive probe measurement that triggers this ADR.
- CLAUDE.md §3, §4, §10.
