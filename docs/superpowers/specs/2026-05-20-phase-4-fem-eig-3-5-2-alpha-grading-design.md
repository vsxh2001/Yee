# Phase 4.fem.eig.3.5.2 — CFS-PML `alpha_alpha(d)` grading + extended H3 thickness

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3.5.2 (Complex-Frequency-Shifted PML `alpha_alpha(d)`
polynomial grading + extended H3 thickness ablation; closes the
5 dB gap to the spec §6 `[-60, -40] dB` strict-absorption-band that
the QQQQQQQQQ Phase 4.fem.eig.3.5.1 R3 retune left open).
**Depends on:** Phase 4.fem.eig.3.5.1 (QQQQQQQQQ R1-R5 shipped: per-axis
`h_α` resolver, `PmlMeshMeta` carrier, R2 ablation binary
`crates/yee-validation/examples/cfs_pml_grading_sweep.rs`, escape-hatch
defaults `(κ_max=5, m=3, thickness=6)` unchanged from OOOOOOOOO).
**Blocks:** retirement of `#[ignore]` on
`fem_eig_003_strict_absorption_floor_gate`,
`fem_eig_003_strict_passive_bound_continuum_limit`, and
`fem_eig_006_magnitude_bounded` — the three gates QQQQQQQQQ R4
left in `#[ignore]` purgatory because the v3.5.1 ablation grid
exhausted ~5 dB short of the strict-band lower bound.

## 1. Goal

Retire the fem-eig-003 strict `[-60, -40] dB` absorption-floor band
and the fem-eig-006 `|S_{11}(30 GHz)| < 0.1` magnitude gate by
adding **two independent knobs** to the QQQQQQQQQ ablation grid:

1. **CFS frequency-shift `alpha_alpha(d)` polynomial grading** per
   Berenger 2002 §VI — replace the constant `alpha_max` carried by
   `PmlConfig` with a depth-dependent
   `alpha_alpha(d) = alpha_max · (1 - d/D)^n` ramping from
   `alpha_max` at the cavity-PML interface (`d = 0`) to `0` at the
   outer truncation surface (`d = D`). Berenger 2002 §VI reports
   ~10-20 dB improvement at the inner boundary on the canonical 2D
   waveguide benchmark when graded `alpha_alpha(d)` replaces the
   constant `alpha_max`.
2. **Extended H3 thickness sweep** beyond the QQQQQQQQQ R2 grid —
   evaluate `thickness_cells ∈ {12, 14, 16}` at the
   `kappa_max = 2, m ∈ {3, 4}` cells of the H3 sub-grid that the
   QQQQQQQQQ probe identified as the local optimum
   (`[-58.13, -35.45] dB` worst-case probe). The QQQQQQQQQ grid
   capped at `thickness = 10`; the extended sweep tests whether
   thickness alone is the binding constraint.

Pick the first parameter combination that retires both strict
gates (`fem_eig_003_s11_max_db < -40 dB` *and*
`fem_eig_006_s11_mag < 0.1`); ship it as the new
`PmlConfig::default()` triple plus the new
`alpha_grading_order: usize` field. If both knobs land but neither
retires both fixtures alone, ship the joint combination (extended
thickness *and* alpha grading) and update the comment block above
`impl Default for PmlConfig` to record the joint measurement.

## 2. Background

### 2.1 QQQQQQQQQ R3 measurement summary

Track QQQQQQQQQ shipped Phase 4.fem.eig.3.5.1 R1-R5 against ADR-0044
on base SHA `56386f1`. R1 introduced the per-axis `h_α` resolver
and the `PmlMeshMeta` carrier; R2 authored the
`cfs_pml_grading_sweep` example binary; R3 captured the partial
ablation table.

The R3 comment block above `impl Default for PmlConfig` in
`crates/yee-fem/src/open_boundary.rs:319-346` records the
measurement:

- **R1 per-axis `h_α` alone** moved fem-eig-003 worst-case from
  the OOOOOOOOO baseline of `-7.48 dB` to `-21.74 dB` — a ~14 dB
  improvement.
- **H2 `kappa_max` sweep** at `(m=3, thickness=6)` showed all
  values in `{1, 1.5, 2, 3, 5}` cluster within ~1 dB at
  `~ -22 dB` worst-case. `kappa_max` is **not** the binding
  constraint at this mesh resolution (consistent with
  Berenger 2002 §V).
- **H3 most-aggressive probe** at `(kappa=2, m=4, thickness=10)`
  reaches `[-58.13, -35.45] dB`. The band minimum is *already
  inside* the `[-60, -40] dB` target window; only the worst-case
  miss-distance of 5 dB above `-40` blocks the strict gate.

The 5 dB shortfall after the v3.5.1 grid is exhausted exactly
matches the magnitude that Berenger 2002 §VI attributes to the
constant-`alpha_max` artefact at the inner boundary — the §VI
canonical sweep shows ~5-10 dB worst-case improvement when the
constant `alpha_max` is replaced by a polynomial ramp
`alpha_alpha(d) = alpha_max · (1 - d/D)^n`.

### 2.2 Berenger 2002 §VI alpha grading

Berenger, J.-P., "Numerical reflection from FDTD-PMLs", *IEEE
TAP* 50(3) (2002), DOI 10.1109/8.999615, §VI ("CFS-PML parameter
optimisation") observes that the worst-case reflection floor of
the CFS-PML is dominated near the **inner** cavity-PML boundary
by spurious reflections off the discontinuity between the
`alpha_alpha = alpha_max` interior boundary value and the
constant-`alpha_max` interior of the PML shell. Grading
`alpha_alpha(d)` from `alpha_max` at the inner boundary to `0`
at the outer truncation surface smooths the discontinuity and
~10-20 dB improvement at low frequencies near the boundary, with
a corresponding ~5-10 dB improvement in the worst-case swept
band.

The grading polynomial mirrors the existing `sigma_alpha(d)`
polynomial in shape but with **opposite sign on the depth axis**:
`sigma_alpha(d)` rises from `0` at the inner boundary to
`sigma_max` at the outer truncation, while `alpha_alpha(d)`
falls from `alpha_max` at the inner boundary to `0` at the outer
truncation. Berenger 2002 §VI canonical sweep uses
`alpha_alpha(d) = alpha_max · (1 - d/D)^n` with `n ∈ {1, 2, 3}`;
`n = 1` is the linear ramp Berenger 2002 §VI defaults to and
the most parsimonious choice.

### 2.3 ADR-0044 §risks (b)

ADR-0044 (Phase 4.fem.eig.3.5.1 scope) §risks (b) explicitly
flagged this gap:

> **(b) `α_α(d)` grading is the other candidate cause.** Berenger
> 2002 §VI shows that for canonical evanescent-mode problems the
> CFS frequency-shift parameter `α_α` is also profitably **graded**
> as `α(d) = α_max · (1 - d/D)^m_α` rather than constant. The
> OOOOOOOOO defaults use a constant `α_max`. This ablation is
> *not* in v3.5.1 scope; deferred to v3.5.2 if H1+H2+H3 retune
> still misses by > 5 dB.

The QQQQQQQQQ R3 measurement triggered exactly this deferral path:
the H3 most-aggressive probe missed `-40 dB` by 5 dB on the
worst-case point. This spec is the v3.5.2 follow-on.

### 2.4 Why two knobs in one phase

`alpha_alpha(d)` and extended thickness are conceptually
independent — alpha grading targets the inner-boundary
discontinuity, thickness targets the absorption-budget magnitude.
Each could be the binding constraint, and the QQQQQQQQQ probe
data is insufficient to decide a priori. Running both ablations
in the same phase (cost: 18 configurations) gives a clean
decision: whichever knob retires the gate first ships; if neither
alone retires but the joint configuration does, both ship.

## 3. Mathematical formulation

### 3.1 Current constant-`alpha_max` formulation (v3.5.1)

The QQQQQQQQQ `pml_stretching_lambda` evaluator in
`crates/yee-fem/src/open_boundary.rs:2444-2506` computes the
per-axis stretching factor as

```text
s_α(d_α, ω) = κ_α(d_α) + σ_α(d_α) / (α_α + j ω ε_0)
```

where

```text
σ_α(d_α) = σ_α_max · (d_α / D_α)^m
κ_α(d_α) = 1 + (κ_max - 1) · (d_α / D_α)^m
α_α      = α_max                       (constant; v3.5)
D_α      = thickness_cells · h_α
```

The constant `alpha_max` enters every `s_α` evaluation identically,
regardless of the depth `d_α` of the quadrature point inside the
PML shell.

### 3.2 New `alpha_alpha(d)` polynomial grading (v3.5.2)

Replace the constant `α_α` with a depth-dependent
`alpha_alpha(d_α)` polynomial:

```text
alpha_alpha(d_α) = alpha_max · (1 - d_α / D_α)^n
```

where `n` is the new `alpha_grading_order: usize` field on
`PmlConfig`. With `n = 0` (the default), `(1 - d/D)^0 = 1` and
the formulation collapses bit-for-bit to the v3.5.1 constant
`alpha_max`. With `n ≥ 1`, the alpha ramps from `alpha_max` at
the cavity-PML interface (`d = 0`) to `0` at the outer truncation
surface (`d = D`). `n = 1` is the linear Berenger 2002 §VI
default; `n = 2, 3` are the higher-order polynomial variants
covered by the §VI sweep.

The updated `s_α(d_α, ω)` formula becomes:

```text
s_α(d_α, ω) = κ_α(d_α) + σ_α(d_α) / (alpha_alpha(d_α) + j ω ε_0)
```

Only the denominator changes; the `σ_α(d_α)` and `κ_α(d_α)`
polynomials are unaffected.

### 3.3 Backward-compatibility default

The new `alpha_grading_order: usize` field defaults to `0`.
`alpha_alpha(d) ≡ alpha_max` reduces to the v3.5.1 constant
formulation bit-for-bit. The v3.5.1 sweep CSV rows and the
QQQQQQQQQ R3 measurement block remain reproducible without
re-running the ablation.

## 4. Public API

`PmlConfig` gains one new field:

```rust
pub struct PmlConfig {
    pub thickness_cells: usize,
    pub sigma_max: f64,
    pub alpha_max: f64,
    pub kappa_max: f64,
    pub m: usize,
    /// Phase 4.fem.eig.3.5.2: `alpha_alpha(d)` polynomial grading
    /// order. `0` (default) recovers v3.5.1 constant `alpha_max`
    /// bit-for-bit. `1, 2, 3` enable Berenger 2002 §VI ramp
    /// `alpha_alpha(d) = alpha_max · (1 - d/D)^alpha_grading_order`.
    pub alpha_grading_order: usize,
}
```

`PmlConfig::default()` keeps `alpha_grading_order = 0` until the
§5 ablation picks a winner. The default updates atomically with
the winning `(kappa_max, m, thickness_cells, alpha_grading_order)`
quadruple in the R3 / S3 comment block.

`ResolvedPmlConfig` (the private carrier consumed by
`pml_stretching_lambda`) gains the same `alpha_grading_order`
field, copied verbatim from the public `PmlConfig` by
`PmlConfig::resolved`. The internal `s_for(axis, d_α)` closure
inside `pml_stretching_lambda` evaluates
`alpha_alpha = alpha_max · (1 - ratio)^alpha_grading_order` in
line with the v3.5.1 `kappa_d` / `sigma_d` evaluations and
substitutes it for the constant `cfg.alpha_max` in the
denominator.

No new public types; no signature change on `with_cfs_pml`. The
Python `yee.fem.solve_open_cavity` binding's `pml_config` kwarg
gains the optional `alpha_grading_order` field with default `0`.

## 5. Ablation grid

Fix `kappa_max = 2` per the QQQQQQQQQ H2 finding that kappa
clusters within 1 dB across `{1, 1.5, 2, 3, 5}` at the v3.5.1
grid resolution. Sweep the three independent axes:

| axis                  | values                |
|-----------------------|-----------------------|
| `m`                   | `{3, 4}`              |
| `thickness_cells`     | `{12, 14, 16}`        |
| `alpha_grading_order` | `{0, 1, 2}`           |

`2 × 3 × 3 = 18 configurations`.

Per spec §6 stopping rule, run fem-eig-003 first on each row
(cheaper at ~70 k extended tets; ~140 s `--release` wall-time per
sweep). Only run fem-eig-006 if fem-eig-003 worst-case
`s11_max_db < -40 dB`. On the first row where *both* fixtures
retire, emit a `WINNER,...` CSV row and exit.

Worst-case wall-time: `18 × 140 s ≈ 42 min` for fem-eig-003 alone
(~30 s additional per fem-eig-006 follow-on row that triggers).
The grid stays well inside the v3.5.1 "1-h design exploration"
budget envelope.

### 5.1 Row layout

The sweep binary `crates/yee-validation/examples/cfs_pml_grading_sweep.rs`
gains a v3.5.2 grid stage after the existing v3.5.1 H1/H2/H3
grid:

```text
H4,kappa=2,m=3,thickness=12,alpha_grading_order=0,...
H4,kappa=2,m=3,thickness=12,alpha_grading_order=1,...
H4,kappa=2,m=3,thickness=12,alpha_grading_order=2,...
H4,kappa=2,m=3,thickness=14,alpha_grading_order=0,...
... (18 rows; pattern continues across m × thickness × alpha_grading_order)
```

The new CSV column `alpha_grading_order` appends to the v3.5.1
header. Pre-v3.5.2 CSV rows record `alpha_grading_order = 0`
(legacy constant `alpha_max`) for reproducibility.

## 6. Stopping rule

Per-row decision flow:

1. Run fem-eig-003 with the row's `(m, thickness_cells, alpha_grading_order)`
   triple (kappa fixed at 2). Record
   `(s11_min_db, s11_max_db, runtime_s)`.
2. If `s11_max_db < -40 dB` (strict band retire), run fem-eig-006
   at the same config. Record `(s11_mag, runtime_s)`.
3. If `s11_mag < 0.1` (magnitude gate retire), emit `WINNER` and
   exit. This is the new `PmlConfig::default()`.
4. Otherwise continue to the next row.

If the full 18-row grid exhausts without a `WINNER`:

- If at least one row retires fem-eig-003 but no row retires
  fem-eig-006: ship the smallest-product
  `(m × thickness_cells)` row that retires fem-eig-003 as
  per-fixture override for fem-eig-003 and queue Phase
  4.fem.eig.3.5.3 for fem-eig-006-specific tuning.
- If no row retires either gate: revert to the QQQQQQQQQ R3
  escape-hatch path. Ship the `alpha_grading_order` field as a
  default-`0` no-op API addition; refresh the
  `impl Default for PmlConfig` comment block with the v3.5.2
  worst-case measurement; queue Phase 4.fem.eig.3.5.3 for further
  work. Do **not** weaken any gate.

## 7. Risks

- **(a) `alpha_alpha(d)` grading destabilises the PML at low
  frequencies.** The CFS `alpha_alpha > 0` term is a causality
  canary — ramping it to `0` at the outer truncation surface
  removes the canary from the outermost cells. If the band-min
  frequency `f_min` is low enough that `j ω ε_0` is small
  compared to the now-zero `alpha_alpha(D)`, the
  `s_α(d, ω) = κ + σ / (0 + j ω ε_0)` evaluation may produce
  `|s_α| → ∞` and NaN-poison the assembled stiffness. Mitigation:
  guard the `alpha_alpha(d)` evaluator with a `1e-12` floor
  identical to the v3.5.1 `kappa_d > 0` guard; if the causality
  canary fails at the band-min frequency (i.e. assembled
  stiffness contains NaN/Inf), the sweep binary records the row
  with `s11_max_db = NaN` and skips fem-eig-006; the R3 decision
  ignores NaN rows.

- **(b) Extended thickness alone may not retire if alpha is the
  binding constraint.** The H3 finding `[-58.13, -35.45] dB` is
  on `thickness = 10`. Linear extrapolation from the v3.5.1 grid
  (`thickness ∈ {6, 8, 10}` worst-case `~ -22, -28, -35 dB`)
  suggests `thickness = 16` would land near `-46 dB`
  worst-case — inside the strict band but only barely. If the
  binding constraint is the inner-boundary alpha discontinuity
  (Berenger 2002 §VI), thickness alone may saturate near
  `-40 dB` without ever crossing. Mitigation: the §5 grid covers
  both axes simultaneously so the data resolves the question
  empirically.

- **(c) Joint configuration ships both knobs; rollback is harder.**
  If the §5 winner is the joint `(thickness ≥ 12, alpha_grading_order ≥ 1)`
  configuration, both knobs land as new defaults simultaneously.
  Rollback to the v3.5.1 escape-hatch defaults requires changing
  two fields. Mitigation: the `alpha_grading_order` default of
  `0` makes the v3.5.1 behaviour the explicit "off" state for the
  alpha grading knob; reverting just the thickness back to `6`
  recovers the v3.5.1 measurement exactly when paired with
  `alpha_grading_order = 0`.

- **(d) PML cell-count growth.** `thickness = 16` doubles the
  PML-shell extended-tet count from the v3.5.1 baseline
  (`thickness = 6`) by ~2.7×; the fem-eig-003 extended mesh grows
  from ~72 k to ~195 k tets. Per-sweep wall-time grows
  proportionally — from ~140 s to ~380 s `--release`. The full
  18-row sweep still fits inside ~2 h. CI-default ungated tests
  only run the **winning** configuration, so the CI wall-time
  impact is one ~380 s per fem-eig-003 frequency sweep (~50
  points = ~5 hours total). If this is judged unacceptable for
  CI, the strict gates can be downsampled to a 10-point sweep
  (~75 min); the spec §6 retire criterion uses the worst-case
  over the sweep, which is robust to downsampling.

## 8. Lane

Spec file:
`docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading-design.md`

Implementation lane (declared here for the S1-S5 plan):

- `crates/yee-fem/src/open_boundary.rs` — add `alpha_grading_order`
  field to `PmlConfig` + `ResolvedPmlConfig`; update
  `PmlConfig::resolved` + `pml_stretching_lambda::s_for` to
  evaluate `alpha_alpha(d) = alpha_max · (1 - d/D)^alpha_grading_order`
  (collapses to constant when `alpha_grading_order == 0`).
- `crates/yee-fem/src/lib.rs` — no re-export change (the field
  rides on the existing `PmlConfig` re-export).
- `crates/yee-fem/tests/pml_open_boundary_assembly.rs` — extend
  to assert `alpha_grading_order = 0` produces bit-for-bit
  identical `Λ(ω)` to the v3.5.1 path.
- `crates/yee-validation/examples/cfs_pml_grading_sweep.rs` —
  append the H4 grid stage (extended thickness + alpha grading);
  add the `alpha_grading_order` CSV column.
- `crates/yee-validation/src/lib.rs` —
  `run_fem_eig_003_*_with_config` + `run_fem_eig_006_*_with_config`
  already pass-through `PmlConfig`; no signature change needed
  since `alpha_grading_order` lives inside `PmlConfig`.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` —
  R4 / S4 un-ignore (or v3.5.3 escape-hatch refresh).
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs` —
  R4 / S4 un-ignore (or v3.5.3 escape-hatch refresh).
- `docs/src/tutorials/07-fem-open-cavity.md` — note the new
  `alpha_grading_order` kwarg + the v3.5.2 measurement.
- `ROADMAP.md` — Phase 4.fem.eig.3.5.2 entry from planned to
  shipped.

Out of lane: `yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`,
`yee-cuda`, `yee-plotters`, `yee-fdtd`, `yee-surrogate`.

## 9. References

- Berenger, J.-P., "Numerical reflection from FDTD-PMLs: a
  comparison of the split PML with the unsplit and CFS PMLs,"
  *IEEE Transactions on Antennas and Propagation* 50(3) (March
  2002), pp. 258-265. DOI 10.1109/8.999615. §VI canonical
  `alpha_alpha(d) = alpha_max · (1 - d/D)^n` sweep with
  `n ∈ {1, 2, 3}`.
- Roden, J. A. and Gedney, S. D., "Convolutional PML (CPML)",
  *IEEE MWCL* 10(5) (May 2000) — Table-I `kappa_max` /
  `sigma_max` defaults inherited via v3.5 / v3.5.1.
- Kuzuoglu, M. and Mittra, R., 1996 *IEEE MWCL* 6(12) — CFS
  causality modification `alpha_alpha > 0` motivating the
  causality canary in §7 (a).
- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune-design.md`
  — v3.5.1 parent spec.
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune.md`
  — v3.5.1 R1-R5 plan; QQQQQQQQQ shipped.
- `docs/src/decisions/0044-phase-4-fem-eig-3-5-1-grading-retune.md`
  — v3.5.1 ADR (§Consequences flagged `alpha_alpha(d)` as
  remaining knob).
- `docs/src/decisions/0045-phase-4-fem-eig-3-5-2-alpha-grading.md`
  — this spec's scope ADR.
- `crates/yee-fem/src/open_boundary.rs:319-346` — QQQQQQQQQ R3
  comment block recording the `[-58.13, -35.45] dB` H3
  most-aggressive probe measurement.
- CLAUDE.md §3, §4, §10.
