# Phase 4.fem.eig.3.5.2 — `alpha_alpha(d)` grading + extended H3 thickness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans`
> to drive this plan step-by-step.

**Companion spec:**
`docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading-design.md`
**Companion ADR:**
`docs/src/decisions/0045-phase-4-fem-eig-3-5-2-alpha-grading.md`
**Base SHA:** `5ec8e90` (main HEAD; QQQQQQQQQ Phase 4.fem.eig.3.5.1
R1-R5 shipped).
**Target phase:** 4.fem.eig.3.5.2 only. 4.fem.eig.3.5.3 (rotated PML,
fem-eig-006-specific tuning) and 4.fem.eig.4 (FEM-BEM hybrid)
explicitly deferred.
**Tech-stack additions:** none. Same `yee-fem` / `yee-validation`
surface area; one new field on `PmlConfig` and an H4 grid stage
appended to the existing `cfs_pml_grading_sweep` example binary.

---

## Goal

Close the 5 dB gap between the QQQQQQQQQ H3 most-aggressive probe
(`[-58.13, -35.45] dB` worst-case at `kappa=2, m=4, thickness=10`)
and the spec §6 strict `[-60, -40] dB` absorption band by adding
two independent knobs to the v3.5.1 ablation grid:

1. `alpha_alpha(d) = alpha_max · (1 - d/D)^n` polynomial grading
   per Berenger 2002 §VI (~5-10 dB worst-case improvement
   expected).
2. Extended thickness sweep `thickness_cells ∈ {12, 14, 16}`
   beyond the v3.5.1 cap of 10 (~5-10 dB worst-case improvement
   expected from absorption-budget growth).

Pick the first (m, thickness, alpha_grading_order) configuration
that retires both `fem-eig-003_strict_absorption_floor_gate` and
`fem-eig-006_magnitude_bounded`, ship it as the new
`PmlConfig::default()` quadruple, and un-ignore the three strict
gates the QQQQQQQQQ R4 left in `#[ignore]` purgatory.

Five-step ladder S1-S5 lands in a single merge train:

1. **S1** — add `alpha_grading_order: usize` field to
   `PmlConfig` + `ResolvedPmlConfig`; thread
   `alpha_alpha(d) = alpha_max · (1 - d/D)^alpha_grading_order`
   through `pml_stretching_lambda::s_for`.
2. **S2** — extend `cfs_pml_grading_sweep` example binary with
   the H4 grid stage (kappa=2 × m ∈ {3, 4} × thickness ∈ {12, 14,
   16} × alpha_grading_order ∈ {0, 1, 2} = 18 configurations).
3. **S3** — run the full sweep; pick winning quadruple; update
   `PmlConfig::default()`.
4. **S4** — un-ignore the three strict gates; verify both
   fixtures pass under the new defaults. If the sweep exhausts
   without a winner, fall through to the v3.5.3 escape-hatch
   path: leave `#[ignore]`'s in place, refresh measurement
   docstrings, queue v3.5.3.
5. **S5** — tutorial note in `docs/src/tutorials/07-fem-open-cavity.md`
   on the new `alpha_grading_order` kwarg + ROADMAP refresh.

CPU-only, single-threaded, scalar FP64 complex. No GPU. No new
dependencies. Same execution model as v3.5.1.

## Pre-flight

Before Step S1, confirm at base SHA `5ec8e90`:

1. `crates/yee-fem/src/open_boundary.rs` exposes
   `AbcOrder::CfsPml(PmlConfig)` and the `pml_stretching_lambda`
   per-axis `s_for(axis, d_alpha)` closure from v3.5.1 R1.
2. `PmlConfig` carries fields `thickness_cells, sigma_max, alpha_max,
   kappa_max, m` (five fields). The new `alpha_grading_order` is
   the sixth.
3. `crates/yee-validation/examples/cfs_pml_grading_sweep.rs` exists
   from v3.5.1 R2 and emits CSV with `hypothesis, kappa_max, m,
   thickness_cells, fem_eig_003_s11_min_db, fem_eig_003_s11_max_db,
   fem_eig_006_s11_mag, runtimes`.
4. `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
   carries the two `#[ignore]`'d strict gates with QQQQQQQQQ R4
   measurement docstrings recording the `[-58.13, -35.45] dB` H3
   probe baseline.
5. `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
   carries the `#[ignore]`'d `fem_eig_006_magnitude_bounded`
   gate with QQQQQQQQQ R4 measurement docstring.

If (1)-(5) blocks, escape-hatch per CLAUDE.md §5 >15-min rule and
surface as a base-SHA drift finding; do **not** weaken any gate.

## File structure

| File | Action | Step | Responsibility |
|------|--------|------|----------------|
| `crates/yee-fem/src/open_boundary.rs` | Modify | S1, S3 | Add `alpha_grading_order` field; update `pml_stretching_lambda::s_for`; new `PmlConfig::default()` after S3 analysis. |
| `crates/yee-fem/tests/pml_open_boundary_assembly.rs` | Modify | S1 | `alpha_grading_order = 0` ≡ v3.5.1 path assertion. |
| `crates/yee-validation/examples/cfs_pml_grading_sweep.rs` | Modify | S2 | Append H4 grid stage; new `alpha_grading_order` CSV column. |
| `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` | Modify | S4 | Remove both `#[ignore]` attributes (common case) or refresh docstrings (escape-hatch case). |
| `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs` | Modify | S4 | Remove the `#[ignore]` attribute (common case) or refresh docstring (escape-hatch case). |
| `docs/src/tutorials/07-fem-open-cavity.md` | Modify | S5 | Note the new `alpha_grading_order` kwarg + measurement. |
| `ROADMAP.md` | Modify | S5 | Phase 4.fem.eig.3.5.2 entry from planned to shipped. |

## Step S1 — `alpha_grading_order` field + `alpha_alpha(d)` closure

**Lane:** `crates/yee-fem/src/open_boundary.rs`,
`crates/yee-fem/tests/pml_open_boundary_assembly.rs`.

Add the new field to `PmlConfig`:

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
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

Update the same field on `ResolvedPmlConfig` (the private carrier
near `crates/yee-fem/src/open_boundary.rs:420`). `PmlConfig::resolved`
copies `alpha_grading_order` verbatim into the carrier; no other
resolver logic changes.

Inside `pml_stretching_lambda::s_for(axis, d_alpha)` (near line
2459 in the v3.5.1 file), replace the constant `cfg.alpha_max`
with a depth-dependent `alpha_alpha` evaluator:

```rust
let s_for = |axis: usize, d_alpha: f64| -> Complex64 {
    let h_alpha = cfg.h_per_axis[axis].max(1.0e-12);
    let d_max = thickness * h_alpha;
    if d_alpha <= 0.0 || d_max <= 0.0 {
        return Complex64::new(1.0, 0.0);
    }
    let ratio = (d_alpha / d_max).clamp(0.0, 1.0);
    let pow_sigma = ratio.powf(m);
    let sigma_d = cfg.sigma_max_per_axis[axis] * pow_sigma;
    let kappa_d = 1.0 + (cfg.kappa_max - 1.0) * pow_sigma;
    // Phase 4.fem.eig.3.5.2: alpha_alpha(d) polynomial grading
    // per Berenger 2002 §VI. With alpha_grading_order = 0 the
    // pow_alpha = 1 collapses bit-for-bit to v3.5.1 constant alpha.
    let pow_alpha = if cfg.alpha_grading_order == 0 {
        1.0
    } else {
        (1.0 - ratio).powi(cfg.alpha_grading_order as i32)
    };
    let alpha_alpha = cfg.alpha_max * pow_alpha;
    let denom = Complex64::new(alpha_alpha, omega * yee_core::units::EPS0);
    if denom.norm_sqr() <= f64::MIN_POSITIVE {
        return Complex64::new(kappa_d, 0.0);
    }
    Complex64::new(kappa_d, 0.0) + Complex64::new(sigma_d, 0.0) / denom
};
```

The `denom.norm_sqr() <= f64::MIN_POSITIVE` guard already covers
the v3.5.2 §7 (a) causality-canary risk: if `alpha_alpha(D) = 0`
and `omega → 0` simultaneously, the closure returns
`Complex64::new(kappa_d, 0.0)` (no NaN poison; the assembled
stiffness stays finite).

**Pattern file:** mirror the v3.5.1 `pml_stretching_lambda`
evaluator directly — only the `pow_alpha` block and the
`alpha_alpha` substitution are new; everything else is verbatim.

**Test update in `pml_open_boundary_assembly.rs`:**

- `alpha_grading_order_zero_matches_v3_5_1` — with
  `alpha_grading_order = 0`, the per-tet assembled blocks at
  every quadrature point match the v3.5.1 path bit-for-bit
  (Frobenius difference < 1e-15).
- `alpha_grading_order_one_smooths_inner_boundary` — with
  `alpha_grading_order = 1` and a single PML cell, the
  `alpha_alpha` at the inner-boundary depth `d = 0` equals
  `alpha_max` exactly, and at the outer truncation `d = D`
  equals `0.0` exactly.

**DoD S1.**
- `cargo check -p yee-fem` exits 0.
- `cargo test -p yee-fem --test pml_open_boundary_assembly` exits 0.
- `grep -q 'alpha_grading_order' crates/yee-fem/src/open_boundary.rs`
  exit 0.

## Step S2 — H4 grid stage in `cfs_pml_grading_sweep`

**Lane:** `crates/yee-validation/examples/cfs_pml_grading_sweep.rs`.

Append a fourth grid stage after the existing v3.5.1 H1/H2/H3
stages. The new stage fixes `kappa_max = 2` and sweeps three
axes:

```rust
// Phase 4.fem.eig.3.5.2 H4: alpha_grading_order + extended thickness.
for &m in &[3_usize, 4] {
    for &thickness in &[12_usize, 14, 16] {
        for &alpha_grading_order in &[0_usize, 1, 2] {
            grid.push(Configuration {
                hypothesis: "H4",
                kappa_max: 2.0,
                m,
                thickness_cells: thickness,
                alpha_grading_order,
            });
        }
    }
}
```

`Configuration` gains the `alpha_grading_order: usize` field;
`as_pml_config()` threads it through.

Append the `alpha_grading_order` column to the CSV header. The
existing v3.5.1 H1/H2/H3 rows record `alpha_grading_order = 0`
(legacy constant `alpha_max`).

Implements the spec §6 stopping rule: run fem-eig-003 first per
row; only run fem-eig-006 if fem-eig-003 worst-case
`s11_max_db < -40 dB`. On the first row where *both* fixtures
retire, emit a final `WINNER,...` row and exit.

**Pattern file:** the existing v3.5.1
`crates/yee-validation/examples/cfs_pml_grading_sweep.rs::build_grid`
fn is the direct template; the H4 stage is a near-copy of the H3
stage with one additional nested loop.

**Smoke test:** the sweep is **not** a CI gate — design
exploration only. The S2 DoD only requires that the binary
compiles and a single-configuration dry run emits one CSV row
with the new column.

**DoD S2.**
- `cargo build -p yee-validation --example cfs_pml_grading_sweep
  --release` exits 0.
- `cargo run -p yee-validation --example cfs_pml_grading_sweep
  --release -- --dry-run` exits 0 and writes one CSV row to
  stdout with the `alpha_grading_order` column populated.
- `grep -q 'alpha_grading_order' crates/yee-validation/examples/cfs_pml_grading_sweep.rs`
  exit 0.

## Step S3 — analyse CSV; pick defaults; update `PmlConfig::default`

**Lane:** `crates/yee-fem/src/open_boundary.rs`.

Run the full v3.5.2 H4 sweep (18 configurations; worst-case
wall-time ~42 min `--release` for fem-eig-003-only, plus ~30 s
per fem-eig-006 trigger). Capture stdout to
`/tmp/phase-4-fem-eig-3-5-2-sweep.csv` (not committed; referenced
in the S5 tutorial).

Apply the spec §6 stopping rule:

1. Pick the first row where `s11_max_db < -40 dB` and
   `s11_mag < 0.1` — this is `WINNER`.
2. If no such row, pick the smallest `(m × thickness)` row that
   retires fem-eig-003; ship it as fem-eig-003-specific override
   and queue v3.5.3 for fem-eig-006.
3. If no row retires either gate, fall through to the
   escape-hatch path: ship `alpha_grading_order` as a default-`0`
   no-op API addition; queue v3.5.3.

Update `PmlConfig::default()` in `open_boundary.rs` with the
winning quadruple `(kappa_max, m, thickness_cells, alpha_grading_order)`.
Annotate the choice in a `// Phase 4.fem.eig.3.5.2 retune
(2026-05-20, sweep CSV row N):` comment block immediately above
the existing v3.5.1 comment block (which stays in place as
historical record).

**Decision-tree exhaustion (escape-hatch):** if no row retires
both fixtures, leave `PmlConfig::default()` at the QQQQQQQQQ
escape-hatch values `(kappa_max=5, m=3, thickness=6, alpha_grading_order=0)`
and add an additional comment block recording the v3.5.2
worst-case measurement. Jump to S4-alternative: do **not**
un-ignore the strict gates; instead update their docstrings with
the new v3.5.2 baseline; queue v3.5.3.

**DoD S3.**
- The CSV exists locally (not committed) and shows the 18
  ablation rows + WINNER (or 18 rows + EXHAUSTED if escape-hatch).
- `grep -q 'Phase 4.fem.eig.3.5.2 retune'
  crates/yee-fem/src/open_boundary.rs` exit 0.
- `cargo check -p yee-fem` exits 0 with the new defaults.

## Step S4 — un-ignore strict gates + verify pass

**Lane:** `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`,
`crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`.

If S3 ended with the spec §6 retire (the common case):

1. Remove `#[ignore = "..."]` from
   `fem_eig_003_strict_absorption_floor_gate` and
   `fem_eig_003_strict_passive_bound_continuum_limit`.
2. Remove `#[ignore = "..."]` from `fem_eig_006_magnitude_bounded`.
3. Run `cargo test -p yee-validation --release` and verify the
   three tests pass under the new `PmlConfig::default()`.
4. Update the docstrings on each gate to record the post-v3.5.2
   measurement (`Phase 4.fem.eig.3.5.2 status: |S_{11}| ∈ [...] ⇒
   s11_db ∈ [...]; retires the strict gate`).

If S3 ended with decision-tree exhaustion (escape-hatch):

1. **Do not** remove any `#[ignore]`.
2. Update each docstring to record the v3.5.2 baseline (best
   non-retiring row from the H4 grid).
3. ROADMAP refresh (S5) marks Phase 4.fem.eig.3.5.2 as "shipped
   `alpha_grading_order` API addition (default-`0` no-op); strict
   gates remain `#[ignore]`'d, queued for Phase 4.fem.eig.3.5.3".

**DoD S4** (common case).
- `cargo test -p yee-validation --release
  --test fem_eig_003_wr90_stub_abc` exits 0 with all
  smoke + un-ignored strict gates passing.
- `cargo test -p yee-validation --release
  --test fem_eig_006_high_aspect_pml` exits 0 with all
  tests passing.
- `grep -c '#\[ignore' crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  prints `0`.
- `grep -c '#\[ignore' crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
  prints `0`.

## Step S5 — tutorial + ROADMAP refresh

**Lane:** `docs/src/tutorials/07-fem-open-cavity.md`, `ROADMAP.md`.

Add an "Alpha-grading order" subsection to
`docs/src/tutorials/07-fem-open-cavity.md` showing:

- The new `alpha_grading_order` `PmlConfig` field with default
  `0`.
- A worked example overriding the default via the
  `pml_config.alpha_grading_order` Python kwarg for a user
  targeting low-frequency evanescent-mode applications.
- A short table of "knob → effect on `|S_{11}|`" derived from
  the S3 sweep CSV (3-4 representative rows from the H4 grid).

Update `ROADMAP.md` Phase 4.fem.eig.3.5.2 entry from "planned /
optional" to "shipped"; link the un-ignored gate references in
`fem_eig_003_wr90_stub_abc.rs` and `fem_eig_006_high_aspect_pml.rs`
(common case) or the v3.5.3 queue note (escape-hatch case).

**DoD S5.**
- `mdbook build docs/` exits 0.
- `grep -q '4.fem.eig.3.5.2' ROADMAP.md` exit 0.
- `grep -q 'Alpha-grading order\|alpha_grading_order'
  docs/src/tutorials/07-fem-open-cavity.md` exit 0.

## Verification roll-up

After S5:

```bash
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
mdbook build docs/
```

All four must exit 0. The `--release` test invocation is required
for the fem-eig-003 strict gate (the budget is ~140-380 s
`--release` per sweep frequency point depending on the winning
`thickness_cells`; debug builds will time out).

## Out of scope

Explicitly deferred:

- **Rotated / non-Cartesian-aligned PML** — Phase 4.fem.eig.3.5.3
  (or later); inherited deferral from v3.5 / v3.5.1.
- **fem-eig-006-specific PML tuning** — if S3 retires
  fem-eig-003 but not fem-eig-006, fem-eig-006 tuning is queued
  for Phase 4.fem.eig.3.5.3.
- **Dispersive interior cavity fills under PML** — Phase
  4.fem.eig.3.6.
- **FEM-BEM hybrid** — Phase 4.fem.eig.4.
- **GPU sparse LU** — open-ended.

## Escape hatches

Per CLAUDE.md §5: any step blocking > 15 minutes → surface and stop.

Step-specific escape hatches:

- **S1 (`alpha_alpha(d)` closure):** if the `pow_alpha` evaluator
  introduces NaN/Inf in the per-tet assembled stiffness at any
  `alpha_grading_order ≥ 1` configuration in the
  `pml_open_boundary_assembly` unit tests, fall back to
  `alpha_grading_order = 0` for that configuration with a
  `// Phase 4.fem.eig.3.5.2: alpha_alpha(d) NaN guard — see §7 (a)`
  comment. The S1 unit tests still pass (the fallback path is
  bit-for-bit v3.5.1).
- **S2 (sweep binary):** if appending the H4 grid stage breaks
  the v3.5.1 CSV header / column count assertions in the
  example's docstring, regenerate the docstring CSV column list
  to include `alpha_grading_order`; downstream consumers (none
  in the workspace at v3.5.1) handle this automatically.
- **S3 (sweep wall-time):** if the full 18-row sweep exceeds 4 h
  wall-time (worst case `thickness = 16` extended mesh ~195 k
  tets), truncate to `thickness ∈ {12, 14}` only (12 rows;
  ~30 min `--release`) and ship the best-of-12 winner; document
  the `thickness = 16` sub-truncation in the S3 comment block
  and queue full H4 sweep for v3.5.3.
- **S3 (decision-tree exhaustion):** see step body. Ship
  `alpha_grading_order` API addition only; leave strict gates
  `#[ignore]`'d; queue Phase 4.fem.eig.3.5.3.
- **S4 (gate still fails after un-ignore):** treat as a base-SHA
  drift finding (S3 chose a row that does not actually retire
  under S4's CI-default invocation; possible cause:
  release-vs-debug numerical drift). Re-`#[ignore]` the gate and
  surface for re-sweep at the current SHA.

## References

- Companion spec
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading-design.md`
  — `alpha_alpha(d)` polynomial grading, 18-row H4 ablation
  grid, decision criteria.
- Companion ADR
  `docs/src/decisions/0045-phase-4-fem-eig-3-5-2-alpha-grading.md`.
- Parent spec
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune-design.md`
  — v3.5.1 grading retune (QQQQQQQQQ shipped).
- Parent plan
  `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune.md`
  — v3.5.1 R1-R5 ladder.
- Parent ADR
  `docs/src/decisions/0044-phase-4-fem-eig-3-5-1-grading-retune.md`
  — v3.5.1 ADR; §risks (b) `alpha_alpha(d)` deferral.
- Berenger 2002 *IEEE TAP* 50(3) (DOI 10.1109/8.999615) §VI —
  `alpha_alpha(d) = alpha_max · (1 - d/D)^n` canonical sweep.
- Roden-Gedney 2000 *IEEE MWCL* 10(5) — `sigma_max` / `kappa_max`
  Table-I defaults inherited via v3.5.1.
- `crates/yee-fem/src/open_boundary.rs:319-346` — QQQQQQQQQ R3
  comment block recording `[-58.13, -35.45] dB` H3 probe.
- CLAUDE.md §3, §4, §10.
