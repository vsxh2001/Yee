# Phase 4.fem.eig.3.5.1 — CFS-PML grading retune — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans`
> to drive this plan step-by-step.

**Companion spec:**
`docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune-design.md`
**Companion ADR:**
`docs/src/decisions/0044-phase-4-fem-eig-3-5-1-grading-retune.md`
**Base SHA:** `56386f1` (main HEAD; OOOOOOOOO P1-P7 CFS-PML shipped).
**Target phase:** 4.fem.eig.3.5.1 only. 4.fem.eig.3.5.2 `α_α(d)`
grading and 4.fem.eig.4 FEM-BEM hybrid explicitly deferred.
**Tech-stack additions:** none. Same `yee-fem` / `yee-validation`
surface area; a new `tools/cfs_pml_grading_sweep.rs` example binary
under `yee-validation`.

---

## Goal

Run the §4 ablation grid (32 configurations across H1/H2/H3) on
fem-eig-003 (WR-90 stub) + fem-eig-006 (100:10:1 high-aspect cavity),
pick a winning `(κ_max, m, thickness_cells)` triple plus the per-axis
`h_α` heuristic, ship it as the new `PmlConfig::default()`, and
un-ignore the three strict gates the OOOOOOOOO P5 measurement left in
`#[ignore]` purgatory.

Five-step ladder R1-R5 lands in a single merge train:

1. **R1** — extend `PmlConfig::resolved` to per-axis `h_α`; add
   `PmlMeshMeta` carrier.
2. **R2** — author `tools/cfs_pml_grading_sweep.rs` running the
   32-configuration ablation grid; emit CSV.
3. **R3** — analyse sweep CSV; pick winning defaults; update
   `PmlConfig::default`.
4. **R4** — un-ignore the three strict gates; verify both fixtures
   pass on the new defaults.
5. **R5** — ROADMAP refresh + tutorial note in
   `docs/src/tutorials/07-fem-open-cavity.md`.

CPU-only, single-threaded, scalar FP64 complex. No GPU. No new
dependencies. Same execution model as v3.5.

## Pre-flight

Before Step R1, confirm at base SHA `56386f1`:

1. `crates/yee-fem/src/open_boundary.rs` exposes
   `AbcOrder::CfsPml(PmlConfig)` (OOOOOOOOO P1) and
   `PmlConfig::resolved(freq_hz, h_cell)` with single-`h_cell`
   signature.
2. `crates/yee-fem/src/pml_mesh.rs` exposes
   `extend_mesh_with_pml(...) -> (TetMesh3D, Vec<PmlClass>, FaceIndexMap)`
   (OOOOOOOOO P2) and `PmlClass::{Interior, PmlX, PmlY, PmlZ, PmlXY,
   PmlYZ, PmlZX, PmlXYZ}`.
3. `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` carries
   exactly two `#[ignore]`'d strict gates with OOOOOOOOO P5
   measurement docstrings recording the `|S_{11}| ∈ [0.281, 0.423]`
   baseline.
4. `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
   carries the `#[ignore]`'d `fem_eig_006_magnitude_bounded` gate
   with OOOOOOOOO P5 measurement docstring `|S_{11}(30 GHz)| = 0.926`.

If (1)-(4) blocks, escape-hatch per CLAUDE.md §5 >15-min rule and
surface as a base-SHA drift finding; do **not** weaken any gate.

## File structure

| File | Action | Step | Responsibility |
|------|--------|------|----------------|
| `crates/yee-fem/src/open_boundary.rs` | Modify | R1, R3 | Per-axis `h_α`; `PmlMeshMeta`; new `PmlConfig::default` after R3 analysis. |
| `crates/yee-fem/src/lib.rs` | Modify | R1 | Re-export `PmlMeshMeta`. |
| `crates/yee-fem/tests/pml_open_boundary_assembly.rs` | Modify | R1 | Per-axis no-op equivalence assertion. |
| `tools/cfs_pml_grading_sweep.rs` | Create | R2 | Ablation sweep binary; emits CSV to stdout. |
| `crates/yee-validation/Cargo.toml` | Modify | R2 | Add `[[example]]` for sweep tool. |
| `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` | Modify | R4 | Remove both `#[ignore]` attributes. |
| `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs` | Modify | R4 | Remove the `#[ignore]` attribute. |
| `docs/src/tutorials/07-fem-open-cavity.md` | Modify | R5 | Note the new defaults; per-fixture override pattern. |
| `ROADMAP.md` | Modify | R5 | Phase 4.fem.eig.3.5.1 entry from planned to shipped. |

## Step R1 — per-axis `h_α` + `PmlMeshMeta`

**Lane:** `crates/yee-fem/src/open_boundary.rs`,
`crates/yee-fem/src/lib.rs`,
`crates/yee-fem/tests/pml_open_boundary_assembly.rs`.

Introduce

```rust
/// Per-axis mesh metadata used by the CFS-PML grading-parameter
/// resolver (Phase 4.fem.eig.3.5.1; replaces single-`h_cell` heuristic).
#[derive(Clone, Copy, Debug)]
pub struct PmlMeshMeta {
    /// Axis-aligned bounding-box extents (m), one per axis [x, y, z].
    pub extents: [f64; 3],
    /// Per-axis cell count of the original cavity mesh.
    pub cell_counts: [usize; 3],
}

impl PmlMeshMeta {
    /// Per-axis characteristic cell length h_α = extents[α] /
    /// cell_counts[α].
    pub fn h_per_axis(&self) -> [f64; 3] { /* ... */ }
}
```

`PmlConfig::resolved` changes signature to
`fn resolved(self, freq_hz: f64, mesh_meta: &PmlMeshMeta) -> Self`.
Internally it materialises a private `ResolvedPmlConfig` carrier that
holds *three* (`σ_α_max`, `α_α_max`) pairs — one per axis — plus the
shared `κ_max`, `m`, `thickness_cells`. The single-`h_cell` overload is
removed.

`pml_stretching_lambda` is reparametrised: its `s_for(d_α)` closure
takes the α-axis `(σ_α_max, κ_max, α_α_max)` triple instead of
referencing `cfg.sigma_max` directly. `Λ(ω) = diag(s_y s_z / s_x,
s_z s_x / s_y, s_x s_y / s_z)` is unchanged in shape; only its
per-axis ingredients change.

`OpenBoundarySolver::with_cfs_pml(config)` derives `PmlMeshMeta` from
the input mesh at builder time:

```rust
let bbox = mesh.bounding_box();
let extents = [bbox.x.length(), bbox.y.length(), bbox.z.length()];
let cell_counts = mesh.kuhn6_cell_counts(); // existing helper if present
let mesh_meta = PmlMeshMeta { extents, cell_counts };
let resolved = config.resolved(freq_hz_center, &mesh_meta);
```

If `kuhn6_cell_counts` does not exist on `TetMesh3D`, infer per-axis
cell counts by binning unique tet-vertex coordinates per axis (an
`itertools::dedup`-style pass on the sorted unique x/y/z coordinates).
Document the inference in a `// Phase 4.fem.eig.3.5.1: ...` comment.

**Pattern file:** mirror the OOOOOOOOO `PmlConfig::resolved`
implementation directly above its new replacement — the existing
single-`h_cell` body documents the analytic formulae; the new body
just dispatches per-axis.

**Test update in `pml_open_boundary_assembly.rs`:**

- `per_axis_resolver_zero_thickness_matches_scalar_path` — with
  `thickness_cells = 0` the per-axis resolver and the v3.5
  single-`h_cell` resolver produce identical `Λ(ω) = I` and identical
  per-tet assembled blocks (Frobenius difference < 1e-12).
- `per_axis_resolver_isotropic_mesh_matches_legacy` — for a
  `(4, 4, 4)` mesh on a 4×4×4 mm cube (`h_x = h_y = h_z`), the
  per-axis path produces `Λ(ω)` entries identical to the legacy
  single-`h_cell` path (the per-axis collapses to single by
  construction).

**DoD R1.**
- `cargo check -p yee-fem` exits 0.
- `cargo test -p yee-fem --test pml_open_boundary_assembly` exits 0.
- `grep -q 'struct PmlMeshMeta' crates/yee-fem/src/open_boundary.rs`
  exit 0.

## Step R2 — `tools/cfs_pml_grading_sweep.rs` ablation binary

**Lane:** `tools/cfs_pml_grading_sweep.rs`,
`crates/yee-validation/Cargo.toml`.

Author a `yee-validation` example binary running the §4 ablation grid.
The binary loops over the 32 configurations, calls
`run_fem_eig_003_wr90_stub_abc_with_config(cfg)` (a new public driver
helper added inline in `crates/yee-validation/src/lib.rs` — pure
parameter pass-through to `with_cfs_pml`, no algorithm change), and
`run_fem_eig_006_high_aspect_pml_with_config(cfg)`. Emits CSV to
stdout:

```text
hypothesis,kappa_max,m,thickness_cells,h_per_axis,
  fem_eig_003_s11_min_db,fem_eig_003_s11_max_db,
  fem_eig_006_s11_mag,fem_eig_003_runtime_s,fem_eig_006_runtime_s
H1,5.0,3,6,per_axis,...,...,...,...,...
H2,1.0,3,6,per_axis,...,...,...,...,...
H2,1.5,3,6,per_axis,...,...,...,...,...
... (32 rows)
```

Implements the §4 stopping rule: on each row, run fem-eig-003 first;
only run fem-eig-006 if fem-eig-003 worst-case `s11_max_db < -40`. On
the first row where both retire, emit a final `WINNER,...` row and
exit.

**Pattern file:** look at `examples/half-wave-dipole/src/main.rs` for
the yee-validation-driver call shape. The sweep binary is a thin loop
wrapping that pattern.

**Smoke test:** the sweep is **not** a CI gate — it is a one-off
exploration tool. The R2 DoD only requires that the binary compiles
and a single-configuration dry run (e.g.
`cargo run -p yee-validation --example cfs_pml_grading_sweep --release
-- --dry-run`) emits one CSV row.

**DoD R2.**
- `cargo build -p yee-validation --example cfs_pml_grading_sweep
  --release` exits 0.
- `cargo run -p yee-validation --example cfs_pml_grading_sweep
  --release -- --dry-run` exits 0 and writes one CSV row to stdout.
- `grep -q 'cfs_pml_grading_sweep' crates/yee-validation/Cargo.toml`
  exit 0.

## Step R3 — analyse CSV; pick defaults; update `PmlConfig::default`

**Lane:** `crates/yee-fem/src/open_boundary.rs`.

Run the full sweep (worst-case 75 min `--release`, see spec §4 stopping
rule). Capture stdout to
`docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-sweep.csv`
(not committed to git; referenced in R5 tutorial). Apply the §3
decision tree:

1. If the `H1`-row alone shows fem-eig-003 max-dB < −25 dB and
   fem-eig-006 magnitude < 0.5 → ship H1 standalone (defaults unchanged
   except for the per-axis resolver landed in R1).
2. Otherwise walk H2 rows; pick the smallest `κ_max ∈ {1.5, 2, 3}` that
   retires fem-eig-003 worst-case ≤ −40 dB.
3. Otherwise walk H3 rows; pick the `(m, thickness)` pair that retires
   *both* fixtures with the smallest `m × thickness` product.

Update `PmlConfig::default()` in `open_boundary.rs` with the winning
triple. Annotate the choice in a `// Phase 4.fem.eig.3.5.1 retune
(2026-05-20, sweep CSV row N):` comment block immediately above the
`impl Default for PmlConfig` block; record the winning row's
`(s11_min_db, s11_max_db, s11_006_mag)` measurements.

**Decision-tree exhaustion (escape-hatch):** if no row retires both
fixtures, leave `PmlConfig::default()` at the OOOOOOOOO values
`(κ_max=5, m=3, thickness=6)` but ship the per-axis `h_α` resolver
from R1, then jump to R4-alternative: do **not** un-ignore the three
strict gates; instead, update the three measurement docstrings with
the new H1-plus-best-H2-or-H3 baseline and leave `#[ignore]`'s in
place. Queue Phase 4.fem.eig.3.5.2 for `α_α(d)` grading per spec §7
(b).

**DoD R3.**
- The CSV exists locally (not committed) and shows the 32 ablation
  rows.
- `grep -q 'Phase 4.fem.eig.3.5.1 retune'
  crates/yee-fem/src/open_boundary.rs` exit 0.
- `cargo check -p yee-fem` exits 0 with the new defaults.

## Step R4 — un-ignore strict gates + verify pass

**Lane:** `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`,
`crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`.

If R3 ended with the §3-tree retire (the common case):

1. Remove `#[ignore = "..."]` from
   `fem_eig_003_strict_absorption_floor_gate` and
   `fem_eig_003_strict_passive_bound_continuum_limit`.
2. Remove `#[ignore = "..."]` from `fem_eig_006_magnitude_bounded`.
3. Run `cargo test -p yee-validation --release` and verify the three
   tests pass under the new `PmlConfig::default`.
4. Update the docstrings on each gate to record the post-retune
   measurement (`Phase 4.fem.eig.3.5.1 status: |S_{11}| ∈ [...] ⇒
   s11_db ∈ [...]; retires the strict gate`).

If R3 ended with decision-tree exhaustion (escape-hatch):

1. **Do not** remove any `#[ignore]`.
2. Update each docstring to record the new H1-on baseline.
3. ROADMAP refresh (R5) marks Phase 4.fem.eig.3.5.1 as "shipped
   per-axis `h_α` resolver; strict gates remain `#[ignore]`'d, queued
   for Phase 4.fem.eig.3.5.2 `α_α(d)` grading".

**DoD R4** (common case).
- `cargo test -p yee-validation --release
  --test fem_eig_003_wr90_stub_abc` exits 0 with all four tests
  passing (smoke + two un-ignored strict + the existing passive
  bound).
- `cargo test -p yee-validation --release
  --test fem_eig_006_high_aspect_pml` exits 0 with all three tests
  passing (magnitude, no-NaN, well-conditioned).
- `grep -c '#\[ignore' crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  prints `0`.
- `grep -c '#\[ignore' crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
  prints `0`.

## Step R5 — tutorial + ROADMAP refresh

**Lane:** `docs/src/tutorials/07-fem-open-cavity.md`, `ROADMAP.md`.

Add a "Grading parameter defaults" subsection to
`docs/src/tutorials/07-fem-open-cavity.md` showing:

- The new `PmlConfig::default()` triple from R3.
- A worked example overriding the default via the `pml_config` Python
  kwarg for a hypothetical user with a different aspect-ratio cavity.
- A short table of "knob → effect on `|S_{11}|`" derived from the R3
  sweep CSV (3-4 rows; pick representative cells).

Update `ROADMAP.md` Phase 4.fem.eig.3.5.1 entry from "planned /
optional" to "shipped"; link the un-ignored gate references in
`fem_eig_003_wr90_stub_abc.rs` and `fem_eig_006_high_aspect_pml.rs`.

**DoD R5.**
- `mdbook build docs/` exits 0.
- `grep -q '4.fem.eig.3.5.1' ROADMAP.md` exit 0.
- `grep -q 'Grading parameter defaults'
  docs/src/tutorials/07-fem-open-cavity.md` exit 0.

## Verification roll-up

After R5:

```bash
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
mdbook build docs/
```

All four must exit 0. The `--release` test invocation is required for
the fem-eig-003 strict gate (the budget is ~140 s `--release` per
sweep frequency point — equivalently ~ 480 s for the full 50-point
sweep; debug builds will time out).

## Out of scope

Explicitly deferred:

- **`α_α(d)` polynomial grading** — Phase 4.fem.eig.3.5.2. Per spec
  §7 (b); ablation only fires if R1+R2+R3 retune misses by > 5 dB on
  either fixture.
- **Rotated / non-Cartesian-aligned PML** — Phase 4.fem.eig.3.5.3 (or
  later); inherited deferral from v3.5.
- **Dispersive interior cavity fills under PML** — Phase
  4.fem.eig.3.6.
- **FEM-BEM hybrid** — Phase 4.fem.eig.4.
- **GPU sparse LU** — open-ended.

## Escape hatches

Per CLAUDE.md §5: any step blocking > 15 minutes → surface and stop.

Step-specific escape hatches:

- **R1 (per-axis resolver):** if `TetMesh3D::kuhn6_cell_counts` (or
  equivalent helper) does not exist and inferring per-axis cell counts
  from sorted unique vertex coordinates blocks > 15 min, fall back to
  re-using the v3.5 single-`h_cell` heuristic for any axis where the
  inference returns `None`, with a `// fallback: per-axis cell count
  inference failed for axis α` comment. The R1 unit tests still pass
  (the fallback path is bit-for-bit the v3.5 path).
- **R2 (sweep wall-time):** if the full 32-row sweep exceeds 2 h
  wall-time, truncate to H1 + H2 only (7 rows, ~20 min `--release`)
  and ship the best H1+H2 winner. Document the H3 sub-truncation in
  the R3 comment block and queue full H3 sweep for Phase
  4.fem.eig.3.5.2.
- **R3 (decision-tree exhaustion):** see step body. Ship per-axis
  `h_α` only; leave strict gates `#[ignore]`'d; queue Phase
  4.fem.eig.3.5.2.
- **R4 (gate still fails after un-ignore):** treat as a base-SHA drift
  finding (R3 chose a row that does not actually retire under R4's
  CI-default invocation; possible cause: release-vs-debug numerical
  drift). Re-`#[ignore]` the gate and surface for re-sweep at the
  current SHA.

## References

- Companion spec
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune-design.md`
  — H1/H2/H3 hypothesis tree, ablation grid, decision criteria.
- Companion ADR
  `docs/src/decisions/0044-phase-4-fem-eig-3-5-1-grading-retune.md`.
- Parent spec
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
  — CFS-PML formulation v3.5.
- Parent plan
  `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-cfs-pml.md`
  — P1-P7 OOOOOOOOO landed.
- Parent ADR `docs/src/decisions/0043-phase-4-fem-eig-3-5-cfs-pml-scope.md`
  §risks "PML grading parameter sensitivity" — the deferral this plan
  fulfils.
- Berenger 2002 *IEEE TAP* 50(3) (DOI 10.1109/8.999615) — CFS-PML
  parameter sweep reference; figures 4-7 the empirical basin source.
- Roden-Gedney 2000 *IEEE MWCL* 10(5) — Table-I defaults the
  OOOOOOOOO baseline inherits.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  §"OOOOOOOOO P5 status" — `|S_{11}| ∈ [0.281, 0.423]` measurement.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
  §"OOOOOOOOO P5 status" — `|S_{11}(30 GHz)| = 0.926` measurement.
