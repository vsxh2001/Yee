# Phase 4.fem.eig.3.5 — CFS-PML — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans`
> to drive this plan step-by-step.

**Companion spec:** `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
**Companion ADR:** `docs/src/decisions/0043-phase-4-fem-eig-3-5-cfs-pml-scope.md`
**Base SHA:** `fc62f82` (main HEAD; NNNNNNNNN mesh refinement merged on
top of Phase 4.fem.eig.3 F1-F7 + LLLLLLLLL dispersive sigma_factor fix).
**Target phase:** 4.fem.eig.3.5 only. 4.fem.eig.3.5.1 grading-parameter
ablation, 4.fem.eig.4 FEM-BEM hybrid explicitly deferred.
**Tech-stack additions:** none new. Same
`faer::sparse::FaerLuSolver<Complex64>` surface.

---

## Goal

Replace the Phase 4.fem.eig.3 Engquist-Majda 2nd-order ABC with a
CFS-PML (Roden-Gedney 2000) volumetric buffer-layer absorber on the
open-boundary FEM solver, retiring both `#[ignore]`'d fem-eig-003
strict gates and adding a new fem-eig-006 high-aspect-ratio stress
fixture.

Seven-step ladder P1-P7 lands in one merge train:

1. **P1** — `PmlConfig` + `PmlRegion` types + `AbcOrder::CfsPml`
   variant. Wire-only; assembly still no-ops.
2. **P2** — PML mesh extension (`extend_mesh_with_pml`) building the
   tet shell outside the cavity.
3. **P3** — anisotropic per-tet `ε_tensor(ω)` assembly path
   (`assemble_tet_element_complex_anisotropic`).
4. **P4** — `with_cfs_pml` builder + OpenBoundarySolver wire-in
   (face classification, per-tet PML class lookup, per-frequency
   stretched-coordinate `Λ(ω)` evaluation).
5. **P5** — un-ignore fem-eig-003 strict gates + new fem-eig-006
   (high-aspect 100:10:1 cavity at 30 GHz).
6. **P6** — Python binding `pml_config` kwarg on
   `yee.fem.solve_open_cavity`.
7. **P7** — tutorial update (`docs/src/tutorials/07-fem-open-cavity.md`
   demonstrates PML mode) + ROADMAP refresh.

CPU-only, single-threaded, scalar FP64 complex, no GPU, single dominant
mode per port, Cartesian-aligned PML only. Same execution model as v3.

## Pre-flight

Before Step P2 starts, confirm at base SHA `fc62f82`:

1. `crates/yee-fem/src/open_boundary.rs` exposes `AbcOrder` with
   variants `First` and `Second` (Phase 4.fem.eig.3 F3+F4). Verify
   the existing 2nd-order ABC integration test
   `crates/yee-fem/tests/abc2_face_block.rs` is green; P1-P4 must not
   regress it.
2. `crates/yee-fem/src/element.rs` exposes
   `assemble_tet_element_complex` taking scalar `ε`, `μ`. P3 adds the
   anisotropic-tensor sibling alongside; the scalar entry point stays.
3. `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` carries
   the `#[ignore]`'d `fem_eig_003_strict_absorption_floor_gate` and
   `fem_eig_003_strict_passive_bound_continuum_limit`. P5 removes both
   `#[ignore]`s; confirm both are present at base SHA.
4. `crates/yee-validation/src/lib.rs` defines
   `run_fem_eig_003_wr90_stub_abc` returning a `FemEig003Result`. P5
   adds `run_fem_eig_006_high_aspect_pml` next to it.

If (1)-(4) blocks, escape-hatch per CLAUDE.md §5 >15-min rule and
surface as a base-SHA drift finding; do **not** weaken the strict
gates.

## File structure

| File | Action | Step | Responsibility |
|------|--------|------|----------------|
| `crates/yee-fem/src/open_boundary.rs` | Modify | P1, P4 | `PmlConfig`, `PmlRegion`, `AbcOrder::CfsPml`, `with_cfs_pml`. |
| `crates/yee-fem/src/pml_mesh.rs` | Create | P2 | `extend_mesh_with_pml`, `PmlClass`, `FaceIndexMap`. |
| `crates/yee-fem/src/element.rs` | Modify | P3 | `assemble_tet_element_complex_anisotropic`. |
| `crates/yee-fem/src/lib.rs` | Modify | P1, P2 | Re-export new types. |
| `crates/yee-fem/tests/pml_mesh_extension.rs` | Create | P2 | Unit test. |
| `crates/yee-fem/tests/anisotropic_tet_assembly.rs` | Create | P3 | Scalar-equivalence + diagonal-Λ unit tests. |
| `crates/yee-fem/tests/pml_open_boundary_assembly.rs` | Create | P4 | PML-end-to-end smoke. |
| `crates/yee-validation/src/lib.rs` | Modify | P5 | Add fem-eig-006 driver; flip fem-eig-003 to CFS-PML. |
| `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` | Modify | P5 | Remove `#[ignore]` on both strict gates. |
| `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs` | Create | P5 | Strict gate driver. |
| `crates/yee-py/src/fem.rs` | Modify | P6 | `pml_config` kwarg. |
| `docs/src/tutorials/07-fem-open-cavity.md` | Modify | P7 | PML demo section. |
| `ROADMAP.md` | Modify | P7 | Phase 4.fem.eig.3.5 entry. |

## Step P1 — `PmlConfig` + `PmlRegion` types + `AbcOrder::CfsPml`

**Lane:** `crates/yee-fem/src/open_boundary.rs`, `crates/yee-fem/src/lib.rs`.

Add `PmlConfig`, `PmlRegion`, and extend `AbcOrder` with the
`CfsPml(PmlConfig)` variant per spec §5. Step P1 is wire-only: the new
variant is recognised by `match` in the assembly path but a placeholder
arm returns `Error::NotEnabled("CFS-PML assembly not wired until P4")`.

`PmlConfig::default()` uses sentinel zeros for `sigma_max` and
`alpha_max`; the solver's `with_cfs_pml(cfg)` builder (P4) recomputes
both from the frequency band and the mesh characteristic length using
Roden-Gedney 2000 §III formulae. Sentinel handling lives in
`PmlConfig::resolved(freq_hz, h_cell)` returning a fully populated
copy.

`Default for AbcOrder` stays `Second` (Phase 4.fem.eig.3 default). Adding
`CfsPml` is additive; existing match arms exhaustive with explicit
`AbcOrder::CfsPml(_)` placeholder branches.

**DoD P1.**
- `cargo check -p yee-fem` exits 0.
- `cargo test -p yee-fem` exits 0 (no new tests yet; existing tests
  untouched).
- `grep -q 'CfsPml' crates/yee-fem/src/open_boundary.rs` exit 0.

## Step P2 — PML mesh extension

**Lane:** `crates/yee-fem/src/pml_mesh.rs`, `crates/yee-fem/src/lib.rs`,
`crates/yee-fem/tests/pml_mesh_extension.rs`.

New module `pml_mesh` implements `extend_mesh_with_pml` per spec §4.2:

1. Identify the original mesh's axis-aligned bounding box and the
   normal direction of every `pml_faces` entry.
2. For each `pml_faces` entry, replicate the boundary face vertex
   layer outward by `thickness_cells` Kuhn-6 brick layers; for edges
   and corners shared between two/three PML axes, generate the wedge
   tets in the standard Berenger 1994 corner-wedge pattern.
3. Tag each new tet with its `PmlClass` and compute the per-axis
   depths `d_x, d_y, d_z` from the original cavity boundary (in
   metres, not cells; the `Λ(ω)` evaluator (P4) uses `d / (thickness
   · h_cell)` to derive the polynomial-grading parameter).
4. The outer-most face of the extended mesh (the truncation surface)
   is tagged `FaceKind::Pec` by default.

For non-Cartesian-aligned PML faces, return
`Error::InvalidArgument("CFS-PML supports Cartesian-aligned faces only;
got face normal {n_vec}; consider rotated-PML extension Phase
4.fem.eig.3.5.1")`.

**Pattern file:** mirror the boundary-vertex extrusion logic in
`yee_mesh::TetMesh3D::from_brick_kuhn` (the existing Kuhn-6 brick
constructor); the PML mesh is a layered re-application of that.

**Unit tests in `tests/pml_mesh_extension.rs`:**

- `pml_shell_tet_count_matches_formula` — for a `(4, 4, 4)` Kuhn-6
  brick cavity with `thickness_cells = 6` on every face, assert the
  extended-mesh tet count matches the analytic formula `n_interior +
  6 · (face_area_cells · 6) + 12 · (edge_length_cells · 6²) + 8 · 6³`
  (interior + face shells + edge wedges + corner wedges; Kuhn-6
  multipliers applied uniformly).
- `pml_inner_boundary_has_continuous_vertex_layer` — every vertex on
  the original cavity boundary is preserved at the inner PML face;
  the PML mesh shares vertices with the cavity (no duplicates).
- `pml_outer_face_tagged_pec` — every face on the outermost PML
  surface is `FaceKind::Pec`.
- `pml_class_depth_monotonic_outward` — for tets along the +x PML
  axis, `d_x` is monotonically increasing with depth.

**DoD P2.**
- `cargo test -p yee-fem --test pml_mesh_extension` exits 0.
- `grep -q 'pub fn extend_mesh_with_pml' crates/yee-fem/src/pml_mesh.rs`
  exit 0.

## Step P3 — anisotropic per-tet ε tensor assembly

**Lane:** `crates/yee-fem/src/element.rs`,
`crates/yee-fem/tests/anisotropic_tet_assembly.rs`.

Implement `assemble_tet_element_complex_anisotropic` per spec §4.3. The
new helper takes `eps_tensor: SMatrix<Complex64, 3, 3>` and
`mu_tensor_inv: SMatrix<Complex64, 3, 3>` instead of the existing
scalar `eps, mu`. The integrand becomes

```text
    K_{ij}  =  ∫_T  ( ∇×N_i )^T · μ_inv · ( ∇×N_j )  dV
            =  V_tet · ( ∇×N_i )^T · μ_inv · ( ∇×N_j )       (curl constant in tet)
    M_{ij}  =  ∫_T  N_i^T · ε · N_j  dV
            =  V_tet · Σ_{a,b} c_{ab}^{ij} · ( unit_a^T · ε · unit_b ),
```

where `c_{ab}^{ij}` are the standard barycentric mass-integration
constants (the integrand `λ_p λ_q` integrates to `V/20 · (1 + δ_{pq})`
on a tet) and `unit_a, unit_b` are the canonical Cartesian unit vectors
selecting the relevant `ε` entry.

For diagonal `ε_tensor = diag(ε_x, ε_y, ε_z)` the inner sum collapses
to three scalar mass contributions; this is the only case that arises
in v3.5 (Cartesian-aligned PML).

**Unit tests in `tests/anisotropic_tet_assembly.rs`:**

- `scalar_equivalence_when_tensor_is_scalar_times_identity` — for
  `ε_tensor = ε · I, μ_inv_tensor = (1/μ) · I` the output matches
  `assemble_tet_element_complex` bit-for-bit (Frobenius difference <
  1e-12).
- `diagonal_anisotropic_block_is_complex_symmetric` — for a diagonal
  complex `ε_tensor = diag(c_x, c_y, c_z)` the resulting `M` is
  complex-symmetric (`M = M^T`, no conjugation).
- `off_diagonal_tensor_rejected_until_v3_5_1` — `ε_tensor` with a
  non-zero off-diagonal entry produces `Error::NotEnabled` (rotated
  PML deferred).

**DoD P3.**
- `cargo test -p yee-fem --test anisotropic_tet_assembly` exits 0.
- `grep -q 'pub fn assemble_tet_element_complex_anisotropic' crates/yee-fem/src/element.rs`
  exit 0.

## Step P4 — `with_cfs_pml` builder + OpenBoundarySolver wire-in

**Lane:** `crates/yee-fem/src/open_boundary.rs`,
`crates/yee-fem/tests/pml_open_boundary_assembly.rs`.

Implement `with_cfs_pml(self, config: PmlConfig) -> Self` per spec §5.
The builder triggers, at construction time:

1. Call `extend_mesh_with_pml` (P2) on the input mesh + ABC faces.
2. Cache the `PmlClass` per extended-mesh tet and the
   `FaceIndexMap` for later port/ABC re-tagging.
3. Resolve the sentinel `PmlConfig` parameters
   (`sigma_max ≈ (m+1) / (150 π h_cell sqrt(ε_r))`, `alpha_max ≈ ω₀ ε_0`)
   from the band centre and mean tet edge length.

At every per-frequency assembly:

1. For each interior tet (`PmlClass::Interior`), call the scalar
   `assemble_tet_element_complex` path bit-for-bit unchanged.
2. For each PML tet, compute `Λ(d_x, d_y, d_z; ω) = diag(s_y s_z /
   s_x, s_z s_x / s_y, s_x s_y / s_z)` with each `s_α(ω) = κ_α(d_α) +
   σ_α(d_α) / (α_α + j ω ε_0)`. Then call
   `assemble_tet_element_complex_anisotropic(ε_iso · Λ, μ_iso · Λ_inv,
   ω)`.
3. The boundary integral on the original cavity-PML interface is
   *zero* (Λ = I at d = 0 by polynomial grading; no surface
   contribution). The outermost PEC face Dirichlet-eliminates as
   usual.

**Unit test in `tests/pml_open_boundary_assembly.rs`:**

- `pml_assembly_matches_scalar_on_zero_thickness` — with
  `thickness_cells = 0`, the PML mesh extension is a no-op and the
  assembled matrix matches the v3 2nd-order-ABC path bit-for-bit.
- `pml_assembly_finite_at_dc` — the assembled matrix at `ω = 0.01 ω_c`
  (well below cutoff) has bounded entries (`max |entry| < 1e6`); this
  is the CFS `α_α > 0` causality canary, distinguishing CFS-PML from
  the original Berenger 1994 PML which would diverge.

**DoD P4.**
- `cargo test -p yee-fem --test pml_open_boundary_assembly` exits 0.
- Existing `cargo test -p yee-fem` suite continues to pass (no
  regression on v3 2nd-order ABC).
- `grep -q 'pub fn with_cfs_pml' crates/yee-fem/src/open_boundary.rs`
  exit 0.

## Step P5 — fem-eig-003 strict un-ignore + new fem-eig-006

**Lane:** `crates/yee-validation/src/lib.rs`,
`crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`,
`crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`.

In `crates/yee-validation/src/lib.rs`, modify
`run_fem_eig_003_wr90_stub_abc` to call
`with_cfs_pml(PmlConfig::default())` instead of the v3
`with_abc_order(AbcOrder::Second)` path. Verify the sweep
runtime (informational) stays under 480 s in `--release` on the
existing `(24, 12, 36)` cavity mesh.

In `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`, remove
both `#[ignore]` attributes per spec §6, and update the absorption
band assertion bounds to `[-60, -40] dB`. Keep the passive-bound,
smoothness, and TE_{101} no-resonance assertions unchanged.

In `crates/yee-validation/src/lib.rs`, add
`run_fem_eig_006_high_aspect_pml` per spec §6: 100 mm × 10 mm × 1 mm
high-aspect cavity, TE-mode drive at `x = 0`, CFS-PML at `x = 100 mm`,
30 GHz single-frequency solve, returns
`FemEig006Result { s11: Complex64, lu_condition_estimate: f64,
notes: String }`.

In `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`,
strict gates per spec §6:

- `fem_eig_006_magnitude_bounded` — `|S_{11}(30 GHz)| < 0.1`.
- `fem_eig_006_no_nan_inf` — `S_{11}` is finite.
- `fem_eig_006_lu_well_conditioned` — `lu_condition_estimate < 1e10`.

**DoD P5.**
- `cargo test -p yee-validation --test fem_eig_003_wr90_stub_abc`
  exits 0 with no `#[ignore]`'d tests passing.
- `cargo test -p yee-validation --test fem_eig_006_high_aspect_pml`
  exits 0.
- `cargo test -p yee-validation` exits 0 (no regression on Phase
  4.fem.eig.{0,1,2,3} drivers).

## Step P6 — Python binding `pml_config` kwarg

**Lane:** `crates/yee-py/src/fem.rs`, `crates/yee-py/tests/test_fem.py`
(extend).

Add a `pml_config: Option<&PyDict>` kwarg to
`yee.fem.solve_open_cavity` per spec §5. When `Some`, parse keys
(`thickness_cells`, `sigma_max`, `alpha_max`, `kappa_max`, `m`) into a
`PmlConfig` and call `with_cfs_pml`; otherwise fall back to the v3
`abc_order` kwarg (default `"second"`).

Extend `crates/yee-py/tests/test_fem.py` with a
`test_solve_open_cavity_with_pml` case driving the fem-eig-003 mini
fixture (4×2×6 mesh, single freq) and asserting `|S_11|` is finite and
sub-unity.

**DoD P6.**
- `maturin develop --release` succeeds (no new dependency).
- `pytest crates/yee-py/tests/test_fem.py -v` exits 0.

## Step P7 — tutorial + ROADMAP refresh

**Lane:** `docs/src/tutorials/07-fem-open-cavity.md`, `ROADMAP.md`.

Add a new "CFS-PML mode" section to
`docs/src/tutorials/07-fem-open-cavity.md` walking through the
`pml_config` kwarg on the fem-eig-003 stub fixture, showing the
absorption-floor improvement vs the 2nd-order ABC default.

Update `ROADMAP.md` Phase 4.fem.eig.3.5 entry from "planned" to
"shipped"; link the fem-eig-003-strict and fem-eig-006 gate
references.

**DoD P7.**
- `mdbook build docs/` exits 0.
- `grep -q '4.fem.eig.3.5' ROADMAP.md` exit 0.

## Verification roll-up

After P7:

```bash
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
mdbook build docs/
```

All four must exit 0. The `--release` test invocation is required for
the fem-eig-003 strict gate (the budget is 480 s in `--release`; debug
builds will time out).

## Out of scope

Explicitly deferred to follow-up phases:

- **Rotated / non-Cartesian-aligned PML** — Phase 4.fem.eig.3.5.1.
  v3.5 rejects with `Error::InvalidArgument`.
- **Grading-parameter ablation sweep** (`thickness_cells`, `m`,
  `α_max`) — Phase 4.fem.eig.3.5.1.
- **PML on dispersive interior cavity fills** — Phase 4.fem.eig.3.6
  (combines v1's `ε(ω)` Newton tracker with v3.5 PML).
- **FEM-BEM hybrid** — Phase 4.fem.eig.4.
- **GPU sparse LU** — open-ended.

## Escape hatches

Per CLAUDE.md §5: any step blocking > 15 minutes → surface and stop.

Step-specific escape hatches:

- **P2 (mesh extension):** if Kuhn-6 corner-wedge tet generation
  blocks, fall back to a simpler "one face only" PML (no edge / corner
  wedges) for fem-eig-003 — the WR-90 stub has a single ABC face so
  edge / corner wedges are not strictly required. fem-eig-006 (3
  faces sharing an edge / corner) would then move to v3.5.1.
- **P3 (anisotropic assembly):** if the mass-matrix barycentric
  formula is uncertain for the off-diagonal anisotropic case, restrict
  the helper to diagonal `ε_tensor` only and gate behind a runtime
  check.
- **P5 (fem-eig-003 strict):** if the un-ignored gates still fail by
  > 5 dB above the `[-60, -40] dB` band after P1-P4 land, **do not
  weaken the bounds**. Surface the measurement, leave the `#[ignore]`s
  in place, and queue a Phase 4.fem.eig.3.5.1 grading-parameter
  retune.

## References

- Companion spec
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
  — math and rationale.
- Companion ADR
  `docs/src/decisions/0043-phase-4-fem-eig-3-5-cfs-pml-scope.md`.
- Berenger 1994; Kuzuoglu-Mittra 1996; Roden-Gedney 2000 — primary
  CFS-PML literature; full citations in companion spec §10.
- Jin §10.8 — PML for frequency-domain FEM.
- Phase 4.fem.eig.3 spec / plan / ADR-0042 — parent.
- NNNNNNNNN mesh refinement
  (`crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  §"NNNNNNNNN status") — the `[-2.22e-2, -2.86e-5] dB` baseline that
  motivates v3.5.
