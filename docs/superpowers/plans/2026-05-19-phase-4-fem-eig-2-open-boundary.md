# Phase 4.fem.eig.2 — Open-Boundary FEM (ABC + Wave Ports) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans`
> to drive this plan step-by-step.

**Companion spec:** `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
**Companion ADR:** `docs/src/decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md`
**Base SHA:** `60ed512` (CLAUDE.md §1 — Phase 4.fem.eig.0 D1–D7 walking skeleton + Phase 4.fem.eig.1 D1–D7 dispersive Newton tracker + fem-eig-002 production gate all shipped to `main`).
**Target phase:** 4.fem.eig.2 only. 4.fem.eig.2.0.1 / 2.1 / 2.5 / 3 are explicitly deferred — see §"Out of scope".
**Tech-stack additions:** none new. `faer` already exposes `Complex64` sparse LU (verified live by Phase 4.fem.eig.1 D2). No new direct dep.

---

## Goal

Phase 4.fem.eig.2 extends the shipped Phase 4.fem.eig.1 closed-cavity
dispersive eigensolver to **open-boundary driven analysis**. Two boundary
kinds are added beyond the Phase 4.fem.eig.0/1 default PEC:

1. **1st-order Engquist–Majda ABC** on tagged exterior faces. Adds a
   `+ j k₀ · area · (n̂ × N_i) · (n̂ × N_j)` block per face into the global
   complex stiffness matrix `K(ω)`.
2. **Modal wave-port faces** with a `NumericalCrossSection`-sourced TE_{10}
   (or other dominant-mode) profile and an incident-amplitude right-hand
   side. Extracts `S_{11}(f)` via a per-frequency complex sparse LU
   back-substitution.

The delivered pipeline: a WR-90 stub mesh (22.86 × 10.16 × 30 mm) is
consumed by a new `yee-fem::open_boundary` module which classifies each
exterior face as `Pec | Abc | WavePort(p)`, assembles a complex-symmetric
driven system per swept frequency, solves once via the Phase 4.fem.eig.1
`faer::sparse::FaerLuSolver<Complex64>` surface, and projects the FEM
solution against the cross-section's modal profile to extract the
S-parameter matrix. Validation gate `fem-eig-003` enforces `|S_{11}(f)|`
within ±0.5 dB of Pozar §3.3 closed-form across 50 sweep points in
8–12 GHz.

CPU-only, single-threaded, scalar FP64 complex, no GPU, single incident
mode per port, scalar isotropic real `ε_r` and `μ_r` on the driven sweep,
PEC + ABC + modal port boundary classification only — same execution model
as Phase 4.fem.eig.1, lifted one axis (closed cavity → open) and wrapped in
an outer driven-frequency sweep.

## Pre-flight — face classification + cross-section modal source

Spec §10 risks and §7 solver-reuse decisions name the load-bearing
dependencies. Before Step E2 starts, confirm at base SHA `60ed512`:

1. `crates/yee-mesh::TetMesh3D` already exposes face iteration with
   stable `FaceId` indices and an outward-normal accessor consistent with
   the Phase 4.fem.eig.0 boundary-edge classification. If face IDs are
   not yet first-class (Phase 4.fem.eig.0 may have only enumerated edges),
   surface this and add a face accessor as Step E0 before E1 / E2 start.
2. `yee_mom::eigensolver::NumericalCrossSection::e_tangential_at(x, y)`
   is a public API (Phase 1.3.1.1) returning a `Vector2<f64>` tangent
   field. The FEM port consumer calls this on per-face Gauss points;
   verify no API drift since 1.3.1.1 ship. If `e_tangential_at` is
   private or renamed, surface and fix at the cross-section side first
   (in the MoM lane, NOT in this lane).
3. `faer::sparse::FaerLuSolver<Complex64>` already factors and
   back-substitutes complex-symmetric (non-Hermitian) matrices in the
   workspace pin — Phase 4.fem.eig.1 D2 exercises this on the lossy-SiO₂
   gate. No new pre-flight needed.
4. The Phase 4.fem.eig.0 boundary classifier walks the tet mesh, finds
   exterior faces, and elimitates the corresponding edges as PEC. v2's
   face classifier extends this: instead of one tag (boundary or not), it
   carries three (`Pec | Abc | WavePort`). The default for an
   unannotated exterior face stays `Pec`, preserving every v0/v1 caller.

If (1)–(2) blocks, escape-hatch per the standard >15-min rule
(CLAUDE.md §5) and surface as a Phase 4.fem.eig.2.0.1 finding; do **not**
weaken the fem-eig-003 gate to compensate.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-fem/src/element.rs` | Modify | Add `assemble_abc_face_block` and `assemble_port_face_block` (3×3 complex face blocks). Extend `assemble_tet_element_complex` with an optional `abc_faces` arg; existing callers pass `&[]`. |
| `crates/yee-fem/src/assembly.rs` | Modify | Face-iteration over per-face `FaceKind`; complex-symmetric driven-system assembly; per-port modal RHS contribution. |
| `crates/yee-fem/src/solve.rs` | Modify | Driven-solve helper `solve_driven_at_frequency(K_driven, rhs)` — single complex sparse LU back-substitution, reusing the Phase 4.fem.eig.1 LU surface. |
| `crates/yee-fem/src/open_boundary.rs` | Create | `OpenBoundarySolver`, `WavePortFace`, `FaceKind`, `SParameterRow`, `SParameters`, `solve_at_frequency`, `sweep`. |
| `crates/yee-fem/src/lib.rs` | Modify | `pub mod open_boundary;`, re-export `OpenBoundarySolver`, `WavePortFace`, `FaceKind`, `SParameterRow`, `SParameters`. |
| `crates/yee-fem/tests/abc_face_block.rs` | Create | E1 unit test — `assemble_abc_face_block` on a hand-built triangular face matches the closed-form `+ j k₀ · area · ...` to `1e-12`. |
| `crates/yee-fem/tests/port_face_block.rs` | Create | E2 unit test — `assemble_port_face_block` recovers TE_{10} modal projection on a synthetic WR-90 cross-section. |
| `crates/yee-fem/tests/driven_solve_pec_stub.rs` | Create | E3 unit test — drive a TE_{10} wave into a fully-PEC-closed WR-90 stub; expect `|S_{11}| = 0 dB ± 0.05` (ABC off, full reflection). |
| `crates/yee-fem/tests/abc_face_eats_wave.rs` | Create | E4 unit test — TE_{10} into a long ABC-terminated stub; assert `|S_{11}|` is below `−35 dB` at mid-band (ABC absorbs the wave). |
| `crates/yee-validation/src/lib.rs` | Modify | E5 — `run_fem_eig_003_wr90_stub_with_abc` driver. |
| `crates/yee-validation/tests/fem_eig_003_wr90_stub_with_abc.rs` | Create | E5 production-gate test. |
| `crates/yee-fem/validation/README.md` | Modify | `fem-eig-003 (WR-90 stub + ABC)` row. |
| `crates/yee-py/src/fem.rs` | Modify (E6, optional) | Python binding `yee.fem.solve_open_cavity(...)`. |
| `crates/yee-py/tests/test_fem_open_boundary.py` | Create (E6, optional) | Python pytest re-running fem-eig-003 from Python. |

No changes to `yee-mom`, `yee-cuda`, `yee-gui`, `yee-plotters`, `yee-io`,
`yee-cli`, `yee-surrogate`, `yee-mesh` (except possibly Step E0 face-ID
surface, if the pre-flight finds it missing). The
`#![forbid(unsafe_code)]` floor is preserved across every touched crate.

## Step ladder

### Step E1 — ABC face-block element helper

- **Brief:** Add `assemble_abc_face_block(face_vertices, outward_normal,
  k0, mu_r_face) -> SMatrix<Complex64, 3, 3>` to
  `crates/yee-fem/src/element.rs`. Compute the three edge tangents
  `t_i = v_{i+1} − v_i` (modulo 3); form `n̂ × N_i = n̂ × t_i` (constant
  per Whitney-1 face basis); the 3×3 block entry is `+ j k₀ · area /
  mu_r_face · (n̂ × N_i) · (n̂ × N_j)`. The face area is `0.5 · ||t_0 ×
  t_1||`. Pattern file: `crates/yee-fem/src/element.rs` itself —
  preserve the existing reference-tet edge-ordering doc. Local-to-global
  orientation sign is applied at the *assembly* layer (E3); the element
  function emits the unsigned block.
- **Lane:** `crates/yee-fem/src/element.rs`,
  `crates/yee-fem/tests/abc_face_block.rs` (create).
- **Base SHA dep:** none — branches off `60ed512` directly.
- **DoD:** unit tests pass — `abc_face_block_matches_closed_form` (a
  unit-area equilateral triangle in the xy-plane with `n̂ = ẑ`, `k₀ =
  2π / 0.03`, `μ_r = 1`: each diagonal entry equals `+ j k₀ · area ·
  ||N_i × ẑ||² = + j k₀ · (1/3)` to `1e-12`);
  `abc_face_block_is_complex_symmetric` (`block^T == block`, not
  Hermitian).
- **Verification:** `cargo test -p yee-fem --release abc_face_block &&
  cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the cross-product sign convention
  → cross-check with the analytic TE_{10} reflection formula on a unit
  test stub before scaling to the gate.
- **LOC:** ~150.

### Step E2 — wave-port face-block element helper + modal RHS

- **Brief:** Two related additions to `crates/yee-fem/src/element.rs`:
  - `assemble_port_face_block(face_vertices, outward_normal, beta_mode,
    mu_r_face) -> SMatrix<Complex64, 3, 3>` — analog to E1 with `β_mode`
    replacing `k₀`. `β_mode = sqrt(k₀² ε_r μ_r − k_c²)` is computed by
    the caller from the cross-section eigensolver.
  - `assemble_port_face_rhs(face_vertices, outward_normal, beta_mode,
    a_inc, e_mode_at_gauss: [Vector2<f64>; 3]) -> SVector<Complex64, 3>`
    — per-face right-hand side encoding `+ 2 j β · a_inc · ∫_face N_i ·
    e_mode dS` computed with 3-point Gauss quadrature on the reference
    triangle. The caller pre-evaluates `e_mode` at the three Gauss
    points via `NumericalCrossSection::e_tangential_at`.
- **Lane:** `crates/yee-fem/src/element.rs`,
  `crates/yee-fem/tests/port_face_block.rs` (create).
- **Base SHA dep:** none — parallel-safe with E1. The two element
  helpers touch disjoint test files.
- **DoD:** unit tests pass — `port_face_block_matches_abc_at_beta_eq_k0`
  (when `β_mode = k₀`, the port face block equals the ABC face block
  E1 emits — sanity check); `port_rhs_te10_matches_analytic` (using the
  analytic TE_{10} profile `e_mode = ŷ sin(π x / a)`, the RHS column for
  a known face on a WR-90 cross-section matches the closed-form
  `∫ N_i · sin(π x / a) ŷ dS` to `1e-8`).
- **Verification:** `cargo test -p yee-fem --release port_face_block &&
  cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on
  `NumericalCrossSection::e_tangential_at` API drift since Phase 1.3.1.1
  → use the analytic TE_{10} profile `ŷ sin(π x / a)` for the test and
  surface a Phase 4.fem.eig.2.0.1 finding to revisit the MoM API in a
  follow-up PR. Do not silently rename the MoM API from this lane.
- **LOC:** ~220.

### Step E3 — `OpenBoundarySolver` + face-kind assembly

- **Brief:** Implement `OpenBoundarySolver::new` + `solve_at_frequency`
  in `crates/yee-fem/src/open_boundary.rs`. Algorithm:
  - At construction, classify every exterior face into `Pec | Abc |
    WavePort(p)` from the caller-supplied `face_kinds` array. Edges
    that lie on a `Pec`-tagged face are added to the PEC Dirichlet set
    (precedence over any `WavePort` face that shares an edge — spec §10).
  - Per swept frequency `ω`:
    - Assemble the closed-cavity complex matrices `K(ω)`, `M(ω)` from
      Phase 4.fem.eig.1's `assemble_complex` path.
    - For each `Abc`-tagged face: call `assemble_abc_face_block` (E1)
      and scatter into `K(ω)` at the corresponding global edges with
      orientation signs.
    - For each `WavePort(p)`-tagged face: call
      `NumericalCrossSection::solve(ω)` once per port per frequency to
      recover the dominant-mode `β_mode` and `e_mode`. Call
      `assemble_port_face_block` (E2) and scatter into `K(ω)`. Call
      `assemble_port_face_rhs` (E2) with the caller's incident amplitude
      and scatter into the global RHS vector.
    - The driven system matrix is `K_driven(ω) = K(ω) − k₀² M(ω) +
      Σ_ABC j k₀ B_ABC + Σ_port j β B_port`; the RHS is `Σ_port
      rhs_port`.
- **Lane:** `crates/yee-fem/src/open_boundary.rs` (create),
  `crates/yee-fem/src/assembly.rs` (face-iteration helper),
  `crates/yee-fem/src/lib.rs` (`pub mod`).
- **Base SHA dep:** E1 + E2 merged.
- **DoD:** unit tests pass — see E4 driven-solve smoke (E3 itself emits no
  numerical result yet, only the driven matrix; the test is "matrix is
  finite, complex-symmetric, has the expected sparsity"). Specifically
  `driven_matrix_is_complex_symmetric` (assemble on a 2-tet fixture and
  assert `K_driven == K_driven^T` to `1e-12`);
  `pec_precedence_over_waveport_at_shared_edges` (a fixture with one
  edge shared between a PEC face and a port face has the PEC tangential
  zero applied, not the modal source).
- **Verification:** `cargo test -p yee-fem --release
  driven_matrix_assembly && cargo clippy -p yee-fem --all-targets --
  -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the local-to-global edge sign
  flip producing a non-symmetric `K_driven` → compare a single-tet
  fixture's `K_driven` against a hand-derived 6×6 reference matrix; the
  orientation flip is exactly the v0 fix from Phase 4.fem.eig.0 D3
  applied to the face blocks. Do not invent a new orientation scheme;
  reuse v0's.
- **LOC:** ~300.

### Step E4 — frequency-sweep driven solve + S-parameter extraction

- **Brief:** Implement `OpenBoundarySolver::sweep(&self, omegas: &[f64])
  -> Result<SParameters, Error>` plus the per-port `S_{11}` extraction.
  Algorithm:
  - For each `omega` in `omegas`:
    - Assemble `K_driven(ω)`, `rhs(ω)` (E3).
    - Factor `K_driven(ω)` once via
      `faer::sparse::FaerLuSolver<Complex64>` and back-substitute against
      `rhs(ω)` (same surface as Phase 4.fem.eig.1 D2 — search-and-
      replace search target).
    - For each port `p`:
      - Project the FEM solution `e` onto port `p`'s modal profile via
        `b_p = 2 ⟨ E_FEM,t , e_mode_p ⟩_port − a_inc_p` (spec §4.3).
      - `S_{p,p}(ω) = b_p / a_inc_p`.
      - Cross-port `S_{p,q}` for multi-port cases lands in 4.fem.eig.2.0.2
        — v0 ships single-port only with stubs returning zero for off-
        diagonal entries.
  - Pack the per-frequency `s` matrices into the returned `SParameters`.
- **Lane:** `crates/yee-fem/src/open_boundary.rs`,
  `crates/yee-fem/src/solve.rs`,
  `crates/yee-fem/tests/driven_solve_pec_stub.rs` (create),
  `crates/yee-fem/tests/abc_face_eats_wave.rs` (create).
- **Base SHA dep:** E3 merged.
- **DoD:** unit tests pass — `driven_solve_pec_stub_returns_full_reflection`
  (drive TE_{10} at 10 GHz into a PEC-closed WR-90 stub, no ABC face;
  expect `|S_{11}| = 0 dB ± 0.05`, phase consistent with the standing-
  wave PEC reflection); `abc_face_eats_wave_at_midband` (drive TE_{10}
  at 10 GHz into a long ABC-terminated stub; expect `|S_{11}| < −35 dB`,
  the documented 1st-order Engquist–Majda floor with margin for
  discretisation). Wall-time < 30 s for both in `--release`.
- **Verification:** `cargo test -p yee-fem --release driven_solve_pec_stub
  abc_face_eats_wave && cargo clippy -p yee-fem --all-targets -- -D
  warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the modal projection returning
  `|S_{11}| > 1` (which is physically impossible) → check the
  `+ j β` vs `− j β` sign and the incident-amplitude normalization
  separately on a hand-built 1-D analog. The standing-wave PEC stub is
  the canary for sign errors; the ABC-terminated stub is the canary for
  the ABC face block sign.
- **LOC:** ~280.

### Step E5 — fem-eig-003 validation gate (WR-90 stub + ABC)

- **Brief:** Implement spec §8 validation. Construct
  `TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, nx, ny, nz)` sized
  to ~25 k tets, ~4 k DoFs after PEC elimination. Classify the face at
  `z = 0` as `Abc`, the face at `z = 30 mm` as `WavePort(0)` with a
  `NumericalCrossSection` solved on the WR-90 cross-section, and the
  four longitudinal sidewalls as `Pec`. Set `a_inc = 1`. Sweep
  `omegas` across `2π · {8.0, 8.08, ..., 12.0} GHz` (50 uniform points,
  80 MHz spacing). Solve the sweep; for each `omega` compute
  `|S_{11}(f)|` in dB and assert: (1) `|S_{11}(f)| ∈ [−45, −35] dB` at
  every swept frequency in 8–12 GHz (the 1st-order Engquist–Majda floor
  is ~ −40 dB ± 5 dB across the band — the ±0.5 dB tolerance in spec §8
  is on the *sweep shape*, not the absolute floor); (2) no swept point
  exceeds `−25 dB` (the would-be PEC resonance at 9.66 GHz is absorbed);
  (3) the phase of `S_{11}(f)` is monotonic across the sweep with no
  unwrapped sign flip; (4) `iterations × omegas.len()` LU factorisations
  complete in `< 180 s` in `--release`. Register
  `run_fem_eig_003_wr90_stub_with_abc` in
  `crates/yee-validation/src/lib.rs` mirroring the existing per-validation
  drivers.
- **Lane:** `crates/yee-validation/src/lib.rs`,
  `crates/yee-validation/tests/fem_eig_003_wr90_stub_with_abc.rs`
  (create), `crates/yee-fem/validation/README.md` (modify).
- **Base SHA dep:** E4 merged (and transitively E1, E2, E3).
- **DoD:** test passes within the spec §8 bounds; `validation/README.md`
  has `fem-eig-003 (WR-90 stub + ABC)` row; CI workflow `ci.yml`
  picks up the new validation test automatically. Wall-time < 180 s in
  `--release`.
- **Verification:** `cargo test -p yee-validation --release
  fem_eig_003_wr90_stub_with_abc` exits 0.
- **Escape hatch:** blocked > 15 min with `|S_{11}|` clipping at the
  `−40 dB` floor across the whole band (which is *expected* for 1st-
  order Engquist–Majda) → green-with-finding and surface
  Phase 4.fem.eig.2.5 as next phase. Do **not** weaken the ±0.5 dB sweep-
  shape tolerance — the *shape* of `|S_{11}(f)|` must still match
  Pozar §3.3 even if the absolute floor is hardware-limited.
- **LOC:** ~320.

### Step E6 (optional) — Python binding `yee.fem.solve_open_cavity(...)`

- **Brief:** Extend `crates/yee-py/src/fem.rs` with
  `solve_open_cavity(mesh, materials, port_faces, abc_faces, omegas) ->
  np.ndarray` returning the swept S-parameter tensor of shape
  `(n_omegas, n_ports, n_ports)`. Mirror the existing
  `yee.fem.solve_cavity_dispersive` binding pattern from Phase 4.fem.eig.1
  D7. Pytest case `crates/yee-py/tests/test_fem_open_boundary.py`
  re-runs a small sweep on the fem-eig-003 geometry and asserts `|S_{11}|
  < −25 dB` across the band.
- **Lane:** `crates/yee-py/src/fem.rs`,
  `crates/yee-py/tests/test_fem_open_boundary.py` (create).
- **Base SHA dep:** E5 merged.
- **DoD:** `maturin develop -p yee-py --release` succeeds; `pytest
  crates/yee-py/tests/test_fem_open_boundary.py` exits 0; the returned
  numpy array has the expected shape and the band-averaged `|S_{11}|`
  is below `−25 dB`.
- **Verification:** `cd crates/yee-py && maturin develop --release &&
  pytest tests/test_fem_open_boundary.py` exits 0.
- **Escape hatch:** blocked > 15 min on PyO3 0.28 returning a 3-D
  complex tensor → ship `(re, im)` paired arrays of shape
  `(n_omegas, n_ports, n_ports, 2)` and surface a
  Phase 4.fem.eig.2.0.3 finding for resolution in a follow-up
  yee-py-lane PR.
- **LOC:** ~200.

## Track sequencing

Critical path: `E1 ‖ E2 → E3 → E4 → E5`.

```
E1 ──┐
     │
     ├── E3 ── E4 ── E5 ──┬── E6 (optional)
     │
E2 ──┘
```

- **E1 and E2 run in parallel** at the start. Both branch off `60ed512`
  and touch disjoint sections of `crates/yee-fem/src/element.rs` (E1: ABC
  face block; E2: port face block + RHS), with disjoint test files.
- **E3 depends on E1 + E2** (consumes both face-block helpers for
  assembly).
- **E4 depends on E3** (frequency sweep + S-parameter extraction
  consumes the driven matrix).
- **E5 depends on E4** (production gate consumes the full pipeline).
- **E6 is optional and depends on E5**.

Within CLAUDE.md §5's "up to 5 parallel agents" envelope: peak
parallelism is at the start (E1 ‖ E2 — two agents). Serial bottleneck is
`E3 → E4 → E5`, ~3 agents-days end-to-end.

## Validation rollup

| Gate | Step | Tolerance | Run-time |
|------|------|-----------|----------|
| **fem-eig-003 sweep shape** — `|S_{11}(f)|` vs Pozar §3.3 across 8–12 GHz | E5 | ±0.5 dB sweep shape; absolute floor in [−45, −35] dB | `< 180 s` `--release` |
| **fem-eig-003 no PEC peak** — would-be 9.66 GHz TE_{101} absorbed | E5 (sub-assertion) | every swept point ≤ −25 dB | covered by same test |
| **fem-eig-003 phase monotonic** — no unwrapped sign flip | E5 (sub-assertion) | phase derivative same sign across band | covered by same test |
| **PEC-closed stub smoke** — ABC off, expect full reflection | E4 unit test | `|S_{11}| = 0 dB ± 0.05` | `< 30 s` `--release` |
| **ABC absorption smoke** — long ABC-terminated stub at mid-band | E4 unit test | `|S_{11}| < −35 dB` | `< 30 s` `--release` |

The sweep-shape and PEC-peak rows land in
`crates/yee-fem/validation/README.md`. Per CLAUDE.md §4 "no solver
feature ships without a published-benchmark validation case" — fem-eig-003
is the published benchmark (Pozar §3.3 dominant-mode characterisation +
Jin §10.4 ABC reflection floor).

Higher-application gates are scoped to later phases:

- **fem-eig-004 (coax-fed dipole inside an ABC box)** — driven Z-input
  validation. Phase 4.fem.eig.2.1.
- **fem-eig-005 (lossy DRA + ABC halo)** — Phase 4.fem.eig.3.

## Lane / file inventory

| Step | Files |
|------|-------|
| E1 | `crates/yee-fem/src/element.rs`, `crates/yee-fem/tests/abc_face_block.rs` (create) |
| E2 | `crates/yee-fem/src/element.rs`, `crates/yee-fem/tests/port_face_block.rs` (create) |
| E3 | `crates/yee-fem/src/{open_boundary,assembly,lib}.rs` |
| E4 | `crates/yee-fem/src/{open_boundary,solve}.rs`, `crates/yee-fem/tests/{driven_solve_pec_stub,abc_face_eats_wave}.rs` (create) |
| E5 | `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fem_eig_003_wr90_stub_with_abc.rs` (create), `crates/yee-fem/validation/README.md` |
| E6 (opt) | `crates/yee-py/src/fem.rs`, `crates/yee-py/tests/test_fem_open_boundary.py` (create) |

Cross-lane consumers (`yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`,
`yee-cuda`) are not touched in 4.fem.eig.2 — `yee-mom`'s
`NumericalCrossSection` is consumed read-only via its existing public API.

## Risk register

Spec §10 risks mapped to steps:

1. **1st-order Engquist–Majda reflection floor (~ −40 dB)** (spec §10).
   This is a *physics* floor, not a bug. **Materialises at Step E5** as
   the gate's absolute-`|S_{11}|` window. The ±0.5 dB tolerance is on
   the *sweep shape* against Pozar §3.3, not the absolute floor; the
   absolute window is the wider `[−45, −35] dB` band. If the gate
   shows the absolute `|S_{11}|` consistently *above* `−35 dB` across
   the band, that's an ABC face-block sign error in E1, not a physics
   floor.
2. **Modal projection on FEM-vs-MoM port mesh mismatch** (spec §10).
   **Materialises at Step E2**. Mitigation: the E2 unit test cross-checks
   the numerical modal RHS against the analytic TE_{10} closed form on
   WR-90; if the projection error there is > 1 %, upgrade interpolation
   to cubic before E5. Do not silently rename the
   `NumericalCrossSection` API.
3. **Port reference plane / phase consistency** (spec §10).
   **Materialises at Step E5** as the phase-monotonicity sub-assertion.
   The reference plane is the port face geometry; if the test fails on
   phase, the mesh has sub-mm port-face placement drift — fix the mesh
   geometry, not the gate.
4. **Complex-symmetric vs Hermitian LU pivoting in `faer`**.
   **Materialises at Step E4** in the per-frequency factor-and-substitute
   loop. Phase 4.fem.eig.1 D2 already exercises this surface on the
   lossy-cavity gate; no new risk, just careful reuse.
5. **PEC corner case at port-face / sidewall edges** (spec §10).
   **Materialises at Step E3** in the face-classifier precedence rule.
   The unit test `pec_precedence_over_waveport_at_shared_edges` is the
   canary.
6. **Sweep frequency density** (spec §10). The 50-point uniform sweep is
   sufficient for the smooth ABC reflection spectrum; adaptive sweeping
   for resonant geometries (fem-eig-004 / 005) is Phase 4.fem.eig.2.1+
   work.

## Out of scope

Explicit non-goals for this plan, per spec §2 and §12:

- **No 2nd-order Engquist–Majda / Higdon ABC.** Phase 4.fem.eig.2.5.
- **No PML (UPML / CFS-PML).** Phase 4.fem.eig.2.5.
- **No multi-mode incident excitation per port.** Single dominant mode
  only.
- **No FEM-BEM hybrid for finite-aperture radiation.** Phase
  4.fem.eig.4+.
- **No driven sweep over the Phase 4.fem.eig.1 dispersive Newton
  tracker.** Phase 4.fem.eig.2.1 — combines both surfaces.
- **No GPU.** CPU-only scalar complex FP64.
- **No higher-order Nedelec on port or ABC faces.** First-order Whitney-1
  only — same as v0/v1.
- **No DRA validation.** Phase 4.fem.eig.3.
- **No CLI / GUI exposure** beyond the optional Python binding (E6).

## Final verification

```bash
cargo build  -p yee-core -p yee-fdtd -p yee-fem -p yee-validation
cargo clippy -p yee-core -p yee-fdtd -p yee-fem -p yee-validation \
  --all-targets -- -D warnings
cargo test   -p yee-core --release
cargo test   -p yee-fdtd --release
cargo test   -p yee-fem  --release
cargo test   -p yee-validation --release fem_eig_003
cargo fmt    --check --all
cargo doc    --no-deps -p yee-core -p yee-fdtd -p yee-fem
mdbook build docs/
```

All nine must exit 0. Every existing `crates/yee-mom/`,
`crates/yee-fdtd/`, `crates/yee-mesh/`, `crates/yee-fem/` test (including
the shipped `fem-eig-001` v0 gate and `fem-eig-002` v1 gate) stays green
— Phase 4.fem.eig.2 is a strict extension, not a refactor.

## Estimated total

- LOC: ~1 470 core (E1 ~150, E2 ~220, E3 ~300, E4 ~280, E5 ~320, E6 ~200).
- Wall-time per agent: 3–5 days end-to-end at one-engineer pace.
  Critical path `E1/E2 → E3 → E4 → E5` is ~3 days; E6 adds ~1 day.
- Risk concentration: Step E5 (fem-eig-003 gate against Pozar §3.3) is
  the load-bearing engineering risk per spec §10; the PEC-closed-stub
  and long-ABC-stub unit tests in E4 are the canaries. Step E2 (modal
  projection against `NumericalCrossSection`) is the load-bearing
  cross-lane interface risk; the analytic TE_{10} cross-check is the
  isolation test.
