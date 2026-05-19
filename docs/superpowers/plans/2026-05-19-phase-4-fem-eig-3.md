# Phase 4.fem.eig.3 — Coupled-Whitney + 2nd-order ABC + Multi-port — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans`
> to drive this plan step-by-step.

**Companion spec:** `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
**Companion ADR:** `docs/src/decisions/0042-phase-4-fem-eig-3-scope.md`
**Base SHA:** `e45692d` (main HEAD; Phase 4.fem.eig.2 E1-E6 + CCCCCCCCC
partial M_pp normalisation shipped).
**Target phase:** 4.fem.eig.3 only. 4.fem.eig.3.0.2 / 3.1 / 3.5 / 4 are
explicitly deferred — see §"Out of scope".
**Tech-stack additions:** none new. Same `faer::sparse::FaerLuSolver<Complex64>`
surface from v1+v2.

---

## Goal

Phase 4.fem.eig.3 retires the BBBBBBBBB `fem-eig-003` strict
absorption-floor and strict passive-bound `#[ignore]`s and adds multi-port
S-parameter matrix extraction to the open-boundary FEM stack. Three
sub-tracks land in one merge train:

1. **F1+F2 coupled exact-Whitney-1.** Lift the modal RHS
   (`assemble_port_modal_rhs`) and the FEM-projection reconstruction
   (`OpenBoundarySolver::e_t_at_face_centroid`) **together** from the
   lumped `N_i(centroid) ≈ t_i / 3` proxy to the exact
   `N_i(ξ) = λ_a(ξ) ∇λ_b − λ_b(ξ) ∇λ_a` identity evaluated at three
   Gauss points. Round-trip cancellation Pozar §3.3 relies on is
   preserved at the exact-basis level.
2. **F3+F4 2nd-order Engquist-Majda ABC.** Add the tangential-curl
   correction `−(1/2k₀)·(n̂×∇×N_i)·(n̂×∇×N_j)` to the ABC face block;
   gate the new term behind an `AbcOrder` enum so v2 behaviour is
   preserved bit-for-bit on `AbcOrder::First`.
3. **F5+F6 multi-port `S_{p,q}` matrix.** New `sweep_matrix` entry
   point + `SParametersMatrix` output type returning the full
   `n_ports × n_ports` complex matrix per swept frequency. Reuses
   per-frequency LU factor across excited ports.

CPU-only, single-threaded, scalar FP64 complex, no GPU, single dominant
mode per port. Same execution model as v2.

## Pre-flight

Before Step F2 starts, confirm at base SHA `e45692d`:

1. `crates/yee-fem/src/element.rs` exposes
   `assemble_abc_face_block`, `assemble_port_face_block`, and
   `assemble_port_modal_rhs` as the v2 surface; the F1 / F3 helpers are
   *additions*, not modifications. Verify the existing helpers' tests
   in `crates/yee-fem/tests/abc_face_block.rs` and
   `crates/yee-fem/tests/port_face_block.rs` are green; F1/F2/F3 must
   not regress them.
2. `OpenBoundarySolver::extract_s11` carries the CCCCCCCCC `M_pp`
   normalisation. Verify by reading the docstring at line ~762 of
   `crates/yee-fem/src/open_boundary.rs`; the F1+F2 change replaces the
   `e_t_at_face_centroid` lumped reconstruction with `e_t_at_face_gauss_pts`
   but **keeps** the `M_pp` normalisation in `extract_s11` unchanged.
   If the M_pp division has reverted, surface as a CCCCCCCCC regression
   and stop.
3. `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` carries
   two `#[ignore]`'d tests (`fem_eig_003_strict_absorption_floor_gate`
   and `fem_eig_003_strict_passive_bound_continuum_limit`). The Step
   F6 deliverable is removing those `#[ignore]`s; confirm both are
   present at base SHA before F6 lands.
4. `crates/yee-validation/src/lib.rs` defines
   `run_fem_eig_003_wr90_stub_abc` returning a
   `FemEig003Result` struct with `gate_a_floor_ok`, `gate_b_passive_ok`,
   `gate_c_smoothness_ok` fields. F6 adds `run_fem_eig_004_wr90_thru_line`
   and `run_fem_eig_005_wr90_t_junction` drivers next to it.

If (1)-(4) blocks, escape-hatch per CLAUDE.md §5 >15-min rule and surface
as a base-SHA drift finding; do **not** weaken the strict gates.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-fem/src/element.rs` | Modify | Add `assemble_port_face_block_gauss_pts`, `assemble_port_face_rhs_gauss_pts`, `assemble_abc2_face_block`. Existing v2 helpers unchanged. |
| `crates/yee-fem/src/open_boundary.rs` | Modify | `AbcOrder` enum, `with_coupled_whitney`, `with_abc_order`, `sweep_matrix`, `SParametersMatrix`, `e_t_at_face_gauss_pts` private helper. |
| `crates/yee-fem/src/lib.rs` | Modify | Re-export `AbcOrder`, `SParametersMatrix`. |
| `crates/yee-fem/tests/port_face_gauss.rs` | Create | F1 unit test — exact Whitney-1 at Gauss points vs lumped centroid on a non-equilateral fixture. |
| `crates/yee-fem/tests/abc2_face_block.rs` | Create | F3 unit test — 2nd-order ABC block recovers 1st-order in the long-wavelength limit and adds the curl correction at the band edges. |
| `crates/yee-fem/tests/open_boundary_matrix.rs` | Create | F5 unit test — 2-port WR-90 thru-line synthetic; `|S_{21}| ≈ 1`, `S_{12} ≈ S_{21}`. |
| `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` | Modify | F6 — remove `#[ignore]` from the two BBBBBBBBB strict gates after F1-F4 land. |
| `crates/yee-validation/src/lib.rs` | Modify | F6 — add `run_fem_eig_004_wr90_thru_line` and `run_fem_eig_005_wr90_t_junction` drivers. |
| `crates/yee-validation/tests/fem_eig_004_wr90_thru_line.rs` | Create | F6 — fem-eig-004 production gate. |
| `crates/yee-validation/tests/fem_eig_005_wr90_t_junction.rs` | Create | F6 — fem-eig-005 production gate. |
| `crates/yee-fem/validation/README.md` | Modify | Append fem-eig-004 / fem-eig-005 rows. |
| `crates/yee-py/src/fem.rs` | Modify (F7, optional) | `coupled_whitney`/`abc_order`/`multi_port` kwargs. |
| `crates/yee-py/tests/test_fem_multi_port.py` | Create (F7, optional) | Python re-run of fem-eig-004. |
| `docs/src/tutorials/08-fem-multi-port.md` | Create (F8, optional) | mdBook tutorial wiring fem-eig-004 from Python. |

No changes to `yee-mom`, `yee-cuda`, `yee-gui`, `yee-plotters`, `yee-io`,
`yee-cli`, `yee-mesh`. `#![forbid(unsafe_code)]` floor preserved.

## Step ladder

### Step F1 — exact-Whitney-1 wave-port face-block + RHS at Gauss points

- **Brief:** Add `assemble_port_face_block_gauss_pts` and
  `assemble_port_face_rhs_gauss_pts` to
  `crates/yee-fem/src/element.rs`. Implementation:
  - Compute the three face-vertex Whitney-1 gradients
    `g_a = ∇λ_a, g_b = ∇λ_b, g_c = ∇λ_c` on the face plane using the
    same identity as `assemble_tet_element_complex` (constant per-face).
  - At each of the three Gauss points
    `ξ_g ∈ {(2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 1/6, 2/3)}`,
    evaluate `N_i(ξ_g) = λ_a · g_b − λ_b · g_a` per directed edge.
  - Sum the per-Gauss-point contributions weighted `w_g = A / 3` to
    form the 3×3 stiffness block and the 3×1 RHS column.
- **Lane:** `crates/yee-fem/src/element.rs`,
  `crates/yee-fem/tests/port_face_gauss.rs` (create).
- **Pattern file:** `crates/yee-fem/src/element.rs::assemble_port_face_block`
  (v2 lumped form) — preserve the comment style and orientation
  convention.
- **Base SHA dep:** none — branches off `e45692d` directly.
- **DoD:** `port_face_block_gauss_matches_lumped_on_equilateral`
  (an equilateral triangle in the xy-plane: the Gauss-rule output equals
  the v2 lumped output to `1e-12` — equilateral triangles are the
  degenerate case where `t_i / 3 = N_i(centroid)`);
  `port_face_block_gauss_differs_on_right_triangle` (a right triangle:
  the Gauss-rule output **differs** from the v2 lumped output by at
  least `1e-3` in at least one entry — this is the F1 fix's whole
  point);
  `port_face_rhs_gauss_te10_matches_analytic_integral`
  (using analytic TE_{10} profile `ŷ sin(π x / a)` on a known WR-90
  cross-section face, Gauss-rule RHS matches closed-form integral to
  `1e-6`).
- **Verification:** `cargo test -p yee-fem --release port_face_gauss &&
  cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on Whitney-1 gradient sign / face
  orientation → cross-check against `assemble_tet_element_complex`
  (interior tet, identical identity). Do not invent a new gradient
  convention.
- **LOC:** ~240.

### Step F2 — wire `OpenBoundarySolver` to the Gauss-point path

- **Brief:** Extend `crates/yee-fem/src/open_boundary.rs`:
  - Add `coupled_whitney: bool` field on `OpenBoundarySolver` (default
    `false` reproduces v2 + CCCCCCCCC bit-for-bit).
  - Add `with_coupled_whitney(coupled: bool) -> Self` builder method.
  - Add a private `e_t_at_face_gauss_pts(face, e_interior,
    interior_dof_of_edge) -> [Vector3<Complex64>; 3]` helper computing
    the FEM-side `E_FEM(ξ_g)` reconstruction at three Gauss points
    using the same exact Whitney-1 identity as F1.
  - Branch the `scatter_port_face` and `extract_s11` paths on
    `self.coupled_whitney`:
    - `false`: call v2 `assemble_port_modal_rhs` and
      `e_t_at_face_centroid` (existing code, no change).
    - `true`: call F1 `assemble_port_face_rhs_gauss_pts` (caller
      pre-evaluates the modal profile at Gauss points by calling
      `port.modal_e_t` three times per face) and `e_t_at_face_gauss_pts`,
      then accumulate the modal projection across the three Gauss
      points before dividing by `M_pp` (which is also recomputed
      via the same Gauss quadrature).
- **Lane:** `crates/yee-fem/src/open_boundary.rs`,
  `crates/yee-fem/src/lib.rs`.
- **Pattern file:**
  `crates/yee-fem/src/open_boundary.rs::scatter_port_face` and
  `::extract_s11` (v2 lumped form).
- **Base SHA dep:** F1 merged.
- **DoD:** `coupled_whitney_default_false_matches_v2`
  (an `OpenBoundarySolver` built without calling `with_coupled_whitney`
  produces a driven matrix + RHS bit-for-bit identical to v2 on a
  2-tet fixture); `coupled_whitney_true_synthesises_matched_port_zero`
  (a synthetic `E_FEM = a_inc · e_mode` at the port face produces
  `|S_{11}| < 1e-10` with `coupled_whitney = true` — the round-trip
  cancellation now holds at the exact-basis level, not just the
  CCCCCCCCC `M_pp` level).
- **Verification:** `cargo test -p yee-fem --release coupled_whitney
  && cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the per-Gauss-point orientation
  sign vs the v2 centroid sign → run the synthetic matched-port
  fixture with `coupled_whitney = false` (v2 path) and confirm it
  *does not* synthesise `|S_{11}| ≈ 0` (because v2 lumped is buggy on
  non-equilateral triangles — that's the whole point of F1+F2). Then
  flip `coupled_whitney = true` and confirm `|S_{11}| ≈ 0` on the same
  fixture. The sign convention is identical to v2; only the basis
  identity changes.
- **LOC:** ~200.

### Step F3 — 2nd-order Engquist-Majda ABC face-block

- **Brief:** Add `assemble_abc2_face_block` to
  `crates/yee-fem/src/element.rs`. Implementation:
  - Compute the 1st-order contribution identical to
    `assemble_abc_face_block` (re-use the helper or inline it).
  - Compute `∇ × N_i = 2 ∇λ_a × ∇λ_b` per directed edge (constant
    per face — curl of a linear vector field).
  - Form the 2nd-order Gram block
    `R_2[i][j] = (n̂ × ∇×N_i) · (n̂ × ∇×N_j)`.
  - Return the composite block `+jk₀ · (A/μ_r) · R_1 + (−1/(2k₀)) ·
    (A/μ_r) · R_2`.
- **Lane:** `crates/yee-fem/src/element.rs`,
  `crates/yee-fem/tests/abc2_face_block.rs` (create).
- **Pattern file:** `crates/yee-fem/src/element.rs::assemble_abc_face_block`.
- **Base SHA dep:** none — parallel-safe with F1 (disjoint test files).
- **DoD:** `abc2_first_order_part_matches_abc1`
  (extract the imaginary part of `abc2 - abc1`'s `R_1` contribution; equals
  zero to `1e-12`); `abc2_second_order_term_is_real_negative`
  (the real part of `abc2 - abc1` is `−(A/(2k₀μ_r)) · R_2`, which is
  real-symmetric and negative-definite on a TE_{10} face);
  `abc2_low_frequency_limit_dominated_by_first_order`
  (at `k₀ → 0` the 2nd-order term `−(1/2k₀)·R_2 → −∞` would dominate;
  the test pins the *normalised* `1/(jk₀) · abc2 → R_1 − (1/(2k₀²))·R_2`
  and asserts the 2nd-order correction scales as `1/k₀²` — the
  Engquist-Majda 1979 frequency-scaling identity).
- **Verification:** `cargo test -p yee-fem --release abc2_face_block
  && cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the curl identity sign → cite
  Engquist-Majda 1979 eq. 9 DOI 10.1109/TAP.1979.1142175 and the Jin
  3rd-ed §10.4 table 10.1 reflection-floor values; cross-check at a
  single grazing-incidence point (`θ = 60°`) that the 2nd-order floor
  is meaningfully better than 1st-order.
- **LOC:** ~220.

### Step F4 — wire `OpenBoundarySolver` ABC order knob

- **Brief:** Extend `crates/yee-fem/src/open_boundary.rs`:
  - Add `AbcOrder { First, Second }` enum with `Default = First`.
  - Add `abc_order: AbcOrder` field on `OpenBoundarySolver` (default
    `First` reproduces v2 bit-for-bit).
  - Add `with_abc_order(order: AbcOrder) -> Self` builder method.
  - Branch `scatter_abc_face` on `self.abc_order`:
    - `First`: call v2 `assemble_abc_face_block` (existing code).
    - `Second`: call F3 `assemble_abc2_face_block`.
- **Lane:** `crates/yee-fem/src/open_boundary.rs`,
  `crates/yee-fem/src/lib.rs`,
  `crates/yee-fem/tests/abc2_face_block.rs` (extend with a small
  integration sub-test).
- **Pattern file:** `crates/yee-fem/src/open_boundary.rs::scatter_abc_face`.
- **Base SHA dep:** F3 merged.
- **DoD:** `abc_order_default_first_matches_v2`
  (an `OpenBoundarySolver` built without calling `with_abc_order`
  produces an ABC scatter bit-for-bit identical to v2 on a 2-tet
  fixture); `abc_order_second_lowers_reflection_on_normal_incidence`
  (on the v2 ABC-eats-wave synthetic fixture from
  `crates/yee-fem/tests/open_boundary_sweep.rs::abc_face_eats_wave`,
  `AbcOrder::Second` produces `|S_{11}|` at least 6 dB lower than
  `AbcOrder::First` at mid-band).
- **Verification:** `cargo test -p yee-fem --release abc_order &&
  cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the 2nd-order term making the
  driven matrix singular near a band edge → fall back to `AbcOrder::First`
  on that frequency and surface a Phase 4.fem.eig.3.5 finding. The
  WR-90 stub fixture should not trigger this on 8-12 GHz with WR-90's
  `f_c = 6.56 GHz`.
- **LOC:** ~140.

### Step F5 — multi-port `S_{p,q}` matrix extraction

- **Brief:** Extend `crates/yee-fem/src/open_boundary.rs`:
  - Add `SParametersMatrix { omegas: Vec<f64>, s: Vec<DMatrix<Complex64>> }`
    output type.
  - Add `sweep_matrix(&self, omegas: &[f64]) -> Result<SParametersMatrix,
    Error>`. Algorithm per spec §7:
    - Per `ω`: assemble the driven system, factor once.
    - For each excited port `p ∈ 0..n_ports`: build a port-specific
      RHS with `a_inc_p = 1` and `a_inc_q = 0` for `q ≠ p`,
      back-substitute against the same LU factor.
    - For each pair `(q, p)`: project the FEM solution onto port `q`'s
      modal profile; pack into `s[k][(q, p)]`.
    - For multi-port modal-overlap correction (spec §10), compute the
      full `M_{pq} = ⟨e_mode_p, e_mode_q⟩_port` matrix at each
      frequency. For geometrically-disjoint ports this is diagonal and
      the per-port extraction reduces to the single-port formula.
      For overlapping ports, solve the small `n_ports × n_ports`
      linear system to recover `b = M^{-1} · ⟨E_FEM, e_mode⟩ − a_inc`.
      Warn if `cond(M) > 1e6`.
  - The single-port `sweep` entry point stays for v2 callers.
- **Lane:** `crates/yee-fem/src/open_boundary.rs`,
  `crates/yee-fem/src/lib.rs`,
  `crates/yee-fem/tests/open_boundary_matrix.rs` (create).
- **Pattern file:** `crates/yee-fem/src/open_boundary.rs::sweep` (v2
  single-port form).
- **Base SHA dep:** F2 + F4 merged.
- **DoD:** `sweep_matrix_two_port_thru_line_reciprocal_and_unitary`
  (synthetic 2-port WR-90 thru-line at 10 GHz: `|S_{21}| ∈ [0.9, 1.05]`,
  `|S_{12} − S_{21}| < 1e-4`, `|S_{11}| + |S_{12}|² ≤ 1.05`);
  `sweep_matrix_lu_factor_reused_across_excited_ports`
  (timing assertion: 2-port `sweep_matrix` runtime is < 1.5× the
  single-port `sweep` runtime at the same `omegas`).
- **Verification:** `cargo test -p yee-fem --release open_boundary_matrix
  && cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on modal-overlap matrix singular
  → log `cond(M)` and fall back to diagonal extraction for the
  fem-eig-004 thru-line case (the two end-face ports are geometrically
  disjoint; `M` is diagonal modulo numerical noise). The
  `cond(M) > 1e6` warning is the canary.
- **LOC:** ~280.

### Step F6 — fem-eig-003 strict un-ignore + fem-eig-004 + fem-eig-005

- **Brief:** Production validation:
  - **fem-eig-003 strict.** In
    `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`:
    update `run_fem_eig_003_wr90_stub_abc` (in
    `crates/yee-validation/src/lib.rs`) to set
    `solver = solver.with_coupled_whitney(true).with_abc_order(AbcOrder::Second)`,
    then remove the two `#[ignore]` attributes from
    `fem_eig_003_strict_absorption_floor_gate` and
    `fem_eig_003_strict_passive_bound_continuum_limit`.
  - **fem-eig-004.** Add
    `run_fem_eig_004_wr90_thru_line` to
    `crates/yee-validation/src/lib.rs`: 60 mm WR-90 section with both
    end faces tagged `WavePort(p)`, four sidewalls PEC,
    `coupled_whitney = true`, `abc_order = First` (no ABC faces).
    At 10 GHz, call `sweep_matrix(&[2π·10e9])` and assert spec §8
    gate criteria.
  - **fem-eig-005.** Add `run_fem_eig_005_wr90_t_junction` to
    `crates/yee-validation/src/lib.rs`: WR-90 H-plane T-junction
    (~50 k tets), three TE_{10} ports, no ABC faces. At 5 GHz,
    assert magnitude conservation and reciprocity per spec §8.
- **Lane:**
  `crates/yee-validation/src/lib.rs`,
  `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`,
  `crates/yee-validation/tests/fem_eig_004_wr90_thru_line.rs` (create),
  `crates/yee-validation/tests/fem_eig_005_wr90_t_junction.rs` (create),
  `crates/yee-fem/validation/README.md`.
- **Pattern file:**
  `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  (BBBBBBBBB form — preserve the gate decomposition style).
- **Base SHA dep:** F5 merged (and transitively F1, F2, F3, F4).
- **DoD:** All three production gates pass:
  - `fem_eig_003_strict_absorption_floor_gate` (un-ignored) passes
    within `[-45, -35] dB`.
  - `fem_eig_003_strict_passive_bound_continuum_limit` (un-ignored)
    passes with `|S_{11}| < 1` strict.
  - `fem_eig_004_thru_line_at_10ghz` passes within the spec §8
    fem-eig-004 gates.
  - `fem_eig_005_t_junction_at_5ghz` passes within the spec §8
    fem-eig-005 gates.
  - `crates/yee-fem/validation/README.md` has fem-eig-004 +
    fem-eig-005 rows.
- **Verification:** `cargo test -p yee-validation --release
  fem_eig_003 fem_eig_004 fem_eig_005` exits 0.
- **Escape hatch:** blocked > 15 min on fem-eig-005 reciprocity gate
  → diagnose via the modal-overlap condition-number log; if
  `cond(M) > 1e6` at the T-junction, the three port profiles have
  non-trivial inner product in the FEM-projected basis. Apply the
  `M^{-1}` correction from F5; if still failing, widen the
  reciprocity tolerance to `1e-2` and surface a Phase 4.fem.eig.3.0.2
  finding for the multi-mode incident-excitation upgrade.
- **LOC:** ~480.

### Step F7 (optional) — Python multi-port binding

- **Brief:** Extend `crates/yee-py/src/fem.rs`:
  `yee.fem.solve_open_cavity(..., coupled_whitney=False,
  abc_order="first", multi_port=False) -> np.ndarray`. When
  `multi_port=True`, return shape `(n_omegas, n_ports, n_ports)`.
  Pytest case re-runs fem-eig-004 from Python and asserts
  `|S_{21}| ≈ 1`.
- **Lane:** `crates/yee-py/src/fem.rs`,
  `crates/yee-py/tests/test_fem_multi_port.py` (create).
- **Base SHA dep:** F6 merged.
- **DoD:** `maturin develop --release` succeeds; pytest passes.
- **Verification:** `cd crates/yee-py && maturin develop --release &&
  pytest tests/test_fem_multi_port.py` exits 0.
- **Escape hatch:** PyO3 0.28 returning a 3-D complex tensor → ship
  `(re, im)` paired 4-D arrays.
- **LOC:** ~180.

### Step F8 (optional) — mdBook tutorial

- **Brief:** Add `docs/src/tutorials/08-fem-multi-port.md` walking
  through fem-eig-004 from Python end-to-end.
- **Lane:** `docs/src/tutorials/08-fem-multi-port.md` (create),
  `docs/src/SUMMARY.md` (link entry).
- **Base SHA dep:** F7 merged.
- **DoD:** `mdbook build docs/` exits 0; tutorial renders cleanly.
- **LOC:** ~180 prose.

## Track sequencing

Critical path: `F1 → F2  ‖  F3 → F4  →  F5 → F6 → F7 → F8`.

```
F1 ─→ F2 ──┐
           │
F3 ─→ F4 ──┼─→ F5 ─→ F6 ─→ F7 (opt) ─→ F8 (opt)
           │
```

- **F1 and F3 run in parallel** at the start. Both branch off `e45692d`
  and touch disjoint sections of `crates/yee-fem/src/element.rs` (F1:
  Gauss-point port helpers; F3: 2nd-order ABC) with disjoint test files.
- **F2 depends on F1** and **F4 depends on F3** — wiring steps.
- **F5 depends on F2 + F4** (matrix sweep consumes both knobs).
- **F6 depends on F5** (production gates consume the full surface).
- **F7 + F8 are optional and depend on F6**.

Peak parallelism: 2 agents (F1 ‖ F3) at start. Serial bottleneck:
`F2 / F4 → F5 → F6`. End-to-end ~4 agent-days at one-engineer pace.

## Validation rollup

| Gate | Step | Tolerance | Run-time |
|------|------|-----------|----------|
| **fem-eig-003 strict absorption floor** — `|S_{11}(f)|` ∈ `[-45, -35] dB` across 8-12 GHz with v3 flags on | F6 | per spec §8 fem-eig-003 strict | `< 240 s` `--release` |
| **fem-eig-003 strict passive bound** — `|S_{11}(f)| < 1` strict | F6 | per spec §8 fem-eig-003 strict | covered by same driver |
| **fem-eig-004 thru-line `|S_{21}|`** — vs lossless WR-90 | F6 | `|S_{21}| ∈ [0.95, 1.05]`, `|S_{12} − S_{21}| < 1e-6`, `|S_{11}|, |S_{22}| < −30 dB` | `< 180 s` `--release` |
| **fem-eig-005 T-junction conservation + reciprocity** | F6 | `Σ_q |S_{q,p}|² ∈ [0.95, 1.05]`, `|S_{p,q} − S_{q,p}| < 1e-3` | `< 300 s` `--release` |
| **F1 coupled Whitney unit** — Gauss-point vs lumped on non-equilateral | F1 | `> 1e-3` entry difference on right triangle | `< 5 s` |
| **F3 2nd-order ABC unit** — 1st + curl correction | F3 | per spec §4.2 / Engquist-Majda 1979 eq. 9 | `< 5 s` |
| **F5 2-port thru-line synthetic** — reciprocity + unitarity | F5 | `|S_{12} − S_{21}| < 1e-4` | `< 30 s` |

The three production rows land in `crates/yee-fem/validation/README.md`.

## Lane / file inventory

| Step | Files |
|------|-------|
| F1 | `crates/yee-fem/src/element.rs`, `crates/yee-fem/tests/port_face_gauss.rs` (create) |
| F2 | `crates/yee-fem/src/open_boundary.rs`, `crates/yee-fem/src/lib.rs` |
| F3 | `crates/yee-fem/src/element.rs`, `crates/yee-fem/tests/abc2_face_block.rs` (create) |
| F4 | `crates/yee-fem/src/open_boundary.rs`, `crates/yee-fem/src/lib.rs`, `crates/yee-fem/tests/abc2_face_block.rs` |
| F5 | `crates/yee-fem/src/open_boundary.rs`, `crates/yee-fem/src/lib.rs`, `crates/yee-fem/tests/open_boundary_matrix.rs` (create) |
| F6 | `crates/yee-validation/{src,tests}/...`, `crates/yee-fem/validation/README.md` |
| F7 (opt) | `crates/yee-py/src/fem.rs`, `crates/yee-py/tests/test_fem_multi_port.py` (create) |
| F8 (opt) | `docs/src/tutorials/08-fem-multi-port.md` (create), `docs/src/SUMMARY.md` |

Cross-lane consumers (`yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`,
`yee-cuda`, `yee-plotters`) are not touched in 4.fem.eig.3.

## Risk register

Spec §10 risks mapped to steps:

1. **3-point Gauss quadrature degree.** Materialises at **Step F1**.
   Mitigation: drop in a 6-point rule if F1 unit-test convergence is
   marginal vs the analytic TE_{10} integral.
2. **2nd-order Mur stability near band edges.** Materialises at **Step
   F4 / F6**. Mitigation: per spec §10, fall back to `AbcOrder::First`
   on per-frequency basis if `|β_mode(ω) − k₀| / k₀ > 0.5`. The
   fem-eig-003 8-12 GHz band is well above the `f_c = 6.56 GHz` WR-90
   cutoff, so this should not trigger.
3. **Multi-port modal-overlap ill-conditioning.** Materialises at **Step
   F5 / F6** (fem-eig-005). Mitigation: compute `cond(M)` per frequency,
   log a warning above `1e6`, apply `M^{-1}` projection correction.
4. **Excited-port LU-factor reuse correctness.** Materialises at **Step
   F5**. Verified via the `sweep_matrix_lu_factor_reused_across_excited_ports`
   timing assertion + the `coupled_whitney_default_false_matches_v2`
   bit-for-bit identity test.
5. **CCCCCCCCC `M_pp` normalisation regression risk.** Materialises at
   **Step F2** when the FEM-projection reconstruction is rewritten.
   Mitigation: F2's `coupled_whitney_default_false_matches_v2` test
   ensures the v2 path is preserved bit-for-bit. The new `M_pp`
   computation under `coupled_whitney = true` uses the same Gauss
   quadrature as the modal-projection numerator, preserving the
   round-trip identity.

## Out of scope

Explicit non-goals for this plan, per spec §2 and §12:

- **No CFS-PML / UPML.** Phase 4.fem.eig.3.5.
- **No multi-mode incident excitation per port.** Phase 4.fem.eig.3.0.2.
- **No higher-order Nedelec.** Phase 4.fem.eig.4+.
- **No driven sweep over Phase 4.fem.eig.1 dispersive Newton tracker.**
  Phase 4.fem.eig.3.1 — combines both surfaces.
- **No GPU.** CPU-only scalar complex FP64.
- **No adaptive sweep / model-order reduction.** Uniform sweep only.
- **No CLI / GUI exposure** beyond the optional Python binding (F7).

## Final verification

```bash
cargo build  -p yee-core -p yee-fdtd -p yee-fem -p yee-validation
cargo clippy -p yee-core -p yee-fdtd -p yee-fem -p yee-validation \
  --all-targets -- -D warnings
cargo test   -p yee-core --release
cargo test   -p yee-fdtd --release
cargo test   -p yee-fem  --release
cargo test   -p yee-validation --release \
             fem_eig_003 fem_eig_004 fem_eig_005
cargo fmt    --check --all
cargo doc    --no-deps -p yee-core -p yee-fdtd -p yee-fem
mdbook build docs/
```

All nine must exit 0. Every existing `yee-mom`, `yee-fdtd`, `yee-mesh`,
`yee-fem` test (including `fem-eig-001` v0, `fem-eig-002` v1, and the
v2 + CCCCCCCCC `fem-eig-003` non-strict smoke + gate B + gate C) stays
green — Phase 4.fem.eig.3 is a strict extension under defaulted-off
config knobs.

## Estimated total

- LOC: ~1 740 core (F1 ~240, F2 ~200, F3 ~220, F4 ~140, F5 ~280, F6
  ~480, F7 ~180 opt, F8 ~180 opt prose).
- Wall-time per agent: 4-6 days end-to-end at one-engineer pace.
  Critical path `F1/F3 → F2/F4 → F5 → F6` is ~4 days; F7+F8 add ~2.
- Risk concentration: Step F2 (coupled Whitney wiring — the actual
  blocker CCCCCCCCC could not retire), Step F4 (2nd-order Mur band-edge
  stability), Step F6 fem-eig-005 (multi-port modal-overlap
  conditioning). The F5 timing assertion is the canary for LU-factor
  reuse correctness.
