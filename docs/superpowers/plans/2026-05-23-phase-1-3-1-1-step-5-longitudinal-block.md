# Phase 1.3.1.1 step 5 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-longitudinal-block-design.md`
**Base SHA:** `c8cd846` (post step-4 LOBPCG merge)
**Lane:** `crates/yee-mom/src/eigensolver/**`,
`crates/yee-mom/src/ports.rs`,
`crates/yee-mom/tests/eigensolver_inhomogeneous.rs` (new),
`ROADMAP.md`, `docs/src/decisions/0051-*.md`.
**Out of lane** (surface as findings, do NOT fix): `crates/yee-fem/**`
(the 3D cavity eigensolver — unrelated), `crates/yee-py/**`, the
existing `eigensolver_wr90.rs` test (run it read-only as a regression
gate, do not edit), `WavePort::rhs` semantics beyond keeping the
`Numerical2D` arm compiling.

## Step ladder

### S1 — `assemble_mixed` + `AssembledMixed`

`crates/yee-mom/src/eigensolver/assembly.rs`:
1. Add interior-vertex map (drop PEC boundary vertices), mirroring the
   interior-edge map in `assemble_transverse`.
2. Add `pub(crate) struct AssembledMixed { a, b, interior_to_global_edges, interior_to_global_verts, n_t, n_z }`.
3. Add `pub(crate) fn assemble_mixed(...)` accumulating `A_tt`,`A_zz`
   into `A` and `B_tt`,`B_zz`,`B_tz=B_ztᵀ` into `B`, using the staged
   `local_a_zz`/`local_b_zz`/`local_b_ze` + the existing transverse
   element matrices. **Read each `local_*` docstring first** to fix the
   block sign/placement (spec §3, §7a). Edge DoFs stacked above vertex
   DoFs.
4. `assemble_transverse` untouched.

**Verification:** `cargo test -p yee-mom --lib eigensolver::assembly`
exit 0 (existing element-matrix unit tests stay green; add a smoke that
`assemble_mixed` on a 2-triangle homogeneous patch yields a `B_tz`
block that is zero when ε_r is uniform — the decoupling sanity check).

### S2 — mixed dense solve

`crates/yee-mom/src/eigensolver/solve.rs`:
1. Add `solve_dense_mixed(&AssembledMixed, freq_hz)` (or generalise
   `solve_dense`): dense generalized eigensolve `A x = β² B x` via the
   existing `nalgebra` path, dominant-quasi-TEM β² selection (largest
   valid β²; energy-ratio spurious-mode filter per spec §4.2), return
   β² + full `[E_t;E_z]` eigenvector.
2. Keep `solve_dense` (transverse) intact for the homogeneous path /
   its tests.

**Verification:** `cargo test -p yee-mom --lib eigensolver::solve`
exit 0.

### S3 — `NumericalCrossSection::solve` wire-in + `mode_profile_ez`

`crates/yee-mom/src/ports.rs`:
1. Add `pub mode_profile_ez: Option<Vec<Complex64>>` (global-vertex
   indexing), `None` in `new`.
2. Switch `solve` to `assemble_mixed` + `solve_dense_mixed`; scatter
   edge DoFs → `mode_profile` (unchanged contract), vertex DoFs →
   `mode_profile_ez`.
3. Preserve `e_tangential_at` + the `Numerical2D` RHS arm behaviour.

**Verification:** `cargo test -p yee-mom --test eigensolver_wr90`
(read-only gate) stays green — homogeneous WR-90 β unchanged within
DoD-V1 0.1%.

### S4 — numerical Z_w

`crates/yee-mom/src/ports.rs`:
1. Replace the `Z_w ≈ η₀k₀/β` line with a numerical extraction
   (spec §4.4): voltage line-integral + power, documented path. Must
   reduce to `η₀k₀/β` within 1% on the homogeneous guide (DoD-V3).

**Verification:** assertion inside the new test (S5).

### S5 — inhomogeneous validation

New `crates/yee-mom/tests/eigensolver_inhomogeneous.rs` (pattern:
`eigensolver_wr90.rs`):
1. DoD-V1: homogeneous WR-90 mixed-solve β reproduces transverse β
   `< 0.1%`.
2. DoD-V2 (published) **or** DoD-V2′ (inequality+regression fallback):
   dielectric-slab-loaded guide. Try the transcendental reference; if
   it trips the >20-min escape-hatch, ship V2′ (`k₀ < β < k₀√ε_r,max`
   + regression value) and note it.
3. DoD-V3: Z_w reduces to TE form on the homogeneous guide; finite +
   regression-tracked on the loaded guide.

**Verification:**
`cargo test -p yee-mom --test eigensolver_inhomogeneous` exit 0.

### S6 — ROADMAP + ADR-0051

1. `ROADMAP.md`: mark Phase 1.3.1.1 step 5 shipped with merge SHA.
2. ADR-0051 `docs/src/decisions/0051-phase-1-3-1-1-step-5-longitudinal-block.md`:
   record the mixed-formulation block placement decision, the Z_w
   definition chosen, and the validation-reference choice (V2 vs V2′).

**Verification:** `mdbook build docs/` exit 0; `grep -n "step 5"
ROADMAP.md` shows the shipped marker.

## Lint floor (every commit)

```bash
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
```

## Full verification (before declaring done; all exit 0)
```
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p yee-mom --lib eigensolver
cargo test -p yee-mom --test eigensolver_wr90          # read-only regression
cargo test -p yee-mom --test eigensolver_inhomogeneous # new
git diff --stat -- '**/Cargo.toml'                     # expect EMPTY
```

## Escape-hatch

Blocked > 20 min on either (a) the block sign/placement (homogeneous
regression β wrong) or (b) the transcendental DoD-V2 reference → for
(a): bisect by zeroing the coupling block `B_tz` and confirming the
mixed solve then equals the transverse solve, then re-introduce
coupling; surface the convention as a finding if still wrong. For (b):
ship the DoD-V2′ inequality+regression gate and queue the transcendental
reference as step-5.1. Do NOT weaken DoD-V1 (homogeneous regression) or
the existing `eigensolver_wr90` gate — those are correctness floors.
If the mixed solve cannot reproduce the homogeneous β at all, commit
S1 (assembly + decoupling smoke) only, leave `NumericalCrossSection::solve`
on the transverse path, mark the inhomogeneous gate `#[ignore]`, and
surface the blocker.

## Out-of-scope (surface as findings, do not fix)

* Sparse mixed solve (LOBPCG/`SparseEigen` for the 2-D cross-section) —
  the dense path is fine at validation DoF counts; sparse is a later step.
* CPW / multi-conductor Z₀ matrix (multi-port modal) — step 5.2+.
* yee-py binding for the E_z profile / numerical Z_w.
* `yee-fem/**` 3D cavity eigensolver (different solver entirely).
