# ADR-0051 — Phase 1.3.1.1 step 5: mixed (E_t, E_z) longitudinal block wire-in

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5 (quasi-TEM microstrip wave-ports)

## Context

`NumericalCrossSection::solve` ships a **transverse-only** Nedelec
eigensolve. That is exact for homogeneous (air-filled) guides but wrong
for the **inhomogeneous** cross-sections the cross-section eigensolver
exists to handle (microstrip, partial dielectric fill, CPW), where the
longitudinal field `E_z` couples through the dielectric interface. The
mixed-formulation longitudinal element matrices (`local_a_zz`,
`local_b_zz`, `local_b_ze`) are already staged and unit-tested in
`eigensolver/assembly.rs` but unused.

Three design questions:

1. **Formulation** — full mixed `(E_t, E_z)` block eigenproblem
   (Lee-Sun-Cendes 1991) vs a scalar quasi-static / transverse-only
   approximation with an ε_eff correction.
2. **Z_w definition** — keep the TE-mode `η₀k₀/β` approximation vs a
   numerical line-integral extraction off the eigenvector.
3. **Inhomogeneous validation reference** — published transcendental
   dielectric-loaded-guide root vs an internal physics inequality.

## Decision

1. **Full mixed `(E_t, E_z)` block formulation** — assemble the staged
   longitudinal blocks into a block generalized eigenproblem with
   `x = [E_t; E_z]`, reusing the existing element matrices. No new
   approximation layer. (As-built: `k_c²`-parameterized and
   symmetric-**indefinite** — see the As-built amendment below; the
   `A x = β² B x` form written here is the design-time prose.)
2. **Numerical Z_w** — voltage line-integral + power, reducing to the
   TE form on the homogeneous guide as a regression guard.
3. **Validation:** published transcendental dielectric-slab-loaded-guide
   β as the primary gate (loose ≤5% tolerance per the
   placeholder-tolerance policy), with a physics-inequality +
   regression fallback (`k₀ < β < k₀√ε_r,max`) mirroring the original
   spec's septum Case B if the transcendental reference trips the
   implementation escape-hatch.

## Rationale

(1) The mixed formulation is the *reason the (E_t, E_z) element
matrices were staged* — steps 2-3 deliberately froze the transverse
path first (walking-skeleton discipline) and deferred the longitudinal
block to here. Using a scalar ε_eff approximation instead would
duplicate Hammerstad-Jensen curve-fits the spec explicitly rejects
(they break for thick conductors / arbitrary stack-ups).

(2) The TE-mode `η₀k₀/β` Z_w is only correct for the air-filled guide;
on a dielectric stack-up the wave impedance is genuinely different and
a line-integral extraction is the standard quasi-TEM definition (Jin
§8.4). Keeping the regression-reduction to the TE form guards the
formula.

(3) CLAUDE.md §4 requires a published benchmark; the dielectric-loaded
guide has one (Pozar §3 / Collin). The inequality fallback keeps the
gate complete and shippable even if the transcendental reference solver
is deferred to step-5.1 — it is the same gate class as the existing
septum Case B.

## Consequences

* New `assemble_mixed` + `AssembledMixed` (block pencil) and
  `solve_dense_mixed` in `eigensolver/`; `assemble_transverse` /
  `solve_dense` retained for the homogeneous path + tests.
* New `mode_profile_ez` field on `NumericalCrossSection`; the
  `e_tangential_at` + `Numerical2D` wave-port RHS contract is preserved
  (the transverse field form is unchanged; its value shifts on
  inhomogeneous guides — the intended effect).
* Highest implementation risk is the block sign/placement convention of
  the staged element matrices. (As-built: the homogeneous-guide β
  regression DoD-V1 does **not** guard the coupling-block sign — the
  coupling decouples globally there; a dedicated horizontal-slab
  `E_z ≠ 0` case guards it. See the As-built amendment.)
* Sparse mixed solve, CPW multi-conductor Z₀ matrix, and the yee-py E_z
  binding are out of scope (later steps).

## As-built amendment (2026-05-23)

Three decisions were refined during implementation; the code is
authoritative where this section and the design prose above differ.

1. **Pencil is `k_c²`-parameterized and symmetric-indefinite.** The
   staged longitudinal element matrices (`local_a_zz`, `local_b_zz`,
   `local_b_ze`) carry no `k₀²` term, so the assembled pencil is
   `A x = k_c² B x` (with `β² = k₀² − k_c²`), matching
   `assemble_transverse`, **not** the `A x = β² B x` of the Decision
   prose. The block `B = [[B_tt, B_tz],[B_zt, B_zz]]` is symmetric
   **indefinite** (the off-diagonal coupling straddles zero even though
   `B_tt`, `B_zz` are individually SPD), so Cholesky / symmetric-
   generalized solvers are invalid. `solve_dense_mixed` forms `B⁻¹A`
   and uses a non-symmetric real-Schur eigensolve with **inverse-iteration**
   eigenvector recovery — acceptable at the `n ≈ 121` validation scale;
   a symmetric-indefinite / LDLᵀ path is the right move when the sparse
   solve lands. (An initial SVD null-space recovery was abandoned: it
   grabbed a spurious E_t-only gradient direction from `A_tt`'s curl
   null-space instead of the physical hybrid mode.)

2. **Coupling-block weight corrected `ε_r → 1/μ_r` — a real bug fix.**
   The originally-staged `local_b_ze` computed `∫ε_r ∇L·N`, a
   divergence-penalty term that the divergence-free curl-null-space mode
   annihilates (Boffi-Brezzi), making the coupling **physically inert on
   every geometry** (`‖E_z‖/‖E_t‖ = 0`). The correct Lee-Sun-Cendes
   coupling is the curl-curl cross term `∫(1/μ_r) ∇L·N` (matching
   `A_zz`); spec §3 prose had the right weight, the staged docstring did
   not. With the fix the coupling is **load-bearing on inhomogeneous
   guides**. The DoD-V1 homogeneous canary cannot guard it — the
   homogeneous dominant mode is pure-TE (`E_z = 0`), weight-independent —
   so the coupling sign/scale is pinned instead by three new guards: a
   **horizontal-slab** `E_z ≠ 0` case (`‖E_z‖/‖E_t‖ = 0.0105`), an
   independent-midpoint-quadrature sign/scale unit test, and a
   zero-`B_tz` β-delta test (β shifts 4.7%). Without the review's demand
   for coupling coverage (P1-1) this no-op coupling would have shipped,
   producing plausible-but-wrong β/Z_w on real microstrip cross-sections.

3. **Validation shipped on DoD-V2′, not the published transcendental.**
   The closed-form dielectric-loaded-guide transcendental reference
   produced no corroborating root in the bring-up window, so the gate
   ships as the rigorous monotonic bracket `β_air < β_loaded < β_full`
   + regression + mesh-convergence (DoD-V2′), with the published
   transcendental reference **queued as step-5.1**. Note the Decision's
   literal band `k₀ < β < k₀√ε_r,max` is wrong for a *closed*
   cutoff-bearing guide (the air-TE10 β already sits below `k₀`); the
   monotonic empty/full bracket is the correct partial-fill statement.

## References

* Lee, Sun, Cendes, IEEE MTT 39(8), 1991.
* Jin, *FEM in EM* 3rd ed. §8.4.
* Pozar, *Microwave Engineering* 4th ed. §3.
* Step 5 spec + plan
  `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-longitudinal-block-design.md`,
  `docs/superpowers/plans/2026-05-23-phase-1-3-1-1-step-5-longitudinal-block.md`.
* ADR-0022 / ADR-0023 — Phase 1.3.1.1 eigensolver spec + TriMesh2D stub.
