# Phase 1.3.1.1 step 5.2 — fix the dielectric β-extraction in the cross-section eigensolver

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.2 (solver correctness fix).
**Depends on:** step 5 (mixed solver; `305d7db`), step 5.1 (verified
reference + the disagreement finding; `981926b`).
**Blocks:** closing the CLAUDE.md §4 published-benchmark gap for the
inhomogeneous cross-section eigensolver.

## 1. Goal

Fix the cross-section eigensolver so its dominant-mode β matches the
step-5.1 verified reference (`reference.rs`) and the trivial analytic β
of a uniformly-filled guide. Then tighten the inhomogeneous gate from
the V2′ bracket to a published-benchmark comparison and close the §4 gap.

Tolerance: match reference ≤5% (tighten if better). No gate weakened.

## 2. Root-cause hypothesis (derived; agent must confirm first)

step-5.1 found the FEM β disagrees with a 3-way-verified reference by
~2.9× (solver ε_eff≈1.35 vs reference 8.17 for a half-ε_r=10.2 fill —
physically impossible). The derivation points at the **β-extraction**,
not the assembly weights or the coupling (which are verified correct).

Current code (`solve.rs:85-235` `solve_dense`, mirrored in
`solve_dense_mixed`): forms `S x = k_c² T_ε x` with
`S = ∫(1/μ)(∇×N)(∇×N)`, `T_ε = ∫ε_r N·N` (ε_r-weighted mass), then
extracts

```
β² = k₀² − k_c²        (solve.rs:193, vacuum k₀)
```

documented (assembly.rs:303-306) as "exact when ε_r is real and
**uniform**". The physical transverse vector-Helmholtz equation is

```
∇×(1/μ_r ∇×E_t) = (k₀² ε_r − β²) E_t
⇒ (k₀² T_ε − S) x = β² T_1 x        with T_1 = ∫ N·N  (UNWEIGHTED mass)
```

i.e. the correct eigenproblem has **eigenvalue β² directly** and the RHS
mass **unweighted** (`T_1`), with `ε_r` appearing only on the
`k₀² T_ε` side. The current `β² = k₀² − k_c²` (vacuum `k₀`, single
scalar subtraction, ε_r-weighted RHS) is algebraically equivalent only
when `ε_r ≡ 1`. For any `ε_r ≠ 1` — **uniform or inhomogeneous** — it
under-counts the dielectric. This is why:
- the homogeneous WR-90 canary (ε_r=1) passes to 4e-14 (the two forms
  coincide), yet
- the slab-loaded β comes out near-air (ε_eff≈1.35).

**Cheapest confirmation (do this first):** add a **uniformly-filled**
guide test (ε_r constant, e.g. 2.55) whose analytic
`β = √(ε_r k₀² − (π/a)²)` is trivial (the step-5.1 reference already
anchors the fully-filled limit at ε_r=2.55 → 305.16 rad/m). The current
solver should get this WRONG too — isolating the β-extraction bug from
inhomogeneity and from the coupling block entirely. Confirm, then fix.

## 3. Approach

1. **Confirm** via the uniform-fill test (§2).
2. **Reformulate** the β-extraction. Two equivalent options — pick the
   one that disturbs the least:
   - **(A)** Solve `(k₀² T_ε − S) x = β² T_1 x` directly (eigenvalue is
     β², RHS = unweighted `T_1`). Cleanest; the eigenvalue is the
     physical quantity, no post-hoc relation.
   - **(B)** Keep `S x = k_c² T_ε x` but extract β² with the correct
     ε_r-weighting (`β² = (k₀² − k_c²)·⟨ε_r⟩_field` is **not** generally
     correct because ⟨ε_r⟩ is mode-dependent — option A avoids this
     trap; prefer A).
   Apply consistently to **both** `solve_dense` (transverse) and
   `solve_dense_mixed` (mixed), and to the mixed block pencil's β
   handling.
3. **Preserve** the homogeneous path exactly (ε_r=1 must stay 4e-14).
4. **Re-validate** the indefinite-pencil solve / mode selection still
   picks the dominant mode under the reformulated pencil (the
   `(k₀²T_ε − S)` operator is itself indefinite; mode selection is now
   "largest β²" directly).

## 4. Validation / DoD

- DoD-1. New **uniformly-filled** guide test: numerical β matches
  analytic `√(ε_r k₀² − (π/a)²)` ≤1% (this is the bug's smoking gun and
  the fix's primary anchor — a published/analytic benchmark).
- DoD-2. Inhomogeneous reconciliation (`eigensolver_inhomogeneous.rs`):
  numerical β now matches `reference.rs` (horizontal slab ε_r=10.2,
  vertical slab ε_r=2.2) ≤5%. The diagnostic becomes a **failing gate**
  (the §4 published-benchmark, replacing the V2′ bracket as primary).
- DoD-3. **No regression:** homogeneous WR-90 canary still 4e-14;
  `eigensolver_wr90` green; coupling guards (horizontal-slab ‖E_z‖/‖E_t‖,
  zero-`B_tz` delta) still green; Z_w still reduces to TE form.
- DoD-4. The step-5 regression β values (180.23, 201.52) are **replaced**
  by the corrected values (document the change — they were wrong).
- DoD-5. ROADMAP + ADR-0053 record the fix; the fem-eig cross-section
  §4 gap is **closed**. No new `Cargo.toml` dependency. Lint floor clean.

## 5. Risks

(a) **Reformulation perturbs the homogeneous path.** Mitigation: the
ε_r=1 canary (4e-14) is the guard; option A reduces to the current form
at ε_r=1 by construction.
(b) **Mode selection under the new pencil.** `(k₀²T_ε − S)` is
indefinite; the dominant mode is now the **largest β²** (not smallest
k_c²). Verify the selection + spurious-mode floor translate correctly.
(c) **Z_w depends on β.** The Z_w extraction uses β; confirm it stays
consistent (it already reduces to the TE form on the homogeneous guide;
re-check the loaded value is now physical).
(d) **The reference itself is wrong** (low probability — 3-way verified,
anchored to analytic limits). Mitigation: the uniform-fill analytic
(DoD-1) is a *fully independent* anchor (closed-form, no transverse
resonance); if the fix matches DoD-1 but not the reference, escalate.

## 6. References

* Jin, *FEM in EM* 3rd ed. §8.2-8.4 (cross-section vector eigenproblem,
  the `(k₀²ε_r − β²)` arrangement).
* Pozar §3.3 (uniformly-filled guide analytic β).
* `crates/yee-mom/src/eigensolver/{solve,assembly}.rs` — the β-extraction
  to fix.
* `crates/yee-mom/src/eigensolver/reference.rs` — the validation oracle.
* ADR-0051 / ADR-0052 (step-5 as-built + the step-5.1 finding).
* Memory: `step5-mixed-solver-dielectric-underweight`.
