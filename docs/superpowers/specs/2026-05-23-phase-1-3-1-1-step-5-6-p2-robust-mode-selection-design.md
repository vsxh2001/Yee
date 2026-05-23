# Phase 1.3.1.1 step 5.6 — p=2-robust mode selection (close the high-contrast gate)

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.6 (close the ε_r=10.2 high-contrast inhomogeneous gate).
**Depends on:** step 5.5 (validated p=2 elements; merge `516acec`),
step 5.3 (β-direct sparse shift-invert; `6eca76a`).
**Blocks:** §4 inhomogeneous gate closure across contrasts (FR-4 done;
ε_r=10.2 open).

## 1. Goal

Make `solve_dense_mixed`'s dominant-mode **selection** robust at p=2 so
the ε_r=10.2 high-contrast case selects the physical dominant mode (not
a spurious gradient-cluster mode), closing the gate to ≤5% vs
`reference.rs`. Plus two step-5.5 review carry-ins: a p=2 uniform-fill
analytic anchor (P1-1) and wiring `ElementOrder::Second` through
`NumericalCrossSection`/`ports.rs` so p=2 is reachable end-to-end.

## 2. Background — the precisely-localised blocker (step 5.5 + review)

step 5.5 validated the p=2 elements (J1 independent-quadrature <1e-12,
J3 analytic-TE10) but ε_r=10.2 did not close: p=2 gives ε_eff≈4.8,
*worse* than p1's ≈5.9 (reference 8.17). Reviewer-confirmed root-cause:
`solve_dense_mixed` selects the **smallest transverse-energy-dominated
k_c²**; at p=2 the curl-free gradient edge functions `∇(λ_aλ_b)` enlarge
the near-null cluster with spurious modes that *pass* the
transverse-energy filter, so the selection locks onto a non-dominant
(low-ε_eff, uniformly-spread) mode. p1 is unaffected (its smaller null
space leaves the dominant mode as the smallest valid k_c²); homogeneous
is unaffected (TE10 cleanly separated).

**Key discriminator:** the physical dominant mode concentrates field in
the high-ε region → **high ε_eff**; the spurious gradient-cluster modes
spread uniformly → **low ε_eff** (below the area-average). So ε_eff (or
equivalently the β-direct Rayleigh quotient) distinguishes them.

## 3. Approach

`crates/yee-mom/src/eigensolver/solve.rs`:

1. **ε_eff-screened selection.** Among the transverse-energy-dominated,
   above-spurious-floor candidates, select the one maximising the
   **β-direct Rayleigh quotient** `β² = R(x) = (xᵀ(k₀²B−A)x)/(xᵀB₁x)`
   (equivalently the highest ε_eff = field-weighted permittivity), rather
   than the smallest cutoff k_c². The physical dominant quasi-TEM mode is
   the slowest wave = largest β² = highest ε_eff; the spurious gradient
   modes have low ε_eff. This directly rejects the gradient-cluster
   contamination.
   - Equivalent / complementary framing: keep the β-direct shift-invert
     but **seed σ₀ from a physics-informed estimate** (e.g. the
     area-average or HJ ε_eff → σ₀ = k₀²·ε_eff_est) so the shift-invert
     converges to the physical mode, then screen the result by ε_eff. Pick
     whichever is more robust across the validation cases; document it.
2. **Preserve p1 + homogeneous + FR-4 exactly.** The new selection must
   reduce to the current behaviour where the current behaviour is correct
   (it already picks the right mode at p1 / homogeneous / FR-4). Guard
   with the existing gates (they must stay bit-identical or improve).

`crates/yee-mom/src/ports.rs`:

3. **Wire `ElementOrder::Second` through `NumericalCrossSection`** (a
   constructor option or a `solve` parameter) so the p=2 path is reachable
   end-to-end (currently `assemble_mixed_p2` is lib-internal, only the
   lib tests reach it). First-order stays the default.

`crates/yee-mom/tests/eigensolver_inhomogeneous.rs` + lib tests:

4. **P1-1 — p=2 uniform-fill anchor.** A p=2 uniformly-filled-guide
   (ε_r=2.55) solve matching the analytic `β=√(ε_r k₀²−(π/a)²)` ≤1% —
   closes the ε_r≠1 `assemble_mixed_p2` coverage gap (a `b_tt`/`b_tt1`
   weighted-vs-unweighted-mass swap would currently hide behind the
   wrong-mode failure).
5. **ε_r=10.2 closure gate.** With the fixed selection, p=2 ε_r=10.2 β
   within ≤5% of `reference.rs` 582.95 → a failing gate. Promote the
   step-5.5 documented-finding study to this gate if it closes.

## 4. Validation / DoD

- DoD-1. ε_eff-screened (or physics-seeded) selection in
  `solve_dense_mixed`; p1/homogeneous/FR-4 unchanged (bit-identical or
  improved — the existing gates green).
- DoD-2. P1-1 p=2 uniform-fill anchor ≤1% vs analytic.
- DoD-3. **§4 high-contrast closure:** p=2 ε_r=10.2 β ≤5% vs reference →
  failing gate. If the fixed selection picks the physical mode but p=2
  *still* exceeds 5% (e.g. needs finer mesh than the dense selection
  allows), document the improved convergence + run a finer release
  example; close if achievable, else a narrower finding (the selection
  fix is still validated by the ε_eff recovery toward 8.17).
- DoD-4. `ElementOrder::Second` reachable through `NumericalCrossSection`;
  a Python or Rust end-to-end smoke optional.
- DoD-5. No regression (FR-4, uniform p1 anchor, ε_r=1 canary, wr90,
  coupling guards, Z_w). No new `Cargo.toml` dependency. Lint clean.
  `reference.rs` untouched.

## 5. Risks

(a) **ε_eff-max selection mis-selects on some geometry** (e.g. a higher
mode with incidentally-high ε_eff). Mitigation: combine with the
transverse-energy filter + the propagating-mode (real, positive β²)
constraint; validate on FR-4 (must stay 1.39%) + homogeneous (must stay
TE10) + the vertical slab — the selection must not regress those.
(b) **p=2 still short of ≤5% after correct selection** (dense mesh cap).
Mitigation: DoD-3's finer-release-example branch; even recovering ε_eff
toward 8.17 (from the wrong-mode 4.8) validates the selection fix and is
reportable progress.
(c) **Selection change perturbs p1 gates.** Mitigation: the new rule must
provably reduce to the old where the old is correct; the FR-4/homogeneous
gates are the guard. If it can't preserve them, scope the new selection
to p=2 only (order-conditional).

## 6. References

* Jin, *FEM in EM* 3rd ed. §9 (spurious modes, edge-element null space).
* Boffi-Brezzi-Demkowicz (the gradient/curl-free null space).
* ADR-0054/0055/0056 (the β-direct solve, the diagnosis chain).
* `crates/yee-mom/src/eigensolver/{solve,assembly,reference}.rs`,
  `crates/yee-mom/src/ports.rs`.
