# Phase 1.3.1.1 step 5.1 — published transcendental reference for the slab-loaded-guide gate

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5.1 (validation hardening for step 5).
**Depends on:** step 5 (mixed longitudinal block; merge `305d7db`).
**Blocks:** closing the CLAUDE.md §4 published-benchmark requirement for
the inhomogeneous cross-section eigensolver.

## 1. Goal

Replace the **inequality + regression** gate on the step-5 inhomogeneous
cases (`eigensolver_inhomogeneous.rs`) with a comparison against a
**published closed-form transcendental dispersion** for the
slab-loaded rectangular waveguide, so the inhomogeneous β rests on a
real benchmark (CLAUDE.md §4), not a self-referential regression value.

Concretely: implement the **transverse-resonance** dispersion of a
dielectric-slab-loaded rectangular waveguide (Pozar §3.6 / Collin §6 /
Harrington), solve it for β at the validation frequency, and assert the
numerical β matches within a loose (≤5%) tolerance — then tighten if
the agreement is better.

## 2. Background — why step 5 shipped without this

Step 5's DoD-V2′ escape-hatch shipped a monotonic bracket
(`β_air < β_loaded < β_full`) + regression because the prior attempt at
the transcendental reference **found no root matching** the
mesh-converged numerical β (vertical slab β=180.23, horizontal slab
β=201.52). Two hypotheses for that non-corroboration:

(a) **Reference-implementation error** (most likely): wrong mode family
(LSE vs LSM), wrong transverse-resonance bracket, or the `m=1`
x-variation not pinned. The numerical solver is now well-validated — the
homogeneous canary reproduces analytic TE10 to 4e-14, and the coupling
block sign/scale is pinned by an independent-quadrature unit test — so
the *machinery* is trustworthy.

(b) **Numerical β wrong** (unlikely given (a)'s evidence, but must be
ruled out): if the implemented reference is provably correct and still
disagrees, that is a real finding about the mixed solver's inhomogeneous
accuracy, escalated, not papered over.

**This step must distinguish (a) from (b)** — implement the reference
*correctly* (verified against a textbook-tabulated case independent of
our solver) before comparing.

## 3. The physics — slab-loaded-guide transverse resonance

Geometry (horizontal slab, the `coupling_block_loadbearing` case):
rectangular guide width `a` (x), height `b` (y), PEC walls; dielectric
`ε_r` fills `0 ≤ y ≤ d₁`, air fills `d₁ ≤ y ≤ b`. Mode `∝ e^{-jβz}` with
`sin(πx/a)` x-dependence (dominant `m=1`).

The dominant mode of a horizontally-stratified guide is an **LSE_{m0}**
(TE-to-y / "H-mode") whose transverse wavenumbers satisfy
`k_{y,i}² = ε_{r,i} k₀² − (π/a)² − β²` in each layer `i`. The
transverse-resonance condition across the interface + PEC at `y=0,b`
gives the transcendental dispersion (LSE form):

```
k_{y2}·tan(k_{y1} d₁) + k_{y1}·tan(k_{y2} d₂) = 0      (d₂ = b − d₁)
```

(verify the exact form + the LSE-vs-LSM distinction against Pozar §3.6 /
Collin §6 — the dominant slab-loaded mode is LSE for the broad-wall
field orientation here; **the implementer must confirm which family
the numerical dominant mode is and use the matching dispersion**). Solve
for the real `β` root in the propagating window
`(π/a) < β/… ` bracketed below by the air-region cutoff and above by the
fully-filled cutoff — i.e. `β_air < β < β_full`, the same bracket step 5
already uses, **which is where the prior attempt failed to find a root,
so the bracket or the mode family was wrong.**

Cross-check the reference implementation against an **independent
textbook-tabulated** slab-loaded-guide value (e.g. a Pozar example or a
published partially-filled-WR-90 number) before comparing to our solver,
so a reference bug cannot masquerade as a solver bug.

## 4. Approach

1. New `crates/yee-mom/src/eigensolver/reference.rs` (or a test-module
   helper if it stays test-only): `slab_loaded_beta(a, b, d1, eps_r,
   freq_hz) -> f64` solving the LSE transverse-resonance transcendental
   by bracketed root-find (bisection/Brent) on the verified dispersion.
2. **Validate the reference itself** against a textbook-tabulated case
   in a unit test (independent of `NumericalCrossSection`).
3. Reconcile: in `eigensolver_inhomogeneous.rs`, assert the numerical β
   (horizontal slab ε_r=10.2, and vertical slab if the geometry maps to
   a solvable LSM form) matches `slab_loaded_beta` within ≤5%.
   - **If they agree:** tighten the gate — keep the bracket as a
     secondary sanity check, make the transcendental comparison the
     primary published-benchmark gate. Update the test docstrings +
     ROADMAP to mark the §4 gap closed.
   - **If they disagree after the reference is independently verified:**
     do NOT relax. Surface as a finding: bisect (mesh refinement, mode
     family, β-extraction) to localise whether the solver or a
     remaining reference subtlety is at fault. Escalate with the
     measurement; keep the V2′ bracket gate; queue step-5.2.
4. The vertical slab (`ε_r=2.2`, x-stratified) is the LSM/LSE dual; cover
   it if the same `slab_loaded_beta` (with axes swapped) applies, else
   scope it to the horizontal case + note the vertical as follow-on.

## 5. Validation / DoD

- DoD-1. `slab_loaded_beta` reference implemented + **independently
  unit-tested** against a textbook-tabulated slab-loaded-guide β (cite
  the source + page).
- DoD-2. Numerical-vs-reference reconciliation in
  `eigensolver_inhomogeneous.rs` for the horizontal-slab case;
  agreement ≤5% (tighten if better) OR a documented, root-caused
  discrepancy with the V2′ bracket retained.
- DoD-3. If DoD-2 agrees: gate upgraded to published-benchmark; test
  docstrings + ROADMAP + ADR-0052 record the §4 gap closed.
- DoD-4. No new `Cargo.toml` dependency (root-find is hand-rolled or via
  existing deps).
- DoD-5. No regression: the homogeneous canary, WR-90 gate, and the
  coupling guards stay green; existing β regression values unchanged
  unless the tighter gate supersedes them (document the swap).
- DoD-6. Lint floor clean.

## 6. Risks

(a) **Reference still won't corroborate.** Then either the solver has an
inhomogeneous-accuracy issue or the reference has a subtlety we still
miss. Mitigation: the independent textbook-tabulated unit test (DoD-1)
isolates which side is wrong; the bracket gate remains as the floor; the
discrepancy is escalated, not hidden. This is an acceptable outcome —
the *finding* is the value.
(b) **LSE/LSM mode-family confusion** (the likely cause of the prior
failure). Mitigation: determine the numerical dominant mode's
field-orientation (the `‖E_z‖/‖E_t‖` and the E_t polarization from the
solved eigenvector) and match the dispersion family to it.
(c) **Mesh-discretisation error** at ε_r=10.2 (7.5 mm λ₀, high contrast).
Mitigation: report numerical β at 2-3 mesh densities; if the
solver↔reference gap is mesh-limited, that is a legitimate ≤5%-band
result, not a failure.

## 7. References

* Pozar, *Microwave Engineering* 4th ed., §3.6 (partially-filled /
  dielectric-loaded rectangular waveguide; LSE/LSM modes).
* Collin, *Field Theory of Guided Waves* 2nd ed., §6 (inhomogeneously
  filled waveguides, transverse resonance).
* Harrington, *Time-Harmonic Electromagnetic Fields* (slab-loaded guide).
* Step 5 spec + ADR-0051 (the as-built mixed solver + the DoD-V2′
  escape-hatch this step closes).
* `crates/yee-mom/tests/eigensolver_inhomogeneous.rs` — the gate to
  upgrade.
