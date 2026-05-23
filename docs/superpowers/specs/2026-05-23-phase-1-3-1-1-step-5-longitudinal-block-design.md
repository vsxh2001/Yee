# Phase 1.3.1.1 step 5 — longitudinal (E_z) block + numerical Z_w for quasi-TEM wave-ports

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 5 (ROADMAP "Pending (high priority)").
**Depends on:** steps 2-3 (transverse-only Nedelec assembly + dense
solve, shipped) and the **already-staged** longitudinal element
matrices `local_a_zz` / `local_b_zz` / `local_b_ze` in
`crates/yee-mom/src/eigensolver/assembly.rs` (with unit tests
`local_a_zz_and_b_zz_symmetric_positive_diagonal`,
`local_b_ze_coupling_smoke_test`).
**Blocks:** inhomogeneous (microstrip / partial-fill / CPW) wave-port
β + Z_w extraction; mom-002 / mom-003 wave-port-driven accuracy.

## 1. Goal

Wire the staged **mixed (E_t, E_z)** Lee-Sun-Cendes (1991) longitudinal
block into `NumericalCrossSection::solve`, and replace the TE-mode
wave-impedance *approximation* with a **numerical Z_w extraction** off
the solved eigenvector. This is the capability the whole cross-section
eigensolver exists for: the dominant quasi-TEM mode of an
**inhomogeneous** cross-section (dielectric stack-up), where the
transverse-only solve is wrong because E_z ≠ 0 couples through the
dielectric interface.

## 2. Background — current state

`NumericalCrossSection::solve` (`crates/yee-mom/src/ports.rs:369`)
currently:
- assembles **transverse-only** `A_tt x = β² B_tt x` via
  `assemble_transverse` (`eigensolver/assembly.rs:293`);
- dense-solves via `solve_dense` (`eigensolver/solve.rs:72`),
  returning `beta_sq` + interior-edge eigenvector;
- caches `Z_w ≈ η₀ k₀ / β` (the TE-mode approximation — `ports.rs:384`),
  exact only for the air-filled rectangular guide;
- the longitudinal blocks `local_a_zz` / `local_b_zz` / `local_b_ze`
  exist and are unit-tested but are **unused** (`ports.rs:364-368`).

On a homogeneous air-filled guide the E_z block decouples (the
coupling `B_ze ∝ (ε_r gradient)` vanishes), so the transverse-only
solve is already correct there — which is why the WR-90 TE10 gate
passes today. The longitudinal block only changes the answer on an
**inhomogeneous** cross-section.

## 3. The mixed formulation

After separating `e^{-jβz}`, the inhomogeneous-guide vector
eigenproblem in the `(E_t, E_z)` unknowns is the standard
Lee-Sun-Cendes / Jin (*FEM in EM* 3rd ed. §8.4) block generalized
eigenproblem

```
[ A_tt    0   ] [E_t]       [ B_tt   B_tz ] [E_t]
[  0     A_zz ] [E_z]  = β²  [ B_zt   B_zz ] [E_z]
```

where (per the staged element matrices' construction):
- `A_tt` = `∫ (1/μ_r)(∇_t×N_i)·(∇_t×N_j) − k₀² ε_r N_i·N_j` (the
  existing transverse stiffness, `assemble_transverse`);
- `B_tt` = `∫ (1/μ_r) N_i·N_j` (transverse mass);
- `A_zz` = `local_a_zz` (`∫ (1/μ_r) ∇_t L_i·∇_t L_j − k₀²ε_r L_iL_j`);
- `B_zz` = `local_b_zz` (`∫ (1/μ_r) L_i L_j`);
- `B_tz` = `B_ztᵀ` = `local_b_ze` coupling (`∫ (1/μ_r) N_i·∇_t L_j`).

**The implementer must confirm the exact sign/placement of each block
against the `local_*` docstrings and Lee-Sun-Cendes 1991 eq. (8)-(12)
before assembling** — the staged element matrices encode a specific
convention; do not re-derive from scratch, read what is there.

DoF ordering: stack interior-edge DoFs (E_t, count `n_t` =
`interior_to_global.len()`) above interior-vertex DoFs (E_z, count
`n_z` = interior vertices). Total `n = n_t + n_z`. PEC walls: Dirichlet
on tangential E_t (already eliminated) **and** on E_z (drop
boundary-vertex DoFs — mirror the edge elimination for vertices).

## 4. Approach

### 4.1 `assemble_mixed` (new, in `eigensolver/assembly.rs`)

A new `pub(crate) fn assemble_mixed(mesh, eps_r, mu_r, table) ->
AssembledMixed` that:
1. builds the interior-edge map (reuse the `assemble_transverse`
   logic) **and** an interior-vertex map (new, drop PEC
   boundary vertices);
2. accumulates the four blocks `A_tt, A_zz` into global `A` and
   `B_tt, B_zz, B_tz(=B_ztᵀ)` into global `B`, both `n × n`
   (`n = n_t + n_z`), using the staged `local_*` element matrices;
3. returns `AssembledMixed { a, b, interior_to_global_edges,
   interior_to_global_verts, n_t, n_z }`.

Keep `assemble_transverse` intact (the homogeneous path and its tests
must not regress).

### 4.2 Mixed solve

`solve_dense` currently assumes the transverse pencil. Add
`solve_dense_mixed(asm: &AssembledMixed, freq_hz) -> DenseModeSolution`
(or generalise `solve_dense`) that dense-solves the `n × n` generalized
eigenproblem `A x = β² B x`, selects the dominant **quasi-TEM** mode
(largest physically valid `β²`, i.e. β closest to `k₀√ε_eff`), and
returns `β²` + the full `[E_t; E_z]` eigenvector. Reuse the existing
β²-selection + physical-validity filtering from `solve_dense`.

**Spurious modes.** The mixed formulation admits spurious / non-physical
solutions (the well-known curl-null-space contamination). The standard
mitigation already implicit in the edge/nodal mixed pair (Nedelec for
E_t kills the gradient null space) plus selecting the largest valid β²
is expected to suffice for the dominant mode; if spurious modes appear
**between** the physical β² and the search target, filter by the
`B`-norm energy ratio (a spurious mode has near-zero transverse energy).
Document whatever filter is used.

### 4.3 `NumericalCrossSection::solve` wire-in

Switch `solve` to call `assemble_mixed` + the mixed solve. Scatter the
`[E_t; E_z]` eigenvector back: edge DoFs → `mode_profile` (unchanged
consumer contract for `e_tangential_at`); E_z DoFs cached in a new
`mode_profile_ez: Option<Vec<Complex64>>` field (per-interior-vertex,
scattered to global-vertex indexing) for Z_w extraction + future field
queries. **`e_tangential_at` and the `Numerical2D` wave-port RHS arm
must keep their current behaviour** (the transverse field they consume
is unchanged in form; only its numerical value shifts on inhomogeneous
guides, which is the point).

### 4.4 Numerical Z_w (line integral)

Replace the `Z_w ≈ η₀ k₀ / β` approximation with the power-voltage or
voltage-current definition for a quasi-TEM mode (Jin §8.4 / Pozar
§3.x):
- `V = −∫_path E_t · dl` along a path from the ground conductor to the
  signal conductor (for microstrip: substrate-normal line under the
  strip);
- `Z_w = V² / (2 P)` (power-voltage) **or** the spec §54 energy-ratio
  form `Z_w = (β / ωε₀) · (∫|E_t|² / |∫ E_t·dl|²)`-style expression —
  pick the form that is well-defined from the cached eigenvector and
  document it.
- For the homogeneous air-filled guide this must reduce to the TE-mode
  `η₀ k₀ / β` to within the validation tolerance (regression guard).

## 5. Validation

Per CLAUDE.md §4 — published-benchmark or defensible-physics gate, no
weakening of existing gates.

- **DoD-V1 (regression, homogeneous):** mixed solve on the air-filled
  WR-90 cross-section reproduces the transverse-only β to `< 0.1%` and
  the existing `eigensolver_wr90` TE10 gate stays green. (E_z block
  contributes zero on a homogeneous guide — proves the wire-in did not
  perturb the working path.)
- **DoD-V2 (capability, inhomogeneous — PUBLISHED reference):**
  a **dielectric-slab-loaded rectangular waveguide** (partial vertical
  fill, e.g. WR-90 with the lower half filled by `ε_r = 2.2`). Pozar
  4th ed. §3 / Collin give the published transcendental dispersion
  relation for the dominant mode. Gate: numerical β within **loose
  tolerance** (`≤ 5%`, per the placeholder-tolerance policy that
  governs mom-002/003 until references tighten) of the transcendental
  root, which the test computes by 1-D root-find. If implementing the
  1-D transcendental reference is the >20-min escape-hatch trigger
  (§Risks), fall back to DoD-V2′.
- **DoD-V2′ (capability, physics inequality — fallback):** on the same
  partial-fill guide, assert `k₀ < β_numerical < k₀√ε_r,max` (the mode
  is slow-wave relative to air, fast-wave relative to the densest
  dielectric) **and** regression-track the numerical β. Mirrors the
  original spec's Case B (septum) inequality gate.
- **DoD-V3 (Z_w):** numerical `Z_w` on the air-filled WR-90 reduces to
  `η₀ k₀ / β` within `1%` (regression guard); on the partial-fill guide
  it is finite, positive-real-dominated, and regression-tracked.

## 6. Definition of done

DoD-1. `assemble_mixed` + `AssembledMixed` in `eigensolver/assembly.rs`;
`assemble_transverse` untouched and still green.
DoD-2. Mixed dense solve selecting the dominant quasi-TEM β².
DoD-3. `NumericalCrossSection::solve` uses the mixed path; new
`mode_profile_ez` field; `e_tangential_at` + `Numerical2D` RHS contract
preserved.
DoD-4. Numerical Z_w extraction replacing the TE approximation, reducing
to the TE form on the homogeneous guide.
DoD-5. DoD-V1…V3 green (V2 or the V2′ fallback). New validation test
file `crates/yee-mom/tests/eigensolver_inhomogeneous.rs` (pattern:
`crates/yee-mom/tests/eigensolver_wr90.rs`).
DoD-6. No new `Cargo.toml` dependency (dense solve via the existing
`nalgebra`).
DoD-7. Lint floor clean (`cargo fmt --check --all`,
`cargo clippy --workspace --all-targets -- -D warnings`).
DoD-8. ROADMAP step-5 line marked shipped; ADR-0051 records the
mixed-formulation wire-in + the Z_w definition chosen + the validation
reference choice.

## 7. Risks

(a) **Block sign/placement convention.** The single highest-risk item —
the staged `local_*` matrices encode a specific Lee-Sun-Cendes
convention; assembling the global pencil with a transposed or
sign-flipped coupling block yields wrong β. Mitigation: the homogeneous
regression gate (DoD-V1) catches gross errors (it must reproduce the
known TE10 β); read the element-matrix docstrings, do not re-derive.
(b) **Spurious modes** (§4.2). Mitigation: largest-valid-β² selection +
energy-ratio filter; the dominant quasi-TEM mode is well-separated for
the validation cases.
(c) **Transcendental reference cost** (DoD-V2). Mitigation: the V2′
inequality+regression fallback is a complete gate on its own.
(d) **Z_w path definition** (§4.4) is geometry-dependent (which path).
Mitigation: for the rectangular validation guides the path is the
mid-cross-section vertical line; document it; the homogeneous-reduction
guard (DoD-V3) validates the formula.

## 8. References

* Lee, Sun, Cendes, "Full-wave analysis of dielectric waveguides using
  tangential vector finite elements", IEEE MTT 39(8), 1991.
* Jin, *The Finite Element Method in Electromagnetics*, 3rd ed.,
  §8.4 (inhomogeneous waveguide, mixed E_t/E_z).
* Pozar, *Microwave Engineering* 4th ed., §3 (dielectric-loaded
  waveguide dispersion; quasi-TEM Z_w definitions).
* `crates/yee-mom/src/eigensolver/assembly.rs` — staged `local_a_zz` /
  `local_b_zz` / `local_b_ze` + `assemble_transverse`.
* `crates/yee-mom/src/eigensolver/solve.rs` — `solve_dense`, β²
  selection.
* `crates/yee-mom/src/ports.rs` — `NumericalCrossSection`.
* Cross-section eigensolver design
  `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md`
  (§2 formulation, §validation Case A/B).
