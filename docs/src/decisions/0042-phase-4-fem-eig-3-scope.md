# ADR-0042: Phase 4.fem.eig.3 scope — coupled-Whitney + 2nd-order ABC + multi-port

## Status

Accepted — 2026-05-19 (spec + plan; implementation deferred to
follow-up tracks).

## Context

Phase 4.fem.eig.2 shipped the open-boundary FEM walking skeleton on
top of the Phase 4.fem.eig.0/1 closed-cavity stack (ADR-0040):
1st-order Engquist-Majda ABC + single-mode wave-port + diagonal
single-port S-parameter extraction. The BBBBBBBBB `fem-eig-003` strict
absorption-floor gate `[-45, -35] dB` and the strict continuum-limit
passive-bound gate `|S_{11}| < 1` both saturated at `|S_{11}| = 1.0`
on the WR-90 stub fixture and remain `#[ignore]`'d under the Phase
4.fem.eig.2 E5 escape hatch.

Track CCCCCCCCC (2026-05-19) identified the modal-projection formula
`b_p = 2 ⟨E_FEM, e_mode⟩ − a_inc` in `OpenBoundarySolver::extract_s11`
silently assumed `M_pp = ⟨e_mode, e_mode⟩_port = 1/2`, while the
standard Pozar §3.3 orthonormalisation used by the
`fem-eig-003` driver carries `M_pp ≈ 1`. CCCCCCCCC normalised the
projection by `M_pp` (replacing `b_p = 2·⟨·,·⟩ − a_inc` with
`b_p = ⟨·,·⟩ / M_pp − a_inc`), which **retired the synthetic
matched-port identity** `E_FEM = a_inc · e_mode  ⇒  S_{11} = 0` but
**did not retire** the empirical `|S_{11}| ≈ 1.0` saturation on the
production fixture. The CCCCCCCCC root-cause analysis traced the
residual saturation to the **lumped Whitney-1 basis-at-centroid
approximation** `N_i(centroid) ≈ t_i / 3` used in
`element::assemble_port_modal_rhs` and
`OpenBoundarySolver::e_t_at_face_centroid`. The exact Whitney-1
identity `N_i(centroid) = (1/3)(∇λ_b − ∇λ_a)` requires per-element
gradient computation — and on every non-equilateral face (every
Kuhn-decomposed face on the WR-90 stub) the lumped form mis-evaluates
`N_i(centroid)`, driving the reconstructed `E_FEM(centroid)` toward zero
on the port face regardless of the FEM solution amplitude.

CCCCCCCCC prototyped a coupled fix (lifting both the RHS and the
projection together to the exact Whitney-1 basis) but found
over-amplification near the closed-stub TE_{10n} resonances at 8 GHz
(`n=1`) and 12 GHz (`n=2`) on the fem-eig-003 stub. The fix needs either
(a) per-Gauss-point modal sampling instead of single-centroid
quadrature, or (b) higher-order Nedelec basis at the port face. (b) is
out of scope for the Whitney-1 floor; (a) is the v3 deliverable.

Independently, the 1st-order Engquist-Majda ABC's `~ −40 dB` reflection
floor was accepted in ADR-0040 §C-3 as the v2 physics limit, with PML
/ 2nd-order ABC reserved as the Phase 4.fem.eig.2.5 upgrade slot. Phase
4.fem.eig.3 takes that slot and merges it with the coupled-Whitney
fix and the multi-port extension into one merge train, because all
three changes touch the same `OpenBoundarySolver` config surface and
the same `element.rs` face-block helpers.

Multi-port `S_{p,q}` extraction is the natural extension of the
single-port `S_{p,p}` v2 path: Sheen, Ali, Abouzahra, Katehi 1990
(*IEEE Trans. MTT* 38(7) p. 849) describes the column extraction via
per-excited-port driven solves with shared LU factor across columns.

## Decision

Phase 4.fem.eig.3 ships **three coupled sub-tracks** behind defaulted-off
config knobs so the v2 + CCCCCCCCC behaviour stays bit-for-bit
identical:

1. **F1+F2 — coupled exact-Whitney-1 RHS + projection at 3-point Gauss
   quadrature.** New `element::assemble_port_face_block_gauss_pts` and
   `assemble_port_face_rhs_gauss_pts` consuming the modal `E_t`
   pre-sampled at three Gauss points on the reference triangle.
   `OpenBoundarySolver::with_coupled_whitney(bool)` toggles the
   exact-Whitney path; default `false` reproduces the v2 + CCCCCCCCC
   lumped-centroid behaviour.
2. **F3+F4 — 2nd-order Engquist-Majda ABC.** New
   `element::assemble_abc2_face_block` carrying both the 1st-order
   Mur term and the 2nd-order tangential-curl correction
   `−(1/2k₀)·(n̂×∇×N_i)·(n̂×∇×N_j)` per Engquist & Majda 1979
   *IEEE T-AP* 27(5) p. 661 eq. 9. Reflection floor for normal
   incidence drops from `~ −40 dB` (1st-order) to `~ −60 dB`
   (2nd-order). `OpenBoundarySolver::with_abc_order(AbcOrder)` selects;
   `AbcOrder::First` default reproduces v2.
3. **F5+F6 — multi-port `S_{p,q}` matrix.** New `sweep_matrix` entry
   point + `SParametersMatrix` output. Per swept ω: assemble + LU
   once, back-substitute against `n_ports` RHS vectors (one per
   excited port), extract every `S_{q,p}(ω)` via Pozar §3.3 / Sheen
   1990 modal projection on every port face.

Six load-bearing decisions:

1. **Default-off knobs preserve v2 + CCCCCCCCC bit-for-bit.**
   `with_coupled_whitney(false)` and `with_abc_order(First)` are the
   defaults. Every v2 caller compiles and produces identical output;
   the change is additive.
2. **3-point Gauss quadrature for F1.** Degree-2-exact on the reference
   triangle. Adequate for `N_i · E_t` (degree 1 × piecewise-linear
   modal profile). 6-point fallback queued for Phase 4.fem.eig.3.0.1
   if the F1 unit-test convergence is marginal against the analytic
   TE_{10}-on-WR-90 integral.
3. **2nd-order Mur, not CFS-PML, for F3.** CFS-PML requires re-doing
   the Whitney-1 basis on stretched coordinates and adding absorbing
   tetrahedra; 2nd-order Mur adds one Gram-matrix term per ABC face
   and stays in the existing assembly path. Deferral of CFS-PML to
   Phase 4.fem.eig.3.5 mirrors ADR-0040's deferral of PML to
   Phase 4.fem.eig.2.5.
4. **`SParametersMatrix` is a new output type alongside `SParameters`.**
   The v2 single-port `SParameters` (with `s_pp: Vec<Vec<Complex64>>`)
   stays for v2 callers; v3 adds `SParametersMatrix` with
   `s: Vec<DMatrix<Complex64>>` per-frequency. Both types live in
   `yee_fem::open_boundary`; both are re-exported from `yee_fem`.
5. **LU-factor reuse across excited ports.** The driven matrix
   `A(ω) = K(ω) − k₀² M(ω) + B_ABC + B_port` is independent of which
   port is driven (every port face contributes its stiffness block;
   only the RHS depends on `a_inc_p`). Per-frequency runtime is
   `O(LU(N) + n_ports · BS(N))`, not `O(n_ports · LU(N))`. Verified
   by the F5 timing assertion.
6. **fem-eig-003 strict gates are un-ignored in F6.** Both the
   `[-45, -35] dB` absorption-floor gate and the strict
   `|S_{11}| < 1` continuum bound become CI-default. fem-eig-004
   (2-port WR-90 thru-line at 10 GHz, `|S_{21}| ≈ 1`) and fem-eig-005
   (3-port WR-90 T-junction at 5 GHz, magnitude conservation +
   reciprocity within `1e-3`) join the validation rollup.

CPU-only, single-threaded, FP64 complex. No GPU. Single dominant mode
per port. Scalar isotropic real `ε_r`, `μ_r` on the driven sweep
(combining v3 with v1's dispersive Newton tracker is Phase
4.fem.eig.3.1). `faer::sparse::FaerLuSolver<Complex64>` continues to
handle the complex-symmetric matrix unchanged.

## Consequences

- **fem-eig-003 strict gates clear without weakening tolerances.** The
  un-ignore is a single attribute removal in F6; the underlying physics
  fix is F1+F2 (coupled Whitney) + F3+F4 (2nd-order Mur).
- **The v2 single-port `OpenBoundarySolver` API stays compile-stable.**
  New methods are builder-style (`with_coupled_whitney`,
  `with_abc_order`); the existing `new`, `solve_at_frequency`, `sweep`,
  `extract_s11` are unchanged.
- **Multi-port S-parameter analysis becomes a first-class FEM
  capability.** Iris-coupled cavity filters, T-junctions, branch-line
  couplers, and Wilkinson-style power dividers all land on the same
  `sweep_matrix` surface. Combining with the v1 dispersive Newton
  tracker (Phase 4.fem.eig.3.1) unlocks lossy filter validation.
- **The 2nd-order ABC's band-edge conditioning becomes a known risk.**
  Near closed-stub TE_{10n} resonances the `−(1/2k₀)·R_2` curl
  correction can over-amplify if the dominant mode goes evanescent.
  Mitigation: per spec §10, the v3 path falls back to `AbcOrder::First`
  on a per-frequency basis if `|β_mode − k₀| / k₀ > 0.5`. The
  WR-90 stub at 8-12 GHz with `f_c = 6.56 GHz` is well above cutoff
  and should not trigger this; the canary is the fem-eig-003 sweep
  smoothness gate (C) inherited from BBBBBBBBB.
- **Multi-port modal-overlap matrix ill-conditioning becomes a new risk
  vector.** When two port modal profiles have non-trivial inner product
  in the FEM-projected basis (geometrically disjoint ports can still
  acquire numerical cross-coupling at the discretisation level), the
  `M^{-1}` correction in extraction is the load-bearing piece. The
  fem-eig-005 reciprocity tolerance `1e-3` (vs fem-eig-004's `1e-6`)
  reflects this risk. `cond(M) > 1e6` is a runtime warning.
- **CFS-PML stays available** as Phase 4.fem.eig.3.5 if 2nd-order Mur
  cannot meet some future tolerance. No code freeze on that path.
- **CCCCCCCCC's `M_pp` normalisation in `extract_s11` is preserved.**
  Under `with_coupled_whitney(true)`, `M_pp` is recomputed via the same
  3-point Gauss quadrature as the modal-projection numerator, keeping
  the round-trip identity at the exact-basis level.

## References

- `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
- `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-3.md`
- ADR-0040 — Phase 4.fem.eig.2 scope (this ADR's parent; §C-3 deferral
  this ADR fulfils).
- ADR-0039 — Phase 4.fem.eig.1 dispersive scope (grandparent).
- ADR-0029 / ADR-0032 — Phase 4.fem.eig.0 scope + plan
  (great-grandparents).
- B. Engquist and A. Majda, "Radiation boundary conditions for acoustic
  and elastic wave calculations", *Comm. Pure Appl. Math.* 32 (1979),
  pp. 313-357 — the 2nd-order ABC derivation (the IEEE T-AP 27(5) p. 661
  variant is the waveguide-mode restatement used in §4.2 of the spec,
  DOI 10.1109/TAP.1979.1142175).
- D. M. Sheen, S. M. Ali, M. D. Abouzahra, P. B. L. Katehi,
  "Application of the three-dimensional finite-difference time-domain
  method to the analysis of planar microstrip circuits",
  *IEEE Trans. Microwave Theory Tech.* 38(7) (1990), pp. 849-857,
  DOI 10.1109/22.55781 — multi-port S-parameter extraction via per-port
  driven solves with shared system matrix.
- J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley 2014, Ch. 10 (driven FEM analysis; §10.4 ABC, §10.5 wave-port,
  §10.7 S-parameters, Table 10.1 reflection floors for 1st- and
  2nd-order ABCs).
- D. M. Pozar, *Microwave Engineering*, 4th ed., 2012, §3.3 (waveguide
  modes), §4.3 (reciprocity for multi-port networks).
- A. Bossavit, "Whitney forms: a class of finite elements for
  three-dimensional computations in electromagnetism", *IEE Proc.*
  135-A (1988), pp. 493-500 — Whitney-1 basis identity.
- J.-P. Berenger, *J. Comput. Phys.* 114 (1994) — PML reference,
  Phase 4.fem.eig.3.5 reserved.
- CLAUDE.md §3, §4.
