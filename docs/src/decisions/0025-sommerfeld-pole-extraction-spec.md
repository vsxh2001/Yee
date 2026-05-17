# ADR-0025: Phase 1.1.1.2 Sommerfeld surface-wave pole extraction — spec, deferred implementation

**Status:** Accepted (spec/plan only; implementation deferred to follow-up tracks)
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0020 shipped multi-image DCIM via GPOF (Phase 1.1.1.0) and
ADR-0024 shipped the edge-clustered strip mesh (Phase 1.1.1.1).
Between the two, `mom-002`'s |Z_in| moved from ~14 kΩ to ~2.1 kΩ
and `Re(Z)` converged to ~-50 Ω. The remaining gap to the
Hammerstad–Jensen `[35, 75] Ω` target is on the **imaginary axis**:
`Im(Z)` plateaus at ~-2.1 kΩ across `nz ∈ {16, 24, 32}` (three
levels of mesh refinement) and refining further would not move it.

The plateau is a **surface-wave-pole** signature, not a mesh-bound
error. On a grounded dielectric slab (FR-4 at 1 GHz, `ε_r = 4.4`,
`h = 1.6 mm`), the slab supports discrete TM_0 and TE_1 modes that
trap energy in the substrate as a guided wave; their dispersion
relations `D_TM(k_ρ) = 0` and `D_TE(k_ρ) = 0` produce **poles in
the spectral reflection coefficient `R(k_ρ)`** that GPOF cannot
fit. GPOF's ansatz is a sum of smooth complex exponentials in the
image domain; near-singular spectral content (a pole sitting just
off the integration contour) does not have a smooth-exponential
representation. The fit converges to "the best smooth
approximation to a near-singular function," which is by definition
wrong by the residue magnitude.

The remedy is well-established in the EM literature (Aksun 1996;
Michalski–Mosig TAP 1997): **subtract the pole contribution
analytically before fitting the smooth remainder with GPOF**. The
pipeline is:

1. **Find the poles.** Newton–Raphson root-finding on `D_TE(k_ρ)
   = 0` and `D_TM(k_ρ) = 0` in the complex `k_ρ` plane, starting
   from analytic seed values (the lossless-slab cutoff
   approximations). A Müller's-method fallback handles the case
   where Newton wanders into a neighbouring basin or onto the
   wrong Riemann sheet.
2. **Extract the residues.** Closed form from
   `r_n = N(k_ρ,n) / D'(k_ρ,n)` (Pozar §3.7; Felsen–Marcuvitz §5).
3. **Subtract the pole contribution.** The spatial-domain
   pole contribution is a Hankel `H_0^(2)(k_ρ,n · ρ)` per pole
   times the residue, evaluated analytically. The remaining smooth
   spectral kernel is what GPOF fits.
4. **Add back analytically.** The full Green's function evaluation
   sums the GPOF image series and the analytic Hankel
   reconstruction.

The decision for **this** ADR is to spec the pipeline and the
public API surface; implementation is staged into the Phase
1.1.1.2 follow-up plan. The spec lives now; the matrix code
lands as a separate sub-project.

The seam on `GreensSpec` (ADR-0015) is one new variant:

```rust
pub enum GreensSpec {
    FreeSpace,
    Microstrip { eps_r: f64, h_m: f64 },
    MicrostripDcim { eps_r: f64, h_m: f64, n_images: usize },
    MicrostripSommerfeld {
        eps_r: f64,
        h_m: f64,
        n_images: usize,
        n_surface_wave_poles: usize,  // default 2 (TM_0 + TE_1)
    },
}
```

`n_surface_wave_poles = 0` is a **bit-for-bit back-compat path** to
the Phase 1.1.1.0 DCIM result — the tripwire test that ADR-0020's
`n_images = 1` test established generalises here.

## Decision

Track DDDDD ships **spec + plan** for Phase 1.1.1.2; implementation
is deferred to a follow-up track. The deliverables are:

- `docs/superpowers/specs/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction-design.md`
  — design spec.
- `docs/superpowers/plans/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction.md`
  — task-by-task plan.

The plan has six steps:

- **Step 0** — re-derive `D_TE(k_ρ)` and `D_TM(k_ρ)` against
  Pozar §3.7 (lossless-grounded-slab dispersion relations); pin
  the analytic expressions in a comment block.
- **Step 1** — Newton–Raphson root finder in the complex `k_ρ`
  plane with closed-form derivative; FR-4 sanity table at
  `f ∈ {1, 2.4, 5} GHz`.
- **Step 2** — residue extraction `r_n = N(k_ρ,n) / D'(k_ρ,n)`;
  `DegeneratePole` escape hatch for the case where `D'` vanishes
  (a degenerate slab thickness right at a mode cutoff).
- **Step 3** — pole-subtracted GPOF: subtract the analytic
  contribution from the sampled kernel, fit the residual with the
  ADR-0020 GPOF pipeline; add the Hankel `H_0^(2)` evaluator
  (small-argument series + large-argument asymptotic).
- **Step 4** — `GreensSpec::MicrostripSommerfeld` variant and
  the new `MultilayerGreens` constructor with the
  `n_surface_wave_poles = 0` bit-for-bit tripwire against
  ADR-0020.
- **Step 5** — swap `mom-002` to the new variant; tighten the
  validation gate from the loose `[1, 100 kΩ]` non-degeneracy
  band to the Hammerstad–Jensen `[35, 75] Ω` target.

The default `n_surface_wave_poles = 2` covers the dominant TM_0
and TE_1 modes for grounded-slab geometries up through low GHz.
Higher modes (TM_1, TE_2, ...) only appear above their respective
cutoff frequencies (above ~10 GHz for FR-4 at `h = 1.6 mm`); the
spec is explicit that the default may need to grow for higher-
frequency cases.

mom-001 (free-space) stays untouched: the new variant is opt-in,
and the `FreeSpace` and `Microstrip` paths are unchanged.

## Alternatives considered

1. **Brute-force GPOF with `n_images` large enough to swallow
   the pole.** Rejected. The mathematics is wrong: a pole is not a
   sum of complex exponentials, and increasing `n_images` does not
   converge to one. Empirically (Track OOOO's |Z_in| sweep at
   `n_images ∈ {1..10}`) the floor at ~kΩ does not budge.
2. **Numerical contour deformation alone.** Possible in principle:
   choose the Sommerfeld integration contour to pass between the
   poles, fit the smooth remainder. Rejected as fragile: the
   contour choice depends on pole locations that are themselves
   the thing being extracted. The Aksun 1996 deformed contour
   (used by ADR-0020) is the right move for the *image* fit but
   does not subtract pole residues from the spatial-domain Green's
   function.
3. **Tabulated Sommerfeld integral evaluation.** Some commercial
   tools precompute `∫ R(k_ρ) J_0(k_ρ ρ) k_ρ dk_ρ` on a dense
   grid and interpolate at query time. Rejected as the wrong
   shape for `yee-mom`: the in-tree solver is impedance-matrix
   based and needs Green's function evaluations at thousands of
   `(ρ, z)` query points per frequency; tabulation cost dominates
   real workloads.
4. **Vector fitting in the frequency domain instead of the
   spectral domain.** Rejected on the same grounds as in ADR-0020:
   rational-function fitting is a conceptual mismatch for the
   complex-exponential ansatz that DCIM requires.

## Consequences

**What becomes easier:**

- **`mom-002`'s `[35, 75] Ω` Hammerstad–Jensen gate is reachable
  in principle.** The spec lays out the algorithm; the plan
  splits it into verifiable steps; the test infrastructure
  (sweep table, headline gate, edge-clustered mesh) is already in
  place from ADR-0024.
- **`mom-003` (2.4 GHz patch), `mom-004` (Wilkinson), `mom-005`
  (branch-line hybrid), `mom-006` (Swanson hairpin BPF) all
  depend on the same Sommerfeld machinery.** Phase 1.1.1.2 is on
  the critical path for tightening all of them.
- **The pole-subtraction abstraction is documented now.** The
  `n_surface_wave_poles` parameter will be visible in
  `GreensSpec::MicrostripSommerfeld` when the implementation
  lands; consumers can write code against the planned variant.

**What becomes harder:**

- **`mom-002`'s tolerance stays loose until 1.1.1.2 implements.**
  CLAUDE.md §4 / §10 still flag `mom-002` and `mom-003` as
  loose-tolerance. The spec narrows the rationale ("waiting on
  surface-wave pole subtraction") but does not change the gate.
- **The Newton root-finder is the riskiest step.** Newton on
  `D(k_ρ) = 0` in the complex plane can wander; the plan's
  Müller's-method fallback exists to handle that, but until
  implementation lands we are taking the literature's word.
- **The default `n_surface_wave_poles = 2` is FR-4-1GHz-specific.**
  Other substrate / frequency combinations may need a different
  default. The spec is explicit; the documentation will need to
  carry it forward.

**What's now closed off:**

- A vector-fitting or rational-function-fitting approach to the
  pole problem. The spec mandates Newton + analytic residue +
  pole-subtracted GPOF + analytic Hankel reconstruction.
- A `MicrostripSommerfeld` variant with a different parameter
  list. The four fields (`eps_r`, `h_m`, `n_images`,
  `n_surface_wave_poles`) are the spec'd surface; changing them
  requires another ADR.
- Quietly tightening `mom-002` / `mom-003` tolerances before
  1.1.1.2 lands. The forward reference in ADR-0024's docstring
  and in this ADR makes the dependency explicit.

## References

- `docs/superpowers/specs/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction-design.md`
  — design spec.
- `docs/superpowers/plans/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction.md`
  — implementation plan.
- Commits 6738d8c (spec), 6fee4f7 (plan), ca4241c (Track DDDDD
  merge).
- D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley, 2011,
  §3.7 (grounded dielectric slab surface-wave dispersion
  relations).
- M. I. Aksun, "A robust approach for the derivation of
  closed-form Green's functions," *IEEE TMTT*, vol. 44, no. 5,
  pp. 651–658, May 1996 (the contour and the multi-image
  framework that pole subtraction sits on top of).
- K. A. Michalski and J. R. Mosig, "Multilayered media Green's
  functions in integral equation formulations," *IEEE TAP*,
  vol. 45, no. 3, pp. 508–519, March 1997 (the
  Michalski–Mosig pole-and-image decomposition).
- L. B. Felsen and N. Marcuvitz, *Radiation and Scattering of
  Waves*, IEEE Press classic reissue, 1994, §5 (canonical
  surface-wave residue source).
- ADR-0015 — `GreensSpec` builder.
- ADR-0020 — Phase 1.1.1.0 multi-image DCIM via GPOF; the
  fit that this ADR's pole subtraction feeds into.
- ADR-0024 — Phase 1.1.1.1 edge-clustered mesh; the `Re(Z)`
  convergence that establishes the surface-wave-pole diagnosis.
- CLAUDE.md §4 — `mom-002` / `mom-003` loose tolerances; §10 —
  `MultilayerGreens` placeholder caveat (still open on the pole
  leg until this ADR's implementation lands).
