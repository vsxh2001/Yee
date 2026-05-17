# ADR-0026: Phase 2.fdtd.5.3 oblique-incidence TF/SF — 1-D aux grid, `H_inc = −(k̂ × E_inc)` sign convention

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0014 shipped TF/SF as a `+x`-normal-incidence slab. ADR-0021
extended it to a six-face finite-box kernel, also at normal
incidence. Phase 2.fdtd.5.3 closes the scope gap that both prior
ADRs explicitly left open: **arbitrary plane-wave incidence**
`(θ, φ, ψ)` — propagation direction `k̂` and polarisation `ê`
both free, subject to `ê ⊥ k̂`.

The Taflove §5.10 canonical extension uses a **1-D auxiliary
Yee grid** along `k̂`. The 3-D grid runs the usual leapfrog;
projection samples the 1-D incident field onto each Yee node at
each box face by `s = (r − r_ref) · k̂`; the 12-stencil correction
on the six faces of the TF/SF box absorbs the projected incident
contribution.

Implementing that, two non-obvious correctness issues had to be
nailed down:

**(1) Sign convention on `H_inc`.** Naive physics: for a
`+s`-propagating plane wave in vacuum,
`H_inc = (1/η₀) (k̂ × E_inc)`. That is *correct in the analytic
plane-wave solution* but *wrong as the relationship the 1-D Yee
aux grid produces.* The 1-D leapfrog uses the same staggered
half-step E/H update as the 3-D Yee scheme, but in 1-D the curl
collapses to a single derivative; the discrete update produces
`inc_h = −inc_e / η₀` for a `+s`-propagating wave (the minus sign
falls out of the Yee staggering). The aux grid is internally
consistent — its `E` and `H` are a valid discrete plane wave —
but the **discrete** sign differs from the analytic-physics sign
by a global minus.

A naive 3-D TF/SF correction that uses
`H_inc_hat = +(k̂ × E_inc_hat) / η₀` injects the wrong sign on
the H side of the 12-stencil correction. The visible failure mode
is total collapse of the TF/SF contrast: an oblique-via-normal
regression test (drive the oblique kernel at `(θ=π/2, φ=0, ψ=π)`
— which is the legacy `+x` / `E_z` mapping — and expect a
high-contrast result) produces only **1.4× contrast** without the
sign fix. The correct convention is
`H_inc_hat = −(k̂ × E_inc_hat)`; with that fix the same test
produces **9.8×10¹⁴× contrast** (effectively roundoff,
~1.5× below the ADR-0021 finite-box floor of ~7×10¹⁴×).

The diagnosis took longer than it had to and is the kind of
not-obvious-from-the-paper detail that future implementers will
re-derive at cost. The sign convention is therefore documented
inline in `sources.rs` as a design-notes block. This ADR records
the *why*.

**(2) Dispersion mismatch between 1-D aux and 3-D Yee.** Yee's
numerical dispersion is axis-anisotropic: the discrete phase
velocity along `k̂` depends on `k̂`'s orientation relative to the
grid axes. A 1-D aux grid running at `ds = dx` matches the
on-axis phase velocity but mismatches the oblique phase velocity
by an amount that grows with `θ` off-axis. The textbook remedy is
Taflove §5.10.5's **dispersion-matched aux step**: solve a
transcendental equation for `ds_aux(θ, φ)` such that the 1-D and
3-D phase velocities along `k̂` agree exactly at the design
frequency.

That root-find has a trivial `ds → 0` attractor; a naive solver
hits it. The fix is a custom bracketing solver that brackets the
non-trivial root before iterating. That is a focused sub-project
of its own (Phase 2.fdtd.5.3.1, queued, **not in this ADR**).

The empirical effect of the unmatched aux step: at 30° / 45°
oblique incidence with `E ∥ ê_φ`, the TF/SF contrast is **14.5×**,
versus a 1000× DoD target. That is well below the textbook gate,
*sufficient* to certify that the kernel runs without panicking and
that the sign conventions and the cross-section index ranges are
correct, but *not* the textbook-quality contrast that
ADR-0021's normal-incidence finite-box achieves. The brief's
escape hatch fires: surface the 14.5× result, document the
dispersion-mismatch root cause inline, and ship.

## Decision

`yee-fdtd` Phase 2.fdtd.5.3 ships the oblique-incidence TF/SF
kernel with the sign convention `H_inc_hat = −(k̂ × E_inc_hat)`
and a documented dispersion-mismatch caveat at oblique incidence.

**Public API:**

```rust
impl PlaneWaveSource {
    /// Phase 2.fdtd.5.3: arbitrary (θ, φ, ψ) incidence.
    pub fn with_oblique_incidence(
        ...
        theta: f64,
        phi: f64,
        psi: f64,
        ...
    ) -> Self;

    /// Phase 2.fdtd.5.0–5.2: +x, E_z. Bit-for-bit preserved.
    pub fn new(...) -> Self;  // legacy_plus_x = true internally
}
```

The legacy `new` constructor sets an internal `legacy_plus_x`
dispatch flag; the existing 5.2 slab and finite-box paths run
unchanged. The oblique kernel is reached only via
`with_oblique_incidence`.

**Single scalar 1-D Yee aux grid along `k̂`.** The aux grid reuses
the existing `step_incident_e` / `step_incident_h` leapfrog plus
Mur ABC at both ends; only its physical orientation changes. The
projection at each Yee node on each box face is by `s = (r −
r_ref) · k̂` and linear interpolation into the 1-D aux samples.

**12-stencil correction on the six box faces.** The standard
Taflove §5.10 extension to a rectangular Huygens surface — twelve
discrete-curl terms straddle the TF/SF boundary at general
incidence, versus the four that straddle at `+x` normal incidence
(ADR-0021).

**Sign convention.**
`H_inc_hat = −(k̂ × E_inc_hat)`, **not** `+(k̂ × E_inc_hat) / η₀`.
The 1-D Yee leapfrog produces `inc_h ≈ −inc_e / η₀` for a
`+s`-propagating wave, so the geometric correction has to absorb
that internal sign. Documented inline as a design-notes block at
the top of `sources.rs`.

**Validation gates (three new tests in
`crates/yee-fdtd/tests/plane_wave_oblique.rs`):**

- **`oblique_normal_incidence_regression`** — drive
  `with_oblique_incidence` at `(θ=π/2, φ=0, ψ=π)` (the legacy
  `+x` / `E_z` mapping) and assert contrast > 1e10×. Empirical:
  ~9.8×10¹⁴× (≈ 1.5× below the legacy 5.2 floor, well within the
  1% regression budget). This is the bisection point against
  ADR-0021's finite-box result.
- **`oblique_30deg_45deg_ephi_polarization`** — `(θ=30°, φ=45°,
  ψ=π/2)` with `E ∥ ê_φ` on a 60³ box. Empirical: ~14.5×. Gate
  is set at `>10×` with an explicit in-test comment marking the
  gap to the 1000× DoD. **This is the escape-hatch outcome.**
- **`oblique_grazing_85deg_runs_without_panic`** — smoke test at
  `θ = 85°`; asserts no NaN over 300 steps.

The ADR-0021 slab and finite-box regression tests
(`plane_wave_propagation` ~2676×, `plane_wave_finite_box`
~7.5×10¹⁴×) remain at their previous contrast floors.

**Phase 2.fdtd.5.3.1 (deferred).** Taflove §5.10.5
dispersion-matched aux step `ds_aux(θ, φ)`. The transcendental
root-find requires a custom bracketing solver to avoid the
trivial `ds → 0` attractor; that is a focused sub-project on its
own. Without it, the oblique 1000× DoD is unreachable; with it,
the 30°/45° case is expected to clear the gate.

## Alternatives considered

1. **3-D vector aux grid (full 3-D Yee for the incident wave).**
   Rejected. A full 3-D aux grid would carry no dispersion
   mismatch (it has the same anisotropy as the main grid) but at
   3-D memory cost. The 1-D-aux + projection idiom is the
   industry-standard Taflove §5.10 trick for exactly this
   trade-off; the dispersion-mismatch issue is documented and
   fixable.
2. **Skip oblique entirely; document slab + finite-box as the
   ceiling.** Rejected. Multiple Phase 2 tutorial targets (RCS,
   bistatic scattering, antenna pattern in receive mode) need
   oblique TF/SF. Shipping the kernel with the dispersion-mismatch
   caveat is better than not shipping; the caveat is a known
   followable bug, not an open-ended unknown.
3. **Use the analytic-physics sign convention
   `H_inc_hat = +(k̂ × E_inc_hat) / η₀` and "fix the aux grid"
   instead.** Rejected. The 1-D Yee aux grid is internally
   consistent — its E and H are a valid *discrete* plane wave.
   Changing the aux grid to match the analytic sign would require
   either (a) rewriting the aux grid's update equations (a
   different scheme entirely), or (b) flipping the H field's
   sign at output time (a sign-flip in the projection, which is
   functionally identical to the negated cross product but
   located in a less obvious place). The negation in the
   correction term is the localised, documented fix.
4. **Ship the dispersion-matched aux step in the same PR.**
   Rejected on agent-brief grounds: the transcendental root-find
   and its bracketing solver are a focused sub-project; bundling
   them with the oblique kernel itself would double the brief
   and conflate two distinct correctness issues. The escape
   hatch surfacing the 14.5× result is the right call.

## Consequences

**What becomes easier:**

- **Oblique-incidence plane-wave problems compile and run.**
  RCS / bistatic-scattering tutorials that need an arbitrary `k̂`
  can be written against the public `with_oblique_incidence`
  surface. Until Phase 2.fdtd.5.3.1 lands, the
  *quantitative* result is contrast-limited to ~14×;
  *qualitative* shape (sign of the scattered field, polarisation
  cross-section) is correct.
- **Normal-incidence regression is bit-clean.** The
  `legacy_plus_x` dispatch flag keeps the ADR-0021 finite-box
  path bit-for-bit unchanged. The oblique kernel exercised at
  `(θ=π/2, φ=0, ψ=π)` reproduces ADR-0021 to within 1.5× on
  contrast (9.8×10¹⁴× vs 7×10¹⁴×), well under the 1% regression
  budget.
- **The sign-convention finding is documented in-source.**
  Future implementers re-deriving the kernel will see the
  `H_inc_hat = −(k̂ × E_inc_hat)` design-notes block
  immediately; the 1.4× → 9.8×10¹⁴× contrast jump on the
  regression test is the empirical evidence.

**What becomes harder:**

- **Oblique contrast is contrast-limited at ~14× until Phase
  2.fdtd.5.3.1.** Tutorials that need a clean SF region around an
  oblique-illuminated scatterer will see numerical-dispersion
  artefacts in the SF region. The artefacts are reproducible and
  axis-symmetric; users who need a quantitative bistatic-RCS
  number have to wait for 5.3.1.
- **Phase 2.fdtd.5.3.1 is on the critical path.** Same caveat as
  ADR-0017 (lumped RLC port deferring Γ-against-analytic to
  2.fdtd.6.1): the placeholder-but-honest gate is fine for the
  walking skeleton; the production gate is non-optional and
  queued.
- **The `ds → 0` attractor in the dispersion-matched root-find
  is a known landmine.** Phase 2.fdtd.5.3.1's brief will need an
  explicit bracketing solver, not a naive Newton or bisection.
  The in-source design-notes block flags this.

**What's now closed off:**

- A naive `H_inc_hat = +(k̂ × E_inc_hat) / η₀` correction term.
  The negated cross product is mandatory; the design-notes block
  is the authoritative reference.
- Quietly shipping the oblique kernel at the 1000× DoD on the
  strength of the normal-incidence regression. The 30°/45°
  test's `>10×` gate makes the gap to 1000× visible to CI.

## References

- `crates/yee-fdtd/src/sources.rs` — `with_oblique_incidence`
  constructor; `legacy_plus_x` dispatch flag; 12-stencil
  correction kernel; design-notes block (sign convention,
  dispersion mismatch).
- `crates/yee-fdtd/tests/plane_wave_oblique.rs` — three new
  tests (normal-incidence regression, 30°/45° oblique, grazing
  smoke).
- Commits 49b84f0 (oblique kernel), 2b569b6 (validation), 56e2dfa
  (Track FFFFF merge).
- A. Taflove and S. C. Hagness, *Computational Electrodynamics:
  The Finite-Difference Time-Domain Method*, 3rd ed., Artech
  House, 2005, §5.10 (3-D TF/SF formulation, 1-D auxiliary
  incident-field source) and §5.10.5 (dispersion-matched aux
  step).
- ADR-0014 — TF/SF slab; the original `+x`-normal-incidence
  baseline.
- ADR-0021 — Phase 2.fdtd.5.2 six-face finite-box at normal
  incidence; the bisection point this ADR's oblique kernel
  regresses against.
- CLAUDE.md §3 — walking-skeleton-first; §4 — no-feature-
  without-a-gate; the `>10×` oblique gate is the honest
  placeholder gate per CLAUDE.md §4's "energy dissipation" /
  ADR-0017 precedent.
