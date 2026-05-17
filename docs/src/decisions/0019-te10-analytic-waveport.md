# ADR-0019: Phase 1.3.1.0 ships an analytic TE10 wave-port, defers arbitrary-cross-section eigensolver to 1.3.1.1

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

Phase 1.3.0 introduced `WavePort` as an API placeholder. Its
`rhs` path matched `DeltaGapPort` bit-for-bit (CLAUDE.md §10 calls
this out explicitly): the modal-distribution machinery existed on the
public surface, but the port itself behaved like a delta-gap
excitation for every cross-section. This was the walking-skeleton
shape — ship the type, defer the modal solver — and it was
sufficient to unblock callers that wanted to *say* "wave-port" in
their geometry definitions.

The follow-up question for Phase 1.3.1 was: **what is the smallest
extension that produces non-trivial modal behaviour?** The textbook
answer for an arbitrary cross-section is a numerical 2-D eigensolver
(Nedelec edge elements + nodal Lagrange for `E_z`, sparse generalised
eigenproblem, shift-and-invert Arnoldi). That is the right end state.
It is also a substantial sub-project: it pulls in `nalgebra-sparse`,
optionally `arpack-rs`, a fresh `TriMesh2D` type in `yee-mesh`, and a
slow-wave validation gate against a loaded-WR-90 reference. Phase
1.3.1 splits along an obvious seam:

- **Phase 1.3.1.0 — analytic rectangular-waveguide TE10.** Closed-form
  Pozar §3.3 quantities: `f_c = c / (2 a √ε_r)`, `β_10 = √(k² − (π/a)²)`,
  `Z_TE10 = η / √(1 − (f_c/f)²)`, `E_y(x) = sin(π x / a)`. No matrix
  assembly, no eigensolve, no new mesh type. Sufficient to validate a
  wave-port-driven simulation on a rectangular waveguide cross-section
  (WR-90 at X-band is the canonical case) and to exercise the
  `ModalDistribution` dispatch surface that the numerical path will
  eventually plug into.
- **Phase 1.3.1.1 — numerical 2-D eigensolver.** Arbitrary
  cross-sections (microstrip, CPW, loaded waveguide). Its own spec,
  its own plan, its own validation gate. The TE10 closed-form is the
  cross-check reference: the numerical eigensolver must reproduce
  TE10 within 0.1% on β and 1% on `Z_w`.

The choice for **this** ADR is whether 1.3.1 ships as a single
all-or-nothing sub-project or as two staged sub-projects. The
walking-skeleton-first principle (CLAUDE.md §3) and the
"each Phase-X.0 placeholder beats a half-finished Phase-X.1" rule
settle it in favour of staging: Phase 1.3.1.0 ships the analytic
TE10, Phase 1.3.1.1 ships the numerical eigensolver behind the same
`ModalDistribution` dispatch.

## Decision

`yee-mom` Phase 1.3.1.0 adds a closed-form rectangular-waveguide
TE10 mode in `RectangularWaveguideTe10`, surfaced through a new
`ModalDistribution::Te10` variant and a `WavePort::with_rectangular_te10`
builder. The numerical 2-D eigensolver for arbitrary cross-sections
is deferred to Phase 1.3.1.1 (separately specced in `docs/superpowers/specs/`,
see ADR-0022).

**Public API:**

```rust
pub struct RectangularWaveguideTe10 {
    pub a: f64,      // broad dimension (m)
    pub b: f64,      // narrow dimension (m)
    pub eps_r: f64,  // relative permittivity of fill
}

impl RectangularWaveguideTe10 {
    pub fn cutoff_hz(&self) -> f64;          // c / (2 a √ε_r)
    pub fn beta(&self, f: f64) -> f64;       // √(k² − (π/a)²), NaN at/below cutoff
    pub fn wave_impedance(&self, f: f64) -> f64; // η / √(1 − (f_c/f)²)
    pub fn e_y_profile(&self, x: f64, y: f64) -> f64; // sin(π x / a)
}

pub enum ModalDistribution {
    Uniform,              // Phase 1.3.0 delta-gap-equivalent default
    Te10(RectangularWaveguideTe10),
}
```

`WavePort::with_rectangular_te10(a, b, eps_r)` is the builder; the
default constructor still produces `ModalDistribution::Uniform`, so
existing callers (mom-001 dipole, mom-002 microstrip) see no
behaviour change. Inside `WavePort::rhs`, the TE10 variant samples
`e_y_profile` at each port-edge midpoint and weights by edge length;
the `Uniform` variant retains the Phase 1.3.0 bit-for-bit equivalence
to `DeltaGapPort`.

The `ports` module is promoted to `pub` so integration tests and
downstream crates can construct ports directly; the `Port` trait
itself remains `pub(crate)` until the numerical eigensolver lands and
the full modal-distribution dispatch surface is finalised.

## Alternatives considered

1. **Ship Phase 1.3.1.1 directly (numerical eigensolver in one step).**
   Rejected. The eigensolver is a roughly 2000-LOC sub-project once
   meshing, assembly, sparse linear algebra, and the validation gate
   are counted. Shipping it without an analytic cross-check is a
   correctness risk: the *only* way to know the eigensolver is
   correct on a rectangular waveguide is to compare against the
   closed-form TE10. Bootstrapping the analytic mode first gives the
   eigensolver work a reference to bisect against.
2. **Keep delta-gap-only (skip 1.3.1.0 entirely).** Rejected. The
   placeholder behaviour of Phase 1.3.0 is documented in CLAUDE.md §10
   as a known limitation; leaving it in place blocks waveguide
   validation cases (and, by extension, the upcoming mom-004 / mom-005
   / mom-006 spec cases that have wave-port-driven variants). The
   analytic TE10 is the smallest concrete step that closes the
   "wave-port is bit-for-bit delta-gap" gap.
3. **Closed-form library other than Pozar (Collin / Harrington /
   Balanis).** Rejected. Pozar §3.3 is the most-cited reference in
   microwave engineering practice and matches what the rest of
   `yee-mom`'s documentation already cites. Using a different
   convention would force a per-formula cross-walk in the docs for no
   technical gain.

## Consequences

**What becomes easier:**

- **WR-90 (and any rectangular-waveguide cross-section) is a
  first-class wave-port geometry**, with a closed-form mode against
  which the numerical eigensolver can be regression-tested when it
  lands. The validation case (a=22.86mm, b=10.16mm, X-band) becomes
  the bisection point for Phase 1.3.1.1.
- **The `ModalDistribution` dispatch surface is exercised
  end-to-end.** Phase 1.3.1.1 plugs in a `Numerical2D(...)` variant
  without touching `WavePort::rhs` or any caller.
- **mom-001 (the NEC-4 dipole gate) is unaffected.** The default
  `Uniform` variant preserves Phase 1.3.0 behaviour bit-for-bit; the
  dipole still uses `DeltaGapPort`, not `WavePort`.

**What becomes harder:**

- **Microstrip, CPW, and arbitrary cross-sections still produce
  delta-gap-equivalent results** until Phase 1.3.1.1. CLAUDE.md §10's
  WavePort caveat narrows but does not disappear: it now reads
  "non-rectangular cross-sections still match delta-gap" instead of
  "all cross-sections match delta-gap."
- **`pub` promotion of the `ports` module increases the public
  surface area** before the modal-distribution dispatch is finalised.
  Phase 1.3.1.1 will likely add fields / variants; semver discipline
  in the run-up to 1.0 will need to budget for that.

**What's now closed off:**

- A wave-port API in which the modal distribution is implicit (e.g.
  hard-coded TE10). The `ModalDistribution` enum is the seam; new
  modes add variants there.
- Bypassing the `Uniform` default. Existing callers (notably mom-001
  and mom-002) cannot accidentally pick up TE10 behaviour by omission;
  it is opt-in via `with_rectangular_te10`.

## References

- `crates/yee-mom/src/ports.rs` — `RectangularWaveguideTe10`,
  `ModalDistribution::Te10`, `WavePort::with_rectangular_te10`.
- `crates/yee-mom/tests/te10_waveport.rs` — WR-90 validation (cutoff
  6.5575 GHz, β finite above cutoff and NaN below, profile zero at
  walls / peak at centre, `Z_TE10 > η_0`).
- Commits 809ae64 (Phase 1.3.1.0 implementation) and 47e6b51 (Track
  IIII merge).
- D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley, 2011, §3.3
  (rectangular-waveguide TE/TM modes).
- ADR-0022 — Phase 1.3.1.1 numerical-eigensolver spec; consumes the
  TE10 closed form as its bisection reference.
- CLAUDE.md §3 — walking-skeleton-first; §10 — `WavePort`
  placeholder caveat (narrowed by this ADR).
