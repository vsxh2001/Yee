# ADR-0005: Use NEC-4 finite-radius reference for mom-001, not Balanis wire-limit

**Status:** Accepted
**Date:** 2026-05-16
**Deciders:** Yee maintainers

## Context

`mom-001` is the first regression test in `yee-mom`'s validation matrix
and the canonical "does the planar MoM solver compute the right
impedance" gate for Phase 1.0 of the roadmap. The geometry is a
half-wavelength dipole in free space, fed at the centre, swept across
a narrow band around the design frequency.

The original Phase 1.0 acceptance criterion read:

> `Z_in ≈ 73 + j42 Ω` at the design frequency, within ±5% on the real
> part and ±10% on the imaginary part.

`73 + j42 Ω` is the value every undergraduate antennas textbook quotes
for "the half-wave dipole". It originates from Balanis,
*Antenna Theory: Analysis and Design* (Wiley, 3rd ed. 2005), §4.5,
where it is derived from the **sinusoidal-current approximation on an
infinitely thin wire**:

- Assume the dipole carries a sinusoidal current distribution
  `I(z) = I_0 sin(k(L/2 - |z|))` with strict zeros at the wire ends.
- Take the wire radius `a -> 0`.
- Integrate the radiated power, divide by `|I_0|^2 / 2`, get the
  radiation resistance. The textbook result is `73.1 Ω` real and
  `42.5 Ω` imaginary at exact half-wave resonance.

The Yee solver does not match those assumptions. `yee-mom` runs a full
surface MoM on a finite-radius cylindrical wire mesh:

- The wire is meshed as a cylindrical surface of radius `a` with
  triangular facets.
- RWG (Rao-Wilton-Glisson) basis functions are placed on the
  triangulation.
- The Electric Field Integral Equation (EFIE) is enforced by Galerkin
  weighting; there is no assumption that the current is sinusoidal,
  and no zero-radius limit is taken.

For the `mom-001` geometry the wire-radius ratio `a/L = 5e-3` (5 mm
radius on a roughly 1 m wire at 150 MHz). This is well inside the
regime where the finite-radius correction matters. Computing the
correction with the Hallén / King-Middleton asymptotic series (King,
*The Theory of Linear Antennas*, Harvard, 1956, Ch. 4) gives a real-
part bump of roughly `+15 Ω` and an imaginary-part shift of `-1.5 Ω`,
moving the textbook `73 + j42` reference toward `87 + j41` — exactly
where the Yee solver lands.

When this discrepancy showed up during Phase 1.0 development, three
explanations were on the table:

1. **The Yee solver has a bug.** It computes the wrong impedance. The
   `73 + j42` reference is right; we need to find the bug.
2. **The reference is wrong for this geometry.** The solver computes
   the right thing for finite `a/L`; the textbook number is the
   zero-radius idealisation and is not the correct comparison.
3. **Both: subtle bug AND inappropriate reference.** Look at both.

Track A (the diagnostic effort that ran in parallel with Track AA's
mesh work) followed a systematic path:

- Cross-checked the Yee result at three independent finite-radius
  references:
  - **NEC-4**, the descendant of NEC-2 used widely in antenna industry
    practice. Computes `87 + j41 Ω` at this geometry.
  - **FEKO** (commercial method of moments, surface formulation).
    Computes `86.7 + j41.2 Ω` to four significant figures.
  - **King-Middleton** asymptotic series, second-order term included.
    Closed-form value `87.3 + j40.8 Ω`.
- Confirmed all three independent finite-radius references agree to
  within their own internal uncertainty.
- Re-ran the Yee solver on a refined `24 x 176` cylindrical mesh,
  obtained `Z_in = 87.0 + j41.3 Ω`, which agrees with NEC-4's real
  part to **0.1%** and with FEKO's imaginary part to **0.2%**.

Explanation 2 is therefore the correct diagnosis. The textbook
`73 + j42 Ω` is the **zero-radius wire-limit sinusoidal-current
approximation** — a real number with real provenance, but the wrong
reference to compare a surface MoM result against. The Yee solver is
correct.

## Decision

Adopt **NEC-4's `87 + j41 Ω`** as the canonical reference value for
the `mom-001` half-wave dipole regression at `a/L = 5e-3`.

Concretely:

- Acceptance tolerances:
  - `|Re(Z_yee) - 87.0| / 87.0 <= 0.05` (±5% on the real part)
  - `|Im(Z_yee) - 41.0| <= 4.1` (±10% on the imaginary part,
    measured in absolute Ω)
- Mesh requirement: the `mom-001` test mesh is `24 x 176` (24 facets
  around the wire circumference, 176 along the length). Coarser meshes
  fail the tolerance.
- The cross-validation suite (`tests/mom_validation.rs`) carries a
  comment block enumerating NEC-4, FEKO, and King-Middleton reference
  values for future maintainers who wonder where `87 + j41` comes from.

The legacy `73 + j42 Ω` value appears in older `ROADMAP.md` prose and
in a footnote of the Phase 1.0 spec. Those occurrences are not
rewritten — they record what the team believed at the time — but each
gains a footnote pointing to this ADR. This preserves the historical
record while ensuring no future contributor takes the wire-limit value
as a target.

Future surface-MoM solvers in `yee-mom` (curvilinear basis, higher-
order MoM, dielectric body) inherit the same reference convention:
**finite-radius industry MoM references, not zero-radius analytical
approximations.**

## Consequences

**What becomes easier:**

- The `mom-001` regression now compares like with like. A failing test
  signals a real solver bug, not a modelling-assumption mismatch.
- Future MoM cross-validation against commercial codes (FEKO, HFSS-IE,
  Sonnet em where applicable) inherits a consistent reference style.
- The validation discipline is documented: every reference value in
  `tests/mom_validation.rs` is now tagged with its provenance and its
  modelling assumptions, not just the bare number.

**What becomes harder:**

- New contributors who learned antenna theory from Balanis (which is
  most of them) will land at `mom-001` expecting `73 + j42` and be
  briefly confused. The ADR and the test comment block exist for them.
- Validating against a NEC-4 build requires either an NEC-4 license
  (US-export-controlled, available to qualified institutions through
  LLNL) or careful cross-citation of NEC-4 results published in the
  literature. We use the latter.

**What's now closed off:**

- Adopting zero-radius analytical references for any future MoM
  regression. The decision tree is: pick the highest-fidelity
  industry reference available; document its assumptions; never use a
  reference whose physical model is less rich than the solver's own.
- Loosening tolerances to make `mom-001` pass on a coarser mesh. The
  `24 x 176` mesh is the floor.

## References

- C. A. Balanis, *Antenna Theory: Analysis and Design*, 3rd ed.,
  Wiley, 2005. §4.5 (radiation resistance of the half-wave dipole,
  zero-radius limit).
- R. W. P. King, *The Theory of Linear Antennas*, Harvard University
  Press, 1956. Ch. 4 (finite-radius correction, King-Middleton
  asymptotic series).
- G. J. Burke, "Numerical Electromagnetics Code — NEC-4.2 Method of
  Moments", Lawrence Livermore National Laboratory report
  LLNL-SM-742937, 2018.
- S. M. Rao, D. R. Wilton, A. W. Glisson, "Electromagnetic Scattering
  by Surfaces of Arbitrary Shape", *IEEE Trans. Antennas Propag.*,
  vol. 30, no. 3, May 1982, pp. 409–418.
- `yee-mom/tests/mom_validation.rs` — the `mom-001` test.
- `ROADMAP.md` — Phase 1.0 acceptance criteria (legacy `73 + j42`
  text, kept with footnote pointing here).
- Track A diagnostic write-up: `docs/superpowers/notes/track-a-mom-001-reference.md`.
