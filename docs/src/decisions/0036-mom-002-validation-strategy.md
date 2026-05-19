# ADR-0036: mom-002 validation strategy — sub-wavelength strip needs reframe

**Status:** Accepted (2026-05-19, Track CCCCCCC investigation).

## Context

The mom-002 validation gate (FR-4 / w = 2.94 mm / h = 1.6 mm / L = 30 mm
at 1 GHz) targets the Hammerstad-Jensen analytic characteristic impedance
`Z_0 ≈ 51 Ω`. The current implementation measures `|Z_in| ≈ 2215 Ω` — a
~43× discrepancy that drove eight diagnostic / fix tracks
(EEEEEE, JJJJJJ, PPPPPP, SSSSSS, TTTTTT, XXXXXX, YYYYYY, CCCCCCC) over
the course of the May 2026 session, each landing a real correctness
improvement but none closing the headline-gate gap.

Track CCCCCCC ran a 10-configuration mesh sensitivity sweep and proved:

1. The 30 mm strip at 1 GHz is electrically `λ_eff / 10` (with
   `λ_eff = c / (f · √ε_eff)` and `ε_eff ≈ 3.32`). At this length the
   line is far short of the electrically-long regime where the input
   impedance approaches the characteristic impedance `Z_0`. The
   sub-wavelength input impedance is instead dominated by the
   open-circuit reactance of a short antenna, `|Z_in| ~ -j · η_0 /
   (k_0 L) · ln(L / a_eff) ≈ -j 2 kΩ`, which is exactly what the solver
   reports.
2. The "Re(Z) = -19 Ω" finding from Track YYYYYY's M2 audit is a
   numerical-precision artifact specific to the end-feed clustered
   Chebyshev mesh (36:1 cell aspect ratios + asymmetric stub). Other
   mesh configurations (centered port, uniform spacing) produce
   `Re(Z) ≥ 0` cleanly. The centered-uniform strip at the 5 GHz
   half-wave resonance gives `Z = +78.7 + j 3.6 Ω` — dipole-like, as
   physics predicts. The radiation resistance at 1 GHz on a 30 mm
   strip is `(k_0 L)² · η_0 / (12π) ≈ 1.3 Ω`, so the observed `-19 Ω`
   is `~0.8 % of |Z|` — well within numerical noise for a stretched
   mesh.
3. The Phase 1.0 mom-001 dipole gate (NEC-4 87 + j 41 Ω at L = 1 m, k_0
   L ≈ π — electrically a half-wave) still passes through the
   identical code paths, confirming no formulation-level bug.

## Decision

The validation gate framing was inappropriate for the chosen geometry.
Either:

1. **Lengthen the strip** to `L ≈ λ_eff / 2 ≈ 82 mm` at 1 GHz so the
   line is a half-wave resonator and `|Z_in|` can be compared to `Z_0`
   via the standard short-circuited / open-circuited line relation; or
2. **Replace the analytic target** with the short-antenna capacitive
   formula `Z_in ≈ -j · η_0 / (k_0 L) · ln(L / a_eff)` and tolerance
   it loosely (the ~10 % accuracy of the Hallen reference is the
   floor on a thin-strip Hallen-type approximation).

Both are open paths; Track GGGGGGG will pick one and implement.

## Consequences

* The kernel-side fixes already landed (EEEEEE Sommerfeld prefactor
  canonical form, TTTTTT residue sign + factor-of-2, and any pending
  DDDDDDD DCIM-TM sign correction) remain **correctness improvements
  independent of the validation reframe**. They fix the Sommerfeld
  surface-wave residue extraction for ALL downstream use cases, not
  only mom-002.
* The six `#[ignore]`-gated diagnostic tests under `crates/yee-mom/tests/`
  (sommerfeld_residue_diagnostic, sommerfeld_synthetic,
  mom_002_extent_sensitivity, mom_002_h2_gpof_diagnostic,
  mom_002_reflection_convention, mom_002_psi_port_audit,
  mom_002_mpie_audit) remain valuable forensic records of what was
  ruled out. They are kept for the historical context and will continue
  to gate on regressions in their respective probes.
* The headline mom-002 gate in `crates/yee-validation/src/lib.rs` will
  shift its measured constant when the strip length / port placement
  changes; the ±5 % tripwire band stays unchanged in spirit but will be
  reapplied around the new measured value.
* The Phase 1.1.0 / 1.1.1.x deferred-tolerance language in CLAUDE.md
  §10 and ROADMAP.md ("loose tolerances until the real multilayer
  Green's function lands") now correctly applies to the **reframed**
  mom-002, not the current 30 mm-strip case.

## References

* Track CCCCCCC investigation report (no committed SHA; surfaced via
  agent escape hatch on 2026-05-19).
* Pozar, *Microwave Engineering*, 4th ed., §2.5 (electrically-long line
  classification), §10.4 (open- and short-circuit input impedances).
* Hammerstad, E. and Jensen, Ø., "Accurate Models for Microstrip
  Computer-Aided Design," *MTT-S Digest*, 1980.
* Hallen, E., "Theoretical Investigations into the Transmitting and
  Receiving Qualities of Antennae," *Nova Acta R. Soc. Sci. Upsala.,*
  ser. IV, vol. 11, 1938 (short-antenna capacitive impedance).
