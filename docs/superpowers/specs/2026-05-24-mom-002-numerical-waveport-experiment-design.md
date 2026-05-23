# mom-002 numerical-microstrip-wave-port — bounded experiment

**Status:** Draft (bounded experiment, not a committed feature)
**Owner:** TBD
**Phase:** MoM beachhead follow-on (leverages the completed Phase 1.3.1.1
cross-section eigensolver).
**Type:** TIGHT BOUNDED EXPERIMENT — the deliverable is the *result*
(a Z_in comparison, or a documented wiring-blocker finding), not a
committed gate change.

## 1. Goal / hypothesis

mom-002 (50 Ω microstrip Z₀ on FR-4, `L = 82 mm`) passes only in a loose
±5% tripwire band at `|Z_in| ≈ 674 Ω` (≈13× the `Z_0 ≈ 51 Ω` target);
10 forensic tracks exonerated the kernel (within 1.83% of HJ ε_eff) and
localised the residual to **delta-gap port-excitation modeling** (the
line is currently fed by the TEM-smoothed delta-gap
`z_in_with_greens_tem`). The cross-section eigensolver — now complete and
**FR-4-validated** (1.39% vs the verified reference) — computes the
microstrip quasi-TEM modal field, and the `WavePort` **`Numerical2D`**
arm injects a cross-section modal `E_t` into the MoM port-edge RHS.

**Hypothesis:** exciting the mom-002 line with a *numerical microstrip
modal* wave-port (the FR-4 cross-section mode) instead of the
TEM-smoothed delta-gap reduces the port-modeling residual and brings
`|Z_in|` closer to `Z_0 ≈ 51 Ω`. This is a fresh angle the 10 forensic
tracks never tried (they worked the delta-gap / kernel / Greens).

## 2. Bounded scope + the explicit STOP conditions

This is an experiment on a heavily-analysed case; it is **bounded**:

- **Phase A — FEASIBILITY (hard 30-min cap).** Can a microstrip
  cross-section (`TriMesh2D`: FR-4 substrate + signal strip + ground +
  air, the transverse plane of the mom-002 line) be built, solved via
  `NumericalCrossSection`, and its modal `E_t` fed to the mom-002 MoM
  line through `WavePort::with_numerical_cross_section` +
  `ModalDistribution::Numerical2D` at all? The `Numerical2D` arm was
  validated for a homogeneous waveguide-TE10, NOT a microstrip into the
  planar-MoM RHS — the cross-section↔strip-current coupling (mapping the
  2-D transverse modal field onto the 2.5-D RWG port edges) may need
  glue that does not exist. **If the coupling cannot be wired cleanly in
  30 min, STOP** — the deliverable is then a documented FINDING: exactly
  what the `Numerical2D` arm lacks to support microstrip-into-planar-MoM
  ports (a scoping result for a future port-infrastructure track). Do
  NOT force it.
- **Phase B — COMPARISON (only if A succeeds).** Extract `|Z_in|` with
  the numerical port; compare to the delta-gap baseline (674 Ω) and the
  HJ target (≈51 Ω). Report the numbers.

## 3. Deliverable (either branch is a success)

- **If wired (Phase B):** a `Z_in`-comparison report — numerical-port
  `|Z_in|` vs delta-gap 674 Ω vs HJ ≈51 Ω — as a non-failing diagnostic
  test (do NOT change the mom-002 gate or its tripwire band; this is an
  experiment, not a re-gate). If the numerical port clearly improves
  `|Z_in|` toward `Z_0`, that is a strong result worth a follow-on track
  to adopt it; if not, the comparison documents that the residual is not
  (only) the port excitation.
- **If not wired (Phase A stop):** a documented finding — the specific
  glue the `Numerical2D` arm needs for microstrip-into-planar-MoM ports.

## 4. NON-NEGOTIABLE constraints

- Do **NOT** re-open the mom-002 kernel / Greens / forensic analysis —
  the kernel is exonerated; this experiment is *only* about the port
  excitation. If the experiment points back at the kernel/Greens, STOP +
  document (do not chase).
- Do **NOT** change the mom-002 gate, its 674 Ω regression, or its
  tripwire band. The comparison ships as a separate non-failing
  diagnostic.
- Do **NOT** edit the cross-section eigensolver (`reference.rs`,
  `assembly.rs`, the validated solve) — consume it read-only.
- No new `Cargo.toml` dependency. Lint floor clean. No regression to
  mom-001/mom-002/mom-003 existing behaviour.

## 5. Definition of done

DoD-1. Phase A feasibility decided within the cap: either the coupling
is wired (→ B) or the blocker is documented (→ stop, finding shipped).
DoD-2. (If B) `|Z_in|` numerical-port vs 674 Ω vs 51 Ω reported as a
non-failing diagnostic test.
DoD-3. No regression (mom-001/002/003 gates + behaviour unchanged); the
mom-002 gate + band untouched; eigensolver untouched. Lint clean.
DoD-4. A clear recommendation: adopt-the-numerical-port (follow-on
track) / residual-is-not-the-port / port-infra-glue-needed.

## 6. References

* ADR-0036 (mom-002 validation reframe), ADR-0037 (R1 retraction),
  the mom-002 forensic-track summary in ROADMAP.
* `crates/yee-validation/src/lib.rs` (`z_in_with_greens_tem`, the
  mom-002 case + constants), `crates/yee-mom/src/ports.rs`
  (`WavePort`, `ModalDistribution::Numerical2D`,
  `with_numerical_cross_section`, `e_tangential_at`), Track GGGGGGGG
  (the Numerical2D arm).
* The completed cross-section eigensolver (ADRs 0050–0058).
