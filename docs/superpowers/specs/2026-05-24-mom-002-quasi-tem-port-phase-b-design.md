# mom-002 quasi-TEM numerical wave-port — Phase B (now unblocked)

**Status:** Draft (experiment continuation + a production wiring step)
**Owner:** TBD
**Phase:** MoM beachhead follow-on — continues ADR-0059's experiment,
now that Phase 1.3.1.2 (ADR-0060) made the quasi-TEM mode selectable.
**Type:** production wiring (Part A — reusable) + bounded experiment
(Part B — the Z_in comparison, deliver either way).

## 1. Goal

The ADR-0059 experiment stopped at Phase A only because the eigensolver
could not find the microstrip quasi-TEM mode. Phase 1.3.1.2 fixed that
(`solve_dense_mixed_quasi_tem`, HJ-validated 1.2%). Phase B:

- **Part A (production wiring, reusable):** make the quasi-TEM mode
  usable for wave-port excitation — add a quasi-TEM solve path to
  `NumericalCrossSection` that caches `mode_profile` (the E_t edge
  amplitudes the `Numerical2D` `WavePort` arm + `e_tangential_at`
  consume). Currently `solve` only runs the closed-guide
  `solve_dense_mixed` (First / Second order); add the quasi-TEM path.
- **Part B (the experiment payoff):** feed the mom-002 microstrip
  cross-section's quasi-TEM modal field to the mom-002 line via the
  `Numerical2D` arm, extract `|Z_in|`, compare to the delta-gap baseline
  (674 Ω) + the HJ target (≈51 Ω). Does a numerical microstrip port
  beat the delta-gap?

## 2. Approach

### Part A — `NumericalCrossSection` quasi-TEM path
`crates/yee-mom/src/ports.rs`: add a quasi-TEM solve path (a
`solve_quasi_tem(freq_hz)` method, or extend the `ElementOrder` /
add a mode selector — mirror the First path's `mode_profile` /
`tri_edges_cache` caching, but call `solve_dense_mixed_quasi_tem` instead
of `solve_dense_mixed`). `solve_dense_mixed_quasi_tem` returns the same
`MixedEigenSolution` (E_t eigenvector), so the existing scatter →
`mode_profile` + the `e_tangential_at` interpolation reuse unchanged.
First-order closed-guide `solve` stays the default + bit-identical.

### Part B — mom-002 numerical port + comparison
Extend `crates/yee-mom/tests/mom_002_numerical_waveport.rs` (the
ADR-0059 diagnostic, which already builds the microstrip cross-section +
attempts the Numerical2D coupling): now the cross-section quasi-TEM solve
SUCCEEDS → feed its `mode_profile` to the mom-002 line via
`WavePort::with_numerical_cross_section` + the Numerical2D RHS, LU-solve,
extract `|Z_in|`. Report vs 674 Ω (delta-gap) + 51 Ω (HJ). NON-FAILING
diagnostic — do NOT re-gate mom-002.

## 3. Bounded framing (the experiment part)

- **Part A feasibility (hard ~30-min cap on the coupling):** can the
  quasi-TEM `mode_profile` be cached + the Numerical2D RHS built for the
  mom-002 line + LU-solved to a finite `|Z_in|`? If the
  2-D-cross-section→2.5-D-RWG modal-RHS mapping for a MICROSTRIP (vs the
  waveguide-TE10 the arm was validated on) needs glue that doesn't
  exist → document the specific blocker + STOP (a finding). Do NOT force.
- **Part B comparison:** report the `|Z_in|` numbers. Either branch is a
  deliverable.

## 4. Definition of done

DoD-1. Part A: `NumericalCrossSection` quasi-TEM solve path caches
`mode_profile` (the closed-guide First/Second paths bit-identical /
unchanged). A smoke that the microstrip cross-section solves quasi-TEM +
populates `mode_profile`.
DoD-2. Part B: `|Z_in|` numerical-quasi-TEM-port vs 674 Ω vs 51 Ω
reported as a non-failing diagnostic, OR a documented coupling-blocker
finding (Part A feasibility stop).
DoD-3. No regression: mom-001/002/003 gates + behaviour unchanged;
mom-002 gate/tripwire untouched; the closed-guide eigensolver paths
bit-identical; the quasi-TEM selection + HJ gate unchanged. Lint clean.
DoD-4. Recommendation: adopt-the-numerical-port (does it beat
delta-gap?) / residual-not-the-port / remaining-glue-needed.

## 5. NON-NEGOTIABLE

- Do NOT re-open the mom-002 kernel / Greens / forensics; do NOT change
  the mom-002 gate or 674 Ω tripwire band.
- Do NOT alter the eigensolver `solve_dense_mixed` / `_quasi_tem` /
  `cutoff_candidates` / the verified `reference.rs` — consume read-only
  (Part A is a `ports.rs` consumer of the existing quasi-TEM solve).
- No new `Cargo.toml` dependency.

## 6. References

* ADR-0059 (the experiment, Phase-A finding), ADR-0060 (the quasi-TEM
  capability this unblocks Part B with).
* `crates/yee-mom/src/ports.rs` (`NumericalCrossSection::solve`, the
  First-path `mode_profile` caching to mirror; `WavePort` Numerical2D
  arm; `e_tangential_at`), `crates/yee-mom/src/eigensolver/solve.rs`
  (`solve_dense_mixed_quasi_tem` → `MixedEigenSolution`),
  `crates/yee-mom/tests/mom_002_numerical_waveport.rs` (the diagnostic
  to extend), `crates/yee-validation/src/lib.rs` (mom-002 constants).
