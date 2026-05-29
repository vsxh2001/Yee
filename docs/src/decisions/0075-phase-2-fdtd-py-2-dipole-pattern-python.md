# ADR-0075: Phase 2.fdtd.py.2 — Python Dipole-Pattern Driver + fdtd-203 Gate

**Date:** 2026-05-27  
**Status:** Accepted  
**Deciders:** orchestrator  
**Supersedes:** —  
**Related:** ADR-0072 (fdtd-202 Python driver), ADR-0073 (fdtd-201 Python driver),
ADR-0074 (ntff-001 aggregator wiring)

---

## Context

`crates/yee-fdtd/src/driver.rs` contains `FdtdDriver` — a complete end-to-end
FDTD driver that combines a vacuum `YeeGrid`, CPML on all six faces, a
Hann-ramped sinusoidal J_z dipole current source, and an `NtffState` DFT
accumulator into a single `run()` call that produces a 37-point θ-cut of
`|E_θ|` normalized to its maximum.  `FdtdDriver` is already exported from
`yee_fdtd` and wrapped as `PyFdtdDriver` / `PyFdtdDriverConfig` /
`PyRadiationPattern` in `crates/yee-py/src/fdtd.rs` (advanced API).

The corresponding integration test `crates/yee-fdtd/tests/dipole_pattern.rs`
(Phase 2.fdtd.4) is `#[ignore]`-gated and validates the sin θ radiation
pattern against the Balanis §4.2 short-dipole reference.  It **passes** when
run with `--release --include-ignored`.

Two gaps remain:
1. No `run_dipole_pattern()` convenience function analogous to `run_cavity_q`
   / `run_cavity_resonance` — Python users can run `FdtdDriver` but must
   assemble it by hand and interpret the raw pattern.
2. No fdtd-203 gate in `Report::run_all()` — the dipole-pattern validation
   milestone is unregistered.

---

## Decision

**Phase 2.fdtd.py.2** closes both gaps following the established pattern
(ADR-0072 / ADR-0073):

1. **`run_dipole_pattern()` pyfunction** in `crates/yee-py/src/fdtd.rs`:
   defaults to 60³ grid, dx=5mm, 800 steps, 1 GHz — the same scenario as
   `dipole_pattern.rs`.  Returns `DipolePatternResult` with `passed: bool`,
   five θ sample points, and numpy-array accessors for the full θ-cut.

2. **`PyDipolePatternResult` struct** (`DipolePatternResult` in Python):
   read-only scalar properties `e_theta_{0,45,90,135,180}`, `passed`, and
   methods `theta_deg_array()` / `e_theta_array()`.

3. **fdtd-203 gate** registered in `Report::run_all()` as `CaseStatus::Skipped`
   (wall-time ~30 s release — same treatment as fdtd-201 / fdtd-201-x).

4. **Tutorial 12** (`docs/src/tutorials/12-fdtd-dipole-pattern-from-python.md`):
   quick-start + full gate run + matplotlib θ-cut overlay.

---

## Gate criteria

| θ    | Expected | Tolerance | Reference              |
|------|----------|-----------|------------------------|
|  0°  | 0        | < 0.15    | Balanis §4.2 null      |
| 45°  | 0.707    | ±0.15     | sin 45° = √2/2 ≈ 0.707 |
| 90°  | 1.0      | ±0.05     | Normalized broadside   |
| 135° | 0.707    | ±0.15     | sin 135° = 0.707       |
| 180° | 0        | < 0.15    | Balanis §4.2 null      |

Tolerances are deliberately loose: the 60³ grid at λ/dx = 60 is a coarse
approximation; NTFF integrates over a finite box.  Tighter tolerances belong
to a finer-grid gate in a future increment.

---

## Consequences

**Positive:**
- `from yee import run_dipole_pattern; assert r.passed` is the canonical
  one-liner for the Phase 2 NTFF radiation-pattern gate — aligns with the
  `run_cavity_q` / `run_cavity_resonance` idiom.
- fdtd-203 is registered in `yee validate all` (Skipped) so the milestone is
  visible in the validation dashboard.
- `DipolePatternResult.theta_deg_array()` / `.e_theta_array()` give Python
  notebooks direct numpy access for plotting θ-cuts — the advanced
  `PyFdtdDriver` API is now complemented by a high-level gate function.

**Negative / neutral:**
- Wall-time ~30 s release for the full 800-step run; kept Skipped in run_all.
- Physics is duplicated in the Python convenience function (by design — same
  reasoning as fdtd-202 / ADR-0071).

---

## Alternatives considered

**A — Promote `dipole_pattern.rs` test directly to the aggregator (inline):**
The test runs 60³ × 800 steps via `FdtdDriver`; in debug mode this is
~5–10 min, too slow for non-gated aggregator use.  Rejected for the same
reason fdtd-201 is Skipped.

**B — Expose only the already-wrapped `PyFdtdDriver` and document it:**
The raw driver forces every caller to assemble `FdtdDriverConfig` and interpret
the `RadiationPattern` themselves.  The convention established by fdtd-202 /
fdtd-201 is a convenience function with a `passed` flag.  Rejected for
consistency.
