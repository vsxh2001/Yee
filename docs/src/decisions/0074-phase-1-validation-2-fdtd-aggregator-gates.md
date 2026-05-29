# ADR-0074: Phase 1.validation.2 — FDTD Gate Integration into Aggregator

**Date:** 2026-05-27  
**Status:** Proposed  
**Deciders:** orchestrator  
**Supersedes:** —  
**Related:** ADR-0008 (validation aggregator), ADR-0062 (fdtd-201), ADR-0071 (fdtd-202)

---

## Context

`Report::run_all()` in `crates/yee-validation/src/lib.rs` has three FDTD gates
hardcoded to `CaseStatus::Skipped` since Phase 1.validation.1:

| ID | Description | Skipped note |
|----|-------------|--------------|
| `cpml-001` | CPML ≥ 30 dB attenuation vs PEC | "deferred to Phase 1.validation.2" |
| `ntff-001` | Broadside/endfire null ≥ 20 dB | "deferred to Phase 1.validation.2" |
| `dispersive-001` | Drude slab Fresnel ≤ 20 % error | "deferred to Phase 1.validation.2" |

The corresponding yee-fdtd integration tests (`cpml_reflection.rs`, `ntff_dipole.rs`,
`dispersive.rs`) all **pass** in `cargo test`. The only gap is that the aggregator
doesn't execute the physics — it just returns a hard-coded Skipped status.

This was deferred because the physics helpers were not yet public or the pattern for
inline-duplicating them was not established. Since `fdtd-202` (ADR-0071) demonstrated
the pattern clearly (duplicate physics inline, no cross-crate test deps), the deferral
is no longer justified.

---

## Decision

Port the three gates from Skipped to actual-running in `yee-validation/src/lib.rs`
following the `fdtd-202` pattern:

1. **cpml-001**: duplicate the PEC-vs-CPML peak-attenuation measurement inline.
2. **ntff-001**: duplicate the broadside/endfire NTFF ratio measurement inline.
3. **dispersive-001**: duplicate the Drude slab Fresnel reflection measurement inline.

All three are estimated < 10 s in debug mode and will NOT be `#[ignore]`-gated
in the aggregator.

Public API used (all already exported from `yee_fdtd`):
- `YeeGrid`, `CpmlParams`, `CpmlState`, `WalkingSkeletonSolver`
- `NtffParams`, `NtffState`
- `Material`, `MaterialMap`, `DispersiveState`
- `yee_fdtd::update::{update_h, update_e}` (pub fn)
- `yee_fdtd::sources::gaussian_pulse_ez` (pub fn)

No new dependencies; no yee-fdtd source changes.

---

## Consequences

**Positive:**
- `yee validate all` goes from 3 Skipped FDTD gates to 3 Passed.
- CI has real execution coverage of CPML absorption quality, NTFF correctness, and
  dispersive material Fresnel accuracy.
- Closes the "deferred to Phase 1.validation.2" note that has been in the codebase
  since Phase 1.validation.1.

**Negative / neutral:**
- yee-validation compile time increases slightly (more code).
- The physics is duplicated (by design — same reasoning as fdtd-202); if the underlying
  tests ever change parameters, the aggregator code needs a matching update.

---

## Alternatives considered

**A — depend on yee-fdtd test helpers**: Rust doesn't allow depending on `#[cfg(test)]`
modules from other crates; this would require making helpers public in yee-fdtd source,
polluting the API surface. Rejected (same reasoning as fdtd-202).

**B — call `cargo test` as a subprocess**: fragile, slow, hard to capture structured
results. Rejected.

**C — keep Skipped forever**: leaves a misleading gap in `yee validate all`. Rejected.
