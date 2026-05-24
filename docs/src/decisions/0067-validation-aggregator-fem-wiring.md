# ADR-0067 — register the shipped FEM-eig gates in the validation aggregator + `yee validate fem`

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** 1.validation (aggregator truth-up)

## Context

A read-only survey of the untouched subsystems (after saturating the
MoM-port, FDTD-cavity, and frontend-plotting areas) found the project's
headline command is **lying**: `yee_validation::Report::run_all` registers
only 6 cases — mom-001/002/003 plus cpml-001/ntff-001/dispersive-001, the
latter three hardcoded `CaseStatus::Skipped` despite their physics gates
passing as `#[cfg(test)]` tests. Meanwhile the **entire FEM eigenmode
suite** — `fem-eig-001`..`006`, the most active subsystem (ADRs 0029-0058),
each with a **public, passing** driver (`run_fem_eig_00N_*` in
`yee-validation`) and a dedicated test — is absent from `run_all`, and the
CLI `ValidateTarget` enum has no `fem` target. Separately,
`validation/README.md` cites the **CLAUDE.md §4-FORBIDDEN** Balanis
`73 + j42 Ω` dipole reference (must be NEC-4 `87 + j41 Ω`).

## Decision

Wire the shipped, public, passing FEM-eig gates into `Report::run_all`
(folding their `FemEig…ValidationResult` into `CaseResult` — the structs
already share `status`/`notes`/`wall_time_seconds`), add a `fem` target to
`yee validate`, and de-stale `validation/README.md` (NEC-4 reference;
reconcile the case list). Walking-skeleton-first: the fast gates in the
default path, heavy ones registered with the `mom-001` wall-time
discipline so the default `yee validate all` stays fast. The FDTD
`Skipped` stubs + the fdtd-201/.x cavity gates need NEW public drivers in
the yee-fdtd lane — a follow-on slice, out of scope here.

## Rationale

(1) **Highest value × dispatchability available.** It makes the headline
command tell the truth (today it hides the most-active subsystem + 2
shipped FDTD milestones, and 3 cases lie as Skipped), and CLAUDE.md §4
calls the aggregator the single source of truth. Yet the drivers are
already public + passing and the result structs already match
`CaseResult`, so it is a fold + two CLI enum arms — no new physics, no new
dependency, no toolchain.

(2) **Genuinely fresh + low-risk.** Validation-aggregator integration
touches none of the three saturated areas and none of the three quagmires
(mom-002/003 accuracy, FDTD subgrid Q6, FEM real-waveguide-port — note the
*passing* fem-eig gates are unrelated to the saturated modal-projection
wave-port). Self-hosting validation contract (the feature IS the harness).

(3) **Fixes a §4 violation in the docs** — the forbidden `73 + j42`
dipole reference, which CLAUDE.md §4 + ADR-0005 explicitly prohibit.

## Consequences

* `yee validate` (and `yee validate fem`) report the FEM eigenmode suite;
  the default path stays fast (heavy gates wall-time-gated like mom-001).
* New aggregator + CLI tests assert fem-eig-001 is registered + Passed.
* `validation/README.md` de-staled (NEC-4 reference; real case list).
* Follow-on (yee-fdtd lane): public drivers to un-`Skip` cpml/ntff/dispersive
  + register fdtd-201/.x.
* No yee-fem / yee-fdtd logic change; no new dependency.

## References

* The fresh-subsystem survey (2026-05-24, read-only): `yee-validation/src/lib.rs`
  (`run_all` ~L130, FEM drivers ~L1481-3875, FDTD Skipped stubs ~L1264),
  `yee-cli/src/main.rs` (`ValidateTarget` ~L316), `tests/integration.rs`
  (`.find()`-based), `validation/README.md` (the forbidden reference).
* ADR-0008 (aggregator shape), ADR-0005 + CLAUDE.md §4 (NEC-4-only mom-001).
* Spec + plan (2026-05-24).
