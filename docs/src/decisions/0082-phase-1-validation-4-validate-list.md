# ADR-0082: `yee validate --list` registered-case inventory (Phase 1.validation.4)

**Status:** Accepted  
**Date:** 2026-05-29  
**Supersedes:** none  
**Related:** ADR-0008 (validation aggregator JSON + PNG), ADR-0067 (FEM-eig
aggregator wiring + `yee validate fem`), ADR-0074 (FDTD aggregator gates)

---

## Context

`yee validate <mom|fdtd|fem|all>` is the single entry point to the validation
aggregator (`yee_validation::Report::run_all`). To discover *which* cases exist
and whether each runs or is deferred, a user today must **execute** the
aggregator — and `run_all()` runs `mom-001`, the 24×176 dipole solve that takes
~7–8 min in `--release` (CLAUDE.md §4/§10). There is no way to see the gate
inventory cheaply: not from the CLI, not from CI, not from a `--help`.

The case set is also currently expressed only as a `vec![run_mom_001(), …]`
literal inside `run_all()`, where each entry is a `fn() -> CaseResult` that
**executes** its case. There is no separate, non-executing description of the
registered cases, so any "list" feature risks drifting out of sync with what
`run_all()` actually runs — exactly the metadata-vs-behaviour drift that bit the
ADR-0067 review (a doc claimed a registration check the code never performed).

## Decision

Ship **`yee validate --list`**: print the registered-case inventory (id, solver,
one-line description, execution policy) and exit 0 **without running any
solver**. Back it with a **single source of truth** so the list cannot drift
from `run_all()`:

- `yee-validation` gains three public, documented items:
  - `enum Solver { Mom, Fdtd, Fem }` — the solver family a case belongs to
    (matches the `yee validate <target>` prefixes `mom-*` / FDTD-family / `fem-*`).
  - `enum ExecutionPolicy { Run, SkippedWallTime, SkippedGateOpen }` — whether
    the case executes in `run_all()` (`Run`), is registered `Skipped` to keep
    the default path fast (`SkippedWallTime`, e.g. `fem-eig-003` ~31 min), or is
    `Skipped` because its gate is open/deferred (`SkippedGateOpen`, e.g.
    `fem-eig-006` `|S11|≈0.955`).
  - `struct CaseDescriptor { id, solver, description, policy }` — non-executing
    metadata for one case.
- A single private `case_registry() -> Vec<(CaseDescriptor, fn() -> CaseResult)>`
  becomes the **one** place the case set is declared. `Report::run_all()` is
  rewritten to map the registry through execution (behaviour unchanged); a new
  `pub fn list_cases() -> Vec<CaseDescriptor>` maps the same registry to its
  descriptors **without** calling any runner. Because both derive from one
  `case_registry()`, the inventory and the executed set cannot diverge.
- `yee-cli`: `yee validate --list` prints a fixed-width table (CASE / SOLVER /
  POLICY / DESCRIPTION) and returns `ExitCode::SUCCESS`. `--list` short-circuits
  before `run_all()`, so it is instant. `--list` and a target are independent:
  `yee validate fem --list` lists only `fem-*`; `yee validate --list` lists all.

### Why a static `policy` label rather than deriving it from a dry execution

The pass/fail of a case is only knowable by running it; `--list` deliberately
does **not** run anything, so it reports *policy* (Run / Skipped-why), not
pass/fail. The case **set** is drift-proof (single `case_registry()`); the
`policy` field is human-maintained descriptive metadata kept truthful by a
doc-comment contract and a unit test that asserts the `Skipped*` cases actually
return `Skipped` (fast cases only — `mom-001` excluded by `#[ignore]`).

## Consequences

**Ships:**
- `Solver`, `ExecutionPolicy`, `CaseDescriptor`, `list_cases()` in `yee-validation`.
- `case_registry()` single-source registry; `run_all()` refactored onto it
  (output byte-for-byte equivalent to before).
- `--list` flag on `yee validate` in `yee-cli`.
- Unit test (registry non-empty; `list_cases()` ids == executed ids for the set;
  `Skipped*` policy matches actual `Skipped` status on fast cases) and a CLI
  smoke test (`yee validate --list` exits 0, prints `mom-001` + `fem-eig-006`,
  completes in < 1 s).

**Not in scope:** changing any gate, tolerance, or `Skipped` decision; JSON
output for `--list` (the existing `--json` runs the aggregator — a non-executing
`--list --json` is a trivial follow-on if wanted); per-case wall-time estimates.

**No new dependency.** Lane: `crates/yee-validation/{src,tests}/**`,
`crates/yee-cli/{src,tests}/**`. No solver, FDTD, or FEM code touched.

---

## References

- ADR-0008 (validation aggregator JSON + PNG), ADR-0067, ADR-0074.
- CLAUDE.md §4 (validation gates), §10 (mom-001 ~7–8 min wall time).
- `docs/superpowers/specs/2026-05-29-validate-list-design.md`
- `docs/superpowers/plans/2026-05-29-validate-list.md`
