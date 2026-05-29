# ADR-0083: `yee validate --list --json` machine-readable inventory (Phase 1.validation.5)

**Status:** Accepted  
**Date:** 2026-05-29  
**Supersedes:** none  
**Related:** ADR-0082 (`yee validate --list` registered-case inventory — this
closes the `--list --json` follow-on it explicitly deferred)

---

## Context

ADR-0082 shipped `yee validate --list`: a human-readable table of the registered
validation cases (id / solver / policy / description) printed without running
any solver. It explicitly deferred a JSON form:

> Not in scope: … JSON output for `--list` (the existing `--json` runs the
> aggregator — a non-executing `--list --json` is a trivial follow-on if wanted).

Today `run_validate_list(target)` ignores the `--json` flag entirely, so
`yee validate --list --json` silently prints the human table. CI and tooling
that want the gate inventory as data must therefore either parse the fixed-width
table (brittle) or run the full aggregator (`--json` without `--list`, which
incurs the ~7–8 min `mom-001` solve). Both are bad options for a machine
consumer that only wants "what cases exist and will each run or skip?".

## Decision

Make `yee validate --list --json` emit the filtered `CaseDescriptor` list as a
pretty-printed JSON array, still running no solver:

- `run_validate_list` takes the `json: bool` already parsed on the `Validate`
  command and, when set, prints `serde_json::to_string_pretty(&cases)` (the
  `Vec<CaseDescriptor>` after the existing `case_matches_target` filter) instead
  of the table, returning `ExitCode::SUCCESS`.
- `CaseDescriptor`, `Solver`, and `ExecutionPolicy` already derive `Serialize`
  (ADR-0082), so each element serializes as
  `{"id":"mom-001","solver":"Mom","description":"…","policy":"Run"}`
  (`policy` ∈ `"Run" | "SkippedWallTime" | "SkippedGateOpen"`). No new
  serialization code, no new dependency (`yee-cli` already depends on
  `serde_json`).
- `yee validate --list` (no `--json`) is unchanged; `yee validate <target>
  --json` (the run path) is unchanged.

### Why a top-level array rather than a wrapped object

The run path emits a `Report` object (`{generated_at, git_sha, cases}`). The
`--list` path has no run metadata (nothing executed), so wrapping it in a
faux-`Report` would be misleading. A bare array of descriptors is the honest,
minimal shape; a consumer that wants metadata runs the aggregator.

## Consequences

**Ships:** `--json` honored by `run_validate_list`; an updated synopsis; a CLI
smoke test `yee_validate_list_json_runs` (not `#[ignore]`'d, < 1 s) asserting
exit 0, output parses as a JSON array, and contains `mom-001` + `fem-eig-006`
with their `policy` strings.

**Not in scope:** changing the run-path `--json` `Report` schema; per-case
wall-time numbers; any gate/tolerance change.

**No new dependency.** Lane: `crates/yee-cli/{src,tests}/**` only.

---

## References

- ADR-0082 (`yee validate --list`); ADR-0008 (validation aggregator JSON + PNG).
- `docs/superpowers/specs/2026-05-29-list-json-design.md`
- `docs/superpowers/plans/2026-05-29-list-json.md`
