# `yee validate --list` — Implementation Plan

**Spec:** `docs/superpowers/specs/2026-05-29-validate-list-design.md`  
**ADR:** ADR-0082  
**Phase:** 1.validation.4  
**Date:** 2026-05-29

---

## Lane

`crates/yee-validation/{src,tests}/**`, `crates/yee-cli/{src,tests}/**`. Out of
lane (docs already committed, any solver/FDTD/FEM/GUI code) → finding, not fix.

## Base

Worktree `worktrees/validation-list`, branch
`feature/phase-1-validation-4-validate-list`, base `origin/main` `2cb0259`.

## Steps

1. **yee-validation types.** Add `pub enum Solver { Mom, Fdtd, Fem }`,
   `pub enum ExecutionPolicy { Run, SkippedWallTime, SkippedGateOpen }`, and
   `pub struct CaseDescriptor { pub id: &'static str, pub solver: Solver,
   pub description: &'static str, pub policy: ExecutionPolicy }`. Derive
   `Debug, Clone, Serialize` (Copy + PartialEq, Eq on the enums). Doc every item
   (`#![warn(missing_docs)]`).

2. **`case_registry()`.** Add private
   `fn case_registry() -> Vec<(CaseDescriptor, fn() -> CaseResult)>` listing
   **every** case currently in `run_all()`, in the same order, each paired with
   its existing `run_*` wrapper coerced `as fn() -> CaseResult`. Author each
   descriptor's `solver`/`policy`/`description` to match the wrapper's real
   behaviour — read each wrapper to set `policy`:
   - a wrapper that calls a driver and propagates its status → `Run`;
   - a wrapper hardcoding `Skipped` for wall-time → `SkippedWallTime`;
   - a wrapper hardcoding `Skipped` for an open gate → `SkippedGateOpen`.
   (As of base `2cb0259`: `mom-001/002/003` Run; the FDTD-family per their
   wrappers; `fem-eig-001/002/004/005` Run; `fem-eig-003` SkippedWallTime;
   `fem-eig-006` SkippedGateOpen. VERIFY against the actual wrapper bodies —
   do not trust this list blindly.)

3. **Refactor `run_all()`** to
   `case_registry().into_iter().map(|(_, f)| f()).collect()` (preserve the
   surrounding `Report { generated_at, git_sha, cases }`). Confirm the case
   ordering and contents are unchanged.

4. **`list_cases()`.** Add
   `pub fn list_cases() -> Vec<CaseDescriptor>` =
   `case_registry().into_iter().map(|(d, _)| d).collect()`. No runner invoked.

5. **CLI `--list`.** In `yee-cli`, add `--list` (bool) to the `validate`
   subcommand. When set, short-circuit before `run_validate`/`run_all`: collect
   `list_cases()`, filter by the existing `case_matches_target`, print a
   fixed-width `CASE | SOLVER | POLICY | DESCRIPTION` table (reuse the column
   sizing idiom from `print_human_report`), return `ExitCode::SUCCESS`. Update
   the crate-level `//! yee validate …` synopsis to mention `--list`.

6. **Tests.**
   - `crates/yee-validation/tests/integration.rs`: `list_cases_matches_registry`
     — `list_cases()` non-empty; contains `mom-001` + `fem-eig-006`; ids match
     the executed-set ids by construction; for every descriptor whose policy is
     `SkippedWallTime`/`SkippedGateOpen` AND whose case is fast (exclude
     `mom-001`), invoking the matching public driver/ wrapper yields
     `CaseStatus::Skipped` (spot-check fem-eig-006 via its wrapper path; do NOT
     call `run_all`).
   - `crates/yee-cli/tests/cli_validate.rs`: `yee_validate_list_runs` (NOT
     ignored) — `yee validate --list` exits 0, stdout has `mom-001` +
     `fem-eig-006`, < 1 s. Add `--list` mention to the existing help test if it
     asserts on synopsis text.

7. **Verify (run all, expect exit 0):**
   ```
   cargo fmt --check --all
   cargo clippy -p yee-validation -p yee-cli --all-targets -- -D warnings
   cargo test -p yee-validation -p yee-cli --release
   cargo run -p yee-cli -- validate --list
   ```
   Keep builds light: `nice -n 19 … --jobs 2`. The CPU box is constrained —
   do NOT run `cargo test --workspace` (pulls the ~8 min mom-001 + ~31 min
   fem-eig-003) and do NOT `--include-ignored`.

## Escape hatch

Blocked > 15 min (e.g. `fn` -pointer coercion fights the existing wrapper
signatures, or `run_all()` output changes) → stop and surface the specific
blocker; do not change gate behaviour or weaken any test to get green.

## Done when

DoD items 1–6 in the spec all pass; lane respected (`git diff --stat
2cb0259..HEAD` shows only the 4 code/test files + the 3 docs already committed);
`run_all()` output unchanged.
