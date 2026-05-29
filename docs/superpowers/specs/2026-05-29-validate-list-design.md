# `yee validate --list` registered-case inventory — Design Spec

**Phase:** 1.validation.4  
**ADR:** ADR-0082  
**Date:** 2026-05-29  
**Status:** Accepted

---

## 1. Goal

Let a user (and CI) see the validation-case inventory in < 1 s, without the
~7–8 min `mom-001` solve that `run_all()` incurs today. Ship
`yee validate --list`: print every registered case (id, solver, policy,
description) and exit 0, running **no** solver.

Secondary goal: make the case set a **single source of truth** so the listed
inventory provably matches what `run_all()` executes (no metadata drift).

## 2. Public API (yee-validation)

```rust
/// Solver family a validation case belongs to. Mirrors the `yee validate
/// <target>` prefix routing (`mom-*`, FDTD-family, `fem-*`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Solver { Mom, Fdtd, Fem }

/// How a case behaves inside [`Report::run_all`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum ExecutionPolicy {
    /// Executed every run; its pass/fail counts toward `has_failures()`.
    Run,
    /// Registered `Skipped` to keep the default path fast; the strict gate
    /// runs via the named `#[ignore]`'d test (e.g. fem-eig-003 ~31 min).
    SkippedWallTime,
    /// Registered `Skipped` because the gate is open / deferred to a future
    /// phase (e.g. fem-eig-006 `|S11|≈0.955`, ADR-0070).
    SkippedGateOpen,
}

/// Non-executing metadata for one registered case.
#[derive(Debug, Clone, Serialize)]
pub struct CaseDescriptor {
    pub id: &'static str,
    pub solver: Solver,
    pub description: &'static str,
    pub policy: ExecutionPolicy,
}

/// All registered cases as descriptors — no solver is executed.
pub fn list_cases() -> Vec<CaseDescriptor>;
```

## 3. Single-source registry

Replace the `vec![run_mom_001(), …]` literal in `run_all()` with one private:

```rust
fn case_registry() -> Vec<(CaseDescriptor, fn() -> CaseResult)> {
    vec![
        (CaseDescriptor { id: "mom-001", solver: Solver::Mom,
            description: "Half-wave dipole impedance vs NEC-4 87 + j41 Ω (~7–8 min)",
            policy: ExecutionPolicy::Run }, run_mom_001 as fn() -> CaseResult),
        // … one row per case …
    ]
}
```

- `Report::run_all()` → `case_registry().into_iter().map(|(_, f)| f()).collect()`.
  Output is **byte-for-byte equivalent** to the current `run_all()` (same order,
  same `CaseResult`s — each runner is the existing wrapper, unchanged).
- `list_cases()` → `case_registry().into_iter().map(|(d, _)| d).collect()` — the
  runners are never invoked.

All existing `fn run_*() -> CaseResult` wrappers keep their current bodies; the
only change is that the registry pairs each with a descriptor. The descriptors'
`policy`/`solver`/`description` are authored to match each wrapper's actual
behaviour (see §6 gate).

## 4. CLI (`yee validate --list`)

- Add a boolean `--list` flag to the `validate` subcommand (clap), independent
  of the existing `<target>` arg and `--json`.
- `--list` handling short-circuits **before** `Report::run_all()`:
  1. `let cases = yee_validation::list_cases();`
  2. filter by `case_matches_target(d.id, target)` (reuse the existing fn — for
     the no-target/`all` case this is all rows).
  3. print a fixed-width table: `CASE | SOLVER | POLICY | DESCRIPTION`
     (`print_human_report`-style column sizing), or skip with a header line.
  4. return `ExitCode::SUCCESS` (listing never "fails").
- Without `--list`, behaviour is exactly as today (`run_validate`).
- Update the crate-level `//! yee validate …` synopsis to mention `--list`.

## 5. Files

- `crates/yee-validation/src/lib.rs` — `Solver`, `ExecutionPolicy`,
  `CaseDescriptor`, `case_registry()`, `list_cases()`, `run_all()` refactor.
- `crates/yee-validation/tests/integration.rs` — registry/list unit tests (§6).
- `crates/yee-cli/src/main.rs` — `--list` flag + handler + synopsis.
- `crates/yee-cli/tests/cli_validate.rs` — `--list` smoke test.

## 6. Definition of Done (machine-checkable)

1. `cargo fmt --check --all` exits 0.
2. `cargo clippy -p yee-validation -p yee-cli --all-targets -- -D warnings` exits 0.
3. `cargo test -p yee-validation -p yee-cli --release` exits 0 (default cases;
   `mom-001` stays `#[ignore]`'d — `--list` and the new unit tests must NOT pull
   it in).
4. New unit test `list_cases_matches_registry`: `list_cases()` is non-empty,
   contains `mom-001` and `fem-eig-006`, and `list_cases()` ids (in order) equal
   the ids `run_all()` would produce — verified WITHOUT running `run_all()` by
   asserting against the same `case_registry()` length/ids, plus a fast spot
   check that every `SkippedWallTime`/`SkippedGateOpen` descriptor's case returns
   `CaseStatus::Skipped` when its (fast) runner is invoked directly (exclude
   `mom-001`).
5. New CLI smoke test `yee_validate_list_runs` (NOT `#[ignore]`'d): runs in < 1 s,
   `yee validate --list` exits 0, stdout contains `mom-001` and `fem-eig-006`.
6. `cargo run -p yee-cli -- validate --list` prints the table and exits 0 in < 1 s.

## 7. Out of scope

Any gate/tolerance/`Skipped` change; `--list --json`; wall-time estimates;
touching solver/FDTD/FEM code. These are explicitly deferred.
