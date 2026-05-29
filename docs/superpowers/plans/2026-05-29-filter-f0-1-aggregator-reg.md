# Filter Phase F0.1 — aggregator registration — Implementation Plan

**Spec:** `2026-05-29-filter-f0-1-aggregator-reg-design.md` · **ADR:** ADR-0085

## Lane
`crates/yee-validation/**`, `crates/yee-cli/**`, root `Cargo.toml` (yee-validation
dep entries only). Out of lane (yee-synth/yee-filter source, other crates) →
finding, not fix.

## Base
Worktree `worktrees/filter-f01`, branch `feature/filter-f0-1-aggregator-reg`,
base `53df105`.

## Pattern files
- `crates/yee-validation/src/lib.rs` — `Solver`/`ExecutionPolicy`/`CaseDescriptor`,
  `case_registry()`, the existing `run_fem_eig_00N()` wrappers (imitate their
  shape for `run_synth_001` etc.), `list_cases()`.
- `crates/yee-cli/src/main.rs` — `ValidateTarget`, `case_matches_target`,
  `run_validate_list` Solver match, the `Validate` synopsis.
- `crates/yee-cli/tests/cli_validate.rs` — smoke-test style.

## Steps
1. `yee-validation/Cargo.toml`: add `yee-synth = { workspace = true }`,
   `yee-filter = { workspace = true }`.
2. `Solver::Synth` variant (doc it; serialize is the variant name).
3. `run_synth_001()/run_synth_002()/run_filt_001()` wrappers per spec — each
   calls into yee-synth/yee-filter, folds the pass check into `CaseResult`
   (id, description, status, notes, wall_time, empty plot_paths). Reuse the
   exact published numbers from ADR-0084 / the F0 crate tests.
4. Append the three to `case_registry()` (after `run_fem_eig_006`), each with
   `CaseDescriptor { id, solver: Solver::Synth, description, policy: Run }`.
5. yee-cli: `ValidateTarget::Synth` + `case_matches_target` arm (`synth-`/`filt-`)
   + `run_validate_list` `Solver::Synth => "Synth"` + synopsis/doc update.
6. Tests: `yee-validation` `synth_filt_cases_registered_and_pass`; `yee-cli`
   `yee_validate_synth_list_runs`.

## Verify (exit 0; nice -n 19, --jobs 2; NO --workspace, NO mom-001/fem-eig-003)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-validation -p yee-cli --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-validation --jobs 2 --test integration synth_filt   # + the new test
nice -n 19 cargo test -p yee-cli --jobs 2 --test cli_validate yee_validate_synth
nice -n 19 cargo run -p yee-cli --jobs 2 -- validate synth --list
```
The synth drivers are ms-scale; do NOT run the full `cargo test -p yee-validation`
unfiltered if it would pull the ~8 min mom-001 / ~31 min fem-eig-003 integration
suites — target by test name.

## Escape hatch
Blocked >15 min (e.g. a yee-validation→yee-synth dep cycle, or run_all ordering
breaks an existing `run_all_includes_*` test) → stop, surface the exact error.
Do NOT weaken any existing test.

## Done when
DoD 1–6 pass; `git diff --stat 53df105..HEAD` shows only yee-validation, yee-cli,
root Cargo.toml + the 3 committed docs; the existing 19 cases unchanged.
