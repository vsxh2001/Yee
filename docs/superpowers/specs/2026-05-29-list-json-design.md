# `yee validate --list --json` machine-readable inventory — Design Spec

**Phase:** 1.validation.5  
**ADR:** ADR-0083  
**Date:** 2026-05-29  
**Status:** Accepted

---

## 1. Goal

Let CI/tooling read the validation-case inventory as data, without running any
solver. Make `yee validate --list --json` emit the filtered `CaseDescriptor`
list as a pretty JSON array (exit 0). Closes the `--list --json` follow-on
deferred by ADR-0082.

## 2. Behaviour

- `yee validate --list --json [target]` → JSON array of the cases matching
  `target` (or all), pretty-printed, exit 0, no solver run. Each element:
  ```json
  { "id": "mom-001", "solver": "Mom", "description": "…", "policy": "Run" }
  ```
  `solver` ∈ `"Mom" | "Fdtd" | "Fem"`; `policy` ∈
  `"Run" | "SkippedWallTime" | "SkippedGateOpen"` (default enum-name
  serialization, already derived).
- `yee validate --list` (no `--json`) → unchanged fixed-width table.
- `yee validate <target> [--json]` (run path) → unchanged.

## 3. Change

`crates/yee-cli/src/main.rs` only:
- `fn run_validate_list(target: ValidateTarget, json: bool) -> ExitCode`
  (add the `json` param). Build the filtered `Vec<CaseDescriptor>` exactly as
  today; then:
  - `if json { println!("{}", serde_json::to_string_pretty(&cases).expect("CaseDescriptor is Serialize")); return ExitCode::SUCCESS; }`
  - else the existing table path.
- `Command::Validate { target, json, list }` arm: call
  `run_validate_list(target, json)` when `list`.
- Update the crate-level `//!` synopsis and the `Validate` doc-comment to note
  `--list --json` emits the inventory as JSON.

`serde_json` is already a `yee-cli` dependency; `CaseDescriptor` already derives
`Serialize` (ADR-0082). No `yee-validation` change.

## 4. Definition of Done (machine-checkable)

1. `cargo fmt --check --all` exits 0.
2. `cargo clippy -p yee-cli --all-targets -- -D warnings` exits 0.
3. `cargo test -p yee-cli --release` exits 0.
4. New CLI test `yee_validate_list_json_runs` (NOT `#[ignore]`'d, < 1 s):
   `yee validate --list --json` exits 0; stdout trims to a string starting `[`
   and ending `]`; contains `"mom-001"`, `"fem-eig-006"`, `"SkippedGateOpen"`,
   and `"Run"`. (Substring assertions — no serde_json dep needed in the test.)
5. `cargo run -p yee-cli -- validate --list --json` prints a JSON array and
   exits 0 in < 1 s; `cargo run -p yee-cli -- validate fem --list --json` prints
   only `fem-*` cases.
6. `yee validate --list` (no `--json`) output is unchanged from ADR-0082.

## 5. Out of scope

Run-path `--json` `Report` schema changes; wall-time fields; any gate/tolerance
change; touching `yee-validation` or any solver crate.
