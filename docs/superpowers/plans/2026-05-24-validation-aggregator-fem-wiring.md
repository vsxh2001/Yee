# Validation aggregator FEM wiring — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-validation-aggregator-fem-wiring-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-validation/src/**`, `crates/yee-validation/tests/**`,
`crates/yee-cli/src/**`, `crates/yee-cli/tests/**`, `validation/README.md`.
NOTHING else.
**Out of lane** (findings, not fixes): `crates/yee-fem/**` (FEM drivers
consumed AS-IS — public + passing; do NOT edit their logic),
`crates/yee-fdtd/**` (the FDTD `Skipped` stubs + cavity gates need
new yee-fdtd-lane drivers — a FOLLOW-ON slice, leave them). No `Cargo.toml`
dependency.

## Step ladder

### S1 — read the registry + a driver + the test pattern
Read `crates/yee-validation/src/lib.rs`: the `run_all` registry (~L130),
ONE FEM driver end-to-end (`run_fem_eig_001_rectangular_cavity` ~L1481 +
its `FemEig…ValidationResult` struct), the hardcoded FDTD `Skipped` stubs
(~L1264, to mirror the `CaseResult` shape). Read `tests/integration.rs`
(`.find()`-based, non-breaking) + `docs/src/decisions/0008-*.md` (aggregator
shape). Measure each FEM driver's mesh (`FEM_EIG_00N_N{X,Y,Z}` constants) +
note which are fast (default path) vs heavy (wall-time-gated like mom-001).

### S2 — register the FEM gates (fold pattern)
For each public FEM-eig driver, add a `run_fem_eig_00N()` wrapper folding
its result into `CaseResult { id: "fem-eig-00N", status, notes, wall_time }`;
push into `run_all`. Fast gates (fem-eig-001 ≈7 s) in the default path;
heavy ones registered with the existing slow-case discipline (do NOT bloat
default `yee validate all` — gate/feature/`#[ignore]`-class as the suite
already does for mom-001; document any deferral in the case `notes`).

### S3 — CLI `fem` target
`crates/yee-cli/src/main.rs`: add `Fem` to `ValidateTarget` + a `fem-*`
prefix arm in `case_matches_target`. `yee validate fem` runs the FEM cases.

### S4 — de-stale the README
`validation/README.md`: replace the forbidden Balanis `73 + j42 Ω` dipole
reference with NEC-4 `87 + j41 Ω` (CLAUDE.md §4 / ADR-0005); reconcile the
documented case list with the actual `run_all` registry (drop fictional
`validation/cases/phase-*/` claims).

### S5 — tests
Extend `tests/integration.rs`: assert `run_all().cases` contains
`fem-eig-001` with `Passed`. Add a `yee-cli` assert_cmd test: `yee validate
fem --json` exits 0 + JSON contains `fem-eig-001`.

## Verification (run in worktree; all exit 0)
```
cargo fmt --check --all
cargo clippy -p yee-validation -p yee-cli --all-targets -- -D warnings
cargo test -p yee-validation        # incl. the new aggregator assertion (fast FEM gates)
cargo test -p yee-cli               # incl. the new `validate fem` test
git diff --stat -- crates/yee-fem crates/yee-fdtd '**/Cargo.toml'   # MUST be empty
grep -c "73 + j42\|73+j42" validation/README.md   # MUST be 0
```
(Do NOT run mom-001 / the heaviest FEM gates in the default test path —
they're wall-time-gated; rely on CI/`--ignored` for those.)

## Escape-hatch
- If a FEM driver's result struct does NOT cleanly fold into `CaseResult`
  (field mismatch), surface it — do NOT change the driver. If `run_all`
  becomes multi-minute, gate the slow gates out of the default path
  (don't bloat CI) + document. Blocked >15 min → surface + stop.
- Do NOT touch yee-fem/yee-fdtd logic. Synchronous; no Monitor/ScheduleWakeup;
  no sub-agents.

## Out-of-scope (findings, not fixes)
* New public yee-fdtd drivers to un-`Skip` cpml/ntff/dispersive + register
  fdtd-201/.x — a FOLLOW-ON slice (yee-fdtd lane).
* Python bindings for the validation runners; Touchstone v2.
