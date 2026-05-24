# Validation aggregator — register the shipped FEM-eig gates + `yee validate fem` + de-stale

**Status:** Draft
**Owner:** TBD
**Phase:** 1.validation (aggregator truth-up)
**Type:** integration + staleness fix (high value × dispatchability)

## 1. Goal

`yee_validation::Report::run_all()` registers only 6 cases (mom-001/002/003
+ cpml-001/ntff-001/dispersive-001, the latter three hardcoded
`CaseStatus::Skipped`), so the project's headline command `yee validate all`
**hides the entire FEM eigenmode suite** — six gates (`fem-eig-001`..`006`)
that have **public, passing drivers** (`run_fem_eig_001_rectangular_cavity`
… `run_fem_eig_006_high_aspect_pml` in `crates/yee-validation/src/lib.rs`)
and dedicated passing tests — and reports no FDTD cavity milestones. The
CLI `ValidateTarget` enum has no `fem` target. Separately,
`validation/README.md` cites the **CLAUDE.md §4-FORBIDDEN** Balanis
`73 + j42 Ω` dipole reference (must be NEC-4 `87 + j41 Ω` only). Make
`yee validate` tell the truth.

## 2. Approach (walking-skeleton-first; single lane)

The FEM drivers already return a `FemEig00NValidationResult` carrying the
same `status: CaseStatus` / `notes: String` / `wall_time_seconds: f64`
triple as `CaseResult` — so registration is a fold, not new physics.

- **`crates/yee-validation/src/lib.rs`:** add thin `run_fem_eig_00N()`
  wrappers that call the existing public drivers and fold their result into
  `CaseResult { id: "fem-eig-00N", … }`; push them into the `run_all` vec.
  **Wall-time discipline:** the fast gates (e.g. fem-eig-001 ≈7 s) go in the
  default `run_all`; for any driver whose mesh (`FEM_EIG_00N_N{X,Y,Z}`
  constants) makes it multi-minute, follow the `mom-001` precedent (it is
  already wall-time-gated in the suite) — register it but keep the default
  path fast (gate/feature/`#[ignore]`-class as the existing slow cases do),
  documenting any deferral. Do NOT bloat the default `yee validate all`.
- **`crates/yee-cli/src/main.rs`:** add a `Fem` variant to `ValidateTarget`
  + a `fem-*` arm in `case_matches_target`, so `yee validate fem` runs the
  FEM cases.
- **`validation/README.md`:** de-stale — replace the forbidden `73 + j42`
  dipole reference with NEC-4 `87 + j41` (per CLAUDE.md §4 / ADR-0005), and
  reconcile the case list with what `run_all` actually registers (drop the
  fictional `validation/cases/phase-{0..4}/` tree claim if present).

## 3. Definition of done

DoD-1. The public FEM-eig gates are registered in `Report::run_all` (the
fast ones in the default path; slow ones registered with documented
wall-time discipline, not bloating default CI). `yee validate` gains a
`fem` target.
DoD-2. New aggregator test: `Report::run_all().cases` contains
`fem-eig-001` with `CaseStatus::Passed` (extend the existing
`tests/integration.rs`, which uses `.find()` — non-breaking).
DoD-3. New CLI test (assert_cmd): `yee validate fem --json` exits 0 + the
JSON contains a `fem-eig-001` entry.
DoD-4. `validation/README.md` no longer cites `73 + j42` (NEC-4 `87 + j41`
only) + its case list matches the registry.
DoD-5. fmt + clippy `-D warnings` clean; `cargo test -p yee-validation -p yee-cli`
green; no new dependency; the existing 6 registered cases + their tests
unchanged in behaviour.

## 4. NON-NEGOTIABLE

- Lane: `crates/yee-validation/src/**`, `crates/yee-cli/src/**`,
  `validation/README.md`, + the relevant `tests/`. Do NOT touch the FEM
  solver (`crates/yee-fem/**`) or its drivers' logic — the drivers are
  consumed AS-IS (they are public + passing). Do NOT touch `crates/yee-fdtd/**`.
- The FDTD `Skipped` stubs (cpml/ntff/dispersive) + the fdtd-201/.x cavity
  gates need NEW public drivers in the yee-fdtd lane — **OUT OF SCOPE**
  here (a follow-on slice); leave them.
- No new `Cargo.toml` dependency (the FEM drivers are already compiled into
  yee-validation).

## 5. References

* `crates/yee-validation/src/lib.rs` — `run_all` registry (~L130), the
  public FEM drivers (`run_fem_eig_001_rectangular_cavity` ~L1481, …006
  ~L3875) + their `FemEig00NValidationResult` structs, the hardcoded FDTD
  `Skipped` stubs (~L1264).
* `crates/yee-cli/src/main.rs` — `ValidateTarget` (~L316), `case_matches_target`
  (~L646), `run_validate`.
* `crates/yee-validation/tests/integration.rs` (the `.find()`-based,
  non-breaking aggregator test to extend).
* `docs/src/decisions/0008-validation-aggregator-json-and-png.md` (the
  aggregator-shape ADR — read before extending). CLAUDE.md §4 + ADR-0005
  (NEC-4-only for mom-001).
