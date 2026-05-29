# Filter Phase F0.1 — register synthesis gates in aggregator — Design Spec

**Phase:** F0.1 · **ADR:** ADR-0085 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
Surface the F0 synthesis gates (`synth-001`, `synth-002`, `filt-001`) in the
`yee-validation` aggregator + `yee validate` CLI, reusing the single-source
`case_registry()` so they appear in `yee validate --list[ --json]` and a new
`yee validate synth` target. Fast, pure-math, no EM.

## yee-validation changes
- Add `Solver::Synth` to the `Solver` enum (serialize → `"Synth"`).
- Add driver wrappers (each returns `CaseResult`, `policy=Run`):
  - `run_synth_001()` — recompute Chebyshev 0.5 dB N=5 g-values via
    `yee_synth::prototype`; `Passed` iff all within 1e-3 of the published
    `[1.7058, 1.2296, 2.5408, 1.2296, 1.7058]` (g6=1.0); else `Failed` with the
    max deviation in `notes`.
  - `run_synth_002()` — `yee_synth::coupling_design` for Chebyshev 0.5 dB N=3,
    FBW=0.10; `Passed` iff `k12==k23` (symmetry, ≤1e-9) and `k`,`Qe` match the
    closed form recomputed in the driver (≤1e-9).
  - `run_filt_001()` — `yee_filter::synthesize` a satisfiable Chebyshev bandpass,
    sweep, `check_mask`; `Passed` iff the mask report passes.
- Register all three in `case_registry()` (append after the fem-eig block), each
  paired with a `CaseDescriptor { id, solver: Solver::Synth, description, policy:
  Run }`. `run_all()`/`list_cases()` pick them up automatically.
- Deps: add `yee-synth`, `yee-filter` to `yee-validation/Cargo.toml`.

## yee-cli changes
- `ValidateTarget::Synth`; `case_matches_target`: `Synth => id.starts_with("synth-")
  || id.starts_with("filt-")`.
- `run_validate_list` Solver match: add `Solver::Synth => "Synth"`.
- Update the `//!` synopsis + `Validate` doc to list the `synth` target.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-validation -p yee-cli --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-validation -p yee-cli` exit 0 (fast; do NOT pull mom-001 /
   fem-eig-003 — use filtered/targeted tests; the 3 synth drivers are ms-scale).
4. yee-validation test `synth_filt_cases_registered_and_pass`: `run_synth_001()`,
   `run_synth_002()`, `run_filt_001()` each return `CaseStatus::Passed`; and
   `list_cases()` contains `synth-001`, `synth-002`, `filt-001` with
   `solver == Solver::Synth`.
5. yee-cli test `yee_validate_synth_list_runs`: `yee validate synth --list`
   exits 0, stdout contains `synth-001`, `synth-002`, `filt-001` and `Synth`.
6. `yee validate --list` now shows the three synth rows; `--list --json` includes
   `"Synth"`.

## Out of scope
New synthesis math; EM; layout; any change to the existing 19 cases' behavior.
