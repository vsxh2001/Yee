# Filter — register coupled/dim/gerber gates in the aggregator — Plan

**Spec:** `2026-05-30-filter-validation-aggregator-gates-design.md` · **ADR:** ADR-0104

## Lane
`crates/yee-validation/**` ONLY (`Cargo.toml`, `src/lib.rs`). Do NOT edit the
filter/layout/export crates or their tests — consume their public API. Out of
lane → finding.

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-validation-aggregator-gates`.

## Pattern files
- `crates/yee-validation/src/lib.rs`:
  - `run_synth_001` (~line 2712) and `run_filt_001` (~line 2876) — the exact
    `CaseResult` construction + pass/measured/reference idiom to MIRROR.
  - `case_registry()` (~line 292) — where the `Solver::Synth` entries are
    registered (`(CaseDescriptor { id, solver, .. }, run_X as fn() -> CaseResult)`);
    add the three new entries here, same shape.
  - the `CaseResult` + `CaseDescriptor` + `Solver` definitions (~line 100) — read
    the field set before constructing.
  - the registry↔`list_cases` invariant test (~line 5852) and any case-count test
    — make sure they pass (update a hard-coded expected count if present).
- Reference data/logic (READ, do not edit): `crates/yee-layout/tests/
  coupled_001_vs_published.rs` (exact Steer numbers + tol), `crates/yee-filter/
  tests/dim_001_inversion_roundtrip.rs` (the round-trip check), the F0 cheb_bpf
  fixture spec values.

## Steps
1. `Cargo.toml`: add `yee-layout` + `yee-export` workspace deps.
2. `src/lib.rs`: add `run_coupled_001`, `run_dim_001`, `run_gerber_001` (mirror
   `run_synth_001`/`run_filt_001`); register all three in `case_registry()` under
   `Solver::Synth` with unique ids.
3. Ensure the registry↔list invariant + any count test pass (adjust expected
   count if hard-coded).

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-validation --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-validation --jobs 2
```
Pure math/text — fast. Do NOT run `cargo test --workspace`, FDTD, mom-001 (other
crates' heavy cases are not exercised by `-p yee-validation`'s own tests).

## Escape hatch
Blocked > 15 min — the `CaseResult`/`CaseDescriptor` shape doesn't fit a
pass/measured/reference filter check, OR the registry has invariants (ordering,
id-prefix conventions, count tests) that fight three additions → STOP and surface
the struct shape + the specific invariant. Do NOT add a new `Solver` variant
(reuse `Synth`). Do NOT weaken any existing case or invariant test. Do NOT edit
other crates.

## Done when
DoD 1–4 pass; `git diff --stat <base>..HEAD` = only `crates/yee-validation/**`
(+ `Cargo.lock`) + the 3 committed docs; the three new ids appear in `list_cases()`.
