# Filter — register coupled/dim/gerber gates in the aggregator — Design Spec

**ADR:** ADR-0104 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Add `coupled-001`, `dim-001`, `gerber-001` to `yee-validation`'s `case_registry()`
(under `Solver::Synth`), so the filter-design pipeline is fully visible in
`yee validate --list` / `Report::run_all`. Pure-math/text re-exercise of shipped
checks; no FDTD; no new `Solver` variant.

## Changes (`crates/yee-validation/**` ONLY)
- `Cargo.toml`: add `yee-layout = { workspace = true }` and `yee-export = {
  workspace = true }` (already has yee-synth + yee-filter).
- `src/lib.rs`: add three `pub fn run_X() -> CaseResult` mirroring `run_synth_001`
  (~line 2712) and `run_filt_001` (~line 2876) — same `CaseResult` construction,
  same pass/measured/reference idiom. Register each in `case_registry()` (~line
  292, beside the existing `Solver::Synth` entries) with a `CaseDescriptor`
  (unique `id` = `"coupled-001"` / `"dim-001"` / `"gerber-001"`, `solver:
  Solver::Synth`, and whatever description/reference fields the descriptor
  carries — copy the shape from the synth entries).

### The three cases
- **`run_coupled_001`** — `yee_layout::coupled_microstrip(w, s, h, eps_r)` for the
  Steer Example 5.6.1 point. Read the EXACT reference numbers + tolerance from the
  shipped `crates/yee-layout/tests/coupled_001_vs_published.rs` (do not invent);
  pass if `Z0e`/`Z0o` within tol. Measured = computed `Z0e`/`Z0o` (or the worse
  relative error); reference = the published values.
- **`run_dim_001`** — build the committed Chebyshev 0.5 dB N=5 spec (the
  `cheb_bpf` fixture: f0=2e9, fbw=0.10, z0=50) via `yee_filter::synthesize`,
  `dimension_edge_coupled` on FR-4 (εr=4.4, h=1.6e-3), then for each section assert
  `coupling_coefficient(coupled_microstrip(width, gap_i, h, εr))` reproduces
  `target_k[i]` within < 1 %. Mirror `crates/yee-filter/tests/dim_001_inversion_roundtrip.rs`.
  Measured = max relative error; reference = `< 0.01`.
- **`run_gerber_001`** — `yee_export::layout_to_gerber` of a small layout (e.g. via
  `dimension_edge_coupled_layout`, or `Polygon::rect` traces). Pass if the output
  has `%FSLAX46Y46*%`, one `G36*`/`G37*` per polygon, and `M02*`. Measured =
  region count; reference = polygon count.

Match whatever `CaseResult` requires (pass bool, measured/reference strings or
f64s, units, notes). Read the struct + a couple of existing `run_*` first.

## DoD (machine-checkable; NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-validation --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-validation` exit 0 — INCLUDING the existing registry↔
   `list_cases` invariant test (~line 5852) and any case-count assertion (update
   the expected count if the test hard-codes one). The three new cases report
   `pass`.
4. `yee validate --list` (or the `list_cases()` API in a test) includes
   `coupled-001`, `dim-001`, `gerber-001`.

## Out of scope
A new `Solver` variant; FDTD/EM cases; modifying the underlying crate tests; the
CLI `validate` subcommand wiring if it already routes `Solver::Synth` (it does —
F0.1). Keep the cases fast (pure math/text).
