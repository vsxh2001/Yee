# ADR-0104: Filter — register coupled/dim/gerber gates in the validation aggregator

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0085 (F0.1 — synth gates in the aggregator, `Solver::Synth`),
ADR-0094 (coupled-line model), ADR-0097 (dimensional synthesis), ADR-0100
(Gerber), `FILTER-DESIGN-ROADMAP.md`

---

## Context

F0.1 (ADR-0085) registered the synthesis gates (`synth-001/002`, `filt-001`) in
the `yee-validation` aggregator under `Solver::Synth`, so they appear in
`yee validate --list` and `Report::run_all`. Since then the filter track has
grown crate-test gates that are NOT in the aggregator: `coupled-001` (ADR-0094
coupled-line model vs Steer), `dim-001` (ADR-0097 dimensional-synthesis inversion
round-trip), and `gerber-001` (ADR-0100 Gerber structure). The aggregator is the
project's single validation dashboard; the filter-design pipeline should be fully
visible there.

## Decision

Register three more filter cases in `crates/yee-validation/src/lib.rs`
`case_registry()`, under the existing `Solver::Synth` (no new `Solver` variant —
avoids exhaustive-match churn; the case IDs distinguish them), each a
`run_X() -> CaseResult` mirroring `run_synth_001` / `run_filt_001`:

- **`coupled-001`** — `yee_layout::coupled_microstrip` for the Steer Example 5.6.1
  point (εr=10, W/h=1, s/h=0.5) reproduces `Z0e≈59`, `Z0o≈37` within tolerance
  (measured Z0e/Z0o vs reference).
- **`dim-001`** — synthesize the Chebyshev N=5 fixture, `dimension_edge_coupled`
  on FR-4, assert each gap re-evaluates to its `target_k` within < 1 % (measured
  max relative error).
- **`gerber-001`** — `layout_to_gerber` of a small layout has the RS-274X header,
  one `G36/G37` region per polygon, and `M02*` (measured region count).

`yee-validation` gains `yee-layout` + `yee-export` deps (it already has
`yee-synth`/`yee-filter`). `list_cases()` derives from `case_registry()`, so it
updates automatically; the registry↔list invariant test stays green.

## Consequences

**Ships:** `yee validate --list` (and `run_all`) now include `coupled-001`,
`dim-001`, `gerber-001` under `Solver::Synth`. One command validates the whole
filter-design pipeline (synthesis → coupled-line model → dimensions → Gerber).

**Gate:** `cargo test -p yee-validation` passes (incl. the registry↔list
invariant + any case-count assertions, updated for the three new IDs); the three
new cases report `pass`. No FDTD; the cases are pure-math/text (fast).

**Not in scope:** a new `Solver` variant; FDTD/EM cases (F1.1b.1+); changing the
underlying crate tests (these aggregator cases re-exercise the same checks, they
do not replace the crate gates).

---

## References
- ADR-0085 (the `Solver::Synth` / `case_registry` pattern this extends).
- `docs/superpowers/specs/2026-05-30-filter-validation-aggregator-gates-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-validation-aggregator-gates.md`.
