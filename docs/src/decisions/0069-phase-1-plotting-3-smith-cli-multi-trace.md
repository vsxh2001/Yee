# ADR-0069 — Phase 1.plotting.3 CLI Smith Chart Multi-trace

**Status:** Accepted  
**Date:** 2026-05-25  
**Context Phase:** 1.plotting.3

## Context

ADR-0068 (Phase 1.gui.5 / 1.plotting.2) shipped `yee-plotters::SmithTrace` +
`plot_smith_chart_multi` and the GUI Smith multi-trace tab, but explicitly
deferred the `yee-cli` multi-trace Smith path:

> **Deferred** (still out of scope): arc impedance labels, VSWR circles,
> admittance chart, `yee-cli --smith --all` multi-trace.

At the time of ADR-0063 (multi-trace S-parameter plotting), `plot_smith_chart_multi`
did not yet exist, so the CLI `run_multi_trace` function rejected
`PlotKind::Smith | PlotKind::Both` with a hard error. After ADR-0068 the
rejection is a dead letter — the plotter already supports what the CLI refuses
to invoke.

## Decision

Wire `plot_smith_chart_multi` into `crates/yee-cli`'s `run_multi_trace`:

* Remove the `PlotKind::Smith | PlotKind::Both` rejection.
* Build `SmithTrace` by converting the already-extracted `SparamTrace` vector
  field-for-field (same label, same values).
* Handle all four `PlotKind` variants in the multi-trace match:
  - `Db` → existing `plot_sparams_db` (unchanged).
  - `Phase` → existing `plot_sparams_phase` (unchanged).
  - `Smith` → `plot_smith_chart_multi`.
  - `Both` → `plot_sparams_db` to `out-db.<ext>` + `plot_smith_chart_multi` to
    `out-smith.<ext>` (mirrors single-trace `Both`).
* Update the now-incorrect doc comments.
* No new dependency; no change to `yee-plotters`.

## Consequences

* `yee plot in.s2p --format smith --entry 11 --entry 21` overlays two Smith
  traces with constant-R/X arc overlays (the Phase 1.plotting.2 arc grid
  comes for free via `plot_smith_chart_multi`).
* `yee plot in.s2p --format smith --all` overlays every entry.
* `yee plot in.s2p --format both --entry 11 --entry 21` emits two files:
  `out-db.png` (dB magnitude overlay) and `out-smith.png` (Smith overlay).
* The single-trace `--port` path is byte-unchanged.
* The integration test `plot_entry_with_smith_errors_cleanly` (which formerly
  verified the rejection) is updated to verify the new success behaviour.

## References

* Spec `docs/superpowers/specs/2026-05-25-phase-1-plotting-3-smith-cli-multi-trace-design.md`
* Plan `docs/superpowers/plans/2026-05-25-phase-1-plotting-3-smith-cli-multi-trace.md`
* ADR-0068 `docs/src/decisions/0068-smith-chart-constant-rx-arcs.md`
* ADR-0063 `docs/src/decisions/0063-multi-trace-sparam-plotting.md`
