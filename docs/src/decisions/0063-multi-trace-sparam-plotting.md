# ADR-0063 — multi-trace S-parameter plotting (off-diagonal overlay)

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** frontend (yee-plotters + yee-cli)

## Context

After two MoM cross-section / microstrip-port tracks (ADRs 0050–0061) and
one FDTD validation rotation (ADR-0062, fdtd-201), the rotation survey's
second candidate is a clean frontend gap: the static plotters expose only
**single-trace** entry points and `yee plot` extracts only the **diagonal**
`S[port][port]`, so there is no way to plot an off-diagonal entry (S21) or
overlay multiple S-parameters — the most common multi-port view.

## Decision

Add a **multi-trace overlay** magnitude-dB plotter (`plot_sparams_db`,
labelled traces + legend) to `yee-plotters`, and wire `yee plot` to
select/overlay arbitrary S-matrix entries (repeated `--entry <ij>` / `--all`),
keeping the existing single-trace fns + the default `--port` diagonal
behaviour unchanged. Bounded first slice = **plotters + CLI only**; the GUI
(egui_plot) overlay is a documented follow-on.

## Rationale

(1) **Clean, certain, dispatchable win** (per the standing value×
dispatchability rule, after a clean FDTD win): additive, self-contained,
disjoint frontend lane, no new dependency, single-pass. It is the
opposite of a grind-risk quagmire.

(2) **Real product value** — multi-port S-parameter overlay (S11+S21+…)
is the standard view a microwave user expects; its absence is a
conspicuous usability gap.

(3) **Breadth** — a different subsystem (frontend) again, continuing to
spread coverage off the MoM-port theme.

(4) **Not a solver gate** — CLAUDE.md §4 (published-benchmark) does not
apply; the appropriate validation is plotter unit tests (the crate's
existing output-non-empty / content-assertion style) + a CLI smoke. This
is acknowledged, not a tolerance dodge.

## Consequences

* `yee-plotters` gains `plot_sparams_db` (+ optional phase) + a
  `SparamTrace`; the single-trace fns are untouched (CLI + tests depend
  on them).
* `yee plot` gains multi-entry / off-diagonal selection; `--port` default
  unchanged.
* GUI overlay + the Smith constant-R/X arc family remain documented
  follow-ons.
* No solver-crate change; no new dependency.

## References

* `crates/yee-plotters/src/lib.rs`, `crates/yee-cli/src/plot.rs`,
  `crates/yee-io/src/touchstone.rs` (the `File`/`SMatrix` data source).
* Spec + plan (2026-05-24). The rotation survey (GUI/plotters gap).
* ADR-0062 (the preceding FDTD rotation).
