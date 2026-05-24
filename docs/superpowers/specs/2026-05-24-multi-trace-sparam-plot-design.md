# Multi-trace S-parameter plotting (off-diagonal S21 overlay)

**Status:** Draft
**Owner:** TBD
**Phase:** frontend (yee-plotters + yee-cli)
**Type:** additive product feature (clean breadth increment)

## 1. Goal

Today the static plotters expose only **single-trace** magnitude/phase/
Smith entry points (`plot_s11_db`, `plot_s11_phase`, `plot_smith_chart`,
each taking one `&[Complex64]` trace), and `yee plot` extracts only the
**diagonal** `S[port][port]` (`yee-cli/src/plot.rs:51`). There is no way
to plot an off-diagonal entry (S21) or overlay multiple traces (S11+S21+…)
— the most common S-parameter view for a multi-port device. Add a
**multi-trace overlay** magnitude-dB (and phase) plot + wire the CLI to
select/overlay arbitrary S-matrix entries.

This is a frontend/usability feature, not a solver gate (CLAUDE.md §4
does not apply); its validation is unit tests on the plotter (the crate's
existing SVG/PNG non-empty + content-assertion style) + a CLI smoke.

## 2. Approach

### Part A — `yee-plotters` multi-trace entry point
Add `plot_sparams_db(freq_hz: &[f64], traces: &[SparamTrace], output, &PlotConfig)`
where `SparamTrace { label: String, values: &[Complex64] }` (or
`(String, Vec<Complex64>)`). Overlay each trace as a labelled magnitude-dB
line with a **legend** + a per-trace colour from a small palette; reuse
the existing axis/`db_clamped` machinery. Optionally a sibling
`plot_sparams_phase`. The existing single-trace fns stay (call them or
keep them as-is; do NOT break their signatures — the CLI + tests use them).

### Part B — `yee-cli` wiring
Extend `yee plot` to select multiple S-matrix entries (e.g. repeated
`--entry <ij>` like `--entry 11 --entry 21`, and/or `--all` for every
entry of a multi-port file), extract each from the row-major `data[k]`
(`n_ports × n_ports`), label them `S<ij>`, and call `plot_sparams_db`.
The existing single-`--port` diagonal behaviour stays the default
(back-compat). Bounds-check requested entries against `n_ports`.

## 3. Definition of done

DoD-1. `plot_sparams_db` (+ optional phase) overlays ≥2 labelled traces
with a legend into a non-empty SVG/PNG; unit-tested in the crate's
existing style (output exists + is non-trivial; a content/string check
where the backend allows). Single-trace fns unchanged.
DoD-2. `yee plot` can plot an off-diagonal entry (S21) and overlay
multiple entries; the default single-`--port` path is unchanged;
out-of-range entries error cleanly. A CLI test covers the multi-entry path.
DoD-3. fmt + clippy `-D warnings` clean; `cargo test -p yee-plotters -p yee-cli`
green; all public items documented (`#![warn(missing_docs)]`).
DoD-4. No solver-crate change; no new dependency (use the existing
plotters palette/colours).

## 4. NON-NEGOTIABLE

- Lane: `crates/yee-plotters/**` + `crates/yee-cli/src/plot.rs` (+ its
  args/tests). Do **not** touch the GUI in this slice (egui_plot overlay
  is a documented follow-on), the solver crates, or `yee-io`'s Touchstone
  parser. Do not break the existing single-trace plotter signatures or
  the default `yee plot --port` behaviour.
- No new `Cargo.toml` dependency.

## 5. References

* `crates/yee-plotters/src/lib.rs` (`plot_s11_db` ~117, `plot_s11_phase`,
  `plot_smith_chart`, `PlotConfig`, `db_clamped`) — the pattern + the
  machinery to reuse.
* `crates/yee-cli/src/plot.rs` (the diagonal-only extraction at :51 to
  generalize).
* `crates/yee-io/src/touchstone.rs` (`File { n_ports, freq_hz, data }`,
  row-major `SMatrix`) — the data source (read-only).
* Rotation survey (the GUI/plotters gap finding) + ADR-0063.
