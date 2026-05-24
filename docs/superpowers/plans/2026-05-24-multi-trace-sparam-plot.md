# Multi-trace S-parameter plotting — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-multi-trace-sparam-plot-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-plotters/**`, `crates/yee-cli/src/plot.rs` (+ its
args + `crates/yee-cli/tests/**` if a CLI test fits there). NOTHING else.
**Out of lane** (findings, not fixes): the GUI (`crates/yee-gui/**` —
egui_plot overlay is a follow-on), the solver crates, `crates/yee-io/**`
(Touchstone read-only). No `Cargo.toml` dependency. Do not change the
existing single-trace plotter signatures or the default `yee plot --port`
behaviour.

## Step ladder

### S1 — read patterns
Read `crates/yee-plotters/src/lib.rs` (`plot_s11_db` body, `PlotConfig`,
`db_clamped`, the existing tests at the file tail) and
`crates/yee-cli/src/plot.rs` (the diagonal extraction + arg parsing) to
match house style + reuse the axis/colour machinery.

### S2 — plotters multi-trace fn
Add `SparamTrace` + `plot_sparams_db(freq_hz, traces, output, &PlotConfig)`:
overlay each trace as a labelled dB line with a legend + distinct colours
from a small fixed palette. Keep the single-trace fns intact. (Optional
`plot_sparams_phase` if cheap.) Unit test: ≥2 traces → output file exists
+ non-trivial size (+ a content assertion where the backend allows),
mirroring the crate's existing plot tests.

### S3 — CLI wiring
Generalize `yee plot`: add repeated `--entry <ij>` (and/or `--all`) to
select S-matrix entries; extract each from row-major `data[k]`; label
`S<ij>`; call `plot_sparams_db`. Keep `--port` (diagonal) as the
unchanged default. Bounds-check entries vs `n_ports` with a clean error.
A CLI test covers a 2-port file → S11+S21 overlay output.

### S4 — docs
Doc all new public items; a one-line `yee plot` help/usage note for the
new flag.

## Verification (run in worktree; all exit 0)
```
cargo fmt --check --all
cargo clippy -p yee-plotters -p yee-cli --all-targets -- -D warnings
cargo test -p yee-plotters
cargo test -p yee-cli
git diff --stat -- crates/yee-gui crates/yee-mom crates/yee-fdtd crates/yee-io '**/Cargo.toml'   # MUST be empty
```

## Escape-hatch
- If `plotters` legend/multi-series rendering needs an unexpected API or a
  new dependency, STOP + surface (do not add a dep). If the CLI arg change
  would break the default `--port` path, STOP + surface.
- Blocked >15 min → surface the specific blocker. Run synchronously; no
  Monitor/ScheduleWakeup; no sub-agents.

## Out-of-scope (findings, not fixes)
* GUI (egui_plot) multi-trace overlay in the S-param panel — a follow-on.
* Smith-chart constant-R/X arc family (a separate plotters increment).
