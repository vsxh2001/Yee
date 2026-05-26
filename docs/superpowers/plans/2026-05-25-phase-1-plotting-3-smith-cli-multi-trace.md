# Phase 1.plotting.3 — CLI Smith Chart Multi-trace Implementation Plan

**Spec:** `docs/superpowers/specs/2026-05-25-phase-1-plotting-3-smith-cli-multi-trace-design.md`  
**ADR:** ADR-0069  
**Lane:** `crates/yee-cli/**`

---

## Steps

### S1 — Update `run_multi_trace` in `crates/yee-cli/src/plot.rs`

1. Remove the `PlotKind::Smith | PlotKind::Both` rejection arm.
2. Refactor to extract `SparamTrace` first (for Db/Phase/Both dB output), then
   derive `SmithTrace` from the same data (for Smith/Both smith output):
   ```rust
   // Convert SparamTrace → SmithTrace field-for-field.
   let smith_traces: Vec<yee_plotters::SmithTrace> = sparam_traces
       .iter()
       .map(|t| yee_plotters::SmithTrace {
           label: t.label.clone(),
           values: t.values.clone(),
       })
       .collect();
   ```
3. Add `PlotKind::Smith` arm:
   ```rust
   PlotKind::Smith => {
       yee_plotters::plot_smith_chart_multi(&smith_traces, &args.output, &config)
           .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
       eprintln!("yee plot: wrote {}", args.output.display());
   }
   ```
4. Add `PlotKind::Both` arm (multi-trace — mirrors single-trace `Both`):
   ```rust
   PlotKind::Both => {
       let db_path = suffixed_path(&args.output, "-db");
       let smith_path = suffixed_path(&args.output, "-smith");
       yee_plotters::plot_sparams_db(&file.freq_hz, &sparam_traces, &db_path, &config)
           .map_err(|e| anyhow::anyhow!("plot (db): {e}"))?;
       yee_plotters::plot_smith_chart_multi(&smith_traces, &smith_path, &config)
           .map_err(|e| anyhow::anyhow!("plot (smith): {e}"))?;
       eprintln!("yee plot: wrote {} and {}", db_path.display(), smith_path.display());
   }
   ```
5. Update the doc-comment on `run_multi_trace` and the module-level doc to remove
   the now-incorrect "smith and both are not supported" note.

### S2 — Update CLI integration tests in `cli_plot_touchstone.rs`

1. Rename `plot_entry_with_smith_errors_cleanly` →
   `plot_entry_with_smith_multi_trace_produces_png` and flip the assertion to
   **success** + file-size check.
2. Add `plot_all_entries_with_smith_produces_png` (uses `--all --format smith` on
   a 2-port `.s2p`; asserts success + size > 1 KB).
3. Add `plot_multi_trace_both_emits_db_and_smith` (uses `--entry 11 --entry 21
   --format both`; asserts two output files each > 1 KB).

### S3 — Verify

```bash
cargo clippy --workspace --all-targets -- -D warnings   # expect exit 0
cargo fmt --check --all                                  # expect exit 0
cargo test --workspace                                   # expect exit 0
```

---

## Definition of Done

- `cargo test --workspace` exits 0 with no skipped or failed tests.
- `cargo clippy --workspace --all-targets -- -D warnings` exits 0.
- `cargo fmt --check --all` exits 0.
- The four integration tests (1 updated + 3 new) all pass.
- `yee plot in.s2p --format smith --all` and `--entry 11 --entry 21` are
  accepted (no error, non-empty PNG produced).
- The `--port` single-trace path is byte-unchanged (existing single-trace tests
  still pass).
