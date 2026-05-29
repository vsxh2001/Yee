# Phase 1.plotting.3 — CLI Smith Chart Multi-trace Design

**Date:** 2026-05-25  
**Phase:** 1.plotting.3  
**ADR:** ADR-0069  
**Status:** Accepted

---

## 1. Context

ADR-0068 (Phase 1.gui.5 / 1.plotting.2) shipped:
- `yee-plotters`: `SmithTrace` + `plot_smith_chart_multi` (multi-trace Smith chart
  with constant-R/X arc overlays).
- `yee-gui`: multi-trace Smith tab with "Show all entries" checkbox.

ADR-0068 explicitly deferred `yee-cli --smith --all` multi-trace:

> **Deferred** (still out of scope): arc impedance labels, VSWR circles,
> admittance chart, **`yee-cli --smith --all` multi-trace**.

The current `yee plot` CLI:
- `--format smith` **single-trace** works (calls `plot_smith_chart` / the
  back-compat wrapper around `plot_smith_chart_multi`).
- `--format smith` + `--entry`/`--all` **is rejected** with a hard error:

  ```
  the `smith` and `both` plot kinds are not supported with
  --entry / --all; use --format db or --format phase
  ```

The rejection exists because `plot_smith_chart_multi` did not exist when multi-trace
was first added. ADR-0068 landed `plot_smith_chart_multi`, making the wiring trivial.

---

## 2. Decision

Wire `plot_smith_chart_multi` into the CLI multi-trace path.

### Scope

- **`crates/yee-cli/src/plot.rs`** — remove the Smith/Both rejection in
  `run_multi_trace`; handle all four `PlotKind` variants.
- **`crates/yee-cli/tests/cli_plot_touchstone.rs`** — update
  `plot_entry_with_smith_errors_cleanly` (now a success test) and add three
  new tests.
- **No change** to `yee-plotters` (already supports multi-trace Smith).
- **No new dependency**.

### Behaviour

| Invocation | Result |
|---|---|
| `yee plot in.s2p --format smith --entry 11` | Single Smith trace (label S11) |
| `yee plot in.s2p --format smith --entry 11 --entry 21` | Two overlaid Smith traces |
| `yee plot in.s2p --format smith --all` | All N² entries overlaid |
| `yee plot in.s2p --format both --all` | Two files: `out-db.png` + `out-smith.png`, each multi-trace |

The `--port` single-trace path is **byte-unchanged**.

---

## 3. Implementation plan

See the companion plan file `2026-05-25-phase-1-plotting-3-smith-cli-multi-trace.md`.

---

## 4. Validation

- `cargo clippy --workspace --all-targets -- -D warnings` exits 0.
- `cargo fmt --check --all` exits 0.
- `cargo test --workspace` exits 0 (all tests including the four new/updated CLI
  integration tests).

---

## 5. References

- ADR-0068 `docs/src/decisions/0068-smith-chart-constant-rx-arcs.md`
- ADR-0063 `docs/src/decisions/0063-multi-trace-sparam-plotting.md`
- `crates/yee-plotters/src/lib.rs` — `SmithTrace`, `plot_smith_chart_multi`
- `crates/yee-cli/src/plot.rs` — `run_multi_trace`, `PlotKind`
