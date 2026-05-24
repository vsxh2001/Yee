# GUI multi-trace S-parameter overlay (egui_plot)

**Status:** Draft
**Owner:** TBD
**Phase:** frontend (yee-gui) — follow-on to ADR-0063 (static multi-trace)
**Type:** additive product feature (GUI capability gain)

## 1. Goal

The static plotters + CLI gained multi-trace / off-diagonal S-parameter
overlay (ADR-0063). The desktop GUI still shows only a **single** trace —
`plots.rs` renders one `Line::new("|S11| (dB)", …)` (`plots.rs:70`) and
the panel is wired to one diagonal entry. So a user inspecting a 2-port
device in the desktop app cannot see S21 or overlay S11+S21. Bring the
GUI to parity: overlay multiple labelled S-parameter traces in the
egui_plot dB panel with a legend.

## 2. Approach

Keep the **data selection / series-building logic separate from the
egui rendering** so the logic is unit-testable (egui rendering needs a
live context and is only smoke-tested).

- **`yee-gui` (`plots.rs` + `app.rs`):**
  - A pure helper (testable, no egui): given the loaded `yee_io::touchstone::File`
    + a selection (set of `(row,col)` entries, or "all"), build a
    `Vec<SparamSeries { label: String, points: Vec<[f64;2]> }>` (freq-GHz
    vs dB), extracting from the row-major `data[k]` exactly like the CLI
    (`flat_idx = r*n + c`). Unit-test this (entry selection, labels `S<ij>`,
    bounds, dB conversion via the existing `db_clamped`-equivalent).
  - The egui_plot dB panel overlays each series as a `Line::new(label, pts)`
    with `egui_plot`'s `Legend` (a colour per series). Single-trace stays
    the default; a small UI control (checkboxes / a "show all entries"
    toggle for multi-port files) selects the overlay set.
- Reuse the magnitude-dB conversion already in `plots.rs`; do NOT
  duplicate it.

## 3. Definition of done

DoD-1. A pure `yee-gui` series-building helper turns a loaded file + an
entry selection into labelled (freq, dB) series; **unit-tested** (entry
selection, `S<ij>` labels, row-major extraction, out-of-range handling).
DoD-2. The egui_plot dB panel overlays the selected traces with a legend
+ per-trace colour; a UI control selects which entries (default = the
existing single-trace behaviour for back-compat / 1-port files).
DoD-3. `cargo build -p yee-gui` + `cargo test -p yee-gui` green; fmt +
clippy `-D warnings` clean; all new public items documented.
DoD-4. No change outside `crates/yee-gui/**`; no new dependency (egui_plot
`Legend`/`Line` already available); the Smith + 3D-viewport panels
unchanged.

## 4. NON-NEGOTIABLE

- Lane: `crates/yee-gui/**` only. Do NOT touch `yee-plotters` (static),
  the solver crates, `yee-io`, or the CLI. No new `Cargo.toml` dependency.
- Do not break the existing single-trace dB panel / Smith panel / viewport.
- egui rendering is hard to unit-test — that is expected; the DoD's
  testable surface is the **series-building logic**, not the egui draw
  calls. Do NOT invent a brittle headless-render test.

## 5. References

* `crates/yee-gui/src/plots.rs` (`plot_s11_db` egui_plot fn `:70`, Smith
  `:86`; the dB-conversion to reuse), `crates/yee-gui/src/app.rs` (the
  loaded `File` + the plot tabs `:47,:231`).
* ADR-0063 + `crates/yee-cli/src/plot.rs` (the CLI multi-entry extraction
  to mirror: `flat_idx = r*n+c`, `S<ij>` labels), `crates/yee-plotters/src/lib.rs`
  (`SparamTrace` — the static analogue).
* egui_plot `Legend` / `Line` API (egui_plot 0.35, already a dep).
