# GUI multi-trace S-parameter overlay — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-gui-multi-trace-overlay-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-gui/**` ONLY.
**Out of lane** (findings, not fixes): `crates/yee-plotters/**` (the static
analogue — read for the pattern, do NOT edit), the solver crates,
`crates/yee-io/**` (Touchstone read-only), `crates/yee-cli/**`. No
`Cargo.toml` dependency.

## Step ladder

### S1 — read patterns
Read `crates/yee-gui/src/plots.rs` (the egui_plot `plot_s11_db` fn + the
dB conversion) + `app.rs` (how the loaded `File` reaches the plot tabs +
any existing entry/port selection state). Read `crates/yee-cli/src/plot.rs`
for the entry-extraction idiom to mirror (`flat_idx = r*n+c`, `S<ij>`).

### S2 — pure series-building helper (the tested surface)
Add a pure fn (no egui) in `plots.rs` (or a small submodule): given
`&yee_io::touchstone::File` + a selection (e.g. `enum Selection { Diagonal(usize), Entries(Vec<(usize,usize)>), All }`),
return `Vec<SparamSeries { label, points: Vec<[f64;2]> }>` (freq-GHz vs dB).
Reuse the existing dB conversion. **Unit-test it**: 2-port `All` → 4 series
labelled S11/S21/S12/S22; a specific entry selection; out-of-range → empty
or error; dB values match the conversion.

### S3 — egui_plot overlay + UI control
The dB panel iterates the series, drawing `Line::new(label, pts)` each,
with `Plot::new(...).legend(egui_plot::Legend::default())`. Add a UI
control (checkbox per entry for small n, or a "show all entries" toggle)
that sets the `Selection`; default to the existing single-trace behaviour
(so 1-port files + the current UX are unchanged). Colours: let egui_plot
auto-assign, or a small palette.

### S4 — docs
Doc all new public items; keep the Smith + viewport panels untouched.

## Verification (run in worktree; all exit 0)
```
cargo fmt --check --all
cargo clippy -p yee-gui --all-targets -- -D warnings
cargo test -p yee-gui
cargo build -p yee-gui
git diff --stat -- crates/yee-plotters crates/yee-cli crates/yee-io crates/yee-mom '**/Cargo.toml'   # MUST be empty
```

## Escape-hatch
- If egui_plot's `Legend`/multi-`Line` needs an API not in the pinned
  egui_plot 0.35, or a new dependency, STOP + surface (do NOT bump egui or
  add a dep — the toolchain is pinned per CLAUDE.md §3).
- Do NOT attempt a brittle headless egui-render test; the tested surface
  is the series-building helper. Blocked >15 min → surface + stop.
- Run synchronously; no Monitor/ScheduleWakeup; no sub-agents.

## Out-of-scope (findings, not fixes)
* Smith-chart multi-trace / constant-R/X arc family (a separate track).
* Per-trace colour/style customization UI beyond a basic palette.
