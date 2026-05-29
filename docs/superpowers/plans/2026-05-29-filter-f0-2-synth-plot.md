# Filter Phase F0.2 — `yee filter synth --plot` — Implementation Plan

**Spec:** `2026-05-29-filter-f0-2-synth-plot-design.md` · **ADR:** ADR-0088

## Lane
`crates/yee-cli/**` only. Base: worktree `worktrees/filter-f02`, branch
`feature/filter-f0-2-synth-plot`, base `5800cb3`.

## Steps
1. `main.rs`: add `plot: Option<PathBuf>` to `FilterCommand::Synth` (doc'd
   `#[arg(long)]`); dispatch `filter::run_synth(&spec, output.as_deref(), plot.as_deref())`.
2. `filter.rs`: extend `run_synth` with a `plot: Option<&Path>` param. After the
   existing sweep + Touchstone write, if `plot` is `Some(path)`: compute
   `s21_db` (20·log10, floored at 1e-12), `let regions = spec_mask_regions(&spec)`,
   build a `PlotConfig { width_px:800, height_px:600, title, format }` (format =
   Svg if path ext is `svg` else Png), call
   `yee_plotters::draw_sparam_with_mask(path, &freqs, &[("S21", &s21_db)], &regions, &cfg)`,
   print `wrote plot: {path}`.
3. Add `fn spec_mask_regions(spec: &FilterSpec) -> Vec<yee_plotters::MaskRegion>`
   per spec §Change (passband Floor + per-stopband Ceiling). Doc it.
4. Tests (in `filter.rs` `#[cfg(test)]` for the adapter; `tests/cli_filter.rs`
   for the CLI): `spec_mask_regions_*` (known spec → expected regions) +
   `yee_filter_synth_plot_writes_png` (`--plot` → exit 0 + PNG > 1024 bytes,
   using CARGO_TARGET_TMPDIR + the cheb_bpf.toml fixture).

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-cli --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-cli --jobs 2
nice -n 19 cargo run -p yee-cli --jobs 2 -- filter synth crates/yee-cli/tests/fixtures/cheb_bpf.toml --plot "$CARGO_TARGET_TMPDIR/cheb.png" || cargo run -p yee-cli -- filter synth crates/yee-cli/tests/fixtures/cheb_bpf.toml --plot /tmp/cheb.png
```
yee-cli pulls yee-plotters → needs libfontconfig1-dev (present). Do NOT run the
workspace test suite (mom-001/fem-eig-003).

## Done when
DoD 1–6 pass; `git diff --stat 5800cb3..HEAD` shows only `crates/yee-cli/**`
+ the 3 committed docs; `filter synth` without `--plot` is byte-unchanged.
