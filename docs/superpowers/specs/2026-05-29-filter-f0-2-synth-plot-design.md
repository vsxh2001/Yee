# Filter Phase F0.2 — `yee filter synth --plot` — Design Spec

**Phase:** F0.2 · **ADR:** ADR-0088 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
Make `yee filter synth` emit a visual |S21| spec-compliance plot: the
synthesized closed-form response with the `SpecMask` forbidden regions shaded.
Connects F0 synthesis ↔ 1.plotting.4 overlay. `yee-cli` only; no EM; no new dep.

## Change (yee-cli)
- `main.rs`: add `#[arg(long)] plot: Option<PathBuf>` to `FilterCommand::Synth`;
  dispatch `filter::run_synth(&spec, output.as_deref(), plot.as_deref())`.
- `filter.rs`:
  - `run_synth(spec_path, output, plot)` — unchanged behaviour when `plot` is
    `None`. When `Some(path)`: build `s21_db[i] = 20·log10(max(|s21[i]|, 1e-12))`
    over the existing sweep; build regions via the adapter; call
    `yee_plotters::draw_sparam_with_mask(path, &freqs, &[("S21", &s21_db)],
    &regions, &cfg)` where `cfg` is a `PlotConfig` titled from the spec
    (`format` from the path extension is NOT used — pass `PlotFormat` chosen by
    extension: `.svg`→Svg else Png; reuse the crate's PlotConfig default size).
    Print `wrote plot: <path>`.
  - `fn spec_mask_regions(spec: &FilterSpec) -> Vec<MaskRegion>`:
    - passband: `f1 = f0·(1 − fbw/2)`, `f2 = f0·(1 + fbw/2)` →
      `MaskRegion { f_lo_hz: f1, f_hi_hz: f2, kind: Floor, limit_db: −mask.passband_ripple_db }`.
    - for each `(f_s, reject_db)` in `mask.stopband`:
      `MaskRegion { f_lo_hz: f_s·0.98, f_hi_hz: f_s·1.02, kind: Ceiling, limit_db: −reject_db }`.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-cli --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-cli` exit 0.
4. Adapter unit test `spec_mask_regions_*`: for a Chebyshev BPF spec with
   `passband_ripple_db=0.5`, `fbw=0.10`, `f0=2e9`, one stopband `(2.4e9, 40.0)`,
   assert the returned vec has a `Floor` region with `f_lo≈1.9e9`, `f_hi≈2.1e9`,
   `limit_db≈−0.5`, and a `Ceiling` region with `f_lo≈2.352e9`, `f_hi≈2.448e9`,
   `limit_db≈−40.0` (tolerances ≤1 Hz-relative / ≤1e-9 dB).
5. CLI test `yee_filter_synth_plot_writes_png`: `yee filter synth <fixture>
   --plot <CARGO_TARGET_TMPDIR/x.png>` exits 0 and the PNG is non-empty
   (> 1024 bytes). Reuse the existing `cheb_bpf.toml` fixture.
6. `cargo run -p yee-cli -- filter synth crates/yee-cli/tests/fixtures/cheb_bpf.toml
   --plot /tmp/cheb.png` writes a PNG and exits 0.

## Out of scope
S11/return-loss mask plotting; GUI; EM; F1.1 FDTD extraction.
