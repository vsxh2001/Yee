# FS.2 implementation plan — far-field products

**Spec:** `docs/superpowers/specs/2026-07-08-fs2-farfield-products-design.md`
**Lane:** `crates/yee-engine/**`, `crates/yee-compute/**` (port records),
`crates/yee-plotters/**`/`crates/yee-io/**` (FS.2c export) + docs.

1. **FS.2a**: `AperturePort` gains an optional per-step `(v_t, i)` record
   (CPU backend; GPU rejects recording ports with `Unsupported` — R.3
   idiom). Engine: `AperturePortSpec::record` (serde default false),
   `JobResult.port_records`. Gate `engine-power-001` (2 release solves):
   lossless-line energy bookkeeping port-A accepted ≈ port-B delivered,
   tolerance measured then pinned.
2. **FS.2b**: audit the `NtffState` |E| normalization against its docs;
   `yee_engine::farfield::gain_dbi`; gate `engine-gain-001` — patch
   broadside 5–8 dBi textbook window + the 2×1-vs-single **array-gain
   differential** (~+2.5–3 dB, cancels modeling bias).
3. **FS.2c**: efficiency (lossless ≈ 1 pinned band; lossy drops —
   direction gate) + full-sphere CSV export, byte-checked.
4. Per-step verification: fmt + workspace clippy floor; unit tests; the
   release gates; ADR-0207; FULL-SUITE-ROADMAP FS.2 row.
