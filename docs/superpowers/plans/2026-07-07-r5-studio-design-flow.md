# Plan — R.5 studio spec→design→export flow

**Spec:** `docs/superpowers/specs/2026-07-07-r5-studio-design-flow-design.md`

1. `studio/src-tauri/src/design.rs`: `FilterDesignRequest`/`Response` +
   `design_filter_impl` (synthesize → dims → layout → coupling-matrix response
   → `.s2p` + Gerber strings); `design_filter` Tauri command; path deps on
   yee-filter/layout/io/export.
2. `studio/src/views.tsx`: `SparamPlot` (SVG, |S21|/|S11| dB, −60 dB floor);
   `studio/src/App.tsx`: `FilterDesignPanel` (spec form, dims readout, plot,
   Blob-download export buttons).
3. Gates: `studio/src-tauri/tests/design_e2e.rs` (byte-checked `.s2p`/Gerber +
   band-pass response + unrealizable-spec error) and
   `studio/src/sparam.test.tsx` (DOM). CI: `cargo test --test design_e2e` step
   in the `studio-build` job.
4. ADR-0198, SUMMARY, RF-TOOL-ROADMAP R.5 row. Commit + push. R.5b (full-wave
   loop streaming in the studio) queued.
