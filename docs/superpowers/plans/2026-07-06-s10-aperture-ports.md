# Plan — S.10 substrate-column ports

**Spec:** `docs/superpowers/specs/2026-07-06-s10-aperture-ports-design.md`

1. `yee_engine::column_port_specs` + unit test.
2. `engine_lpf_verify.rs`: drive + load become column ports; release re-run; record.
3. Adopt/revert by measurement; ADR-0187; roadmap S.10 row + footer; SUMMARY.
4. fmt/clippy/tests; commit + push; continue to the next queued task without stopping.
