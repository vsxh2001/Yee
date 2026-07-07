# Plan — S.11 / F1.2.1.0 EM-in-the-loop refinement

**Spec:** `docs/superpowers/specs/2026-07-06-s11-em-in-loop-refine-design.md`

1. Gate `crates/yee-filter/tests/engine_lpf_refine.rs` (`#[ignore]`, release):
   self-contained N=3 scenario (margins 30, ~8500 steps, CPML-xy, aperture ports);
   `verify_cutoff(f_c_synth) -> f_3db` helper runs the two-job measurement; the test
   does seed → correction → refined and asserts per the spec.
2. CI: one step in `compute-engine-gates` (≈ 4 release solves, ~10 min).
3. ADR-0188, ENGINE-STUDIO-ROADMAP S.11 row + footer, FILTER-DESIGN-ROADMAP F1.2.1.0
   note, SUMMARY.md.
4. fmt/clippy/fast tests; release gate run; record measured numbers; commit + push;
   continue to the next queued task.
