# Plan: FS.5b.1 — ASM with the engine as fine model

**Spec:** `docs/superpowers/specs/2026-07-12-fs5b1-asm-em-fine-design.md`

1. Gate `crates/yee-filter/tests/sm_em_001.rs` (`#[ignore]`, release):
   stub layout parameterized by stub length; fine closure = graded
   two-port measure + parabolic notch refine, logging every eval;
   coarse = TL formula; `space_map` per spec; asserts per spec.
2. Run boxed/background, pin measured numbers.
3. CI: step in the yee-filter release job. ADR-0218, SUMMARY, roadmap
   FS.5 row.
