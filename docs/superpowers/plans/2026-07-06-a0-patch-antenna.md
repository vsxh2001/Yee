# Plan — A.0 patch antenna walking skeleton

**Spec:** `docs/superpowers/specs/2026-07-06-a0-patch-antenna-design.md`

1. `yee-layout`: `PatchDims` + `patch_antenna_dims` + `edge_fed_patch`; unit test vs
   hand-computed Balanis values (2.45 GHz FR-4: W ≈ 37.2 mm, L ≈ 28.8 mm).
2. Gate `yee-engine/tests/antenna_patch_s11.rs` (`#[ignore]`, release): dip position
   ±10 % of f₀, dip depth ≥ 2 dB below band median; print the S11 table.
3. CI: add to the yee-engine `--include-ignored` step (automatic).
4. ENGINE-STUDIO-ROADMAP Part 3 (A.*) section; ADR-0190; SUMMARY.
5. fmt/clippy/tests; release gate; record numbers; commit + push.
