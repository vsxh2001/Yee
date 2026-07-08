# FS.1b implementation plan — 2×1 corporate-fed patch array

**Spec:** `docs/superpowers/specs/2026-07-08-fs1b-patch-array-design.md`
**Lane:** `crates/yee-layout/**`, `crates/yee-engine/tests/**` (+ docs).

1. `PatchArrayDims` + `patch_array_2x1(f0, substrate, z0)` in yee-layout:
   spine (λg/2, probe room) → junction at x_v → λg/4 70.7 Ω transformers
   ±y → 50 Ω verticals to y = ±d/2 (d = 0.5 λ₀) → horizontal 50 Ω feeds →
   the 4-rect inset-patch construction (pattern file:
   `inset_fed_patch_with_depth`) translated to ±d/2, inset 0.25·L.
   Unit tests: symmetry, connectivity, patch dims, transformer length.
2. Gate `engine-antenna-007` (`antenna_patch_array.rs`, `#[ignore]`,
   A.1 fixture idiom, auto_dx-seeded, classic stack + A.2 open-top
   boundary): directional |S11| dip within ±10 % of 2.45 GHz; depth
   pinned after first measurement. Add to the antenna CI job when green.
3. Gate `engine-antenna-008` (`antenna_patch_array_pattern.rs`): y-z-cut
   NTFF beam narrowing + x-z patch-like cut; asserts pinned from the
   first instrumented run.
4. Verification per step: fmt + workspace clippy floor,
   `cargo test -p yee-layout --lib`, release gate runs, ADR-0206,
   FULL-SUITE-ROADMAP FS.1 row update.
