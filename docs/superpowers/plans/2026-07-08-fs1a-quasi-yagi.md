# FS.1a implementation plan — quasi-Yagi

**Spec:** `docs/superpowers/specs/2026-07-08-fs1a-quasi-yagi-design.md`
**Lane:** `crates/yee-voxel/**`, `crates/yee-layout/**`,
`crates/yee-engine/tests/**` (+ docs lane).

1. **FS.1a.0**: add `ground_x_max_m: Option<f64>` to `VoxelOptions`
   (`..Default` construction sites unaffected; field defaults `None`).
   In `voxelize_microstrip`, restrict the k = 0 `Ex`/`Ey` PEC rows to
   `x ≤ ground_x_max_m` when set. Unit gate `voxel_002`:
   (a) `None` produces masks equal to a hand-built full-ground
   expectation on the reference layout (regression pin);
   (b) truncation cuts exactly at the requested column, traces unaffected.
2. **FS.1a.1**: `QuasiYagiDims` + `quasi_yagi(...)` in yee-layout
   (scaling-rule seeds, all `Polygon::rect`s); pattern file:
   `inset_fed_patch`. Gate `engine-antenna-005` mirroring
   `antenna_patch_inset.rs` (one solve, directional |S11|).
3. **FS.1a.2**: pattern gate via the A.2 NTFF machinery
   (`antenna_patch_pattern.rs` pattern file), F/B floor measured first,
   then pinned.
4. Verification at each step: fmt + clippy floor; `cargo test -p
   yee-voxel`; the release gates via the blanket engine CI step; ADR-0205
   at FS.1a.0+1, updated at FS.1a.2; FULL-SUITE-ROADMAP FS.1 row.
