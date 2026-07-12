# Plan: FS.3.2a — diagonal-edge geometry under full-wave test

**Spec:** `docs/superpowers/specs/2026-07-12-fs32a-diagonal-geometry-design.md`

1. `yee_layout::double_jog(w, jog_dy, run_x, MiterStyle)` + `MiterStyle`
   enum; serde-free pure generator returning `Layout` with two x-facing
   ports at equal y. Unit tests: port symmetry, mitered corner polygon
   has 5 verts with one 45° edge, square corners are rects.
2. `voxel-poly-001` in `crates/yee-voxel/tests/` (instant).
3. `engine-miter-001` in `crates/yee-engine/tests/engine_miter.rs`
   (`#[ignore]`, release): graded fixture, measured asserts per spec;
   run boxed, pin numbers from the first green run.
4. CI: add to the blanket yee-engine release gates step (it does not
   match the `--skip graded_`/antenna filters — verify) or a named step
   in the graded job. Verify: fmt, clippy, unit tests, gate.
5. ADR-0217, SUMMARY line, roadmap FS.3 row (+ FS.4.1 pointer note).
