# FS.0b.2a — graded two-port board fixture (plan)

**Spec:** `docs/superpowers/specs/2026-07-12-fs0b2a-graded-fixture-design.md`

1. `board.rs`: `GradedBoardOptions`, `GradedTwoPortBoardJob`,
   `two_port_board_jobs_graded` (logic lifted from
   `tests/engine_graded_notch.rs`); structural unit test.
2. Refactor `engine_graded_notch.rs` onto the fixture; asserts and pinned
   numbers unchanged.
3. Verify: `cargo test -p yee-engine` (fast), clippy floor, fmt; the
   release gate re-runs in CI (and locally once, boxed, to confirm the
   numbers did not move).
4. ADR-0210 addendum note + roadmap FS.0 row note.
