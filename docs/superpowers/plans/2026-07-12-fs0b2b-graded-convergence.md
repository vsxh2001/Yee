# Plan: FS.0b.2b — graded convergence loop

**Spec:** `docs/superpowers/specs/2026-07-12-fs0b2b-graded-convergence-design.md`

1. `GradedMeshOptions.scale` (default 1.0 in `for_board`; validated
   `(0, 1]` by `auto_spacings`), applied to coarse + fine before band
   construction. Fix up every in-repo literal constructor.
2. `ConvergencePass.cells`; the uniform loop fills it from the job's
   `nx·ny·nz`, the graded loop from `GradedTwoPortBoardJob.cells`.
3. `converge_two_port_graded` per the spec (npml/spacing_cells rescaled
   per pass in coarse cells; scale/√2 per pass; linear ΔS criterion).
4. Unit tests (fast): scale=0.5 halves coarse and fine exactly; invalid
   scale rejected; npml/spacing rescaling arithmetic.
5. Gate `crates/yee-engine/tests/automesh_graded.rs` (`#[ignore]`,
   release): stub board per the spec's asserts. Run boxed locally, pin
   measured numbers.
6. Verify: fmt, clippy -D warnings, `cargo test -p yee-engine`, gate in
   release. ADR-0216, SUMMARY line, roadmap FS.0 row → FS.0b COMPLETE.
