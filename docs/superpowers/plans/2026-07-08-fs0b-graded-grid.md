# FS.0b.0 implementation plan — graded grid in the yee-compute CPU kernel

**Spec:** `docs/superpowers/specs/2026-07-08-fs0b-graded-grid-design.md`.

## Steps

1. **`yee-compute/src/spec.rs`** — public `GradedSpacings` (`dx/dy/dz:
   Vec<f64>` primal widths) with `validate(&FdtdSpec)`,
   `validate_cpml_layers(npml, faces)`, `courant_limit()` (identical
   expression shape to `FdtdSpec::courant_limit`, over per-axis minima);
   crate-private `SpacingArrays` (`primal` len n, `dual` len n+1 per axis)
   with `uniform(&FdtdSpec)` and `graded(&FdtdSpec, &GradedSpacings)`
   constructors. Unit tests for dual arithmetic and validation errors.
2. **`yee-compute/src/cpu.rs`** — `CpuFdtd` carries `SpacingArrays`
   (uniform-filled by every existing constructor); new
   `set_spacings(&GradedSpacings)` (asserts validity, CPML-layer uniformity,
   and no dispersive map; `set_dispersive` asserts not graded).
   `update_h`/`update_e` divide by `primal`/`dual` array entries instead of
   the scalar `s.dx/dy/dz` — literal division kept so the uniform fill is
   bit-exact by construction. Resistive/aperture port updates use local
   spacings.
3. **`yee-compute/src/cpml.rs`** — `update_e`/`update_h` take
   `&SpacingArrays`; same primal/dual mapping as the bulk kernel.
4. **`yee-compute/src/lib.rs`** — export `GradedSpacings`; module docs.
5. **Gate compute-018** — `tests/graded_uniform_bitexact.rs` (fast): CPML +
   soft source + resistive port + aperture port + probes, and a PEC variant;
   scalar vs constant-array runs bit-identical (probes + all six field
   components, exact equality).
6. **`yee-engine/src/lib.rs`** — serde `GradedSpacings` + `JobSpec.spacings`
   (`serde(default)`); `run_job` validates (lengths/positivity, CPML layers,
   graded Courant dt incl. `dt_s` check), rejects `gpu`+spacings with the
   `ComputeError::Unsupported` message (auto → CPU) and `ntff`+spacings;
   `build_drive` computes aperture height/area as spacing sums when graded;
   CPU path attaches spacings. Tests: JSON round-trip + legacy default,
   protocol bit-exactness of constant arrays vs scalar, error cases.
   Update `JobSpec` struct literals across existing tests
   (`spacings: None`) — the permitted obvious-default exception
   (yee-engine tests, board.rs, yee-server + yee-filter test fixtures).
7. **Gate compute-019** — `tests/graded_interface_reflection.rs`
   (`#[ignore]`, release): measure the grading reflection floor via the
   upstream-difference method (spec §gates), print the number, then pin the
   assert just above the measured value.
8. **Docs** — ADR-0208 (+ `SUMMARY.md` line) recording the implementation
   choice (one kernel, literal division by spacing arrays) and the measured
   compute-019 number; honest negatives if any.

## Verification

```sh
cargo fmt --check --all \
  && cargo clippy -p yee-compute -p yee-engine --all-targets -- -D warnings \
  && cargo test -p yee-compute \
  && cargo test -p yee-engine --lib \
  && cargo test -p yee-compute --release --test graded_interface_reflection -- --ignored --nocapture
```
