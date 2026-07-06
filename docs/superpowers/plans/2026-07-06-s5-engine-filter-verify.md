# Plan — S.5 engine-powered filter verify (walking skeleton)

**Spec:** `docs/superpowers/specs/2026-07-06-s5-engine-filter-verify-design.md`

1. **`yee-engine` protocol** (`crates/yee-engine/src/lib.rs`):
   - Add `MaterialsSpec` (serde; mirrors `yee_compute::Materials` field-for-field).
   - `JobSpec` gains `#[serde(default)] materials: Option<MaterialsSpec>` and
     `#[serde(default)] dt_s: Option<f64>`.
   - `run_job`: validate lengths (`(nx+1)(ny+1)(nz+1)` maps; per-component staggered
     masks) and `dt_s > 0` **before** constructing the stepper; on failure emit
     `JobEvent::Error` and return. Apply `spec.dt = dt_s` when present; hand
     `Materials` to both the CPU and GPU constructors.
2. **Fast tests** (in-crate `mod tests`): serde round-trip incl. materials/dt;
   small heterogeneous-job bit-parity vs direct `CpuFdtd`; three error paths.
3. **Gate** `crates/yee-engine/tests/verify_line_eeff.rs` — `engine-verify-001`,
   `#[ignore]`'d, release-only; dev-deps `yee-voxel` (path) + `yee-layout`
   (workspace), same as `yee-compute`'s. Scenario copied from `compute-008`
   with the drive/probes/materials expressed as a `JobSpec`; ε_eff vs
   Hammerstad–Jensen ≤ 15 %.
4. **CI**: extend the `compute-engine-gates` job with
   `cargo test -p yee-engine --release -- --include-ignored --nocapture`
   (fast tests re-run in seconds; the gate is the payload).
5. **Verify**: `cargo fmt`, `cargo clippy --workspace --all-targets -- -D warnings`,
   `cargo test -p yee-engine`, then the release gate with `--ignored --nocapture`
   (expect ~90 s / ~2.9 M cells, matching compute-008).
6. **Ship**: ADR-0182, `ENGINE-STUDIO-ROADMAP.md` S.5 row + footer, SUMMARY.md,
   commit + push.
