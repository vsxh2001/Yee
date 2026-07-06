# E.1 — CPML + per-cell materials on `yee-compute` (implementation plan)

**Spec:** `docs/superpowers/specs/2026-07-06-e1-cpml-materials-design.md`

1. `src/materials.rs` — `Materials` (+ shape validation against `FdtdSpec`), `Boundary`,
   `CpmlConfig` (with `for_spec`, mirroring `CpmlParams::for_grid` including
   `sigma_max_optimal`).
2. `src/cpml.rs` — flat-buffer `CpuCpmlState`: profile construction (`make_profiles` port,
   shared with the GPU upload path), `pml_depth`, rayon-slab `update_e` / `update_h` passes
   with per-cell arithmetic identical to `yee_fdtd::cpml::CpmlState`.
3. `src/cpu.rs` — per-cell arms in `update_h`/`update_e` (same match structure as the
   reference), `apply_pec_box`, `apply_pec_mask`, step orchestration (boundary phases in
   reference order), step counter, `with_config`, `step_with_gaussian_ez`.
4. Gate `tests/cpu_e1_reference_parity.rs` (compute-003) — bit-exact vs
   `WalkingSkeletonSolver` for both Cpml and PecBox scenarios. Run it before touching the GPU.
5. Gate `tests/cpml_reflection.rs` (compute-004) — ≥ 30 dB via `CpuFdtd`.
6. `src/gpu.rs` + `src/shaders/fdtd.wgsl` — arena-buffer refactor per spec §3.3; fused
   bulk+CPML kernels; mask clamp kernels; `with_config` (PecBox = host-side face zeroing).
   compute-002 must still pass unchanged.
7. Gate `tests/gpu_e1_parity.rs` (compute-005) — run for real on llvmpipe locally.
8. Docs: ADR-0176 with measured outcome; SUMMARY.md entry; ENGINE-STUDIO-ROADMAP E.1 →
   SHIPPED.
9. Verify: `cargo fmt --check --all`; `cargo clippy -p yee-compute --all-targets -- -D
   warnings`; `cargo test -p yee-compute` (all five compute gates green-or-skipped);
   `cargo check -p yee-compute --no-default-features`; `cargo check --workspace`.

**DoD:** step 9 green; compute-001/002 unchanged and green; compute-003 bit-exact;
compute-004 ≥ 30 dB; compute-005 green on llvmpipe; pushed to the feature branch.

**Lane:** `crates/yee-compute/**` + `docs/**` + `ENGINE-STUDIO-ROADMAP.md` (yee-fdtd is
read-only reference). Escape hatch: blocked > 15 min → surface and stop.
