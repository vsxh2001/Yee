# FS.0b.1 implementation plan — graded mesh rules + graded voxelization

**Spec:** `docs/superpowers/specs/2026-07-11-fs0b1-graded-rules-design.md`.

## Steps

1. **`yee-voxel/src/lib.rs`** — `GradedVoxelGrid` (per-axis primal spacing
   vectors + `x0_m`/`y0_m` origin + `k_gnd`/`k_top`),
   `GradedMicrostripModel` (raw `Array3` eps/PEC masks + node-coordinate
   vectors + port cells — no `YeeGrid`, whose scalar dx/dt are meaningless
   graded), `voxelize_microstrip_graded` mirroring `voxelize_inner`'s loop
   structure and z-stack against true cell centres/nodes from cumulative
   sums; port lookup by `partition_point` over nodes. Uniform entry points
   untouched.
2. **Gate `voxel-graded-001`** —
   `yee-voxel/tests/voxel_graded_001_uniform_bitexact.rs` (fast): constant
   arrays equal to dx vs `voxelize_microstrip` — exact equality on
   eps/PEC arrays, dims, port cells. Plus graded-z substrate fill and
   nonuniform-x port-lookup sanity tests.
3. **`yee-engine/src/automesh.rs`** — extract `trace_boxes` (shared with
   `min_feature_m`); `GradedMeshOptions` (+ `for_board`), `AutoSpacings`
   (+ `to_spacings`), `auto_spacings` with the per-axis rules (coarse =
   `auto_dx`; fine = `min(min_feature/2, coarse/2)`; edge±guard + gap fine
   intervals; geometric ladder ≤ growth; z substrate `ceil(h/(coarse/2))`
   exact-fit cells then growth into air; absorbers uniform coarse).
   Private `mesh_axis` marcher. Unit tests per the spec's fast-gate list.
4. **Gate `engine-graded-001`** —
   `yee-engine/tests/engine_graded_notch.rs` (`#[ignore]`, release): stub
   fixture from `board_automesh.rs`; `auto_spacings` grid; DUT + reference
   voxelized on the SAME grid; JobSpec built directly (CPML-xy npml 10,
   `dx_m = coarse`, `spacings` attached, dt from the engine); probe triples
   on scanned uniform-coarse stretches (12 coarse cells); double-ratio
   |S21| over 3.5–6 GHz / 50 MHz; assert notch ∈ 4.850 GHz ± 2 %, depth
   ≤ −20 dB, and `cells_graded/cells_uniform` under a measured-then-pinned
   ceiling (uniform = the dx0/2 = 0.267 mm pass-2 grid from
   `two_port_board_job` with the loop's rescaled options, built not
   solved). Print cells, ratio, runtime. Iterate rules (guard, fine) if
   tolerance missed; every iteration's numbers go in ADR-0210.
5. **Docs** — ADR-0210 + `SUMMARY.md` line: rules as shipped, iteration 0
   (brief-literal rule ≡ uniform pass-0 mesh ⇒ 5.100 GHz, 5.2 % — fails by
   ADR-0204's own measurement, no new solve), every solved iteration's
   notch/depth/cells/runtime, honest negatives.

## Verification

```sh
cargo fmt --check --all \
  && cargo clippy --workspace --all-targets -- -D warnings \
  && cargo test -p yee-voxel \
  && cargo test -p yee-engine --lib \
  && cargo test -p yee-engine --release --test engine_graded_notch -- --ignored --nocapture
```
