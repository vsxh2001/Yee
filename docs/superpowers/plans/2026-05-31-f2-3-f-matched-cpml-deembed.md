# Filter Phase F2.3-f — matched-CPML board de-embed — Plan

**Spec:** `2026-05-31-f2-3-f-matched-cpml-deembed-design.md` · **ADR:** ADR-0131

## Lane
`crates/yee-voxel/**` ONLY (`src/lumped_sim.rs` + module doc). Do NOT edit
yee-fdtd/yee-filter. Out of lane → finding. (Cargo.lock/ci.yml from the merge.)

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(tip `9c57e5f` — aperture port + CW drive + the F2.3-e finer-grid commit). Merge
current `main` first (Cargo.lock `--theirs`, keep all CI jobs). `cargo check
-p yee-voxel` green in container.

## Pattern files (READ FIRST)
- `crates/yee-voxel/src/lumped_sim.rs` — `run_board_solve` (the CW drive + the
  lumped-resistor load + the load-cell voltage you replace with a matched-CPML
  output + transmitted-wave reference-plane measurement). Keep the aperture-port
  placement + the CW per-frequency drive.
- `crates/yee-fdtd/tests/cpml_per_axis_001.rs` (READ-ONLY) — how to drive an
  x-only-CPML guide: `with_axes([true,false,false])` + the custom step
  (update_h_only → cpml.update_h → source → update_e_only → cpml.update_e →
  transverse-PEC clamp → advance_clock). Mirror it for the F2.3 board.
- `crates/yee-fdtd/tests/cap_cw_001.rs` (READ-ONLY) — the CW steady-state amplitude
  measurement idiom.
- ADR-0131 (the matched-CPML de-embed) + ADR-0129 Outcome (the over-unity/collapse
  symptom) + ADR-0123 (why CPML≠matched-at-DC — but CW sidesteps it).

## Steps
1. Merge `main` into the branch; `cargo check -p yee-voxel` green (container).
2. In `run_board_solve`: terminate the microstrip with x-only CPML at the output
   (and input, behind the source) + transverse-PEC; lengthen the board so the
   reference planes clear the discontinuities + the CPML has room. Replace the
   lumped-resistor load + load-cell voltage with a transmitted-wave amplitude at an
   output reference plane (steady-state CW). Keep the aperture ports + the CW drive.
3. DUT/thru: `S21(f) = |V_out,ss| / |V_thru,ss|`.
4. Verify the |S21| is PHYSICAL (≤~1, no over-unity) + dx-stable (re-run at dx =
   0.4 & 0.2 mm — should now converge). Then re-run `fdtd_lumped_001`.

## Verify (bounded container — heavy FDTD + CPML)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- `... scripts/yee-box.sh cargo test -p yee-voxel --release --test fdtd_lumped_001
  -- --ignored --nocapture` → REPORT the |S21| (physical? dx-stable? notch depth) +
  GREEN-or-how-close. (cargo direct or `bash -c`, NEVER `bash -lc`.)

## Escape hatch
Do NOT weaken `fdtd_lumped_001`. If the matched-CPML de-embed gives a physical,
dx-stable |S21| but the notch is still shallow (< 20 dB), that is a CLEAN, valuable
result — record the |S21| (how close) → the residual is now isolated to the
aperture-port accuracy → the sub-cell reactance correction is the next sub-increment
(separate ADR). If the CPML termination destabilizes the CW board or the over-unity
persists, surface the precise behavior (the de-embed still isn't converging — quote
the |S21| vs dx). Blocked > 90 min → surface. Do NOT touch yee-fdtd/yee-filter.

## Done when
Either: `fdtd_lumped_001` GREEN at the strict 20 dB (physical + dx-stable de-embed +
the proven port reach it) → branch ready for review + the F2.3 merge, EM-sim ships
6/6. OR a precise "physical + dx-stable now, notch caps at X dB → port-accuracy
residual, sub-cell correction next" finding. diff = `crates/yee-voxel/**` (+ merge
artifacts).
