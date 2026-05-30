# Phase 2.fdtd.6.7 — per-axis CPML face selection — Plan

**Spec:** `2026-05-30-fdtd-6-7-per-axis-cpml-design.md` · **ADR:** ADR-0122

## Lane
`crates/yee-fdtd/**` ONLY (`src/cpml.rs` + `tests/cpml_per_axis_001.rs`; may touch
`ci.yml` for a release-gate job). Edit nothing else. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first). Branch `feature/fdtd-6-7-per-axis-cpml`.

## Pattern files (READ FIRST)
- `crates/yee-fdtd/src/cpml.rs` — the CURRENT symmetric all-faces CPML:
  `CpmlParams` (npml, sigma_max, …), `CpmlState::new`, `pml_depth(i, n) ->
  Option<(depth, side)>`, `update_e`, `update_h`. This is what you generalize. Note
  the module docs already flag "per-face slabs is a future optimization."
- `crates/yee-fdtd/tests/cpml_reflection.rs` — the ≥30 dB reflection-vs-PEC idiom
  (build a grid, drive a pulse, measure peak E attenuation PEC vs CPML). MIRROR it
  for the x-only case; this all-faces gate must stay green with the default mask.
- `crates/yee-fdtd/src/boundary.rs` — `apply_pec` (the transverse-wall clamp the
  disabled axes rely on).

## Steps
1. Add `axes: [bool; 3]` to `CpmlParams`, **default `[true; 3]`** (a `Default`
   impl and/or a `with_axes([bool;3])` builder; keep `for_grid` behaviour
   identical — it sets all-true). Thread the mask into `CpmlState`.
2. In `update_e`/`update_h`/`pml_depth`, skip any axis whose flag is `false`
   (no CPML stretch / `pml_depth` returns `None` there). Verify `axes=[true;3]` is
   numerically identical to today (the existing tests pin this).
3. `tests/cpml_per_axis_001.rs` (`#[ignore]`'d): a guide with **x-only** CPML
   (`axes=[true,false,false]`) + PEC on y/z (apply_pec); drive an x-travelling
   Gaussian; measure interior-probe E reflection vs an all-PEC control → assert
   **≥30 dB reduction**. ALSO assert the transverse PEC walls intact (tangential E
   on y/z faces ≈ 0 → the guide mode survives, isn't absorbed).
4. `ci.yml`: a `fdtd-per-axis-cpml-gate` release job (mirror `fdtd-lumped-rlc-gate`).
5. Iterate IN THE CONTAINER:
   `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
   --release --test cpml_per_axis_001 -- --ignored --nocapture`
   (cargo direct or `bash -c '…'`, NEVER `bash -lc`).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd` + `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container, `bash -c`) → exit 0.
- No regression: `cargo test -p yee-fdtd --release --test cpml_reflection
  -- --include-ignored` GREEN (the all-faces gate, default mask). Also the
  `cpml_per_axis_001` gate GREEN.
- The unit tests in `cpml.rs` (`pml_depth_lookup_*`) stay green.

## Escape hatch
Blocked > 60 min, OR the x-only CPML cannot reach ≥30 dB without disturbing the
transverse PEC walls → STOP and surface: the mask plumbing you implemented, the
measured x-only reduction (dB), and the transverse-wall field levels. Do NOT
weaken the ≥30 dB target or the all-faces gate; do NOT touch other crates. A
precise partial is acceptable.

## Done when
DoD 1–4: `cpml_per_axis_001` GREEN (≥30 dB x-only + PEC walls intact), `cpml_reflection`
+ FDTD line/coupling gates non-regressed, fmt/clippy clean; diff =
`crates/yee-fdtd/**` (+ optional `ci.yml`). Enables the matched-line bench
(next sub-increment of ADR-0121's increment 3).
