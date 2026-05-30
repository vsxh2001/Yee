# Filter Phase F2.3-c — wire F2.3 onto the aperture lumped port — Plan

**Spec:** `2026-05-30-f2-3-c-aperture-port-driver-design.md` · **ADR:** ADR-0126

## Lane
`crates/yee-voxel/**` ONLY (`src/lumped_sim.rs` + module doc). Do NOT edit
yee-fdtd/yee-filter. Out of lane → finding. (`Cargo.lock`/`ci.yml` from the
main-merge are expected.)

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(tip `bbc7e26` — has the air-gap fix + the ADR-0124 single-edge sheet). **FIRST:
merge current `main` into the branch** (it needs the aperture port from 6.9 +
all the 6.x work). `Cargo.lock` → `--theirs` + `cargo check` in container; `ci.yml`
→ keep ALL gate jobs (CLAUDE.md §5). Commit the merge. Then also `rm` the untracked
`crates/yee-voxel/tests/scratch_*.rs` left from the investigation (harmless, but
clean them).

## Pattern files (READ FIRST)
- `crates/yee-voxel/src/lumped_sim.rs` — the CURRENT driver (the ADR-0124 sheet
  placement + the air-gap-fixed line-band detection). You replace the per-branch
  placement with one aperture port per branch.
- `crates/yee-fdtd/src/lumped.rs` `LumpedRlcPort::aperture(...)` + `ApertureSpec`
  (READ-ONLY) — the constructor signature + what the aperture spec needs (the
  `(y,z)` cells, `w`, `h`, the aggregate R/L/C). Mirror how `aperture_port_001.rs`
  builds an `ApertureSpec`.
- `docs/src/decisions/0126-...md` (the decision + the capacitor/window caveat),
  ADR-0125 (the aperture port + the cap CW caveat), ADR-0124 (the air-gap fix).

## Steps
1. Merge `main` into the F2.3 branch (Cargo.lock --theirs, keep all ci.yml jobs);
   `cargo check -p yee-voxel` green in the container. Remove the scratch_*.rs.
2. Build the `(y,z)` aperture per element: line band `[j_lo,j_hi)` (trace width,
   already detected) × substrate height `k=0..n_sub` at the element's x-column.
   Construct `ApertureSpec` accordingly.
3. Replace the per-branch placement: series branch → one `aperture(spec, ESR, L,
   C, None)`; shunt branch → `aperture(spec, ESR, L, ∞, None)` ‖ `aperture(spec,
   ESR, 0, C, None)` at the same aperture. Aggregate R/L/C (NO C/N split — the
   aperture port handles it). Keep drive/load ports.
4. Raise `n_steps` (or `LumpedSimConfig` default) for the capacitor steady state
   (a generous record; the band-pass needs the slow tail).
5. Re-run `fdtd_lumped_001` in the container; capture the FULL |S21| sweep.
6. Update the `lumped_sim.rs` module doc (aperture placement).

## Verify (bounded container — heavy FDTD)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- `... scripts/yee-box.sh cargo test -p yee-voxel --release --test fdtd_lumped_001
  -- --ignored --nocapture` → REPORT the |S21| sweep + GREEN-or-how-close.
- (cargo direct or `bash -c`, NEVER `bash -lc`.)

## Escape hatch
Do NOT weaken `fdtd_lumped_001` to force GREEN. If, on the aperture port + a long
window, the band-pass still doesn't form (the capacitor-under-transient limit),
that is an HONEST result: record the achieved |S21| (how close — quote the peak +
stopband dB) → the CW single-frequency drive is the next increment. Blocked > 60
min on the merge/build → surface. Do NOT touch yee-fdtd/yee-filter.

## Done when
fmt/clippy clean; `fdtd_lumped_001` re-run on the aperture port + |S21| sweep
reported. Either GREEN (→ branch ready for review + the F2.3 merge, EM-sim ships)
OR a precise "how close + CW drive needed" finding. diff = `crates/yee-voxel/**`
(+ main-merge artifacts).
