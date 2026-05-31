# Filter Phase F2.3-h — clean forward-wave launch — Plan

**Spec:** `2026-05-31-f2-3-h-clean-forward-launch-design.md` · **ADR:** ADR-0133

## Lane
`crates/yee-voxel/**` ONLY (`src/lumped_sim.rs` + module doc). Do NOT edit
yee-fdtd/yee-filter. Out of lane → finding. (Cargo.lock/ci.yml from the merge.)

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(tip `4bde9dd` — the F2.3-g PEC-box 2-point de-embed). Merge current `main` first
(Cargo.lock `--theirs`, keep all CI jobs). `cargo check -p yee-voxel` green.

## Pattern files (READ FIRST)
- `crates/yee-voxel/src/lib.rs` `run_line_eeff` (ADR-0108) — the TIME-GATED
  incident-wave launch on a PEC line: how it launches a pulse + time-gates the
  forward incident wave before reflections. THIS is the clean `a₁` reference to
  adopt.
- `crates/yee-voxel/src/lumped_sim.rs` `run_board_solve` — the F2.3-g PEC-box +
  2-point de-embed (KEEP it) + the soft CW source (the weak-forward-launch problem)
  + the output probe region (β_out=0 issue). Fix the launch + output-probe placement.
- ADR-0133 (the clean-launch decision) + ADR-0132 Outcome (the floor diagnosis:
  source reflects, β_out=0, |b₂| degenerate).
- `crates/yee-fdtd` TF/SF source (ADR-0014/0021/0026) IF a directional launch is
  pursued (optional; the time-gated incident reference is the simpler first try).

## Steps
1. Merge `main`; `cargo check -p yee-voxel` green (container).
2. Clean forward `a₁`: adopt `run_line_eeff`'s time-gated incident-wave reference
   (launch into a long-enough lead-in, time-gate the forward incident before the
   first reflection). And/or a directional source (PEC/absorbing backing) to reduce
   the near-pure-standing-wave at the input.
3. Lengthen the line + place the output reference region clear of the PEC end wall /
   evanescent zone so a propagating forward wave is clean there (β_out>0, `b₂` above
   the floor).
4. DUT response: keep the CW steady-state (tanks ring up), referenced to the
   trustworthy `a₁` (hybrid time-gated-`a₁` + CW-`b₂` is fine; document it).
5. Re-run `fdtd_lumped_001`; report whether `a₁`/`b₂` are now well-resolved + the
   |S21| sweep + the disambiguation verdict.

## Verify (bounded container — minutes, dx=0.4mm, few freqs; NO finer-dx/multi-hour)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- `... scripts/yee-box.sh cargo test -p yee-voxel --release --test fdtd_lumped_001
  -- --ignored --nocapture` → REPORT a₁/b₂ resolution + |S21| + verdict.
  (cargo direct or `bash -c`, NEVER `bash -lc`. A single run should be MINUTES.)

## Escape hatch
Do NOT weaken `fdtd_lumped_001`. The deliverable is a TRUSTWORTHY launch (well-
resolved a₁/b₂, β>0) + a definitive verdict: band-pass ≥20 dB (ship), clean
inverted/notch-at-f0 (→ real topology inversion, cheap fix next), or
still-degenerate (→ surface the measurement-research wall). If a single run exceeds
~45 min or the launch still can't resolve b₂ above the floor, STOP + surface the
precise diagnostics. Blocked > 90 min → surface. Do NOT touch yee-fdtd/yee-filter;
do NOT merge/push.

## Done when
fmt/clippy clean; `fdtd_lumped_001` re-run with a trustworthy launch + the |S21| +
verdict reported. Either GREEN@20dB (→ review + F2.3 merge, EM-sim ships 6/6) OR a
definitive classification (topology inversion / measurement wall). diff =
`crates/yee-voxel/**` (+ merge artifacts).
