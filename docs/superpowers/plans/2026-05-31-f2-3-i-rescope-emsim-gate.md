# Filter Phase F2.3-i — re-scope the lumped EM-sim gate — Plan

**Spec:** `2026-05-31-f2-3-i-rescope-emsim-gate-design.md` · **ADR:** ADR-0134

## Lane
`crates/yee-voxel/**` ONLY (`tests/fdtd_lumped_001.rs` + `src/lumped_sim.rs` if a
small metric helper is needed). Do NOT edit yee-fdtd/yee-filter; do NOT touch
aperture_port_001 / cap_cw_001 / ladder_s21 (the isolation + circuit validations
stand). Out of lane → finding.

## Base / worktree
Existing worktree `worktrees/lumped-fdtd`, branch `feature/filter-f2-3-lumped-fdtd`
(tip `23a52e0` — the F2.3-h clean-launch + 2-point de-embed). Merge current `main`
first (Cargo.lock `--theirs`, keep all CI jobs). `cargo check -p yee-voxel` green.

## Pattern files (READ FIRST)
- `docs/src/decisions/0134-...md` (the principled-re-scope mandate + the integrity
  guardrail: the re-scoped assertion MUST be non-vacuous — fail for inert/broken).
- `crates/yee-voxel/tests/fdtd_lumped_001.rs` — the current gate (the ≥20 dB
  assertions to replace; the EM-sim pipeline + the de-embed to keep).
- ADR-0133 Outcome (the cavity wall — the F2.3-h |S21| numbers: thru over-unity at
  a box mode, the loaded sweep) + ADR-0124 Outcome (the INERT single-cell flat≈1
  response — the negative reference the re-scoped gate must distinguish from).
- ADR-0111 (`ladder_s21` — the sharp cross-validation to delegate to);
  ADR-0125/0127 (the isolation port gates).

## Steps
1. Merge `main`; `cargo check -p yee-voxel` green (container).
2. Determine, from the actual F2.3 board sweep (run it in the container), what is
   RELIABLY true + meaningful: the elements load the line (the DUT/thru deviates
   meaningfully + is frequency-dependent) vs the inert flat≈1 (ADR-0124). Pick a
   metric + a threshold set ABOVE the inert-noise floor (so it's non-vacuous) and
   BELOW what the cavity-limited measurement reliably delivers (so it's achievable).
3. Replace the ≥20 dB assertions with: finite/non-trivial sweep + the
   elements-load-the-line metric. Keep `ladder_s21` PASS-vs-Pozar referenced in the
   docstring as the sharp cross-validation; document the cavity wall (ADR-0133).
4. Confirm NON-VACUITY: show the asserted metric clears for the loaded board but
   would fail for the inert response (a comment with the inert≈0 dB margin vs the
   threshold, or a quick inert-control check).
5. Re-run `fdtd_lumped_001` → GREEN at the achievable bar.

## Verify (bounded container)
- `YEE_BOX_DIR=/home/hadassi/Code/Yee/worktrees/lumped-fdtd ... scripts/yee-box.sh
  bash -c 'cargo fmt --check -p yee-voxel && cargo clippy -p yee-voxel --all-targets
  -- -D warnings'` → exit 0.
- `... scripts/yee-box.sh cargo test -p yee-voxel --release --test fdtd_lumped_001
  -- --ignored --nocapture` → GREEN; REPORT the sweep + the asserted metric value
  + the inert-control margin (proving non-vacuity). (cargo direct or `bash -c`,
  NEVER `bash -lc`; bounded — minutes.)
- No-regression: aperture_port_001, cap_cw_001 still green (--include-ignored).

## Escape hatch
The re-scope MUST be principled + non-vacuous (the maintainer authorized re-scoping
the BAR, NOT faking a pass). If you cannot find a meaningful achievable assertion
that the loaded board passes AND the inert response fails, do NOT ship a vacuous
gate — STOP and surface that (the EM-sim genuinely doesn't demonstrate a
distinguishable real effect). Do NOT weaken to always-pass; do NOT touch
yee-fdtd/yee-filter/the isolation gates. Blocked > 60 min → surface.

## Done when
fmt/clippy clean; `fdtd_lumped_001` GREEN at the re-scoped achievable bar with a
demonstrated non-vacuity margin (loaded passes, inert fails); docstring delegates
the sharp cross-validation + documents the wall; no regression. Then I (dispatcher)
run a code-reviewer to CONFIRM the re-scoped gate is honest/non-vacuous → merge F2.3
→ EM-sim ships 6/6. diff = `crates/yee-voxel/**` (+ merge artifacts).
