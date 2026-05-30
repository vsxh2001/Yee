# Phase 2.fdtd.6.5 — reactive lumped-port V+I de-embedding bench — Plan

**Spec:** `2026-05-30-fdtd-6-5-reactive-port-deembed-bench-design.md` · **ADR:** ADR-0119

## Lane
`crates/yee-fdtd/**` ONLY (`src/lumped.rs` + `tests/reactive_deembed_001.rs`; may
update `ci.yml` for a release-gate job). May READ yee-voxel for context; edit
nothing else. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first). Branch `feature/fdtd-6-5-deembed-bench`.
First action: bring the canonical per-element updates from branch
`feature/fdtd-6-4-canonical-lc` (`021bed2`) onto this branch (cherry-pick or merge
that single commit) so `lumped.rs` has the per-edge-verified canonical port. Then
build the bench against it.

## Pattern files (READ FIRST)
- `docs/src/decisions/0117-fdtd-6-3-reactive-magnitude.md` + `0118-...` Outcomes —
  the contradiction (port-local proxy correct vs line-reflection wrong). This
  bench resolves it. Do NOT re-attempt coefficient fixes.
- `crates/yee-fdtd/tests/lumped_rlc_twoway_001.rs` — the parallel-plate TEM line
  harness (grid, full-width source/load sheet, two-run difference, DFT bins,
  calibration). REUSE its line + stepping; the new bench adds a **current**
  measurement (Ampère loop) alongside the existing voltage, to form Z = V/I.
- `crates/yee-fdtd/tests/cpml_reflection.rs` — the FDTD integration-test idiom.
- The canonical updates on `021bed2` (`crates/yee-fdtd/src/lumped.rs`).

## Steps
1. Cherry-pick `021bed2`'s `lumped.rs` change onto the branch (canonical port).
2. Add a current probe: `I(ω) = ∮ H·dl` (discrete Ampère loop) at the reference
   plane, single-bin-DFT'd like the existing voltage bin.
3. Measure `Z₀(ω) = V_inc/I_inc` from a matched/open run (no fitting).
4. For pure-L, pure-C, series-RLC and a known resistor: measure `Z_in = V/I`,
   de-embed `Z_L(ω)`, compare to `R + jωL + 1/(jωC)` at 4/6/9 GHz (+ one more).
5. ASSERT the resistor anchor `Z_L → R` (loose tol). Record the reactive `Z_L(ω)`
   table; assert the reactive arms IF they match, else assert the anchor + a
   `// VERDICT: port-{correct|wrong} — <numbers>` note.
6. Iterate IN THE CONTAINER (fast):
   `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
   --release --test reactive_deembed_001 -- --ignored --nocapture`
   (cargo direct or `bash -c '…'`, NEVER `bash -lc`).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd` + `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container, `bash -c`) → exit 0.
- No regression: `cargo test -p yee-fdtd --release --test lumped_lc_resonance
  --test lumped_resistor --test lumped_rlc_twoway_001 -- --include-ignored` GREEN.
- Gate: `reactive_deembed_001` GREEN (resistor anchor asserted).

## Escape hatch
Blocked > 60 min OR the current/Z₀ measurement itself won't validate (the resistor
anchor `Z_L → R` won't converge) → STOP and surface the bench you built, the
measured `Z₀(ω)`, and the resistor + reactive `Z_L(ω)` tables. Do NOT fake the
anchor, do NOT weaken it, do NOT touch yee-voxel/yee-filter. A working bench + an
honest verdict (even "port is wrong, here are the numbers") is the success
condition — this increment is a DECISION, not a port fix.

## Done when
DoD 1–4; resistor anchor asserted GREEN; the reactive `Z_L(ω)` verdict recorded;
diff = `crates/yee-fdtd/**` (+ optional `ci.yml`). The verdict drives increment 2.
