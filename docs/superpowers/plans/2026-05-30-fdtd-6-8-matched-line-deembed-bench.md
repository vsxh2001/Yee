# Phase 2.fdtd.6.8 — matched-line reactive de-embed bench — Plan

**Spec:** `2026-05-30-fdtd-6-8-matched-line-deembed-bench-design.md` · **ADR:** ADR-0123

## Lane
`crates/yee-fdtd/**` ONLY (`tests/reactive_deembed_matched_001.rs` NEW; `ci.yml`
for the gate job; `src/cpml.rs` only if a tiny accessor is unavoidable — prefer
not). May READ anything. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first; main has per-axis CPML `da66616`+).
Branch `feature/fdtd-6-8-matched-bench`.

## Pattern files (READ FIRST)
- `crates/yee-fdtd/tests/reactive_deembed_001.rs` — the EXISTING PEC-source bench:
  the parallel-plate guide, full-width source/load sheets, V=∫E·dz, modal I=∮H·dl,
  single-bin DFT, the resistor anchor, the per-load Z_L tables. REUSE the geometry
  + V/I measurement; CHANGE the boundaries to matched (x-only CPML both ends) and
  the de-embed to clean incident/reflected (no κ/A hack).
- `crates/yee-fdtd/tests/cpml_per_axis_001.rs` — how to drive an x-only-CPML guide:
  `CpmlParams::for_grid(&grid, NPML).with_axes([true,false,false])`, the custom
  step (update_h_only → cpml.update_h → source → update_e_only → cpml.update_e →
  transverse PEC clamp → advance_clock), and the transverse-PEC clamp helper.
- `crates/yee-fdtd/tests/cpml_reflection.rs` — the CPML reflection idiom.
- `docs/src/decisions/0123-...md` (the ADR — the de-embed math + the outcome gate)
  and the 0121 Outcome (why the PEC-source bench couldn't do this).

## Steps
1. Build the matched guide: parallel-plate, x-only CPML both ends + PEC transverse
   walls, soft `E_z` source near low-x (inside the line, past the source-end CPML),
   full-width lumped load at a load plane, a reference plane between them. Run long
   enough for the full reactive tail (no echo to truncate against — set the window
   generously; the source-end CPML absorbs the reflection after one pass).
2. Incident run (no/matched load): DFT `V_inc(ω)`, `I_inc(ω)` at the reference
   plane; `Z₀(ω)=V_inc/I_inc`.
3. Load runs (resistor anchor, pure-L, pure-C, series-RLC): total `V`,`I`;
   reflected = total − incident; `Z_in=V/I`; `Γ`; de-embed `Z_L(ω)`.
4. ASSERT the resistor anchor (`Z_L→R`, loose tol) — if it won't converge, the
   bench is wrong; fix it (the incident/reflected separation, the reference plane,
   the modal I loop), do NOT proceed to a reactive verdict.
5. Measure the reactive `Z_L(ω)`. If the well-conditioned capacitor is within
   `react_tol` of `1/(jωC)` → assert PORT-CORRECT (flip the verdict honestly,
   tell the dispatcher to re-run F2.3 + update ADR-0121/0123). Else pin the
   confirmed residual with a `// VERDICT: single-cell limit, brick 3 needed` note.
6. `ci.yml`: a `fdtd-matched-deembed-gate` release job (mirror
   `fdtd-per-axis-cpml-gate`).
- Container loop:
  `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
  --release --test reactive_deembed_matched_001 -- --ignored --nocapture`
  (cargo direct or `bash -c '…'`, NEVER `bash -lc`).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd` + `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container, `bash -c`) → exit 0.
- No regression: `cargo test -p yee-fdtd --release --test cpml_per_axis_001
  --test cpml_reflection --test reactive_deembed_001 --test lumped_lc_resonance
  --test lumped_resistor --test lumped_rlc_twoway_001 -- --include-ignored` GREEN.
- Gate: `reactive_deembed_matched_001` GREEN (resistor anchor + the reactive
  verdict).

## Escape hatch
Blocked > 60 min, OR the resistor anchor won't converge on the matched line (the
incident/reflected separation or the modal-I loop is off) → STOP and surface: the
bench you built, the measured `Z₀(ω)`, the resistor `Z_L(ω)`, and the reactive
`Z_L(ω)` tables. Do NOT fake/weaken the anchor; do NOT touch yee-voxel/yee-filter.
A working matched bench + an honest pinned verdict (PORT-CORRECT or single-cell
limit confirmed) IS the success condition — this brick is a DECISION.

## Done when
DoD 1–4: resistor anchor asserted GREEN on the matched line; the reactive `Z_L(ω)`
pinned with an explicit verdict; no regression; diff = `crates/yee-fdtd/**`
(+ `ci.yml`). The verdict drives either the F2.3 re-run or brick 3.
