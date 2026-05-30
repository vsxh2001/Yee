# Phase 2.fdtd.6.9 — multi-cell aperture lumped port — Plan

**Spec:** `2026-05-30-fdtd-6-9-aperture-lumped-port-design.md` · **ADR:** ADR-0125

## Lane
`crates/yee-fdtd/**` ONLY (`src/lumped.rs` + `tests/aperture_port_001.rs`; `ci.yml`
for the gate job). Do NOT edit yee-voxel/yee-filter. Out of lane → finding.

## Base
New worktree off `main` (re-fetch first; main has the two-way port, per-axis CPML,
the de-embed bench). Branch `feature/fdtd-6-9-aperture-port`.

## Pattern files (READ FIRST)
- `docs/src/decisions/0125-fdtd-6-9-aperture-lumped-port.md` (the ADR — the
  mechanism + the 4-item formulation; esp. item 1 = aperture-area back-action is
  the ROOT fix for the O(dx²) inductor collapse) and ADR-0124 Outcome (the dx
  sweep + the O(dx²) diagnostic: inductor back-action `∝ dt²·dz/(ε₀·dA·L_cell)`
  with single-cell `dA=dx²`; capacitor `ε_eff = ε₀ + C_cell·dz/dA` frozen).
- `crates/yee-fdtd/src/lumped.rs` — the CURRENT two-way `correct_e` (the single-edge
  `V = E_z·dz`, the `(dt/(ε₀·dA))·I` single-cell back-action, `K`/`β`). You add an
  APERTURE coupling: modal `V = ∫E_z·dz` over the substrate height + back-action
  referenced to the aperture area `A = w·h`, distributed over `(y,z)`. Keep the
  single-edge path + resistor-exact intact.
- `crates/yee-fdtd/tests/reactive_deembed_001.rs` — the TRUSTWORTHY PEC-source
  de-embed harness (V=∫E·dz, modal I=∮H·dl, resistor anchor κ). Reuse it; the new
  gate runs it at TWO dx and adds the dx-stability assertion. (Do NOT use the
  matched-line bench — ADR-0123 showed it's corrupted.)
- `crates/yee-fdtd/tests/cpml_per_axis_001.rs` — `with_axes` usage if a matched
  guide helps.

## Steps
1. Add an aperture-port coupling to `LumpedRlcPort` (additive constructor/spec):
   the element spans `(y,z)` aperture cells; branch voltage = modal `∫E_z·dz` over
   the substrate height; the two-way update's back-action + `K`/`β` reference the
   aperture area `A = w·h` (NOT single-cell `dx²`), so the realized `Z_L` is
   dx-independent. Derive the per-cell distribution from the aperture normalization
   (not ad-hoc `C/N`). Keep `K+β>0` + the exact resistor limit.
2. `tests/aperture_port_001.rs`: de-embed pure-L / pure-C / series-RLC aperture
   ports at dx = 0.4 & 0.2 mm. Assert (a) resistor anchor `Z_L→R`, (b) reactive
   `Z_L` within a loose tol of `R+jωL+1/(jωC)`, (c) **dx-stability** — reactive
   `Z_L` agrees across the two dx (the O(dx²) collapse is GONE). This (c) is the
   decisive new check; prioritize killing the collapse.
3. `ci.yml`: a `fdtd-aperture-port-gate` release job (mirror `fdtd-per-axis-cpml-gate`).
4. Iterate IN THE CONTAINER:
   `YEE_BOX_DIR=<abs worktree path> scripts/yee-box.sh cargo test -p yee-fdtd
   --release --test aperture_port_001 -- --ignored --nocapture`
   (cargo direct or `bash -c`, NEVER `bash -lc`).

## Verify
- LOCAL light: `cargo fmt --check -p yee-fdtd` + `cargo clippy -p yee-fdtd
  --all-targets -- -D warnings` (container) → exit 0.
- No regression: `cargo test -p yee-fdtd --release --test reactive_deembed_001
  --test lumped_lc_resonance --test lumped_resistor --test lumped_rlc_twoway_001
  --test cpml_per_axis_001 --test cpml_reflection -- --include-ignored` GREEN.
- Gate: `aperture_port_001` GREEN (anchor + reactive tol + dx-stability).

## Escape hatch
This is a substantial core formulation. If the reactive tol isn't fully met but the
**O(dx²) collapse is killed** (reactive `Z_L` now dx-stable), that is real,
shippable progress — assert the dx-stability + record the residual; the remaining
accuracy is a follow-on. If blocked > 90 min OR the aperture back-action
destabilizes (NaN) / can't be made dx-stable → STOP and surface the update you
implemented, the per-dx `Z_L` tables, and where it diverges. Do NOT weaken the
resistor anchor, fake, or touch other crates. A precise partial is a good outcome.

## Done when
DoD 1–4: resistor anchor exact; reactive `Z_L` dx-stable (O(dx²) gone) + within
loose tol (or the dx-stability achieved + residual recorded); no regression; diff =
`crates/yee-fdtd/**` (+ ci.yml). Then F2.3 re-runs on the aperture port (next
sub-increment).
