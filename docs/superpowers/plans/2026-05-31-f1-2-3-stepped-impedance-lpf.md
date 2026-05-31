# Filter F1.2.3 — Stepped-Impedance LPF — Plan

**Spec:** `2026-05-31-f1-2-3-stepped-impedance-lpf-design.md` · **ADR:** ADR-0137

## Lane
`crates/yee-filter/**` (the synthesis + dimensions + gate). `crates/yee-layout/**`
ONLY if a missing helper is genuinely needed (prefer composing existing helpers).
Out of lane → finding, do NOT fix in place. NO studio changes this increment.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/f1-2-3-stepped-impedance`.

## Pattern files (READ FIRST)
- `crates/yee-filter/src/dimension.rs` — `dimension_edge_coupled` /
  `dimension_edge_coupled_layout` / `EdgeCoupledDimensions` and the `dimension_hairpin`
  pair. **Mirror these exactly** for style: the `yee_layout` imports (`microstrip_width`,
  `eps_eff`, `Layout`, `Substrate`), the `DimError` handling, the doc density,
  `c = 299_792_458`, the `_layout` companion. Add the stepped-Z functions alongside.
- `crates/yee-synth/src/lib.rs` — `Prototype { g: [g0..g_{N+1}] }`, `prototype(approx,
  order)`, `min_order`. g[1..=N] are the reactive elements; `proto.order()` = N.
- `crates/yee-filter/src/lib.rs` — the crate-root re-export block for the dimension
  types/functions (add the stepped-Z ones).
- The spec §Method (the exact βl formulas) and the ADR.

## Steps
1. `SteppedSection` + `SteppedImpedanceDimensions` structs (spec §Types). Document all.
2. `dimension_stepped_impedance(proto, f_c, z0, z_high, z_low, sub)`: for k=1..=N, map
   `g_k` to a section — **section 1 (k=1) is shunt-cap / low-Z** (`high_z=false`),
   alternating. `βl = g_k·z_low/z0` (low-Z) or `g_k·z0/z_high` (high-Z). Width via
   `microstrip_width(z, εr, h)`; `ε_eff` via `eps_eff(width, h, εr)`; `λg = c/(f_c·√ε_eff)`;
   `length = βl/(2π)·λg`. Return `DimError` on a non-physical input.
3. `dimension_stepped_impedance_layout(...)`: place the sections in-line (mirror
   `dimension_edge_coupled_layout`'s `Layout`/`Substrate` construction).
4. Re-export the new items from the crate root.
5. Gate `crates/yee-filter/tests/dim_stepped_001.rs` (spec DoD §1): Pozar Example 8.6
   (Butterworth N=6, f_c=2.5e9, z0=50, z_high=120, z_low=20) → assert the six electrical
   lengths in degrees within ±1.0° of `[11.85, 33.76, 44.28, 46.12, 32.41, 12.34]`,
   assert section 0 is low-Z, assert all physical lengths > 0 and finite. Derive βl from
   the function output (do NOT hardcode the expected as the computed).

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-filter --test dim_stepped_001` — quote the "test result: ok" line.
- `cargo test -p yee-filter` (full crate, no regressions).
- `cargo clippy -p yee-filter --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
(yee-filter is light/pure Rust — host is fine; NO Docker box, NO FDTD here.)

Commit on the branch: `yee-filter: stepped-impedance LPF synthesis + dimensions + Pozar
§8.6 gate (F1.2.3, ADR-0137)` + the Co-Authored-By trailer.

## Escape hatch
If the Pozar gate values don't reproduce within ±1°, **surface it** — recheck the βl
formula and the prototype-element/alternation convention (Pozar starts with a shunt
capacitor); do NOT widen the tolerance or hardcode the answer to force a pass. If a
`yee_layout` helper is missing, surface it rather than editing yee-layout speculatively.
Blocked > 30 min → surface and stop.

## Done when
`dimension_stepped_impedance` + the non-vacuous Pozar §8.6 gate are green; clippy/fmt/
check clean; diff = `crates/yee-filter/**` (+ `crates/yee-layout/**` only if justified).
Then I (dispatcher) verify + review + merge. Studio low-pass wiring is a follow-on.
