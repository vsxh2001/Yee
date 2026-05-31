# App.2.2 — Low-pass stepped-impedance flow in the studio — Plan

**Spec:** `2026-05-31-app-2-2-studio-lowpass-stepped-design.md` · **ADR:** ADR-0139

## Lane
`crates/yee-filter/src/lib.rs` (the `ideal_response_lowpass` API + its gate) and
`crates/yee-studio-web/**` (the low-pass stepped-Z flow). NO yee-layout / engine-physics
edits (`dimension_stepped_impedance` already exists). Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-2-lowpass-stepped`.

## Pattern files (READ FIRST)
- **The lumped parallel flow is the TEMPLATE — mirror it for stepped-Z low-pass:**
  `crates/yee-studio-web/src/engine.rs` — `LumpedDesigned`, `design_lumped` /
  `design_lumped_from`, and the in-crate tests. `crates/yee-studio-web/src/stages.rs` —
  `lumped_synthesis_stage`, `lumped_layout_stage`, `Stage::LUMPED` rail, `Topology`,
  `technique_status`, `topology_label`/`topology_name`/`length_label`, `spec_stage`.
  `crates/yee-studio-web/src/main.rs` — the `lumped` signal + its re-derivation effect,
  `Stage::rail`, `StageCanvas` `lumped_flow` branching.
- `crates/yee-filter/src/lib.rs` — `ideal_response` (the band-pass analogue to mirror),
  the private `lowpass_s21_squared` (reuse at `Ω = f/f_c`), `dimension_stepped_impedance`
  + `SteppedImpedanceDimensions` / `SteppedSection` (the F1.2.3 engine you surface).
- `crates/yee-studio-web/src/svg.rs` — `board_svg(&Layout)` (generic; reuse).
- The spec §Method (the 2 parts) + ADR-0139.

## Steps
1. **yee-filter:** add public `ideal_response_lowpass(approx, order, cutoff_hz, freqs_hz)`
   reusing `lowpass_s21_squared` at `Ω = f/f_c`. Document it. Gate test (spec DoD §1).
2. **studio engine.rs:** `SteppedLowpassDesigned` (mirror `LumpedDesigned`) +
   `design_stepped_from(spec)` / `design_stepped()`: synth g-values, `dimension_stepped_impedance`
   sections, low-pass `|S21|` sweep via `ideal_response_lowpass`, low-pass mask bands +
   PASS/FAIL, the board `Layout`, board size, `dim_error`. The non-vacuous host test (DoD §2).
3. **studio stages.rs:** `Topology::SteppedImpedance`; `Stage::rail(SteppedImpedance)` =
   the stepped rail; `stepped_synthesis_stage` + `stepped_layout_stage`; `technique_status`
   `=> Live`; card `selects: Some(SteppedImpedance)`; `topology_label`/`topology_name`/
   `length_label` arms; `spec_stage` low-pass awareness (Cutoff label, hide FBW, set
   `Response::Lowpass`) when the active topology is SteppedImpedance.
4. **studio main.rs:** a `stepped` signal recomputed on spec edit (mirror `lumped`);
   `StageCanvas` `stepped_flow` branch routes Synthesis/Layout to the stepped renderers.

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-filter` — the new `ideal_response_lowpass` gate passes (quote it).
- `cargo test -p yee-studio-web` — the new non-vacuous routing test passes; existing
  band-pass + lumped tests unregressed.
- `cargo clippy -p yee-filter -p yee-studio-web --all-targets -- -D warnings` ;
  `cargo fmt --check`.
- `cargo check --workspace`.
- `dx build --platform web --release` (in `crates/yee-studio-web`, dx 0.6.3 + wasm32)
  → EXIT 0; confirm the served bundle built.
(All light/pure or wasm — host is fine; NO Docker box, NO FDTD.)

Commits (two ok): `yee-filter: ideal_response_lowpass + gate (App.2.2, ADR-0139)` and
`yee-studio-web: low-pass stepped-impedance flow (App.2.2)` + the Co-Authored-By trailer.

## Escape hatch (IMPORTANT — guarantees a clean partial)
This is a larger increment. If the full studio flow blocks > 45 min (the Spec low-pass
mode, the `stepped_flow` routing, or a Dioxus signal issue you cannot resolve), **land
the MINIMUM SHIPPABLE clean subset and surface the rest**: the `yee-filter`
`ideal_response_lowpass` + its gate AND the `SteppedLowpassDesigned` engine +
`design_stepped_from` + the non-vacuous engine test (all gateable on host) — and stop,
surfacing the remaining stage-UI wiring as a follow-on. Do NOT half-wire the UI, do NOT
fake a low-pass response or stub the sections, do NOT weaken either gate. A lit
SteppedImpedance card MUST route to the real `dimension_stepped_impedance` +
`ideal_response_lowpass`. NEVER edit yee-layout / the physics engines to force it.

## Done when
`ideal_response_lowpass` + its strong gate are green; the studio low-pass stepped-Z flow
is live (SteppedImpedance card routes to the real low-pass engine, non-vacuous test
passes); dx build EXIT 0; band-pass + lumped flows unregressed; clippy/fmt/check clean;
diff = `crates/yee-filter/src/lib.rs` + `crates/yee-studio-web/**`. Then I verify + review.
