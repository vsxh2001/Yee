# App.2.3 — TopBar shows the active flow's verdict — Plan

**Spec:** `2026-05-31-app-2-3-topbar-active-flow-verdict-design.md` · **ADR:** ADR-0140

## Lane
`crates/yee-studio-web/src/{engine.rs, main.rs}` ONLY. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-3-topbar-verdict`.

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/main.rs` — the current `TopBar(designed)` (the summary +
  PASS/FAIL chip to make topology-aware) and the `App` signal decls / `TopBar { designed }`
  call site (`topology`, `designed`, `lumped`, `stepped` are all already in scope).
- `crates/yee-studio-web/src/engine.rs` — `Designed` (`.spec`, `.order()`, `.report.pass`),
  `LumpedDesigned` (`.order()`, `.verdict.pass`, shares the band-pass spec),
  `SteppedLowpassDesigned` (`.spec`, `.order`, `.pass`, `.cutoff_hz()`); the existing
  `#[cfg(test)]` tests (add the new one there).
- The spec §Method (the `topbar_view` signature + the three branches) + ADR-0140.

## Steps
1. `engine.rs`: add `pub fn topbar_view(topology, designed: &Designed, lumped:
   Option<&LumpedDesigned>, stepped: &SteppedLowpassDesigned) -> (String, Option<bool>)`
   per the spec branches. Factor the approximation label (currently inlined in TopBar)
   as needed. Document it.
2. `engine.rs`: the non-vacuous test (spec DoD §1) — `SteppedImpedance` summary contains
   "cutoff" and not "%"; `EdgeCoupled` summary contains "%"; verdicts come from the
   right flow.
3. `main.rs`: `TopBar(topology, designed, lumped, stepped)` calls `topbar_view`, renders
   the summary chip + PASS / FAIL chip (verdict `None` → a muted "geometry not
   realizable" chip). Update the `App` call site to pass the four.

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — the new non-vacuous `topbar_view` test passes;
  existing tests unregressed. Quote the "test result: ok" line.
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0 (dx 0.6.3 +
  wasm32 already installed). (Studio is light; no Docker box, no FDTD.)

Commit on the branch: `yee-studio-web: TopBar shows the active flow's verdict (App.2.3,
ADR-0140)` + the Co-Authored-By trailer.

## Escape hatch
If threading the three signals into `TopBar` causes a Dioxus borrow/lifetime issue you
cannot resolve in 30 min, surface the specific blocker. Do NOT fake the verdict, do NOT
leave the band-pass-only TopBar with a misleading note. Keep `topbar_view` pure (no
signal reads inside it — the component reads the signals and passes references). NEVER
edit yee-filter / engine physics.

## Done when
`topbar_view` + its non-vacuous test are green; `TopBar` shows the active flow's summary
+ verdict (lumped ladder verdict; stepped low-pass cutoff + verdict); dx build EXIT 0;
existing flows unregressed; clippy/fmt/check clean; diff =
`crates/yee-studio-web/src/{engine.rs, main.rs}`. Then I verify + review + merge.
