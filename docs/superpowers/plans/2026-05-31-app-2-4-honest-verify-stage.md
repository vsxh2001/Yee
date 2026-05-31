# App.2.4 — Honest Verify stage — Plan

**Spec:** `2026-05-31-app-2-4-honest-verify-stage-design.md` · **ADR:** ADR-0141

## Lane
`crates/yee-studio-web/src/{engine.rs, stages.rs, main.rs}` ONLY. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-2-4-honest-verify`.

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/engine.rs` — the App.2.3 `topbar_view` pure helper + its
  test (MIRROR this pattern for `verify_view`); the flow structs: `Designed.report`
  (`MaskReport`: `pass`, `worst_passband_ripple_db`, `worst_return_loss_db`, `stopband:
  Vec<(f,ach,req,met)>`), `LumpedDesigned.verdict` (`MaskVerdict`: `pass`,
  `worst_passband_ripple_db`, `worst_return_loss_db`, `worst_stopband_rej_db`),
  `SteppedLowpassDesigned` (`pass`, `worst_passband_ripple_db`, `worst_return_loss_db`,
  `stopband: Vec<(f,ach,req,met)>`).
- `crates/yee-studio-web/src/stages.rs` — the current no-arg `verify_stage()` stub (the
  fake "—" cards to replace) + a real stage (e.g. `lumped_synthesis_stage`) for the
  card / chip / note RSX idiom.
- `crates/yee-studio-web/src/main.rs` — the `StageCanvas` `Stage::Verify =>
  stages::verify_stage()` arm (pass the signals) + how the other arms pass `designed`/
  `lumped`/`stepped`/`topology`.
- The spec §Method (the `VerifyView` shape + branches) + ADR-0141.

## Steps
1. `engine.rs`: `VerifyLevel` { RealizedLadder, SynthesizedIdeal } + `VerifyView` +
   `verify_view(topology, &Designed, Option<&LumpedDesigned>, &SteppedLowpassDesigned)`
   per the spec branches (stopband rej = min achieved over the `stopband` Vec for
   band-pass/stepped; `worst_stopband_rej_db` for lumped — `None` when absent/∞). Pure,
   documented.
2. `engine.rs`: the non-vacuous test (spec DoD §1) — per-flow metrics equal the source
   structs' fields, level differs (lumped RealizedLadder vs distributed SynthesizedIdeal),
   the lumped `None` case.
3. `stages.rs`: `verify_stage(topology, designed, lumped, stepped)` — render the three
   real metrics (ripple / return loss / stopband rejection; `"—"` only when genuinely
   `None`), the PASS / FAIL / "not realizable" chip, a `level` label ("Realized LC
   ladder" / "Synthesized ideal response"), and an honest EM note (full-wave EM of the
   board is a native step — the deferred ADR-0133 frontier — not run in the browser).
   Remove the fake "FDTD realized response" `"—"` cards + the "2.4 GHz" hard-code + the
   stale "App.D.5" claim.
4. `main.rs`: the `Stage::Verify` arm passes `topology()`, `designed`, `lumped`, `stepped`.

## Verify (run these; expected EXIT 0; quote output)
- `cargo test -p yee-studio-web` — the new non-vacuous `verify_view` test passes;
  existing tests unregressed. Quote "test result: ok".
- `cargo clippy -p yee-studio-web --all-targets -- -D warnings` ; `cargo fmt --check`.
- `cargo check --workspace`.
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0. Confirm the
  removed strings ("FDTD realized response" stat cards) are gone from `stages.rs`.

Commit on the branch: `yee-studio-web: honest topology-aware Verify stage (App.2.4,
ADR-0141)` + the Co-Authored-By trailer.

## Escape hatch
If a flow's metric field is not where expected, surface it — do NOT fabricate a metric
or leave a "—" where a real value exists. Keep `verify_view` pure (the component reads
signals, passes refs — bind guards to named locals as in App.2.3's TopBar). Do NOT add
any EM run, do NOT reopen ADR-0133. NEVER edit yee-filter/physics. Blocked > 30 min →
surface.

## Done when
`verify_view` + its non-vacuous test are green; the Verify stage shows the active flow's
real metrics + honest EM framing (no fake "—"/"2.4 GHz"/"App.D.5"); dx build EXIT 0;
existing flows unregressed; clippy/fmt/check clean; diff = `crates/yee-studio-web/src/
{engine.rs, stages.rs, main.rs}`. Then I verify + review + merge.
