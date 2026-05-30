# ADR-0130: App.D.2 — merge the Dioxus studio + retire the eframe view

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0110 (Dioxus redesign, POC-first — the maintainer's merge gate),
ADR-0120 (the lumped-LC studio stages), ADR-0107/0096/0090 (the eframe yee-studio +
its wasm/CI), the lumped-LC → PCB goal ("polished UI" component),
[[project-lumped-lc-and-studio-redesign]]

---

## Context

The Dioxus studio (`crates/yee-studio-web`, branch `feature/app-d0-dioxus-poc`) has
the full lumped-LC flow rendered from real engine data (ladder, E24/E96 BOM,
Monte-Carlo yield, dimensioned board) plus the POC distributed Synthesis/Layout;
Spec/Technique/Export are styled stubs. The maintainer reviewed it (served at
:8088) and chose (AskUserQuestion, 2026-05-31) **"Approve — build out + merge":**
fan out the remaining stages, retire the eframe studio, merge to main. This ships
the goal's **polished-UI** component.

## Decision

Merge the Dioxus studio to `main` as the studio, and retire the eframe view:

1. **Bring the distributed-flow stages to a shippable state** — Spec (the spec
   input form driving the synthesis), Technique (the topology gallery; distributed
   entries selectable or honestly "Soon"), Export (the param-sheet + downloads) —
   from styled stubs to real (the lumped flow is already real). Match the design
   system.
2. **Wire CI** — a `yee-studio-web` wasm build + `fmt`/`clippy` job (or extend the
   existing `wasm-build` job from `yee-studio` to `yee-studio-web`); ensure the
   workspace + the docs build stay green.
3. **Retire the eframe view** — remove the eframe app (`crates/yee-studio`'s
   `app.rs` + the eframe `main`/`desktop` feature path), keeping `StudioState` + the
   engine crates (which `yee-studio-web` reuses). If a clean delete is risky in one
   pass, **deprecate** (feature-gate off / mark deprecated) and delete in a tight
   follow-up — but the Dioxus studio is the studio on `main`.
4. **Merge to main** (`--no-ff`) after a code-review.

## Consequences

**Ships:** the goal's **polished-UI** component — a pure-Rust, web-first Dioxus
studio on `main` with the real lumped-LC journey (Spec → Technique → Synthesis →
Components+BOM → Tolerance → Layout) + the distributed Synthesis/Layout, replacing
the chunky eframe tool. With F2.0/F2.1/F2.2/F2.4 (and EM-sim pending Track A), the
lumped-LC goal's UI + component-choosing + BOM + tolerance + PCB are all on `main`
behind a polished UI.

**Gate:** the `yee-studio-web` wasm build + fmt/clippy green in CI; the workspace
+ docs build green; the eframe-retirement leaves `StudioState`/engine intact (the
non-studio crates unaffected). Code-review before merge.

**Risk / care:** retiring eframe touches `crates/yee-studio` + CI (the `wasm-build`
job referenced `yee-studio`). Reviewer verifies no orphaned refs, the workspace
builds, and the lumped flow renders the real engine data post-merge.

**Not in scope:** the EM-Verify stage (rides on Track A's EM-sim); desktop/webview
packaging; mobile; a fine distributed-flow polish pass (incremental after merge).

---

## References
- ADR-0110 (the POC-first gate the maintainer just cleared); ADR-0120 (the lumped
  studio stages); the maintainer's AskUserQuestion choice (2026-05-31).
- `docs/superpowers/specs/2026-05-31-app-d2-merge-dioxus-retire-eframe-design.md`;
  `docs/superpowers/plans/2026-05-31-app-d2-merge-dioxus-retire-eframe.md`.
