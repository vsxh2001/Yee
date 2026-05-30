# App.D.2 — merge the Dioxus studio + retire the eframe view — Plan

**Spec:** `2026-05-31-app-d2-merge-dioxus-retire-eframe-design.md` · **ADR:** ADR-0130

## Lane
`crates/yee-studio-web/**`, `crates/yee-studio/**` (the eframe retirement),
`.github/workflows/ci.yml` (the wasm job), and the workspace `Cargo.toml` if the
member set changes. Do NOT change engine/physics crates (`yee-filter`, `yee-synth`,
`yee-layout`, `yee-export`, `yee-fdtd`, `yee-voxel`) or `StudioState`'s public
behaviour — `yee-studio-web` calls them read-only. Out of lane → finding.

## Base / worktree
Existing worktree `worktrees/studio-dx`, branch `feature/app-d0-dioxus-poc` (has
the lumped stages, ADR-0120). Merge current `main` first (Cargo.lock `--theirs`,
keep all CI jobs). Build/serve must work (dx 0.6.3 + wasm32 installed).

## Pattern files (READ FIRST)
- `crates/yee-studio-web/src/{main,stages,engine,svg}.rs`, `assets/studio.css` —
  the shell + the real lumped stages + the styled distributed stubs. Build out
  Spec/Technique/Export to match.
- `crates/yee-studio/` — the eframe app (`app.rs`, the eframe `main`, the `desktop`
  feature) to retire; and `StudioState` / the egui-free core to KEEP.
- `.github/workflows/ci.yml` — the `wasm-build` job (targets `yee-studio` today);
  repoint/add `yee-studio-web`.
- ADR-0130 + ADR-0110 (the POC-first gate, now cleared).

## Steps
1. Merge `main` into the branch; `cargo check --workspace` green.
2. Build out Spec / Technique / Export stages to real-enough (Spec input form →
   `synthesize`; Technique gallery; Export param-sheet + downloads). Match the
   design system. Distributed *fine* polish deferred; honest "Soon" labels OK.
3. Wire CI: a `yee-studio-web` wasm32 build + fmt/clippy job (repoint/extend
   `wasm-build`). Keep workspace + docs jobs green.
4. Retire the eframe view: remove `crates/yee-studio`'s eframe `app.rs` + eframe
   `main`/`desktop` render path, KEEP `StudioState` + the egui-free core + engine
   reuse. `cargo check --workspace` → no orphaned refs. (If a clean delete is risky
   in one pass, feature-gate it off + note the follow-up delete.)
5. Serve locally + smoke-check the real lumped flow renders real engine data.

## Verify
- `cargo fmt --check --all` + `cargo clippy --workspace --all-targets
  --no-default-features -- -D warnings` → exit 0 (the CI variant).
- `cargo build -p yee-studio-web --target wasm32-unknown-unknown` (+ the dx/trunk
  bundle as the POC did) → succeeds.
- `cargo check --workspace` → green, no orphaned `yee-studio` refs post-retirement.
- `mdbook build docs/` (or the docs CI step) → green.
- Serve the bundle + confirm the lumped stages render real values (smoke).

## Escape hatch
If retiring eframe in one pass orphans refs or risks the workspace build, FEATURE-
GATE the eframe app off (keep it compiling) + note the delete as a tight follow-up
rather than a messy big-bang delete — surface that choice. Do NOT change engine
crates / `StudioState` behaviour. Blocked > 60 min → surface. This branch is
maintainer-APPROVED to merge (ADR-0110 gate cleared) — but a code-review runs first.

## Done when
DoD 1–5: workspace + wasm + docs green, eframe retired (or cleanly gated), the
Dioxus studio is the studio. Then I (dispatcher) run a code-reviewer → merge
`--no-ff` to main → the polished-UI component ships. Report: what built out, the
CI change, how eframe was retired, the smoke-check, and any out-of-lane findings.
