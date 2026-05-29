# App.0 — `yee-studio` desktop skeleton — Implementation Plan

**Spec:** `2026-05-29-app-0-yee-studio-desktop-design.md` · **ADR:** ADR-0090

## Lane
`crates/yee-studio/**` (new), root `Cargo.toml` (add member). Out of lane
(yee-gui or any other crate, docs already committed) → finding, not fix.

## Base
Worktree `worktrees/app0`, branch `feature/app-0-yee-studio`, base `06a5ff6`.

## Pattern files
- `crates/yee-gui/` — the eframe app shape: `Cargo.toml` (eframe `["wgpu"]` +
  egui + egui_plot deps, `[[bin]]`), `src/main.rs` (`eframe::run_native`),
  `src/app.rs` (`impl eframe::App`), `src/lib.rs`, `src/plots.rs` (egui_plot use).
  Imitate its module layout + the `[lints.rust]` manifest form.
- `crates/yee-cli/src/filter.rs` — the `spec_mask_regions` mapping + the
  401-pt sweep + `s21_db` computation to mirror (passband Floor at −ripple,
  per-stopband Ceiling at −reject over ±2%).

## Steps
1. `crates/yee-studio/Cargo.toml` per spec (deps: yee-synth, yee-filter, egui,
   eframe `["wgpu"]`, egui_plot, num-complex). bin `yee-studio` → `src/main.rs`.
2. `src/lib.rs`: `StudioState`, `MaskRegionView`, `from_spec`, `recompute` per
   spec §lib. NO egui types here (keep WASM-safe + testable). Doc every public item.
3. `src/app.rs`: `StudioApp { state: StudioState }` `impl eframe::App` with the
   left spec-editor `SidePanel`, central synthesis panel, and the `egui_plot`
   |S21|-vs-mask view. On spec edits call `state.recompute()`.
4. `src/main.rs`: thin `eframe::run_native` with a default Chebyshev BPF spec.
5. Root `Cargo.toml`: add `"crates/yee-studio"` to `members`.
6. Tests: `tests/studio_state.rs` (or `#[cfg(test)]` in lib) — `studio_state_
   recompute_pass` + `_fail` per spec §DoD 4. Headless (no egui rendering).

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-studio --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo build -p yee-studio --jobs 2          # windowed bin builds; NOT run
nice -n 19 cargo test -p yee-studio --jobs 2
```
egui/eframe/wgpu are already compiled as `yee-gui` deps (shared workspace
`target/`), so this is mostly link + the new crate's own code. Do NOT run the
windowed app (no display in CI). Do NOT run `cargo test --workspace`.

## Escape hatch
Blocked >15 min — eframe 0.34 API mismatch (run_native signature, `App::update`),
or egui_plot box-shading API → STOP, commit what compiles, surface the exact
error. Do NOT add a new dependency or pull in yee-gui internals.

## Done when
DoD 1–5 pass; `git diff --stat 06a5ff6..HEAD` shows only `crates/yee-studio/**`
+ root `Cargo.toml`/`Cargo.lock` + the 3 committed docs.
