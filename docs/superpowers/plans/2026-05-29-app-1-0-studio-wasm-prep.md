# App.1.0 — `yee-studio` WASM-prep — Implementation Plan

**Spec:** `2026-05-29-app-1-0-studio-wasm-prep-design.md` · **ADR:** ADR-0092

## Lane
`crates/yee-studio/**` ONLY. Out of lane (any other crate, root Cargo.toml) →
finding, not fix.

## Base
Worktree `worktrees/wasmprep`, branch `feature/app-1-0-studio-wasm-prep`,
base `3b63e87`.

## Steps
1. `crates/yee-studio/Cargo.toml`: set `eframe`/`egui`/`egui_plot` `optional = true`;
   add `[features] default = ["desktop"]`, `desktop = ["dep:eframe","dep:egui","dep:egui_plot"]`.
2. `src/lib.rs`: `#[cfg(feature = "desktop")] pub mod app;`. Leave `StudioState`,
   `MaskRegionView`, `recompute`/`apply_derived`, sweep/mask helpers un-gated.
   Update the `TODO(App.1)` doc comment to reflect that the gate now exists and
   only the wasm32/trunk build remains.
3. `src/main.rs`: `#[cfg(feature = "desktop")]` on the eframe `main` (and its
   `use eframe`/`egui` imports); add a `#[cfg(not(feature = "desktop"))] fn
   main() { println!("yee-studio built without the `desktop` feature (no GUI); …"); }`
   so the bin target links under `--no-default-features`.
4. Confirm the tests (`studio_state_recompute_*`) are NOT under `mod app` and use
   no egui types, so they compile + pass under `--no-default-features`.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-studio --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo clippy -p yee-studio --no-default-features --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo build -p yee-studio --jobs 2
nice -n 19 cargo build -p yee-studio --no-default-features --jobs 2
nice -n 19 cargo tree -p yee-studio --no-default-features | grep -E "eframe|egui" && echo "FAIL: egui present" || echo "OK: egui absent"
nice -n 19 cargo test -p yee-studio --jobs 2
nice -n 19 cargo test -p yee-studio --no-default-features --jobs 2
```
(The `grep … && echo FAIL || echo OK` line must print `OK: egui absent`.)

## Escape hatch
Blocked >15 min — a non-egui item in lib.rs unexpectedly depends on egui (so
`--no-default-features` won't compile), or the bin can't link without the feature
→ STOP, surface the exact item. Do NOT move logic into `app` to dodge it.

## Done when
DoD 1–6 pass (esp. egui absent from the no-default-features dep tree);
`git diff --stat 3b63e87..HEAD` shows only `crates/yee-studio/**` + the 3 docs.
