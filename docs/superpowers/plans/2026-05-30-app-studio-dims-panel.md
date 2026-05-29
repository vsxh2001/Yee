# App — yee-studio physical-dimensions panel — Plan

**Spec:** `2026-05-30-app-studio-dims-panel-design.md` · **ADR:** ADR-0099

## Lane
`crates/yee-studio/**` ONLY (`Cargo.toml`, `src/lib.rs`, `src/app.rs`, tests). Do
NOT edit `yee-filter`/`yee-layout`/any other crate. Out of lane → finding. Keep
`StudioState` egui-free + WASM-safe (ADR-0089) — the new derived `dims` must be
`Result<_, String>`, no `egui`/`eframe`/native type in `StudioState`.

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/app-studio-dims-panel`.

## Pattern files
- `crates/yee-studio/src/lib.rs` — `StudioState`, `recompute`/`apply_derived`,
  the existing derived fields (sweep/s21/mask) — add `dims` the same way; mirror
  the doc style. Note the `#[cfg(any(feature="desktop", feature="web"))] pub mod app;`.
- `crates/yee-studio/src/app.rs` — the existing panels (spec editor / synthesis /
  plot) — add the "Physical dimensions" section in the same idiom.
- `crates/yee-filter/src/dimension.rs` — `dimension_edge_coupled` /
  `EdgeCoupledDimensions` / `DimError` (call + read fields).
- `crates/yee-layout/src/lib.rs` — `Substrate` field names/units.

## Steps
1. `Cargo.toml`: add `yee-layout = { workspace = true }`.
2. `lib.rs`: add `eps_r`/`h_m` fields (FR-4 defaults) + derived
   `dims: Result<EdgeCoupledDimensions, String>`; compute it in recompute.
3. `app.rs`: the "Physical dimensions" section (editable εr/h + readout / error).
4. tests: `studio_state_dims` per DoD 4.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-studio --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-studio --jobs 2
nice -n 19 cargo check -p yee-studio --no-default-features --target wasm32-unknown-unknown --jobs 2
nice -n 19 cargo tree -p yee-studio --no-default-features --target wasm32-unknown-unknown -i egui   # egui MUST be absent (errors)
```
Do NOT run `cargo test --workspace`, FDTD, or mom-001. The native eframe build is
cached from App.1.2a; the wasm check is a fast incremental.

## Escape hatch
Blocked > 15 min — `StudioState` cannot hold `dims` without pulling an egui/native
type (breaking WASM-safety), or the default fixture's couplings can't be realized
on FR-4 (so `dims` is always `Err`) → STOP and surface the blocker (+ the
target_k vs achievable range if a fixture issue). Do NOT weaken the WASM-safety
gate (DoD 5) or fabricate dims. Do NOT edit `yee-filter`/`yee-layout`.

## Done when
DoD 1–5 pass; `git diff --stat <base>..HEAD` shows only `crates/yee-studio/**`
(+ `Cargo.lock`) + the 3 committed docs; `StudioState` still compiles egui-free
for `--no-default-features` wasm32.
