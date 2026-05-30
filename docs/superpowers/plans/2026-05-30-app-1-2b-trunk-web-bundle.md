# App.1.2b — `trunk` web bundle + CI build — Plan

**Spec:** `2026-05-30-app-1-2b-trunk-web-bundle-design.md` · **ADR:** ADR-0107

## Lane
`crates/yee-studio/**` (new `index.html`, `Trunk.toml` — do NOT touch `src/` or
`Cargo.toml`; the `web` feature + `start_web` entry already exist) AND
`.github/workflows/ci.yml` (extend the wasm job). NOTHING ELSE. Out of lane →
finding. Do NOT add deps, do NOT edit `src/main.rs`/`src/app.rs`/`src/lib.rs`.

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/app-1-2b-trunk-web`.

## Pattern / context files (READ first)
- `crates/yee-studio/src/main.rs` (~lines 75–100) — `start_web()` mounts on
  `web_sys` canvas `the_canvas_id`; confirm the exact id string. The `[[bin]]`
  name is `yee-studio` (Cargo.toml).
- `crates/yee-studio/Cargo.toml` — the `web` feature list (already complete);
  do NOT modify it.
- `.github/workflows/ci.yml` (~lines 75–96) — the existing `wasm-build` job
  (`runs-on: ubuntu-latest`, `dtolnay/rust-toolchain@stable` toolchain 1.92 +
  `targets: wasm32-unknown-unknown`, `Swatinem/rust-cache@v2`, then `cargo build
  -p yee-studio --no-default-features --target wasm32`). MIRROR its style.

## Steps
1. `crates/yee-studio/index.html`: `<!DOCTYPE html>` … `<canvas
   id="the_canvas_id"></canvas>` + `<link data-trunk rel="rust"
   data-bin="yee-studio" data-cargo-no-default-features data-cargo-features="web"
   data-wasm-opt="2" />` + minimal full-viewport CSS + `<title>`. Match the id
   string EXACTLY to what `start_web` queries.
2. `crates/yee-studio/Trunk.toml`: `[build]` with `target = "index.html"` and
   `dist = "dist"`.
3. `.github/workflows/ci.yml`: in the `wasm-build` job, AFTER the existing
   `--no-default-features` compile step, add steps: install trunk (prefer
   `jetli/trunk-action@v0.5.0`, or `cargo install --locked trunk` as fallback) +
   the matching `wasm-bindgen-cli`; `trunk build --release
   crates/yee-studio/index.html`; `actions/upload-artifact@v4` with `path:
   crates/yee-studio/dist`. Keep the existing step intact. Do NOT add a Pages
   deploy.

## Verify (light only — MEMORY-constrained box; do NOT run trunk / wasm build)
```
test -f crates/yee-studio/index.html && grep -q 'the_canvas_id' crates/yee-studio/index.html && grep -q 'data-cargo-features="web"' crates/yee-studio/index.html && echo HTML_OK
test -f crates/yee-studio/Trunk.toml && grep -q '\[build\]' crates/yee-studio/Trunk.toml && echo TOML_OK
python3 -c "import yaml,sys; yaml.safe_load(open('.github/workflows/ci.yml')); print('YAML_OK')"
```
If `actionlint` is installed, also run it on `ci.yml`. Do NOT run `trunk build`,
`cargo build --features web --target wasm32`, or any workspace build — those are
heavy and OOM-kill this box; CI proves the bundle builds. App.1.2a already
proved the `--features web` wasm compile is green.

## Escape hatch
Blocked > 15 min (uncertain trunk `data-trunk` directive syntax for selecting a
crate feature + bin; whether trunk needs an explicit `wasm-bindgen-cli` version
pin; YAML won't validate) → STOP and surface the exact uncertainty. Do NOT run
the heavy wasm/trunk build locally to "see if it works"; do NOT add a Pages
deploy; do NOT edit yee-studio `src/` or `Cargo.toml`.

## Done when
DoD 1–4 pass; `git diff --stat <base>..HEAD` = only `crates/yee-studio/{index.html,Trunk.toml}`
+ `.github/workflows/ci.yml` (+ the 3 committed docs); no `src/`/`Cargo.toml`
change.
