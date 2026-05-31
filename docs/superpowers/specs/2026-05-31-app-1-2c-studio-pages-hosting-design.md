# App.1.2c — deploy the Dioxus studio live via GitHub Pages — Design Spec

**ADR:** ADR-0135 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

The shipped `yee-studio-web` Dioxus studio (ADR-0130) is not deployed — the broader
app goal's web surface isn't accessible. `docs.yml` already deploys the mdBook to
Pages (works). Extend Pages to publish the studio at a `/studio/` subpath.

## Goal

The studio is live on GitHub Pages at `<owner>.github.io/Yee/studio/` (docs stay at
root), built + deployed in CI, with the base-path correct so assets/wasm resolve.

## Method

A unified Pages workflow (extend `docs.yml` → rename intent, or a new `pages.yml`
that supersedes the docs-only one — keep ONE Pages deploy, ONE `concurrency: pages`
group, to avoid two workflows racing the same Pages environment):

1. Build mdBook → `docs/book` (as today).
2. Build the studio: install wasm32 + dx 0.6.3, set the dx web `base_path` to the
   project-site subpath (`Yee/studio` — so URLs are `/Yee/studio/…`), `dx build
   --platform web --release` (working-dir `crates/yee-studio-web`), copy
   `target/dx/yee-studio-web/release/web/public/*` → `docs/book/studio/`.
3. Upload `docs/book` (now docs + `studio/`) as the Pages artifact; deploy via
   `actions/upload-pages-artifact@v3` + `actions/deploy-pages@v4` (the docs.yml
   pattern). `permissions: pages: write, id-token: write`.
4. Trigger: push to `main` on `docs/**`, `crates/yee-studio-web/**`, the workflow.

**Base-path is the crux:** set `[web.app] base_path = "Yee/studio"` in
`crates/yee-studio-web/Dioxus.toml` (or the dx-equivalent) so the built
`index.html` references `/Yee/studio/assets/…` and the wasm loader path is correct.
If dx 0.6.3's base_path key differs, find the right one (dx docs) — the test is the
built HTML's URLs.

## Changes

- `.github/workflows/` — the unified Pages workflow (build mdBook + studio,
  combine, deploy). If extending `docs.yml`, update its name/comment; if a new
  `pages.yml`, retire docs.yml's deploy to avoid two Pages deploys.
- `crates/yee-studio-web/Dioxus.toml` — the `base_path` for the `/Yee/studio/`
  subpath.
- Do NOT change studio source / engine crates.

## DoD (machine-checkable, locally verifiable parts)

1. `dx build --platform web --release` (with the base-path set) succeeds; the built
   `target/dx/yee-studio-web/release/web/public/index.html` references its
   assets/wasm under `/Yee/studio/` (grep the built HTML — the live-correctness gate).
2. The workflow YAML is valid (yamllint or a parse) + mirrors the proven docs.yml
   deploy (upload-pages-artifact@v3 + deploy-pages@v4, pages permissions, single
   `concurrency: pages`); the combine step places the studio at `docs/book/studio/`.
3. `mdbook build docs/` still works (docs unaffected); no studio-source/engine
   change; `cargo check --workspace` green (no Dioxus.toml break).
4. Document (in the workflow + ADR) the Pages-source requirement + that the live
   deploy runs on the next push to main.

## Out of scope

Desktop/webview; custom domain; SSR. The actual live deploy (runs in CI on push).

## Why

It makes the shipped, polished studio **accessible** — the broader app goal's
deployed web surface — extending the proven docs Pages deploy, with the one risk
(base-path) locally verified via the built HTML URLs.
