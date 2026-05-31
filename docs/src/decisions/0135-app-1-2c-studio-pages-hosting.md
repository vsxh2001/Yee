# ADR-0135: App.1.2c — deploy the Dioxus studio live via GitHub Pages

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0130 (the Dioxus studio merged + eframe retired — now needs to go
live), ADR-0107 (the studio web bundle), the desktop+web app goal (ADR-0089),
`docs.yml` (the working mdBook → Pages deploy), [[project-filter-design-final-goal]]

---

## Context

The Dioxus filter studio (`yee-studio-web`) is shipped on `main` (ADR-0130) and CI
builds its `dx build --platform web --release` bundle (ci.yml `wasm-build`), but it
is **not deployed** — the polished, pure-Rust web studio that is the broader app
goal's whole point is not yet accessible to anyone. `docs.yml` already deploys the
mdBook to GitHub Pages successfully (CLAUDE.md §8), so the Pages infra works; this
extends it to publish the studio.

## Decision

Publish the studio bundle to a **`/studio/` subpath** of the existing GitHub Pages
site (one Pages site per repo), alongside the mdBook docs at root:

- A unified Pages workflow (extend `docs.yml` or a new `pages.yml`) that:
  (1) builds the mdBook → `docs/book`; (2) builds the studio
  (`dx build --platform web --release`, with the **base-path** set for the project
  site so asset/wasm URLs resolve under the subpath) → into `docs/book/studio/`;
  (3) uploads the combined `docs/book` as the Pages artifact; (4) deploys via the
  standard `actions/deploy-pages@v4` (single `concurrency: pages` group).
- Trigger on `docs/**`, `crates/yee-studio-web/**`, and the workflow file.
- **Base-path:** a GitHub *project* Pages site serves at
  `https://<owner>.github.io/<repo>/` (here `/Yee/`), so the studio lives at
  `/Yee/studio/`. Set the dx web `base_path` (Dioxus.toml `[web.app] base_path`) so
  the built `index.html` references `/Yee/studio/assets/…` + the wasm loader resolves
  there (the verifiable crux — grep the built HTML).

## Consequences

**Ships:** the shipped Dioxus studio goes **live** on GitHub Pages at
`<owner>.github.io/Yee/studio/` (docs stay at root) — the broader app goal's
deployed, accessible web surface. Makes the whole filter-design journey usable in a
browser, no local build.

**Gate / verification (what's locally verifiable):** the studio bundle builds with
the project-site base-path; the **built `index.html` asset/wasm URLs resolve under
`/Yee/studio/`** (grep-checkable); the combined `docs/book` artifact contains the
docs at root + the studio at `/studio/`; the workflow YAML is valid + mirrors the
proven `docs.yml` deploy pattern. The **live deploy** runs on the next push to
`main` (a CI/CD step, like the existing docs deploy) — not locally reproducible.

**Requirement (documented):** GitHub repo Settings → Pages → Source = "GitHub
Actions" (already set for the docs deploy per CLAUDE.md §8). If the studio base-path
is wrong, the live page loads blank (404 on wasm/assets) — hence the built-HTML
URL check is the gate.

**Not in scope:** desktop/webview packaging; a custom domain; SSR/hydration.

---

## References
- `docs.yml` (the proven mdBook → Pages deploy); ci.yml `wasm-build` (the studio
  bundle build); ADR-0130 (the studio); ADR-0107 (the web bundle).
- `docs/superpowers/specs/2026-05-31-app-1-2c-studio-pages-hosting-design.md`;
  `docs/superpowers/plans/2026-05-31-app-1-2c-studio-pages-hosting.md`.
