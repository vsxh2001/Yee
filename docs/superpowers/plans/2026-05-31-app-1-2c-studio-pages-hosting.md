# App.1.2c — deploy the Dioxus studio live via GitHub Pages — Plan

**Spec:** `2026-05-31-app-1-2c-studio-pages-hosting-design.md` · **ADR:** ADR-0135

## Lane
`.github/workflows/**` (the Pages workflow) + `crates/yee-studio-web/Dioxus.toml`
(the base-path). Do NOT change studio source / engine crates. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first). Branch `feature/app-1-2c-pages`.

## Pattern files (READ FIRST)
- `.github/workflows/docs.yml` — the PROVEN mdBook → Pages deploy (build job →
  `upload-pages-artifact@v3` path `docs/book`; deploy job `deploy-pages@v4`;
  `permissions: pages: write, id-token: write`; `concurrency: pages`; the
  "Settings → Pages → GitHub Actions" requirement comment). EXTEND this (or
  supersede it with a unified `pages.yml`) — keep ONE Pages deploy.
- `.github/workflows/ci.yml` `wasm-build` job — the studio bundle build
  (wasm32 target, `cargo install dioxus-cli --version 0.6.3 --locked`,
  `dx build --platform web --release` in `crates/yee-studio-web`, bundle at
  `target/dx/yee-studio-web/release/web/public`).
- dx 0.6.3 web base-path config (find the right Dioxus.toml key for a subpath).
- ADR-0135 (the decision + the base-path crux + the verifiable gate).

## Steps
1. Set the dx web base-path in `crates/yee-studio-web/Dioxus.toml` for `/Yee/studio/`
   (e.g. `[web.app] base_path = "Yee/studio"` — verify the exact dx 0.6.3 key).
2. Unified Pages workflow: build mdBook → `docs/book`; build the studio (wasm32 +
   dx, base-path set) → copy `…/web/public/*` to `docs/book/studio/`; upload
   `docs/book` as the Pages artifact; deploy (deploy-pages@v4). Trigger on
   `docs/**` + `crates/yee-studio-web/**` + the workflow. If you add `pages.yml`,
   remove docs.yml's deploy job (keep its build logic folded in) so only ONE
   workflow owns the `pages` concurrency group + the Pages environment.
3. Document the Settings→Pages requirement + the live-on-push note in the workflow
   header (mirror docs.yml's comment).

## Verify (locally verifiable parts — the live deploy runs in CI on push)
- `cargo build`-free: install dx 0.6.3 + wasm32, `dx build --platform web --release`
  in `crates/yee-studio-web` with the base-path set → succeeds; **grep the built
  `target/dx/yee-studio-web/release/web/public/index.html`** → its asset/JS/wasm
  URLs are under `/Yee/studio/` (the live-correctness gate). Quote the URLs.
- `mdbook build docs/` → succeeds (docs unaffected).
- `cargo check --workspace` → green (Dioxus.toml change doesn't break the build).
- Validate the workflow YAML (parse / a yaml linter); confirm it uses
  `upload-pages-artifact@v3` + `deploy-pages@v4` + the pages permissions + a single
  `concurrency: pages`, and the combine step puts the studio at `docs/book/studio/`.
Commit on the branch: `ci: deploy the Dioxus studio to GitHub Pages at /studio/
(App.1.2c, ADR-0135)` + the Co-Authored-By trailer. (CI lane — the live deploy is
verified by the CI run on merge to main, the docs.yml precedent.)

## Escape hatch
If dx 0.6.3 has no working base-path for a subpath (the built HTML keeps absolute
`/assets` URLs), surface it — options: serve the studio at the Pages ROOT (and docs
under `/docs/`), or a post-build URL-rewrite step. Do NOT ship a workflow whose
built HTML URLs are wrong (the live page would be blank). Do NOT change studio
source / engine crates. Blocked > 45 min → surface.

## Done when
The Pages workflow builds mdBook + the studio (correct `/Yee/studio/` base-path,
verified in the built HTML) into one Pages artifact + deploys via the standard
pattern; YAML valid; docs + workspace unaffected; diff = `.github/workflows/**` +
`Dioxus.toml`. Then I (dispatcher) review + merge; the live deploy runs on the push.
