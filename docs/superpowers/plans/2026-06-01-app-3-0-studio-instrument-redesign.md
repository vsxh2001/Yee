# App.3.0 — Studio "Instrument" redesign — Plan

**Spec:** `2026-06-01-app-3-0-studio-instrument-redesign-design.md` · **ADR:** ADR-0152

## Lane
`crates/yee-studio-web/{assets/studio.css, assets/*.woff2 (new fonts if self-hosted),
index.html, src/stages.rs, src/main.rs, src/svg.rs}`. NO `yee-filter`/engine edits.
Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first; main 5caaa38). Branch `feature/app-3-0-instrument`.

## Pattern files (READ FIRST — edit ONLY in the worktree)
- `crates/yee-studio-web/assets/studio.css` — the FULL design system you're elevating: tokens
  (`:root` — copper `--copper #e6b24d` / substrate `--substrate #3f9e72` / teal `--accent
  #2dd4bf` / `--muted #828c98` AA), `.topbar`/`.brand`, `.rail`/`.item.on`, `.card`,
  `.plot`/`.board-frame`, `.tcard`/`.tcard.sel`, the `.rec-*` recommender, the `@media`
  blocks (App.2.9 — KEEP + extend). It currently has `0 @keyframes` — you ADD the motion.
- `crates/yee-studio-web/src/stages.rs` + `src/main.rs` — the rsx! shell (topbar/rail/canvas),
  the stage renderers, the board-frame + plot wrappers, the wordmark. Add entrance-stagger
  classes, the board-hero wrapper, draw-in SVG attrs, keep all a11y attrs (App.2.9).
- `crates/yee-studio-web/src/svg.rs` — the SVG plot/board emitters (for the trace/board
  draw-in: add `stroke-dasharray`/`stroke-dashoffset` + a CSS-animated class, or an opacity
  reveal — without breaking the existing geometry).
- `crates/yee-studio-web/index.html` — the dx template (has `lang="en"`); add the font
  `<link>`/preload here if using a CDN, or `@font-face` in studio.css if self-hosting.
- The frontend-design skill principles (already in the dispatcher's context — embedded in the
  brief): bold, intentional, distinctive, anti-AI-slop.

## Steps (commit to the "Instrument" vision — see spec)
1. **Fonts**: pick an OFL display face (characterful technical grotesque/display-mono — NOT
   Inter/Roboto/Arial/system, NOT Space Grotesk) + a refined mono for data. Self-host woff2 in
   `assets/` + `@font-face` (preferred), or a `fonts.googleapis.com` `<link>` in index.html.
   Wire `--font-display` / `--mono` tokens; apply display to `.brand`/`h1`/`.stat .v`/numerals.
2. **Identity**: promote copper/substrate to brand/active/CTA/focus; demote teal to secondary.
   Do NOT recolor the board SVG. Re-verify contrast on any changed text token (≥4.5:1).
3. **Depth/atmosphere**: layered card shadows; a faint dot-grid/grain background (CSS
   gradient/SVG data-uri, cheap); a radial glow behind the active hero; subtle vignette.
4. **Board-as-hero**: elevate the Layout/Synthesis board/plot canvas (shadow + glow + size).
5. **Motion** (`@keyframes`): a staggered stage-content entrance (fade+translateY,
   `animation-delay` stagger across cards) on stage render; |S21| traces + board draw-in
   (`stroke-dashoffset` animation) on synthesis. Wrap ALL of it in `@media
   (prefers-reduced-motion: reduce)` → no/instant motion.
6. **Keep App.2.9**: responsive `@media` (no 390px overflow), AA contrast, aria/lang. Re-verify.

## Verify (FROM THE WORKTREE; `bash -c`; EXIT 0; quote output)
- `cd crates/yee-studio-web && dx build --platform web --release` → EXIT 0.
- `cargo test -p yee-studio-web` green; `cargo clippy -p yee-studio-web --all-targets -- -D
  warnings`; `cargo fmt --check`; `cargo check --workspace`.
- **Playwright** (load MCP; serve the built bundle via a bg static server over a `Yee/studio/`
  path so base_path resolves — the App.2.9 agent's approach): screenshot redesigned desktop
  (Synthesis + Layout board-hero + Technique, ≥1280px) + 390px mobile to
  `/home/hadassi/Code/Yee/.playwright-mcp/app-3-0/`; assert 390px `scrollWidth<=innerWidth`;
  confirm WASM loads + a technique click re-renders the board. Quote the overflow result +
  list the screenshot paths.
- `git -C /home/hadassi/Code/Yee status --porcelain crates/` EMPTY (main untouched).

Commit (worktree): `yee-studio-web: instrument redesign — copper identity, display type,
depth, motion (App.3.0, ADR-0152)` + Co-Authored-By trailer.

## Escape hatch
COMMIT to the distinctive vision (don't produce timid flat-dark-with-a-new-accent — the point
is a real aesthetic step-change). BUT the hard gates are non-negotiable: dx build EXIT 0, 15
tests green, NO functional regression, responsive (390px) + AA contrast + a11y + reduced-motion
all HELD. Do NOT touch `yee-filter`/engine. If a self-hosted font can't be fetched, use the
CDN `<link>` fallback (don't ship a broken font). If Playwright can't serve/load, ship + state
the gap honestly (don't fake a screenshot). Blocked > 60 min → surface + stop.

## Done when
The studio is visibly a distinctive "instrument" (copper identity, display type, depth,
entrance + draw-in motion) — NOT flat GitHub-dark — with the full pipeline working, App.2.9
responsive/AA/a11y held, reduced-motion honored, dx build EXIT 0, tests/clippy/fmt green,
screenshots captured; diff = `crates/yee-studio-web/**`. Then I (dispatcher) verify +
adversarial-review + judge the look + merge, then surface to the maintainer for final judgment.
