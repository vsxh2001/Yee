# App.3.0 — Studio "Instrument" distinctiveness redesign — Design Spec

**ADR:** ADR-0152 · **Date:** 2026-06-01 · **Status:** Accepted
**Follows:** the 2026-06-01 frontend review + App.2.9 (ADR-0151, responsive/AA/a11y base).
Increment B of the maintainer-chosen "P0 fixes then redesign" pair. Maintainer approved
"commit + redesign directly" on the **Instrument/Workbench** direction.

## Problem

The review's verdict: the studio's *substance* is distinctive (live WASM EDA pipeline, real
Gerber/KiCad, the copper-on-substrate board, honest ideal-vs-realized), but its *chrome* is
"AI-slop-adjacent" — flat GitHub-dark surfaces, default Tailwind teal, **system fonts**, zero
motion, zero depth/atmosphere. Highest-leverage fix: stop treating the board as one panel
among many; make the copper/substrate identity + the board the hero, give it a typographic
voice, depth, and one orchestrated motion moment.

## Aesthetic direction (commit fully — "Instrument, not dashboard")

Yee is a precision RF **instrument**, not a generic dark SaaS dashboard. Execute with
intentionality (the frontend-design skill: distinctive, production-grade, memorable — never
generic). Concretely:

1. **Identity from the board, not the teal.** Promote **copper `#e6b24d` + substrate-green
   `#3f9e72`** to the brand signature — the wordmark, active stage, primary CTAs, hero
   accents, focus glows. Demote the generic Tailwind-teal `#2dd4bf` to a quiet secondary/info
   role (or retune it). **Do NOT recolor the board SVG itself** (copper traces / green
   substrate stay — the chrome now *echoes* the board, unifying the identity).
2. **Typographic voice.** Introduce a distinctive **display face** (OFL/open-licensed) for the
   wordmark + H1 + the stat/readout numerals — a characterful technical grotesque or display
   mono. **NOT** Inter/Roboto/Arial/system, and **NOT** Space Grotesk (overused). Keep/upgrade
   a refined **mono** for engineering data (the "instrument readout"). Self-host the woff2 in
   `assets/` (`@font-face`) for an offline/self-contained deploy; a `fonts.googleapis.com`
   `<link>` in `index.html` is an acceptable fallback. Body may stay the system stack or
   upgrade to a refined sans.
3. **Depth + atmosphere.** Move off flat panels: layered elevation shadows on cards, a faint
   **dot-grid or grain/noise background texture**, a subtle **radial glow** behind the active
   stage's hero element, a gentle vignette. A "lit instrument panel," not flat fill.
4. **Board-as-hero.** On the Layout (and Synthesis) stages, the board/plot canvas becomes a
   larger, **elevated** workbench surface (depth shadow + glow + a touch more prominence) —
   the visual anchor, not an equal panel.
5. **One orchestrated motion moment** (NOT scattered micro-interactions): a staggered
   stage-content **entrance** (~150–250ms fade + small translate, `@keyframes` +
   `animation-delay` stagger) on stage change, and the **|S21| plot traces + the board
   drawing in** (SVG `stroke-dashoffset` / opacity animation) when synthesis updates. Fill the
   data-light stages' empty canvas with atmosphere rather than void.

## Hard constraints (non-negotiable — gates)

- **Keep every App.2.9 gain:** responsive `@media` (NO horizontal overflow at 390px), the
  **AA-contrast** floor (the NEW palette/text must STILL compute ≥ 4.5:1 for body/caption text
  — re-verify; copper/substrate on dark can be low-contrast, so check), the a11y attrs
  (`aria-current`, `aria-label`, `aria-hidden`, `lang`).
- **`prefers-reduced-motion: reduce`** — disable/curtail all the new motion (a11y requirement).
- **No functional regression:** all 6 techniques, recommender, compare/overlay, board, BOM,
  tolerance, verify, export must work exactly as before; `cargo test -p yee-studio-web` (15)
  green; the engine (`yee-filter`) is NOT touched.
- **`dx build --platform web --release` EXIT 0**; the live app still loads (WASM mounts).
- Self-hosted fonts must be bundled so the Pages deploy serves them (no broken/blocked font).

## DoD (machine-checkable + judged)

1. `dx build --platform web --release` EXIT 0; `cargo test -p yee-studio-web` green;
   `cargo clippy … -D warnings` + `cargo fmt --check`; `cargo check --workspace`.
2. **Playwright visual capture + regression re-check:** load the built bundle; screenshot the
   redesigned **desktop** (≥1280px) Synthesis + Layout (board-hero) + Technique stages AND the
   **390px** mobile to `/home/hadassi/Code/Yee/.playwright-mcp/app-3-0/`; assert 390px still
   has NO horizontal overflow (`scrollWidth<=innerWidth`); confirm the WASM app loads + a
   technique click still re-renders the board.
3. **Contrast re-verify**: body/caption text colors compute ≥ 4.5:1 on their backgrounds
   (state the ratios for any changed text token).
4. **Reduced-motion**: a `@media (prefers-reduced-motion: reduce)` block neutralizes the new
   animations.
5. Distinctiveness (judged by dispatcher + reviewer + maintainer from the screenshots): the
   chrome is visibly NOT flat-GitHub-dark — the copper/instrument identity, the display font,
   depth, and the entrance/draw-in motion are present and cohesive.

## Changes

- `crates/yee-studio-web/assets/studio.css` (the redesign: tokens, depth, texture, motion
  `@keyframes`, board-hero, reduced-motion), `assets/*.woff2` (self-hosted fonts, if used),
  `index.html` (font `<link>`/preload + keep `lang`), `src/{stages.rs, main.rs, svg.rs}`
  (markup for the entrance stagger classes / board-hero wrapper / draw-in SVG attrs / wordmark).
  NO `yee-filter` / engine edits.

## Out of scope

Engine/verdict logic; new features/stages; a full motion system (one orchestrated moment, not
many). The EM-verify wall (ADR-0133/0147).

## Why

Turns the studio from "competent but generic dark tool" into something *memorable* + true to
what it is (a precision RF instrument that outputs real fab files) — the review's
highest-leverage change — without regressing the functional pipeline or the just-won
responsive/AA/a11y base.
