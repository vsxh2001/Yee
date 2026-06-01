# App.2.9 — Studio responsive + accessibility + P1 polish — Design Spec

**ADR:** ADR-0151 · **Date:** 2026-06-01 · **Status:** Accepted
**Follows:** the frontend review (fresh-agent live walkthrough + source verification, 2026-06-01)
of the deployed studio (`yee-studio-web`, ADR-0110 design system). Increment A of the
maintainer-chosen "P0 fixes then redesign" pair (the distinctiveness redesign is App.3.0 /
ADR-0152, a follow-on).

## Problem (verified defects in the deployed studio)

The review found, and source confirms:
- **P0 — no responsive layout** (`0 @media` queries in `studio.css`): at ≤420px the content
  overflows horizontally — clipped number inputs/`<select>`, the brand wraps to 3 lines, the
  fixed 78px rail eats the viewport. Broken, not degraded.
- **P0 — WCAG-AA contrast failures**: `--muted #6b7480` on `--bg #0b0d11` ≈ 3.6:1 (AA needs
  4.5:1 for normal text) and it styles *every* caption at 13px; rail labels are 9px.
- **P1s**: the recommender "Use this" CTA is unstyled (`.rec-use { cursor: pointer }` only →
  default light browser button, breaks the dark theme); no `lang` on `<html>`; the active
  stage is conveyed by color only (no `aria-current`); decorative rail glyphs are unlabeled
  (read aloud as "◈ Spec"); the selected technique card is under-differentiated (border-only,
  same teal badge as others); the Export "Design summary" shows "Spec verdict PASS" directly
  above "yield 0.0%" with no nominal-vs-yield distinction (misleading to a novice).

## Scope (Increment A — the verified P0/P1 fixes; NOT the redesign)

Fix the defects within the existing design system (no aesthetic overhaul — that's App.3.0):

1. **Responsive (`studio.css` `@media`)** — add breakpoints so nothing overflows at 390px:
   - `≤ 760px`: stage rail becomes a **horizontal top tab-strip** (or wraps) instead of a
     fixed 78px column; `.grid-2` / `.grid-3` / `.tgrid` / `.compare-table` wrapper collapse
     to single column / horizontal-scroll where a table must stay wide; `.canvas` padding
     reduced; the top bar wraps gracefully (brand stays one line; spec-chip drops if needed).
   - `≤ 420px`: number inputs / `<select>` shrink to fit (`max-width: 100%`); the export-row
     buttons wrap; no element exceeds the viewport width (no horizontal scrollbar on the page).
   - Desktop layout unchanged above the breakpoint.
2. **Contrast (tokens)** — raise `--muted` to a value that clears **AA 4.5:1** on `--bg`
   (e.g. `#8b95a1` is ~5.0:1 — verify the exact ratio and pick the lightest that still reads
   as "muted"); reuse `--muted-2` where appropriate. Raise the rail `.lab` font-size from 9px
   to **≥ 11px**. Keep the palette identity (teal/copper/substrate) intact.
3. **P1 fixes**:
   - Style `.rec-use` as a proper accent CTA (mirror `.export-row .btn.dl` hover or the
     `.rec-primary` accent treatment) — no default white button.
   - `lang="en"` on the document root (Dioxus index template / `Dioxus.toml`, or the app's
     root element — whichever the dx build honors).
   - `aria-current="step"` on the active rail item; an accessible name on each rail button
     (`aria-label` = the stage name); `aria-hidden="true"` on the decorative glyph spans.
   - Selected technique card: add a filled-accent badge / inset glow (not border-only) so the
     active card is unmistakable.
   - Export "Design summary": label the spec verdict as **nominal** (spec-mask) and surface a
     yield caveat chip when MC yield is low, so "PASS" + "0% yield" no longer read as a
     contradiction. (Wording/markup only — no engine change.)

## DoD (machine-checkable)

1. **Responsive**: `studio.css` contains `@media` breakpoints; a Playwright re-check at 390px
   shows **no horizontal page overflow** and no clipped inputs (before/after screenshots);
   desktop (≥1000px) visually unchanged.
2. **Contrast**: the new `--muted` computes to **≥ 4.5:1** on `#0b0d11` (state the ratio); rail
   `.lab` ≥ 11px.
3. **A11y**: `lang` present on the root; `aria-current="step"` on the active rail item;
   rail buttons have accessible names; decorative glyphs `aria-hidden`.
4. **Build/tests**: `dx build --platform web --release` EXIT 0; `cargo test -p yee-studio-web`
   green; `cargo clippy … -D warnings` + `cargo fmt --check` clean; `cargo check --workspace`.
5. No engine/`yee-filter` change; diff confined to `crates/yee-studio-web/**`.

## Changes

- `crates/yee-studio-web/assets/studio.css` (responsive `@media` + contrast tokens + `.rec-use`
  + selected-card + the rail-label size),
- `crates/yee-studio-web/src/{stages.rs, main.rs}` (aria attrs, `lang`, the yield-caveat
  wording, the selected-card markup),
- `crates/yee-studio-web/Dioxus.toml` only if that's where `lang` must be set.
  NO `yee-filter` edits.

## Out of scope

The distinctiveness redesign (App.3.0 / ADR-0152 — board-as-hero, copper brand identity,
display font, depth/atmosphere, motion). A full a11y audit beyond the review's findings.
Engine/verdict logic (the yield-vs-PASS fix is wording only).

## Why

Fixes real, verified, user-facing defects in the *deployed* studio — a broken mobile layout
and AA-failing contrast are correctness issues, not polish. Clean, low-risk, within the
existing design system; ships fast and unblocks the redesign on a sound, accessible base.
