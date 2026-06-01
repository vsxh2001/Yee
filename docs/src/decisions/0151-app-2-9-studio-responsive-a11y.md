# ADR-0151: App.2.9 — Studio responsive + accessibility + P1 polish

**Status:** Accepted
**Date:** 2026-06-01
**Related:** ADR-0110 (the studio design system this fixes within), ADR-0130/0135 (the Dioxus
studio, deployed on Pages), the 2026-06-01 frontend review; ADR-0152 (App.3.0 distinctiveness
redesign — the follow-on), [[lumped-lc-and-studio-redesign]]

---

## Context

A fresh-agent live walkthrough of the deployed studio (`yee-studio-web`) + source verification
found real, user-facing defects: **no responsive layout** (`0 @media` → broken at ≤420px),
**WCAG-AA contrast failures** (`--muted #6b7480` ≈ 3.6:1, used for every caption; 9px rail
labels), and several P1s (unstyled recommender "Use this" CTA, no `lang`, color-only active
stage, unlabeled decorative glyphs, under-differentiated selected card, a misleading
"Spec PASS" shown above "yield 0.0%"). The maintainer chose to fix these **then** do a
distinctiveness redesign — this ADR is the fixes (Increment A); the redesign is ADR-0152.

## Decision

Fix the verified P0/P1 defects **within the existing design system** (no aesthetic overhaul):

- **Responsive** — add `@media` breakpoints (≤760px, ≤420px): the fixed stage rail becomes a
  horizontal top strip, multi-column grids/tables collapse or scroll, inputs/buttons fit the
  viewport — **no horizontal overflow at 390px**; desktop unchanged.
- **Contrast** — raise `--muted` to clear **AA 4.5:1** on `#0b0d11`; rail labels 9px → ≥ 11px;
  palette identity (teal/copper/substrate) kept.
- **P1** — `.rec-use` styled as an accent CTA; `lang="en"` on the root; `aria-current="step"`
  + accessible names on the rail + `aria-hidden` decorative glyphs; selected technique card
  filled-accent (not border-only); the Export verdict labelled **nominal** with a low-yield
  caveat (wording/markup only, no engine change).

## Consequences

**Ships:** an accessible, mobile-usable studio on a sound base — broken-mobile + AA-contrast
are correctness defects, not polish, so this is a genuine quality fix to the deployed app.
Gated by `dx build` EXIT 0 + a Playwright 390px no-overflow re-check (before/after) + the
contrast ratio stated + the a11y attrs present + desktop unchanged + tests/clippy/fmt green.
Low-risk, design-system-internal; unblocks the App.3.0 redesign on an accessible foundation.

**Not in scope:** the distinctiveness redesign (ADR-0152 — board-as-hero, copper brand
identity, display font, depth/atmosphere, motion); engine/verdict logic (the yield-vs-PASS
fix is wording only); a full a11y audit beyond the review's findings.

---

## References
- `crates/yee-studio-web/{assets/studio.css, src/stages.rs, src/main.rs}`.
- `docs/superpowers/specs/2026-06-01-app-2-9-studio-responsive-a11y-design.md`;
  `docs/superpowers/plans/2026-06-01-app-2-9-studio-responsive-a11y.md`.
- Frontend review (fresh-agent Playwright walkthrough + source verification), 2026-06-01.
