# ADR-0152: App.3.0 — Studio "Instrument" distinctiveness redesign

**Status:** Accepted
**Date:** 2026-06-01
**Related:** ADR-0151 (App.2.9 responsive/AA/a11y base this builds on), ADR-0110 (the original
design system this elevates), the 2026-06-01 frontend review, [[lumped-lc-and-studio-redesign]]

---

## Context

The 2026-06-01 frontend review (fresh-agent live walkthrough + source verification) found the
studio's substance distinctive but its chrome "AI-slop-adjacent": flat GitHub-dark surfaces,
default Tailwind teal, system fonts, zero motion, zero depth. The maintainer chose to fix the
P0/P1 defects first (App.2.9, shipped) **then** do a distinctiveness redesign, and approved
committing directly to the **Instrument/Workbench** direction.

## Decision

Redesign the studio chrome (`yee-studio-web`) to a distinctive "precision RF instrument"
aesthetic, committing fully (per the frontend-design skill): 

- **Identity from the board, not the teal** — promote copper `#e6b24d` + substrate-green to the
  brand signature (wordmark, active states, CTAs, focus/hero glows); demote the generic teal to
  a quiet secondary. The board SVG itself is unchanged — the chrome now echoes it.
- **Typographic voice** — a distinctive OFL display face (technical grotesque / display mono,
  NOT system / Inter / Roboto / Space Grotesk) for wordmark + headings + readout numerals,
  with a refined mono for data; self-hosted woff2 (CDN `<link>` fallback).
- **Depth + atmosphere** — layered elevation, a faint dot-grid/grain texture, a radial glow
  behind the active hero — off flat panels.
- **Board-as-hero** — the layout/synthesis board canvas becomes an elevated workbench surface.
- **One orchestrated motion moment** — a staggered stage-content entrance + |S21|/board
  draw-in on synthesis; gated behind `prefers-reduced-motion: reduce`.

## Consequences

**Ships:** a memorable, identity-driven studio true to what it is (a precision instrument that
outputs real fab files) — the review's highest-leverage change. **Hard gates held
(non-negotiable):** every App.2.9 gain preserved (responsive — no 390px overflow; AA contrast
re-verified on the new palette; a11y attrs + `lang`), `prefers-reduced-motion` honored, NO
functional regression (all 6 techniques / recommender / compare / board / BOM / tolerance /
verify / export; 15 tests green), engine untouched, `dx build` EXIT 0, fonts bundled so the
Pages deploy isn't broken. Distinctiveness is judged from captured screenshots by the
dispatcher + reviewer + ultimately the maintainer; if the look misses, it iterates.

**Not in scope:** engine/verdict logic, new features, a full motion system (one moment), the
EM-verify wall (ADR-0133/0147).

---

## References
- `crates/yee-studio-web/{assets/studio.css, assets/*.woff2, index.html, src/{stages.rs,
  main.rs, svg.rs}}`.
- `docs/superpowers/specs/2026-06-01-app-3-0-studio-instrument-redesign-design.md`;
  `docs/superpowers/plans/2026-06-01-app-3-0-studio-instrument-redesign.md`.
- Frontend review (2026-06-01); App.2.9 / ADR-0151.
