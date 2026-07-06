# ADR-0181: S.3b 3-D field surface + S.4 Dioxus parity audit (retirement deferred)

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0179/0180 (studio track), ADR-0175 (Dioxus feature-freeze), ADR-0110/0130
(the Dioxus studio lineage).

---

## S.3b — 3-D field surface (three.js)

The engine's field slice now renders as a height-mapped, vertex-colored 3-D mesh with orbit
controls (`FieldSurface3D`). Design points:

- **Geometry is a pure function** (`buildSurface`: positions/colors/indices from a slice),
  vitest-gated against hand-computable values (vertex/triangle counts, ±heightScale extremes,
  grid centering, index validity) — the same keep-the-math-testable discipline as
  `analysis.ts`.
- **WebGL fallback**: without a WebGL context the component renders a text notice — which is
  exactly what jsdom provides, so the DOM gate exercises the fallback path for real.
- **Code-split**: three.js rides a lazy chunk (133.3 kB gz) fetched only when a result is on
  screen; the initial bundle stays **49.4 kB gzipped**. 11 vitest tests green.

## S.4 — parity audit: `yee-studio-web` (Dioxus) vs the Tauri studio

| Capability | Dioxus (`yee-studio-web`) | Tauri studio |
|---|---|---|
| Filter spec → synthesis → response/mask | ✔ (5 topologies incl. lumped-LC) | ✖ |
| Dimensional synthesis + SVG layout | ✔ | ✖ |
| Export: Gerber, KiCad, `.s2p`, JLCPCB set | ✔ | ✖ |
| Deployed to GitHub Pages | ✔ (`/Yee/studio/`) | ✖ (desktop app) |
| Full-wave EM simulation (GPU/CPU engine) | ✖ (EM deliberately out of dep graph, ADR-0089) | ✔ (yee-engine jobs, progress, cancel) |
| Field visualization (probe/spectrum/slice/3-D) | ✖ | ✔ |
| Remote execution path | ✖ | ✔ (same protocol over `yee-server` WS) |

**Verdict: not at parity, and parity is the wrong frame today.** The two studios serve
disjoint roles: Dioxus is the shipped *filter-design* product; the Tauri studio is the
*engine* product. Retiring Dioxus now would delete a working, deployed designer with no
replacement.

**Decision.**
1. The ADR-0175 freeze stands: `yee-studio-web` stays deployed and feature-frozen.
2. **Convergence path** (the point of the engine track all along): the filter studio's Verify
   stage is the ideal coupling-matrix response because ADR-0089 kept EM out of the WASM dep
   graph — and the missing "native server the web client calls" now exists (`yee-server`,
   ADR-0180). The planned convergence is *engine-powered verify* — the filter flow (in
   whichever frontend) submits full-wave jobs over the S.0 protocol — rather than a
   line-by-line port of 7.3 k lines of Dioxus UI.
3. Retirement is re-decided (own ADR) when the Tauri/React studio hosts the filter flow or
   the filter flow consumes engine jobs, whichever lands first.

S.4 closes as **AUDITED — retirement deferred with a defined convergence path**. The studio
track has no other queued phases; next studio-track work items are engine-powered verify
(S.5, new) and interactive/live-streamed visualization on top of the WS transport.
