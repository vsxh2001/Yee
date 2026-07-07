# A.0 — Patch antenna on the engine (walking skeleton): synthesis + S11 resonance

**Date:** 2026-07-06
**Phase:** A.0 — first phase of the antenna track (Part 3 of ENGINE-STUDIO-ROADMAP),
opened per the project goal "design an antenna or filter with the engine". Inherits the
S.9–S.12 measurement stack.
**Plan:** `docs/superpowers/plans/2026-07-06-a0-patch-antenna.md`

## Design

- **Synthesis (closed form, Balanis §14.2)** in `yee-layout` (the geometry crate — pure
  math, WASM-safe): `patch_antenna_dims(f0, ε_r, h)` → width `W = c/(2f₀)·√(2/(ε_r+1))`,
  `ε_eff(W, h)` (existing Hammerstad form), open-end extension ΔL (Hammerstad — same
  formula the stub gate uses), length `L = c/(2f₀√ε_eff) − 2ΔL`. Plus
  `edge_fed_patch(f0, substrate, z0_feed)` → a two-rect `Layout` (feed line + patch,
  edge-fed at the patch's radiating edge centre, one 50 Ω `PortRef`).
- **Verify (engine)**: gate `engine-antenna-001` in `yee-engine/tests` (dev-deps
  already present): voxelize, S.5 materials job, CPML-xy walls, S.10 aperture drive
  port, S.7 two-run |S11| (reference = bare feed on the same bbox). The patch resonates
  where the cavity model says: assert the |S11| **dip frequency** within **±10 %** of
  the designed f₀ = 2.45 GHz. The dip position is robust even though an edge-fed patch
  is badly mismatched (edge resistance ~200–300 Ω → shallow but localized dip); depth
  is asserted loosely (≥ 2 dB below the band median) — matching/inset feed is A.1.

## Non-goals

Radiation pattern over the protocol (needs NTFF exposure + an open z-top — per-face
CPML, its own increment); input-impedance matching (inset feed, A.1); gain/efficiency.
