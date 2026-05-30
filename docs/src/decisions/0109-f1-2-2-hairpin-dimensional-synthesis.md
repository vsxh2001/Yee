# ADR-0109: Filter Phase F1.2.2 — hairpin dimensional synthesis

**Status:** Accepted (design; implementation queued)
**Date:** 2026-05-30
**Related:** ADR-0097 (F1.2.0 edge-coupled dimensional synthesis — the pattern),
ADR-0094 (coupled-line model), `FILTER-DESIGN-ROADMAP.md`,
[[project-filter-design-final-goal]]

---

## Context

The product goal names **topologies** (plural), but only the edge-coupled
half-wave topology is dimensioned today (`dimension_edge_coupled`, ADR-0097). The
`yee_layout::hairpin_bpf` geometry generator already exists (Hong & Lancaster
ch. 6), with no synthesis driving it. A hairpin is a half-wave line folded into a
U; adjacent hairpins couple through the edge gap between their adjacent arms —
the **same** edge-coupled mechanism `dimension_edge_coupled` already inverts.

## Decision

Add closed-form **hairpin dimensional synthesis** to `yee-filter`
(`dimension_hairpin` / `dimension_hairpin_layout` + `HairpinDimensions`),
mirroring F1.2.0:

- Line width = `microstrip_width(z0)`; ε_eff from `eps_eff`.
- Resonator arm length = `λ_g/4 = c/(4·f0·√ε_eff)` (the U-folded half-wave
  resonator is two ≈λ/4 arms — the factor-4 vs edge-coupled's factor-2 λ/2).
- Inter-resonator gaps = the SAME gap-bisection as edge-coupled: solve
  `coupling_coefficient(coupled_microstrip(w, s, h, εr)) == target_k[i]` with
  `target_k[i] = fbw·m[i][i+1]`, bisection (monotone gap→k, `coupled_002`).
- Layout via `hairpin_bpf`; the single-`coupling_gap_m` vs per-section-gap
  mismatch is resolved either by a minimal backward-compatible per-section
  extension to `HairpinParams` or a documented uniform-gap walking-skeleton
  limitation (implementation picks, surfaced in the report).

Pure `f64`, WASM-safe, NO FDTD/surrogate — the initial dimensioning F1.2.1 BO
refines. `dimension_edge_coupled` is unchanged.

## Consequences

**Ships:** a SECOND filter topology in the synthesis→dimensions→layout pipeline —
the goal's "topologies" deepened beyond edge-coupled.

**Gate:** `hairpin_dim_001` (mirror of `dim_001`): synthesize the Chebyshev N=5
fixture, `dimension_hairpin` on FR-4, assert each realized gap re-evaluates to its
`target_k` within < 1 %, `arm_length ≈ λ_g/4`, `gaps.len() == N-1`. Pure-math,
sub-second; `cargo test -p yee-filter` green.

**Not in scope:** FDTD validation of the hairpin; tapped-feed Qe synthesis;
combline/interdigital (need shorted-resonator + via models); CLI/studio
`--topology hairpin` wiring (a cross-lane follow-on).

---

## References
- ADR-0097 (the edge-coupled synthesis this mirrors); ADR-0094 (coupled-line k).
- `docs/superpowers/specs/2026-05-30-f1-2-2-hairpin-dimensional-synthesis-design.md`;
  `docs/superpowers/plans/2026-05-30-f1-2-2-hairpin-dimensional-synthesis.md`.
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, ch. 6
  (hairpin) / ch. 8 (coupling-matrix synthesis).
