# ADR-0111: Filter Phase F2.0 — lumped-element LC ladder synthesis

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0084 (synthesis core / g-values), ADR-0097/0109 (distributed
dimensional synthesis — the pattern), `FILTER-DESIGN-ROADMAP.md`,
[[project-filter-design-final-goal]]

---

## Context

New product goal (2026-05-30): a **full lumped-LC filter → PCB** flow with a
polished UI, component choosing, EM simulation, BOM, and tolerance analysis. The
lumped-LC technology track is greenfield — `yee-synth`/`yee-filter` produce the
abstract prototype (g-values, coupling matrix) but nothing realizes a **lumped
L/C** filter. (Distributed edge-coupled/hairpin are done; `yee-fdtd` already has
`LumpedRlcPort::series_rlc` for the later EM step.)

## Decision

Add closed-form **LC ladder synthesis** to `yee-filter` (`src/lumped.rs`):
`synthesize_lumped(&FilterProject) -> LumpedLadder` applying the textbook
low-pass-prototype → band-pass ladder transform (Pozar §8.3) to map each g-value
to a series or shunt **LC resonator** tuned to f0 (`L_k`, `C_k` from g_k, ω0, Δ,
Z0; alternating series/shunt). A private ABCD-cascade `ladder_s21` computes the
realized response. Pure `f64`, WASM-safe, no FDTD/parts/PCB.

Gate `lumped_001`: synthesize the committed Chebyshev N=5 fixture, assert N
resonators each tuned (`L·C·ω0² ≈ 1`), and the ladder `|S21|` (ABCD cascade)
**meets the same spec mask** as the synthesized prototype — the LC realization
reproduces the design (self-consistent + textbook-transform benchmark).

## Consequences

**Ships:** the lumped-track foundation — ideal L/C element values + a reference
realized response. Feeds F2.1 (map ideal → real E-series/vendor parts + ESR/SRF
→ BOM), F2.3 (voxelize + `LumpedRlcPort` → FDTD S-params), F2.4 (Monte-Carlo over
part tolerances → yield), and the Dioxus lumped UI track.

**Gate:** `cargo test -p yee-filter` green incl. `lumped_001` (≤ the spec mask).
Pure-math, sub-second.

**Not in scope:** parts/parasitics (F2.1), PCB/footprints (F2.2), FDTD lumped sim
(F2.3), tolerance (F2.4), UI, and non-band-pass transforms.

---

## References
- `docs/superpowers/specs/2026-05-30-f2-0-lumped-lc-synthesis-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-0-lumped-lc-synthesis.md`.
- Pozar, *Microwave Engineering* §8.3 (lumped-element filter transforms);
  Hong & Lancaster ch. 3.
