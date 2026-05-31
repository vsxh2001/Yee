# ADR-0144: Filter F1.2.5 — Combline dimensional synthesis

**Status:** Accepted
**Date:** 2026-05-31
**Related:** the maintainer pick (AskUserQuestion 2026-05-31 — combline, with the
proper-gate caveat), the combline research sweep (workflow `waws53n82`, confidence 0.9),
ADR-0109/0138 (hairpin — the coupled-resonator dimensioning pattern reused), ADR-0136
(the recommender already recommends combline), ADR-0142/0143 (compare/overlay — combline
will join them once lit), [[lumped-lc-and-studio-redesign]]

---

## Context

The maintainer chose combline as the next technique, explicitly conditioning it on a
PROPER published-design gate (not the shallow λg/4 mirror that left interdigital
gate-blocked). A 3-source research sweep (Hong & Lancaster §5.2.5, Pozar, Matthaei) both
confirmed the synthesis and — crucially — identified the **non-tautological** gate:

- The loading cap `C_L = cot(θ0)/(2π·f0·Z0)` (short-circuited θ0<90° line + shunt cap →
  resonance at f0) is correct, **but gating on it is tautological** (it is the engine's
  own emit formula).
- The proper gate is **H&L eq (5.46)**: a *published* 5-pole 0.1 dB Chebyshev combline
  design (g=[1.1468,1.3712,1.9750,1.3712,1.1468], FBW=0.1) → external Q **Qe=11.468** and
  inter-resonator couplings **M₁₂=0.07975, M₂₃=0.06077** (reproduced to machine precision
  from the g-values; a second FBW=0.15 pseudocombline point Qe=7.645/M₁₂=0.11962/
  M₂₃=0.09115). Plus a first-principles resonance check (independently root-find the
  loaded-stub susceptance zero → f0), which catches a wrong cap/length/dispersion without
  inverting the cap formula.

## Decision

Add `dimension_combline` to `yee-filter::dimension`, mirroring `dimension_hairpin`:

- **Reuse** the validated coupling realization (`target_k = FBW·m[i][i+1]` → `solve_gap`
  over `coupled_microstrip`/`coupling_coefficient`) — identical to edge-coupled/hairpin.
- **Combline-distinct:** a short-circuited resonator of electrical length `θ0` (default
  45° = λg/8 for compactness), physical length `θ0/β(f0)`, and a loading cap
  `C_L = cot(θ0)/(2π·f0·Z0)` at the open end (the other end is a via to ground).
- Gate `dim_combline_001`: the H&L eq 5.46 published Qe/M benchmark (non-tautological) +
  the first-principles resonance consistency check + coupling-solved/positive-dims.

Engine first; the **studio lighting** of combline is a follow-on (the stepped-Z /
hairpin engine→studio pattern).

## Consequences

**Ships:** combline dimensional synthesis — the compact, high-Q, narrow-band band-pass
realization the recommender already points at — with a PROPER published-benchmark gate
(H&L eq 5.46), honestly satisfying the maintainer's condition. Reuses the proven coupling
machinery; adds the combline-distinct short-circuited θ0 resonator + loading cap.

**Gate (non-vacuous, published):** `dim_combline_001` — synthesized Qe/M match H&L's
published combline-design numbers (a wrong/constant synthesis fails by ≫tol); the
dimensioned resonator+cap resonate at f0 by independent root-find. NOT the cap-formula
tautology (the research explicitly flagged that trap).

**Honest scope:** this first-order engine reuses the `solve_gap` coupling realization
(like hairpin) rather than the rigorous Getsinger/Cristal self-/mutual-capacitance
coupled-bar synthesis (H&L eq 5.44) — documented. Discrete E-series `C_L` selection,
via/short-circuit 3-D modelling, and the studio lighting are out of scope. The EM-verify
wall (ADR-0133) is untouched.

---

## References
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications* §5.2.5–5.2.6
  (eqs 5.42–5.48); Matthaei-Young-Jones Ch. 8.
- `crates/yee-filter/src/dimension.rs` (`dimension_hairpin`); `crates/yee-synth` (`prototype`).
- `docs/superpowers/specs/2026-05-31-f1-2-5-combline-synthesis-design.md`;
  `docs/superpowers/plans/2026-05-31-f1-2-5-combline-synthesis.md`; research workflow `waws53n82`.
