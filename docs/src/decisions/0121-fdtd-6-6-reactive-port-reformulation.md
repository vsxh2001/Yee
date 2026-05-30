# ADR-0121: Phase 2.fdtd.6.6 — reactive lumped-port reformulation (sheet→mode coupling), bench-validated

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0119 (the de-embed bench + PORT-WRONG verdict — this is the
research track's increment 2), ADR-0118 (canonical per-element updates, verified
per-edge), ADR-0115 (F2.3, blocked), the lumped-LC → PCB goal (maintainer chose
"open the reactive-port research track"), [[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0119's V+I de-embedding bench (`reactive_deembed_001`) is shipped and honest
(asserted resistor anchor). It returned **PORT-WRONG**, with the **well-conditioned
capacitor** as the load-bearing evidence: the canonical shunt capacitor presents
`Z_in ≈ 94 Ω` where `κ/(ωC) ≈ 3175 Ω` is expected — **over-coupled** by ~N×. Yet
the per-element constitutive is verified **correct per-edge** (ADR-0118: a forced
single edge gives `L → +488j`, `C → −496j`). The discrepancy is therefore in the
**sheet → guide-mode coupling**, not the constitutive law:

The full-width lumped *sheet* places the element on every transverse interior
`E_z` edge. If each of the `N` cells carries the **full** element value, the sheet
sums to `N ×` the intended admittance (a shunt capacitor → `N·C` → over-short; a
shunt inductor → `L/N`), so the port loads the line by the wrong amount. The
modal voltage `V = ∫E·dz` (a line integral) and modal current `I = ∮H·dl` (the
whole cross-section) relate to the per-cell injected current with a geometric /
normalization factor that the resistor's instantaneous `κ`-calibration absorbs but
the reactive arms do not.

Now that increment 1 gives a fast, honest `Z_L(ω)` read-out, the reformulation is
**bench-iterable**, not open-ended.

## Decision

Reformulate the reactive lumped-port → guide-mode coupling so the de-embedded
`Z_L(ω)` matches `R + jωL + 1/(jωC)` (the per-edge-correct constitutive), within a
loose tol, on `reactive_deembed_001`. Concrete hypotheses, in order, each
validated against the bench:

1. **Sheet value-normalization** — distribute the lumped element across the `N`
   sheet cells (a shunt `C` → `C/N` per cell, a shunt `L` → `N·L` per cell; series
   topologies dual), so the sheet sums to the intended `Z_L`. Cheapest; the
   well-conditioned capacitor over-coupling points straight at it.
2. **Modal coupling factor** — if (1) is insufficient, the per-cell field
   back-action `(dt/(ε₀·dA))·I` must use the modal cross-section / a port-face
   normalization consistent with `V = ∫E·dz` and `I = ∮H·dl`, not the bare
   single-cell `dA`. Derive the factor from the resistor's measured `κ`.
3. **Modal lumped port** — if needed, a proper 1-port termination that enforces
   `V = Z_L·I` on the *measured modal* V and I (the bench's own quantities) with
   the current distributed over the port face ∝ the mode.

Strengthen `reactive_deembed_001` to **assert** the reactive arms (shunt-C first —
well-conditioned; shunt-L and series-RLC as the de-embed conditioning allows) once
they match. The resistor anchor stays asserted; never weakened, never faked.

## Consequences

**Ships (target):** a reactive lumped-port whose de-embedded `Z_L(ω)` is correct —
the foundation that **unblocks F2.3** (the lumped-filter selectivity, the goal's
EM-sim component). The bench flips from PORT-WRONG to PORT-CORRECT, asserted.

**Gate:** `reactive_deembed_001` GREEN with the reactive-arm assertions (shunt-C at
minimum); resistor anchor + one-way + fdtd-206 non-regressed.

**Escape hatch (honest partial):** if hypotheses (1)→(3) do not bring the
well-conditioned capacitor within tol while staying stable, surface the bench
`Z_L(ω)` after each attempt and the precise coupling tried; the bench quantifies
the residual. A measured partial (e.g. "C within tol, L still ill-conditioned") is
acceptable and recorded — never weaken the anchor or fake a pass. If genuinely
intractable after the increment, the multi-week modal-port (3) is the remaining
scope.

**Not in scope:** F2.3's board sim (rides on this once the port is correct — a
follow-on); F2.3's own element placement (re-checked after the port is fixed);
SRF/ESR parasitics.

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-6-reactive-port-reformulation-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-6-reactive-port-reformulation.md`.
- ADR-0119 bench + verdict; ADR-0118 canonical per-edge constitutive.
- Taflove & Hagness, lumped-element FDTD; modal/lumped-port de-embedding.
