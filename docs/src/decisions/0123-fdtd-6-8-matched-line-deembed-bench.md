# ADR-0123: Phase 2.fdtd.6.8 — matched-line reactive de-embed bench

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0122 (per-axis CPML — the enabling brick 1, just shipped),
ADR-0121/0119 (the PEC-source de-embed bench + its ≈0.37 PORT-WRONG verdict and
its inability to do a long-window + echo-free + clean-anchor reactive measurement),
the lumped-LC → PCB goal (maintainer-approved reactive-port research track),
[[project-lumped-lc-and-studio-redesign]]

---

## Context

The increment-1/2 de-embed bench (ADR-0119/0121) measured the reactive lumped
port at ≈0.37 off, but on a **PEC-terminated** line: its `Z₀`/κ calibration is
tied to the *reflecting* source, so it cannot give a simultaneously long-window
(full dispersive tail), echo-free, and clean-anchor reactive measurement —
absorbing the source end broke the anchor. That ambiguity is why the ≈0.37
residual is only *bracketed*, not pinned, and why "PORT-CORRECT vs single-cell
limit" is unresolved.

Brick 1 (ADR-0122, per-axis CPML `with_axes`) now lets us build a **matched
line**: a parallel-plate guide with **x-only CPML at both x-ends** (absorbing — no
source-end echo, no far-wall echo) and **PEC on the transverse walls** (guide mode
preserved). On a matched line the standard incident/reflected de-embed applies
with **no multiple bounces**, so the PEC-source κ hack is unnecessary.

## Decision

Add a matched-line reactive de-embed bench (gate `reactive_deembed_matched_001`,
yee-fdtd) that pins the reactive port's `Z_L(ω)` definitively:

- **Guide:** parallel-plate, x-only CPML both ends (`CpmlParams::for_grid(..).
  with_axes([true,false,false])`), PEC transverse walls. A soft `E_z` source
  launches +x; a full-width lumped load sits at a load plane; a reference plane
  sits between source and load.
- **De-embed (matched, standard VNA):** a **reference run** (no load / matched)
  gives the incident `V_inc(ω)`, `I_inc(ω)` at the reference plane (V = ∫E·dz,
  I = ∮H·dl) and the measured line `Z₀(ω) = V_inc/I_inc`. A **load run** gives the
  total `V`, `I`; the reflected wave `= total − incident` (clean two-run
  difference — the CPML absorbs it after one pass, no bounces); `Γ = (Z_in−Z₀)/
  (Z_in+Z₀)`, `Z_in = V/I`, back out `Z_L(ω)`. No scalar-`A` / κ calibration tied
  to a reflecting source — the matched line makes incident/reflected separation
  honest on its own.
- **Honesty anchor (asserted):** a known resistor de-embeds to `Z_L → R` within a
  loose tol. The reactive arms are then measured with a **long** window (full
  dispersive tail, since nothing bounces) AND **echo-free** AND **clean anchor** —
  the measurement the PEC-source bench could not give.

**Outcome gate (decision):** if the well-conditioned capacitor `Z_L(ω)` matches
`1/(jωC)` within the loose tol (sub-`react_tol`) → **PORT-CORRECT**: the reactive
lumped port is fine; the EM-sim blocker was the *measurement*, and F2.3 should be
re-run (its flatness then traces to element placement, a smaller fix) — assert it.
If it still sits at ≈0.37 with the clean bench → the **single-cell ε_eff limit is
confirmed real** → brick 3 (multi-cell aperture port). Either way, record the
pinned number; never weaken the anchor, never fake.

## Consequences

**Ships:** a trustworthy matched-line reactive de-embed bench + a *pinned* (not
bracketed) verdict on the reactive port — resolving the ADR-0121 ambiguity and
deciding whether brick 3 (the multi-cell port) is needed or F2.3 is nearly
unblocked.

**Gate:** `reactive_deembed_matched_001` GREEN (resistor anchor asserted; reactive
arms asserted to the pinned result); the existing lumped + CPML gates non-regressed.

**Not in scope:** the multi-cell aperture port (brick 3, only if this confirms the
limit); F2.3 re-run (follow-on, once the verdict is PORT-CORRECT or the port is
fixed); the studio UI (Track B).

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-8-matched-line-deembed-bench-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-8-matched-line-deembed-bench.md`.
- ADR-0122 per-axis CPML; ADR-0119/0121 the PEC-source bench + its limits.
- Standard FDTD matched-line VNA de-embedding (incident/reflected two-run).
