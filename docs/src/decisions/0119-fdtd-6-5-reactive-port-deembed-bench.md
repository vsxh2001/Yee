# ADR-0119: Phase 2.fdtd.6.5 — reactive lumped-port V+I de-embedding bench (research-track increment 1)

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0116/0117/0118 (the two-way port + the two disproven reactive
attempts), ADR-0115 (F2.3, blocked), the lumped-LC → PCB goal (maintainer chose
"open the reactive-port research track"), [[project-lumped-lc-and-studio-redesign]]

---

## Context

The maintainer green-lit the multi-week **reactive lumped-port research track** to
unblock the goal's "EM simulation". Before committing weeks to a new port
formulation, increment 1 must resolve a **contradiction in the prior findings**:

- ADR-0117's investigation found the port-local reflection proxy `g·I/E*`
  (current and field measured *at the element*) gives the **correct** inductor
  magnitudes (≈0.51/0.37/0.25 across 4/6/9 GHz) — i.e. the port *does* present
  ~`jωL`;
- yet the **line-reflection measurement** (two-run difference + gated DFT +
  scalar-`A` / `z0_eff` calibration) reports the inductor as transparent (≈0.013)
  and the capacitor as a near-short.

These cannot both be the whole truth. Either (a) the port is correct and the
gate's *single-load reflection de-embed* mis-measures reactive (phase-shifting,
energy-storing) loads — in which case the fix is a better measurement and a
correct F2.3 element placement, **not** a multi-week port rewrite; or (b) the port
genuinely mis-loads the line and a new formulation is required. Increment 1
**decides which**.

## Decision

Build a **clean 1-port VNA-style de-embedding bench** (a new `#[ignore]`'d gate
`reactive_deembed_001`, yee-fdtd) that extracts the load impedance `Z_L(ω)`
*directly* from the measured **voltage and current** at a reference plane on the
parallel-plate TEM line (not from a calibrated reflection magnitude):

- measure the line's own `Z₀(ω)` from the incident-wave V/I ratio (a property of
  the discretised line, no fitting);
- with a single canonical lumped load (pure-L, pure-C, series-RLC) terminating /
  shunting the line, measure `V(ω)` and `I(ω)` at the load reference plane via
  single-bin DFT, form `Γ = (Z_in − Z₀)/(Z_in + Z₀)` and back out `Z_L(ω)`;
- compare `Z_L(ω)` to the intended `R + jωL + 1/(jωC)`.

**Outcome gate (decision, not pass/fail of the port):**
- if `Z_L(ω)` matches `R + jωL + 1/(jωC)` within a loose tol → **the port is
  correct**; the EM-sim blocker is the *measurement + F2.3 placement*, which
  increment 2 fixes (no port rewrite). Record this and re-scope the track.
- if `Z_L(ω)` does **not** match → the port truly mis-loads; increment 2 is the
  multi-cell-aperture / TL-Z₀ port reformulation. Record the measured `Z_L(ω)`.

Either way the bench itself is the deliverable: a trustworthy V+I reactive-load
measurement that the rest of the track is validated against. The bench asserts the
**known-good resistor** case (`Z_L → R`) so it cannot silently lie; the reactive
arms are asserted only to the extent the decision above warrants (never weakened
to a no-op; never a fake pass).

## Consequences

**Ships:** a validated V+I de-embedding bench + a recorded, evidence-based verdict
on whether the reactive port is correct — the de-risking foundation for the rest
of the research track, and possibly a shortcut (if the port is fine, EM-sim is
much closer than "multi-week").

**Gate:** `reactive_deembed_001` GREEN in CI (resistor `Z_L → R` asserted);
container-iterated; the reactive verdict recorded in this ADR's follow-up.

**Not in scope:** the port reformulation itself (increment 2, gated on this
verdict); F2.3's element placement (increment 2+); the studio UI (Track B).

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-5-reactive-port-deembed-bench-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-5-reactive-port-deembed-bench.md`.
- ADR-0117/0118 outcomes (the contradiction this resolves).
- Taflove & Hagness, lumped-element FDTD; standard VNA 1-port de-embedding.
