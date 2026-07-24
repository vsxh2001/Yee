# FS.4.2c — automesh stackup integration (the b ≥ 16 lesson becomes a rule)

**Date:** 2026-07-24 · **Track:** FS.4/FS.0 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-engine/**` (+ docs)
**Predecessors:** FS.0a auto_dx rulebook (ADR-0204); FS.4.0 stackup + the measured
lesson "confined lidded modes need ≥ ~16 cells across b — at 8 cells β reads 7 % high
from discrete transverse-operator error — a future automesh rule, FS.4.2" (ADR-0215);
FS.4.2a/b stripline gates (ADR-0225/0226).

## Gap

`yee_engine::automesh::auto_dx(layout, f_max)` knows one substrate (`layout.substrate`).
Stackup boards (FS.4.0) have N layers + optional lid; nothing turns a `Stackup` into a
rulebook dx, and the ADR-0215 b-resolution lesson lives only in prose + hand-set
fixtures. Push-button meshing (the FS.0 wedge) must extend to multilayer.

## Deliverables

1. **`auto_dx_stackup(layout, stackup, f_max_hz) -> f64`** in
   `yee_engine::automesh` — largest dx satisfying:
   - wavelength: dx ≤ λ_min/20 with λ_min = c/(f_max·√ε_r_max) over the layers;
   - per-layer resolution: dx ≤ h_i/3 for EVERY layer (generalizes the h/3 rule;
     a buried interface under-resolved is the same silent failure);
   - feature: dx ≤ min_feature/2 (reuse `min_feature_m`);
   - **lid rule (the ADR-0215 lesson): if `stackup.lid`, dx ≤ b/16** where
     b = Σ h_i (total ground→lid dielectric height);
   - same [1 µm, 1 mm] clamp, same doc style as `auto_dx` (cite ADR-0215 for the
     16 in the doc comment).
   Unit tests: each rule binding in turn (construct stackups where each term is the
   min); single-layer no-lid case degenerates to `auto_dx` exactly (consistency).
2. **Gate `engine-automesh-stackup-001`**: the FS.4.0 stripline ε_eff fixture built
   with **no hand-set dx anywhere** — `auto_dx_stackup` seeds the grid (assert which
   rule binds: expected the lid rule b/16 for the standard fixture; print it) — and
   the measured ε_eff vs the exact TEM ε_r within the same ≤ 2 % bar as
   `engine-stripline-eeff-001` (pin measured + margin). The point: the rulebook
   alone lands inside the certified-fixture tolerance.
3. **ADR-0227** + FS.4 roadmap row (FS.4.2 remainder after this: MoM cross-check).

## Constraints

- Existing gates unmodified/green (both stripline gates, alpha gate, automesh gates,
  bit-exact suite). `auto_dx` itself untouched (new function, no behavior change).
- Honest pin; > 5 % ε_eff error from the rulebook grid → STOP and root-cause (is a
  rule too loose?) rather than widening — a loose rulebook is the product defect
  this gate exists to catch.

## Non-goals

Graded `auto_spacings` stackup variant (uniform-dx rulebook first — walking
skeleton; graded multilayer is a follow-on); MoM cross-check; automesh via/arc rules.
