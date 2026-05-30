# ADR-0125: Phase 2.fdtd.6.9 — multi-cell aperture lumped port

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0124 (F2.3 sheet placement — necessary but insufficient; the
grid-vs-port investigation lives in its Outcome + here), ADR-0121 (the ≈0.37
single-cell floor), ADR-0119/0123 (the de-embed benches), ADR-0116 (the two-way
port), ADR-0122 (per-axis CPML), the lumped-LC → PCB goal (maintainer-approved
reactive-port research track), [[project-lumped-lc-and-studio-redesign]]

---

## Context — the grid-vs-port investigation (decisive)

A design-investigation (de-risking before implementing) ran F2.3's full-width-sheet
sim at dx = 0.4 / 0.2 / 0.1 mm (n_steps ∝ 1/dx) and a zero-FDTD placement
diagnostic. Result: **PORT-LIMITED, not grid-limited — finer dx makes it strictly
worse.** As dx shrinks the FDTD |S21| converges to a **transparent through-line**
(→ 1.0 flat across 1.6–2.8 GHz); the coarse 0.4 mm run has *more* residual
structure (0.5 dB) than the fine runs (0.0 dB). No notch ever forms.

**Mechanism (with a scaling law):** the shunt resonators are extreme low-impedance
tanks (`|Z_L|=|Z_C|≈2.93 Ω` ≪ Z₀=50 Ω). With the ADR-0124 `C/N`, `N·L` sheet
(N = trace-width cells ∝ 1/dx):

- the **capacitor arm is dx-invariant** — `ε_eff = ε₀ + C_cell·dz/dA` with the
  cubic `dz/dA = 1/dx` cancels the `C/N` scaling, so `ε_eff/ε₀ ≈ 1000` at every dx
  (a frozen per-cell near-short);
- the **inductor arm collapses as O(dx²)** — its two-way field back-action is
  `dt²·dz/(ε₀·dA·L_cell)` with `dt ∝ dx`, `dz/dA = 1/dx`, `L_cell = N·L ∝ 1/dx` ⇒
  `∝ dx²` (exactly 4× weaker per 2× refinement, measured).

A parallel L‖C resonance needs both arms balanced. The inductor (half the tank)
goes inert (`O(dx²)`) while the capacitor stays a fixed short ⇒ the sheet
degenerates to a frequency-flat shunt capacitance, which the DUT/thru
normalization divides out to `|S21|≈1.0`. This is the same single-cell reactance
floor ADR-0121 measured (≈0.37; `ε_eff≈4.6` lumped-C-as-dielectric-scatterer) —
here the tanks sit far past it (`ε_eff/ε₀≈1000`). **No dx meets the 20 dB gate; the
floor is fundamental to the single-cell port formulation.**

## Decision

Implement a **multi-cell aperture lumped port** in `yee-fdtd` whose field coupling
is referenced to the **modal port-face**, not one Yee cell — the formulation the
investigation pins (in priority order):

1. **Port-face-area back-action.** Inject the lumped branch current into the field
   with a back-action referenced to the **modal aperture area `w·h`** (trace width
   × substrate height) the mode occupies — NOT the bare single-cell `dA = dx²`.
   This is the root fix: it removes the `O(dx²)` inductor collapse (the back-action
   becomes dx-stable). (= ADR-0121's flagged-but-unvalidated "modal coupling
   factor".)
2. **Modal branch voltage.** Tie the branch voltage to the modal `V = ∫E·dz` over
   the **full substrate height** (all `n_sub` `E_z` edges), not a single `E_z`
   edge.
3. **(y, z) aperture distribution.** Spread the lumped element over the full
   `(y, z)` port face (trace width × substrate height) with a value scaling that
   holds the **aggregate `Z_L` fixed independent of cell count** — so neither arm
   scales away with dx.
4. **De-embedded reference plane.** Read the port `S`/`Z` from line currents at a
   reference plane (TL-style), not a single sense edge.

Keep the existing single-edge `series_rlc`/`pure_resistor` API + the resistor-exact
path (additive: a new aperture-port constructor / type, or an aperture spec on the
port). Validate on a de-embed bench: the **capacitor reactance must improve from
≈0.37 toward the loose tol** AND the **`O(dx²)` inductor collapse must be gone**
(reactive `Z_L` dx-stable). Then F2.3's `fdtd_lumped_001` re-runs on top.

## Consequences

**Ships (target):** a physically-correct, dx-stable reactive lumped FDTD port — the
EM-sim unblocker. With it, F2.3's L‖C tanks resonate and the band-pass forms →
the goal's EM-simulation component. The investigation has **de-risked** it from
"uncertain" to a specified formulation with a known failure mechanism to fix.

**Gate:** a de-embed bench shows the aperture port's reactive `Z_L` accurate +
dx-stable (resistor anchor still exact, never weakened); then `fdtd_lumped_001`
GREEN within its loose tol. The existing lumped/CPML/de-embed gates non-regressed.

**Scope/risk:** a substantial yee-fdtd-core formulation + re-validation (likely
multiple sub-increments: the port + its bench validation, then F2.3). The
investigation's mechanism makes the FIRST sub-increment concrete: fix the
back-action normalization (item 1) + the (y,z) aperture (items 2–3) and show the
dx² collapse is gone on the bench.

**Not in scope:** tight-tol EM; SRF/ESR; the studio UI (Track B).

---

## References
- ADR-0124 Outcome (the dx sweep + the O(dx²) diagnostic); ADR-0121 (≈0.37 floor);
  ADR-0119/0123 (de-embed benches); ADR-0122 (per-axis CPML).
- `docs/superpowers/specs/2026-05-30-fdtd-6-9-aperture-lumped-port-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-9-aperture-lumped-port.md`.
- Taflove & Hagness lumped-element FDTD; CST/HFSS modal/aperture lumped-port
  formulations (integration line over the port face).
