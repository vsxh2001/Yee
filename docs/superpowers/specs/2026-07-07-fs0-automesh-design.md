# FS.0 — automatic meshing + convergence loop

**Date:** 2026-07-07
**Track:** FULL-SUITE-ROADMAP FS.0 (top priority per the market research:
manual meshing is the #1 practitioner-cited barrier to open-EM adoption).

## Decomposition

- **FS.0a (this spec's walking skeleton): auto-uniform dx + convergence.**
  No kernel change. Two pieces:
  1. `yee_engine::automesh::auto_dx(layout, f_max)` — the rulebook a novice
     can't be expected to know, as code: dx ≤ λ_min/20 (in the dielectric),
     dx ≤ h_substrate/3 (vertical field resolution), dx ≤ min_feature/2
     (smallest gap/width in the layout), floored to avoid absurd grids.
     Pure, unit-gated against hand-computed values.
  2. `yee_engine::automesh::converge_two_port(...)` — the HFSS-style
     adaptive-pass loop, FDTD flavour: solve the S.12 directional |S21| at
     dx₀ = auto_dx, refine dx by 1/√2 per pass, stop when the max |ΔS| dB
     between consecutive passes drops below tolerance (or the pass budget
     runs out — reported honestly in the result, never silently). Reuses
     `yee_engine::board::two_port_board_job`, so every design flow gets it
     for free. GPU-aware (each pass is one job; `backend` passes through) —
     the loop is exactly why the GPU matters: re-solves are the price of
     push-button accuracy.
- **FS.0b (follow-on): graded/nonuniform grid kernel** (Taflove ch. 11
  dual-step formulation) in yee-compute, gated bit-exact-on-uniform, then
  graded mesh *rules* (refine at edges/gaps, growth ratio ≤ 1.3) replacing
  the uniform ladder — the full commercial-style answer. Own spec when
  FS.0a's loop is proven.

## Gates

- Unit (`automesh` module tests): each rule binds on a scenario constructed
  to make it the binding constraint; the combined `auto_dx` picks the min.
- **`engine-automesh-001`** (release, up to 6 solves): the S.6 λ/4
  open-stub notch board with **no hand-set dx anywhere** — the loop starts
  from `auto_dx`, converges within the pass budget, and the converged notch
  frequency lands within 5 % of transmission-line theory
  (f = c/(4·(l+ΔL)·√ε_eff)) at ≤ −20 dB depth. Chosen over the plain-line
  ε_eff scenario because a resonant feature is the harder, more
  convergence-sensitive observable — it is exactly where a novice's
  hand-meshed run goes wrong. This is the "novice gets a trustworthy answer
  push-button" criterion made machine-checkable.

## Implementation lessons (2026-07-08, recorded in ADR-0204)

1. Convergence must be judged on **linear** |ΔS| (HFSS's ΔS convention),
   not dB — a converged deep notch still swings tens of dB per bin.
2. The loop must hold the fixture's physical sizes (CPML margin, air
   height, absorber depth) constant in metres as dx shrinks — the
   first version scaled them in cells and the fine pass read a
   non-physical +10.7 dB broadband |S21|.

## Consequences

The single biggest usability gap starts closing with zero risk to the
certified kernel. The convergence-loop infrastructure (solve → compare →
refine) is exactly what FS.0b's graded passes will reuse.
