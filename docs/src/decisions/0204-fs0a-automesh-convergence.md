# ADR-0204: FS.0a — push-button meshing: the rulebook + the convergence loop

**Status:** Accepted
**Date:** 2026-07-08
**Related:** FULL-SUITE-ROADMAP FS.0 (market research: manual meshing is the
#1 practitioner-cited open-EM adoption barrier; AMR is the core of HFSS's
commercial accuracy positioning), ADR-0199 (the shared board fixture this
rides), S.12/ADR-0189 (the directional |S21| observable).
**Spec:** `docs/superpowers/specs/2026-07-07-fs0-automesh-design.md`

## Decision

`yee_engine::automesh` ships the two pieces a novice cannot be expected to
know, as code:

1. **The rulebook** — `auto_dx(layout, f_max)` picks the largest dx
   satisfying dx ≤ λ_min/20 (in-dielectric), dx ≤ h_substrate/3, and
   dx ≤ min_feature/2 (`min_feature_m`: AABB widths + axis-aligned gaps),
   clamped to [1 µm, 1 mm]. Unit-gated per rule.
2. **The convergence loop** — `converge_two_port(layout, reference, opts,
   freqs, tol, max_passes)`: solve the two-port at dx₀, refine dx → dx/√2,
   re-solve, stop when the max per-frequency Δ|S21| between consecutive
   passes is ≤ tol. Unconverged results are reported honestly
   (`Converged::converged = false`), never hidden. Each pass holds the
   **physics constant and varies only the discretization** — see lesson 2.

## Two measured lessons (instrumented gate runs, 2026-07-08)

1. **The convergence criterion is LINEAR |ΔS|, not dB.** The first gate run
   converged the notch physics (4.900 GHz across two consecutive passes) yet
   measured **Δ = 15.35 dB** at the notch bin: near a deep null, a tiny
   frequency/depth shift produces tens of dB of per-bin delta while the
   linear |S| change is milliunits. Commercial adaptive refinement (HFSS's
   ΔS) uses the linear metric for exactly this reason. Tolerance for uniform
   staircased FDTD at walking-skeleton fidelity: **0.10** (HFSS's reference
   point is ~0.02; FS.0b's graded grid is the path toward it).
2. **Every pass must solve the same physical problem.** The board fixture
   sizes its CPML margin, air-under-lid height, and CPML absorber depth in
   *cells*; the first loop version left those counts fixed, so refining
   dx 0.533 → 0.267 mm silently halved the physical margin/lid/absorber
   (18.1 → 9.1 mm; CPML 5.3 → 2.7 mm). The per-bin dump showed the fine
   pass reading a non-physical broadband |S21| up to **+10.7 dB**: the stub
   junction scatters into the lowered lid's parallel-plate environment and
   the thinned absorber, the clean reference line doesn't, and the
   DUT/reference ratio explodes. Fix: the loop rescales `margin_cells`,
   `air_above_cells`, and the new `TwoPortBoardOptions::npml` each pass to
   hold their **metre** sizes constant; only dx/dt/n_steps vary.

## Gate `engine-automesh-001` (release, up to 6 solves)

The S.6 λ/4 open-stub notch board with **no hand-set dx anywhere**:
`auto_dx` seeds (measured **0.533 mm**, the h/3 substrate rule binding —
λ/20 = 1.19 mm, feature/2 = 1.50 mm), the loop refines
0.533 → 0.377 → 0.267 mm, and the converged curve is held against
transmission-line theory (f = c/(4·(l_stub+ΔL)·√ε_eff)):

- notch frequency within **5 %** of theory (the diagnostic runs before the
  constant-physics fix already read 4.900 GHz / 2.0 %)
- depth ≤ **−20 dB**
- `converged = true` within the 3-pass budget (final Δ|S| linear ≤ 0.10)

## Consequences

The single biggest usability gap starts closing with zero kernel risk: every
flow that consumes the shared board fixture (gates, studio, Python, WS) can
now be seeded and converged push-button. The loop is the exact workload the
GPU backend exists for (re-solves at shrinking dx; `opts.backend` passes
through). FS.0b (graded/nonuniform kernel, Taflove ch. 11) reuses this
solve→compare→refine scaffolding and tightens the tolerance toward the
commercial 0.02.
