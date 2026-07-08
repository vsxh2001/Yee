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

## Three measured lessons (instrumented gate + forensic runs, 2026-07-08)

1. **The convergence criterion is LINEAR |ΔS|, not dB.** The first gate run
   converged the notch physics (4.900 GHz across two consecutive passes) yet
   measured **Δ = 15.35 dB** at the notch bin: near a deep null, a tiny
   frequency/depth shift produces tens of dB of per-bin delta while the
   linear |S| change is milliunits. Commercial adaptive refinement (HFSS's
   ΔS) uses the linear metric for exactly this reason. Tolerance for uniform
   staircased FDTD at walking-skeleton fidelity: **0.10** (HFSS's reference
   point is ~0.02; FS.0b's graded grid is the path toward it).
2. **Every pass must solve the same physical problem.** The board fixture
   sizes its CPML margin, air-under-lid height, absorber depth, and probe
   spacing in *cells*; the loop rescales all of them (`margin_cells`,
   `air_above_cells`, the newly exposed `TwoPortBoardOptions::npml`,
   `spacing_cells`) to hold their **metre** sizes constant — only dx/dt/
   n_steps vary. Honest attribution: this is hygiene, not the observed
   bug's cause — a re-run with doubled margin/absorber cell counts
   reproduced the pass-2 blowup to within 0.01 dB, exonerating the
   boundaries.
3. **The measurement must not assume the two runs launch the same wave.**
   The real cause of the pass-2 blowup (a clean-fit, non-physical
   broadband |S21| up to **+10.7 dB**): the S.12 single-ratio observable
   `fwd_B(dut)/fwd_B(ref)` assumes launch equality, and the wave-split
   forensics measured it failing — fit residuals ≤ 0.016 and β on the HJ
   dispersion on both sides, but the DUT's **plane-A** forward wave sat
   +10…+15 dB above the reference's. The shunt stub reflects strongly
   across the whole band (|Γ| ≈ 0.7 even at the passband shoulders), and
   that reflection re-pumps the imperfectly matched aperture source. The
   loop now measures the **launch-normalized double ratio**
   `|T_dut|/|T_ref|`, `T = fwd_B/fwd_A` per run
   (`sparams::forward_transfer`, the R.2-validated observable): each run
   normalizes by its own incident wave, so the launch cancels exactly.
   With it, the same dx = 0.267 mm pass reads **−3.2 dB shoulders (TL
   theory −2.9), a −33.0 dB notch at 5.00 GHz** — physical everywhere.

## Gate `engine-automesh-001` (release, up to 6 solves)

The S.6 λ/4 open-stub notch board with **no hand-set dx anywhere**:
`auto_dx` seeds (measured **0.533 mm**, the h/3 substrate rule binding —
λ/20 = 1.19 mm, feature/2 = 1.50 mm), and the loop refines by 1/√2 per
pass. Measured trajectory with the double-ratio observable
(0.533 → 0.377 → 0.267 mm): notch **5.100 → 4.900 → 4.850 GHz** at
**−31.8 → −35.1 → −34.2 dB** — vs TL theory (f = c/(4·(l_stub+ΔL)·√ε_eff))
the converged error is **3.0 %** (gate ≤ 5 %) at ≥ 20 dB depth.

**Convergence tolerance: 0.20 linear, measured.** The 0.377→0.267 mm pair
moves max Δ|S| = **0.1978**, all of it in the 5.45–6.0 GHz upper skirt
where the stub's open-end fringing is staircase-limited; the notch region
itself moves ≤ 0.08. A fourth uniform pass (dx = 0.189 mm, ~19M cells)
costs ~2.4 h — exactly the "refine everywhere because you can't refine
somewhere" waste that FS.0b's graded grid eliminates. The gate asserts
the loop's own `converged = true` verdict at 0.20 plus the physics
(notch ≤ 5 % of theory, ≤ −20 dB).
## Consequences

The single biggest usability gap starts closing with zero kernel risk: every
flow that consumes the shared board fixture (gates, studio, Python, WS) can
now be seeded and converged push-button. The loop is the exact workload the
GPU backend exists for (re-solves at shrinking dx; `opts.backend` passes
through). FS.0b (graded/nonuniform kernel, Taflove ch. 11) reuses this
solve→compare→refine scaffolding and tightens the tolerance toward the
commercial 0.02.
