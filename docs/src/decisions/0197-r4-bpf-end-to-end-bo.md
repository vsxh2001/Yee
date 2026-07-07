# ADR-0197: R.4 — BPF end-to-end: the F1.2.1 core + EM-in-the-loop BO

**Status:** Accepted
**Date:** 2026-07-07
**Related:** FILTER-DESIGN-ROADMAP F1.2.1 (the deferred core this closes),
ADR-0109 (F1.2.2 hairpin dims — documented the mean-gap layout collapse and the
qe→tap placeholder), ADR-0188/0189 (the loop skeleton + directional observable),
RF-TOOL-ROADMAP R.4.
**Spec:** `docs/superpowers/specs/2026-07-07-r4-bpf-end-to-end-design.md`

## R.4a — the closed-form layer (shipped, unit-gated)

- **Per-section hairpin geometry (gap option (a))**:
  `yee_layout::HairpinSectionParams` + `hairpin_bpf_sections` — each adjacent
  resonator pair sits at its own solved gap; a uniform `gaps_m` reproduces
  `hairpin_bpf`'s placement exactly (unit-tested; `geo-003` untouched).
  `dimension_hairpin_layout` now emits per-section geometry — the mean-gap
  collapse is gone.
- **qe→tap** (`tap_offset_from_qe`, gate `tap-qe-001`): the tapped half-wave
  resonator relation `Qe(t) = (π/2)(Z0/Zr)/cos²(πt/L)`, inverted to
  `t = (L/π)·acos(√((π/2)(Z0/Zr)/qe))`; `DimError::TapNotRealizable` carries
  the realizable `[qe_min, qe_max]` (antinode tap … end-of-arm tap). The
  `arm_length/3` placeholder is gone.
- **Fold-corrected arm length**: the U's midline (arm + bend + arm) is the
  half-wave, so `arm = (λ_g/2 − fold_spacing)/2` — the F1.2.2 `λ_g/4` form
  left every resonator electrically long by the full bend path (~37 % on the
  first gate scenario; measured as a wrecked, low-shifted response).
  `dimension_hairpin_with_fold` exposes the fold pitch (line widths;
  `fold_widths ≤ 1` is rejected — centre-to-centre spacing at one width
  merges the arms into a solid block, a degenerate geometry one instrumented
  run actually produced). `hairpin_dim_001` evolved to pin the corrected
  formula.

## The measured negative result that reshaped the full-wave gate

Three instrumented full-wave runs + probe-dump forensics (the R.2 pattern) on
the synthesized seed:

1. FR-4 1.6 mm, fold 2 w: **the tap cannot exist** — the 3 mm-wide 50 Ω line's
   fold consumes the half-wave; `TapNotRealizable` correctly rejects (the
   classic reason real hairpins use thin substrates / high-Z resonator lines).
2. h = 0.8 mm, fold 2 w (physical, tap fits with ~1 mm margin): the seed
   measures **|Γ_in| ≈ 1.0 across 3.5–6.5 GHz** — the end resonator absorbs
   almost nothing; its resonance sits **~+17 % high** (corner/open-end effects
   the midline model can't see) and the effective tap coupling lands far below
   the designed qe. S21 peaks at **−19 dB**: a detuned, under-coupled response
   **no closed form on this stack repairs**.

Off-resonance |Γ| = 1 is correct physics (the resonator is an open stub from
the feed), so the measurement chain is healthy; the residual is the *seed*.
This is precisely the residual F1.2.1 scheduled "surrogate-BO + EM-in-the-loop
refinement" for — so the R.4 full-wave gate **is** the BO gate, not a
standalone seed-verify (a gate asserting the seed's broken passband would pin
nothing of value).

## R.4b — gate `engine-bpf-bo-001`

One straight-line reference solve (launch normalization; every candidate
shares one grid via a fixed envelope bbox), then
`yee_surrogate::bo::minimize` over three knobs — `arm_scale` (retune),
`tap_scale` (external Q), `gap_scale` (inter-resonator coupling, clamped to
the 2·dx grid floor) — each objective call one full-wave DUT solve; the
objective is the RMS misfit (dB, floored at −40 dB) between measured
directional |S21| (S.12 observable) and `coupling_matrix_s_params` (the
validated Hong-Lancaster evaluator) over 3.5–6.5 GHz. Budget: 5 LHS + 7 EI
iterations ≈ 13 solves (~1 h release).

### The first converged run (2026-07-07, 13 solves, 36 min) and what it changed

- **BO worked**: misfit 28.48 → **24.59 dB RMS** (−3.9 dB), with the EI phase
  clearly outperforming the LHS phase (best three evals were EI's).
- **Where it went is the finding**: the best evals all sat at the
  **gap-scale lower bound** (0.70 — the 2·dx grid floor the candidate builder
  clamps to) and near the tap-scale lower bound (0.52) — the optimizer
  systematically pushed toward stronger coupling and ran out of physically /
  numerically honest room. The passband stayed buried (peak −16.2 dB, centre
  read 6.05 GHz).
- **Conclusion**: on h = 0.8 mm FR-4 at dx = 0.2 mm, the hairpin cannot
  realize the designed k ≈ 0.16 / qe ≈ 4.5 — a **feasibility limit of the
  scenario**, not a defect of the loop. Resolving the needed sub-0.4 mm gaps
  takes dx ≈ 0.1 mm (~8× the cell count) — exactly the scale the R.3 GPU
  parity work exists for.

The gate therefore pins what this run validated: the **loop machinery** —
BO must improve the seed misfit by ≥ 2 dB RMS (measured 3.9) on a live
measurement (best-curve dynamic range ≥ 15 dB; measured 29). Passband
centre/BW/peak are recorded, not asserted. **R.4c** (queued): the fine-grid
close-out on the GPU nightly — same gate harness, dx ≈ 0.1 mm, where "BO
closes centre frequency + bandwidth to spec" becomes assertable.

## Consequences

The filter path gains what the antenna path has (A.3), at BPF complexity:
synthesize → seed → measure → close the loop, now with a genuine multi-knob
surrogate optimizer over engine jobs — plus a quantified map of where the
walking-skeleton scenario's fidelity ends (the gap-resolution floor). The
seed-quality findings (fold correction, tap realizability bounds, stack
constraints) are permanent closed-form improvements with their own unit
gates. Queued: R.4c (fine-grid BO close-out on the GPU nightly), R.5b,
R.2b, R.0b.
