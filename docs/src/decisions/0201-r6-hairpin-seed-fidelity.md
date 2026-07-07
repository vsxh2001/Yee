# ADR-0201: R.6 — hairpin seed fidelity: the corner correction lands on-frequency

**Status:** Accepted
**Date:** 2026-07-07
**Related:** ADR-0197 (R.4 — the +17 % measured detune and thick-stack tap
wall), ADR-0109 (F1.2.2 dims).
**Spec:** `docs/superpowers/specs/2026-07-07-r6-hairpin-seed-fidelity-design.md`

## Decision

`yee_filter::HairpinOptions { fold_widths, resonator_z_ohm, corner_widths }`
+ `dimension_hairpin_opts`; `dimension_hairpin`/`_with_fold` delegate.

- **Corner correction (default ON, κ = 0.85)**: each of the U's two 90°
  corners shortens the electrical path by ≈ κ·w relative to the midline, so
  `arm = (λ_g/2 − fold)/2 + κ·w`. κ was **calibrated from the single R.4
  instrumented data point** (the seed's |Γ_in| dip at 5.95 vs designed
  5.0 GHz → a 2.62 mm electrical deficit over two corners of a 1.529 mm
  line, net of open-end lengthening).
- **Resonator impedance** (`resonator_z_ohm`, default spec-Z0): resonator
  width / ε_eff / λ_g / fold / gaps computed at Zr; the tap solve's
  `(Z0/Zr)` factor is where thinner lines buy tap room. The feed stays a Z0
  line — `HairpinDimensions.feed_width_m` separates the two widths, and the
  layout builders / studio flows stop conflating them. Unit gate: Zr = 70 Ω
  dimensions the previously-`TapNotRealizable` h = 1.6 mm stack (thinner
  resonator line, Z0 feed, tap + feed half-width on the arm).

`hairpin_dim_001` evolved to pin the corrected formula; the studio
`design_e2e` unrealizable-spec case re-triggers via a wide fold
(`fold_widths = 3.5`) because the corner correction made the old trigger
realizable — itself evidence the correction adds real arm length.

## The blind full-wave check (re-run of `engine-bpf-bo-001`)

The κ calibration came from the Γ-dip of one run; the verification is the
independent S21 measurement of the corrected seed:

- **Seed passband centre: 5.000 GHz vs designed 5.00 GHz** (pre-R.6:
  5.95 GHz, +19 %). The corner correction lands on-frequency on its first
  blind test.
- Seed peak −33 dB: properly tuned, the response sits deeper in the
  under-coupled regime at f0 — the dx = 0.2 mm coupling floor (ADR-0197)
  is unchanged and remains R.4c's job.
- **The BO phase then went 2.5× further than pre-R.6**: misfit 29.87 →
  **19.92 dB RMS** (−9.95 dB; the R.4-era run managed −3.9), and the
  optimized response shows a **real emerging passband — peak −6.08 dB**
  (pre-R.6 best: −16.2), reading −7.4 dB at 5.5 GHz, centre +8 %, BW
  200 MHz vs 1100 designed. A better seed multiplies what the same solve
  budget buys — and the R.4c close-out criteria (peak ≥ −6 dB, centre
  ±5 %) are now within reach of the fine grid rather than fantasy.

## R.4c wiring (shipped with this ADR)

The BO gate is now **fidelity/backend-parameterized**
(`YEE_BPF_BO_DX_MM`, `YEE_BPF_BO_BACKEND`; steps scale inversely with dx;
the candidate builder's 2·dx gap floor relaxes automatically), and the GPU
nightly gains the R.4c job: dx = 0.1 mm on `backend: gpu`, where the gate's
close-out asserts activate (passband peak ≥ −6 dB, centre ±5 %, BW ±30 %).
It runs when `YEE_GPU_RUNNER_ENABLED` is set; its first hardware run may
retune numbers (the usual nightly-landing pattern).

## Consequences

The seed's frequency placement is solved; the remaining seed→spec residual
is coupling amplitude, which is a resolution problem (R.4c), not a formula
problem. κ = 0.85 is a single-stack calibration — revisit when a second
stack is measured (the constant is one line in `HairpinOptions`).
