# ADR-0189: S.12 directional S-parameters — the refine loop converges

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0188 (the blocker this removes), ADR-0129/0131/0133 (the F2.3
standing-wave artifact and `fit_standing_wave`, ported here), ADR-0183 (sparams).
**Spec:** covered by `2026-07-06-s11-em-in-loop-refine-design.md` + ADR-0188's fix plan.

## Decision

`yee_engine::sparams` gains the repo's established directional-measurement machinery,
as pure protocol-side post-processing (no engine change):

- `fit_standing_wave(v0, v1, v2, spacing)` — verbatim port of `yee-voxel`'s three-probe
  fit: `cos(βd) = (V₀+V₂)/(2V₁)` recovers β, a linear solve splits
  `V(x) = a·e^{−jβx} + b·e^{+jβx}` into forward/backward phasors (degenerate βd → 0/π
  guarded, consistency residual exposed).
- `directional_transmission_db(dut_triple, ref_triple, dt, spacing, freqs)` — the
  **forward-wave** |S21|: immune to the port-to-port reflected wave that rippled the
  single-probe ratio.

Unit gates: the fit recovers known (a, b, β) to 1e-9; a DUT carrying 0.5× the forward
wave **plus a strong backward wave** still reads −6.02 dB.

## Measured result — the ADR-0188 blocker is gone

The refine loop (`engine-refine-001`), identical except for the observable (three
probes, 5.1 mm spacing, on the output feed):

| iter | synth f_c | measured cutoff | err |
|---|---|---|---|
| 0 (seed) | 2.000 GHz | 1.460 GHz | −27.0 % |
| 1 | 2.740 GHz | 2.510 GHz | +25.5 % |
| 2 | 2.380 GHz | 2.210 GHz | +10.5 % |
| 3 | 2.129 GHz | **2.020 GHz** | **+1.0 %** |

The map is smooth and monotone where the single-probe version oscillated
(−27 → +26.5 → −13 → +16.5), and the secant converges: **a filter synthesized by the
pipeline, measured by the engine, and corrected until its measured cutoff sits 1.0 %
from the design target**. The convergence asserts (final ≤ half seed, ≤ 10 %) are
restored and the gate joins CI.

## Consequences

The engine now *designs*: synthesize → full-wave verify → correct → converge, all over
the job protocol. Follow-ons: the multi-knob/BPF slice of F1.2.1 (per-section knobs,
`yee-surrogate` BO when the knob count grows), directional S11 at the input feed by the
same fit, and the antenna track (which inherits every fidelity lever built in
S.9–S.12).
