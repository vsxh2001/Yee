# ADR-0207: FS.2a — aperture-port power records: the gain denominator

**Status:** Accepted
**Date:** 2026-07-08
**Related:** ADR-0206 (FS.1b — flagged that every pattern gate is
relative-|E| only), S.10/ADR-0187 (the aperture port whose state this
records), R.3/ADR-0196 (the GPU `Unsupported` idiom reused).
**Spec:** `docs/superpowers/specs/2026-07-08-fs2-farfield-products-design.md`

## Decision

`AperturePortSpec::record` (serde default `false`) → per-step
`(v_src, v_terminal, i_branch)` triples in `JobResult::port_records`. The
port loop already computes all three every step (the LumpedRlcPort
correction), so recording is free of new field probes and modal-
normalization guesswork. CPU-only; the GPU backend rejects recording
ports with `Unsupported` (Auto falls back) until a readback buffer lands.

## Measured lesson: account on the circuit side, not the aperture side

The first gate run summed the naive aperture-side `v_term·i` and read a
non-physical B/A energy ratio of **1.596**: the port's implicit β term
(`β = dt·h/2ε₀A` ≈ **14.5 Ω** on this fixture — comparable to the 50 Ω
branch!) sits between the terminal voltage and the resistor, so `v_term·i`
mixes real transfer with reversible discretization storage. The honest
identity uses **circuit-side quantities only**: EMF supply `Σ −v_src·i·dt`
vs branch-resistor dissipation `Σ i²R·dt`.

## Gate `engine-power-001` (release, 1 solve ≈ 55 s) — GREEN

Lossless matched through line, both ports recording:

- **closure = 0.9917** (pinned [0.95, 1.0]): E_emf = 5.557e-13 J supplied;
  A's resistor 2.708e-13 + B's resistor 2.802e-13 dissipated; 0.8 % left
  in CPML leakage + residual ring.
- accepted-by-field `E_emf − E_R(A)` = **51.3 %** of the supply — the
  textbook matched-source halving, and the FS.2b gain denominator.
- passive-branch sample-wise sanity `v_term·i ≥ 0` (a load only
  dissipates).

## FS.2b (same ADR): the absolute-scale forensic — three measured results

1. **`farfield::gain_dbi` shipped** (normalization chain audited and
   documented in the module: NTFF = continuous-transform pattern
   amplitude; accepted density = per-frequency circuit-side identity;
   every 1/π and DFT scale cancels in G = 4π|F|²/(η₀·p_acc)).
2. **Gate `engine-scale-001` GREEN — the NTFF absolute scale is right.**
   Soft `E_z += s` is an exactly-known Hertzian moment
   (I·dl = ε₀·dx³·S(ω)/dt); measured NTFF/analytic = **1.048 (θ = 90°) /
   1.029 (45°)**, reproducible to 3 decimals across dx ∈ {1.5, 2, 2.5} mm
   × f ∈ {1.8, 2.45, 3.2} GHz. Lesson: the first attempt used a BASEBAND
   Gaussian and measured ±40 % direction-dependent scatter — near-DC
   energy survives the CPML and leaks into the single-bin DFT; the
   zero-DC `SourceSpec::GaussianPulseEz` was added for it.
3. **Gate `engine-gain-001` RED, root-cause hypothesis.** Patch read
   22.15 dBi (physics caps it at ~3–6), array 23.92, while the
   differential (1.77 dB) and all relative patterns stay healthy. With
   the transform and the port power both independently certified, the
   excess isolates to the fixture: the voxelizer's substrate slab spans
   the whole domain, so the equivalence box intersects dielectric where
   the strongest guided fields live and free-space η₀ misprices them.

## Queued

**FS.2b.1**: finite-extent substrate in the voxelizer (real boards end;
the NTFF box then passes through air — the openEMS practice), then
re-measure engine-gain-001. FS.2c: efficiency + full-sphere export. GPU
port-record readback.
