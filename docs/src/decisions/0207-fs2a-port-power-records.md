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

## FS.2b.1 — the finite board: root cause CONFIRMED, gain gate GREEN

`yee_voxel::voxelize_finite_board` bounds the dielectric slab and the
ground sheet to bbox + margin (gate voxel_003: cell-exact bounds,
huge-margin ≡ infinite bit-identically). Re-measured on the finite board
(bbox + 15 mm, lifted stack, all-face CPML, box fully in air):

- single patch **5.42 dBi** (textbook 5–7); the pattern amplitude
  dropped **16.7 dB** while p_acc stayed identical (2.517e-23 →
  2.535e-23) — the 22 dBi excess was entirely the dielectric-crossing
  equivalence box, as hypothesized;
- 2×1 array **7.63 dBi**; differential **+2.21 dB** (ideal 3 minus
  mutual coupling). Asserts pinned: single [4, 7.5] dBi, differential
  [1.5, 3.5] dB.

FS.2b is complete: the pipeline now produces absolute commercial-grade
far-field numbers — gain in dBi certified by an analytic Hertzian pin, a
textbook patch window, and a bias-cancelling array differential.

## FS.2c — efficiency + full-sphere export: SHIPPED, GREEN first run

`farfield::sphere_grid` / `radiation_efficiency` (the gain theorem
`∮G dΩ = 4π·η` by midpoint quadrature; unit-gated: isotropic → 1 to
< 0.1 %) / `pattern_csv` (byte-stable export). Gate `engine-eff-001`
(A.1 patch, finite board, 12×16 sphere):

- lossless η = **0.806** (pinned [0.65, 1.0] — the certified NTFF scale,
  quadrature, and absorber leakage each shave a few percent);
- tan δ = 0.02 → η = **0.294**, squarely in the 30–60 % FR-4-patch
  literature range (substrate loss only); and the lossy antenna's
  accepted power ROSE (2.53 → 3.23e-23) — loss broadens the match, a
  free physics cross-check.

**FS.2 is complete**: the pipeline delivers the absolute far-field
products commercial tools quote — gain in dBi, radiation efficiency,
and a full-sphere export — each behind an analytic or textbook gate.

## Queued

GPU port-record readback. Consider migrating the A-track pattern gates
to the finite-board fixture (their relative asserts are unaffected but
the fixture is more physical).
