# FS.4.2b — per-layer stackup loss: σ map + stripline attenuation gate

**Date:** 2026-07-23 · **Track:** FS.4 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-voxel/**`, `crates/yee-engine/**` (+ docs)
**Predecessors:** FS.4.0 stackup (ADR-0215), FS.4.2a Z₀ gate (ADR-0225), FS.2c loss/efficiency (single-layer `substrate_sigma_cells`).

## Gap

`StackupLayer.loss_tangent` exists but nothing consumes it: `substrate_sigma_cells`
(FS.2c) maps ONE tan δ over the whole substrate ε-map. Multilayer boards mix
materials (FR-4 core tan δ ≈ 0.02 over low-loss prepreg ≈ 0.004); per-layer loss is
table stakes for the FS.4 "every real board" story.

## Deliverables

1. **`yee_voxel::stackup_sigma_cells(model, stackup, f_ref_hz) -> Vec<f64>`** —
   per-cell σ = 2π f_ref ε₀ ε_r(layer) tan δ(layer), assigned by each cell's layer
   (same k-band bookkeeping `voxelize_stackup` used to fill ε; air/metal cells σ = 0).
   Unit gate: two-layer stackup with distinct tan δ → exact per-band σ values,
   boundaries at the right k; single-layer case must agree with
   `substrate_sigma_cells` at the same tan δ (consistency pin, not bit-equality if
   the code paths differ — assert exact equality only if both compute the identical
   expression).
2. **Gate `engine-stripline-alpha-001`** (yee-engine): symmetric stripline
   (FS.4.2a fixture idiom, ε_r 2.2, b ≥ 16 cells), tan δ = 0.02 via
   `stackup_sigma_cells`. Stripline is pure-TEM entirely inside the dielectric, so
   the closed form is exact: α_d = (π f √ε_r / c)·tan δ nepers/m (× 8.686 → dB/m);
   conductor loss = 0 (PEC). Measure attenuation from the gated forward-wave ratio
   at two measurement planes a known distance apart (two V-columns per FS.4.2a's
   extraction, or the launch-normalized double-ratio idiom — ADR-0204 lesson:
   per-plane normalization, no absolute single ratios). Evaluate at a few bins
   around f₀; note σ is frequency-independent (constant-σ model) while true tan δ
   loss scales ∝ f — grade at f_ref where the model is exact, and document the
   off-f_ref deviation rather than asserting on it. Target ≤ 10 % on α at f_ref
   (loss extraction is noisier than phase); pin measured + margin; > 20 % off →
   STOP and root-cause, never widen.
3. **ADR-0226** + FS.4 roadmap row update.

## Constraints

- Existing gates unmodified/green: `engine-stripline-eeff-001`, `engine-stripline-z0-001`
  (lossless fixtures — must not be touched by the loss plumbing), bit-exact suite,
  `voxel-stackup-001/002`, the FS.2c σ tests.
- Loss OFF (all tan δ = 0) must be a provable no-op: σ all-zero → existing lossless
  results bit-identical (cheap assert in the unit gate).

## Non-goals

Frequency-dependent σ(f) / Debye fitting (dispersive ADE lane exists for that);
conductor (finite-σ metal) loss; MoM cross-check; automesh stackup integration.
