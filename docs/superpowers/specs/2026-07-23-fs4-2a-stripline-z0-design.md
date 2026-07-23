# FS.4.2a — stripline Z₀ vs closed form (H-probes + V/I extraction gate)

**Date:** 2026-07-23 · **Track:** FS.4 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-compute/**`, `crates/yee-engine/**` (+ docs)
**Predecessor:** FS.4.0 (ADR-0215) — `voxelize_stackup` + `engine-stripline-eeff-001` (ε_eff 0.065 % vs exact TEM); FS.4.1 (ADR-0221) vias.

## Goal

The FS.4 roadmap row's next named gate: **stripline Z₀ vs closed form**. For a
zero-thickness symmetric stripline the exact conformal-mapping result is

  Z₀ = (η₀ / 4√ε_r) · K(k′)/K(k),  k = sech(πw/2b),  k′ = tanh(πw/2b)

(Pozar §3.8 quotes the Wheeler-fit approximation; use the exact K-ratio form and
cite both — the fit is within ~1 % and serves as a cross-check on the reference
itself). This is the first *impedance* (not propagation-constant) full-wave
validation in the engine flow — it requires measuring V and I on the line.

## Physics / extraction method

- **V(t)** at a measurement plane: line integral of E from ground plane to trace
  along z — computable TODAY as a column of point `Ez` probes summed ·Δz host-side.
  (Stripline: integrate from the k=0 ground to the trace plane; symmetry makes the
  lower half sufficient, but integrate the actual path used and document it.)
- **I(t)** at the same plane: Ampère loop ∮H·dl on a rectangular contour around the
  trace in the (y,z) cross-section plane — **requires H-component probes, which the
  engine does not have** (`Probe.component` is `EComponent`-only). That is
  deliverable 1.
- **Z₀ = V/I of the forward traveling wave**: drive the line (existing stripline
  fixture idiom from `engine-stripline-eeff-001`), time-gate the first pass of the
  pulse at the measurement plane (before any end reflection returns — same
  windowing discipline as the ε_eff gates), and take Z₀ from the gated V/I ratio
  (time-domain ratio at pulse peak, or band-averaged |V(f)/I(f)| over the fixture's
  validated band — implementer picks per what the fixture supports, documents why,
  reports both if cheap).
- Yee-grid staggering: V(Ez) and I(H) live at different half-step times/positions;
  handle the ½-cell/½-step offsets explicitly (document the choice; at TEM
  frequencies the correction is small but state it rather than ignore it).

## Deliverables

1. **H-component probes** in `yee-compute`: `Probe` gains H sampling (e.g. a
   `FieldComponent` enum or parallel `HComponent` field — pick what fits the
   existing plumbing with the least churn; CPU exact, GPU via the existing
   record_probes path with a CPU↔GPU parity gate on H probes, or a named
   `Unsupported` rejection if the GPU path is disproportionate — walking-skeleton
   rules apply, rejection must be tested).
2. **Gate `engine-stripline-z0-001`** (yee-engine tests): symmetric stripline via
   `voxelize_stackup` (b ≥ 16 cells per the ADR-0215/0221 lesson), V-column + H-loop
   probes at a source-far plane, gated V/I → Z₀, vs the exact K-ratio closed form.
   Target ≤ 5 % error; pin the measured value with honest margin (repo convention).
   If the first measurement misses by > 10 %, STOP and root-cause (staggering,
   loop placement, window) before pinning anything — do not widen to pass.
3. **ADR-0225** + FS.4 roadmap row update.

## Constraints

- Existing gates unmodified/green (bit-exact suite; `engine-stripline-eeff-001` in
  particular — same fixture territory).
- New probe machinery: no change to existing probe semantics (E probes bit-identical
  results before/after — the recorded stream layout may not silently reorder).
- Closed form implemented with a documented complete-elliptic-K (AGM iteration is
  ~10 lines, no new dependency).

## Non-goals

Per-layer tan δ, MoM cross-check, automesh stackup integration (rest of FS.4.2);
microstrip Z₀ (harder reference — quasi-TEM, not exact); port-impedance
renormalization of S-parameters to the measured Z₀ (future).
