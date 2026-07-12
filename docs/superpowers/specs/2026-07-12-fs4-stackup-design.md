# FS.4.0 — multilayer stackup walking skeleton (design)

**Date:** 2026-07-12 · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`) · **ADR:** 0215

## Problem

Yee voxelizes exactly one stackup: ground + single dielectric + top-metal
trace (+ air). Real boards are N-layer; the protocol already carries
per-cell ε and 3-D masks, so the gap is purely the description + the
voxelizer. FS.4.0 is the walking skeleton: an N-layer `Stackup` type, a
stackup voxelizer with a buried trace, and the buried-line analog of
S.5's ε_eff gate.

## Design

- **`yee_layout::Stackup { layers: Vec<StackupLayer>, lid: bool }`**,
  bottom-up, ground plane below layer 0; `StackupLayer { eps_r, height_m,
  loss_tangent }`. `lid: true` puts a PEC sheet directly on top of the
  last layer (stripline/shielded boards); `lid: false` leaves the top
  open (air above, the microstrip case).
- **`yee_voxel::voxelize_stackup(layout, &stackup, trace_layer, opts)`**:
  uniform grid; layer i quantizes to `round(h_i/dx).max(1)` cells (height
  error ≤ dx/2 — acceptable at walking-skeleton tolerance, recorded in
  the ADR); ε fills each layer's E_z edges contiguously from k = 0
  (ground) — **no air gaps anywhere in the stack** (the ADR-0108 lesson
  generalized: any accidental series air gap between layers poisons
  every downstream result); trace PEC at the TOP of `layers[trace_layer]`;
  lid PEC (whole plane) at the stack top when `lid`. Returns the existing
  `MicrostripModel` (k_gnd = 0, ports at the trace plane) so every
  downstream fixture works unchanged.
- Single dielectric + `trace_layer` = last + `lid: false` reproduces
  `voxelize_microstrip` — pinned bit-identical (ε and masks) by a unit
  gate, the FS.0b idiom.

## Gates

1. **`voxel-stackup-001`** (instant): single-layer stackup ≡
   `voxelize_microstrip` bit-identical; a 3-layer stack puts each ε_r in
   the right k-band, trace at the right interface, lid masked, and **no
   air gap between layers**.
2. **`engine-stripline-eeff-001`** (release, ignored): symmetric
   stripline — 2 identical layers, trace at the mid interface, lid on —
   time-gated two-probe phase velocity (the S.5/verify_line_eeff method)
   → ε_eff vs the **exact TEM value ε_eff = ε_r** (homogeneous
   dielectric; no Hammerstad–Jensen approximation involved). Measured
   number pinned; walking-skeleton tolerance ≤ 5 %.

## Non-goals (FS.4.1+)

Stripline Z₀ vs closed form (needs port-impedance extraction), buried
traces on multiple layers simultaneously, through/blind vias, per-layer
loss tangents in the engine materials (single tan δ today), MoM
multilayer cross-check.

## Lane

`crates/yee-layout/**`, `crates/yee-voxel/**`, `crates/yee-engine/tests/**`,
`docs/**`. (Disjoint from the FS.0b.2-GPU worktree track in
`crates/yee-compute/**`.)
