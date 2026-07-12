# ADR-0215: FS.4.0 — multilayer stackups: Stackup type + voxelizer + the stripline gate

**Date:** 2026-07-12 · **Status:** accepted · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-12-fs4-stackup-design.md`

## Context

The engine protocol already carries per-cell ε and 3-D masks; the
single-substrate assumption lived only in the layout description and the
voxelizer. FS.4.0 is the walking skeleton for N-layer boards.

## Decision

1. **`yee_layout::Stackup { layers, lid }`** (bottom-up, ground below
   layer 0; `lid: true` = PEC sheet directly on the last layer) +
   `StackupLayer { eps_r, height_m, loss_tangent }` +
   `Stackup::symmetric_stripline`. Serde like the rest of the layout
   types.
2. **`yee_voxel::voxelize_stackup(layout, &stackup, trace_layer, opts)`**:
   uniform grid; layers quantize to `round(h/dx).max(1)` cells and fill
   **contiguously from k = 0** — no air gap anywhere in the stack (the
   ADR-0108 lesson generalized); trace PEC at the top of
   `layers[trace_layer]`; open tops get one guaranteed cell above the
   stack (mirroring the microstrip voxelizer, which makes the
   single-layer case **bit-identical** — gate `voxel-stackup-001`);
   lidded domains end exactly at the lid plane. Returns the existing
   `MicrostripModel`, so every downstream fixture works unchanged.

## Measured gates

- **`voxel-stackup-001`** (instant, GREEN): single-layer open stackup ≡
  `voxelize_microstrip` bit-identical (ε + both PEC masks); 3-layer
  ε-bands land exactly (2.2/4.4/3.0 at their k-ranges, no gaps), buried
  trace at the chosen interface, lid masked whole-plane.
- **`engine-stripline-eeff-001`** (release, ignored, GREEN): symmetric
  stripline (ε_r = 4.4, b = 3.2 mm, lid on) — time-gated two-probe phase
  velocity vs the **exact** TEM value ε_eff = ε_r. **Measured 0.065 %**
  (grid 1184×48×16, dx = 0.2 mm, ~89 s release); assert pinned ≤ 2 %.
  Runs in the blanket yee-engine CI release step.

## The three measured lessons (each cost one red run)

1. **Box modes:** a lidded homogeneous box is itself a waveguide; at the
   first margin (w = 22.3 mm) its TE₁₀ cutoff was 3.2 GHz — inside the
   drive band — and mode mixing read ε_eff 13.8 % high. Keep the box
   width's cutoff above the band (margin 10 cells ⇒ f_c ≈ 7.5 GHz).
2. **Window clipping:** narrowing the bandwidth lengthens the pulse; at
   L = 6 λ_g the reflection gate clipped the tail at probe B (14.5 %
   high). L = 8 λ_g restores headroom. Both hygiene bounds must hold
   simultaneously.
3. **Confined-mode resolution (the big one):** with b = 8 cells the
   stripline mode's transverse confinement (decay scale b/π ≈ 2.5 cells)
   is under-resolved; the discrete transverse operators break the TEM
   Laplacian cancellation and β reads 7 % high — probe-separation
   doubling scaled the phase excess 2×, proving a real β error, not a
   measurement artifact. b = 16 cells collapses the error to 0.065 %.
   **Lidded/confined modes need ≥ ~16 cells across the confinement
   dimension** — a meshing rule the FS.0 automesh rulebook should learn
   when stackups integrate (FS.4.2).

## Non-goals / queued (FS.4.1+)

Stripline Z₀ vs closed form (port-impedance extraction), simultaneous
multi-layer traces, through/blind vias, per-layer tan δ in the engine
materials, MoM cross-check, automesh integration (FS.4.2).
