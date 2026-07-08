# ADR-0206: FS.1b — the 2×1 corporate-fed patch array

**Status:** Accepted
**Date:** 2026-07-08
**Related:** ADR-0205 (FS.1a quasi-Yagi — the antenna-catalog track this
continues), ADR-0191/0193 (A.1/A.3 — the inset element and its measured
0.25·L depth), ADR-0201 (R.6 — the corner-path lesson the tree's symmetry
sidesteps), ADR-0204 (FS.0a auto_dx, which seeds every FS.1 gate).
**Spec:** `docs/superpowers/specs/2026-07-08-fs1b-patch-array-design.md`

## Decision

`yee_layout::patch_array_2x1(f0, substrate, z0)`: two A.0/A.3-certified
inset-fed patches side by side along y (an H-plane pair) at **0.5 λ₀**
centre spacing, fed in phase through a symmetric corporate tree — 50 Ω
spine (λg/2, probe room) → junction → two **λg/4 70.7 Ω transformers**
along ±y (each transforms its branch's 50 Ω to 100 Ω; the pair in
parallel presents 50 Ω at the junction) → 50 Ω branches → each element's
inset at the A.3-measured 0.25·L depth. All axis-aligned rects. **Phase
balance is exact by mirror symmetry** — both branches see identical
corners, so the R.6 corner-path error can only detune the *common* match,
never split the element phases; the S11 gate measures exactly that.

Unit gates: every trace AABB has its exact mirror partner about y = 0
(this caught a w50/2 offset in the −y transformer during development);
single connected component under edge-adjacency (the certified inset
construction's bands share exact arithmetic edges — contiguous metal under
rasterization; corner-point contact does not count); seed dims equal the
closed forms.

## Gate `engine-antenna-007` — GREEN first run

auto_dx-seeded (0.533 mm, h/3 binding), classic floor-ground stack, A.2
open-top boundary, single-run directional |S11|:

- dip at **2.450 GHz / −21.1 dB** vs the designed 2.45 GHz — **0.0 %**
  error (the 25 MHz raster hit f₀ exactly); the 0.5 λ₀ mutual-coupling
  detune is below the raster.
- Asserts pinned: frequency error ≤ 5 %, depth ≤ −10 dB.

The corporate tree is the real DUT here — with two matched elements, only
the transformer pair's junction match can fail, and it measured textbook.

## Gate `engine-antenna-008` — pattern multiplication

The reason arrays exist: in the array plane (y-z) the 2-element factor
AF(θ) = cos(π/2·sin θ) multiplies the element pattern (−3 dB at 30°,
−13.6 dB at 60°, horizon null), while the E-plane cut stays patch-like.
Asserts: array-plane θ = 60°/75° ≥ 10 dB below broadside on both sides;
E-plane broadside beats θ = 60° (the A.2 idiom). Measured numbers recorded
on the first green run.

## Consequences and queued follow-ons

The antenna catalog now spans broadside single element (A.0–A.3),
end-fire traveling-wave (FS.1a), and broadside array with feed-network
synthesis (FS.1b) — three of the four workhorse planar topologies.
Queued: FS.1b N×1 generalization (the tree recursion is mechanical once
the 2×1 junction is certified), FS.1c thin-wire subcell, FS.2 gain in
dBi / efficiency (the NTFF cut is |E|-only today — gain needs the
input-power normalization).
