# FS.1b — 2×1 corporate-fed patch array

**Date:** 2026-07-08
**Track:** FULL-SUITE-ROADMAP FS.1 (antenna catalog), after FS.1a
(quasi-Yagi) shipped complete. Arrays are the first "more gain by
replication" step every commercial antenna tool ships (pattern
multiplication + corporate feed synthesis).

## Decomposition

- **FS.1b.0 (walking skeleton, this spec):** `yee_layout::patch_array_2x1`
  — two A.0/A.3-certified inset-fed patches side by side along y (H-plane
  pair) at 0.5 λ₀ centre spacing, fed in phase by a symmetric corporate
  tree: 50 Ω spine → junction → two **λg/4 70.7 Ω transformers** along ±y
  (each transforms 50 → 100 Ω; the two in parallel present 50 Ω at the
  junction) → 50 Ω branches → each patch's inset (the A.3-measured
  0.25·L depth). All axis-aligned rects with overlapped joints; the
  classic full-ground stack (no lifted stack needed — patches radiate up).
  Unit gates: mirror symmetry about y = 0, connectivity (single component,
  union-find), patch rects reproduce `patch_antenna_dims`, transformer
  length = λg(70.7 Ω)/4.
- **FS.1b.1 gate `engine-antenna-007`** (release, 1 solve, A.1 machinery):
  the array's directional |S11| dips within ±10 % of the designed
  2.45 GHz with a depth tripwire → pinned from measurement. The corporate
  tree is the DUT here: a wrong transformer leaves both patches
  mismatched at the junction.
- **FS.1b.2 gate `engine-antenna-008`** (release, 1 solve + NTFF): array
  physics — in the **array plane** (the y-z cut) the 2-element factor
  narrows the beam: broadside beats θ ≥ 45° by a measured-then-pinned
  margin, and the E-plane cut (x-z) stays patch-like. Asserts chosen
  after the first instrumented run (A-track pattern).

## Risks

- Corner path errors in the tree (R.6 κ): symmetric for the pair (both
  branches see identical corners), so phase BALANCE is safe by symmetry —
  only the common match can detune. The S11 gate catches it.
- Mutual coupling at 0.5 λ₀ shifts the resonance a little vs the single
  patch (literature: a few %); the ±10 % gate absorbs it — record the
  measured shift.
