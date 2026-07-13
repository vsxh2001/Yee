# ADR-0221: FS.4.1 — vias through multilayer stackups

**Date:** 2026-07-13 · **Status:** accepted · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-13-fs41-stackup-vias-design.md`

## Context

FS.4.0 (ADR-0215) gave the FDTD flow N-layer stackups, buried traces and
lids; the only via primitive was the single-layer ground-to-trace
`with_via_at_cell` (R.1, ADR-0194). The engine protocol already carries
`pec_mask_ez` end-to-end, so — like FS.4.0 itself — the gap was purely
the voxel-side helper.

## Decision

1. **`yee_voxel::with_via_between(model, i, j, k_lo, k_hi)`** — blind
   via: `E_z` edges `k_lo..k_hi` at grid column `(i, j)` become PEC (a
   post from node-plane `k_lo` to node-plane `k_hi`). Grid **cell**
   indices, not stackup layer indices — the caller quantizes layer
   heights the same `round(h/dx).max(1)` way `voxelize_stackup` did;
   the cell-index post-processing idiom is kept.
2. **`with_through_via_at_cell(model, i, j)`** — through-via: the full
   stack, `with_via_between(…, 0, nz)`; on a lidded stackup node `nz`
   *is* the lid, so this is a ground-to-lid barrel through every layer
   and any trace it passes.
3. **`with_via_at_cell`** is re-expressed as
   `with_via_between(…, 0, k_top)` — bit-identical (gated), so R.1's
   `engine-via-001` and unit test are untouched.

## Measured gates

- **`voxel-stackup-002`** (instant, GREEN): on the FS.4.0 3-layer lidded
  stack (2+3+2 cells, trace k = 5, lid k = 7) a through-via masks
  exactly `E_z` `k = 0..7` and a blind via `with_via_between(…, 2, 5)`
  exactly `k = 2..5`; the **whole-mask set count equals 7 + 3** (nothing
  else touched anywhere); neighbour columns clear; `with_via_at_cell` ≡
  `with_via_between(…, 0, k_top)` bit-identical.
- **`engine-stackup-via-001`** (release, ignored, GREEN): symmetric
  stripline (ε_r = 4.4, b = 3.2 mm = 16 cells at dx = 0.2 mm, lid on)
  with a mid-line λ/4 stub, three runs on one 477×88×16 grid (12000
  steps each, **~5.5 min total** release CPU):
  - open stub: |S21| notch **−39.81 dB at 5.075 GHz** (design 5 GHz
    after pre-compensating the stub by the stripline open-end
    Δl = b·ln2/π — 1.5 % off);
  - same stub shorted by a **through-via** at its far end: **+0.62 dB**
    at 5.075 GHz and a whole-band (4.0–5.6 GHz) min of **−1.18 dB** —
    the notch is gone everywhere, exactly the shorted-λ/4-stub
    (input-open, own notch near 2 f₀) behaviour.
  - Pinned: control ≤ −15 dB; via variant ≥ −3 dB at the control notch
    **and** everywhere in the band. (The small positive dB is the
    ADR-0204 single-ratio launch artifact — fine for notch-shaped
    asserts, not for absolute |S21|.)

## Fixture lessons (mapping the ADR-0215 hygiene to a two-port stripline)

- **b ≥ 16 cells** carries over unchanged (confined lidded mode).
- The **box-mode cutoff rule** dissolves when the lateral walls are CPML
  instead of PEC: with absorbing side walls there is no lateral cavity,
  and parallel-plate waves radiated by the tee/via are absorbed. The
  stub end sits 14 working cells (2.8 mm) inside the absorber
  (exp(−πx/b) fringing ≈ 6 %).
- The **pulse-tail/time-gate rule** is replaced by absorbing
  terminations (resistive aperture ports + CPML beyond) and a
  ring-down-length record (12000 steps ≈ 4.6 ns).

## Non-goals / queued

Finite via barrel diameter + antipads/pads; a layer-index convenience
over `Stackup`; via inductance vs closed form; multi-trace-layer
voxelization; automesh awareness of vias (FS.4.2). Roadmap FS.4 row
update deferred to the merge (out of the track's lane).
