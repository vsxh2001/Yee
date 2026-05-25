# ADR-0070 — Phase 4.fem.eig.3.5.6: Lee-Mittra First-Order Absorbing-Mode Wave-Port

**Status:** Accepted  
**Date:** 2026-05-25  
**Context Phase:** 4.fem.eig.3.5.6 (queued by ADR-0049)

## Context

`fem-eig-006` (high-aspect 100×10×1 mm rectangular cavity, WavePort at +x face, 40 GHz)
has `|S₁₁| = 0.955` under the current modal-projection wave-port. ADR-0049 queued the
Lee-Mittra 1997 §IV absorbing-mode port as the next fix.

### Diagnosis of the current wave-port stiffness

The existing `scatter_port_face_gauss` adds, for each mode m in the port:

```
K += jβ_m × assemble_port_face_block_gauss_pts(face, normal, β_m, 1.0)
```

`assemble_port_face_block_gauss_pts` computes the **full face Gram matrix** `B_face`:

```
B_face[i,j] = ∑_g w_g (n̂ × N_i(ξ_g)) · (n̂ × N_j(ξ_g))
```

This is independent of the modal shape — mode shapes enter only through the RHS.
For {TE₁₀, TE₂₀, TE₀₁} at 40 GHz (TE₀₁ evanescent → β=0 → block=0):

```
K_port = j(β₁₀ + β₂₀) B_face  ≈ j·1330 B_face
```

The port acts as a **scalar ABC** with β_eff ≈ 1330 rad/m for ALL modal content.
This is a poor absorber for modal content outside the {TE₁₀, TE₂₀} basis.

### Lee-Mittra first-order absorbing BC (§IV)

The correct Lee-Mittra port imposes mode-specific impedance matching:

```
K_absorbing[i,j] = jk₀ B_face[i,j] + ∑_m j(β_m − k₀) R_m[i,j]
```

where `R_m[i,j] = ∫_Γ [(n̂×N_i)·e_t^m] [(e_t^m·n̂×N_j)] dS` is the rank-1 modal
projection for mode m, and k₀ = ω/c₀.

- For modes in the basis (TE₁₀, TE₂₀): mode-specific β_m absorption
- For modal content in the complement: first-order k₀ ABC absorption
- For evanescent modes (β_m = 0): `β_m − k₀ = −k₀` (complement handles them)

## Decision

Implement the Lee-Mittra first-order absorbing-mode port as `PortDefinition::with_absorbing_complement()`:

1. New element-layer function `assemble_port_face_block_projected_gauss_pts` (exact Whitney-1
   Gauss quadrature; rank-1 outer product of modal projections).
2. New `assemble_port_face_block_projected` (centroid approximation for non-Gauss path).
3. `PortDefinition::absorbing_complement: bool` (default false → backward-compat).
4. New branch in `scatter_port_face_gauss` and `scatter_port_face` when absorbing_complement=true.
5. fem-eig-006 port_1 gains `absorbing_complement: true`.

## Measurement (to be filled after implementation)

| Configuration | `|S₁₁|(40 GHz)` |
|---|---|
| v3.5.5 multi-mode TE₁₀+TE₂₀+TE₀₁ (baseline) | 0.955397 |
| v3.5.6 with Lee-Mittra absorbing complement | TBD |

## Consequences

- `fem_eig_006_magnitude_bounded` is un-ignored if `|S₁₁| < 0.1`; kept `#[ignore]`
  with updated measurement if not.
- `PortDefinition::absorbing_complement = false` by default → all existing port
  definitions unchanged; all existing gates must remain green.
- The new `assemble_port_face_block_projected_gauss_pts` is a reusable higher-order
  absorbing BC primitive for future Phase 4 ports.
- If the gate does not close, the finding and new measurement are recorded here;
  the next step (Phase 4.fem.eig.3.5.7) would investigate higher-order absorbing BC.

## References

- Lee, M.-F., and R. Mittra. "Absorbing Boundary Conditions for Wave-port
  Excitation." *IEEE Trans. Microw. Theory Tech.* 45(7), 1997, §IV.
- Jin, J.-M. *The Finite Element Method in Electromagnetics*, 3rd ed. §10.5–10.6.
- ADR-0046 — Phase 4.fem.eig.3.5.3 wave-port termination.
- ADR-0048 — Phase 4.fem.eig.3.5.5 frequency-retune disposition.
- ADR-0049 — queues this ADR; records v3.5.5 baseline measurement.
- Spec: `docs/superpowers/specs/2026-05-25-phase-4-fem-eig-3-5-6-absorbing-mode-wave-port-design.md`
- Plan: `docs/superpowers/plans/2026-05-25-phase-4-fem-eig-3-5-6-absorbing-mode-wave-port.md`
