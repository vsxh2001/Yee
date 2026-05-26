# Phase 4.fem.eig.3.5.6 — Lee-Mittra Absorbing-Mode Wave-Port — Design Spec

**Status:** Draft  
**Owner:** TBD  
**Phase:** 4.fem.eig.3.5.6  
**ADR:** ADR-0070  
**Depends on:** Phase 4.fem.eig.3.5.5 (shipped; fem-eig-006 at 40 GHz, TE₂₀ propagating)  
**Closes:** `fem_eig_006_magnitude_bounded` gate (`|S₁₁| < 0.1`; `#[ignore]`'d since v3.5.3)

---

## 1. Context

`fem-eig-006` (high-aspect 100×10×1 mm cavity, 40 GHz) has `|S₁₁| = 0.955` under
the current modal-projection wave-port (WavePort). ADR-0049 queues the Lee-Mittra
first-order absorbing-mode wave-port (Lee-Mittra 1997 §IV) as the next fix.

The residual is confirmed NOT to be discretisation-limited (v3.5.5 refinement probe)
or modal-degeneracy-limited (v3.5.4 cutoff analysis): it is a genuine limitation of
the current port stiffness block, which adds a **scalar** (non-projected) `jβ_m B_face`
term for each mode.

### Diagnosis of the existing wave-port stiffness

The current `scatter_port_face_gauss` loop adds, for each mode m:

```
K += jβ_m × assemble_port_face_block_gauss_pts(face, normal, β_m, 1.0)
```

`assemble_port_face_block_gauss_pts` computes:

```
B_face[i,j] = ∑_g w_g (n̂ × N_i(ξ_g)) · (n̂ × N_j(ξ_g))
```

This is the **full face Gram matrix** — it does **not** project onto the modal shape e_t^m.
For a 3-mode basis {TE₁₀, TE₂₀, TE₀₁} at 40 GHz (TE₀₁ evanescent → β=0 → block=0):

```
K_port = j(β₁₀ + β₂₀) B_face
```

The port acts as a **scalar ABC** with β_eff = β₁₀ + β₂₀ ≈ 1330 rad/m for ALL modal
content at the face — including modes NOT in the {TE₁₀, TE₂₀} basis.

### Lee-Mittra first-order absorbing BC (§IV)

The correct absorbing-mode wave-port (Lee-Mittra 1997 §IV) imposes mode-specific
impedance matching on the port face Γ_port:

```
K_absorbing[i,j] = ∑_m jβ_m ∫_Γ (n̂ × N_i) · P_m · (n̂ × N_j) dS
                 + jk₀  ∫_Γ (n̂ × N_i) · (I - P) · (n̂ × N_j) dS
```

where P_m = e_t^m ⊗ e_t^m (dyadic modal projection), P = ∑_m P_m, and k₀ = ω/c₀.

Equivalently:

```
K_absorbing[i,j] = jk₀ B_face[i,j] + ∑_m j(β_m - k₀) R_m[i,j]
```

where `R_m[i,j] = ∫_Γ [(n̂×N_i)·e_t^m] [(e_t^m·n̂×N_j)] dS` is the rank-1 modal
projection block for mode m.

**Effect:**
- Modes in basis: absorbed with mode-specific β_m (correct impedance matching)
- Modes NOT in basis (complement): absorbed with k₀ (first-order ABC)
- Modes with β_m > k₀ (strong propagation): complement reduces over-absorption of
  non-projected content
- Evanescent modes (β_m = 0): term vanishes, complement absorbs them via k₀

---

## 2. Scope

**In scope:**

- New function `assemble_port_face_block_projected_gauss_pts` in `crates/yee-fem/src/element.rs`
  computing the rank-1 projected block R_m (exact Whitney-1 / 3-pt Gauss).
- New function `assemble_port_face_block_projected` in `element.rs` computing R_m via
  the centroid approximation (for the non-gauss path in `scatter_port_face`).
- `PortDefinition::absorbing_complement: bool` field (default: `false`); backward-compat.
- `PortDefinition::with_absorbing_complement()` builder.
- New branch in `scatter_port_face` and `scatter_port_face_gauss` when `absorbing_complement = true`:
  - REPLACES the scalar `jβ_m B_face` stiffness accumulation with the Lee-Mittra formula.
  - RHS accumulation (a_inc × 2jβ_m × ∫N_i·e_t dS) is UNCHANGED.
- fem-eig-006 fixture (`run_fem_eig_006_high_aspect_pml_with_config` in
  `crates/yee-validation/src/lib.rs`): port_1 gains `absorbing_complement: true`.
- Gate disposition: if `|S₁₁|(40 GHz) < 0.1`, un-ignore `fem_eig_006_magnitude_bounded`.
  If not, surface the measured value as a finding and keep `#[ignore]` with updated notes.
- Unit tests for `assemble_port_face_block_projected_gauss_pts`:
  - Rank-1 structure check (non-zero entry forces outer-product structure)
  - Uniform e_t = (1,0,0): projected block matches scaled scalar block when entire face
    normal-cross aligns with e_t (degenerate case)
  - Orthogonal e_t = (0,0,1) to face: projected block is zero
  - `absorbing_complement = false` path: backward-compat canary — fem-eig-004 thru-line
    `|S₂₁|` and fem-eig-001 cavity frequency unchanged (regression guard)

**Out of scope:**

- Higher-order absorbing BC (Lee-Mittra 1997 §V rational-function extension)
- Updating `scatter_absorbing_complement` to handle CFS-PML combinator
- Changing S-parameter extraction logic (M_pp normalisation, `extract_s11`)
- Modifying the RHS assembly in any way
- Multi-port S-matrix (only `S₁₁` extracted for fem-eig-006)
- Python bindings change (no new public surface on `yee_py`)

---

## 3. Physics and Implementation Detail

### 3.1 Existing stiffness path (non-absorbing_complement)

For a port face with N modes:

```
// for each mode m:
K += sign_ij × jβ_m × B_face[i,j]
// ← unchanged, backward-compat
```

### 3.2 New absorbing_complement stiffness path

When `port.absorbing_complement = true`:

```
// Step A: compute k₀ × full face Gram
let k0 = omega / C0;
let K_full = assemble_port_face_block_gauss_pts(verts, normal, Complex64::new(k0, 0.0), 1.0);

// Step B: for each mode m, compute mode-projected rank-1 block scaled by (β_m - k₀)
let mut K_correction = SMatrix::<Complex64, 3, 3>::zeros();
for mode in port.modes {
    let beta = (mode.beta_mode)(omega);
    let e_t_gauss = [mode.modal_e_t(p_0), mode.modal_e_t(p_1), mode.modal_e_t(p_2)];
    let R_m = assemble_port_face_block_projected_gauss_pts(verts, normal, beta, &e_t_gauss, 1.0);
    K_correction += Complex64::new(beta - k0, 0.0) * R_m; // (β_m - k₀) correction
}
// Net Lee-Mittra block = jk₀ B_face + j∑_m (β_m - k₀) R_m = K_full + j*K_correction
let K_lee_mittra = K_full + Complex64::i() * K_correction;
scatter K_lee_mittra into triplets;
```

Wait — to be precise, `assemble_port_face_block_gauss_pts` already includes the `j` factor
(it returns `jβ × B_face`). So the exact formulation:

```
K_lee_mittra = assemble_port_face_block_gauss_pts(verts, normal, Complex64::new(k0, 0.0), 1.0)
             + ∑_m assemble_port_face_block_projected_gauss_pts(verts, normal, β_m - k₀, e_t_m, 1.0)
```

where `assemble_port_face_block_projected_gauss_pts(verts, normal, β_eff, e_t_gauss, mu_r)` returns:
`jβ_eff × ∑_g w_g [(n̂ × N_i(ξ_g)) · e_t_g] [(e_t_g · n̂ × N_j(ξ_g))]`

Note: when β_eff = 0 (evanescent mode), the term vanishes identically (same as existing).

### 3.3 `assemble_port_face_block_projected_gauss_pts` function

Signature:
```rust
pub fn assemble_port_face_block_projected_gauss_pts(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_eff: f64,               // may be β_m, (β_m - k₀), etc.
    modal_e_t_at_gauss_pts: [Vector3<f64>; 3],
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3>
```

Formula (exact Whitney-1, 3-point Gauss):
```
block[i,j] = j · β_eff / μ_r × ∑_g w_g · [(n̂ × N_i(ξ_g)) · e_t_g] · [(e_t_g · n̂ × N_j(ξ_g))]
```

`w_g = face_area / 3.0` (equal Gauss weight for 3-pt midpoint quadrature on triangle).

Outer product = rank-1: `block[i,j] = j β_eff μ_r⁻¹ × (∑_g w_g a_i(ξ_g) a_j(ξ_g))`
where `a_i(ξ_g) = (n̂ × N_i(ξ_g)) · e_t_g`.

### 3.4 Centroid-approximation variant (non-gauss path `scatter_port_face`)

Mirror the above for `assemble_port_face_block_projected` (centroid quadrature):

```
block[i,j] = j · β_eff / μ_r × face_area · [(n̂ × t_i) · e_t_c] · [(e_t_c · n̂ × t_j)]
```

where `e_t_c = mode.modal_e_t(face_centroid)`.

---

## 4. Validation Gates

### Gate A (primary, un-ignore target)
```
fem_eig_006_magnitude_bounded: |S₁₁(40 GHz)| < 0.1
```
Must pass after the absorbing-complement change on port_1.

### Gate B (existing, must remain green)
```
fem_eig_006_no_nan_inf: S₁₁ is finite
```
Must pass.

### Gate C (backward-compat regression guards, must remain green)
```
fem_eig_001: TE₁₀₁ WR-90 cavity freq ≤ 0.09% rel err
fem_eig_002: lossy SiO₂ cavity Re(f) ≤ 1.3e-3, Im(f) ≤ 2.96e-3
fem_eig_004: |S₂₁| ≈ -0.045 dB, reciprocity 2e-15
fem_eig_005: passivity + reciprocity
```

---

## 5. Escape Hatch

If `|S₁₁|` does not reach 0.1 after the Lee-Mittra complement:
- Record the measured `|S₁₁|` in the test docstring and in ADR-0070
- Keep `fem_eig_006_magnitude_bounded` under `#[ignore]` with updated notes
- Surface the finding to the orchestrator with: measured value, % improvement from 0.955,
  and a root-cause hypothesis (e.g., centroid vs Gauss sampling of e_t, sign convention,
  or genuinely a higher-order absorbing BC needed)
- Do NOT weaken the `< 0.1` tolerance

**Blocked > 15 min on any single step → surface and stop.**

---

## 6. Files Touched (lane)

Lane: `crates/yee-fem/src/**`, `crates/yee-validation/src/lib.rs`,
      `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`,
      `docs/superpowers/specs/`, `docs/superpowers/plans/`, `docs/src/decisions/`,
      `docs/src/SUMMARY.md`, `ROADMAP.md`, `CLAUDE.md`

Out of lane: yee-py, yee-gui, yee-cli, yee-fdtd, yee-mom, CI workflows — surface if found.

---

## 7. References

- Lee, M.-F., and R. Mittra. "Absorbing Boundary Conditions for Wave-port
  Excitation." *IEEE Trans. Microw. Theory Tech.* 45(7), 1997, §IV.
- Jin, J.-M. *The Finite Element Method in Electromagnetics*, 3rd ed.
  Wiley 2014, §10.5–10.6 (modal wave-port + termination).
- ADR-0046 — Phase 4.fem.eig.3.5.3 W1 wave-port termination.
- ADR-0048 — Phase 4.fem.eig.3.5.5 frequency-retune disposition.
- ADR-0049 — Phase 4.fem.eig.3.5.6 absorbing-mode port queued.
