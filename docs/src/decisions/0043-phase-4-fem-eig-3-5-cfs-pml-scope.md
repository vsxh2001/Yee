# ADR-0043: Phase 4.fem.eig.3.5 scope — CFS-PML retires fem-eig-003 strict gate

## Status

Accepted — 2026-05-20 (spec + plan; implementation deferred to
follow-up tracks).

## Context

Phase 4.fem.eig.3 (ADR-0042) shipped the 2nd-order Engquist-Majda ABC
on `OpenBoundarySolver`, lowering the normal-incidence reflection floor
from `~ −40 dB` (1st-order) toward the `~ −60 dB` Jin §10.4 Table 10.1
figure. The fem-eig-003 strict gates (`[-45, -35] dB` absorption floor
and strict `|S_{11}| < 1` continuum bound) remained `#[ignore]`'d under
the Phase 4.fem.eig.2 E5 escape hatch.

Track NNNNNNNNN (2026-05-20, commit `9cf4b3a`, merge `fc62f82`) refined
the fem-eig-003 cavity mesh from `(16, 8, 24) = 18 432 tets` to
`(24, 12, 36) = 62 208 tets` and re-ran the swept driver with v3's
F1+F2 coupled exact-Whitney-1 modal RHS + projection plus F3+F4
2nd-order Engquist-Majda ABC. The measurement: `|S_{11}(f)| ∈ [0.9976,
0.99997]` ⇒ `s11_db ∈ [-2.22e-2, -2.86e-5] dB` across 8-12 GHz. The
refined mesh is ~2× better in dB than the JJJJJJJJJ
`(16, 8, 24)` baseline (`s11_db ∈ [-5.0e-2, -8.1e-5] dB`) but **still
~35 dB above** the spec §8 `[-45, -35] dB` window. Per the Track
NNNNNNNNN brief escape hatch ("strict gate still fails > 5 dB above
−35 dB → fundamental limit reached; queue Phase 4.fem.eig.3.5 PML"),
both strict gates remain `#[ignore]`'d. The binding constraint at the
`(24, 12, 36)` mesh tier is no longer modal-sampling discretisation
but the 2nd-order Engquist-Majda ABC's **intrinsic floor for
off-normal modal content scattered by the truncation surface** —
exactly the deferral path ADR-0042 §risks queued for v3.5.

The CFS-PML upgrade was historically queued at three different deferral
points in Phase 4.fem.eig.{2,3} ADRs (0040 §C-3, 0042 §risks). Phase
4.fem.eig.3.5 takes that slot and ships the CFS-PML kernel that
finally retires the strict gates without weakening tolerances.

## Decision

Phase 4.fem.eig.3.5 ships CFS-PML (Complex Frequency Shifted Perfectly
Matched Layer, Roden-Gedney 2000) as a new `AbcOrder::CfsPml(PmlConfig)`
variant on `OpenBoundarySolver`. The PML is a thin (6 tet-layer
default) volumetric buffer outside the original cavity volume, with
stretched-coordinate complex permittivity tensor `ε_eff(ω) = ε · Λ(ω)`
where

```text
    Λ(ω)    =  diag( s_y s_z / s_x,  s_z s_x / s_y,  s_x s_y / s_z ),
    s_α(ω)  =  κ_α(d_α)  +  σ_α(d_α) / ( α_α  +  j ω ε_0 ),
```

per Roden-Gedney 2000 §II. The CFS modification (`α_α > 0`,
Kuzuoglu-Mittra 1996) retires the evanescent-mode divergence of the
original Berenger 1994 PML. The original `FaceKind::Abc` truncation
surface is replaced by a volumetric PML buffer; the new outermost
extended-mesh face is tagged `FaceKind::Pec` and absorbs via PML
round-trip attenuation `R(θ) ≈ exp(−2 cos(θ) · ∫_0^D σ_eff dd /
(ε_0 c)) ≈ exp(−12) ≈ −105 dB` for the default `σ_max · D / (ε_0 c)
≈ 6` rule of thumb.

Six load-bearing decisions:

1. **CFS-PML, not the original Berenger 1994 PML.** The CFS
   modification `α_α > 0` retires the evanescent-mode divergence
   (Roden-Gedney 2000 §IV). The original Berenger PML is not an option
   for the fem-eig-006 high-aspect-ratio fixture (which excites
   exactly the evanescent-mode regime CFS was designed to absorb).
2. **Stretched-coordinate anisotropic-`ε` formulation, not
   split-field.** Roden-Gedney 2000 §II reframes PML as a complex
   anisotropic permittivity tensor `Λ(ω)`, avoiding the 12-component
   split-field unknowns of Berenger 1994. For frequency-domain FEM,
   no convolution is required — the `ε · Λ(ω)` form is direct, and
   per-tet assembly is a straightforward extension of
   `assemble_tet_element_complex` to per-tet 3×3 complex `ε_tensor`
   matrices.
3. **6-cell PML thickness default.** Roden-Gedney 2000 §III's
   recommended "6 to 10 cells for microwave applications". Configurable
   via `PmlConfig::thickness_cells`; ablation is deferred to Phase
   4.fem.eig.3.5.1.
4. **Cartesian-aligned PML only.** The Cartesian-aligned case keeps
   `Λ(ω)` diagonal in the global frame, preserving the
   complex-symmetric structure of `K(ω) − k₀² M(ω)` and the
   `faer::sparse::FaerLuSolver<Complex64>` complex-LDLᵀ factorisation
   path. Non-Cartesian-aligned faces are rejected with
   `Error::InvalidArgument` until Phase 4.fem.eig.3.5.1.
5. **Existing `AbcOrder::{First, Second}` paths stay.** CFS-PML is a
   new `AbcOrder::CfsPml` variant alongside the 1st/2nd-order surface
   ABC kernels, not a replacement. `AbcOrder::Second` remains the
   default for v3.5 callers; only fem-eig-003-strict and fem-eig-006
   flip to `CfsPml`. The change is additive; v0/v1/v2/v3 paths stay
   bit-for-bit identical on `AbcOrder::{First, Second}`.
6. **fem-eig-003 strict gates un-ignore + fem-eig-006 new fixture.**
   Both `#[ignore]`'d gates in
   `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` become
   CI-default. New fem-eig-006 fixture (100 mm × 10 mm × 1 mm
   high-aspect cavity, TE-mode drive at 30 GHz) stress-tests PML
   stability on off-normal modal content — the regime where v3's
   2nd-order ABC saturates at `|S_{11}| ≈ 0.95`. Gate: `|S_{11}(30
   GHz)| < 0.1`.

CPU-only, single-threaded, FP64 complex. No GPU. Cartesian-aligned PML
only. `faer::sparse::FaerLuSolver<Complex64>` continues to handle the
complex-symmetric matrix unchanged (diagonal-`Λ` preserves complex
symmetry).

## Consequences

- **fem-eig-003 strict gates clear without weakening tolerances.** The
  un-ignore in P5 is a single attribute removal; the underlying physics
  fix is the CFS-PML volumetric absorption replacing the
  Engquist-Majda surface integral.
- **`AbcOrder` enum gains a `CfsPml(PmlConfig)` variant.** The existing
  `First` and `Second` variants are unchanged. v0-v3 callers continue
  to use the surface-integral ABC kernels bit-for-bit; only
  fem-eig-003-strict and fem-eig-006 (and future PML callers) flip to
  `CfsPml`.
- **`assemble_tet_element_complex` gains an anisotropic-tensor
  sibling.** D1+D3+D4 (Phase 4.fem.eig.{0,1,2} scalar-`ε` paths) all
  remain reachable bit-for-bit; the new
  `assemble_tet_element_complex_anisotropic` is called only on PML
  tets. The scalar entry point stays.
- **High-aspect-ratio cavity analysis becomes a first-class FEM
  capability.** fem-eig-006 stress-tests the regime where v3's
  2nd-order ABC fails; future microwave / millimetre-wave thin-film
  filter validation lands on the same CFS-PML surface.
- **Rotated-PML / non-Cartesian-aligned faces deferred to v3.5.1.**
  The implementation rejects non-axis-aligned PML face normals with
  `Error::InvalidArgument`. This is a known restriction; the
  cuboid-cavity fem-eig-003 and the axis-aligned high-aspect
  fem-eig-006 both satisfy it.
- **Complex-symmetric stiffness is preserved under Cartesian-aligned
  PML.** Diagonal `Λ(ω) = diag(λ_x, λ_y, λ_z)` keeps the assembled
  `K(ω) − k₀² M(ω)` complex-symmetric, so `faer`'s complex LDLᵀ path
  is unchanged. Off-diagonal `Λ` (deferred) would force complex LU
  with ~2× factorisation cost.
- **PML grading parameters become a tunable surface.** `PmlConfig`
  exposes `thickness_cells`, `sigma_max`, `alpha_max`, `kappa_max`,
  and `m`. Sensible Roden-Gedney 2000 Table-I defaults apply; an
  ablation sweep is queued for Phase 4.fem.eig.3.5.1.

## References

- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
  — Phase 4.fem.eig.3.5 design spec.
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-cfs-pml.md` —
  P1-P7 implementation plan.
- ADR-0042 — Phase 4.fem.eig.3 scope (this ADR's parent; §risks
  deferral path this ADR fulfils).
- ADR-0040 — Phase 4.fem.eig.2 open-boundary scope (grandparent;
  §C-3 originally queued PML to Phase 4.fem.eig.2.5).
- J.-P. Berenger, "A perfectly matched layer for the absorption of
  electromagnetic waves", *J. Comput. Phys.* 114 (1994), pp. 185-200,
  DOI 10.1006/jcph.1994.1159 — the original split-field PML
  reference.
- M. Kuzuoglu and R. Mittra, "Frequency dependence of the constitutive
  parameters of causal perfectly matched anisotropic absorbers",
  *IEEE MWCL* 6(12) (1996), pp. 447-449, DOI 10.1109/75.541428 — the
  CFS modification.
- J. A. Roden and S. D. Gedney, "Convolutional PML (CPML): An
  efficient FDTD implementation of the CFS-PML for arbitrary media",
  *IEEE MWCL* 10(5) (May 2000), pp. 27-29 — the stretched-coordinate
  anisotropic-material formulation this spec implements directly in
  frequency domain.
- J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley 2014, §10.4 (ABC reflection floors), §10.8 (PML for FEM).
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` —
  §"NNNNNNNNN status" measurement that motivates this ADR.
- CLAUDE.md §3, §4.
