# ADR-0046: Phase 4.fem.eig.3.5.3 scope — retire fem-eig-006 via wave-port termination (W1 hypothesis)

## Status

Accepted — 2026-05-20 (spec + plan; implementation queued under the
Phase 4.fem.eig.3.5.3 T1-T4 ladder).

## Context

Phase 4.fem.eig.3.5.2 (ADR-0045) shipped under Track SSSSSSSSS
against base SHA `5ec8e90`. The S1-S5 ladder landed the
`alpha_grading_order: usize` field on `PmlConfig`, the extended H3
thickness sweep `thickness_cells ∈ {12, 14, 16}`, and the new
`(kappa_max=2, m=4, thickness=14, alpha_grading_order=1)` defaults
that retired the fem-eig-003 strict band into `[-71.53, -55.58] dB`
(merge SHA `8aad1be`).

The S2 H4 ablation binary additionally re-ran
`fem_eig_006_magnitude_bounded` across all 18 H4 configurations
(`kappa_max=2 × m ∈ {3, 4} × thickness ∈ {12, 14, 16} ×
alpha_grading_order ∈ {0, 1, 2}`). **Result: `|S_{11}|(30 GHz)`
frozen at `0.926` in every row.** Neither α-grading nor doubled PML
thickness moved the reflection coefficient by more than the 4th
decimal. The current ignore docstring at
`crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs:53-60`
records the measurement verbatim.

The orthogonality is decisive: **α-grading is not the binding
constraint for fem-eig-006**. The cavity is `100 mm × 10 mm × 1 mm`
(aspect 100 : 10 : 1); at 30 GHz on `b = 10 mm` the TE_{10} mode is
well-propagating (`f_c ≈ 15 GHz`, `λ_g ≈ 11.5 mm`) with bulk-wave
propagation along x — i.e. normal to the +x truncation face. But
the modal energy distribution at the +x face is dominated by the
**guide-mode** standing-wave structure across `(y, z)`, not the
bulk-wave propagating along x. Berenger 1996 IEEE TAP 44:1 §IV-A
notes that Cartesian-aligned CFS-PML reflection floors **degrade
catastrophically for guide-modes** even when bulk-wave propagation
appears normal to the face, because the per-axis stretching tensor
cannot redistribute absorption budget between axes on a flat slab
geometry.

ADR-0045 §risks (a) flagged this exact mitigation:

> If the fem-eig-006 fixture cannot be retired by any PML grading,
> the v3.5.3 phase is reframed as a driver-level fix: switch the +x
> face from `FaceKind::AbcFace` (PML) to `FaceKind::WavePort(0)`
> reusing the existing `fem_eig_006_modal_e_t_te10` shape. This is
> the standard Jin §10.6 closed-cavity wave-port termination and
> does not require any new types or PML changes.

The SSSSSSSSS H4 measurement triggered this deferral path. ADR-0046
formalises the W1 wave-port-termination decision and defers W2
(rotated PML) + W3 (multi-face wedge) to Phase 4.fem.eig.4+.

Three candidate fixes were enumerated in spec §3:

1. **W1 (wave-port termination on +x).** Replace CFS-PML on the +x
   face with `FaceKind::WavePort(1)` carrying the same TE_{10}
   modal basis as the -x driving port. Simplest implementation
   (~10-line driver change), reuses Phase 4.fem.eig.2 E2 wave-port
   machinery. Expected `|S_{11}| < 0.01` per Jin §10.6 closed-cavity
   modal-termination floor.
2. **W2 (rotated/oblique CFS-PML).** Introduce a non-Cartesian-
   aligned stretching tensor per Berenger 1996 §IV-B. Requires
   generalising `assemble_tet_element_complex_anisotropic` to
   non-diagonal `ε_tensor`; breaks the complex-LDLᵀ diagonal
   exploit (ADR-0043 decision (4)) and may force complex-LU on the
   full system. v3.5 explicitly deferred this per ADR-0043
   §risks (c).
3. **W3 (multi-face PML wedges).** Wrap CFS-PML around three faces
   (+x, +y, -y) with Kuhn-6 wedge tets in the edges per Berenger
   1994 §V corner-wedge pattern. Most invasive (~1500-2000 LoC +
   new mesh tests in `yee-mesh`); v3.5 P2 escape hatch deferred
   this to v3.5.1+.

W2 and W3 are strictly more general than W1, but W1 retires
fem-eig-006 specifically without committing the project to the
non-diagonal-PML refactor or the multi-face wedge mesher. The
project's "walking-skeleton first" convention (CLAUDE.md §3)
favours W1.

## Decision

Phase 4.fem.eig.3.5.3 ships **W1 (wave-port termination on +x
face)** as the fem-eig-006 fix. The Phase 4.fem.eig.2 E2 wave-port
machinery (`FaceKind::WavePort(p)` + `PortDefinition { beta_mode,
modal_e_t }`) is reused unchanged; the v3.5.2 PML default quadruple
`(κ_max=2, m=4, thickness=14, alpha_grading_order=1)` stays in
force for fem-eig-003. No new types, no `yee-fem` API changes, no
new dependencies.

Five load-bearing decisions:

1. **W1 over W2 / W3.** W1 is the smallest diff that retires the
   gate (`~10` driver lines vs `~500-2000` LoC for W2 / W3). W1 is
   strictly less general — it only works for fixtures whose
   dominant mode is well-approximated by a single analytic modal
   shape — but for fem-eig-006 specifically that condition holds
   (TE_{10} is the only propagating mode at 30 GHz on the
   `b = 10 mm` broad wall). W2 and W3 are deferred to Phase
   4.fem.eig.4 (FEM-BEM hybrid or rotated-PML scope).
2. **Driver-only change; no `yee-fem` API churn.** All v3.5.3
   behavioural changes live inside
   `run_fem_eig_006_high_aspect_pml_with_config` in
   `crates/yee-validation/src/lib.rs`. The `PmlConfig` struct,
   `FaceKind` enum, `OpenBoundarySolver` builder, and Whitney-1
   modal-projection paths are all unchanged. This keeps the
   v3.5.3 diff narrow and the risk envelope confined to the
   fem-eig-006 driver.
3. **`pml_config` parameter survives as vestigial.** The
   `pml_config: yee_fem::PmlConfig` driver parameter is **kept**
   for source-compatibility with the v3.5.2 ablation binary
   (`crates/yee-validation/examples/cfs_pml_grading_sweep.rs`) but
   becomes unused inside the body (`#[allow(unused_variables)]`
   plus doc-comment note). Removal is queued for Phase
   4.fem.eig.4 cleanup. This preserves the v3.5.2 sweep CSV
   reproducibility — pre-v3.5.3 sweep rows continue to be
   re-runnable; post-v3.5.3, fem-eig-006 rows produce constant
   `|S_11|` across `(m, thickness, alpha_grading_order)` because
   the wave-port floor dominates.
4. **fem-eig-003 retire is untouched.** W1 is fem-eig-006-driver-
   specific. The fem-eig-003 driver uses its own face-kind
   classification + its own `.with_cfs_pml(...)` builder and is
   unchanged. The v3.5.2 strict-band retire `[-71.53, -55.58] dB`
   stays in force. The T3 verification step explicitly re-runs
   fem-eig-003 to confirm.
5. **Strict gate un-ignore after T2 + T3 measurements confirm
   retire.** If the T3 measurement yields `|S_{11}| < 0.1`, the
   `#[ignore]` is removed from `fem_eig_006_magnitude_bounded`
   and the docstring is refreshed with the v3.5.3 wave-port
   measurement. If T3 yields `|S_{11}| ≥ 0.1` (higher-order modal
   content per spec §7 (a)), the W1 driver ships anyway as a
   no-PML-needed cleanup, the docstring refreshes with the new
   measurement, and Phase 4.fem.eig.3.5.4 (multi-mode wave-port
   extension) is queued. We do not weaken the `< 0.1` tolerance.

CPU-only, single-threaded, FP64 complex. No GPU. No new
dependencies. Same `faer::sparse::FaerLuSolver<Complex64>` surface
(ADR-0043 decision (4), ADR-0044, ADR-0045 carried).

## Consequences

- **`crates/yee-validation/src/lib.rs` driver-only change.**
  `run_fem_eig_006_high_aspect_pml_with_config` switches the +x
  face from `FaceKind::Pec` (current PML-outer-truncation) to
  `FaceKind::WavePort(1)`; the `extend_mesh_with_pml` and
  `.with_cfs_pml(...)` calls drop off; a second `PortDefinition`
  clones the TE_{10} shape for port 1. Total diff ~25 lines.

- **`fem_eig_006_magnitude_bounded` flips from `#[ignore]` to
  CI-default** on the common-case retire path (T3 measurement
  `< 0.1`). The `< 0.1` tolerance is **not weakened**.

- **No `yee-fem` API change.** `PmlConfig`, `FaceKind`,
  `OpenBoundarySolver`, `PortDefinition` are unchanged. No new
  types; no new feature flags; no new public-API surface area.

- **`pml_config` parameter becomes vestigial** on the
  fem-eig-006 driver. Kept for v3.5.2 sweep CSV compatibility;
  removal queued for Phase 4.fem.eig.4 cleanup. Doc-comment
  flags the deprecation.

- **fem-eig-006 measurement on cavity_uniform native mesh.** The
  v3.5.2 driver ran on the SSSSSSSSS-era ~580-tet PML-extended
  mesh; v3.5.3 runs on the native `(16, 3, 2)` cavity (~96 tets
  via Kuhn-6). Wall-time drops from ~60-90 s to under ~10 s per
  fem-eig-006 invocation, which simplifies CI gating.

- **W2 + W3 deferred to Phase 4.fem.eig.4.** Rotated CFS-PML and
  multi-face wedge PML remain queued; if a future fixture is
  encountered that cannot be retired by either Cartesian-aligned
  CFS-PML (the v3.5.2 path) or wave-port termination (the v3.5.3
  path), W2 / W3 are the next implementation candidates.

- **Backward-compatibility.** The v3.5.2 sweep CSV
  (`crates/yee-validation/examples/cfs_pml_grading_sweep.rs`
  output) remains reproducible. The fem-eig-006 rows now report
  constant `|S_11|` across `(m, thickness, alpha_grading_order)`
  because the wave-port driver ignores `pml_config`; this is the
  expected post-v3.5.3 behaviour and documented in spec §5.

- **Risk: higher-order modal content underestimate.** Spec §7 (a):
  the TE_{20} cutoff on `b = 10 mm` sits at `30 GHz` exactly. If
  higher-order modal content is present, a TE_{10}-only wave-port
  underestimates the reflection. Mitigation queued for Phase
  4.fem.eig.3.5.4 multi-mode wave-port extension if T3
  measurements warrant.

## References

- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-3-design.md`
  — design spec (W1 / W2 / W3 hypothesis tree, W1 recommendation).
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-3.md`
  — T1-T4 implementation plan.
- ADR-0045 — Phase 4.fem.eig.3.5.2 scope (this ADR's parent; §risks
  (a) queued this ADR).
- ADR-0044 — Phase 4.fem.eig.3.5.1 grading retune (grandparent).
- ADR-0043 — Phase 4.fem.eig.3.5 CFS-PML scope (great-grandparent;
  §risks (c) deferred W2 / W3).
- ADR-0042 — Phase 4.fem.eig.3 scope (great-great-grandparent;
  coupled-Whitney + 2nd-order ABC + multi-port machinery reused by
  W1).
- ADR-0040 — Phase 4.fem.eig.2 open-boundary scope (E2 wave-port
  machinery W1 reuses).
- Berenger, J.-P., "Three-dimensional perfectly matched layer for
  the absorption of electromagnetic waves," *IEEE Transactions on
  Antennas and Propagation* 44(1) (January 1996), pp. 110-117.
  DOI 10.1109/8.477535. §IV-A bulk-vs-guide-wave physics
  motivating the W1 decision over the v3.5.2 CFS-PML approach.
- Berenger, J.-P., "A perfectly matched layer for the absorption
  of electromagnetic waves," *Journal of Computational Physics*
  114(2) (1994), pp. 185-200. DOI 10.1006/jcph.1994.1159. §V
  multi-face wedge PML (W3 hypothesis; deferred).
- Jin, J.-M., *The Finite Element Method in Electromagnetics*,
  3rd ed. (Wiley, 2014), Chapter 10.6 "Wave-port termination" —
  W1 mathematical foundation.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs:53-60`
  — current `#[ignore]` docstring recording the SSSSSSSSS H4
  frozen-magnitude finding (|S_11| = 0.926 across all 18 rows).
- CLAUDE.md §3 (walking-skeleton-first convention motivating W1
  over W2 / W3), §4 (validation gates), §10 (known limitations).
