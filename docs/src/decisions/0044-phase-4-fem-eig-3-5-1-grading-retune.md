# ADR-0044: Phase 4.fem.eig.3.5.1 scope — CFS-PML grading retune via ablation sweep, retires fem-eig-003 + fem-eig-006 strict gates

## Status

Accepted — 2026-05-20 (spec + plan; implementation queued under the
Phase 4.fem.eig.3.5.1 R1-R5 ladder).

## Context

Phase 4.fem.eig.3.5 (ADR-0043) shipped CFS-PML
(`AbcOrder::CfsPml(PmlConfig)`) on `OpenBoundarySolver` with the
literature-recommended defaults: `κ_max = 5` (Roden-Gedney 2000 Table
I), polynomial grading order `m = 3`, `thickness_cells = 6`, and the
analytic
`σ_max = (m+1) / (150 π h_cell √ε_r)`,
`α_max = 2 π f_centre ε_0` calibrations from Roden-Gedney 2000 §III/IV.

Track OOOOOOOOO's P5 measurement on the fem-eig-003 WR-90 stub
(`(24, 12, 36)` cavity + 6-cell PML shell, ~72 k extended tets,
~140 s `--release` per-frequency point) recorded
`|S_{11}(f)| ∈ [0.281, 0.423]` ⇒ `s11_db ∈ [-11.0, -7.48] dB` across
8-12 GHz. This is ~10 dB better than the Phase 4.fem.eig.3 2nd-order
Engquist-Majda baseline (`[-2.22e-2, -2.86e-5] dB`) — confirming the
PML kernel is functional — but **~30 dB above** the v3.5 spec §6
`[-60, -40] dB` target window. The companion fem-eig-006 high-aspect
(100 mm × 10 mm × 1 mm) fixture's strict `|S_{11}(30 GHz)| < 0.1` gate
measured `0.926`: CFS `α_α > 0` causality canary passes (finite, not
`NaN`), magnitude gate misses by ~9.3×.

Both production strict gates (`fem_eig_003_strict_absorption_floor_gate`,
`fem_eig_003_strict_passive_bound_continuum_limit`,
`fem_eig_006_magnitude_bounded`) consequently remain `#[ignore]`'d at
the v3.5 merge. The ~30 dB miss is not a mesh-density artefact —
NNNNNNNNN's earlier refinement showed ~2 dB/level saturation, well
short of the 30 dB gap. The binding constraint is the **grading
parameter set**.

Three candidate root causes, each independent:

1. **H1 (h_cell heuristic).** The single-`h_cell` resolver in
   `PmlConfig::resolved(freq_hz, h_cell)` mis-predicts the per-axis
   optimal `σ_α_max` on aspect-ratio cells. WR-90: 0.952 mm broad-wall
   cells, 0.847 mm narrow-wall cells, 0.833 mm axial — ~14 % spread.
2. **H2 (`κ_max` calibration).** The Roden-Gedney 2000 Table-I
   `κ_max = 5` is FDTD-calibrated. Berenger 2002 parameter sweeps
   suggest the frequency-domain FEM regime sits at `κ_max ∈ [1.5, 3]`.
3. **H3 (polynomial order).** `m = 3` on a 6-cell shell concentrates
   90 %+ of absorption in the outer 3 cells; Berenger 2002 §IV
   suggests `m = 2` with `thickness = 8` distributes more uniformly.

ADR-0043 §risks "PML grading parameter sensitivity" queued exactly
this retune as Phase 4.fem.eig.3.5.1.

## Decision

Phase 4.fem.eig.3.5.1 runs a 32-configuration ablation grid across the
H1/H2/H3 hypothesis tree against both production fixtures, picks the
winning `(κ_max, m, thickness_cells)` triple plus the per-axis `h_α`
heuristic, and locks it in as the new `PmlConfig::default()`.

Five load-bearing decisions:

1. **Hypothesis tree is depth-first, cheapest knob first.** H1
   (per-axis `h_α`) is essentially free — one resolver-signature
   change, no new sweep cost. H2 (`κ_max`) is six runs. H3
   (`m, thickness`) is nine runs. Ship on first retire; only fall
   through if the prior leaf misses.
2. **Decision criterion per leaf is concrete and shell-checkable.**
   H1 ships if fem-eig-003 worst-case ≤ −25 dB *and* fem-eig-006
   magnitude < 0.5. H2 ships the smallest `κ_max ∈ {1.5, 2, 3}` that
   retires the fem-eig-003 `[-60, -40] dB` band. H3 ships the
   `(m, thickness)` pair retiring *both* fixtures with the smallest
   `m × thickness` product.
3. **Ablation tool lives in `tools/cfs_pml_grading_sweep.rs` as a
   yee-validation example binary.** Not a CI gate — one-off design
   exploration. Runs locally; emits CSV. Only the winning defaults
   land in CI as the un-ignored strict gates.
4. **Per-axis `h_α` replaces single-`h_cell` heuristic
   unconditionally**, even if H1 alone misses both decision criteria.
   The single-`h_cell` heuristic was a known shortcut documented in
   `pml_stretching_lambda`'s "h_cell back-formula" comment; the
   per-axis path is strictly more correct and the unit tests prove
   bit-for-bit equivalence on isotropic meshes (`(4, 4, 4)` cube). The
   `PmlMeshMeta` carrier is the new public type.
5. **Strict gates un-ignore only after retire.** If the §3 decision
   tree exhausts without both-fixture retire (spec §7 risk a),
   per-fixture override values land in
   `docs/src/tutorials/07-fem-open-cavity.md`, the three strict gates
   stay `#[ignore]`'d with refreshed measurement docstrings, and Phase
   4.fem.eig.3.5.2 (CFS `α_α(d)` polynomial grading) is queued. We do
   not weaken any gate tolerance.

CPU-only, single-threaded, FP64 complex. No GPU. No new dependencies.
Same `faer::sparse::FaerLuSolver<Complex64>` surface; per-axis `Λ(ω)`
stays diagonal in the global frame so complex-LDLᵀ is preserved
(ADR-0043 decision (4) carried).

## Consequences

- **`PmlConfig::resolved` signature changes** from
  `(freq_hz: f64, h_cell: f64)` to
  `(freq_hz: f64, mesh_meta: &PmlMeshMeta)`. The single-`h_cell`
  overload is removed. Callers go through
  `OpenBoundarySolver::with_cfs_pml`, which derives `PmlMeshMeta`
  internally; the public signature on `with_cfs_pml` is unchanged. No
  downstream caller-facing break.
- **`PmlConfig::default()` ships with new `(κ_max, m,
  thickness_cells)` defaults selected by the R3 ablation analysis.**
  The new defaults are annotated with the winning sweep-CSV row
  number and the post-retune `|S_{11}|` measurement in a comment
  block immediately above the `impl Default for PmlConfig`. Callers
  that explicitly pass a `PmlConfig` with all fields set are
  bit-for-bit unaffected; callers using `PmlConfig::default()` see
  the new measurement.
- **fem-eig-003 + fem-eig-006 strict gates flip from `#[ignore]`'d
  to CI-default**, without weakening any tolerance. Three gates total:
  `fem_eig_003_strict_absorption_floor_gate`,
  `fem_eig_003_strict_passive_bound_continuum_limit`,
  `fem_eig_006_magnitude_bounded`. Both fixtures become first-class
  CI tests.
- **New `PmlMeshMeta` public type** in `crates/yee-fem/src/lib.rs`.
  Carries per-axis extents and cell counts; computed by
  `OpenBoundarySolver` at builder time. Python users do not see it —
  the Rust builder owns the derivation.
- **`tools/cfs_pml_grading_sweep.rs` lands as a permanent ablation
  artefact** under the `yee-validation` crate's `[[example]]` list.
  Future regrades (Phase 4.fem.eig.3.5.2, future cavity geometries)
  re-use it without re-deriving the loop.
- **`α_α(d)` polynomial grading is the only remaining CFS-PML knob
  not exhausted by v3.5.1.** Deferred to Phase 4.fem.eig.3.5.2 if the
  R3 ablation exhausts without retire (spec §7 risk b) — otherwise
  remains optional/open-ended.

## References

- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune-design.md`
  — design spec.
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-1-grading-retune.md`
  — R1-R5 implementation plan.
- ADR-0043 — Phase 4.fem.eig.3.5 scope (this ADR's parent; §risks
  "PML grading parameter sensitivity" deferral path).
- ADR-0042 — Phase 4.fem.eig.3 scope (grandparent).
- Berenger, J.-P., "Numerical reflection from FDTD-PMLs: a comparison
  of the split PML with the unsplit and CFS PMLs," *IEEE Transactions
  on Antennas and Propagation* 50(3) (March 2002), pp. 258-265,
  DOI 10.1109/8.999615 — the canonical CFS-PML parameter sweep study;
  the empirical basin `κ_max ∈ [1.5, 3], m ∈ {2, 3}` is from figures
  4-7 of this paper.
- Roden, J. A. and Gedney, S. D., "Convolutional PML (CPML)",
  *IEEE MWCL* 10(5) (May 2000) — Table-I defaults the OOOOOOOOO
  baseline inherits.
- Kuzuoglu-Mittra 1996 — CFS modification `α_α > 0`.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  §"OOOOOOOOO P5 status" — `|S_{11}| ∈ [0.281, 0.423]` measurement.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
  §"OOOOOOOOO P5 status" — `|S_{11}(30 GHz)| = 0.926` measurement.
- CLAUDE.md §3, §4, §10.
