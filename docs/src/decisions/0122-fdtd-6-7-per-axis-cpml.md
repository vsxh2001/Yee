# ADR-0122: Phase 2.fdtd.6.7 — per-axis CPML face selection

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0121 (the reactive-port research track — increment 3 needs this),
ADR-0119 (the de-embed bench whose PEC-source calibration is the blocker), the
lumped-LC → PCB goal (maintainer chose "open the reactive-port research track"),
[[project-lumped-lc-and-studio-redesign]]

---

## Context

The reactive-port research track's increment 3 (ADR-0121) needs a **matched-line**
de-embed bench: a parallel-plate guide that **absorbs at both x-ends** (no
source-end echo, no far-wall echo — the artifacts that made the increment-1/2
bench unable to give a long-window + echo-free + clean-anchor reactive
measurement) while keeping **PEC on the transverse (y, z) faces** so the guide
mode is preserved.

The repo's CPML (`crates/yee-fdtd/src/cpml.rs`) is **symmetric on all six faces**
(its own docs: "per-face slabs is a future optimization"). Applied to a guide it
absorbs the transverse PEC walls and destroys the mode — confirmed in 2.fdtd.6.6,
where every attempt to absorb the source end (CPML / σ-sponge / Mur) broke the
bench's `Z₀`/κ calibration. So a matched line is **not buildable** with the
current all-faces CPML. This is the enabling first brick of increment 3.

## Decision

Add **per-axis CPML selection** to `cpml.rs`: a `[bool; 3]` axis mask on
`CpmlParams` (default `[true; 3]` → existing all-axes behaviour, byte-compatible),
honoured in `CpmlState::update_e`/`update_h` and the `pml_depth` lookup so a
disabled axis is skipped (its faces then behave as plain PEC via the existing
`apply_pec`). A caller can then request **x-only** CPML and keep PEC on y/z — the
matched parallel-plate line increment 3 needs.

Validation gate `cpml_per_axis_001` (`#[ignore]`'d, release): on a 1-D-style guide
with **x-only** CPML at both x-ends and PEC transverse walls, an x-travelling
pulse reflects ≥ the existing CPML target (**≥30 dB reduction vs an all-PEC
control**, mirroring `cpml_reflection.rs`), AND the transverse PEC walls are
verified intact (the guide mode survives — a y/z field-symmetry check). The
existing `cpml_reflection` gate (all-faces) must stay green (default `[true;3]`).

## Consequences

**Ships:** a general per-axis CPML capability (useful beyond this track — any
waveguide / matched-line FDTD setup) with its own validation gate. **Unblocks**
increment 3's matched-line bench (the next sub-increment): an x-matched guide on
which the reactive port's `Z_L(ω)` can be measured with a long window AND no
echoes AND a clean anchor — the measurement the 2.fdtd.6.6 bench could not give.

**Gate:** `cpml_per_axis_001` GREEN in CI; `cpml_reflection` (all-faces) and the
FDTD line/coupling gates non-regressed (default axis mask = all-true).

**Not in scope:** the matched-line bench itself (next sub-increment); the
multi-cell aperture port; F2.3 (all ride on the completed increment 3); per-face
(as opposed to per-axis) asymmetric slabs — per-axis is sufficient here.

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-7-per-axis-cpml-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-7-per-axis-cpml.md`.
- `crates/yee-fdtd/src/cpml.rs` (the symmetric all-faces CPML being generalized);
  `crates/yee-fdtd/tests/cpml_reflection.rs` (the ≥30 dB idiom to mirror).
- Roden & Gedney CPML; Taflove & Hagness Ch. 7.
