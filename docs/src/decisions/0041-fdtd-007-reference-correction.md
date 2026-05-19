# ADR-0041: fdtd-007 Maloney-Smith reference correction (open question)

## Status

Open question — 2026-05-19. Citation in the Phase 2.fdtd.7
subgridding spec is **definitively wrong** for the cited geometry;
the correct reference is **not yet identified**. Surfacing this as an
ADR routes the open question to an `fdtd-007.1` follow-up track and
prevents the wrong citation from being silently re-quoted in
downstream specs, plots, or commit messages.

## Context

Track UUUUUUUU (commit `d56c460`, 2026-05-19) rewired the `fdtd-007`
Maloney-Smith dielectric-loaded thin-slot validation driver onto the
per-cell `ε` map and per-component PEC mask infrastructure
(MMMMMMMM `cb6f8ed`, PPPPPPPP `c57592f`). The rewired uniform-fine
driver measures

```
f_res ≈ 5.30 GHz,   |S_11(f_res)| ≈ −6.4 dB
```

on a 32 × 80 × 25 grid at `dx = 0.5 mm`, 2000 steps — `|df|/f ≈ 0.40`
against the digitised `8.9 GHz` reference encoded in
`yee_validation::FDTD_007_FRES_REF_HZ`, far outside the `±5 %`
digitisation envelope the original Phase 2.fdtd.7 Q7 escape hatch
allowed for.

UUUUUUUU surfaced a finding (verbatim, finding #1):

> The Phase 2.fdtd.7 spec
> (`docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`)
> cites Maloney & Smith 1993 IEEE T-AP 41(5) as the reference for the
> slot, but that paper title — "A study of transient radiation from
> the Wu-King resistive monopole" — is a *cylindrical monopole*
> paper, not a slot. The `FDTD_007_FRES_REF_HZ = 8.9 GHz` value may
> therefore be misattributed; an `fdtd-007.1` follow-up should verify
> the correct reference (or replace the geometry to match
> Maloney-Smith's actual published case).

Track XXXXXXXX verified the finding against IEEE Xplore document
[222286](https://ieeexplore.ieee.org/document/222286). The verbatim
paper record is:

> J. G. Maloney and G. S. Smith, "A study of transient radiation
> from the Wu-King resistive monopole — FDTD analysis and
> experimental measurements", *IEEE Trans. Antennas Propag.*,
> vol. 41, no. 5, pp. 668–676, May 1993.

The abstract describes "a cylindrical monopole antenna with
continuous resistive loading ... using a resistance variation
proposed by Wu and King (1965)". There is no slot, no dielectric
substrate, no `Fig. 9` published `S_11` curve for the geometry the
spec describes. The citation is therefore **wrong for the cited
geometry**, regardless of whether `f_res ≈ 5.3 GHz` or `8.9 GHz` is
the physically correct answer.

## Decision

1. **Treat the existing `FDTD_007_FRES_REF_HZ = 8.9 GHz` and
   `FDTD_007_S11_DB_REF = -22 dB` constants as unverified.** The
   physics gates against these constants stay `#[ignore]`'d per the
   LLLLLLLL escape hatch ("Do NOT relax to > 5 %"). The constants are
   not touched in this ADR — Track XXXXXXXX's lane is documentation
   only; changing the constants is an `fdtd-007.1` deliverable that
   must follow the *resolved* reference.

2. **Mark the spec citation as `[TBD verify]`** with an inline
   pointer to this ADR. Keep the Wu-King reference visible in the
   References block so the next reader can see what was previously
   cited and what is now disputed; do not silently delete it.

3. **Open `fdtd-007.1`** to (a) identify the correct published source
   for the dielectric-loaded thin-slot geometry, or (b) replace the
   geometry to match a Maloney-Smith case that *is* in the published
   record. The plausible candidate identified in the Track XXXXXXXX
   brief — Maloney, Smith, Scott, "Accurate Computation of the
   Radiation from Simple Antennas Using the Finite-Difference
   Time-Domain Method", IEEE T-AP, vol. 38, no. 7, pp. 1059–1068,
   July 1990 — analyses an open-ended parallel-plate waveguide, a
   cylindrical monopole, and a conical monopole. It is **also not a
   slot antenna paper**, so it is not a drop-in fix either.

## Chain of evidence

- **Citation as listed in the spec
  (`docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`,
  pre-this-ADR, §References, line 195):**
  > Maloney, J. G., Smith, G. S., "A study of transient radiation
  > from the Wu-King resistive monopole — FDTD analysis and
  > experimental measurements", *IEEE Trans. Antennas Propag.*
  > 41(5), 1993, pp. 668–676.

- **Actual paper, per the journal record (IEEE Xplore document
  222286, confirmed via web search 2026-05-19):**
  > J. G. Maloney and G. S. Smith, "A study of transient radiation
  > from the Wu-King resistive monopole — FDTD analysis and
  > experimental measurements", *IEEE Trans. Antennas Propag.*,
  > vol. 41, no. 5, pp. 668–676, May 1993.
  >
  > Abstract: cylindrical monopole antenna with continuous resistive
  > loading per Wu and King (1965). No slot, no dielectric substrate.

- **Mismatch:** the spec quotes the journal title correctly but
  attributes a dielectric-loaded thin-slot geometry and a `Fig. 9`
  `S_11` curve that this paper does not contain.

- **Resolved reference:** **open question.** No Maloney-Smith
  publication exhibiting the cited geometry (`w = 0.5 mm`,
  `L = 30 mm` slot, `ε_r = 2.2`, `h = 1.524 mm` substrate,
  delta-gap fed) was located in three web searches (2026-05-19,
  Track XXXXXXXX). The candidate Maloney-Smith-Scott 1990 paper
  (T-AP 38(7), 1059–1068) covers monopoles and waveguides, not
  slots.

## Candidates the follow-up track should check

These are *plausible* sources for a `Fig. 9` `S_11` curve on a
dielectric-loaded thin slot; none is verified.

1. **Original spec author's notes.** The cleanest resolution is the
   spec author resurfacing the source they were reading when they
   wrote `f_res = 8.9 GHz`. The figure number `Fig. 9` and the
   specific dimensions (`L = 30 mm`, `w = 0.5 mm`, `ε_r = 2.2`,
   `h = 1.524 mm`) are precise enough that they almost certainly
   come from one paper, not a synthesis.

2. **Other Maloney/Smith papers, 1990–1995.** A bibliography sweep
   on J. G. Maloney's Google Scholar profile (Georgia Tech) would
   take an hour and is the next step if (1) fails.

3. **Different first authors.** The geometry parameters are common
   enough that the source may be a different research group entirely
   (e.g. Sheen, Ali, Abouzahra, Katehi 1990 IEEE T-MTT 38(7) on
   microstrip patches and CPW). A `±20 %` envelope on
   `f_res = 5.3 GHz` against published slot-antenna catalogues is
   enough to triangulate the source.

4. **Replace the geometry.** If the source cannot be located,
   `fdtd-007.1` replaces the validation case with a geometry that
   does have a published reference — e.g. the Maloney-Smith-Scott
   1990 cylindrical monopole, which has a clean `Z_in(f)` curve and
   is already a Phase-2 candidate gate.

## Consequences

- **Spec citation marked as `[TBD verify]`** prevents the wrong
  citation propagating into downstream specs, commit messages, plot
  captions, and notebook docstrings.
- **Constants left untouched** keeps the physics-gate `#[ignore]`
  state intact; UUUUUUUU's `fdtd_007_uniform_fine_smoke` (passivity-
  only) gate continues to run on default CI.
- **`fdtd-007.1` follow-up gated** by this ADR — the implementing
  agent must resolve the reference (or replace the geometry) and
  update both `FDTD_007_FRES_REF_HZ` / `FDTD_007_S11_DB_REF` and the
  spec citation in lockstep.

## References

- Track UUUUUUUU commit body, `d56c46037449fd4a41eeeb88e7859d317f6b306e`
  (the original finding).
- IEEE Xplore document 222286
  (https://ieeexplore.ieee.org/document/222286) — Wu-King resistive
  monopole paper record, confirming title / abstract / page range.
- IEEE Xplore document 55618
  (https://ieeexplore.ieee.org/document/55618) — Maloney-Smith-Scott
  1990 monopoles/waveguides paper, candidate but not a slot.
- `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`,
  §Validation gate (`fdtd-007`) and §References.
- `crates/yee-validation/src/lib.rs` lines 1648–1690
  (`FDTD_007_*` constants, including the existing TBD flag on
  `FDTD_007_FRES_REF_HZ`).
- `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs`
  (LLLLLLLL gate tests; physics gates `#[ignore]`'d).
