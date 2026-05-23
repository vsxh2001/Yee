# ADR-0049 — Phase 4.fem.eig.3.5.6 queues the absorbing-mode wave-port

**Status:** Accepted (queues v3.5.6; ADR-0048 Option (b) selected)
**Date:** 2026-05-23
**Context Phase:** 4.fem.eig.3.5.5 (frequency-retune escape-hatch).

## Context

ADR-0048 left two candidate dispositions for the `fem-eig-006`
`|S_{11}| < 0.1` gate, which had been pinned at ~0.926 since v3.5.2:

* **(a)** retune `FEM_EIG_006_F_HZ` off the TE_{20} cutoff edge so the
  v3.5.4 multi-mode wave-port basis carries real propagating content;
* **(b)** replace the modal-projection wave-port with an
  evanescent-mode-absorbing port per Lee-Mittra 1997.

Phase 4.fem.eig.3.5.5 executed Option (a): `FEM_EIG_006_F_HZ` was
retuned from 30 GHz to 40 GHz, where TE_{20} propagates with
`β ≈ 554 rad/m` (33% above its `f_c = c / B = 30 GHz` cutoff). The
cavity, mesh, and 3-mode `[TE_{10}, TE_{20}, TE_{01}]` driver were
otherwise unchanged.

### Measurement

| Configuration                                       | `|S_{11}|`        |
|-----------------------------------------------------|-------------------|
| v3.5.3 W1 single-mode TE_{10} (30 GHz)              | 0.925644          |
| v3.5.4 multi-mode (30 GHz, cutoff-degenerate)       | 0.925637          |
| **v3.5.5 multi-mode (40 GHz, TE_{20} propagating)** | **0.955397 (-0.40 dB)** |
| v3.5.5 refinement probe (40 GHz, NY 3→9 NZ 2→6, 5184 tets) | 0.913956 (-0.78 dB) |

The retune did **not** retire the gate. At 40 GHz `|S_{11}|` is
marginally **above** the cutoff-degenerate 30 GHz value, so the
v3.5.4 modal-degeneracy was not the binding constraint either.

### Discretisation excluded

ADR-0048 / spec §4(a) flagged that the native (16, 3, 2) mesh resolves
the transverse cross-section at only ~2.3 cells/λ at 40 GHz, raising the
possibility that the residual was discretisation-limited rather than a
modal-mismatch. A one-shot refinement probe (NY 3→9, NZ 2→6; a 9×
transverse element-count increase to 5184 tets) moved `|S_{11}|` from
0.955 to 0.914 — a ~0.04 shift, nowhere near the 0.1 gate. The residual
is **not** discretisation-limited. The probe was reverted; the native
(16, 3, 2) mesh stands.

## Decision

**Adopt ADR-0048 Option (b): land an absorbing-mode wave-port in Phase
4.fem.eig.3.5.6.** With both modal degeneracy (v3.5.4) and
discretisation (v3.5.5 probe) excluded, the residual `~0.95` is a
genuine limitation of the modal-projection wave-port: projecting onto
a finite TE_{mn} basis cannot fully match the field at the truncation
face of a strongly off-square (100 : 10 : 1) cavity, where the local
field carries content beyond the named propagating modes. Adding more
named modes does not close the gap (v3.5.4 already demonstrated the
collapse-to-single-mode failure mode; v3.5.5 demonstrates that even with
TE_{20} propagating the projection saturates near 0.95).

The Lee-Mittra 1997 §IV absorbing-mode port models the +x face as a
half-space with frequency-dependent absorption matched to the local
impedance of the propagating modal content, rather than projecting onto
a finite basis. This handles unknown / continuum modal content directly
and is the physically correct termination for this fixture.

## Consequences

* `fem_eig_006_magnitude_bounded` remains `#[ignore]`'d through v3.5.6
  design and implementation. The v3.5.5 escape-hatch measurement
  (`|S_{11}|(40 GHz) = 0.955397`) and refinement probe are recorded in
  the test docstring and the `#[ignore]` reason.
* `FEM_EIG_006_F_HZ = 40.0e9` is **retained** — the 40 GHz operating
  point with TE_{20} propagating is the correct regime for an
  absorbing-mode port to be exercised against, and reverting to 30 GHz
  would reintroduce the cutoff degeneracy. v3.5.6 builds on the 40 GHz
  fixture, not the 30 GHz one.
* The v3.5.4 multi-mode API (`PortMode`, `Vec<PortMode>`,
  `PortDefinition::single_mode`) is unaffected and remains the canonical
  modal-projection shape. The absorbing-mode port lands alongside it as
  a new `FaceKind` variant + assembly path (ADR-0048 Option (b) cost
  note), not as a replacement.
* Per ADR-0048's consequences, if v3.5.6 retires the gate the
  published-benchmark provenance for `fem-eig-006` at 40 GHz must be
  established (CLAUDE.md §4). fem-eig-006 is a synthetic
  matched-modal-termination self-consistency fixture (spec §2.3), not an
  external published benchmark like mom-001; the v3.5.6 design must
  state which class of gate the retired value certifies.
* fem-eig-003 strict band `[-71.53, -55.58] dB` retire from v3.5.2 is
  unaffected.

## References

* ADR-0048 `docs/src/decisions/0048-phase-4-fem-eig-3-5-5-disposition.md`
  — the (a)-vs-(b) disposition this ADR resolves.
* ADR-0047 — multi-mode wave-port API
  `docs/src/decisions/0047-phase-4-fem-eig-3-5-4-multi-mode-wave-port.md`.
* Phase 4.fem.eig.3.5.5 spec §2.3, §4
  `docs/superpowers/specs/2026-05-21-phase-4-fem-eig-3-5-5-design.md`.
* Phase 4.fem.eig.3.5.5 plan
  `docs/superpowers/plans/2026-05-21-phase-4-fem-eig-3-5-5.md`.
* Lee, S. T. and Mittra, R., 1997, *IEEE Trans. Antennas Propag.*
  45(4), 671–678 — absorbing-mode wave-port termination.
* Pozar, *Microwave Engineering*, 4th ed., §3.3 TE_{mn} cutoff
  frequencies + field patterns.
