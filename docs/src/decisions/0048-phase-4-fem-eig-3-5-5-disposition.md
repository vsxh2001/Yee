# ADR-0048 — Phase 4.fem.eig.3.5.5 disposition of fem-eig-006

**Status:** Accepted (queued; v3.5.5 design will pick (a) vs (b))
**Date:** 2026-05-20
**Context Phase:** 4.fem.eig.3.5.5 (post-v3.5.4 multi-mode wave-port
landing, ADR-0047; cutoff-degeneracy finding closed the modal-basis
approach at the current test frequency).

## Context

ADR-0046 §Decision (5) blueprinted v3.5.4's multi-mode wave-port
extension. Track XXXXXXXXX shipped it (commits `ef48b5d` → `0399ab9`
→ `bfdd065` → `5653d05`): `PortDefinition` now carries
`Vec<PortMode>` and the fem-eig-006 +x face is populated with
`[TE_{10} (a_inc=1), TE_{20} (a_inc=0), TE_{01} (a_inc=0)]`.

The measurement, however, **did not move the needle**:

| Configuration                                  | `|S_{11}|(30 GHz)`  |
|------------------------------------------------|---------------------|
| v3.5.2 CFS-PML (best H4 row)                   | 0.926               |
| v3.5.3 W1 single-mode TE_{10}                  | 0.925644            |
| v3.5.4 multi-mode (TE_{10}, TE_{20}, TE_{01})  | 0.925637            |

The cause was a **geometry-convention mis-derivation in the v3.5.4
spec §2.2**: the spec assumed `a = 100 mm, b = 10 mm` predicting
TE_{20} cutoff at 3 GHz and TE_{01} at 15 GHz (both propagating at
30 GHz). The real geometry is `A = 100 mm` (cavity propagation
length along +x), and the port-face cross-section is `B = 10 mm`
(broad wall) × `D = 1 mm` (narrow wall). Under the real geometry:

* `TE_{10}` cutoff = `c / (2 B) = 15.0 GHz` — propagating at 30 GHz ✓
* `TE_{20}` cutoff = `c / B = 30.0 GHz` exactly — `β = 0` at the
  test frequency; multi-mode stiffness block contribution vanishes
  identically.
* `TE_{01}` cutoff = `c / (2 D) = 150.0 GHz` — evanescent at 30
  GHz; carries no propagating modal content.

The multi-mode basis therefore collapses to single-mode at 30 GHz
exactly, regardless of how many additional modes the basis names.

## Decision

Defer the v3.5.5 design to a dedicated brainstorming pass. Two
candidates are on the table, both physically motivated:

### Option (a) — Retune `FEM_EIG_006_F_HZ`

Move the test frequency to a value well above the TE_{20} cutoff
(e.g. `40 GHz`, giving `TE_{20}` a margin of 33% above cutoff and
`β` real with non-trivial magnitude). The multi-mode basis then
carries real propagating modal content and the v3.5.4 projection
step has something non-trivial to do. The cavity, mesh, and driver
are otherwise unchanged.

**Risk:** Re-derives the published-benchmark provenance. fem-eig-006
was originally chosen as a stress test of CFS-PML at off-normal
incidence; retuning the frequency changes the physical regime. A
new benchmark target (e.g. published rectangular-cavity Q at the
retuned frequency) is needed to validate the gate value `< 0.1`.

### Option (b) — Absorbing-mode wave-port (Lee-Mittra 1997)

Replace the modal-projection wave-port with an
**evanescent-mode-absorbing** wave-port per Lee-Mittra 1997 §IV.
Rather than projecting onto a finite modal basis, the absorbing-mode
port models the +x face as a half-space with frequency-dependent
absorption matching the local impedance of the propagating modal
content. This handles cutoff-degenerate test frequencies and
unknown modal content alike.

**Risk:** Higher implementation cost than (a) — the absorbing-mode
operator is a new variant of `FaceKind` and a new assembly path,
not a driver-only edit. Lee-Mittra's derivation also assumes a
homogeneous medium on the absorbing face; verify cavity materials
allow it before committing.

## Consequences

* `fem_eig_006_magnitude_bounded` remains `#[ignore]`'d through
  v3.5.5 design and implementation.
* v3.5.4 multi-mode API (`PortMode`, `Vec<PortMode>`,
  `PortDefinition::single_mode`) is permanent — neither (a) nor (b)
  reverts it. Multi-mode wave-port becomes the canonical
  modal-projection shape for any future driver; absorbing-mode (b)
  plugs in alongside as a separate `FaceKind` variant if chosen.
* fem-eig-003 strict band `[-71.53, -55.58] dB` retire from v3.5.2
  is unaffected by both options.
* If v3.5.5 picks (a), the published-benchmark reference for
  fem-eig-006 at the retuned frequency must land before the gate
  retires. CLAUDE.md §4 ("No solver feature ships without a
  published-benchmark validation case") applies.

## References

* ADR-0046 §Decision (5)
  `docs/src/decisions/0046-phase-4-fem-eig-3-5-3-fem-eig-006-retire.md`.
* ADR-0047 — multi-mode wave-port API
  `docs/src/decisions/0047-phase-4-fem-eig-3-5-4-multi-mode-wave-port.md`.
* Phase 4.fem.eig.3.5.4 spec §2.2 (geometry mis-derivation)
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-4-design.md`.
* Lee, S. T. and Mittra, R., 1997, *IEEE Trans. Antennas Propag.*
  45(4), 671–678 — absorbing-mode wave-port termination.
* Pozar, *Microwave Engineering*, 4th ed., §3.3 TE_{mn} cutoff
  frequencies.
