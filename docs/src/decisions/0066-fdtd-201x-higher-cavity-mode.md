# ADR-0066 — fdtd-201.x higher-order cavity mode resonance gate

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** 2.fdtd validation (sibling to fdtd-201, ADR-0062)

## Context

fdtd-201 (ADR-0062) validated the dominant TE₁₀₁ cavity mode. A read-only
scoping of the FDTD validation landscape (verified against source) found
FDTD already has ~11 live gates (CPML reflection, NTFF dipole pattern,
Fresnel point-source, TF/SF quiet-zone/oblique, dispersive, lumped), and
that the genuinely-clean, infra-ready next increment is a **higher-order
cavity mode** gate — every alternative is either already shipped (dipole
pattern), blocked on missing measurement infra the codebase itself flags
(tight Fresnel VSWR, coax/TEM Z₀), needs materially more build (waveguide
cutoff: PEC-mask + CPML + multi-freq dispersion), or is the lowest-novelty
re-run (strict-±0.5% same-mode refinement). The quagmires (subgrid Q6
energy-balance at 75-79% drift; fdtd-007 wrong-reference) are confirmed
and avoided.

## Decision

Add a tests-only **higher-order rectangular-cavity mode** gate (TE₂₀₁
resonance) that clones the fdtd-201 harness (PEC cavity, Gaussian E_y
injection, single-bin-DFT scan, peak-find), using **a ≠ d** to break the
TE₂₀₁/TE₁₀₂ degeneracy so a *named* mode is validated. Assert against
analytic Pozar §6.3 within a documented loose ±2.5% (grid dispersion
floor, worse at higher f), `#[ignore]`-gated. No `src/` change.

## Rationale

(1) **Maximally dispatchable / near-zero risk** — same source, observable,
reference family, PEC-cavity setup, and DFT scanner as a gate that landed
cleanly days ago (fdtd-201). The harness is proven on this exact geometry.

(2) **Genuine new §4 coverage** — validates mode *selectivity* and the
solver's higher-frequency grid-dispersion behaviour, a distinct claim
from "the dominant mode is right" (so it is not a re-run of fdtd-201).

(3) **Value × dispatchability** — a marginal-but-certain clean win, the
right autonomous-loop move now that the high-value tracks are
blocked/quagmire/major (mom-port intrinsic per ADR-0064; FDTD Q6; FEM
real-port) and the cleaner alternatives are shipped or infra-blocked.

(4) **Avoids all quagmires** — touches none of the subgrid / fdtd-007
surface; tests-only, so it cannot regress the solver.

## Consequences

* A new `#[ignore]`-gated higher-mode cavity test (in `cavity_higher_mode.rs`
  or `cavity_resonance.rs`) + a `validation/README.md` row naming the
  validated mode.
* a≠d geometry chosen so the target mode is non-degenerate + cleanly
  separated in the scan band.
* The strict ±0.5% refinement + Q-factor extraction remain documented
  follow-ons. No `src/` change; no new dependency.

## References

* Pozar §6.3 (`f_mnp`). The FDTD-gate scoping (2026-05-24, read-only):
  `cavity_resonance.rs` (the harness), `heterogeneous_substrate.rs:263` +
  `lumped_resistor.rs:9` (why the Fresnel/Z₀ candidates are infra-blocked),
  `subgrid_energy_balance.rs` + `fdtd_007_*` (the quagmires).
* ADR-0062 (fdtd-201, the gate this extends). Spec + plan (2026-05-24).
