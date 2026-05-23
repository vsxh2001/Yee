# ADR-0060 — Phase 1.3.1.2: quasi-TEM mode selection for microstrip wave-ports

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** 1.3.1.2 (open-microstrip quasi-TEM capability)

## Context

The mom-002 numerical-wave-port experiment (ADR-0059) found that the
cross-section eigensolver — though built for "quasi-TEM microstrip
wave-ports" — currently serves only **cutoff-bearing** closed /
slab-loaded guides: its dominant-mode selection takes the largest-β²
mode with cutoff `k_c²` above the spurious floor, and a **microstrip
quasi-TEM mode has `k_c² ≈ 0`** (no cutoff), so it is rejected as a
curl-free spurious mode. The §4 FR-4 case validated a *closed*
slab-loaded rectangular guide, not an open microstrip. mom-002's port
residual is therefore a missing **capability**, not a wiring bug.

## Decision

Add a **quasi-TEM-aware selection path** in `solve_dense_mixed`: gather
the `k_c² ≈ 0` candidates too (not floored out), and discriminate the
genuine quasi-TEM mode (transverse-energy-dominated, `β² > 0`) from the
gradient null (`E_t ≈ 0`) by the **converged-eigenvector transverse
screen** (the step-5.6-proven discriminator), keeping the highest-β²
survivor. Validate against the **Hammerstad-Jensen `ε_eff`** for a
canonical (shielded) microstrip — the published microstrip benchmark
that has been absent (the §4 cases were closed guides). Scope it
**feasibility-first** with a bounded cap: surface the quasi-TEM mode, or
document the gathering blocker and stop — no multi-step chase.

## Rationale

(1) **On the PRIMARY mission.** The planar-MoM beachhead needs accurate
microstrip wave-ports; this is the genuinely-needed capability the
experiment surfaced (the real microstrip port + the mom-002
numerical-port path), distinct from the closed-guide work already done.

(2) **Reuses proven machinery + a clean discriminator.** The
quasi-TEM-vs-gradient-null separation is by FIELD (transverse energy),
which step 5.6 established is reliable on the converged eigenvector; the
near-zero gathering mirrors step-5.7's σ-ladder idea. So it is
well-targeted, not greenfield.

(3) **Feasibility-first bounds the risk.** Separating the quasi-TEM mode
from the gradient cluster (both at `k_c²≈0`) has real numerical depth; a
bounded feasibility cap + "document the blocker = a deliverable" prevents
an ε_r=10.2-style multi-step chase, per the standing
bounded-experiment / don't-grind rule.

## Consequences

* `solve_dense_mixed` gains a quasi-TEM selection path; the closed-guide
  selection stays bit-identical (the WR-90/FR-4/homogeneous gates guard
  it; if a unified relaxation regresses them, the quasi-TEM path becomes
  a separate entry-point).
* A new microstrip-vs-HJ validation (the missing open-line benchmark).
* If feasible + validated, unblocks the mom-002 numerical port (a
  follow-on adoption track) + true quasi-TEM wave-ports.
* If the near-zero separation proves hard, a documented blocker finding
  scopes the follow-on; no grind.

## References

* Hammerstad & Jensen 1980 (microstrip ε_eff / Z₀). Pozar §3.8.
* ADR-0059 (the experiment scope finding). ADR-0056/0057/0058 (the
  selection + sparse machinery this builds on).
* `crates/yee-mom/src/eigensolver/solve.rs`.
