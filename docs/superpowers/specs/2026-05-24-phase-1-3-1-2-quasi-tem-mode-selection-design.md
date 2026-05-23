# Phase 1.3.1.2 — quasi-TEM mode selection (open/shielded microstrip wave-ports)

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.2 (the open-microstrip quasi-TEM capability the
cross-section eigensolver's wave-port purpose ultimately needs).
**Depends on:** the complete cross-section eigensolver (steps 4→5.8,
ADRs 0050–0058) + the mom-002 experiment scope finding (ADR-0059).
**Blocks:** numerical microstrip wave-ports (the mom-002 numerical-port
path; true quasi-TEM ports).

## 1. Goal

Extend `solve_dense_mixed`'s dominant-mode selection so it lands the
**quasi-TEM mode** of a (shielded) microstrip cross-section — the mode
with **`k_c² ≈ 0`** (no low-frequency cutoff; propagates to DC),
`β² ≈ ε_eff k₀²`. The mom-002 experiment (ADR-0059) found the current
selector rejects `k_c² ≈ 0` modes wholesale (the spurious/gradient
floor), so it cannot find the microstrip dominant mode at all
("no propagating cutoff candidate"). Close that with a quasi-TEM-aware
selection, validated against the Hammerstad-Jensen `ε_eff` for a
canonical microstrip.

This is a NEW capability (open quasi-TEM), distinct from the
closed/slab-loaded cutoff-bearing modes steps 4→5.8 validated. **Not** a
re-grind of the ε_r=10.2 mode-family question.

## 2. The challenge — quasi-TEM vs the gradient null, both at k_c² ≈ 0

The cutoff pencil `A x = k_c² B x` has the curl-free **gradient null
space** at `k_c² ≈ 0` (the spurious floor exists to reject it). The
genuine **quasi-TEM** mode is ALSO at `k_c² ≈ 0` — but it is
**transverse-energy-dominated and propagating** (`β² = (k₀²−k_c²)⟨ε_r⟩ >
0`, large), whereas the gradient nulls are `E_t ≈ 0` (curl-free). So
they are **indistinguishable by `k_c²` alone** but **separable by the
field** — exactly the converged-eigenvector transverse screen step 5.6
established as reliable.

The hard part is *gathering* the quasi-TEM candidate (it's buried in the
gradient cluster at `k_c²≈0`, so a naive shift-invert near 0 is dominated
by the nulls) — analogous to step-5.7's σ-ladder for the β-direct pencil,
applied here to surfacing the near-zero transverse mode.

## 3. Approach (feasibility-first — see §5 escape-hatch)

`crates/yee-mom/src/eigensolver/solve.rs`:

1. **Gather near-zero candidates too.** Relax the candidate-gathering so
   `k_c² ≈ 0` modes are included (not floored out), e.g. a shift-invert
   rung at/just-above `k_c² = 0` (and/or seed from a TEM-like uniform-E_t
   vector), surfacing the quasi-TEM mode alongside the gradient cluster.
2. **Discriminate by field (the existing screen).** Shift-invert each
   candidate's β-direct RQ, screen the **converged** β-direct eigenvector
   for transverse-dominance — the gradient nulls (E_t≈0) fail, the
   quasi-TEM (E_t-dominant) passes — and keep the highest-β² survivor
   (the quasi-TEM dominant). This reuses the step-5.6 machinery.
3. **Preserve the closed-guide path exactly.** A cutoff-bearing guide
   (WR-90, FR-4 slab-loaded) must still select the same mode (its
   dominant `k_c² > 0` mode is unaffected by also-gathering near-zero
   candidates that then lose the highest-β² or transverse comparison).
   The existing gates (WR-90, FR-4, homogeneous, uniform) are the guard —
   bit-identical or unchanged.

## 4. Validation / DoD

- DoD-1. **Quasi-TEM selection lands the microstrip mode:** a (shielded)
  canonical microstrip cross-section (strip width / substrate height /
  ε_r with a published HJ `ε_eff`) solves to a quasi-TEM mode with
  `k_c² ≈ 0`, transverse-dominated, `β > 0`.
- DoD-2. **HJ validation (PUBLISHED benchmark, loose):** the quasi-TEM
  `ε_eff = (β/k₀)²` matches the Hammerstad-Jensen `ε_eff` for that
  canonical microstrip within a loose tolerance (≤5–10%, per the
  placeholder-tolerance policy — the box truncation perturbs the open
  HJ value; a larger box tightens it). A new
  `crates/yee-mom/tests/eigensolver_microstrip_quasi_tem.rs`.
- DoD-3. **No regression to the closed-guide path:** WR-90 TE10, FR-4
  §4 gate (1.39%), homogeneous canary, uniform anchor, vertical-slab,
  coupling guards — all bit-identical / unchanged.
- DoD-4. No new `Cargo.toml` dependency (faer + the existing sparse /
  selection machinery). Lint clean. `reference.rs` untouched (HJ is a
  closed-form formula, add it as a test helper or a small reference fn).

## 5. Risks + escape-hatch (feasibility-first)

(a) **Separating quasi-TEM from the gradient cluster at k_c²≈0** may be
delicate (both at ~0; the near-zero shift-invert is null-dominated).
Mitigation: the transverse-energy screen on the converged eigenvector is
the discriminator (step-5.6-proven); the gathering uses a near-zero
rung / TEM seed. **Feasibility-first:** the FIRST milestone is "can the
selection surface a transverse-dominated quasi-TEM mode on a microstrip
at all?" — if not within a bounded budget (~40 min), document the
specific blocker (what the near-zero gathering needs) + STOP, queue a
follow-on. Do NOT grind into a multi-step chase.
(b) **Closed-guide regression** from relaxing the floor. Mitigation: the
WR-90/FR-4/homogeneous gates are the non-negotiable guard; if a unified
relaxation regresses them, scope the quasi-TEM path as a separate
entry-point (a `solve_quasi_tem`-style method or a flag), keeping the
closed-guide selection bit-identical.
(c) **HJ-vs-shielded mismatch** (box perturbs ε_eff). Mitigation: loose
tol + a large box; document the box size; HJ is the canonical open-line
reference (the defensible benchmark for microstrip ε_eff).

## 6. References

* Hammerstad & Jensen, "Accurate Models for Microstrip Computer-Aided
  Design", IEEE MTT-S 1980 (the ε_eff / Z₀ closed forms — already used
  in yee-design's Balanis calculator + the mom-002 HJ comparison).
* Pozar §3.8 (microstrip quasi-TEM). Jin §8 (quasi-TEM cross-section).
* `crates/yee-mom/src/eigensolver/solve.rs` (selection + the step-5.6/5.7
  machinery), the mom-002 experiment diagnostic (ADR-0059), ADRs
  0056/0057/0058.
