# ADR-0156: F1.1b Qe extraction deferred — the FEM driven-sweep has a numerical unloaded-Q floor

**Status:** Accepted
**Date:** 2026-06-03
**Related:** ADR-0155 (F1.1b coupling-k via the FEM driven-sweep — K1/K2 shipped; this defers its K3
sibling), ADR-0108 (the abandoned FDTD resonant route), ADR-0097 (F1.2.0 dimensional synthesis,
which already deferred the Qe→I/O-feed dimensioning), `FILTER-DESIGN-ROADMAP.md`,
[[fem-driven-sweep-s21-viable]], [[project-filter-design-final-goal]].

---

## Context

F1.1b extracts the coupled-resonator design parameters from a full-wave sim. The coupling
coefficient **k is done** (ADR-0155 K1/K2: validated, monotonic in gap, graded vs the like-for-like
ε_eff-split). The remaining piece was **Qe** (the external quality factor — how strongly the I/O
feed loads the end resonators). K3 was to extract Qe via the wall-free FEM driven-sweep.

**A de-risk probe ran first** (the maintainer's standing "scope/probe before a research-open
build"; spike `spike/fem-qe-probe`, `8bc5161`): a single λ_g/2 open-open microstrip resonator with
**one** weakly gap-coupled `microstrip_port_numerical_at` feed, `sweep_matrix` → S11(ω), Qe from the
S11 group delay at resonance (`Qe = ω₀·τ_g/2`, the lossless one-port all-pass relation), at two feed
gaps. Measured:

| feed gap | f₀ | τ_g,peak | **\|S11\|@res** | **Qe (naive)** |
|---|---|---|---|---|
| 1.0 mm | 2.291 GHz | 0.982 ns | **0.836** | **7.07** |
| 2.0 mm | 2.299 GHz | 0.142 ns | **0.838** | **1.03** |

Three problems, and a decisive diagnosis:

1. **Qe is implausibly low.** A weakly-coupled microstrip resonator should have Qe in the
   tens-to-hundreds; 1–7 corresponds to near-critical/over-coupling, not the weak coupling the wide
   gaps should give.
2. **The trend is inverted.** Qe *falls* (7.07 → 1.03) as the gap widens; weaker coupling should
   *raise* Qe.
3. **The smoking gun — `|S11|@res` is coupling-invariant** (0.836 vs 0.838 across a 2× gap change).
   ~30 % of the resonant power vanishes *independent of coupling strength*. That is not a
   coupling-dependent absorption (which would track the gap) — it is a **fixed numerical / model
   unloaded-Q floor** in the FEM driven sweep, which caps the measurable Q at single digits.

**Why this differs from k (which worked).** k is a peak-*location* ratio
`(f_hi²−f_lo²)/(f_hi²+f_lo²)` — robust to amplitude/loss, like B4's ε_eff. **Qe is a Q** — it
requires the model to store energy loss-free over ~Qe cycles. The FEM driven sweep's numerical
dissipation (≈30 % per resonance here, coupling-invariant) sets an unloaded-Q floor far below the
Qe filters need, so the naive group-delay extraction reads the floor, not the coupling.

## Decision

**Defer Qe.** Ship the validated coupling-k (ADR-0155 K1/K2) as F1.1b's coupling-extraction
deliverable; do **not** build a K3 Qe gate now. Record the numerical-Q-floor finding so FEM-Qe is
not re-attempted naively.

This is consistent with F1.2.0 (ADR-0097), which **already deferred** the Qe→I/O-feed dimensioning
("do NOT invent a qe→gap formula"). F1.2.1's surrogate-BO EM-in-loop can begin refining on **k**
(the validated, loss-robust quantity); Qe re-enters only when a feed-dimensioning loop genuinely
needs it.

**If/when Qe is revisited, the naive group-delay is NOT the method.** The candidate is a
**Kajfez / S11-circle fit** that explicitly separates the external Qe from the unloaded Q₀ (= the
numerical floor here) — it *might* recover Qe despite the floor, *if* the numerical Q₀ is high
enough relative to the target Qe to resolve (uncertain on the current mesh; a finer mesh / a
lower-dissipation port would raise the floor). An FDTD ring-down in the matched-CPML box is a
distinct back-end with a different loss mechanism. Both are research-open sub-efforts to be scoped
when the need is real — not now.

## Consequences

**F1.1b coupling extraction = k-validated, Qe-deferred.** The app's design engine can EM-validate
coupling (k) but not yet the feed loading (Qe) via FEM. This is an honest scope line, not a hidden
gap: the probe is committed as evidence (non-gating), and the limitation is documented here so the
back-end's capability map is accurate (FEM driven-sweep: GO for k / ε_eff / S21-shape — all
loss-robust phase/location quantities; NOT-GO for Qe / any absolute-Q without a loss-separating fit).

**Not in scope / do NOT re-attempt naively:** the `Qe=ω₀τ/2` group-delay extraction on the current
FEM mesh (numerical-Q-floored). Do not reopen the FDTD resonant-split (ADR-0108), the MoM port
(ADR-0064), or `fem-eig-006`.

---

## References
- De-risk probe (committed, non-gating evidence): `spike/fem-qe-probe` (`8bc5161`),
  `crates/yee-fem/tests/qe_probe.rs` — the one-port resonator + S11 group-delay measurement.
- The shipped coupling-k it complements: ADR-0155, `yee_fem::coupled_resonator_k`,
  `fem-coupling-001/002` gates.
- The prior Qe-dimensioning deferral: ADR-0097 (F1.2.0).
