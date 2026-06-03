# ADR-0155: F1.1b coupling-coefficient extraction via the FEM driven-sweep (not FDTD-resonant)

**Status:** Accepted
**Date:** 2026-06-03
**Related:** ADR-0108 (F1.1b.1 — the FDTD coupled-resonator driver whose **resonant-split method
was abandoned**; the propagation ε_eff gate shipped instead), ADR-0153/0154 (the FEM driven-sweep
track + the numerical-eigenmode port this reuses), ADR-0093/0094 (the `extract_coupling` DSP
primitive + the analytic coupled-line k reference), `FILTER-DESIGN-ROADMAP.md`,
[[project-filter-design-final-goal]], [[fem-driven-sweep-s21-viable]].

---

## Context

The filter-design app's pipeline is closed-form end-to-end (synthesis → dimensioning → layout →
manufacturing) and the **EM-verification leg is the gap**: extract the inter-resonator coupling
coefficient `k` from a full-wave simulation of a dimensioned coupled-resonator pair, cross-checked
against the analytic reference, so the app's design engine can be EM-validated. The maintainer
chose to **pivot to the app** after the FEM-EM microstrip-S21 track (ADR-0153/0154) reached its
endpoint, and delegated the sub-item pick; F1.1b (coupling extraction) is the highest-value,
most-dispatchable next step and is not blocked on the eframe-vs-Dioxus frontend fork.

**The back-end is the crux.** F1.1b.1 (ADR-0108) already tried the obvious route — an **FDTD
resonant coupled-resonator split** (`run_coupled_pair`) — and **abandoned it** after a multi-iter
saga: *"no box is simultaneously high-Q and non-confining"* (a small PEC box confines the fringing
fields that set the split; a large PEC box rings as a cavity that swamps the spectrum; an open CPML
box collapses the resonator Q so there are no detectable peaks). What shipped instead was a
**propagation** ε_eff gate (`run_coupled_line_eeff` → even/odd ε_eff + a `k_split` ratio). But
`k_split = (ε_eff_e − ε_eff_o)/(ε_eff_e + ε_eff_o)` is **not** the filter's resonant coupling `k`.

A scoping pass recommended retrying the FDTD-resonant `run_coupled_pair` — but that **re-opens
ADR-0108's documented wall with no fundamentally new strategy**, only geometry-tuning optimism.
This ADR does **not** do that.

Instead it uses the **FEM frequency-domain driven sweep** (ADR-0153/0154, just shipped): it is
**inherently wall-free** — one linear solve per frequency, no time-stepping, no cavity ring-up, no
box-mode swamping. B7 already showed FEM microstrip resonances resolve (a 3-pole bandpass with
**correct peak locations**), and `k` depends only on peak *locations*, which are robust to the FEM
analytic-port absolute-level floor (exactly as B4's ε_eff was robust to the same floor).

## Decision

Extract F1.1b's coupling `k` from a **FEM driven-sweep of a coupled-microstrip-resonator pair**,
reusing the ADR-0154 machinery (`layered_microstrip_filter_mesh` + `TraceRect`,
`microstrip_port_numerical_at`, `OpenBoundarySolver::sweep_matrix` with `with_coupled_whitney(true)`)
and the shipped `yee_filter::extract_coupling` peak-split primitive, validated against the analytic
`yee_layout::coupling_coefficient`.

**A decisive de-risk probe ran this before committing** (spike `spike/fem-coupled-k-probe`,
`933940f`): two identical λ_g/2 FR-4 resonators (W = 1 mm, S = 2 mm, h = 1 mm, ε_r = 4.4, f0 =
2.4 GHz), weakly gap-coupled feeds, 63 k tets, swept 2.10–2.70 GHz. Result — **GO**:

- **Two cleanly resolvable peaks**: f_lo = 2.230 GHz, f_hi = 2.340 GHz, with a **−61.7 dB valley
  19.4 dB below the shallower peak** (unambiguous separation, not one smeared bump).
- **k_fem = 0.0481** vs **k_eps = 0.0581** (the like-for-like even/odd ε_eff-split reference) =
  **17.2 %**; vs **k_imp = 0.0646** (`coupling_coefficient`, the impedance-based synthesis-side
  reference) = 25.6 % — both inside the loose ≲30 % de-risk band.
- The low absolute |S21| (~−35 dB, weak coupling × the analytic-port floor) does **not** matter:
  `k` is a peak-*location* ratio. Both peaks sit slightly low (mesh dispersion + feed-gap loading)
  but the *ratio* is barely affected.

Decompose into **ordered N-bricks** (lettered to continue the FEM lineage), each machine-checkable:

| # | Brick | Gate (machine-checkable) | Risk | Deps |
|---|-------|--------------------------|------|------|
| **K1** | Productionize the probe into a `yee-fem` coupled-resonator-k API + a `fem-coupling-001` gate | builds clean (clippy `-D warnings`, fmt); `#[ignore]`'d + `--release`; **two resolvable peaks** (valley a real margin below both) **AND** `|k_fem − k_ref| / k_ref ≤ 0.30` vs `coupling_coefficient` (the synthesis-side k), with the ε_eff-split k also reported | **eng** (de-risked: probe = GO, 17–26 %) | ADR-0154 |
| K2 | Sweep k vs gap S (monotonicity) — a small S-sweep showing k decreases with S, tracking the analytic curve | `k_fem(S)` monotonic-decreasing across ≥3 gaps + each within the K1 tolerance of `coupling_coefficient(S)` | eng | K1 |
| K3 (later) | Qe extraction (external Q) via the FEM driven-sweep group-delay / singly-loaded resonance | a `fem-qe-001` gate vs a synthesis-side Qe reference | research-open | K1 |

**Critical path** = K1 (the k gate — the F1.1b walking skeleton). K2 strengthens it (a curve, not a
point — harder to pass by coincidence). K3 (Qe) is a later increment.

**The gate is non-circular.** The reference is the Kirschning-Jansen quasi-static closed-form
(`coupling_coefficient`/`coupled_microstrip`); the FEM is a full-wave Maxwell solve on the meshed
geometry — the analytic model does not set up the FEM mesh. The probe was a *measurement* (asserted
only pipeline non-degeneracy); K1 turns it into a real PASS/FAIL with a **two-peaks-resolvable**
tripwire + the k-tolerance, so a future regression that smears the peaks or floors k cannot pass.

## Consequences

**K1 achievable in ~1-2 days** (promote the validated probe + a tripwire gate), **high** confidence
(the physics is de-risked). The gate tolerance is a **loose walking-skeleton 30 %** vs the
synthesis-side `coupling_coefficient` — the probe measured 17 % vs the like-for-like ε_eff-split and
26 % vs the impedance k, and there is a **known systematic** (both FEM peaks pulled low by finite
mesh dispersion + the feed-gap capacitive load) that a finer mesh / a Qe-de-embedded peak would
tighten. K1 reports BOTH reference comparisons for traceability; tightening the tolerance is a K2+
lever, **not** something to force now.

**This validates the FEM driven-sweep as the filter app's EM-verification back-end** for coupling —
wall-free, reusing the entire ADR-0154 investment. The FDTD path stays at its ADR-0108 endpoint
(propagation ε_eff) and is **not** reopened.

**Resource discipline (standing):** every heavy cargo run through the bounded box
(`scripts/yee-box.sh`, ≤14 g / 3 cpu); the coupling gate is `#[ignore]`'d + a dedicated `--release`
CI job (the `mom-001`/`fem-eigen` pattern); never the debug workspace test. **Honesty:** the
reviewer enforces `gate_is_real` (the two-peaks tripwire + the k-tolerance are real measured
quantities, not a tautology); a NO-GO geometry would be a printed measurement, not a forced pass.

**Not in scope / do NOT reopen:** the FDTD resonant-split (ADR-0108, abandoned), the FDTD cavity
wall (ADR-0133), the planar-MoM port (ADR-0064), `fem-eig-006`, the FEM strict-in-mask filter
follow-on (the maintainer-deferred ADR-0154 N3 continuation). The Layout→FEM-mesh app integration
(turning a filter `Layout` into the coupled-pair mesh automatically) is a later F1.2.x step — K1/K2
validate the *physics* on a directly-built geometry first.

---

## Update (2026-06-03) — K2 finding: grade vs the ε_eff-split, not the impedance-k

**K2 (k-vs-gap monotonicity, `fem-coupling-002`) confirmed the FEM coupling-`k`
is physically sound.** Across S = 1.5 / 2.0 / 3.0 mm (W = 1 mm, h = 1 mm,
ε_r = 4.4, f0 = 2.4 GHz; only the gap moves) `k_fem` is **strictly
monotone-decreasing — 0.0611 > 0.0481 > 0.0321** — tracking the analytic
fall-off, and **two transmission peaks resolve at every gap** (valley 14–21 dB
below the shallower peak, all far past the 6 dB re-smearing tripwire). The
extraction is reliable over the whole gap range, which is exactly what K2 was
built to stress (a curve is far harder to fluke than the K1 single point).

**But K2 exposed that the two analytic "coupling coefficients" diverge at strong
coupling.** The impedance-`k` `coupling_coefficient` (`k_imp = (Z0e−Z0o)/(Z0e+Z0o)`)
and the ε_eff-split (`k_eps = (f_e²−f_o²)/(f_e²+f_o²)`) are *not* the same number
as the gap closes: at S = 1.5 mm `k_imp/k_eps = 1.375`, **outside the `[1.0,1.3]`
"comparable" band the src already encodes** (`analytic_k_references_finite_positive_and_agree`).
The FEM measures a **resonant split**, so grading it against `k_imp` is
apples-to-oranges at tight gaps:

| S (mm) | k_fem | vs k_eps (like-for-like) | vs k_imp (impedance) |
|--------|-------|--------------------------|----------------------|
| 1.5 | 0.0611 | **10.5 %** | 34.9 % (> 30 % gate) |
| 2.0 | 0.0481 | 17.2 % | 25.6 % |
| 3.0 | 0.0321 | 21.7 % | 7.4 % |

The S = 1.5 mm row is the **best** fit vs the like-for-like `k_eps` (10.5 %) yet
*fails* a 30 % gate measured against `k_imp` (34.9 %) — not because `k` is floored
or smeared (its valley is the deepest of the three, −63.7 dB) but purely because
`k_imp` has drifted away from the resonant split it is being compared to.

**Resolution (maintainer-endorsed).** Both `fem-coupling-001` (K1) and
`fem-coupling-002` (K2) now grade `k_fem` against the **like-for-like ε_eff-split
`k_eps`** (≤ 30 %, the **same** tolerance — a *reference correction*, not a
weakening), reporting `k_imp` for traceability with the strong-coupling
divergence noted inline. `k_eps` is also a Kirschning-Jansen closed-form (the
even/odd ε_eff of `coupled_microstrip`), so the gate stays non-circular (KJ
closed-form vs full-wave FEM). Both gates are green: K1 17.2 % vs `k_eps`; K2
10.5 / 17.2 / 21.7 % across the three gaps.

**Design-loop implication (for F1.2.1).** The synthesis side
(`dimension_edge_coupled`) targets the *impedance*-`k`, but the EM realizes the
*resonant-split* `k` — and the two diverge at strong coupling. The F1.2.1 EM-in-
the-loop (BO) refinement should therefore target the **resonant-`k`** (or
explicitly account for the `k_imp`↔`k_eps` divergence at tight gaps) rather than
assume the impedance-`k` is what the simulator returns. This is the canonical
reference choice now — **do not reopen it**.

---

## References
- De-risk probe: `spike/fem-coupled-k-probe` (`933940f`),
  `crates/yee-fem/tests/coupled_k_probe.rs` (the K1 seed).
- Reused machinery: `yee_fem::{layered_microstrip_filter_mesh, TraceRect, microstrip_port_numerical_at}`,
  `OpenBoundarySolver::sweep_matrix` (ADR-0154); `yee_filter::extract_coupling`
  (`crates/yee-filter/src/extract.rs:57`); `yee_layout::{coupled_microstrip, coupling_coefficient}`
  (`crates/yee-layout/src/coupled.rs:144/207`).
- The abandoned FDTD path: ADR-0108 Update (resonant → propagation).
- Spec: `docs/superpowers/specs/2026-06-03-f1-1b-fem-coupled-resonator-k-design.md`;
  plan: `docs/superpowers/plans/2026-06-03-f1-1b-fem-coupled-resonator-k.md`.
