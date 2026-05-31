# The Ideal Filter-Design Web App — Product Vision & Competitor Analysis

**Date:** 2026-05-31 · **Status:** Product north-star (not an implementation spec)
**Author:** autonomous product pass (maintainer brief: "think what the ideal web app
for designing full filters would look like, from choosing the technique to seeing the
parameters at the end; browse other tools for comparison; do the product work")

This is a **steering document**, not a build spec. It exists so that future increments
are chosen against a competitor-grounded north star instead of ad-hoc. It maps the
landscape, states where Yee's studio already wins, and prioritizes the gaps.

---

## 1. The ideal end-to-end flow (technique → final parameters)

A filter-design app should carry one editable intent from "what do you want" to a
manufacturable board, with every stage live-derived from a single source of truth and
honest about confidence. The ideal flow has **seven stages** plus a **dual entry**:

**Dual entry (the Nuhertz pattern):**
- **Guided / novice** — "I want a band-pass at 2.4 GHz, 100 MHz wide, 20 dB rejection
  at 2.7 GHz." The app *recommends* a technique (and says why) and pre-fills the spec.
- **Expert / gallery** — pick the topology directly; full control.

**The seven stages:**
1. **Technique** — choose realization (lumped LC, edge-coupled, hairpin, combline,
   interdigital, stepped-impedance, stubs…), with guidance on what each is good for
   (size, Q, frequency range, manufacturability).
2. **Spec** — response class (LP/HP/BP/BS), approximation (Butterworth / Chebyshev /
   Elliptic-Cauer / Bessel), order or auto-order-from-mask, f0/BW, ripple, return
   loss, source/load Z. A **draggable spec mask** (ADI's signature affordance).
3. **Synthesis** — ideal transfer function vs the spec mask, PASS/FAIL, with the
   prototype g-values / coupling matrix exposed.
4. **Realization** — turn the math into reality: lumped → real **E24/E96 component
   values + BOM**; distributed → **physical microstrip dimensions** + per-resonator
   table + material stackup.
5. **Tolerance / yield** — Monte-Carlo over component & fab tolerances → yield %, the
   single most underserved stage in free tools (only ADI active + Nuhertz commercial
   have it).
6. **Verify (EM)** — a full-wave EM check of the *physical* realization (not just the
   circuit), because distributed filters always need it ("first approximation, verify
   in 3D EM" — every free tool's disclaimer).
7. **Export** — parameter sheet, BOM CSV, and **manufacturable PCB**: Gerber + KiCad
   (+ ideally STEP/3D). This is where every free tool quits.

Cross-cutting: live re-derivation on every edit; an honest confidence chip per stage;
instant tweak/optimize (ADI-style) rather than re-run-from-scratch.

---

## 2. Competitor landscape (browsed 2026-05-31)

Three tiers, each leaving a distinct gap.

### Tier A — Commercial heavyweights
**Ansys Nuhertz FilterSolutions, Keysight Genesys (S/Filter + M/Filter), Cadence AWR
Microwave Office iFilter.**
- Full lumped **and** distributed synthesis; ~200 topological transforms (Genesys);
  pole/zero shaping, duplexers/multiplexers, custom transfer functions.
- **Parameterized physical layout** → seamless hand-off to a 3D EM solver (HFSS /
  AXIEM) with ports/materials/optimization auto-set-up.
- **Dual UI** — Nuhertz "FilterQuick" (novice) vs "Filter Advanced" (expert).
- Export: SPICE netlist, C-code (digital), schematic into the larger EDA flow.
- **Gap they leave:** closed-source, expensive, desktop-only, and the EM step needs a
  *separate* heavyweight license. Not accessible to a hobbyist / student / quick check.

### Tier B — Free web calculators
**Marki Microwave (microstrip + LC tools), rftools.io, rf-tools.com, changpuak
interdigital, RF Wireless World.**
- One-click synthesis from a spec; Marki's microstrip tool is the closest analogue to
  Yee's distributed flow (low-pass stubs/stepped-Z, band-pass interdigital, arbitrary
  Z, substrate pick → exact widths/lengths).
- **Gap they leave (verbatim from Marki's own tool):** *no* Gerber/DXF/PCB export,
  *no* EM verification, *no* BOM or tolerancing, *no* S-parameter performance plots.
  "The design is a first approximation and should be verified in a 3D EM simulation."
  **The flow terminates at dimension output.**

### Tier C — Best-in-class web UX (adjacent domain)
**ADI Analog Filter Wizard** (active op-amp filters, not distributed RF).
- Spec → real-op-amp circuit "in minutes"; draggable response, instant tweak for
  noise/power/voltage, and **live component-tolerance visualization**.
- **Gap they leave:** wrong domain (active analog, not passive RF/microwave) — but the
  *interaction bar* (instant, draggable, tolerance-visualized) is what to match.

---

## 3. Where Yee already wins

Yee's studio (`yee-studio-web`, live at `<owner>.github.io/Yee/studio/`) **already
does what no free tool does and what only Tier-A commercial tools do** — and is the
only one that is open, browser-native (WASM), and end-to-end to a manufacturable PCB:

| Capability | Tier A ($$$) | Tier B (free web) | **Yee** |
|---|---|---|---|
| Lumped LC synthesis | ✅ | ✅ | ✅ |
| Distributed synthesis | ✅ | ✅ (dims only) | ✅ (edge-coupled) |
| S-param vs spec mask + PASS/FAIL | ✅ | ❌ | ✅ |
| E24/E96 component selection + **BOM** | ✅ | ❌ | ✅ (lumped) |
| **Monte-Carlo tolerance / yield** | ✅ | ❌ | ✅ (lumped) |
| Real dimensioned **board + stackup** | ✅ | partial | ✅ |
| **Gerber / KiCad export** | via EDA | ❌ | ✅ |
| **Own** full-wave EM (no extra license) | ❌ (needs HFSS) | ❌ | ⏳ (FDTD; native) |
| Open-source, browser, zero-install | ❌ | ✅ | ✅ |

**Yee's unique position: the only open, browser-native, *end-to-end-to-PCB* filter
designer with its own EM engine.** Free tools stop at dimensions; commercial tools cost
money + a separate EM license + a desktop. Yee closes the whole loop for free.

---

## 4. Gaps vs the ideal — prioritized

Ranked by value × dispatchability (and honesty about the known hard wall).

**P1 — Technique breadth (the most *visible* gap).** The Technique gallery shows 2 live
(edge-coupled, lumped) and **4 greyed "Soon"** (hairpin, combline, interdigital,
stepped-impedance). Every one is pure-math synthesis + dimensional layout —
WASM-safe, validatable against a published design equation, and each lights a gallery
card → live. Order of impact:
- **Stepped-impedance low-pass** — adds an entirely new *response class* (LOW-PASS;
  Yee is band-pass-centric today) and is the simplest distributed topology (alternating
  hi-/lo-Z line sections from g-values, Pozar §8.6). Highest capability-per-effort.
- **Interdigital / hairpin band-pass** — broadens band-pass realization; both have
  canonical published references to gate against (Hong & Lancaster).

**P1 — Guided entry (the most *product-distinctive* gap).** A "recommend a technique
from my spec" wizard (f0/BW/order/response → recommended topology + the *why*). Matches
Nuhertz FilterQuick, differentiates from every free calculator, and is pure decision
logic — WASM-safe and dispatchable. Turns the expert gallery into a dual-UI.

**P2 — Response breadth.** Lumped path is band-pass-centric; add LP/HP/BS so the
guided entry has somewhere to route. Elliptic/Cauer approximation for sharper skirts.

**P2 — Interactivity / optimize.** ADI-style instant tweak + tolerance-effect overlay
(Yee has Monte-Carlo; surface it as a live "drag a tolerance, watch the yield" view).

**P2 — Manufacturing depth.** KiCad STEP/3D export; design-rule hints (min trace/gap
vs the chosen fab class).

**Deferred — full-board EM verify (the honest hard wall).** The Verify stage stays
"Soon" because high-Q microstrip full-board CW S21 in a stable PEC box is
cavity-dominated and matched CPML is unstable into substrate (ADR-0133). The shipped
EM-sim *loads the line* and is cross-validated at the circuit level (F2.3); a true
≥20 dB full-board cross-validation needs a finer-grid / different-measurement research
track (multi-week, deferred). **Do not reopen without a new measurement strategy.**

---

## 5. Recommended next increments (the steer)

1. **Stepped-impedance low-pass topology** — lights the gallery, adds the low-pass
   response class, simplest distributed synthesis, gate-validatable vs Pozar §8.6.
   *Best value-per-effort; recommended first.*
2. **Guided technique-recommender entry** — the product-distinctive dual-UI move;
   pure logic, WASM-safe.
3. **Interdigital / hairpin band-pass** — band-pass breadth, published-reference gates.
4. Then: response breadth (LP/HP/BS), optimize/tune affordance, KiCad STEP.

Each is a normal Yee increment (spec+plan+ADR → worktree+agent → reviewer → gate →
merge), WASM-safe, and validatable — no gate weakening, no faking, no reopening the EM
wall. The EM-verify frontier remains a separately-scoped research track.

---

## 6. Why this matters

Yee is one or two breadth increments away from being **the** obvious free, open,
browser-based answer to "design me a real, manufacturable RF filter end-to-end" — a
niche currently split between dimension-only free calculators and expensive closed
desktop suites. The product work says: **fill the technique gallery and add a guided
entry next**; defer only the EM-verify wall, which is genuinely hard and already
honestly labeled.
