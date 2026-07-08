# Full-Suite Roadmap (FS.*) — from RF design tool to commercial-class EM suite

**Date opened:** 2026-07-07
**Companion to:** `RF-TOOL-ROADMAP.md` (R.* — complete except R.4c/hardware),
`ROADMAP.md`, `FILTER-DESIGN-ROADMAP.md`, `ENGINE-STUDIO-ROADMAP.md`.

Goal, as set by the maintainer: **a full-suite EM solver capable of designing
antennas, filters, and general RF hardware, comparable to commercial products.**
This document records the market research that scopes that goal, the gap
analysis against Yee's shipped state, and the prioritized phase plan.

---

## 1. Market research

### 1.1 Sourced findings (web research, 2026-07-07)

> Method note: a 5-angle research sweep (solver technologies, workflows,
> deliverables, meshing/GPU, pricing/practitioner-minimal-set) gathered and
> extracted the claims below with sources. A partial adversarial verification
> pass (2026-07-08) confirmed the two gprMax GPU claims **3-0 with exact
> source quotes**; the remaining claims are **sourced but not panel-verified**
> (further verification deliberately stopped — cost). They cohere with each
> other and with practitioner folklore.

**The open-source GPU gap is real and unfilled:**
- Mainline openEMS has no upstream GPU path; the only visible attempt is an
  abandoned third-party fork whose last commits (Aug 2022) read "CUDA compiles
  but doesn't work" ([openEMS-CUDA](https://github.com/aWZHY0yQH81uOYvH/openEMS-CUDA)).
- Users report openEMS is memory-bandwidth-bound and under-utilizes high-end
  CPUs — precisely the profile GPU FDTD fixes
  ([openEMS discussion #94](https://github.com/thliebig/openEMS-Project/discussions/94)).
- The one mature open-source GPU FDTD, gprMax, is a ground-penetrating-radar
  code, not an RF design tool; its published numbers — **1194 Mcells/s
  (Kepler), 3405 Mcells/s (Pascal), ~30× over an i7-4790K OpenMP build** — are
  the c. 2019 open-source bar, and it found consumer GeForce cards
  cost-effective because **FDTD is memory-bandwidth-bound**
  ([Comput. Phys. Commun. 2019](https://www.sciencedirect.com/science/article/pii/S0010465518303990)).
  **Both claims panel-VERIFIED 3-0** against the paper's own text ("up to
  1194 Mcells/s and 3405 Mcells/s on NVIDIA Kepler and Pascal architectures
  … up to 30 times faster than the parallelised (OpenMP) CPU solver"; the
  GeForce "cost–performance benefit … is especially notable" while the
  Tesla P100 wins on absolute performance "due to its use of high-bandwidth
  memory").

**Usability, not physics, is the #1 adoption barrier for open EM tools:**
- Manual mesh creation and boundary setup are called the biggest barrier to
  openEMS adoption vs commercial automated meshing (#94); a novice
  "fine-tuning λ/30 vs λ/20 regions" cannot obtain trustworthy results, and
  results are highly sensitive to small manual mesh changes
  ([#130](https://github.com/thliebig/openEMS-Project/discussions/130)).
- Curved 3-D geometry (helical/discone) on the rectilinear grid is an
  acknowledged open gap — the maintainer confirms (Dec 2024) a robust
  curved-surface mapping is still needed; the community's workaround is
  "switch to a NEC MoM code" (#130).
- Naive FDTD cannot serve the RFIC/MMIC segment: no efficient lossy
  thick-metal model, and sub-µm cells blow up runtime (#94).
- The second sweep (2026-07-08) strengthened this axis from both sides.
  Open side: an EM-consulting review states openEMS reaches **numerical
  accuracy comparable to the commercial tools** for antenna problems at zero
  cost, with its main deficit being **less versatile/advanced meshing**
  ([EpsilonForge](https://www.epsilonforge.com/post/commercial-electromagnetic-software/))
  — solver accuracy is not the gap; meshing is. Commercial side: **adaptive
  mesh refinement is branded as the source of HFSS's "gold-standard
  accuracy"** (solve → refine → repeat until converged), it rides a *suite*
  of initial meshers (Classic, TAU, Phi) each specialized per geometry
  class, and the flagship TAU mesher is known to degrade on high-aspect-ratio
  **planar/PCB** geometry where Classic does better
  ([SemiWiki: the HFSS mesh evolution](https://semiwiki.com/eda/306866-a-mesh-by-any-other-name-the-hfss-mesh-evolution/)).
  Push-button convergence is *the* accuracy product commercially — and even
  HFSS finds planar RF stackups a hard meshing case. FS.0 is aimed exactly
  here.

**The licensing wedge:**
- Commercial pricing (dated forum data, order-of-magnitude only): HFSS
  ~$50–70k + per-core add-ons; CST suite ~€100k; academic licenses ~2–3k/yr
  subscription; prices heavily negotiation-dependent
  ([edaboard thread](https://www.edaboard.com/threads/price-of-hfss-cst.166061/)).
- License *terms* push practitioners out: a consultant reports an Ansys
  license contractually barred from consulting use, evaluating openEMS for
  IoT antenna work (#94).
- openEMS's scriptability already enables GDSII-driven flows benchmarked to
  110 GHz (#94) — scriptability is a strength to preserve, not replace.
- Firmer price anchors (2026-07-08 sweep, reseller-published): CST Studio
  Suite HF **perpetual from ~$62,500**, **quarterly lease from ~$3,500**
  (~$14k/yr — the cheapest published commercial entry)
  ([Fidelis/Dassault reseller](https://www.fidelisfea.com/post/how-much-does-cst-studio-suite-cost-and-whats-included));
  and the freemium pattern: **Sonnet Lite** ships the full engine
  capacity-capped at 16 MB solver memory — full accuracy, toy problems only
  ([edaboard](https://www.edaboard.com/threads/free-electromagnetic-simulators-rather-than-commercial-ones.180440/)).

### 1.2 Commercial landscape (domain knowledge — re-verify before quoting)

> 2026-07-08: the second sweep independently sourced (not yet panel-verified)
> the solver positioning of the three flagship rows — HFSS as FEM-first,
> CST differentiating on FIT time-domain broadband solves (hence its
> transient/EMC niche), FEKO as MoM/hybrid for antennas-on-platforms and RCS
> ([EpsilonForge](https://www.epsilonforge.com/post/commercial-electromagnetic-software/)).

| Tool | Solvers | Known for | The workflow moat |
|---|---|---|---|
| **Ansys HFSS** | FEM (freq-domain) + FEM transient, SBR+ asymptotic, eigenmode | The accuracy reference for 3-D structures; **adaptive mesh refinement** (solve → error estimate → refine → repeat until ΔS < tol) is the core of its "push-button accuracy" claim | Optimetrics (sweeps/optimization/yield), wave/lumped ports with de-embedding, HPC licensing |
| **CST Studio** | FIT/TLM time-domain (its heart), FEM, MoM/MLFMM, asymptotic — one UI, many solvers | Broadband time-domain solves of electrically large structures; GPU acceleration of the T-solver is mature | Schematic co-sim, Filter Designer 3D (coupling-matrix assisted tuning), EMC/SAR outputs |
| **Keysight ADS (+ Momentum/RFPro)** | Circuit + harmonic balance + planar MoM + FEM | The RF *board/MMIC designer's* daily driver — circuit-EM co-simulation | Design guides (filter/matching synthesis), foundry PDKs, tuning cockpit |
| **Cadence AWR (MWO + AXIEM/Analyst)** | Circuit + planar MoM + FEM | iFilter synthesis assistant, MMIC flows | Same co-sim moat as ADS |
| **Sonnet** | Shielded planar MoM | The planar accuracy gold standard; small shop pricing | Narrow but deep: planar S-params people trust |
| **Altair FEKO** | MoM + MLFMM + FEM + PO/UTD hybrid | Antennas on platforms (ships/aircraft), RCS, EMC | Solver hybridization for electrically huge problems |
| **Open source** | openEMS (EC-FDTD), MEEP (photonics FDTD), gprMax (GPR FDTD, GPU), scikit-rf (network post-processing, de facto standard), QucsStudio/Qucs-S (circuit) | Scriptable, free | No integrated design→verify→export loop anywhere; no GPU RF FDTD; manual meshing everywhere |

**What practitioners actually use day-to-day** (folklore + the sourced
threads): S11/S21 vs frequency, one antenna pattern cut + gain + efficiency,
Touchstone in/out, a tune loop, and enough meshing automation to trust the
answer. The long tail (SAR, thermal co-sim, RCS, EMC suites) wins enterprise
deals but is not what makes a tool *usable*.

### 1.3 The strategic wedge for Yee

1. **GPU-first FDTD for RF design is an empty niche** — openEMS never got
   there, gprMax serves GPR. Yee already has a certified wgpu backend
   (compute-015/016 parity gates) and a nightly perf harness. Beating the
   3405 Mcells/s Pascal-era bar on modern hardware is table stakes to claim
   leadership; publish benches.
2. **The spec→design→verify→export loop is the differentiator nobody open
   has.** HFSS/CST sell solvers; ADS/AWR sell workflows. Yee's R-track
   (synthesis → EM-in-the-loop BO → byte-checked .s2p/Gerber/JLCPCB) is an
   ADS-style workflow on an open GPU solver — keep leading with it.
3. **Meshing automation is the credibility gate.** Every practitioner thread
   says manual meshing is why open tools stay niche. An FDTD analog of
   adaptive refinement (graded mesh from geometry + solve-refine-resolve
   convergence loop) is the single highest-leverage usability feature.
4. **License pain is the market opening** — five-figure seats, consulting
   bans, annual subscriptions. Open + GPU + workflow hits all three.

---

## 2. Gap analysis — Yee today vs the credible-alternative bar

**Already shipped and gated** (see R.*/A.*/F.*/E.* roadmaps): GPU/CPU FDTD
with parity gates; planar MoM with NEC-4-validated dipole; CPML (per-face),
dispersive ADE, dielectric + strip-conductor loss, vias; aperture ports;
directional S-parameter extraction with complex Γ/T; NTFF patterns; filter +
patch antenna synthesis with EM-in-the-loop BO; Touchstone/Gerber/KiCad/
JLCPCB export; Tauri studio with design/verify panels; Python bindings; WS
server; surrogate GP/BO/NSGA-II.

**The gaps, ranked by the research:**

| Gap | Evidence | Phase |
|---|---|---|
| Automatic graded meshing + convergence loop | #1 barrier in every practitioner thread | FS.0 |
| Antenna catalog beyond one patch (quasi-Yagi, arrays, wire) | this repo's own R.5c assessment | FS.1 |
| Commercial-grade far-field outputs (gain dBi, efficiency, 3-D pattern export) | day-to-day practitioner set | FS.2 |
| Layout **import** (Gerber/DXF in, not just out) | GDSII flows are how real users arrive | FS.3 |
| Multilayer stackups | every real board; MoM side has it, FDTD flow doesn't | FS.4 |
| Yield/tolerance + space-mapping optimization | commercial Optimetrics parity; surrogate crate is ready | FS.5 |
| Network algebra / circuit co-sim (cascade .sNp, matching synthesis) | the ADS/AWR moat; scikit-rf proves the open appetite | FS.6 |
| Published GPU performance leadership | the 3405 Mcells/s bar | FS.7 |
| MMIC support | **deferred** — needs thick-metal loss + sub-µm meshing (#94 says naive FDTD can't); revisit after FS.0/FS.4 | FS.8 |

Explicitly out of scope: parabolic reflectors / electrically-huge asymptotic
solving (FEKO's PO/UTD class — different solver physics), SAR/thermal/EMC
suites (enterprise long tail).

---

## 3. Phase plan

Conventions unchanged: walking skeleton first; every solver-adjacent phase
ships behind a machine-checkable gate against a strong reference; ADRs for
decisions; specs+plans in `docs/superpowers/`.

| Phase | Scope (walking skeleton first) | Gate sketch | Status |
|---|---|---|---|
| **FS.0** | **Auto-mesh + convergence**. **FS.0a (walking skeleton) SHIPPED** (ADR-0204): `yee_engine::automesh` — `auto_dx` rulebook (λ/20-in-dielectric, h/3, min_feature/2, clamped) + `converge_two_port` adaptive-pass loop (dx/√2 per pass, everything cell-denominated held constant in metres, **linear** ΔS criterion per HFSS's ΔS convention, unconverged reported honestly). Three measured lessons: dB criteria blow up at deep notches (15.35 dB at a converged notch); constant-physics rescaling is necessary hygiene but wasn't the bug; the single-ratio observable's launch-equality assumption fails (+10…15 dB plane-A inequality from stub-reflection source re-pumping) — the loop measures the launch-normalized double ratio `\|T_dut\|/\|T_ref\|`, `T = fwd_B/fwd_A`. **FS.0b (queued)**: graded/nonuniform kernel (Taflove ch. 11, bit-exact-on-uniform gate) + refine-where-flagged rules — the measured motivation: the residual 0.198 pass-to-pass movement is all in the stub's staircase-limited open-end skirt, and the next uniform pass costs ~2.4 h | `engine-automesh-001` (release, in the blanket yee-engine CI gates step): the S.6 stub-notch board with **no hand-set dx anywhere** — auto_dx seeds 0.533 mm (h/3 binding), notch trajectory 5.100→4.900→4.850 GHz / −31.8→−35.1→−34.2 dB, converged err **3.0 %** (≤ 5 %) at ≥ 20 dB depth, loop verdict asserted at tol 0.20 (measured 0.1978) | **FS.0a SHIPPED**, FS.0b queued — **top priority** |
| **FS.1** | **Antenna catalog**. **FS.1a.0+1 SHIPPED** (ADR-0205): `truncate_ground_at_cell` (exact-edge unit gate voxel_002) + `yee_layout::quasi_yagi` (Kaneda/Deal topology, scaling-rule seeds, FDTD-calibrated dipole ε = 1+0.18(ε_r−1) — the half-space (ε_r+1)/2 measured 29 % high on thin FR-4) + the **lifted stack** `voxelize_microstrip_open`/`AperturePortSpec::k_lo` (measured root cause: the domain floor's PEC face was an image plane no mask truncation removes; no compute-kernel change — AperturePort was already cell-list based). Gate `engine-antenna-005` GREEN: **dip 5.950 GHz / −20.9 dB vs designed 5.8 → 2.6 %**, in the antenna CI job. **FS.1a.2 SHIPPED**: end-fire NTFF gate `engine-antenna-006` GREEN first run — **F/B 12.3 dB** (pinned ≥ 6), main lobe toward the director, minimum over the reflector; the balun verified by radiation physics. **FS.1b SHIPPED** (ADR-0206): `patch_array_2x1` — 2×1 corporate-fed H-plane pair (λg/4 70.7 Ω transformer junction, exact mirror symmetry); gates GREEN first run: S11 **2.450 GHz / −21.1 dB (0.0 %)**, pattern multiplication within **~0.6 dB of AF theory** (θ = 60°: −14.2 vs −13.6 predicted). FS.1c thin-wire subcell (Holland) queued; N×1 tree recursion mechanical | per-topology: closed-form seed + full-wave S11 + pattern gate (the A-track template); thin-wire vs the MoM NEC-4 dipole | **FS.1a + FS.1b COMPLETE**, FS.1c queued |
| **FS.2** | **Far-field products**. **FS.2a SHIPPED** (ADR-0207): `AperturePortSpec::record` → per-step `(v_src, v_term, i)` in `JobResult::port_records` (the port already computes all three; GPU rejects recording ports, R.3 idiom). Measured lesson: account on the **circuit side** — the naive aperture-side v·i read a non-physical 1.596 ratio because β = dt·h/2ε₀A ≈ 14.5 Ω rivals R. Gate `engine-power-001` GREEN: **closure 0.9917** (EMF supply vs two-resistor dissipation), accepted-by-field 51.3 % (textbook matched-source halving). **FS.2b partial**: `farfield::gain_dbi` shipped with the audited normalization chain; gate `engine-scale-001` GREEN — the NTFF absolute scale certified vs the analytic Hertzian dipole (**1.048/1.029 across 3 (dx, f) configs**; lesson: baseband-Gaussian near-DC leakage caused ±40 % scatter → zero-DC `GaussianPulseEz` added). Gate `engine-gain-001` RED with root-cause hypothesis: patch read 22 dBi (differential healthy at 1.77 dB) because the whole-domain substrate slab forces the equivalence box through dielectric — **FS.2b.1 queued**: finite-extent substrate in the voxelizer, then re-measure. FS.2c efficiency + full-sphere export queued | `engine-power-001` (in the blanket engine CI gates step); gain of the validated dipole vs 2.15 dBi; efficiency = 1 lossless sanity; pattern export byte-checked | **FS.2a SHIPPED**, FS.2b/c queued |
| **FS.3** | **Layout import**: Gerber (RS-274X subset) and DXF → `Layout` polygons; round-trip gate with our own writer; then "import → verify → export" studio flow | import(export(L)) ≡ L byte-semantics; an imported reference board measures within tolerance of its native-built twin | queued |
| **FS.4** | **Multilayer stackups**: N-layer `Stackup` in yee-layout, voxelizer builds multi-ε_r z-stacks + buried traces + through/blind vias (protocol already carries per-cell ε and 3-D masks) | stripline Z₀ vs closed form (the buried-line analog of S.5's ε_eff gate); MoM multilayer cross-check | queued |
| **FS.5** | **Optimization maturity**: Monte-Carlo yield over dimension tolerances on the surrogate (cheap — GP already fits the BO history); space mapping (coarse = closed forms, fine = EM — our pair is exactly the textbook setup); expose in studio | yield estimate vs brute-force MC on a closed-form testcase; space-mapping converges in fewer EM solves than direct BO on the R.4 scenario | queued |
| **FS.6** | **Network algebra**: cascade/de-embed .sNp blocks (S↔T conversion), matching-network synthesis (L-section/stub from complex Γ — R.2's output), renormalization; CLI + studio | textbook cascade identities; match synthesized from a measured antenna Γ improves its measured S11 full-wave | queued |
| **FS.7** | **Performance leadership**: publish Mcells/s on the GPU nightly across grid sizes; kernel fusion/occupancy passes; beat the gprMax Pascal bar (3405 Mcells/s) on current hardware and say so with a reproducible bench | yee-bench numbers in CI artifacts; the README claim backed by the nightly | queued (starts when the GPU runner lands) |
| **FS.8** | **MMIC** (deferred): thick-metal/multi-sheet conductor model + sub-µm graded mesh + GDSII import; only after FS.0 + FS.4 + R.0b follow-ons prove out | vs foundry-published line data | deferred |

**Sequencing rationale**: FS.0 is first because the research is unambiguous
that meshing automation is the adoption gate, and because graded meshing
multiplies every other phase (finer cells only where gaps/edges need them —
the R.4c coupling-floor problem becomes cheap instead of 8× cells). FS.1/FS.2
ride the existing uniform grid and can proceed in parallel lanes with FS.0's
kernel work. FS.7 needs the user-side GPU runner.

*Last updated: 2026-07-08 (later) — **FS.1a COMPLETE** (ADR-0205: truncated ground voxel_002; quasi_yagi generator with FDTD-calibrated dipole ε; the lifted stack + AperturePortSpec::k_lo after the measured floor-is-a-ground negative result; S11 gate 5.950 GHz/−20.9 dB/2.6 %; pattern gate F/B 12.3 dB — both in the antenna CI job). Earlier same day: FS.0a SHIPPED (ADR-0204: auto_dx rulebook +
convergence loop, gate engine-automesh-001 green with three measured
lessons — linear ΔS criterion, constant-physics rescaling, and the
launch-normalized double-ratio observable). §1.1 upgraded with the partial
verification results (the two gprMax GPU claims panel-verified 3-0; CST
pricing anchors, HFSS AMR-centrality, openEMS meshing-gap claims added as
sourced); further verification deliberately stopped (session budget) — the
remaining claims stay "sourced but not panel-verified".*
