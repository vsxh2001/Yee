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
| **FS.0** | **Auto-mesh + convergence**. **FS.0a (walking skeleton) SHIPPED** (ADR-0204): `yee_engine::automesh` — `auto_dx` rulebook (λ/20-in-dielectric, h/3, min_feature/2, clamped) + `converge_two_port` adaptive-pass loop (dx/√2 per pass, everything cell-denominated held constant in metres, **linear** ΔS criterion per HFSS's ΔS convention, unconverged reported honestly). Three measured lessons: dB criteria blow up at deep notches (15.35 dB at a converged notch); constant-physics rescaling is necessary hygiene but wasn't the bug; the single-ratio observable's launch-equality assumption fails (+10…15 dB plane-A inequality from stub-reflection source re-pumping) — the loop measures the launch-normalized double ratio `\|T_dut\|/\|T_ref\|`, `T = fwd_B/fwd_A`. **FS.0b.0 SHIPPED** (ADR-0208, worktree-agent track): `JobSpec.spacings` + one CPU kernel dividing by per-axis primal/dual spacing arrays (uniform fill → every divisor bit-equal, so bit-exactness is by construction); gates `compute-018` (bit-exact-on-uniform: **max Δ = 0 exactly**) and `compute-019` (ratio-1.122/cell taper reflection **−52.68 dB**, pinned −48); graded-inside-CPML/dt/NTFF/dispersive rejections; GPU rejects with Unsupported (Auto→CPU). **FS.0b.1 SHIPPED** (ADR-0210, worktree-agent track): `auto_spacings` rulebook (coarse ceiling **without** the feature/2 term — features refine locally, the graded payoff; fine bands at trace edges/gaps ± guard, growth ≤ 1.3 = the compute-019-certified regime; bit-equal-coarse absorber layers; h/n_sub substrate stack) + `yee_voxel::voxelize_microstrip_graded` (bit-identical to the uniform voxelizer on a uniform grid, gate `voxel-graded-001`); gate **`engine-graded-001`** GREEN: the S.6 stub board on the graded grid reproduces the uniform-converged notch — **4.900 GHz @ −37.2 dB, err 1.03 %, at 0.190× the cells** (1.27 M vs 6.68 M; pinned < 0.25); dedicated CI release-gate job. **FS.0b.2a SHIPPED**: the graded fixture is a library API — `board::two_port_board_jobs_graded` returns the (DUT, reference) pair on one DUT-derived grid (the ADR-0204 same-grid lesson in the API shape); `engine-graded-001` now certifies the builder end-to-end + instant structural gate. **FS.0b.2-GPU SHIPPED** (ADR-0214, worktree-agent track): graded spacings on the wgpu backend — one inverse-spacing buffer (binding 8), kernels multiply per-cell; uniform fill **bit-equal** to the scalar GPU path by construction (gate `compute-020`: 0 differing elements on llvmpipe); taper CPU↔GPU probe parity 4.7e-6, GPU taper reflection **−52.68 dB** = the CPU ADR-0208 figure (gate `compute-021`); engine rejection lifted — `BackendChoice::Gpu`/`Auto` run graded via `GpuFdtd::set_spacings` (NTFF-on-graded / z-taper apertures reject, Auto falls back); hardware certification on the GPU nightly. **FS.0b.2b SHIPPED** (ADR-0216): `GradedMeshOptions.scale` (one refinement knob; coarse+fine move together) + `ConvergencePass.cells` + **`converge_two_port_graded`** — the FS.0a loop on rulebook grids, npml/probe-span rescaled per pass in coarse cells. Gate **`engine-automesh-002`** GREEN (heavy-weekly.yml — ~90 min release): stub board, no hand-set dx — trajectory 4.900→5.100→5.050 GHz (converged err **1.0 %**, depth −37.7 dB), final linear Δ|S| **0.1351** (tol 0.20), every pass **≤ 0.19×** the equivalent-resolution uniform cells (1.27/3.52/9.87 M). Measured band-edge lesson: convergence-criterion bins must stop at ~0.96·f_max — the full-band first run failed at Δ = 0.84 entirely from 5.85–6.0 GHz artifacts (a wandering dip + a non-physical +1.05 dB bin at f_max); root-cause queued (FS.0b.2c, low priority) | `engine-automesh-001` (release, in the blanket yee-engine CI gates step): the S.6 stub-notch board with **no hand-set dx anywhere** — auto_dx seeds 0.533 mm (h/3 binding), notch trajectory 5.100→4.900→4.850 GHz / −31.8→−35.1→−34.2 dB, converged err **3.0 %** (≤ 5 %) at ≥ 20 dB depth, loop verdict asserted at tol 0.20 (measured 0.1978) | **FS.0 COMPLETE** (FS.0a + FS.0b.0/1/2a/2-GPU/2b; FS.0b.2c band-edge root-cause queued, low priority) |
| **FS.1** | **Antenna catalog**. **FS.1a.0+1 SHIPPED** (ADR-0205): `truncate_ground_at_cell` (exact-edge unit gate voxel_002) + `yee_layout::quasi_yagi` (Kaneda/Deal topology, scaling-rule seeds, FDTD-calibrated dipole ε = 1+0.18(ε_r−1) — the half-space (ε_r+1)/2 measured 29 % high on thin FR-4) + the **lifted stack** `voxelize_microstrip_open`/`AperturePortSpec::k_lo` (measured root cause: the domain floor's PEC face was an image plane no mask truncation removes; no compute-kernel change — AperturePort was already cell-list based). Gate `engine-antenna-005` GREEN: **dip 5.950 GHz / −20.9 dB vs designed 5.8 → 2.6 %**, in the antenna CI job. **FS.1a.2 SHIPPED**: end-fire NTFF gate `engine-antenna-006` GREEN first run — **F/B 12.3 dB** (pinned ≥ 6), main lobe toward the director, minimum over the reflector; the balun verified by radiation physics. **FS.1b SHIPPED** (ADR-0206): `patch_array_2x1` — 2×1 corporate-fed H-plane pair (λg/4 70.7 Ω transformer junction, exact mirror symmetry); gates GREEN first run: S11 **2.450 GHz / −21.1 dB (0.0 %)**, pattern multiplication within **~0.6 dB of AF theory** (θ = 60°: −14.2 vs −13.6 predicted). FS.1c thin-wire subcell (Holland) queued; N×1 tree recursion mechanical | per-topology: closed-form seed + full-wave S11 + pattern gate (the A-track template); thin-wire vs the MoM NEC-4 dipole | **FS.1a + FS.1b COMPLETE**, FS.1c queued |
| **FS.2** | **Far-field products**. **FS.2a SHIPPED** (ADR-0207): `AperturePortSpec::record` → per-step `(v_src, v_term, i)` in `JobResult::port_records` (the port already computes all three; GPU rejects recording ports, R.3 idiom). Measured lesson: account on the **circuit side** — the naive aperture-side v·i read a non-physical 1.596 ratio because β = dt·h/2ε₀A ≈ 14.5 Ω rivals R. Gate `engine-power-001` GREEN: **closure 0.9917** (EMF supply vs two-resistor dissipation), accepted-by-field 51.3 % (textbook matched-source halving). **FS.2b SHIPPED**: `farfield::gain_dbi` (audited normalization chain); `engine-scale-001` GREEN — NTFF absolute scale certified vs the analytic Hertzian (**1.048/1.029 across 3 (dx, f) configs**; lesson: baseband near-DC CPML leakage caused ±40 % scatter → zero-DC `GaussianPulseEz`). **FS.2b.1 SHIPPED**: `voxelize_finite_board` (real boards end — gate voxel_003) fixed the measured 22 dBi excess (|F| −16.7 dB, p_acc unchanged — the whole-domain slab had forced the box through dielectric); `engine-gain-001` GREEN: patch **5.42 dBi** (textbook 5–7), array **7.63 dBi**, differential **+2.21 dB**. **FS.2c SHIPPED**: `radiation_efficiency` (gain theorem, quadrature) + byte-stable `pattern_csv`; gate `engine-eff-001` GREEN first run — lossless η **0.806**, tan δ = 0.02 → **0.294** (the 30–60 % FR-4 literature range), lossy p_acc rose (loss broadens the match). **FS.2 COMPLETE** | `engine-power-001` (in the blanket engine CI gates step); gain of the validated dipole vs 2.15 dBi; efficiency = 1 lossless sanity; pattern export byte-checked | **FS.2 COMPLETE** |
| **FS.3** | **Layout import**. **FS.3.0 SHIPPED** (ADR-0209): `yee_export::import` — the RS-274X region subset our writer emits, exact coordinates, named rejections for everything else (inches/polarity/arcs/stroked/flashes). Gate `gerber-rt-001` GREEN first run: vertex-exact + `export∘import∘export` **byte-identical** on 4 real generator layouts. **FS.3.1a+b SHIPPED**: `gerber_to_outline` (gate `gerber-rt-002` corner-exact) + the studio `import_gerber` command with a **byte-provable echo** (gate `studio-import-e2e-001`, in CI); **FS.3.1c SHIPPED**: studio `ImportPanel` (pickers + paste, stackup/port fields, SVG preview, **echo badge** = strict byte equality, vitest DOM gates); **FS.3.2a SHIPPED** (ADR-0217): geometry generality under full-wave test — `yee_layout::double_jog` (four-bend through line; `MiterStyle::Square` vs `Mitered{f}` 45° outer-corner cuts, per-corner polygons so the automesh AABBs refine every bend) + gate `voxel-poly-001` (45°-cut staircase exact; mitered masks fewer cells by the right area) + gate **`engine-miter-001`** GREEN first run (graded fixture, ADR-0216 criterion band, 621 s): mitered band-mean |S21| **0.9738 vs square 0.9665**, mitered better at **every bin** (worst −0.27 dB) — the repo's first non-axis-aligned edges in a measurement; measured lesson: the miter advantage is U-shaped (four-bend interference), so the gate asserts bin-wise dominance, not a frequency trend. **FS.3.2b SHIPPED** (ADR-0220, workflow-agent track): import-side G02/G03 arcs (tessellated at pinned 1 µm chord tolerance, endpoints exact) + D03 C/R flashes with per-D-code %AD bookkeeping; named rejections keep the subset boundary explicit (FlashInRegion, UnknownAperture, UnsupportedAperture, BadArc; G74). Gate `gerber-rt-003` GREEN (pinned 18-seg quarter-arc, sagitta ≤ 1 µm, CW+CCW, rejection matrix); KiCad-style pad/rounded-corner files now import (incl. through the studio panel). Remaining FS.3.2: DXF + the imported-reference-board-vs-native-twin measurement gate | import(export(L)) ≡ L byte-semantics ✓; an imported reference board measures within tolerance of its native-built twin (FS.3.2) | **FS.3.0+3.1 SHIPPED** |
| **FS.4** | **Multilayer stackups**. **FS.4.0 SHIPPED** (ADR-0215): `yee_layout::Stackup` (N layers + lid) + `yee_voxel::voxelize_stackup` (contiguous ε-bands from k = 0 — the ADR-0108 no-gap lesson generalized; buried trace at any interface; single-layer case **bit-identical** to `voxelize_microstrip`, gate `voxel-stackup-001`); gate **`engine-stripline-eeff-001`** GREEN: symmetric stripline vs the **exact** TEM ε_eff = ε_r — **measured 0.065 %** (pinned ≤ 2 %). Three measured lessons: box-mode cutoff must clear the band; gate window must clear the pulse tail; **confined lidded modes need ≥ ~16 cells across b** (at 8 cells β reads 7 % high from discrete transverse-operator error — a future automesh rule, FS.4.2). **FS.4.1 SHIPPED** (ADR-0221, workflow-agent track): `with_via_between` (blind, node-plane to node-plane) + `with_through_via_at_cell` (ground→trace→lid barrel; `with_via_at_cell` now a pure delegation) — structural gate `voxel-stackup-002` pins the exact masked-cell set (7+3, neighbours clear); full-wave gate **`engine-stackup-via-001`** GREEN (~5.5 min release, bit-reproducible): open λ/4 stripline stub notches **−39.81 dB @ 5.075 GHz** (1.5 % off design via the b·ln2/π open-end correction), the through-via-shorted stub erases it band-wide (min −1.18 dB). Measured lesson: the ADR-0215 box-mode rule dissolves under CPML side walls; b ≥ 16 cells carries over. FS.4.2+ queued: stripline Z₀, per-layer tan δ, MoM cross-check, automesh stackup integration | stripline ε_eff vs exact TEM ✓ (`engine-stripline-eeff-001`); stripline Z₀ vs closed form; MoM multilayer cross-check | **FS.4.0 SHIPPED** |
| **FS.5** | **Optimization maturity**. **FS.5a SHIPPED** (ADR-0211): `yee_surrogate::yield_mc` — model-agnostic `yield_estimate` (the pass closure wraps a closed form, a GP, or the engine), deterministic in-crate splitmix64 + Box-Muller (no `rand` dep; same seed ⇒ bit-identical, gate `yield-mc-002`), Wilson 95 % CI (non-collapsing at yield → 1). Gates GREEN first run: `yield-mc-001` MC brackets Φ(z) at n = 1e5; **`surrogate-yield-001`** (the roadmap gate) patch-resonance closed form, FR-4 L ± 0.1 mm / ε_r ± 0.05, spec ±40 MHz — **brute-force 0.9721 vs GP-surrogate 0.9720 (Δ = 1e-4)**, both on the analytic 2Φ(2.2)−1 ≈ 0.972. **FS.5b.0 SHIPPED** (ADR-0213): `yee_surrogate::spacemap` — Broyden ASM + Gauss–Newton parameter extraction, fully deterministic; gate `surrogate-sm-001` GREEN first run: patch two-mode HJ-warp testcase, **ASM 0.00143 % spec error in 4 fine evals vs direct BO 44.8 % at the same 5-eval budget** (~31 000×; asserts pinned ≤0.1 % / ≥5×). **FS.5b.1 SHIPPED** (ADR-0218): ASM with the **engine as the fine model** — gate `sm-em-001` GREEN: the measured stub notch driven onto an off-design 5.3 GHz target in **2 fine evals** (seed 0.988 % → **0.121 %**, 531 s). The increment's real finding: **`GradedMeshOptions.snap_edges`** — un-snapped rasterization quantizes geometry to fine-cell multiples (three lengths spanning 34 µm read the identical notch; Broyden oscillated inside one ~2.5 % step), and snapping the nearest node onto every trace edge makes the response continuous; it also halved the apparent coarse-model bias (2.22 % → 0.99 % at the same seed — half was rasterization error). **FS.5c SHIPPED** (ADR-0222, workflow-agent track): studio `YieldPanel` + `yield_estimate` Tauri command over `yee_surrogate::yield_estimate` (ADR-0211 patch-resonance testcase, deterministic seeded MC + Wilson CI) — reproduces the `surrogate-yield-001` brute-force number **exactly** (0.9721 at the gate seed); vitest gate `studio-yield-dom-001` (the suite's first invoke-mocking test). FS.5b.2 (multi-knob R.4 BPF scenario) + FS.5c.1 (space-mapping panel) queued | yield estimate vs brute-force MC on a closed-form testcase ✓ (`surrogate-yield-001`); space-mapping converges in fewer EM solves than direct BO on the R.4 scenario | **FS.5a SHIPPED** |
| **FS.6** | **Network algebra**. **FS.6.0 SHIPPED** (ADR-0212): `yee_io::network` — 2-port S↔T (chain convention `T_cas = T_A·T_B`, derived in module docs), `cascade`, `deembed_left`, strict `cascade_files` (identical grid + z₀, named rejections, no silent resample); singular conversions (`s21 = 0`) are `Error::Network`, never NaN. Gate `net-001` GREEN first run: S↔T 1e-15, thru identity, 3+3 = 6.000 dB with phases summed, associativity, de-embed recovers DUT, **series impedances cascade to their sum** (ABCD identity in S-form, non-zero reflections), File-level happy + 3 rejections. **FS.6.1 SHIPPED**: `renormalize` (Kurokawa reduced to Möbius `S′ = (S−rI)(I−rS)⁻¹` for identical real-z₀ ports), `deembed_right`, `renormalize_file` (the strict `cascade_files` stays strict — explicit renormalize-then-cascade); gate `net-002` GREEN first run (bit-exact identity, 75 Ω closed form 1e-13, round-trip 1e-14, File-layer unblock). **FS.6.2a SHIPPED** (ADR-0219): `yee_layout::single_stub_match` — the Smith-chart shunt-open-stub construction in closed form; gate `stub-match-001` GREEN first run (Pozar Ex 5.2 position d = 0.1104 λ; 96-load machine contract: synthesized pair nulls combined reflection < 1e-9). **FS.6.2b SHIPPED**: gate **`match-em-001`** GREEN — the measured edge-fed-patch Γ (−4.50 dB) matched to **−11.00 dB** by the synthesized stub (**6.49 dB** improvement, bar ≥ 6; antenna CI job). Three-run lesson pinned in ADR-0219: near resonant radiators, a free-β standing-wave fit locks onto the bulk velocity and corrupts Γ — synthesize from a patch-far plane with `fit_standing_wave_known_beta` (the new constrained split; residual flags contamination); FS.6.3 CLI/studio queued | textbook cascade identities ✓ (`net-001`/`net-002`); match synthesized from a measured antenna Γ improves its measured S11 full-wave | **FS.6.0+6.1 SHIPPED** |
| **FS.7** | **Performance leadership**: publish Mcells/s on the GPU nightly across grid sizes; kernel fusion/occupancy passes; beat the gprMax Pascal bar (3405 Mcells/s) on current hardware and say so with a reproducible bench. **FS.7.0 SHIPPED** (ADR-0223, workflow-agent track): `GpuFdtd::sync()` (public blocking device-wait, the benchmark seam) + `crates/yee-compute/examples/bench.rs` (`cargo run -p yee-compute --release --example bench`; vacuum sweep 64³–224³, CPU vs GPU Mcells/s + GB/s + readback ms, 64³ rel-L2 sanity, `--json`). Two measured negative results: workgroup-shape tuning (`(64,1,1)`/`(32,2,2)`/`(8,8,4)` all lost to the existing `(4,4,4)`, root-caused to `gid.x` mapping to the shader's slowest-varying axis) and `STEPS_PER_SUBMIT` chunking (flat across a 32× sweep) — both reverted, kept as documented negatives. The 192³+ throughput dip is root-caused (not fixed) as a memory-roofline effect (working set outgrowing cache reuse), not chunking/bind-group overhead/the >128 MiB binding/thermal. **Verdict: bar NOT MET on this RTX 5060 Ti** — peak 2864 Mcells/s @ 96³ (84.1 % of 3405), while already at 91.9 % of the card's 448 GB/s bus at that grid size; the reference card has 1.63× more bandwidth. FS.7.1 candidates (not started): E/H pass fusion, a uniform-coefficient fast path, `gid.x → k` kernel index remap. **FS.7.1 SHIPPED (ADR-0224, main-lane track) — BAR MET, reversing the FS.7.0 NO-GO:** both queued levers landed, bit-exact gates unmodified throughout. `gid.x → k` index remap (all 9 volume kernels; adjacent threads → adjacent k-fastest memory) erased the 128³→224³ decline outright (flat ~2980–3170 Mc/s vs the old 2410→934 fall), then a from-scratch post-remap workgroup-shape re-tune picked `(32,2,2)` (FS.7.0's shape table was invalidated by the remap — old flat-x 4–5× penalty was a `gid.x`-maps-to-slowest-axis artifact, now gone). E/H pass fusion 6→2 dispatches (`update_h`/`update_e`, union-extent `cell_dims()`, per-component WGSL bodies kept byte-identical, no reassociation) measured positive at every grid size (+29–54 %), no revert needed. **Definitive final sweep: peak 12665 Mcells/s @ 96³ = 372.0 % of the 3405 bar; every grid size in the 64³–224³ sweep clears the bar, worst case (128³) at 130.5 %.** GB/s discussion: the bench's fixed 144 B/cell/step traffic model now reads above the 448 GB/s nameplate at every grid size (not just the small/cache-resident ones as after Task 1 alone) — fusion's cache reuse means the model's streaming-traffic assumption no longer matches what the fused kernel actually moves; an honest overstatement, not a hardware claim. FS.7.2 (shared-memory tiling to close the remaining 96³→128³ cliff) is optional upside, not required — the spec's own conditional ("only if remap+fusion both land and the bar is within ~10 %") is moot since the bar is cleared outright. **FS.7.2 wrap SHIPPED**: `gpu-nightly.yml` now runs the bench after the existing GPU test steps and uploads `bench-results.json`/`bench-results.txt` as the `gpu-bench-results` artifact (90-day retention) on every enabled nightly run; README gained a Performance section publishing the ADR-0224 definitive medians with the reproduce command and the ADR-0223/0224 links | yee-bench numbers in CI artifacts (`cargo run -p yee-compute --release --example bench`); the README claim backed by the nightly | **FS.7.0+FS.7.1+FS.7.2 SHIPPED — bar MET, honest GO recorded and published (ADR-0223, ADR-0224)** |
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
