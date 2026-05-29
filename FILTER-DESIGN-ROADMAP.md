# Yee — End-to-End RF Filter Design Roadmap

**Status:** Proposed (strategic — the project's stated *final goal*)
**Date:** 2026-05-29
**Owner:** (TBD)
**Companion to:** `ROADMAP.md` (the engine roadmap). This file is the
*application* roadmap that sits on top of it.

---

## 0. Vision

> A designer states a filter specification; Yee guides them — stage by stage,
> with the designer approving each step — from that spec to a **manufacturable
> RF filter**: synthesized prototype → coupling matrix → physical layout →
> full-wave-verified S-parameters → fabrication files (KiCad/Gerber for planar
> & lumped; STEP/mechanical CAD for waveguide).

Today Yee is a strong EM **analysis + optimization back-end**: given a geometry
it returns S-parameters (MoM / FDTD / FEM), and `yee-surrogate` (GP + Bayesian
optimization + NSGA-II + active learning) can tune parameters. The filter
*front-end* — synthesis, dimensional mapping, parametric layout, manufacturing
export, and the interactive design flow that ties them together — does not yet
exist. This roadmap builds it.

**Scope decisions (locked 2026-05-29):**
- **Technologies:** planar (microstrip/stripline), waveguide/cavity, and
  lumped-element LC — a technology-agnostic core with three back-ends.
- **Automation:** *synthesis-assisted interactive* — the tool proposes each
  stage; the designer inspects, tunes, and approves before the next stage.
- **Manufacturing output:** KiCad/Gerber (planar + lumped PCB) **and**
  STEP/mechanical CAD (waveguide). (Not GDSII in this roadmap.)
- **Home:** new crates inside the Yee monorepo (`yee-synth`, `yee-filter`,
  `yee-layout`, `yee-export`), reusing the existing workspace + CI + the
  multi-track orchestration pattern (`ROADMAP.md` §5 / CLAUDE.md §5).

---

## 1. The end-to-end pipeline (technology-agnostic)

Ten stages. Each row maps to an existing Yee capability (**HAVE**) or new work
(**NEW**). The interactive flow gates between stages: the tool produces a
stage artifact, the designer reviews/edits in the GUI, then proceeds.

| # | Stage | What it does | Status |
|---|-------|--------------|--------|
| 0 | **Spec capture** | Designer enters response type, f₀, FBW, ripple/return-loss, rejection mask, Z₀, technology + substrate/band/E-series constraints → typed `FilterSpec`. | NEW (`yee-filter`); can reuse the Phase-3 NL design surface to parse natural-language specs. |
| 1 | **Approximation / prototype** | `FilterSpec` → order N + lowpass-prototype g-values (Butterworth, Chebyshev, elliptic, Bessel). | NEW (`yee-synth`). Closed-form (Matthaei-Young-Jones; Hong-Lancaster). |
| 2 | **Network synthesis** | Prototype → bandpass transform → **coupling matrix M** (coupled-resonator) / J-K inverters (distributed) / LC ladder (lumped). Cross-coupled & elliptic via Cameron synthesis. | NEW (`yee-synth`). |
| 3 | **Circuit realization** | M / inverters / ladder → ideal network; compute its S-parameters = the **target response**. | NEW (`yee-filter`). ABCD/S cascade. |
| 4 | **Dimensional synthesis** | Map the abstract circuit to **physical dimensions** for the technology (resonator lengths, coupling gaps, iris widths, tap points, E-series parts). Initial values from closed-form curves, then refine by **parameter/coupling-matrix extraction + surrogate-BO with the EM solver in the loop**. | NEW driver (`yee-layout`) + **HAVE** `yee-surrogate` (BO) + **HAVE** EM engines. |
| 5 | **Layout generation** | Physical dimensions → full parametric 2-D/3-D geometry, meshable by the EM engine and exportable. | NEW (`yee-layout`). |
| 6 | **Full-wave verification** | EM-simulate the *complete* layout → S-params → compare to the spec mask; loop 4–6 if off. | **HAVE** MoM/FDTD/FEM + Touchstone + plotting; NEW: spec-mask overlay in GUI. |
| 7 | **Tolerance / yield** | Monte-Carlo over etch/εr/machining tolerances (cheap via the surrogate) → yield, sensitivity, tuning suggestions. | **HAVE** `yee-surrogate` (cheap MC) + NEW analysis. |
| 8 | **Manufacturing export** | Final geometry → **KiCad/Gerber** (planar/lumped) or **STEP/mechanical** (waveguide). | NEW (`yee-export`). |
| 9 | **Design report** | Bundle spec, synthesis, dimensions, EM-verified response, yield, and fab files into a report. | NEW (small) in `yee-cli`/`yee-gui`. |

**Cross-cutting — the flow orchestrator:** a stage-gated state machine + a
persisted `FilterProject` document, driven by a **GUI wizard** (the
synthesis-assisted interactive experience). NEW (`yee-filter` engine +
`yee-gui` wizard + `yee-cli filter` + `yee-py` design API).

---

## 2. New crates

```
crates/
  yee-synth/    — NEW. Pure-math synthesis: approximation (g-values), bandpass
                  transform, coupling-matrix & inverter synthesis, Cameron
                  cross-coupled synthesis. No EM, no I/O. Heavily unit-testable
                  against published tables. (depends on: yee-core, nalgebra)
  yee-filter/   — NEW. Filter-domain types (FilterSpec, Prototype, CouplingMatrix,
                  Topology, FilterProject), circuit realization (ABCD/S cascade,
                  ideal response), and the stage-gated design-flow orchestrator.
                  (depends on: yee-core, yee-synth, yee-io)
  yee-layout/   — NEW. Technology-specific parametric geometry generators
                  (planar / waveguide / lumped) + the dimensional-synthesis
                  driver (coupling extraction + surrogate-BO + EM-in-the-loop).
                  (depends on: yee-core, yee-mesh, yee-mom, yee-fdtd, yee-fem,
                  yee-surrogate)
  yee-export/   — NEW. Manufacturing writers: KiCad S-expr + Gerber (RS-274X)
                  for planar/lumped; STEP / B-rep for waveguide mechanical.
                  (depends on: yee-core, yee-layout)
```

Extended existing crates: `yee-gui` (design wizard, spec-mask plot), `yee-cli`
(`yee filter` subcommand), `yee-py` (Python design API), `yee-io` (spec-mask,
report), `yee-surrogate` (filter-aware acquisition if needed).

---

## 3. New external dependencies / tools (to evaluate)

| Need | Candidate | Notes / decision |
|------|-----------|------------------|
| STEP / B-rep CAD export (waveguide) | **`truck`** (pure-Rust CAD kernel) vs **OpenCASCADE** via `bindgen` FFI | Prefer pure-Rust `truck` to honor `#![forbid(unsafe_code)]` default; fall back to OCCT FFI (precedent: Gmsh FFI in `yee-mesh`) if `truck` STEP-write is immature. **ADR required.** |
| KiCad output | in-house S-expression writer (no dep) | KiCad `.kicad_pcb`/footprint is documented S-expr; write directly. |
| Gerber output | in-house RS-274X writer (no dep) | Gerber is a simple aperture/D-code text format. |
| Synthesis linear algebra | **`nalgebra`** (already in tree) | Coupling-matrix eigen-decomposition, similarity transforms. |
| Circuit cascade | in-house (small) | 2-port ABCD ↔ S; N-port for junctions. |

No new heavy runtime dependency is expected except the CAD kernel for STEP.

---

## 4. Hard technical realities that shape the plan (READ FIRST)

These come from the engine's current state and are non-negotiable constraints:

1. **The planar EM-analysis engine must be FDTD, not MoM — for now.** The MoM
   microstrip *port* is **proven ill-posed for planar MoM** (ADR-0064 / CLAUDE.md
   §10): the microstrip quasi-TEM mode's dominant field is substrate-normal
   `E_z`, orthogonal to the in-plane RWG basis, so MoM microstrip S-params are
   port-limited (loose tolerance only). FDTD handles microstrip excitation
   correctly and is well-validated (cpml/ntff/dispersive/lumped/skin-depth gates
   all pass). **Planar filter verification (Stage 6) and dimensional synthesis
   (Stage 4) run on FDTD.** A principled MoM microstrip port (aperture/frill
   reciprocity, or TL-based Z₀ de-embedding) is a separate multi-week engine
   track; if it lands, MoM becomes a fast alternative back-end. Do **not** block
   the filter roadmap on it.
2. **The waveguide track gates on FEM wave-port maturation.** FEM multi-port
   S-matrices work for thru-line (fem-eig-004) and 3-port junction
   (fem-eig-005), but the high-aspect wave-port termination (fem-eig-006,
   `|S11|≈0.955`) is **gate-open**, pending the higher-order absorbing-mode
   wave-port (Phase 4.fem.eig.3.5.7, ADR-0070/0049). Waveguide filter S-param
   accuracy depends on closing this. **Phase F3 is sequenced after that engine
   work** (or uses FDTD for waveguide as an interim back-end).
3. **EM-in-the-loop optimization is expensive.** A single FDTD filter solve is
   seconds-to-minutes; dimensional synthesis needs many. `yee-surrogate`
   (GP + BO + active learning) is exactly the tool — build the dimensional
   synthesizer *around* the surrogate, never raw grid search.
4. **Do not reopen the deferred quagmires** (CLAUDE.md §10): MoM microstrip
   port, FDTD subgrid Q6/Q7, FEM real-port beyond the queued increment, full
   Sommerfeld tail. The filter roadmap routes *around* them by engine choice.
5. **Every shipped stage needs a published-benchmark validation gate**
   (CLAUDE.md §4). Synthesis stages validate against textbook g-tables and
   Cameron coupling matrices (exact); end-to-end stages reproduce *published
   filters* (Hong-Lancaster, Matthaei, Cameron, Pozar).

---

## 5. Phased plan (walking-skeleton first)

Phase IDs follow the `ROADMAP.md` convention. Each phase = spec + plan + ADR
(lockstep) before code, dispatched on disjoint lanes, reviewed before merge.

### Phase F0 — Synthesis walking skeleton (spec → ideal response)
*The minimal end-to-end pipe: pure math, no EM, no layout, no new heavy deps.*
- `yee-synth`: Butterworth + Chebyshev lowpass-prototype g-values; lowpass→
  bandpass transform; J/K-inverter & coupling-matrix synthesis (all-pole).
- `yee-filter`: `FilterSpec`, `Prototype`, `CouplingMatrix`, `Topology` types;
  ideal-circuit S-parameters (ABCD/S cascade); `FilterProject` document.
- `yee-cli`: `yee filter synth <spec.toml>` → prototype + coupling matrix +
  ideal S-params (Touchstone) + a spec-mask pass/fail.
- **Gates:** `synth-001` g-values vs Matthaei-Young-Jones Table 4.05-2(a)
  (exact, ≤1e-6); `synth-002` all-pole coupling matrix vs a published example;
  `filt-001` ideal Chebyshev response meets its own ripple/RL/rejection mask.
- **Why first:** establishes the whole data model + the spec→response contract
  with zero EM cost; everything downstream plugs into it.

### Phase F1 — Planar track to first manufacturable filter (FDTD-backed)
*The headline end-to-end demonstration.*
- `yee-layout`: parametric **edge-coupled** and **hairpin** microstrip
  generators (substrate stack, resonators, coupling gaps, tapped feed).
- Dimensional synthesis: closed-form initial gaps/lengths (Hong-Lancaster
  coupling curves) → **coupling-matrix extraction + `yee-surrogate` BO with
  FDTD in the loop** to hit the Stage-2 targets.
- Stage-6 full-wave verify on FDTD; GUI spec-mask overlay (S-params vs mask).
- `yee-export`: **KiCad + Gerber** writer; export round-trip gate.
- GUI wizard MVP: drive Stages 0→8 with stage-gate approvals.
- **Headline gate:** reproduce the **published Swanson 5-pole hairpin BPF**
  (already targeted as `v1-001` in `validation/README.md`) end-to-end —
  spec → synthesis → layout → FDTD S-params within **±1 dB to 4 GHz** of the
  published reference — and emit Gerber that re-imports to matching geometry.

### Phase F2 — Lumped-element LC track
- `yee-synth` ladder element values; `yee-layout` PCB with component pads +
  parasitic-aware placement; **E-series rounding** + re-verification.
- Stage-6 verify via FDTD lumped-RLC port (engine `LumpedRlcPort`, ADR-0017/0080).
- `yee-export` KiCad PCB with footprints + BOM.
- **Gate:** a published lumped Chebyshev/elliptic LC filter vs analytic +
  measured response.

### Phase F3 — Waveguide / cavity track (STEP output)
*Sequenced after FEM wave-port maturation (Phase 4.fem.eig.3.5.7) — or FDTD interim.*
- `yee-synth` iris/aperture coupling design; `yee-layout` 3-D parametric
  iris-coupled rectangular-cavity generator.
- Stage-6 verify via FEM multi-port S-matrix (or FDTD interim).
- `yee-export` **STEP / mechanical CAD** (CAD-kernel ADR from §3).
- **Gate:** a published X-band (WR-90) iris-coupled 4-pole cavity filter vs
  reference S-params; STEP re-imports to matching solid.

### Phase F4 — Advanced synthesis (elliptic / cross-coupled / multiband)
- Cameron general coupling-matrix synthesis (prescribed transmission zeros),
  cross-coupled topologies, dual-band, diplexers/multiplexers.
- **Gates:** Cameron's published coupling matrices (exact); a cross-coupled
  quasi-elliptic filter reproduced end-to-end on the F1 planar back-end.

### Phase F5 — Tolerance, yield & tuning
- Monte-Carlo over manufacturing tolerances via the surrogate (cheap);
  sensitivity ranking; tuning-screw / trim-pad suggestions; re-tune loop.
- **Gate:** yield estimate for the F1 filter matches a Monte-Carlo EM ground
  truth within tolerance; sensitivity ranking matches analytic expectation.

### Phase F6 — Interactive Filter Design Studio (the product)
- Polished `yee-gui` wizard: spec entry, per-stage review/edit/approve,
  live spec-mask, project save/load, one-click report + fab export.
- `yee-py` scripting API for the whole flow; `yee filter` CLI parity.
- **Gate:** a new user designs, verifies, and exports a spec-compliant filter
  end-to-end through the GUI without touching code (recorded walkthrough).

---

## 6. Dependency graph (what unlocks what)

```
F0 (synth core) ──┬─→ F1 (planar, FDTD) ──┬─→ F4 (advanced synth)
                  │                        └─→ F5 (yield) ─→ F6 (studio)
                  ├─→ F2 (lumped, FDTD)
                  └─→ F3 (waveguide, FEM*) ──────────────────┘
        * F3 gated on Phase 4.fem.eig.3.5.7 (FEM wave-port) or FDTD interim
```

F0 is the unconditional prerequisite. F1/F2 can proceed immediately after F0 on
the proven FDTD back-end. F3 waits on (or routes around) the FEM port. F4–F6
build on F1.

---

## 7. Validation gates summary (the §4 contract)

| Gate | Phase | Reference |
|------|-------|-----------|
| `synth-001` g-values | F0 | Matthaei-Young-Jones Tbl 4.05-2(a) |
| `synth-002` coupling matrix (all-pole) | F0 | published example |
| `filt-001` ideal response meets mask | F0 | self-consistent (spec mask) |
| `filt-planar-001` hairpin BPF end-to-end | F1 | Swanson 5-pole (`v1-001`), ±1 dB |
| `export-001` Gerber/KiCad round-trip | F1 | geometry equivalence |
| `filt-lumped-001` LC Chebyshev | F2 | published LC + analytic |
| `filt-wg-001` iris cavity BPF | F3 | published WR-90 4-pole |
| `export-002` STEP round-trip | F3 | solid equivalence |
| `synth-cameron-001` cross-coupled matrix | F4 | Cameron published tables |
| `filt-yield-001` Monte-Carlo yield | F5 | EM ground-truth MC |
| `studio-001` GUI end-to-end walkthrough | F6 | spec-compliant design, no code |

---

## 8. Risks & open questions

- **CAD kernel for STEP** (`truck` pure-Rust vs OCCT FFI) — ADR needed before F3;
  affects the `#![forbid(unsafe_code)]` posture.
- **FDTD throughput for in-loop synthesis** — may need the GPU path (cuSOLVER /
  CUDA) or surrogate-heavy strategies to keep F1 dimensional synthesis tractable.
- **Coupling-matrix extraction robustness** from noisy EM S-params (the
  group-delay / least-squares extraction must be stable enough to drive BO).
- **MoM microstrip port** — if/when a principled port lands it adds a fast
  planar back-end; tracked in the engine `ROADMAP.md`, not here.
- **FEM wave-port** — F3 timeline is coupled to Phase 4.fem.eig.3.5.7.
- **Scope creep** — three technology tracks is broad; F0→F1 (planar) is the
  single most important proof; resist starting F3 before F1 ships.

---

## 9. Immediate next step

Stand up **Phase F0**: write the spec + plan + ADR for `yee-synth` +
`yee-filter` (the synthesis core and data model), then dispatch the build on the
`crates/yee-synth/**` + `crates/yee-filter/**` lane with the `synth-001` /
`filt-001` published-table gates as the DoD. No EM, no new heavy deps — a clean,
fully-validatable walking skeleton that every later phase plugs into.
