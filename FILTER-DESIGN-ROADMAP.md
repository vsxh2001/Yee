# Yee — End-to-End RF Filter Design Roadmap

**Status:** Proposed (strategic — the project's stated *final goal*)
**Date:** 2026-05-29
**Owner:** (TBD)
**Companion to:** `ROADMAP.md` (the engine roadmap). This file is the
*application* roadmap that sits on top of it.

---

## 0. Vision

> A designer opens an **app** (desktop **or** in the browser), states a filter
> specification, and is guided — stage by stage, approving each step — from that
> spec to a **manufacturable RF filter**: synthesized prototype → coupling matrix
> → physical layout → full-wave-verified S-parameters → fabrication files
> (KiCad/Gerber for planar & lumped; STEP/mechanical CAD for waveguide).

**The final deliverable is an interactive filter-design application, shipped as
both a desktop app and a web app** (clarified 2026-05-29). The synthesis/design
flow built by Phases F0–F4 is the *engine* of that app; the **App/Studio track**
(§5a) is the product surface that makes it usable end-to-end without the CLI.

Today Yee is a strong EM **analysis + optimization back-end**: given a geometry
it returns S-parameters (MoM / FDTD / FEM), and `yee-surrogate` (GP + Bayesian
optimization + NSGA-II + active learning) can tune parameters. The filter
*front-end* — synthesis, dimensional mapping, parametric layout, manufacturing
export, and the interactive design flow that ties them together — is what F0–F4
build; the app/web-app wraps it.

**Scope decisions (locked 2026-05-29):**
- **Delivery: a desktop + web app.** One Rust/`egui` codebase via `eframe`,
  which targets **native** (desktop) and **WASM** (browser, WebGL/WebGPU). No
  separate JS rewrite — the existing `yee-gui` (egui 0.34 / wgpu 29) is the seed.
  See ADR-0089 for the architecture (light flow client-side in WASM; heavy EM +
  surrogate optimization on a native `yee-server` the web client calls).
- **Technologies:** planar (microstrip/stripline), waveguide/cavity, and
  lumped-element LC — a technology-agnostic core with three back-ends.
- **Automation:** *synthesis-assisted interactive* — the app proposes each
  stage; the designer inspects, tunes, and approves before the next stage.
- **Manufacturing output:** KiCad/Gerber (planar + lumped PCB) **and**
  STEP/mechanical CAD (waveguide). (Not GDSII in this roadmap.)
- **Home:** new crates inside the Yee monorepo (`yee-synth`, `yee-filter`,
  `yee-layout`, `yee-export`, plus `yee-studio` + `yee-server` for the app),
  reusing the existing workspace + CI + the multi-track orchestration pattern.

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

### Phase F0 — Synthesis walking skeleton (spec → ideal response) — **SHIPPED** (ADR-0084, merge `dbfc5c5`)
*The minimal end-to-end pipe: pure math, no EM, no layout, no new heavy deps.*
Shipped 2026-05-29: `yee-synth` + `yee-filter` crates + `yee filter synth` CLI;
gates `synth-001`/`synth-002`/`filt-001` pass as crate tests; closed-form ideal
response used for `filt-001` (coupling-matrix→S realization is F1+). The
`yee-validation` aggregator registration of the three gates is the small
follow-on **Phase F0.1**.
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
*Superseded/expanded by the App/Studio track (§5a) — F6 is the desktop-app
milestone within it. Retained here as the capstone of the F-series flow.*
- Polished filter-design app: spec entry, per-stage review/edit/approve,
  live spec-mask, project save/load, one-click report + fab export.
- `yee-py` scripting API for the whole flow; `yee filter` CLI parity.
- **Gate:** a new user designs, verifies, and exports a spec-compliant filter
  end-to-end through the app without touching code (recorded walkthrough).

---

## 5a. App / Web-app track (the final deliverable)

The F-series above builds the *flow*; this track wraps it in the shipped
**desktop + web app**. One `egui`/`eframe` codebase, two build targets (native +
WASM). Architecture per **ADR-0089**: the *light* flow (spec → synthesis →
coupling matrix → layout preview → spec-mask plot — all pure-Rust, WASM-safe,
already shipped as F0/F0.1/F0.2/F1.0) runs **client-side**; the *heavy* steps
(FDTD/FEM verification, surrogate dimensional synthesis, mesh/export) run on a
native **`yee-server`** the web client calls over HTTP, and in-process for the
desktop app. New crates: `yee-studio` (the egui app) + `yee-server` (axum API).

- **App.0 — `yee-studio` desktop skeleton. ✅ SHIPPED** (ADR-0090, merge `338a35c`).
  A native `eframe` app: spec-editor panel → synthesis panel (g-values, coupling
  matrix, Qe, PASS/FAIL) → `egui_plot` |S21|-vs-spec-mask view, recomputed live.
  `StudioState` (the flow logic) is egui-free + WASM-safe per ADR-0089 (App.1
  reuses it); only `app.rs`/`main.rs` are native. Gate `studio_state_recompute`
  (headless, pass+fail). Layout preview deferred (needs the F1.2 dims mapping).
  TODO(App.1): cfg-gate `mod app` before the WASM build.
- **App.1 — WASM web build of the light flow.** Compile the App.0 light path to
  `wasm32-unknown-unknown` via `eframe` web; deploy as a static site (CI →
  Pages). Everything through the ideal-response spec-mask view runs fully in the
  browser, no server. **Gate:** `trunk build` / `wasm-pack` produces a loadable
  bundle; a headless WASM smoke test (wasm-bindgen-test) exercises the flow.
- **App.2 — `yee-server` EM/optimization backend.** An axum service exposing the
  heavy steps (FDTD/FEM verify, surrogate dimensional synthesis from F1.1+,
  mesh, KiCad/Gerber/STEP export) as JSON/artifact endpoints. The web client
  calls it for the F1.1+ stages; the desktop app links the engine directly.
  **Gate:** a round-trip — web client POSTs a `FilterProject`, server returns an
  EM-verified S-parameter Touchstone + a pass/fail against the spec mask.
- **App.3 — full end-to-end in the app + deploy.** Spec → … → fab-file download,
  in both desktop and browser; project save/load; design report. **Gate:** the
  F6 walkthrough, performed in the deployed web app, end-to-end without code.

**Sequencing:** App.0 can start now (it consumes only shipped light crates) and
proceeds in parallel with the F1.1+ EM work; App.1 follows App.0; App.2 lands
once F1.1–F1.4 give the server something to verify/optimize; App.3 is the
capstone. The light client (App.0/App.1) is **not** blocked on the EM loop.

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

## 9. Status & next step

**SHIPPED so far (2026-05-29):**
- **F0** (ADR-0084, merge `dbfc5c5`): `yee-synth` + `yee-filter` +
  `yee filter synth`; `synth-001`/`synth-002`/`filt-001` green.
- **F0.1** (ADR-0085, merge `e71e400`): the three synthesis gates registered in
  the `yee-validation` aggregator under a new `Solver::Synth` / `yee validate
  synth` target — they now appear in `yee validate --list[ --json]`.
- **F1.0** (ADR-0086, merge `9a51655`): `yee-layout` crate — parametric
  microstrip geometry (edge-coupled + hairpin generators), Hammerstad-Jensen
  width/ε_eff synthesis, dependency-free SVG preview; gates `geo-001/002/003`.
  Geometry-only (no EM yet); consumes explicit dims (the coupling-matrix→dims
  mapping is F1.2).
- **1.plotting.4** (ADR-0087, merge `8d6e81f`): `yee-plotters` spec-mask overlay
  (`draw_sparam_with_mask` + `mask_violations`) for the Stage-6 verification view.
- **F0.2** (ADR-0088, merge `4de1a28`): `yee filter synth --plot` — renders the
  synthesized |S21| with the spec mask overlaid (the spec→visual pipe).
- **Filter-synthesis theory chapter** (merge `fe45016`): `docs/src/theory/filter-synthesis.md`.
- **App.0** (ADR-0090, merge `338a35c`): `yee-studio` eframe desktop app — the
  first product surface; spec → synthesis → spec-mask plot, live.
- **App.1.0** (ADR-0092, merge `ead2819`): `yee-studio` eframe shell gated behind a
  default `desktop` feature so `StudioState` builds egui-free (web-ready).
- **App.1.1** (ADR-0095, merge `d901a2c`): `yee-studio --no-default-features` PROVEN
  to compile to `wasm32-unknown-unknown` (egui absent); CI `wasm-build` job gates it.
- **F1.1a** (ADR-0091, merge `c4f3af4`): `yee-voxel` crate —
  `voxelize_microstrip(&Layout) → YeeGrid` (tangential Ex+Ey PEC masks); gate voxel_001.
- **F1.1b.0** (ADR-0093, merge `be2d2bc`): `yee-filter::extract` — `extract_coupling`
  + `extract_q_ringdown`, validated vs analytic signals (no FDTD).
- **F1.1b.gate** (ADR-0094, merge `e676c42`): `yee-layout::coupled_microstrip`
  (Kirschning-Jansen 1984 even/odd model + coupler k); gates coupled_001 vs Steer
  Ex 5.6.1 (≤0.21%) + coupled_002 monotonic k. The validatable `k` reference for F1.1b.1.
- **App.1.2a** (ADR-0096, merge `92f1696`): the full `yee-studio` eframe `StudioApp`
  now compiles for `wasm32-unknown-unknown` behind a `web` feature + a
  `#[wasm_bindgen(start)]` `eframe::WebRunner` browser entry (eframe split
  per-target; WebRunner takes a DOM `HtmlCanvasElement`; no RUSTFLAGS — wgpu 29
  WebGL2 fallback). Gate `cargo check -p yee-studio --target wasm32 --features web`
  exit 0; native + headless-wasm builds unregressed. (App.1.2b = trunk bundle/deploy.)

**Final goal: a desktop + web APP** (ADR-0089) — one `egui`/`eframe` codebase,
native + WASM. The shipped light flow (F0/F0.1/F0.2/F1.0) is WASM-safe and is the
in-browser front-end; heavy EM goes behind a native `yee-server`. See §5a.

**Two parallel fronts next:**
- *Product:* **App.1.0/1.1 ✅ SHIPPED** (ADR-0092 `ead2819` / ADR-0095 `d901a2c`):
  `yee-studio` eframe shell behind a default `desktop` feature, and the
  `--no-default-features` light path PROVEN to compile to `wasm32-unknown-unknown`
  (egui absent from the dep tree) under a CI `wasm-build` job. App.1.2 split in two:
  **App.1.2a ✅ SHIPPED** (ADR-0096, merge `92f1696`) — the full eframe `StudioApp`
  compiles for wasm32 behind a `web` feature + a `#[wasm_bindgen(start)]`
  `eframe::WebRunner` entry (eframe split per-target with `default-features=false`
  + `["wgpu","default_fonts"]` on wasm; `WebRunner::start` takes a DOM
  `HtmlCanvasElement` fetched via `web-sys`, not a string id — eframe 0.34.2; no
  RUSTFLAGS needed, wgpu 29 WebGL2 fallback compiles clean). Gate
  `cargo check -p yee-studio --target wasm32 --features web` exit 0; native +
  headless-wasm unregressed. **NEXT = App.1.2b** = `trunk` bundle +
  `index.html`/`Trunk.toml` + static deploy (needs `cargo install trunk` — HEAVIER;
  riskiest = wgpu-WebGL2 at runtime; the web entry already expects
  `<canvas id="the_canvas_id">`). Delivers the loadable in-browser app.
- *Engine (toward the Swanson-hairpin FDTD gate):* **F1.1a `yee-voxel` ✅ SHIPPED**
  (ADR-0091, merge `c4f3af4`): `voxelize_microstrip(&Layout) → YeeGrid` (ground
  PEC + substrate ε_r slab + trace PEC, point-in-polygon rasterized; tangential
  Ex+Ey masks per the review P0 fix); gate voxel_001 (no FDTD run). `yee-layout`
  untouched (WASM-safe). **F1.1b.0 `extract` ✅ SHIPPED** (ADR-0093, merge
  `be2d2bc`): `yee-filter::extract_coupling` (k from the two split peaks
  `(f2²−f1²)/(f2²+f1²)`) + `extract_q_ringdown` (Qe = π f0 τ decay-fit), validated
  vs analytic signals — the extraction the FDTD driver feeds. **F1.1b.gate
  coupled-line model ✅ SHIPPED** (ADR-0094, merge `e676c42`): `yee-layout::
  coupled_microstrip` (Kirschning-Jansen 1984 quasi-static even/odd z0e/z0o/εeff +
  coupler `k=(z0e−z0o)/(z0e+z0o)`); gates coupled_001 vs Steer Ex 5.6.1 (≤0.21%) +
  coupled_002 5-gap monotonic k; pure f64, WASM-safe. This is the validatable `k`
  reference for F1.1b.1 and the initial-dimensioning model for F1.2.
  **NEXT = F1.1b.1** — the FDTD coupled-resonator DRIVER: `yee-voxel` voxelize a
  coupled pair → `LumpedRlcPort`s → run `yee-fdtd` → single-bin DFT → `extract_*`,
  gated against the F1.1b.gate even/odd `k`/split-frequency reference. HEAVY
  (multi-min FDTD). Then **F1.2** surrogate-BO dimensional synthesis; **F1.3**
  verify + mask gate; **F1.4** `yee-export`. **App.2** (`yee-server`) once F1.1+ exist.
  (Tutorial 17 — filter design via CLI + Studio — shipped, merge `c6e477c`.)
