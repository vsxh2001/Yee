# Engine + Studio Roadmap (GPU/CPU compute engine → web studio)

Direction set by **ADR-0175** (2026-07-05). This is the third top-level roadmap, alongside
`ROADMAP.md` (core EM solvers, Phases 0–4) and `FILTER-DESIGN-ROADMAP.md` (filter application).
It tracks the two-part re-centering of the project:

- **Part 1 — the engine (E.\*):** a fast Rust simulation engine that uses GPU *and* CPU.
  Portable wgpu/WGSL compute + rayon CPU in the new `crates/yee-compute`, with `yee-fdtd`'s
  scalar kernels kept as the validated reference and `yee-cuda`'s cuSOLVER LU lane unchanged.
- **Part 2 — the studio (S.\*):** an engine-service protocol (`yee-engine` → `yee-server`) and a
  modern web-technology studio (Tauri 2 shell + React/TypeScript frontend) that drives the
  engine in the background, in-process on desktop and over WebSocket in the browser.

Spec: `docs/superpowers/specs/2026-07-05-gpu-engine-web-studio-design.md`
Plan: `docs/superpowers/plans/2026-07-05-gpu-engine-web-studio.md`

Conventions match the other roadmaps: every phase ships behind a machine-checkable validation
gate; walking-skeleton first; phases get ADRs when they make a decision worth recording.

---

## Part 1 — Engine track (E.*)

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| **E.0** | `yee-compute` walking skeleton: `FdtdSpec`/`Fields`/`FdtdEngine`, rayon FP64 `CpuFdtd`, wgpu/WGSL FP32 `GpuFdtd`, uniform lossless vacuum + PEC box | `compute-001` (CPU **bit-exact** vs `yee-fdtd` scalar reference, 25 steps, non-cubic grid); `compute-002` (GPU vs CPU, rel-L2 < 1e-4 / L∞ < 1e-3, 100 steps; self-skips without adapter, real on GPU nightly) | **SHIPPED** (ADR-0175, this branch) |
| **E.1** | CPML + per-cell ε_r/μ_r/σ + interior PEC masks + legacy PEC box on both backends; GPU arena-buffer layout (5 storage bindings — inside WebGPU browser limits) | `compute-003` (CPU **bit-exact** vs reference, heterogeneous + CPML + masks, both boundary modes); `compute-004` (CPML reflection: **69.3 dB** measured vs ≥ 30 dB target); `compute-005` (GPU vs CPU on the full E.1 scenario: ~2e-7 E / ~3e-6 H family-rel on llvmpipe; CPML holds 210× less ‖H‖ than PEC) | **SHIPPED** (ADR-0176) |
| **E.2** | Drive layer: `SoftSource`/`ResistivePort`/`Probe`/`Drive` on both backends (GPU: whole-run f64-precomputed tables + on-GPU step counter → zero per-step host round-trips) | `compute-007` driven step **bit-exact** vs reference; `compute-006` cavity TE₁₀₁ vs **analytic Pozar**: CPU −0.063 %, GPU −0.063 %, CPU↔GPU 0.0000 %; `compute-008` line-eeff on the engine vs **Hammerstad–Jensen**: 0.132 % (≤ 15 % gate), 88.6 s release | **SHIPPED** (ADR-0177) |
| **E.3** | Precision policy: FP32-GPU/FP64-CPU characterized (WGSL has no f64 — SHADER_F64 unreachable without SPIR-V passthrough; noted) | `compute-009` drift over 10⁴ energy-conserving steps: 3e-6…2e-5 family-rel (√N random-walk), 100× inside the 1e-3 gate | **SHIPPED** (ADR-0177) |
| **E.4** | Performance: `yee-bench` `compute_step` (scalar vs rayon CPU vs GPU) landed; container numbers recorded | 4-core container: rayon scales 2.2× internally but nets **0.78×** vs scalar (flat-buffer kernel ~2.8× slower single-thread — bounds-checked idx arithmetic). Row-sliced kernels landed (ADR-0179): single-thread −27 %, 4-thread ≈ scalar (bandwidth-bound container); bit-exact gates unchanged. Real-hardware numbers via the GPU nightly bench; the 20×-dGPU target remains to be certified there | **CLOSED** (ADR-0179; hardware numbers pending nightly) |
| **E.5a** | Far-field on the engine: engine steps, reference `NtffState` consumes fields via host adapter | `compute-010` vs **analytic sin θ**: broadside/endfire 327.9 dB (≥ 20 dB gate) | **SHIPPED** (ADR-0177) |
| **E.5b** | On-GPU full-field DFT phasor accumulation (`accumulate_dft` kernel, psi-arena tail, on-GPU step counter — zero per-step readback); reference `NtffState` projects via two synthetic samples | `compute-013`: GPU-resident dipole — **315.4 dB** analytic null, broadside matches the CPU path to **2.9e-7** | **SHIPPED** (ADR-0179) |
| **E.5c** | Dispersive ADE (Drude/Lorentz/Debye) on both backends: verbatim CPU port; unified-ADE GPU form folded into the coeff/psi arenas | `compute-011` **bit-exact** vs `yee_fdtd::dispersive` (four-arm scenario); `compute-012` differential GPU gate (ADE ≤ 20× standard-pair error, measured ≤ 6×; drift-class backstop) | **SHIPPED** (ADR-0179) |

Non-goals for E.*: replacing `yee-cuda`'s cuSOLVER LU lane (stays as-is); MoM/FEM assembly on
wgpu (revisit after E.4 with data).

## Part 2 — Studio track (S.*)

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| **S.0** | `yee-engine` crate: serde `JobSpec`/`JobEvent`/`JobResult` protocol + threaded chunked executor with progress streaming, cooperative cancel, cpu/gpu/auto backend selection | 4 unit tests + doctest: serde round-trip, progress stream, cancellation, auto-backend | **SHIPPED** (ADR-0179) |
| **S.1** | `yee-server` (axum 0.8): `/healthz` + WS `/v1/jobs` streaming live `JobEvent`s; disconnect cancels via `JobCanceller`; `yee serve` CLI subcommand | end-to-end tokio-tungstenite gates in the workspace suite (round trip incl. probe series + field slice; invalid-spec error event); `/healthz` verified live | **SHIPPED** (ADR-0180) |
| **S.2** | Tauri 2 + React/TS/Vite studio shell (`studio/`, outside the root workspace) speaking S.0 in-process: `run_job` command + `job://progress` events + probe SVG plot. Frontend 47.9 kB gzipped | walking skeleton verified in-container: `cargo check` (webkit2gtk) + `npm run build` green; interactive run + CI wiring are the S.2 follow-on | **SKELETON SHIPPED** (ADR-0179) |
| **S.3** | Visualization walking skeleton: engine `slice` option (final E-plane in `JobResult`) → canvas heatmap + single-bin-DFT spectrum plot in the studio; `studio-build` CI job (typecheck + vite + vitest + Tauri cargo check) | vitest gates: DFT recovers a **known sinusoid** to one bin; color-map extremes; DOM smoke renders of both views (7 tests) | **SKELETON SHIPPED** (ADR-0180) |
| **S.3b** | three.js 3-D field surface: height-mapped vertex-colored mesh + orbit controls, lazy-chunked (initial bundle stays 49.4 kB gz; three rides a 133 kB on-demand chunk); WebGL fallback | geometry is a pure function gated against hand-computable values; the fallback path DOM-renders under jsdom (11 vitest tests total) | **SHIPPED** (ADR-0181) |
| **S.4** | Parity audit done (ADR-0181 capability table): the studios serve disjoint roles — Dioxus = shipped filter designer, Tauri = engine studio. Freeze stands; retirement re-decided when the filter flow consumes engine jobs (the defined convergence path via `yee-server`) | audit table + decision recorded | **AUDITED** (ADR-0181; retirement deferred) |
| **S.5** | Engine-powered verify, walking skeleton: `JobSpec` gains per-cell materials + interior PEC masks (`MaterialsSpec`) and an explicit `dt_s` (both `#[serde(default)]` — the protocol stays backward-compatible), validated at submission (`Error` events, no worker panics); voxelized layouts now run over the S.0/S.1 protocol on both backends | `engine-verify-001` (`#[ignore]`, release CI): the fdtd-line-eeff-001 FR-4 microstrip expressed **as a `JobSpec`** through `submit()`/events → ε_eff vs **Hammerstad–Jensen** ≤ 15 %; fast gates: serde round-trip + legacy specs, heterogeneous job **bit-exact** vs direct `CpuFdtd`, 4 malformed-spec error paths | **SHIPPED** (ADR-0182) |
| **S.6** | S-parameters on the engine, walking skeleton: `yee_engine::sparams` (`single_bin_dft`, `transmission_db` — pure post-processing over `JobResult` probe series) + the two-run reference/DUT transmission method with a passive resistive-port termination (`v0 = 0`); no protocol or `yee-compute` changes | `engine-sparams-001` (`#[ignore]`, release CI): λ/4 open-stub bandstop over the job protocol — notch **4.850 GHz / −36.8 dB**, **3.0 %** from the closed-form TL-theory prediction (±15 % / ≥ 8 dB gate; band-edge standing-wave ripple bounded at |12| dB, measured +8.7/+5.2 dB); fast gates: known-sinusoid DFT, −6.02 dB scaled copy | **SHIPPED** (ADR-0183) |
| **S.7** | |S11| via incident/reflected separation: `sparams::reflection_db` — the reference run's port-1 probe is the incident wave, `dut − ref` isolates the device reflection; zero extra solve cost (one more probe on the same two jobs) | `engine-sparams-001` extended: at the stub notch **\|S11\| = −0.93 dB** (≥ −4 dB gate — a λ/4 open stub reflects ~everything at resonance) and **\|S11\|²+\|S21\|² = 0.807** (physical band [0.5, 1.3]; the deficit is the documented second-order re-reflection + ripple); fast gate: synthetic 0.25 reflection → −12.04 dB | **SHIPPED** (ADR-0184) |
| **S.8** | = **F1.3.0**: the first filter **synthesized by the pipeline and verified by the engine** — N=5 Butterworth stepped-impedance LPF (f_c = 2 GHz, FR-4) from `yee_synth::prototype` → `dimension_stepped_impedance_layout` → voxelize → two engine jobs → measured response vs the `ideal_response_lowpass` design targets; gate lives in `yee-filter/tests` (dev-deps on the engine keep the lib WASM-safe) | `engine-filter-verify-001` (`#[ignore]`, own release CI step): cutoff **1.900 GHz vs designed 2.0 GHz — 5.0 %** (±20 % gate, sustained-crossing scan); **rejection 30.6 dB** (ideal 30.1; ≥ 20 dB gate); passband-mean ripple bound ±6 dB (measured +3.4 dB; PEC-box single-probe ripple documented up to +17.8 dB at the band edge). Boundary finding recorded: all-face CPML collapsed the passband — root-caused in S.9 | **SHIPPED** (ADR-0185) |
| **S.9** | Per-axis CPML on the protocol: `BoundarySpec::Cpml` gains serde-defaulted `axes: [x, y, z]` (pre-S.9 JSON still parses) → `CpmlConfig::with_axes`. **ADR-0185 collapse root-caused**: the ~5-cell substrate sat *inside* the 10-layer z-min absorber — a scenario error, not a CPML defect. Board-level open boundary = `[true, true, false]` (absorbing side walls, PEC ground/lid), adopted by the LPF gate | Re-measured `engine-filter-verify-001` under CPML-xy: passband mean **+1.32 dB** (was +3.42 PEC-box), stopband **−32.9 dB**, rejection **34.2 dB**, cutoff unchanged 1.900 GHz — better on every aggregate; residual ripple attributed to the lumped-port mismatch (the next fidelity lever). Fast gates: legacy-JSON default + axes round-trip | **SHIPPED** (ADR-0186) |

Standing decision during S.*: **`yee-studio-web` (Dioxus) is feature-frozen but stays deployed**
until S.4 concludes (ADR-0175). `yee-gui` (egui EM-analysis shell) is unaffected by this track.

---

*Last updated: 2026-07-06 (latest) — S.9 SHIPPED (ADR-0186): per-axis CPML on the protocol
(`BoundarySpec::Cpml.axes`, serde-defaulted) and the ADR-0185 collapse root-caused — the
~5-cell substrate sat inside the 10-layer z-min absorber. Board-level open boundary =
CPML x/y + PEC ground/lid; the LPF gate re-measured better on every aggregate (passband
mean +1.32 dB, rejection 34.2 dB, cutoff unchanged); residual ripple = lumped-port
mismatch, the next fidelity lever. Before that, S.8/F1.3.0 SHIPPED (ADR-0185): the
design→verify loop closed for the first time — an N=5 Butterworth stepped-impedance LPF
synthesized by the pipeline, run through the engine, measured cutoff 1.900 GHz vs designed
2.0 GHz (5.0 %) with 30.6 dB passband/stopband rejection (ideal 30.1 dB). Follow-ons:
matched terminations / de-embedding, F1.3 spec-mask API, F1.2.1 EM-in-loop refinement.
Before that, S.7 SHIPPED (ADR-0184): |S11| via incident/reflected
separation (`sparams::reflection_db`, zero extra solve cost); at the stub notch
|S11| = −0.93 dB and |S11|²+|S21|² = 0.807 — both halves of a filter response now come out
of two engine jobs. Before that, S.6 SHIPPED (ADR-0183): `yee_engine::sparams`
(single-bin DFT + transmission ratio, pure post-processing) and gate `engine-sparams-001`:
a λ/4 open-stub bandstop run twice over the job protocol (reference line / DUT, passive
resistive-port termination) notches at **4.850 GHz / −36.8 dB — 3.0 %** from the closed-form
TL-theory prediction. The full filter-verify chain now exists on the engine: layout →
voxelize → JobSpec → FDTD → probes → |S21|(f) + |S11|(f). Earlier: S.5 SHIPPED (ADR-0182): per-cell
materials + PEC masks + explicit dt on the job protocol (`MaterialsSpec`/`dt_s`,
serde-defaulted, validated at submission with Error events instead of worker panics); gate
`engine-verify-001` runs the fdtd-line-eeff-001 FR-4 microstrip **as a JobSpec** through
`submit()`/events and recovers ε_eff to **0.132 %** of Hammerstad–Jensen (identical to
compute-008 — the protocol adds no physics), plus a bit-exact heterogeneous-job parity gate
vs direct `CpuFdtd`. Remaining candidates: spec-mask overlay (F1.3 proper), Touchstone
export of engine-measured responses, complex/de-embedded S-parameters, live-streamed
visualization over WS, real-GPU nightly numbers.
Earlier same day — engine track COMPLETE through E.5 (ADR-0179): E.4 closed
(row-sliced kernels), E.5b shipped (on-GPU NTFF accumulation, 315.4 dB / 2.9e-7 cross-backend),
E.5c shipped (dispersive ADE, bit-exact CPU + differential GPU gate). Python bindings
`yee.compute` shipped (ADR-0178). Studio track underway: S.0 `yee-engine` job API SHIPPED,
S.2 Tauri 2 + React skeleton SHIPPED (47.9 kB gzipped frontend; cargo check + vite build green
in-container). Later same day (ADR-0180): S.1 `yee-server` SHIPPED (WS job streaming +
cancel-on-disconnect + `yee serve`; e2e WS gates) and the S.3 visualization skeleton SHIPPED
(engine field-slice → heatmap + DFT spectrum, vitest/DOM gates, `studio-build` CI job).
Still later (ADR-0181): S.3b SHIPPED (three.js field surface, lazy-chunked, pure-function
geometry gates) and S.4 AUDITED (capability table; Dioxus retirement deferred with a defined
convergence path — engine-powered verify over the S.0/S.1 protocol). The S.* track as
originally scoped is complete; next candidates: engine-powered filter verify ("S.5"),
live-streamed visualization over WS, real-GPU perf numbers from the nightly.
Earlier: E.2/E.3/E.5a (ADR-0177), E.1 (ADR-0176), E.0 (ADR-0175).*
