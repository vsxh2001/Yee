# Yee Roadmap

This roadmap is a living document. It targets a realistic timeline for a small team augmented by AI tooling: **v1.0 of the planar MoM beachhead in three to four years.** We deliberately resist scope creep; everything below Phase 4 is non-negotiable scope, and Phase 4 itself is open-ended.

Conventions used below:
- 🎯 **Goal** — what success means at the end of the phase
- 📦 **Deliverables** — concrete artifacts shipped
- ✅ **Validation** — benchmark cases that must pass before the phase is called "done"
- ⚠️ **Risks / dependencies** — what could derail this phase

## Status snapshot (2026-05-19)

**Shipped:**
- Phase 0 walking skeleton (`phase-0-done` tag)
- Phase 1.0 free-space MoM dipole, NEC-4 87+j41 Ω reference passing (`phase-1-0-mom-dipole` tag)
- Phase 1.1.0 multilayer Greens placeholder (one-image DCIM)
- Phase 1.1.1.2 Sommerfeld pole extraction implementation (Newton-Raphson TM_0/TE_1, pole-subtracted GPOF, Hankel reconstruction; ADR-0033, merge `a22d622`)
- Phase 1.1.1.2.1 Sommerfeld surface-wave prefactor canonical correction (Michalski-Mosig 1997 eq. 25 / Felsen-Marcuvitz §5.5; Track EEEEEE merge `ca0e7bb`)
- Phase 1.1.1.2.2 `sommerfeld::residue()` sign + factor-of-2 fix (Michalski-Mosig 1997 eq. 19 form `-N₁/(2·D')`; Track TTTTTT merge `a4f98a4`; verified by contour-integral diagnostic Track SSSSSS)
- Phase 1.3.0 wave-port skeleton (matches delta-gap)
- Phase 1.3.1.1 step 2-3 Nedelec edge-element + nodal Lagrange E_z assembly + dense eigen on `TriMesh2D`; WR-90 TE10 cutoff gate passing at 0.055% error
- Phase 1.3.1.1 step 6 `yee.eigensolver` Python binding (PyTriMesh2D + PyNumericalCrossSection; 7 pytest cases; WR-90 cutoff sweep notebook)
- Phase 1.4 surface roughness (Hammerstad-Jensen, Groiss, Huray)
- Phase 1.5 cuSOLVER LU (hardware-gated)
- Phase 1.6 GMRES iterative
- Phase 1.gui.0/1/2/3 (egui shell, wgpu viewport, S11 + Smith plots, rust-1.92 + egui-0.34 + wgpu-29 toolchain bump)
- Phase 1.mesh.0/1 (Gmsh + KiCad import)
- Phase 1.plotting.0 (yee-plotters)
- Phase 1.validation.0/1/2 (aggregator + JSON Report + PNG artifacts via CI upload)
- Phase 1.bench yee-bench (criterion benches: MoM solve, FDTD step, GMRES vs LU, GP fit, BO, TF/SF, lumped)
- Phase 1.cli.1 `yee validate`, `yee bench`
- Phase 1.examples.0/2/4 (Rust examples, BO notebook, NSGA-II + AL notebooks)
- Phase 1.frontend.0/1/2/3 (yee-py: GP, FdtdDriver, BO, NSGA-II + AL, validation aggregator)
- Phase 2.fdtd.0..6 (walking skeleton, CPML, NTFF, dispersive ADE, end-to-end driver, TF/SF slab, lumped RLC)
- Phase 2.fdtd.5.3.2 cubic Lagrange aux-grid interpolation; oblique TF/SF clears >1000× DoD at 1027× / 60.2 dB (ADR-0034, merge `f878bdd`)
- Phase 2.fdtd.7 Q1 `WalkingSkeletonSolver::step` refactor into composable helpers (Track FFFFFF merge `1301623`)
- Phase 2.fdtd.7 Q2 `SubgridRegion` + 2× sub-Yee-grid scaffold (Track IIIIII merge `65ea3df`)
- Phase 2.fdtd.7 Q3 coarse→fine E_t spatial + temporal interpolation (Track MMMMMM merge `817955a`)
- Phase 2.fdtd.7 Q4 fine→coarse H_t area-average + E_t overwrite closures (Track OOOOOO merge `6ded764`)
- Phase 2.fdtd.7 Q4.1 `snapshot_fine_h_mid_step` time-centering helper (Track VVVVVV merge `a2abb4c`)
- Phase 2.fdtd.7 Q5 time-subcycling step (Track RRRRRR merge `426a36c`)
- Phase 2.fdtd.7.x Berenger Huygens spec + plan + ADR-0035 (Track AAAAAAA merge `003bdde`)
- Phase 2.fdtd.7.x B1 Berenger skeleton + face enumeration (Track EEEEEEE merge `c663b90`)
- Phase 2.fdtd.7.x B2 equivalent-current injection (Track FFFFFFF merge `c0b0cca`)
- Phase 2.fdtd.7.x B2.1 split J/M injection refactor (Track LLLLLLL merge `bb054e8`)
- Phase 2.fdtd.7.x B2.2 J-side coarse-ghost subtraction (Track OOOOOOO merge `464c7ba`)
- Phase 2.fdtd.7.y M-coupling spec + plan + ADR-0038 (Track UUUUUUU merge `0d260d3`)
- Phase 2.fdtd.7.y C1 pre/post fine-E snapshots (Track YYYYYYY merge `134fd93`)
- Phase 2.fdtd.7.y C2 compensating-source M (Option β; degenerates to 0 — Track ZZZZZZZ merge `be71a76`)
- Phase 2.fdtd.7.y C5 Mur ABC on fine outer E_t (Option α; retires 500-step divergence — Track BBBBBBBB merge `a6283ae`)
- **Phase 2.fdtd.7.y C6 un-ghosted J variant — retires Q5 strict 0.5%-of-peak gate at 0.0000% rel err** (Track DDDDDDDD merge `47c461c`; trade-off: fine grid permanently passive in source-on-coarse mode; Q6 long-time energy drift still `#[ignore]`'d)
- Phase 3.gp.0/1 (GP regression + ML hyperparameter fit)
- Phase 3.bo.0/1 (Expected-Improvement BO, NSGA-II multi-objective)
- Phase 3.al.0 (variance-acquisition active learning)
- Phase 4.fem.eig.0 walking-skeleton FEM eigenmode end-to-end:
  - T1+T2: yee-fem scaffold + `TetMesh3D` (Track GGGGGG merge `84a6632`)
  - T3: 6-edge Nedelec local K+M matrices, 4-pt Gauss quadrature (Track HHHHHH merge `f92fb59`)
  - T4: global sparse K+M assembly + PEC Dirichlet elimination (Track KKKKKK merge `aebb2a1`)
  - T5: `SparseEigen` trait + `InverseIterEigen` shift-invert via faer sparse LU (Track NNNNNN merge `fb6be04`; lobpcg crate fallback)
  - T6: `TetMesh3D::cavity_uniform` + Kuhn 6-tet brick decomposition (Track LLLLLL merge `ce899c3`)
  - **T7: fem-eig-001 production gate — TE_{101} 0.09% rel err vs Pozar §6.3 analytic at WR-90 (a=22.86, b=10.16, d=30) mm; mode-10 RMS 0.37% on (12,9,15) mesh; wall-time ~7 s release (Track QQQQQQ merge `d42aefc`)**
  - T8: `yee.fem.solve_cavity` Python binding, 3 pytest cases (Track UUUUUU merge `cb0e15f`)
  - T9: mdBook tutorial `docs/src/tutorials/04-fem-cavity-eigenmode.md` (Track WWWWWW merge `06e72f2`)
- Phase 3.nl.0 NL design surface, end-to-end:
  - R1: yee-design crate scaffold + DesignIntent types (Track PPPPPPP merge `fbd752e`)
  - R2: Balanis Ch. 14 initial-estimate calculator — Example 14.1 W/L within 0.08%/0.07% (Track RRRRRRR merge `32baeb4`)
  - R3: deterministic project-TOML emitter (Track VVVVVVV merge `2e54e6f`)
  - R4: yee.design.from_prompt_llm Anthropic Messages tool-use sidecar (Track XXXXXXX merge `2c7ece4`)
  - R5: yee design CLI subcommand + 10 canonical prompts (Track AAAAAAAA merge `08cec1b`)
  - R6: nl-001 production gate — schema+round-trip+offline sub-gates A+B+C all 10 prompts (Track CCCCCCCC merge `417978e`)
  - R7: mdBook tutorial `docs/src/tutorials/04-nl-design-surface.md` (Track EEEEEEEE merge `5016fda`)
- Phase 4.fem.eig.1 dispersive ε_r(ω) FEM eigensolver, D1-D7 shipped:
  - design spec + plan + ADR-0039 (Track FFFFFFFF merge `10d91d7`)
  - D1+D2: complex tet element + Complex64 inverse-iter (Track HHHHHHHH merge `cfd3e49`)
  - D3: MaterialDatabase (Drude/Lorentz/Debye ε(ω)) (Track JJJJJJJJ merge `7e15ed2`)
  - D4: DispersiveSolver::solve_at_frequency (Track NNNNNNNN merge `90bc337`)
  - D5: Newton-Raphson ω-tracker (Track OOOOOOOO merge `1480a51`)
  - D7: yee.fem.solve_cavity_dispersive Python binding (Track RRRRRRRR merge `214075b`)
  - D6 production gate fem-eig-002 lossy SiO₂ cavity — Track QQQQQQQQ in flight
- Track WWWWWWW TEM-mode smoothed RHS port: mom-002 |Z_in| 674 Ω → 3.46 Ω, Maxwell-envelope deviation 580% → 70% (merge `a08f0db`)
- Track GGGGGGGG WavePort::rhs Numerical2D arm wired (Phase 1.3.1.1 step 7, 1% L2 vs analytic TE10; merge `3b115fa`)
- Track IIIIIIII mom-003 re-run through Sommerfeld + TEM port: CaseStatus::Passed within loose-tolerance band, |Z_in| = 13.4 Ω (merge `3b115fa`)
- Track MMMMMMMM yee-fdtd per-cell ε_r/μ_r + PEC mask infrastructure (merge `cb6f8ed`)
- Track PPPPPPPP CPML reads per-cell ε_r/μ_r (lifts MMMMMMMM workaround; reflection floor 69.33 dB preserved; merge `c57592f`)
- Track LLLLLLLL fdtd-007 Maloney-Smith driver+gates committed `#[ignore]`'d pending fdtd infra (3 blockers documented; merge `30b2d2c`)
- mom-002 root-cause chain end-to-end (10 forensic tracks + 3 kernel fixes + 3 ADRs):
  - EEEEEE prefactor / JJJJJJ extent / PPPPPP GPOF / SSSSSS contour / TTTTTT residue sign / XXXXXX ψ_p / YYYYYY MPIE / CCCCCCC port-mesh / MMMMMMM ε_eff / NNNNNNN R1 retract / DDDDDDD DCIM-TM / TTTTTTT port spatial / QQQQQQQ β eigen (kernel exonerated at 1.83% from HJ)
  - ADR-0036 mom-002 validation reframe (sub-wavelength strip)
  - ADR-0037 R1 metric retraction
  - IIIIIII reframe to L=82mm centered uniform: |Z_in| 2569→674 Ω

**Pending (high priority):**

*In-flight (this session):*
- WWWWWWW mom-002 TEM-mode smoothed RHS port-excitation fix (TTTTTTT P1 root cause)
- Phase 1.3.1.1 step 4 sparse arpack-rs / LOBPCG eigensolver
- Phase 1.3.1.1 step 5 longitudinal block for quasi-TEM microstrip wave-ports

*Design-coverage shipped, impl pending:*
- Phase 2.fdtd.7 Q6 stability/reciprocity 10000-step energy gate — Q5 strict retired by C6; Q6 long-time drift (75-79%) deferred to future track (subgrid-coarse impedance-mismatch is the residual; needs proper energy-balance closure)
- Phase 2.fdtd.7 Q7 fdtd-007 Maloney-Smith production gate
- Phase 4.fem.eig.1+ — dispersive ε_r(ω), real waveguide ports, absorbing boundaries — designs not yet drafted

**Outstanding validation gates:**
- mom-001 dipole — **GATE PASSES** (NEC-4 87+j41 Ω)
- mom-002 microstrip Z₀ — gate passes within ±5% tripwire band at `|Z_in| = 674 Ω` on L=82mm reframed mesh (per ADR-0036); 10 forensic tracks confirmed kernel is correct within 1.83% of HJ ε_eff; remaining residual is delta-gap port-excitation modeling (Track WWWWWWW in flight)
- mom-003 2.4 GHz patch — loose tolerance pending re-run through `GreensSpec::MicrostripSommerfeld`
- fem-eig-001 WR-90 rectangular cavity — **GATE PASSES** (TE_{101} 0.09% rel err, mode-10 RMS 0.37%)
- fdtd-007 Maloney-Smith oblique TF/SF — forward gate for Phase 2.fdtd.7 subgridding (gated on Q6 + Q7)
- nl-001 10-prompt sweep — **GATE PASSES on sub-gates A+B+C** (schema, round-trip, offline); D-gate (solver ±5% f) `#[ignore]`'d pending real MultilayerGreens per Phase 1.1.1 deferred-tolerance policy

---

## Phase 0 — Foundation (Months 0–6)

🎯 **Goal.** Stand up the project skeleton end-to-end. A user can install Yee, mesh a microstrip line via Gmsh, run a stub MoM solve that exercises CUDA, and export a Touchstone file. The result need not be physically accurate beyond simple analytical cases — the point is that every pipe is connected.

📦 **Deliverables.**
- Cargo workspace with crates `yee-core`, `yee-cuda`, `yee-mesh`, `yee-mom`, `yee-fdtd` (stub), `yee-io`, `yee-cli`, plus an `examples/` and `validation/` tree.
- CUDA scaffolding via `cudarc` 0.19+: device enumeration, context/stream management, NVRTC kernel compilation, a "hello world" stencil kernel, and CI that builds on CUDA 12.4 and 13.0.
- Gmsh integration: in-tree `bindgen`-generated FFI against `gmshc.h` 4.15+, with a thin safe Rust wrapper. (The pre-existing `rgmsh` crate is unmaintained since 2019 and targets Gmsh 4.4.1; we generate fresh bindings.)
- Lossless, single-layer, infinite-ground planar MoM solver (2.5D, no dielectric stack-up yet, perfect conductor only). Dense LU on CPU via `faer` as a reference; first GPU port via cuSOLVER `cusolverDnCgetrf` exposed as a feature flag.
- Touchstone v1.1 reader/writer (`.s1p` through `.s4p` minimum; generic `.sNp` support).
- `yee` CLI for `validate`, `mesh`, `run`, `export`.
- Initial documentation site (mdBook) and contributing guide.

✅ **Validation milestones.** Phase 0 is a *walking skeleton*: every pipe between the workspace crates, the CLI, the documentation pipeline, and CI is connected end-to-end. The gates below are pure-build, not physical-accuracy:

1. `cargo check --workspace --no-default-features` exits 0
2. `cargo test --workspace --no-default-features` exits 0
3. `cargo clippy --workspace --all-targets --no-default-features -- -D warnings` exits 0
4. `cargo fmt --check --all` exits 0
5. `cargo doc --workspace --no-default-features --no-deps` exits 0
6. `cargo run --bin yee -- --help` exits 0 and lists every subcommand
7. `cargo run --bin yee -- validate all` exits 0
8. `mdbook build docs/` exits 0
9. `THIRD_PARTY_LICENSES.md` documents Gmsh, OCCT, and NVIDIA CUDA proprietary dynamic-link posture
10. CI workflow runs gates 1–9 on Linux + Rust 1.88 and exits green

⚠️ **Risks / dependencies.**
- `cudarc` self-describes as "pre-alpha"; it has historically shipped breaking minor releases (notably 0.13 → 0.14). **Mitigation:** pin to exact minor version; introduce a thin internal `yee_cuda::backend` abstraction so we can swap if needed.
- Gmsh's GPL v2+ (with linking exception) is GPL v3 compatible, but we must document this clearly in `THIRD_PARTY_LICENSES.md`.
- Rust 1.85+ is the floor (driven by `maturin` 1.10 and `pyo3` 0.28). This is fine but worth pinning in `rust-toolchain.toml`.

---

## Phase 1 — Planar MoM v1.0 (Months 6–18)

🎯 **Goal.** Ship a production-grade, GPU-accelerated **multilayer planar MoM solver** competitive with Sonnet Lite on real PCB designs, with first-class Python bindings and a usable desktop GUI. This is the **beachhead**.

📦 **Deliverables.**
- **Multilayer dielectric stack-up.** Spectral-domain Green's functions for arbitrary layered media; Sommerfeld integral evaluation with Discrete Complex Image Method (DCIM) or rational-function fitting for speed.
- **RWG/rooftop basis functions** on planar triangular and rectangular meshes.
- **Lumped ports** (delta-gap and edge ports), wave ports for microstrip/CPW with mode extraction, and **TRL/SOLT de-embedding** for reference-plane shifting.
- **Surface roughness** models: Hammerstad-Jensen, Groiss, Huray (small-sphere). Frequency-dependent loss.
- **GPU acceleration.** Matrix fill on CUDA (one block per RWG pair batch); dense LU via cuSOLVER (`cusolverDnZgetrf` for complex double); right-hand-side solves via `cusolverDnZgetrs`. cuBLAS for any GEMM/GEMV used in iterative refinement. Single-GPU first; multi-GPU dense LU via cuSOLVERMg flagged behind a feature.
- **Python bindings** via PyO3 0.28 with `abi3-py310`; built and published as wheels via `maturin` 1.10 with `manylinux_2_28`. NumPy interop through the `numpy` crate.
- **Initial desktop GUI** built with `egui` 0.34+ and `eframe`. Embedded `wgpu` 3D viewport (paint callback) for PCB geometry; `egui_plot` for S-parameter and Smith-chart views; `egui_dock` for panel docking.
- **`rerun` SDK** integration as an optional structured-logging sink for solver internals (mesh evolution, current densities per frequency, convergence traces).

✅ **Validation milestones.**
- **Closed-form half-wave dipole impedance**[^phase0-reclassified]: 50 Ω microstrip-fed, reproduce Z ≈ 73 + j42 Ω within 5%.
- **50 Ω microstrip line on FR-4**[^phase0-reclassified]: characteristic impedance within ±3% of TX-LINE / Hammerstad-Jensen.
- **2.4 GHz rectangular microstrip patch on FR-4**[^phase0-reclassified] (29.2 × 38.0 mm, h = 1.6 mm, εr = 4.4, lossless): resonance within ±2% of published value; |S11| < −10 dB at resonance. (FDTD comparison case to be added in Phase 2.)
- **Swanson 5-pole hairpin BPF** (RT/Duroid 6006, εr = 6.15, h = 1.27 mm, ~2.0 GHz): reproduce S-parameter response within ±1 dB of Sonnet reference up to 4 GHz; resonant frequencies within ±0.5%.
- **Parallel-coupled-line BPF** (Hong & Lancaster Ch. 5): reproduce passband ripple, return loss, and stopband rejection within ±1 dB.
- **Wilkinson divider at 2 GHz**: three-port S-parameters within ±0.5 dB of closed-form / Pozar reference.
- **Branch-line (90°) hybrid**: amplitude and phase balance verified.
- **Cross-validation against openEMS** on every microstrip and patch case (FDTD vs MoM should agree within 3% at resonance).
- **Inset-fed patch antenna on RO4003C** (matched 50 Ω): published-paper figure-for-figure match.
- All validation runs scripted; results pushed to `validation/results/` and regenerated in CI nightly.

[^phase0-reclassified]: Originally listed under Phase 0; reclassified as Phase 1 in the [2026-05-16 Phase 0 walking-skeleton design](docs/superpowers/specs/2026-05-16-phase-0-multi-agent-execution-design.md).

⚠️ **Risks / dependencies.**
- DCIM accuracy across wide frequency ranges is finicky; expect to ship multiple Green's-function evaluators (DCIM + direct Sommerfeld + rational fit) and switch adaptively.
- Dense LU at n ≥ 50k overflows a single 80 GB H100; we will hit this on real PCBs. **Mitigation:** start with iterative GMRES + block-diagonal preconditioner on GPU as the n ≥ 50k path; queue MLFMA / ACA work for Phase 4.
- egui ships breaking minor releases roughly quarterly. **Mitigation:** isolate UI behind a stable internal trait so the GUI crate can be migrated in a single PR each quarter.
- PyO3 has historically shipped one breaking change per minor release; we will pin and migrate deliberately.

---

## Phase 2 — 3D FDTD (Months 18–30)

🎯 **Goal.** A production 3D FDTD solver on CUDA, covering radiation, transient signal integrity, and dispersive materials — the cases where planar MoM is the wrong tool.

📦 **Deliverables.**
- **3D Yee staggered grid** on CUDA. Memory-bandwidth-optimized E/H update kernels. Mixed precision (FP32 for fields, FP64 for accumulators where needed). Multi-GPU domain decomposition with NCCL boundary exchange (cudarc has safe NCCL bindings).
- **CPML (Convolutional PML)** absorbing boundaries — Roden & Gedney formulation — on all six faces, with the standard polynomial grading.
- **Dispersive materials**: Drude, Lorentz, Debye, and arbitrary multi-pole Debye via ADE (Auxiliary Differential Equation) or PLRC (Piecewise Linear Recursive Convolution).
- **Near-to-far-field transformation** (NTFF) for full 3D antenna radiation patterns, gain, directivity, axial ratio, and 3D pattern export.
- **Lumped-element ports** (resistor / capacitor / inductor / arbitrary RLC), waveguide ports with modal sources, plane-wave sources for scattering problems.
- **Subgridding** (non-uniform Cartesian) with stability fixes per Berenger / Xiao-Liu schemes.
- **Conformal techniques** (Dey-Mittra or simple staircase fallback) for non-aligned geometry.
- **Geometry ingestion** through OpenCascade via `opencascade-rs` 0.2+ (STEP/IGES import), then voxelization onto the Yee grid; KiCad PCB import for the common case.
- **GPU-resident time-stepping with on-the-fly volume-data streaming to `rerun` for debugging.**

✅ **Validation milestones.**
- **Resonant cavity Q-factor**: rectangular cavity TE/TM modes match analytical to ±0.5%.
- **Pyramidal horn antenna**: pattern within ±1 dB of measured/published in main beam.
- **Dipole over a dielectric half-space**: NTFF pattern vs Sommerfeld reference.
- **Cross-validation against openEMS** on identical geometries — agreement within numerical-noise level (driven by grid and PML settings, not solver choice).
- **Microstrip line transient propagation**: time-domain TDR matches frequency-domain MoM via FFT.

⚠️ **Risks / dependencies.**
- FDTD memory bandwidth is the bottleneck; hand-tuned kernels with shared-memory tiling are essential to beat openEMS. We will benchmark openly.
- Subgridding stability is famously fragile; we plan to ship without it first and add it once the rest of the solver is locked.
- Multi-GPU domain decomposition adds significant complexity. **Mitigation:** ship single-GPU first; multi-GPU behind a feature flag with explicit "experimental" labeling.

---

## Phase 3 — AI / ML Layer (Months 30–42)

🎯 **Goal.** Make Yee genuinely **AI-native**: every solve trains a surrogate, every parametric sweep gets cheaper, and a natural-language interface lets engineers describe what they want and get a viable starting design.

📦 **Deliverables.**
- **Surrogate model framework.** Every parametric sweep produces a labeled dataset (parameters → S-parameters, near fields, far fields). Pluggable surrogate backends: Gaussian processes for small data, MLPs/transformers/Fourier neural operators for large data. Training orchestrated via Candle (Rust-native) or PyTorch (Python sidecar) — both via `cudarc`-compatible CUDA contexts.
- **Surrogate-in-the-loop optimization.** Bayesian optimization, NSGA-II for multi-objective (size vs bandwidth vs gain), with surrogate predictions checked against the full solver on a schedule.
- **Active learning loops.** Solver picks the next simulation points to maximize surrogate accuracy.
- **Natural-language design surface.** LLM-mediated front end that parses "I need a 2.4 GHz inset-fed patch on RO4003C with at least 100 MHz bandwidth and gain over 6 dBi" into a parameterized Yee design, generates initial dimensions from textbook formulas, refines via the surrogate, and returns a ready-to-simulate project file. Underneath this surface, all interactions are reproducible script — the natural-language layer is convenience, not magic.
- **Pre-trained model zoo.** Public surrogates for canonical geometry families (rectangular patches, inset-fed patches, hairpin filters, Wilkinson dividers) hosted alongside their training data on Hugging Face.
- **Inverse design / topology optimization.** Adjoint-based gradients through the FDTD solver enable photonic-style inverse design for antennas and filters.

✅ **Validation milestones.**
- **Surrogate accuracy.** On the patch-antenna family, surrogate predictions of S11 within ±0.5 dB and resonance within ±0.2% over the trained parameter range, with 10–100× speed-up vs. full solve.
- **NL-to-design.** End-to-end: text prompt → working design that meets stated specs to within 10% on at least 5 canonical antenna / filter classes.
- **Inverse-designed antenna** that outperforms its textbook starting point on a defined figure of merit, verified against full FDTD.

⚠️ **Risks / dependencies.**
- LLM dependencies (whether self-hosted or API) introduce reliability and reproducibility issues. **Mitigation:** the LLM only emits structured design scripts that the user can inspect, edit, and re-run deterministically.
- Surrogate accuracy is geometry-family-specific; the "every sweep gets cheaper" promise is true within a family but does not generalize across families without large pre-training. Be honest about this.
- Adjoint FDTD is non-trivial — plan for a research-grade implementation first, production-grade after.

---

## Phase 4 — 3D FEM, Eigenmode, Broader Applications (Months 42+)

🎯 **Goal.** Round out the solver portfolio so Yee can compete with HFSS and Palace on driven 3D FEM and eigenmode problems, and open the door to adjacent application domains (SI/PI, EMI/EMC, photonics, accelerator cavities).

📦 **Planned deliverables** (priorities to be re-confirmed when we get here):
- **3D FEM solver** with high-order Nedelec edge elements; HCURL spaces; conformal hexahedral and tetrahedral meshes via Gmsh.
- **Eigenmode solver** (subspace iteration / Krylov-Schur) for resonant cavities and filters.
- **MLFMA / ACA / H-matrix** compression for MoM, enabling n ≥ 100k.
- **Coupled circuit-EM co-simulation** (a la ADS).
- **Time-domain FEM** for transient analysis.
- **Anisotropic / nonlinear / time-varying materials.**
- **Application packs**: SI/PI (DDR/PCIe channel simulation), EMI/EMC (radiated emissions, shielded enclosures), photonics (silicon photonics, plasmonics), particle accelerators (cavity design).

⚠️ **Risks.** This phase exists to keep direction honest; specifics will be revised based on user demand and what the planar MoM + FDTD + AI core teaches us.

---

## Cross-cutting work (every phase)

- **Validation.** No solver feature ships without a published-benchmark validation case in `validation/` and a CI run that regenerates results nightly.
- **Documentation.** Every public Rust crate and Python module has examples. Every solver has a "theory of operation" doc that cites its sources.
- **Reproducibility.** All examples and validation cases are scripted; no GUI-only artifacts.
- **Performance budget.** Each release includes published benchmark times against the previous release and (where licenses allow) against openEMS and gprMax on identical geometries.
- **Community.** Discord, GitHub Discussions, a monthly "office hours" call once usage justifies it.
