# Yee Roadmap

This roadmap is a living document. It targets a realistic timeline for a small team augmented by AI tooling: **v1.0 of the planar MoM beachhead in three to four years.** We deliberately resist scope creep; everything below Phase 4 is non-negotiable scope, and Phase 4 itself is open-ended.

Conventions used below:
- 🎯 **Goal** — what success means at the end of the phase
- 📦 **Deliverables** — concrete artifacts shipped
- ✅ **Validation** — benchmark cases that must pass before the phase is called "done"
- ⚠️ **Risks / dependencies** — what could derail this phase

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

[^phase0-reclassified]: Originally listed under Phase 0; reclassified as Phase 1 in the 2026-05-16 Phase 0 walking-skeleton design (`docs/superpowers/specs/2026-05-16-phase-0-multi-agent-execution-design.md`).

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
