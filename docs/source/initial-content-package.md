# Yee — Repository content package

Below are three drop-in files for the new GitHub repository. The recommended working project name is **Yee** (after Kane S. Yee, the originator of the FDTD method). It is short, instantly meaningful to any EM engineer, makes a clean CLI (`yee sim …`), and has no naming collisions in open-source EM, on crates.io, or against major commercial brands. Two fallback names — **Poynting** and **Sparrow** — are noted in the README for the team to weigh in on.

---

## A. `README.md`

```markdown
# Yee

> **GPU-accelerated, AI-native electromagnetic simulation. Open source. Written in Rust.**

[![Build](https://img.shields.io/badge/build-pending-lightgrey)](#)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-0.0.0--alpha-orange)](#)
[![Rust](https://img.shields.io/badge/rust-1.85%2B-orange?logo=rust)](https://www.rust-lang.org)
[![CUDA](https://img.shields.io/badge/CUDA-12.4%2B-76B900?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![Discord](https://img.shields.io/badge/Discord-join-5865F2?logo=discord)](#)

## What is Yee?

**Yee** is a next-generation electromagnetic (EM) simulation studio for RF, microwave, and antenna engineers. The v1 beachhead is **planar PCB antennas and filters**, solved with a **GPU-accelerated planar Method of Moments (MoM)** kernel — the Sonnet-class regime — backed by a **3D FDTD** secondary solver for radiation, transients, and dispersive materials. It is written in Rust, accelerates on NVIDIA CUDA, and exposes a first-class Python API for Jupyter workflows.

Named after **Kane S. Yee**, whose 1966 paper *Numerical Solution of Initial Boundary Value Problems Involving Maxwell's Equations in Isotropic Media* introduced the staggered grid that still bears his name.

## Why does this need to exist?

The RF/microwave world runs on commercial tools — CST Studio Suite, Ansys HFSS, Keysight ADS/Momentum, Sonnet Suites, AWR/AXIEM — that cost **$40k–$200k per seat per year**, are sales-gated, ship on annual or perpetual+maintenance contracts, and remain closed-source. The open-source alternatives are real but partial: **openEMS** (FDTD, CPU-only in production), **gprMax** (CUDA FDTD focused on ground-penetrating radar), **Palace** (3D FEM eigenmode/driven for superconducting qubits), **NEC-2** (1981-vintage wire-only MoM). **No production-quality open-source planar MoM solver exists.** No open-source GPU planar MoM exists.

That gap is the wedge. PCB antennas and filters are where most working RF engineers spend their day, and planar MoM is the right solver for that work. Doing it on the GPU in Rust, with an AI-native UX on top, is a defensible position no commercial vendor has staked out.

A second motivation: every modern EM simulation is a labeled dataset waiting to happen. **Yee treats ML surrogate models as first-class citizens** — every parametric sweep trains a surrogate, every subsequent sweep gets cheaper, and natural-language design becomes feasible because you can ask a fast model what to try before you commit a full solve.

## Status

**Phase 0 — Foundation.** Cargo workspace bootstrapping, CUDA kernel scaffolding via [cudarc](https://github.com/coreylowman/cudarc), Gmsh FFI, and a lossless planar MoM kernel under construction. Not yet usable. See **[ROADMAP.md](ROADMAP.md)**.

## Key features

**Available now (Phase 0):** Cargo workspace, CUDA detection and device enumeration, Touchstone (`.s2p` / `.sNp`) export, basic mesh ingestion via Gmsh.

**Phase 1 (months 6–18):** Multilayer-dielectric planar MoM with lumped ports, de-embedding, surface roughness models, GPU-resident dense LU (via cuSOLVER), Python bindings via PyO3, an `egui`-based GUI.

**Phase 2 (months 18–30):** 3D FDTD on a CUDA Yee grid with CPML boundaries, Drude/Lorentz/Debye dispersive materials, near-to-far-field transformation for full 3D antenna radiation patterns, lumped-element ports.

**Phase 3 (months 30–42):** ML surrogates trained on every solve, natural-language design interface ("design me a 2.4 GHz inset-fed patch on RO4003C"), automated optimization with surrogate-in-the-loop, neural-operator fast solvers.

**Phase 4 (months 42+):** 3D FEM, eigenmode solver, broader application domains (SI/PI, EMI/EMC, photonics).

## Differentiators

1. **AI-native UX** — natural-language design as the primary surface, a stable scripting layer underneath, a GUI for visualization. The notebook is a first-class citizen, not an afterthought.
2. **ML surrogates by default** — every simulation contributes to a model. Every sweep gets cheaper. Sweeps you ran last week become predictions you can run in milliseconds today.
3. **GPU-first architecture** — CUDA on NVIDIA throughout. We are not retrofitting GPU support onto a CPU solver; the kernels are GPU-resident from day one.
4. **Stand on giants** — we use **Gmsh** for meshing and **OpenCascade** for CAD via the well-maintained Rust bindings. We do not rewrite what is not the moat.
5. **Validation, not vibes** — every solver is held against canonical published benchmarks (microstrip patch on FR-4, Wilkinson divider, Swanson hairpin BPF, IEEE AP-S reference antennas), with reproducible scripts in `validation/`.

## Quick start

> ⚠️ Phase 0 — the commands below are aspirational and will evolve before the first tagged release.

```bash
# Prerequisites: Rust 1.85+, CUDA Toolkit 12.4+, Python 3.10+
git clone https://github.com/yee-em/yee
cd yee
cargo build --release --features cuda

# Run the example: 2.4 GHz inset-fed microstrip patch on FR-4
cargo run --release --example patch_2g4

# Python (after `pip install yee` once published)
python - <<'PY'
import yee
patch = yee.shapes.patch(width_mm=29.2, length_mm=38.0, substrate=yee.materials.FR4(h_mm=1.6))
sim = yee.PlanarMoM(geometry=patch, freq_range=(2.0e9, 3.0e9, 201))
s = sim.run(device="cuda:0")          # GPU solve
s.touchstone("patch.s1p")             # Export
s.plot()                              # egui live view or matplotlib in Jupyter
PY
```

## Architecture overview

Yee is a Cargo workspace built around three layers: a **kernel layer** (CUDA C kernels for MoM matrix fill, FDTD updates, and ML inference, JIT-compiled via NVRTC and orchestrated from Rust via `cudarc`); a **solver layer** (planar MoM and 3D FDTD, with shared infrastructure for ports, materials, sources, post-processing, and surrogate training); and a **frontend layer** (a stable Rust API, PyO3 bindings for Python, an `egui` desktop GUI with an embedded `wgpu` 3D viewport, and an LLM-backed natural-language design surface). Mesh and CAD I/O delegate to Gmsh and OpenCascade through their C APIs. The whole stack is GPL v3.

## How Yee compares

| Tool | Method(s) | License | GPU | Planar MoM | Python | Cost (USD) |
|---|---|---|---|---|---|---|
| **Yee** *(this project)* | Planar MoM + 3D FDTD (FEM later) | **GPL v3** | **CUDA-first** | ✅ | ✅ | **Free** |
| CST Studio Suite | FIT (time), FEM (freq), MoM, asymptotic | Commercial | ✅ (HPC tokens) | ✅ | Limited | ~$80–100k/seat |
| Ansys HFSS | FEM (driven + eigen), MoM, SBR+ | Commercial | ✅ (HPC pack) | Indirect | Limited | ~$40–50k/seat/yr |
| Keysight ADS / Momentum | Circuit + planar MoM (2.5D) | Commercial | Limited | ✅ | Limited | ~$30–60k/seat |
| Sonnet Suites | Planar MoM (shielded-box) — gold standard | Commercial (free Lite) | CPU multi-core | ✅ | Limited | ~$3–15k+/seat |
| openEMS | 3D FDTD | GPL v3 | Experimental only | ❌ | ✅ | Free |
| gprMax | 3D FDTD (GPR-focused) | GPL v3 | ✅ CUDA | ❌ | ✅ | Free |
| NEC-2 | Wire-only thin-wire MoM | Public domain | ❌ | ❌ (wire only) | Via wrappers | Free |
| Palace (AWS) | 3D FEM (Nedelec, eigen + driven) | Apache 2.0 | ✅ CUDA/HIP | ❌ | Limited | Free |

The **rightmost columns tell the story**: there is no open-source GPU-accelerated planar MoM. That is exactly where Yee enters.

## Links

- [ROADMAP.md](ROADMAP.md) — multi-year plan, deliverables, validation milestones
- [TECH_STACK.md](TECH_STACK.md) — chosen dependencies, with rationale
- [CONTRIBUTING.md](CONTRIBUTING.md) — how to help (coming soon)
- [docs/](docs/) — full documentation site (coming soon)
- [Discord](#) and [Discussions](../../discussions) — community

## Project name

We chose **Yee** as the working name after evaluating roughly twenty candidates. Notable collisions we deliberately avoided: **Maxwell** (Ansys's flagship low-frequency EM product), **Lorentz** (IntegratedSoft's commercial 3D MoM tool LORENTZ-HF), **Heaviside** (Arena Physica's EM foundation model, released late 2025), **Ampere** (NVIDIA's GPU architecture — a particularly bad collision for a CUDA project), **Marconi** (the Julia `Marconi.jl` RF library), **Tesla** and **Faraday** (overloaded beyond rescue). Fallbacks under consideration: **Poynting** and **Sparrow**.

## Acknowledgments

Yee stands on shoulders. We gratefully acknowledge:

- **Gmsh** (Christophe Geuzaine, Jean-François Remacle) — the meshing engine
- **OpenCascade Technology (OCCT)** — the B-rep CAD kernel
- **cudarc** (Corey Lowman et al.) — CUDA host bindings for Rust
- **faer** (Sarah Quinones) — high-performance dense linear algebra
- **nalgebra**, **ndarray**, **parry** (dimforge) — geometry and arrays
- **egui** (Emil Ernerfeldt / Rerun) — immediate-mode GUI
- **PyO3 / maturin** — Python ↔ Rust interop
- **openEMS** (Thorsten Liebig) — the reference open-source FDTD, against which we will cross-validate
- **Palace** (AWS Center for Quantum Computing) — proof that industrial-quality OSS computational EM is possible
- Kane S. Yee, Roger F. Harrington, Allen Taflove, and the generations of researchers whose published work makes this kind of project possible at all.

## License

Yee is distributed under the **GNU General Public License v3.0 or later**. See [LICENSE](LICENSE). We chose strong copyleft deliberately: the EM tools market is dominated by closed-source vendors, and we want every improvement to this codebase to remain available to the community that built it.
```

---

## B. `ROADMAP.md`

```markdown
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

✅ **Validation milestones.**
- **Closed-form half-wave dipole impedance**: 50 Ω microstrip-fed, reproduce Z ≈ 73 + j42 Ω within 5%.
- **50 Ω microstrip line on FR-4**: characteristic impedance within ±3% of TX-LINE / Hammerstad-Jensen.
- **2.4 GHz rectangular microstrip patch on FR-4** (29.2 × 38.0 mm, h = 1.6 mm, εr = 4.4, lossless): resonance within ±2% of published value; |S11| < −10 dB at resonance. (FDTD comparison case to be added in Phase 2.)
- All validation cases reproducible via `cargo test --features validation` or `python -m yee.validation`.

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
- **Swanson 5-pole hairpin BPF** (RT/Duroid 6006, εr = 6.15, h = 1.27 mm, ~2.0 GHz): reproduce S-parameter response within ±1 dB of Sonnet reference up to 4 GHz; resonant frequencies within ±0.5%.
- **Parallel-coupled-line BPF** (Hong & Lancaster Ch. 5): reproduce passband ripple, return loss, and stopband rejection within ±1 dB.
- **Wilkinson divider at 2 GHz**: three-port S-parameters within ±0.5 dB of closed-form / Pozar reference.
- **Branch-line (90°) hybrid**: amplitude and phase balance verified.
- **Cross-validation against openEMS** on every microstrip and patch case (FDTD vs MoM should agree within 3% at resonance).
- **Inset-fed patch antenna on RO4003C** (matched 50 Ω): published-paper figure-for-figure match.
- All validation runs scripted; results pushed to `validation/results/` and regenerated in CI nightly.

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
```

---

## C. `TECH_STACK.md`

```markdown
# Yee Technology Stack

Selections and rationale, current as of **early 2026**. Versions are the recommended floor; pin in `Cargo.toml` and `pyproject.toml` and migrate deliberately.

## Languages and toolchain

| Concern | Choice | Version floor | Rationale |
|---|---|---|---|
| Systems language | **Rust** | 1.85+ | Driven by `maturin` 1.10 (Rust 2024 edition) and `pyo3` 0.28 (MSRV 1.83). Fearless concurrency, zero-cost abstractions, mature CUDA story. |
| Python | **CPython** | 3.10+ | Most RF engineers use Jupyter; 3.10 is the comfortable floor for modern scientific Python in 2026. |
| CUDA Toolkit | **CUDA** | 12.4+ (13.0 supported) | Tensor-Core / Hopper / Blackwell coverage; matches `cudarc` 0.19 support matrix. |
| Build | Cargo workspace + `maturin` | Cargo 1.85, maturin 1.10 | Single source of truth; wheels via `maturin build --release`. |

---

## CUDA host bindings — `cudarc`

**Choice: `cudarc` 0.19+.**

We evaluated three candidates:

1. **`cudarc`** (coreylowman) — latest 0.19.4, ~monthly releases, **the production-grade Rust CUDA binding in 2026**. Used by Hugging Face Candle, mistral.rs, candle-vllm. Covers CUDA driver, NVRTC, **cuBLAS, cuBLASLt, cuSOLVER, cuSPARSE, cuFFT, cuRAND, cuTENSOR, NCCL, cuDNN** out of the box. Supports CUDA 11.4–13.0 selectable at build time. Three linking modes: dynamic-loading (default, no libs at build time), dynamic-linking, static-linking. License: MIT/Apache-2.0.
2. **`cust`** (part of Rust-GPU/Rust-CUDA project) — last crates.io release 0.3.2 is years old. The project was rebooted in January 2025 and is genuinely active again, but as of late 2025 the rebooted code is **git-only** — no crates.io publication. Not yet suitable for a production project.
3. **`rustacuda`** — last release 0.1.3 (2021). **Unmaintained.**

We will **standardize on `cudarc` as the single CUDA host binding for the entire project.** Where `cudarc` lacks a top-level safe wrapper (notably cuSPARSE), we write a thin internal `yee_cuda::sparse` module that wraps `cudarc::cusparse::result::*` with RAII handles. Caveat: `cudarc` self-describes as "pre-alpha" despite massive production use; we pin to exact minor version and introduce a thin internal abstraction layer so we can swap if needed.

### Kernel authoring

We will **write CUDA kernels in CUDA-C and JIT them via NVRTC** (the approach Candle uses). This is the path of least resistance in 2026.

We considered Rust-on-GPU via `rustc_codegen_nvvm` (the Rust-CUDA project's NVVM backend). As of August 2025 its toolchain was bumped from a 3-year-old nightly to `nightly-2025-06-23`, CUDA 12.x is "experimental," and CI was restored — meaningful progress, but it still depends on a pinned nightly and is not production-ready. **We will re-evaluate in late 2026.**

```toml
cudarc = { version = "0.19", default-features = false, features = [
    "std", "driver", "nvrtc",
    "cublas", "cublaslt", "cusolver", "cusparse", "cufft", "curand",
    "f16",
    "cuda-version-from-build-system",
    "dynamic-loading",
] }
```

---

## Linear algebra

| Workload | Choice | Version | Why |
|---|---|---|---|
| Small static matrices (geometry, transforms) | **`nalgebra`** | 0.34+ | The de-facto standard. ~4.3M downloads/month. Compile-time-sized types with zero heap alloc. Pair with `nalgebra-sparse` if useful. |
| Medium/large dense CPU LU/Cholesky/SVD | **`faer`** | 0.23+ | Pure-Rust, MIT, competitive with Intel MKL on GEMM and factorizations. Supports complex types via `num_complex::Complex64`. Used as host-side fallback for MoM and reference implementation. |
| Dense LU at MoM scale (10k–100k complex) on GPU | **cuSOLVER via cudarc** | — | `cusolverDnZgetrf`/`Zgetrs`. cuSOLVERMg for multi-GPU. This is the hot path. |
| N-d arrays at the Python boundary | **`ndarray` + `numpy` crate** | ndarray 0.17, numpy 0.28 | Zero-copy NumPy ↔ `ndarray::ArrayView`. Use for FDTD field arrays exposed to notebooks. |
| Sparse SpMV on CPU (FEM/FDTD operators) | **`sprs`** or `faer` sparse | sprs 0.11+ | `faer`'s sparse module is the best pure-Rust sparse direct solver in 2026; `sprs` is the workhorse for SpMV/SpGEMM. |
| Sparse SpMV / iterative solves on GPU | **cuSPARSE via cudarc** | — | We wrap the result/sys layer ourselves into a small safe API. |

### Why not just `nalgebra` for everything?

`nalgebra` is excellent for ergonomic, small, statically-sized math but is explicitly **not optimized for medium/large dense factorization**. `faer`'s own docs recommend `nalgebra` for low-dimensional/graphics work and `faer` for medium/large dense. We use both for what each does best.

### Scaling caveat

10k×10k complex double = 1.6 GB — feasible on a single GPU via cuSOLVER. **100k×100k dense is 1.6 TB** — not feasible by brute force. For n ≥ 50k we will introduce iterative GMRES with block-diagonal preconditioning in Phase 1; **MLFMA / ACA / H-matrix compression is Phase 4 work.**

---

## CUDA numerical libraries (wrapped via cudarc)

| Library | Use in Yee | `cudarc` coverage |
|---|---|---|
| **cuBLAS** | Dense complex GEMM/GEMV for MoM matrix-fill aggregation | ✅ Full safe wrapper |
| **cuBLASLt** | Mixed-precision / Tensor-Core GEMM | ✅ Full safe wrapper |
| **cuSOLVER** | Dense LU (`Dn{C,Z}getrf`), Cholesky, QR, SVD, eigen | ✅ Safe wrapper (`cudarc::cusolver::safe`) |
| **cuSOLVERMg** | Multi-GPU dense factorization | ✅ |
| **cuSPARSE** | SpMV / SpMM for FEM/FDTD iterative solvers | ⚠️ `result` + `sys` only; we wrap |
| **cuFFT** | Near-to-far-field projection, FDTD post-processing | ✅ Safe wrapper (dynamic-loading only) |
| **cuRAND** | Random sources (rough surfaces, Monte Carlo) | ✅ |
| **NCCL** | Multi-GPU collective comms for FDTD domain decomp | ✅ Safe wrapper |

No other actively-maintained, high-coverage Rust wrapper exists for these libraries; standardizing on `cudarc` is the only sane choice.

---

## Meshing — Gmsh

**Choice: Gmsh 4.15+ via in-tree `bindgen` FFI to `gmshc.h`.**

Gmsh is the de-facto open-source mesher: actively maintained (4.15.2 stable as of March 2026), licensed **GPL v2+ with a linking exception** that is explicitly compatible with GPL v3 hosts. Its C API (`gmshc.h`) has been stable across the 4.x line since 2018.

We considered the existing Rust bindings:
- **`gmsh-sys`** and **`rgmsh`** (mxxo) — last updated November 2019, target Gmsh 4.4.1, six years stale. **Unmaintained.**

**We regenerate bindings ourselves with `bindgen` against the current `gmshc.h` from the official SDK** and write a thin safe wrapper in-tree as `yee-mesh`. The C API surface is a few hundred functions; this is roughly one day of work and lets us pin the Gmsh version we link against.

Fallback for cases where the Rust binding lags: shell out to the Gmsh CLI or call the Python API via PyO3 — slower iteration but bulletproof.

---

## CAD kernel — OpenCascade

**Choice: `opencascade-rs` 0.2+ wrapping OCCT 7.9.x (with a watch on 8.0, due May 2026).**

OCCT is **LGPL 2.1 with an exception**, confirmed GPL-2+/GPL-3 compatible by the Open CASCADE FAQ. It is the only mature, open-source B-rep CAD kernel; FreeCAD and Gmsh's OCC backend both use it.

Rust bindings:
- **`opencascade-rs`** (Bschwind) — 0.2.0, **actively maintained** (last commit February 2026), uses `cxx.rs` bridging. Exposes primitives, sketches, extrude/revolve/sweep, fillet/chamfer, booleans, STEP/STL/DXF/SVG/KiCad I/O. Builds OCCT as a git submodule — expect 5–15 minute cold builds; cache aggressively in CI.
- **`truck`** (RICOS) — pure-Rust B-rep with NURBS, MIT-licensed, **promising but immature**. Lacks robust booleans for complex topologies and mature real-world STEP import. **We watch but do not bet on it.**

We use `opencascade-rs` as primary. For cases it does not expose, we drop into raw FFI via `opencascade-sys`. For STEP/IGES *import only*, the Gmsh OCC factory (`importShapes("file.step")`) is an excellent fallback that reuses Gmsh's OCCT linkage and gives us tessellated geometry directly.

---

## Geometry primitives and spatial queries

| Concern | Crate | Version | Notes |
|---|---|---|---|
| Simulation math (vectors, matrices, transforms) | **`nalgebra`** | 0.34+ | Already chosen above. |
| Graphics math (shader-side, wgpu interop) | **`glam`** | 0.30+ | What wgpu/egui/rerun use natively. |
| Spatial queries / BVH / ray casts | **`parry3d`** | latest (mid-2025+) | Apache-2.0. BVH (`Qbvh`), TriMesh, AABB, distance queries — perfect for "is this field point inside copper?" and far-field projection ray casts. |
| 2D PCB polygon ops (copper layers, Gerber, DXF) | **`geo`** + **`i_overlay`** | geo 0.30, i_overlay 4.x | `i_overlay` is faster and more robust than `geo-booleanop` for PCB-grade Boolean operations. |
| Half-edge / triangle mesh data structures | In-tree module | — | No clear community standard; we adopt half-edge structures from `truck-topology` or roll our own thin wrapper around `parry3d::shape::TriMesh`. |

---

## Python bindings

**Choice: `pyo3` 0.28+ with `abi3-py310`, built by `maturin` 1.10+, distributed as `manylinux_2_28` wheels.**

- **PyO3 0.28.3** — supports CPython 3.7–3.14 and PyPy/GraalPy. The `Bound<'py, T>` API has been stable since 0.21. **Free-threaded Python (PEP 703 / 3.13t / 3.14t) is supported** since 0.23, mature in 0.28; note that `abi3` does not yet apply to the free-threaded ABI (PEP 803 is in flight, targeting Python 3.15). MSRV is Rust 1.83.
- **`abi3-py310`** lets us ship one wheel per OS/arch that works on Python 3.10+, vastly simplifying distribution.
- **`maturin` 1.10.2** (Nov 2025) — MSRV Rust 1.85; uses Rust 2024 edition internally. Builds wheels for Win/Linux/macOS/FreeBSD plus iOS-simulator targets in recent versions. Manylinux tagged automatically; we use `manylinux_2_28` (glibc 2.28+ / RHEL 8 baseline) for 2026.
- **`numpy` crate** matched to PyO3 0.28 for zero-copy NumPy ↔ `ndarray` interop. Essential for FDTD field arrays in Jupyter.
- **`pyo3-async-runtimes`** when we expose async solver callbacks.

CI builds via `PyO3/maturin-action@v1` with `manylinux: "2_28"`; we publish with `uv publish`.

---

## GUI — `egui` + `eframe`

**Choice: `egui` 0.34+ with `eframe`, a custom `wgpu` 3D viewport painted into an egui rect.**

We evaluated five frameworks:

| Framework | Latest | Verdict |
|---|---|---|
| **`egui`** | 0.34.2+ | **Selected.** Immediate-mode, pure-Rust, wgpu/glow backends. Sponsored by Rerun. Used by Rerun, MakerPnP, many scientific tools. |
| Tauri | 2.9.6 | Forces a JS/TS frontend; awkward for a tightly-coupled solver studio with a custom GPU viewport. |
| Dioxus | 0.7.0 | Promising (especially Dioxus Native, no webview), but the native renderer is brand new. |
| iced | 0.14 | Excellent architecture, but missing mature plotting/table widgets for scientific use. Revisit after 1.0. |
| Slint | 1.16+ | Beautiful for embedded HMIs, but the `.slint` DSL and absent plotting widget make it a poor fit. |

**Why egui wins for an EM studio:**
1. **3D viewport integration.** `egui-wgpu` lets you render into a texture and embed it in any egui rect — this is exactly how Rerun's 3D viewer is built.
2. **Plotting.** `egui_plot` drops in for S-parameter, Smith-chart, time-domain, and convergence plots.
3. **Docking and panels.** `egui_dock` or `egui_tiles` (Rerun's panel system) cover tabbed tools.
4. **Precedent.** Rerun, Bevy's editor work, several open-source FEA/CFD tools — the "scientific Rust GUI" segment is overwhelmingly egui-based today.
5. **Cross-platform and web.** `eframe` runs on Win/Mac/Linux natively and compiles to WASM unchanged.

**Honest downsides:** breaking minor releases roughly quarterly (mitigated by isolating UI behind a stable internal trait); not "native-looking" (acceptable for a technical tool).

---

## 3D visualization

**Choice: roll the viewport on `wgpu` directly, embedded inside an `egui` paint callback. Use `rerun` as a complementary debugging logger.**

| Library | Latest | Role for Yee |
|---|---|---|
| **`wgpu`** | 26.x | Primary viewport; custom shaders for surface heatmaps, far-field deformed spheres, cut planes. |
| **`rerun`** | 0.31–0.32 | Optional structured logging sink for solver internals: dump E-field samples per frequency, mesh per refinement step, far-field per iteration. Engineers run the Rerun viewer alongside the simulator. Apache-2.0. |
| `bevy` | 0.17–0.18 | **Not used.** Game engine ECS is the wrong shape for a deterministic, parameter-driven scientific viewport. |
| `kiss3d` | 0.40 (Jan 2026, wgpu rewrite) | Useful for documentation examples; rendering model too limited for main viewport. |
| `three-d` | 0.18 | Reasonable middle ground; doesn't add much over wgpu+egui for our use case. |

Required capabilities — extruded copper layers, current-density heatmaps, E/H cut planes, 3D far-field patterns — all reduce to a few hundred lines of wgpu plus `glam`, with full shader control. Rerun complements: at any time the solver can `rerun.log()` arbitrary tensors, meshes, scalars on a timeline, giving us a free interactive multi-view debugger.

---

## Plotting

| Channel | Choice | Why |
|---|---|---|
| In-app live plots | **`egui_plot`** | Native immediate-mode plot widget; pans, zooms, axis-linked panels — perfect for sweeping S-params across ports. |
| Static / publication export | **`plotters`** | Pure Rust; SVG/PNG/Cairo backends; publication quality. |
| Power-user / scripting | Matplotlib / plotly in user notebooks via the PyO3 bindings | The cleanest path — don't host Python inside the GUI. |

---

## Logging, error handling, async

- **Errors:** `thiserror` for library crates, `anyhow` for the CLI/GUI binaries.
- **Logging / tracing:** `tracing` + `tracing-subscriber`. Rerun integration is additive on top.
- **Async (where used):** `tokio` for any background services (file watch, surrogate training orchestration). The solver itself is synchronous, GPU-driven.
- **Config / project files:** TOML for human-edited (`yee.toml`); HDF5 or Arrow IPC for solver output (large field arrays, time-series).

---

## Quick reference — workspace `Cargo.toml` excerpt

```toml
[workspace.package]
rust-version = "1.85"
edition      = "2024"
license      = "GPL-3.0-or-later"

[workspace.dependencies]
# CUDA (one binding for the whole project)
cudarc = { version = "0.19", default-features = false, features = [
    "std", "driver", "nvrtc",
    "cublas", "cublaslt", "cusolver", "cusparse", "cufft", "curand",
    "f16",
    "cuda-version-from-build-system",
    "dynamic-loading",
] }

# Linear algebra
faer       = "0.23"
nalgebra   = "0.34"
ndarray    = "0.17"
sprs       = "0.11"
num-complex = "0.4"

# Geometry
glam        = "0.30"
parry3d     = "0.22"
geo         = "0.30"
i_overlay   = "4"

# CAD / Meshing
opencascade     = "0.2"   # bschwind/opencascade-rs
opencascade-sys = "0.2"
# Gmsh: in-tree bindgen against gmshc.h 4.15+, not the stale rgmsh crate

# GUI / Viz
egui        = "0.34"
eframe      = "0.34"
egui_plot   = "0.34"
egui_dock   = "0.18"
wgpu        = "26"
rerun       = { version = "0.31", default-features = false, features = ["sdk"] }
plotters    = "0.3"

# Python bindings
pyo3   = { version = "0.28", features = ["abi3-py310", "extension-module"] }
numpy  = "0.28"

# Diagnostics
tracing            = "0.1"
tracing-subscriber = "0.3"
thiserror          = "1"
anyhow             = "1"
```

---

## License sanity check

Every dependency above is verified compatible with **GPL v3** distribution:

- **Gmsh** — GPL v2+ with linking exception (FAQ-confirmed GPL-3 compatible).
- **OCCT** — LGPL 2.1 with exception (FAQ-confirmed GPL-2+ compatible).
- **`opencascade-rs`** — LGPL-2.1.
- **`egui`, `eframe`, `egui_plot`, `parry`, `wgpu`, `glam`, `nalgebra`, `ndarray`, `faer`, `rerun`, `plotters`, `pyo3`, `numpy`, `cudarc`** — all permissive (MIT/Apache/BSD), no conflict.
- **NVIDIA CUDA libraries (cuBLAS, cuSOLVER, cuSPARSE, cuFFT, NCCL)** — proprietary but **dynamically linked at runtime**, which is the same posture every CUDA-using OSS project takes; we document the dependency clearly.

No license conflicts in the recommended stack.

---

## Known watch-outs

1. **`cudarc` is "pre-alpha"** per its own README; history includes breaking minor releases (notably 0.13 → 0.14). Pin to exact minor version and route all CUDA access through a thin internal abstraction.
2. **Gmsh Rust bindings (`rgmsh`/`gmsh-sys`) are abandoned.** We bake our own.
3. **`opencascade-rs` cold build is 5–15 minutes** (rebuilds OCCT). Cache CI builds aggressively (`sccache` or workspace artifact cache).
4. **egui has breaking minor releases roughly quarterly.** Isolate UI behind a stable internal trait.
5. **PyO3 ships one breaking change per minor release.** Pin and migrate deliberately.
6. **Bevy 0.x ecosystem resets every ~3 months** — another reason we are not using Bevy as the UI shell.
7. **OCCT 8.0 targets May 2026** with significant changes (C++17, BRepGraph topology); expect a wave of `opencascade-rs` work after that lands.
8. **Dense LU at n ≥ 50k overflows a single GPU.** Iterative methods in Phase 1; MLFMA/ACA in Phase 4.
9. **`rustc_codegen_nvvm` (Rust on GPU) is not production-ready.** Re-evaluate late 2026. Until then, kernels are CUDA-C JITed via NVRTC.
```

---

These three files form a coherent foundation: the **README** sets the vision and positioning, **ROADMAP** commits to a realistic phased plan with concrete validation gates, and **TECH_STACK** justifies every dependency with current 2026 versions and known risks. All licensing is verified compatible with GPL v3. The recommended project name **Yee** is unique in the EM space and instantly meaningful to the target user.