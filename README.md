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

## Performance

FDTD throughput, vacuum cubic grids, FP32, sync-correct (device-wait timed, not submit-only), 3-repetition medians:

| Grid | Mcells·step/s |
|---|---|
| 64³  | 11537 |
| 96³  | 12665 |
| 128³ | 4442 |
| 160³ | 4519 |
| 192³ | 4685 |
| 224³ | 4685 |

Measured on a single RTX 5060 Ti 16 GB. These are hardware-specific numbers from one consumer card, not a cross-vendor or "fastest" claim. They clear the gprMax CUDA-FDTD reference bar (3405 Mcells/s) at every grid size measured, peaking at 3.7× the bar. See [ADR-0223](docs/src/decisions/0223-fs70-gpu-bench.md) (bench methodology, first honest bar miss) and [ADR-0224](docs/src/decisions/0224-fs71-kernel-opt.md) (kernel-fusion optimization that cleared it).

Reproduce locally:

```bash
cargo run -p yee-compute --release --example bench
```

The GPU nightly workflow (`.github/workflows/gpu-nightly.yml`) re-runs this bench and publishes fresh JSON + table artifacts on every scheduled run, but only when a GPU runner is registered (repo variable `YEE_GPU_RUNNER_ENABLED`) — without one, the job no-ops and no numbers are published.

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
